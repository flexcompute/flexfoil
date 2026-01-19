//! Global Newton-Raphson viscous-inviscid coupling solver.
//!
//! This module implements XFOIL-style simultaneous solution of the inviscid
//! panel equations and boundary layer equations using Newton's method.
//!
//! The key advantage over sequential iteration is quadratic convergence
//! and better handling of separated flows.
//!
//! # State Vector
//!
//! The solver solves for:
//! - gamma[0..N]: Panel vorticity at N nodes
//! - theta[0..M]: Momentum thickness at M BL stations
//! - h[0..M]: Shape factor H at M BL stations
//!
//! # Residual Equations
//!
//! - R_inviscid: Panel stream function BC with transpiration
//! - R_momentum: von Kármán momentum integral equation
//! - R_shape: Shape factor closure (Head or Cτ lag)
//!
//! # Reference
//!
//! Drela & Giles (1987) "Viscous-Inviscid Analysis of Transonic and Low
//! Reynolds Number Airfoils", AIAA J. 25(10)

use crate::boundary_layer::{BLConfig, TurbulentModel};
use crate::boundary_layer::closure::{ludwieg_tillmann_cf, head_h1, head_entrainment};
use crate::inviscid::{FlowConditions, FactorizedSolution, InviscidSolution};
use crate::viscous::compute_transpiration;
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

/// Configuration for the Newton VII solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NewtonConfig {
    /// Reynolds number
    pub reynolds: f64,
    /// Critical N-factor for transition
    pub n_crit: f64,
    /// Turbulent model
    pub turbulent_model: TurbulentModel,
    /// Maximum Newton iterations
    pub max_iterations: usize,
    /// Convergence tolerance (relative residual norm)
    pub tolerance: f64,
    /// Finite difference step for Jacobian
    pub fd_epsilon: f64,
    /// Line search parameters
    pub line_search_max_iter: usize,
    /// Armijo condition parameter
    pub armijo_c: f64,
    /// Initial damping factor (0 < damping <= 1)
    /// Start conservative and increase on success
    pub initial_damping: f64,
    /// Minimum damping factor
    pub min_damping: f64,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            reynolds: 1e6,
            n_crit: 9.0,
            turbulent_model: TurbulentModel::default(),
            max_iterations: 50,
            tolerance: 1e-6,
            fd_epsilon: 1e-7,
            line_search_max_iter: 10,
            armijo_c: 1e-4,
            initial_damping: 0.3,  // Start conservative
            min_damping: 0.05,     // Minimum step size
        }
    }
}

impl NewtonConfig {
    /// Create config from BLConfig
    pub fn from_bl_config(bl_config: &BLConfig) -> Self {
        Self {
            reynolds: bl_config.reynolds,
            n_crit: bl_config.n_crit,
            turbulent_model: bl_config.turbulent_model,
            ..Default::default()
        }
    }
}

/// Geometry information needed for Newton solver.
#[derive(Debug, Clone)]
pub struct NewtonGeometry {
    /// Arc-length coordinates at each panel node
    pub s_coords: Vec<f64>,
    /// X-coordinates at each panel node
    pub x_coords: Vec<f64>,
    /// Y-coordinates at each panel node
    pub y_coords: Vec<f64>,
    /// Chord length
    pub chord: f64,
    /// Stagnation point index
    pub stag_idx: usize,
    /// Number of panel nodes
    pub n_panels: usize,
    /// Upper surface indices (from stag to TE)
    pub upper_indices: Vec<usize>,
    /// Lower surface indices (from stag to TE)
    pub lower_indices: Vec<usize>,
}

impl NewtonGeometry {
    /// Create geometry from coordinates.
    pub fn new(
        s_coords: Vec<f64>,
        x_coords: Vec<f64>,
        y_coords: Vec<f64>,
        chord: f64,
        stag_idx: usize,
    ) -> Self {
        let n_panels = s_coords.len();
        
        // Upper surface: from stag_idx down to 0 (marching from LE to TE)
        let upper_indices: Vec<usize> = (0..=stag_idx).rev().collect();
        
        // Lower surface: from stag_idx to n-1 (marching from LE to TE)
        // Note: stag_idx appears in both surfaces (like in XFOIL)
        let lower_indices: Vec<usize> = (stag_idx..n_panels).collect();
        
        Self {
            s_coords,
            x_coords,
            y_coords,
            chord,
            stag_idx,
            n_panels,
            upper_indices,
            lower_indices,
        }
    }
    
    /// Get total number of BL stations (both surfaces).
    /// This equals n_panels + 1 because stagnation point is in both surfaces.
    pub fn n_bl_stations(&self) -> usize {
        self.upper_indices.len() + self.lower_indices.len()
    }
}

/// State vector for Newton iteration.
///
/// Contains all unknowns that are solved simultaneously.
/// 
/// # XFOIL-Style Variables
/// 
/// The BL unknowns per station are `[Ctau/Ampl, theta, mass]` where:
/// - `mass = Ue * delta_star = Ue * theta * H`
/// - `Ue` is NOT an unknown - it's computed from: `Ue = Ue_inviscid + sum_j(DIJ[i][j] * mass[j])`
/// - `H` is a derived quantity: `H = mass / (Ue * theta)`
/// 
/// This matches XFOIL's state vector structure for proper V-I coupling.
#[derive(Debug, Clone)]
pub struct NewtonState {
    /// Panel vorticity at each node (N values)
    pub gamma: Vec<f64>,
    /// BL momentum thickness at each station
    /// Ordered: upper surface (LE to TE), lower surface (LE to TE)
    pub theta: Vec<f64>,
    /// BL mass defect at each station: mass = Ue * delta_star
    /// This is the XFOIL-style 3rd variable that couples to inviscid through DIJ
    pub mass: Vec<f64>,
    /// Whether each station is turbulent
    pub is_turbulent: Vec<bool>,
    /// Transition index on upper surface (index into upper stations)
    pub transition_upper: Option<usize>,
    /// Transition index on lower surface (index into lower stations)
    pub transition_lower: Option<usize>,
    /// Ctau (shear stress coefficient) for turbulent stations
    /// In laminar regions, this stores amplification factor
    pub ctau: Vec<f64>,
    /// Whether station is in inverse mode (separated)
    pub inverse_mode: Vec<bool>,
    /// Target Hk for inverse mode stations
    pub hk_target: Vec<f64>,
}

impl NewtonState {
    /// Create a new state with given sizes.
    pub fn new(n_panels: usize, n_bl_stations: usize) -> Self {
        // Initialize with reasonable defaults
        // Initial theta from Blasius: θ ~ 0.664*sqrt(ν*x/U) ~ 1e-4 for Re~1e6
        let theta_init = 1e-5;
        // Initial H ~ 2.6 (Blasius), so mass = Ue * θ * H ~ 1.0 * 1e-5 * 2.6
        let mass_init = theta_init * 2.6;
        
        Self {
            gamma: vec![0.0; n_panels],
            theta: vec![theta_init; n_bl_stations],
            mass: vec![mass_init; n_bl_stations],
            is_turbulent: vec![false; n_bl_stations],
            transition_upper: None,
            transition_lower: None,
            ctau: vec![0.03; n_bl_stations], // Typical initial Ctau
            inverse_mode: vec![false; n_bl_stations],
            hk_target: vec![2.5; n_bl_stations],
        }
    }
    
    /// Initialize from inviscid solution.
    pub fn from_inviscid(inviscid: &InviscidSolution, n_bl_stations: usize) -> Self {
        let n_panels = inviscid.gamma.len();
        let mut state = Self::new(n_panels, n_bl_stations);
        state.gamma = inviscid.gamma.clone();
        state
    }
    
    /// Total number of unknowns: gamma + theta + mass
    pub fn total_size(&self) -> usize {
        self.gamma.len() + self.theta.len() + self.mass.len()
    }
    
    /// Get value at flat index.
    pub fn get(&self, idx: usize) -> f64 {
        let n_gamma = self.gamma.len();
        let n_theta = self.theta.len();
        
        if idx < n_gamma {
            self.gamma[idx]
        } else if idx < n_gamma + n_theta {
            self.theta[idx - n_gamma]
        } else {
            self.mass[idx - n_gamma - n_theta]
        }
    }
    
    /// Set value at flat index.
    pub fn set(&mut self, idx: usize, value: f64) {
        let n_gamma = self.gamma.len();
        let n_theta = self.theta.len();
        
        if idx < n_gamma {
            self.gamma[idx] = value;
        } else if idx < n_gamma + n_theta {
            self.theta[idx - n_gamma] = value;
        } else {
            self.mass[idx - n_gamma - n_theta] = value;
        }
    }
    
    /// Perturb a single variable by epsilon.
    pub fn perturb(&mut self, idx: usize, eps: f64) {
        let current = self.get(idx);
        self.set(idx, current + eps);
    }
    
    /// Convert to flat vector.
    pub fn to_vec(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.total_size());
        v.extend(&self.gamma);
        v.extend(&self.theta);
        v.extend(&self.mass);
        v
    }
    
    /// Update from flat vector.
    pub fn from_vec(&mut self, v: &[f64]) {
        let n_gamma = self.gamma.len();
        let n_theta = self.theta.len();
        
        self.gamma.copy_from_slice(&v[0..n_gamma]);
        self.theta.copy_from_slice(&v[n_gamma..n_gamma + n_theta]);
        self.mass.copy_from_slice(&v[n_gamma + n_theta..]);
    }
    
    /// Apply a step with given scale factor.
    pub fn apply_step(&self, step: &NewtonStep, alpha: f64) -> Self {
        let mut new_state = self.clone();
        
        for i in 0..self.gamma.len() {
            new_state.gamma[i] += alpha * step.d_gamma[i];
        }
        for i in 0..self.theta.len() {
            new_state.theta[i] = (self.theta[i] + alpha * step.d_theta[i]).max(1e-10);
        }
        for i in 0..self.mass.len() {
            // Mass must stay positive
            new_state.mass[i] = (self.mass[i] + alpha * step.d_mass[i]).max(1e-12);
        }
        
        // Copy auxiliary state
        new_state.ctau = self.ctau.clone();
        new_state.inverse_mode = self.inverse_mode.clone();
        new_state.hk_target = self.hk_target.clone();
        
        new_state
    }
    
    /// Compute displacement thickness: δ* = mass / Ue
    /// 
    /// Since mass = Ue * δ*, we have δ* = mass / Ue
    pub fn delta_star(&self) -> Vec<f64> {
        // For now, return mass as δ* approximation (assuming Ue ~ 1)
        // The caller should use compute_delta_star_with_ue for accuracy
        self.mass.clone()
    }
    
    /// Compute displacement thickness given edge velocities.
    /// δ* = mass / Ue
    pub fn delta_star_with_ue(&self, ue: &[f64]) -> Vec<f64> {
        self.mass.iter()
            .zip(ue.iter())
            .map(|(&m, &u)| m / u.max(1e-10))
            .collect()
    }
    
    /// Compute shape factor H given edge velocities.
    /// H = δ* / θ = mass / (Ue * θ)
    pub fn h_with_ue(&self, ue: &[f64]) -> Vec<f64> {
        self.mass.iter()
            .zip(self.theta.iter())
            .zip(ue.iter())
            .map(|((&m, &th), &u)| m / (u.max(1e-10) * th.max(1e-10)))
            .collect()
    }
    
    /// Compute edge velocity from inviscid Ue plus DIJ mass coupling.
    /// 
    /// This is the key XFOIL relationship:
    /// `Ue[i] = Ue_inviscid[i] + sum_j(DIJ[i][j] * mass[j])`
    /// 
    /// # Arguments
    /// * `ue_inviscid` - Edge velocity from inviscid solution (|gamma|)
    /// * `dij` - Source influence matrix
    /// * `bl_to_panel` - Mapping from BL station index to panel index
    /// 
    /// # Returns
    /// Edge velocity at each panel node
    pub fn compute_ue_from_mass(
        &self,
        ue_inviscid: &[f64],
        dij: &[Vec<f64>],
        bl_to_panel: &[usize],
    ) -> Vec<f64> {
        let n_panels = ue_inviscid.len();
        let mut ue = ue_inviscid.to_vec();
        
        // Build mass at each panel from BL stations
        let mut mass_at_panel = vec![0.0; n_panels];
        for (bl_idx, &panel_idx) in bl_to_panel.iter().enumerate() {
            if panel_idx < n_panels && bl_idx < self.mass.len() {
                mass_at_panel[panel_idx] = self.mass[bl_idx];
            }
        }
        
        // Apply DIJ coupling: dUe[i] = sum_j DIJ[i][j] * mass[j]
        for i in 0..n_panels {
            if i < dij.len() {
                for j in 0..n_panels.min(dij[i].len()) {
                    ue[i] += dij[i][j] * mass_at_panel[j];
                }
            }
        }
        
        // Ensure positive
        for u in &mut ue {
            *u = u.abs().max(1e-10);
        }
        
        ue
    }
    
    /// Initialize mass from H and Ue: mass = Ue * θ * H
    pub fn init_mass_from_h(&mut self, h: &[f64], ue: &[f64], bl_to_panel: &[usize]) {
        for (bl_idx, &panel_idx) in bl_to_panel.iter().enumerate() {
            if bl_idx < self.mass.len() && bl_idx < h.len() {
                let ue_val = if panel_idx < ue.len() { ue[panel_idx].abs().max(1e-10) } else { 1.0 };
                self.mass[bl_idx] = ue_val * self.theta[bl_idx] * h[bl_idx];
            }
        }
    }
}

/// Newton step (dx in J*dx = -R).
#[derive(Debug, Clone)]
pub struct NewtonStep {
    /// Change in gamma
    pub d_gamma: Vec<f64>,
    /// Change in theta
    pub d_theta: Vec<f64>,
    /// Change in mass defect
    pub d_mass: Vec<f64>,
}

impl NewtonStep {
    /// Create from flat vector.
    pub fn from_vec(v: &[f64], n_gamma: usize, n_bl: usize) -> Self {
        Self {
            d_gamma: v[0..n_gamma].to_vec(),
            d_theta: v[n_gamma..n_gamma + n_bl].to_vec(),
            d_mass: v[n_gamma + n_bl..].to_vec(),
        }
    }
    
    /// Create a zero step.
    pub fn zero(n_gamma: usize, n_bl: usize) -> Self {
        Self {
            d_gamma: vec![0.0; n_gamma],
            d_theta: vec![0.0; n_bl],
            d_mass: vec![0.0; n_bl],
        }
    }
    
    /// Compute L2 norm.
    pub fn norm(&self) -> f64 {
        let sum: f64 = self.d_gamma.iter().map(|x| x * x).sum::<f64>()
            + self.d_theta.iter().map(|x| x * x).sum::<f64>()
            + self.d_mass.iter().map(|x| x * x).sum::<f64>();
        sum.sqrt()
    }
}

/// Residual evaluation result.
#[derive(Debug, Clone)]
pub struct Residuals {
    /// Inviscid residuals (N values) - panel BC with transpiration
    pub r_inviscid: Vec<f64>,
    /// Momentum equation residuals (M values)
    pub r_momentum: Vec<f64>,
    /// Shape equation residuals (M values)
    pub r_shape: Vec<f64>,
}

impl Residuals {
    /// Create new residuals with given sizes.
    pub fn new(n_panels: usize, n_bl: usize) -> Self {
        Self {
            r_inviscid: vec![0.0; n_panels],
            r_momentum: vec![0.0; n_bl],
            r_shape: vec![0.0; n_bl],
        }
    }
    
    /// Total size.
    pub fn total_size(&self) -> usize {
        self.r_inviscid.len() + self.r_momentum.len() + self.r_shape.len()
    }
    
    /// Get value at flat index.
    pub fn get(&self, idx: usize) -> f64 {
        let n_inv = self.r_inviscid.len();
        let n_mom = self.r_momentum.len();
        
        if idx < n_inv {
            self.r_inviscid[idx]
        } else if idx < n_inv + n_mom {
            self.r_momentum[idx - n_inv]
        } else {
            self.r_shape[idx - n_inv - n_mom]
        }
    }
    
    /// Convert to flat vector.
    pub fn to_vec(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.total_size());
        v.extend(&self.r_inviscid);
        v.extend(&self.r_momentum);
        v.extend(&self.r_shape);
        v
    }
    
    /// Compute L2 norm.
    pub fn norm(&self) -> f64 {
        let sum: f64 = self.r_inviscid.iter().map(|x| x * x).sum::<f64>()
            + self.r_momentum.iter().map(|x| x * x).sum::<f64>()
            + self.r_shape.iter().map(|x| x * x).sum::<f64>();
        sum.sqrt()
    }
    
    /// Compute max absolute residual.
    pub fn max_abs(&self) -> f64 {
        let max_inv = self.r_inviscid.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let max_mom = self.r_momentum.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        let max_shape = self.r_shape.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        max_inv.max(max_mom).max(max_shape)
    }
}

/// Result of Newton solver iteration.
#[derive(Debug, Clone)]
pub struct NewtonResult {
    /// Final state
    pub state: NewtonState,
    /// Final residual norm
    pub residual_norm: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Whether converged
    pub converged: bool,
    /// Residual history
    pub residual_history: Vec<f64>,
}

// ============================================================================
// XFOIL-Style Block-Tridiagonal Structures
// ============================================================================

/// 3x3 block for BL Jacobian (matches XFOIL VA/VB blocks).
/// 
/// Variables are: [0] = Ctau/Ampl, [1] = Theta, [2] = mass defect m = Ue*delta*
#[derive(Debug, Clone, Copy, Default)]
pub struct BLBlock {
    /// Block data [row][col]
    pub data: [[f64; 3]; 3],
}

impl BLBlock {
    /// Create a zero block.
    pub fn zero() -> Self {
        Self { data: [[0.0; 3]; 3] }
    }
    
    /// Create an identity block.
    pub fn identity() -> Self {
        Self {
            data: [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        }
    }
    
    /// Access element (i, j).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i][j]
    }
    
    /// Set element (i, j).
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        self.data[i][j] = value;
    }
    
    /// Multiply block by vector: result = A * v
    pub fn mul_vec(&self, v: &[f64; 3]) -> [f64; 3] {
        [
            self.data[0][0] * v[0] + self.data[0][1] * v[1] + self.data[0][2] * v[2],
            self.data[1][0] * v[0] + self.data[1][1] * v[1] + self.data[1][2] * v[2],
            self.data[2][0] * v[0] + self.data[2][1] * v[1] + self.data[2][2] * v[2],
        ]
    }
    
    /// Solve 3x3 system: A * x = b using Gaussian elimination with partial pivoting.
    /// Returns None if singular.
    pub fn solve(&self, b: &[f64; 3]) -> Option<[f64; 3]> {
        let mut a = self.data;
        let mut x = *b;
        
        // Forward elimination with partial pivoting
        for k in 0..3 {
            // Find pivot
            let mut max_idx = k;
            let mut max_val = a[k][k].abs();
            for i in k + 1..3 {
                if a[i][k].abs() > max_val {
                    max_val = a[i][k].abs();
                    max_idx = i;
                }
            }
            
            if max_val < 1e-15 {
                return None; // Singular
            }
            
            // Swap rows
            if max_idx != k {
                a.swap(k, max_idx);
                x.swap(k, max_idx);
            }
            
            // Eliminate
            let pivot = a[k][k];
            for i in k + 1..3 {
                let factor = a[i][k] / pivot;
                for j in k + 1..3 {
                    a[i][j] -= factor * a[k][j];
                }
                x[i] -= factor * x[k];
            }
        }
        
        // Back substitution
        for i in (0..3).rev() {
            for j in i + 1..3 {
                x[i] -= a[i][j] * x[j];
            }
            x[i] /= a[i][i];
        }
        
        Some(x)
    }
    
    /// Compute Frobenius norm of block.
    pub fn norm(&self) -> f64 {
        let mut sum = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                sum += self.data[i][j] * self.data[i][j];
            }
        }
        sum.sqrt()
    }
}

/// Scaling factors for Newton variables to improve conditioning.
/// 
/// XFOIL uses mass defect m = Ue*delta* instead of delta* directly,
/// which varies more smoothly and improves the Jacobian condition number.
/// 
/// The scaled system is: D_r * J * D_v^{-1} * (D_v * dx) = D_r * R
/// where D_r scales residuals and D_v scales variables.
#[derive(Debug, Clone)]
pub struct NewtonScaling {
    /// Scale for first BL variable (Ctau for turbulent, Ampl for laminar)
    pub s_ctau: f64,
    /// Scale for momentum thickness theta
    pub s_theta: f64,
    /// Scale for mass defect m = Ue*delta*
    pub s_mass: f64,
    /// Reference edge velocity for scaling
    pub ue_ref: f64,
    /// Reference length (chord)
    pub length_ref: f64,
    /// Residual scale for first equation (shear lag / amplification)
    pub r_eq1: f64,
    /// Residual scale for momentum equation
    pub r_mom: f64,
    /// Residual scale for shape equation
    pub r_shape: f64,
}

impl NewtonScaling {
    /// Create scaling factors from flow conditions.
    pub fn from_flow(reynolds: f64, chord: f64, ue_max: f64) -> Self {
        // Typical BL thickness ~ chord / sqrt(Re)
        let bl_scale = chord / reynolds.sqrt();
        
        Self {
            s_ctau: 1.0,                    // Ctau ~ O(1)
            s_theta: bl_scale,               // theta ~ chord/sqrt(Re)
            s_mass: ue_max * bl_scale,       // m ~ Ue * delta*
            ue_ref: ue_max.max(1.0),
            length_ref: chord,
            // Residual scales (targeting unit magnitude)
            r_eq1: 1.0,              // Shear lag eqn residual
            r_mom: bl_scale,          // Momentum eqn has theta terms
            r_shape: 1.0,            // Shape eqn is ~O(1)
        }
    }
    
    /// Default scaling (no scaling applied).
    pub fn none() -> Self {
        Self {
            s_ctau: 1.0,
            s_theta: 1.0,
            s_mass: 1.0,
            ue_ref: 1.0,
            length_ref: 1.0,
            r_eq1: 1.0,
            r_mom: 1.0,
            r_shape: 1.0,
        }
    }
    
    /// Scale a station state vector [ctau, theta, mass] for Newton iteration.
    pub fn scale(&self, v: &mut [f64; 3]) {
        v[0] /= self.s_ctau;
        v[1] /= self.s_theta;
        v[2] /= self.s_mass;
    }
    
    /// Unscale a station state vector.
    pub fn unscale(&self, v: &mut [f64; 3]) {
        v[0] *= self.s_ctau;
        v[1] *= self.s_theta;
        v[2] *= self.s_mass;
    }
    
    /// Scale residual vector [eq1, momentum, shape].
    pub fn scale_residual(&self, r: &mut [f64; 3]) {
        r[0] /= self.r_eq1;
        r[1] /= self.r_mom;
        r[2] /= self.r_shape;
    }
    
    /// Scale a Jacobian block: J_scaled = D_r * J * D_v^{-1}
    pub fn scale_block(&self, block: &mut BLBlock) {
        // Row scaling (residuals)
        let r_scales = [self.r_eq1, self.r_mom, self.r_shape];
        // Column scaling (variables) - inverse
        let v_inv_scales = [self.s_ctau, self.s_theta, self.s_mass];
        
        for i in 0..3 {
            for j in 0..3 {
                block.data[i][j] *= v_inv_scales[j] / r_scales[i];
            }
        }
    }
    
    /// Apply scaling to entire block-tridiagonal system.
    pub fn scale_system(&self, jacobian: &mut BlockTridiagJacobian) {
        for iv in 0..jacobian.n_stations {
            // Scale diagonal block
            self.scale_block(&mut jacobian.diag[iv]);
            
            // Scale sub-diagonal block
            if iv > 0 {
                self.scale_block(&mut jacobian.sub_diag[iv]);
            }
            
            // Scale RHS
            self.scale_residual(&mut jacobian.rhs[iv]);
        }
    }
    
    /// Unscale solution from block system.
    pub fn unscale_solution(&self, jacobian: &mut BlockTridiagJacobian) {
        for iv in 0..jacobian.n_stations {
            self.unscale(&mut jacobian.solution[iv]);
        }
    }
}

/// VACCEL acceleration parameters for BLSOLV (from XFOIL).
/// 
/// When |VM[k,iv]| < vacc, skip the elimination step.
/// This speeds up the solve and often improves stability.
#[derive(Debug, Clone)]
pub struct VaccelParams {
    /// Base acceleration parameter (XFOIL default: 0.01)
    pub vaccel: f64,
    /// Threshold for row 1 (Ctau) elimination
    pub vacc1: f64,
    /// Threshold for row 2 (theta) elimination
    pub vacc2: f64,
    /// Threshold for row 3 (mass) elimination
    pub vacc3: f64,
}

impl Default for VaccelParams {
    fn default() -> Self {
        Self {
            vaccel: 0.01,
            vacc1: 0.01,
            vacc2: 0.01,
            vacc3: 0.01,
        }
    }
}

impl VaccelParams {
    /// Create params scaled by arc length (as in XFOIL BLSOLV).
    pub fn from_arc_length(vaccel: f64, s_total: f64) -> Self {
        let vacc_scaled = vaccel * 2.0 / s_total.max(1e-10);
        Self {
            vaccel,
            vacc1: vaccel,
            vacc2: vacc_scaled,
            vacc3: vacc_scaled,
        }
    }
}

/// Block-tridiagonal Jacobian storage for BL Newton system.
/// 
/// Matches XFOIL's BLSOLV structure:
/// ```text
///     A  |  |  .  |  |  .  |    d       R
///     B  A  |  .  |  |  .  |    d       R
///     |  B  A  .  |  |  .  |    d   =   R
///     .  .  .  .  |  |  .  |    d       R
///     |  Z  |  |  B  A  .  |    d       R
///     |  |  |  |  |  |  B  A    d       R
/// ```
/// 
/// - A: 3x3 diagonal blocks (derivatives wrt current station)
/// - B: 3x3 sub-diagonal blocks (derivatives wrt previous station)
/// - Z: Wake jump block at TE
/// - VM: Mass influence matrix (sparse coupling to inviscid via Ue)
///   VM[i][j][k] = d(residual[i][k]) / d(Ue[j]) where k is equation index
#[derive(Debug, Clone)]
pub struct BlockTridiagJacobian {
    /// Number of BL stations
    pub n_stations: usize,
    /// Diagonal blocks VA[station] - 3x3 each
    pub diag: Vec<BLBlock>,
    /// Sub-diagonal blocks VB[station] - 3x3 each (VB[0] unused)
    pub sub_diag: Vec<BLBlock>,
    /// Mass influence: VM[iv][jv] = 3-element vector for Ue coupling
    /// d(res[iv][k]) / d(Ue[jv]) for k=0,1,2
    pub mass_influence: Vec<Vec<[f64; 3]>>,
    /// Right-hand side (residuals)
    pub rhs: Vec<[f64; 3]>,
    /// Solution vector
    pub solution: Vec<[f64; 3]>,
    /// Wake jump block VZ at TE (optional, for wake coupling)
    pub wake_jump: Option<[[f64; 3]; 2]>,
    /// Index of TE station on side 1
    pub te_idx_1: Option<usize>,
    /// VACCEL parameters
    pub vaccel: VaccelParams,
}

impl BlockTridiagJacobian {
    /// Create a new block-tridiagonal Jacobian for n_stations.
    pub fn new(n_stations: usize, _n_panels: usize) -> Self {
        Self {
            n_stations,
            diag: vec![BLBlock::zero(); n_stations],
            sub_diag: vec![BLBlock::zero(); n_stations],
            mass_influence: vec![vec![[0.0; 3]; n_stations]; n_stations],
            rhs: vec![[0.0; 3]; n_stations],
            solution: vec![[0.0; 3]; n_stations],
            wake_jump: None,
            te_idx_1: None,
            vaccel: VaccelParams::default(),
        }
    }
    
    /// Set diagonal block at station iv.
    pub fn set_diag(&mut self, iv: usize, block: BLBlock) {
        if iv < self.n_stations {
            self.diag[iv] = block;
        }
    }
    
    /// Set sub-diagonal block at station iv (relates iv to iv-1).
    pub fn set_subdiag(&mut self, iv: usize, block: BLBlock) {
        if iv < self.n_stations && iv > 0 {
            self.sub_diag[iv] = block;
        }
    }
    
    /// Set RHS residual at station iv.
    pub fn set_rhs(&mut self, iv: usize, residual: [f64; 3]) {
        if iv < self.n_stations {
            self.rhs[iv] = residual;
        }
    }
    
    /// Set mass influence: d(residual[iv][k]) / d(Ue[jv]) for each equation k.
    pub fn set_mass_influence_vec(&mut self, iv: usize, jv: usize, values: [f64; 3]) {
        if iv < self.n_stations && jv < self.mass_influence[iv].len() {
            self.mass_influence[iv][jv] = values;
        }
    }
    
    /// Set mass influence for a single equation.
    pub fn set_mass_influence(&mut self, iv: usize, jv: usize, eqn: usize, value: f64) {
        if iv < self.n_stations && jv < self.mass_influence[iv].len() && eqn < 3 {
            self.mass_influence[iv][jv][eqn] = value;
        }
    }
    
    /// Solve the block-tridiagonal system using XFOIL-style BLSOLV.
    /// 
    /// This is a custom solver that exploits the block structure:
    /// 1. Forward elimination with 3x3 block inversions
    /// 2. VACCEL skipping for small off-diagonal entries
    /// 3. Back substitution
    pub fn solve(&mut self) -> bool {
        if self.n_stations == 0 {
            return true;
        }
        
        // Copy RHS to solution (will be modified in place)
        for iv in 0..self.n_stations {
            self.solution[iv] = self.rhs[iv];
        }
        
        // Forward elimination
        for iv in 0..self.n_stations {
            // Invert diagonal block VA[iv]
            if !self.invert_diagonal_block(iv) {
                return false; // Singular
            }
            
            // Eliminate sub-diagonal VB[iv+1] if exists
            if iv + 1 < self.n_stations {
                self.eliminate_subdiag(iv);
            }
            
            // Handle wake jump block at TE
            if let (Some(te_idx), Some(_vz)) = (self.te_idx_1, &self.wake_jump) {
                if iv == te_idx {
                    self.eliminate_wake_jump(iv);
                }
            }
            
            // Eliminate lower VM column entries (with VACCEL skipping)
            if iv + 2 < self.n_stations {
                self.eliminate_mass_influence_column(iv);
            }
        }
        
        // Back substitution
        for iv in (1..self.n_stations).rev() {
            self.back_substitute_column(iv);
        }
        
        true
    }
    
    /// Invert diagonal block at station iv (in-place Gauss elimination).
    fn invert_diagonal_block(&mut self, iv: usize) -> bool {
        let va = &mut self.diag[iv];
        let sol = &mut self.solution[iv];
        
        // Row 1: normalize
        let pivot = va.data[0][0];
        if pivot.abs() < 1e-15 {
            return false;
        }
        let inv_pivot = 1.0 / pivot;
        va.data[0][1] *= inv_pivot;
        va.data[0][2] *= inv_pivot;
        sol[0] *= inv_pivot;
        
        // Eliminate column 1 in rows 2, 3
        for k in 1..3 {
            let factor = va.data[k][0];
            va.data[k][1] -= factor * va.data[0][1];
            va.data[k][2] -= factor * va.data[0][2];
            sol[k] -= factor * sol[0];
        }
        
        // Row 2: normalize
        let pivot = va.data[1][1];
        if pivot.abs() < 1e-15 {
            return false;
        }
        let inv_pivot = 1.0 / pivot;
        va.data[1][2] *= inv_pivot;
        sol[1] *= inv_pivot;
        
        // Eliminate column 2 in row 3
        let factor = va.data[2][1];
        va.data[2][2] -= factor * va.data[1][2];
        sol[2] -= factor * sol[1];
        
        // Row 3: normalize (this is now just the (3,3) element after elimination)
        let pivot = va.data[2][2];
        if pivot.abs() < 1e-15 {
            return false;
        }
        sol[2] /= pivot;
        
        // Back substitute within block to get final solution components
        sol[1] -= va.data[1][2] * sol[2];
        sol[0] -= va.data[0][1] * sol[1] + va.data[0][2] * sol[2];
        
        true
    }
    
    /// Eliminate sub-diagonal block VB[iv+1] using inverted VA[iv].
    fn eliminate_subdiag(&mut self, iv: usize) {
        let ivp = iv + 1;
        if ivp >= self.n_stations {
            return;
        }
        
        // VB[ivp] * VA[iv]^{-1} elimination
        // In XFOIL, this is done row by row
        for k in 0..3 {
            let vb_k0 = self.sub_diag[ivp].data[k][0];
            let vb_k1 = self.sub_diag[ivp].data[k][1];
            let vb_k2 = self.sub_diag[ivp].data[k][2];
            
            // Update RHS
            self.solution[ivp][k] -= vb_k0 * self.solution[iv][0]
                                   + vb_k1 * self.solution[iv][1]
                                   + vb_k2 * self.solution[iv][2];
            
            // Update diagonal block at ivp (already handled by forward sweep)
        }
    }
    
    /// Eliminate wake jump block (VZ) at TE.
    fn eliminate_wake_jump(&mut self, iv: usize) {
        // Wake jump connects upper TE to first wake station
        // Implementation follows XFOIL's VZ elimination in BLSOLV
        if let Some(vz) = &self.wake_jump {
            // VZ is 2x3: rows for theta, mass coupling
            for k in 0..2 {
                let row = k + 1; // Skip ctau row
                self.solution[iv + 1][row] -= vz[k][0] * self.solution[iv][0]
                                            + vz[k][1] * self.solution[iv][1]
                                            + vz[k][2] * self.solution[iv][2];
            }
        }
    }
    
    /// Eliminate lower VM column entries with VACCEL skipping.
    /// 
    /// This is the key VACCEL acceleration: skip elimination steps
    /// when the mass influence coefficients are small.
    fn eliminate_mass_influence_column(&mut self, iv: usize) {
        // Get mass variable (third component) from current station
        let sol_mass = self.solution[iv][2];
        
        for kv in iv + 2..self.n_stations {
            // VM[kv][iv] contains d(res[kv]) / d(mass[iv])
            let vm = self.mass_influence[kv][iv];
            
            // Skip if below threshold (XFOIL's VACCEL acceleration)
            if vm[0].abs() > self.vaccel.vacc1 {
                self.solution[kv][0] -= vm[0] * sol_mass;
            }
            if vm[1].abs() > self.vaccel.vacc2 {
                self.solution[kv][1] -= vm[1] * sol_mass;
            }
            if vm[2].abs() > self.vaccel.vacc3 {
                self.solution[kv][2] -= vm[2] * sol_mass;
            }
        }
    }
    
    /// Back substitute upper VM columns.
    /// 
    /// Propagates the mass variable influence back to upstream stations.
    fn back_substitute_column(&mut self, iv: usize) {
        let sol_mass = self.solution[iv][2];
        
        for kv in 0..iv {
            // VM[kv][iv] = d(res[kv]) / d(mass[iv])
            let vm = self.mass_influence[kv][iv];
            
            // VACCEL skipping in back substitution
            if vm[0].abs() > self.vaccel.vacc1 {
                self.solution[kv][0] -= vm[0] * sol_mass;
            }
            if vm[1].abs() > self.vaccel.vacc2 {
                self.solution[kv][1] -= vm[1] * sol_mass;
            }
            if vm[2].abs() > self.vaccel.vacc3 {
                self.solution[kv][2] -= vm[2] * sol_mass;
            }
        }
    }
    
    /// Get solution at station iv.
    pub fn get_solution(&self, iv: usize) -> Option<[f64; 3]> {
        if iv < self.n_stations {
            Some(self.solution[iv])
        } else {
            None
        }
    }
    
    /// Compute residual norm.
    pub fn rhs_norm(&self) -> f64 {
        let mut sum = 0.0;
        for iv in 0..self.n_stations {
            for k in 0..3 {
                sum += self.rhs[iv][k] * self.rhs[iv][k];
            }
        }
        sum.sqrt()
    }
    
    /// Set up mass influence from DIJ source influence matrix.
    /// 
    /// This is the key V-I coupling: DIJ tells us how source strength (mass defect)
    /// at station J affects surface velocity at station I.
    /// 
    /// The mass influence VM is:
    ///   VM[iv][jv][k] = d(residual[iv][k]) / d(Ue[jv]) * dUe[jv]/dm[jv]
    /// 
    /// where dUe[jv]/dm[jv] comes from DIJ:
    ///   dUe[i] = Σ_j DIJ[panel_i][panel_j] * dm[j]
    /// 
    /// For BL equations, dR/dUe appears in the momentum and shape equations.
    /// 
    /// # Arguments
    /// * `dij` - Source influence matrix (N x N) where N is panel count
    /// * `bl_to_panel` - Maps BL station index to panel node index
    /// * `ue` - Current edge velocity at each BL station
    /// * `theta` - Current momentum thickness at each BL station
    /// * `h` - Current shape factor at each BL station
    /// * `reynolds` - Reynolds number
    pub fn set_mass_influence_from_dij(
        &mut self,
        dij: &[Vec<f64>],
        bl_to_panel: &[usize],
        ue: &[f64],
        theta: &[f64],
        h: &[f64],
        reynolds: f64,
    ) {
        let n_bl = self.n_stations;
        
        for iv in 0..n_bl {
            // Panel index for this BL station
            let panel_i = bl_to_panel.get(iv).copied().unwrap_or(iv);
            let ue_i = ue.get(iv).copied().unwrap_or(1.0).max(1e-10);
            let theta_i = theta.get(iv).copied().unwrap_or(1e-4);
            let h_i = h.get(iv).copied().unwrap_or(2.0);
            
            for jv in 0..n_bl {
                // Panel index for source station
                let panel_j = bl_to_panel.get(jv).copied().unwrap_or(jv);
                
                // Get DIJ influence: dGamma/dSigma
                let dij_val = dij.get(panel_i)
                    .and_then(|row| row.get(panel_j))
                    .copied()
                    .unwrap_or(0.0);
                
                // dUe/dm ≈ DIJ / Ue_j (since sigma = dm/ds, and gamma ≈ Ue)
                let ue_j = ue.get(jv).copied().unwrap_or(1.0).max(1e-10);
                let due_dm = dij_val / ue_j;
                
                // Now compute dR/dUe for each BL equation at station iv
                // 
                // Momentum equation: dθ/ds + (H+2)(θ/Ue)(dUe/ds) = Cf/2
                // Residual form: R_mom = θ_2 - θ_1 - [(H+2)(θ/Ue)(dUe)] - [Cf/2 * ds]
                // dR_mom/dUe has contributions from:
                //   - the (θ/Ue)(dUe/ds) term
                //   - Cf depends on Re_θ = Ue*θ*Re, so dCf/dUe exists
                //
                // Shape equation: similar dependencies
                //
                // For now, use simplified linearization: dR/dUe ≈ scaling factor
                
                // Equation 0: Ctau/Amplification - weak Ue dependence
                let dr0_due = 0.0; // Ctau equation has minimal direct Ue dependence
                
                // Equation 1: Momentum - strong Ue dependence
                // dR_mom/dUe ≈ -(H+2)*θ/Ue² * dUe approximation
                let dr1_due = -(h_i + 2.0) * theta_i / (ue_i * ue_i);
                
                // Equation 2: Mass defect m = Ue*δ* 
                // dm/dUe = δ* = θ*H
                let dr2_due = theta_i * h_i;
                
                // Combined: d(R)/d(m_j) = d(R)/d(Ue_i) * d(Ue_i)/d(m_j)
                self.mass_influence[iv][jv][0] = dr0_due * due_dm;
                self.mass_influence[iv][jv][1] = dr1_due * due_dm;
                self.mass_influence[iv][jv][2] = dr2_due * due_dm;
            }
        }
    }
    
    /// Update edge velocities from mass defect changes using DIJ.
    /// 
    /// After solving the Newton system, the mass defect changes dm need to
    /// feed back to update Ue through the inviscid coupling.
    /// 
    /// # Arguments
    /// * `dij` - Source influence matrix
    /// * `bl_to_panel` - Maps BL station to panel node
    /// 
    /// # Returns
    /// Change in edge velocity at each BL station
    pub fn compute_ue_update(&self, dij: &[Vec<f64>], bl_to_panel: &[usize]) -> Vec<f64> {
        let n_bl = self.n_stations;
        let mut due = vec![0.0; n_bl];
        
        for iv in 0..n_bl {
            let panel_i = bl_to_panel.get(iv).copied().unwrap_or(iv);
            
            for jv in 0..n_bl {
                let panel_j = bl_to_panel.get(jv).copied().unwrap_or(jv);
                
                // Get DIJ influence
                let dij_val = dij.get(panel_i)
                    .and_then(|row| row.get(panel_j))
                    .copied()
                    .unwrap_or(0.0);
                
                // Get mass defect change from solution (third component)
                let dm_j = self.solution[jv][2];
                
                // Accumulate Ue change: dUe[i] += DIJ[i][j] * dm[j]
                due[iv] += dij_val * dm_j;
            }
        }
        
        due
    }
}

// ============================================================================
// XFOIL-Style Newton System for Viscous-Inviscid Coupling
// ============================================================================

/// Configuration for XFOIL-style Newton solver.
#[derive(Debug, Clone)]
pub struct XfoilNewtonConfig {
    /// Maximum Newton iterations for global system
    pub max_iterations: usize,
    /// Convergence tolerance (residual norm)
    pub tolerance: f64,
    /// Reynolds number
    pub reynolds: f64,
    /// Critical N-factor for transition
    pub n_crit: f64,
    /// Initial damping factor
    pub initial_damping: f64,
    /// Minimum damping factor
    pub min_damping: f64,
}

impl Default for XfoilNewtonConfig {
    fn default() -> Self {
        Self {
            max_iterations: 25,
            tolerance: 1e-5,
            reynolds: 1e6,
            n_crit: 9.0,
            initial_damping: 0.5,
            min_damping: 0.1,
        }
    }
}

/// XFOIL-style global Newton system for viscous-inviscid coupling.
/// 
/// This combines:
/// - Block-tridiagonal BL Jacobian (A/B blocks)
/// - Mass influence columns (from DIJ coupling)
/// - Inverse mode handling for separated flows
/// 
/// The system structure matches XFOIL's SETBL/BLSOLV:
/// ```text
/// | A  |     .  |     .  |     |   d       R
/// | B  A     .  |     .  |     |   d   =   R
/// | .  B  A  .  |     .  |     |   d       R
/// | |  |  |  B  A     .  |     |   d       R
/// | .  .  .  .  .  .  .  |     |   d       R
/// ```
/// 
/// Where:
/// - A, B are 3x3 blocks for BL equations
/// - | represents mass influence columns (from DIJ)
/// - d is the solution update [dCtau, dθ, dm]
/// - R is the residual vector
#[derive(Debug)]
pub struct XfoilNewtonSystem {
    /// Configuration
    pub config: XfoilNewtonConfig,
    /// Block-tridiagonal Jacobian for upper surface
    pub upper_jac: BlockTridiagJacobian,
    /// Block-tridiagonal Jacobian for lower surface  
    pub lower_jac: BlockTridiagJacobian,
    /// DIJ source influence matrix (for Ue coupling)
    pub dij: Vec<Vec<f64>>,
    /// Mapping from BL station index to panel node index (upper surface)
    pub upper_bl_to_panel: Vec<usize>,
    /// Mapping from BL station index to panel node index (lower surface)
    pub lower_bl_to_panel: Vec<usize>,
    /// Current edge velocities at each panel node
    pub ue: Vec<f64>,
    /// Inviscid edge velocities (from panel method)
    pub ue_inviscid: Vec<f64>,
    /// Converged flag
    pub converged: bool,
    /// Final residual norm
    pub residual_norm: f64,
    /// Number of iterations taken
    pub iterations: usize,
}

impl XfoilNewtonSystem {
    /// Create a new XFOIL Newton system.
    pub fn new(
        config: XfoilNewtonConfig,
        n_upper: usize,
        n_lower: usize,
        n_panels: usize,
        dij: Vec<Vec<f64>>,
        upper_indices: Vec<usize>,
        lower_indices: Vec<usize>,
        ue_inviscid: Vec<f64>,
    ) -> Self {
        Self {
            config,
            upper_jac: BlockTridiagJacobian::new(n_upper, n_panels),
            lower_jac: BlockTridiagJacobian::new(n_lower, n_panels),
            dij,
            upper_bl_to_panel: upper_indices,
            lower_bl_to_panel: lower_indices,
            ue: ue_inviscid.clone(),
            ue_inviscid,
            converged: false,
            residual_norm: f64::MAX,
            iterations: 0,
        }
    }
    
    /// Update edge velocities from mass defect using DIJ coupling.
    /// 
    /// Ue[i] = Ue_inviscid[i] + Σ_j DIJ[i][j] * mass[j]
    pub fn update_ue_from_mass(&mut self, mass_upper: &[f64], mass_lower: &[f64]) {
        let n_panels = self.ue_inviscid.len();
        self.ue = self.ue_inviscid.clone();
        
        // Build mass at each panel node
        let mut mass_at_panel = vec![0.0; n_panels];
        for (j, &panel_idx) in self.upper_bl_to_panel.iter().enumerate() {
            if panel_idx < n_panels && j < mass_upper.len() {
                mass_at_panel[panel_idx] = mass_upper[j];
            }
        }
        for (j, &panel_idx) in self.lower_bl_to_panel.iter().enumerate() {
            if panel_idx < n_panels && j < mass_lower.len() {
                mass_at_panel[panel_idx] = mass_lower[j];
            }
        }
        
        // Apply DIJ coupling: dUe[i] = Σ_j DIJ[i][j] * mass[j]
        for i in 0..n_panels {
            if i < self.dij.len() {
                for j in 0..n_panels.min(self.dij[i].len()) {
                    self.ue[i] += self.dij[i][j] * mass_at_panel[j];
                }
            }
            self.ue[i] = self.ue[i].abs().max(1e-10);
        }
    }
    
    /// Set up mass influence columns from DIJ matrix.
    /// 
    /// This computes VM[iv][jv][k] = dR_k/dUe * dUe/dm_j for the coupling.
    pub fn setup_mass_influence(
        &mut self,
        upper_states: &[StationState],
        lower_states: &[StationState],
    ) {
        // Upper surface mass influence
        let upper_theta: Vec<f64> = upper_states.iter().map(|s| s.theta).collect();
        let upper_h: Vec<f64> = upper_states.iter().map(|s| s.h).collect();
        let upper_ue: Vec<f64> = self.upper_bl_to_panel.iter()
            .map(|&i| self.ue.get(i).copied().unwrap_or(1.0))
            .collect();
        
        self.upper_jac.set_mass_influence_from_dij(
            &self.dij,
            &self.upper_bl_to_panel,
            &upper_ue,
            &upper_theta,
            &upper_h,
            self.config.reynolds,
        );
        
        // Lower surface mass influence
        let lower_theta: Vec<f64> = lower_states.iter().map(|s| s.theta).collect();
        let lower_h: Vec<f64> = lower_states.iter().map(|s| s.h).collect();
        let lower_ue: Vec<f64> = self.lower_bl_to_panel.iter()
            .map(|&i| self.ue.get(i).copied().unwrap_or(1.0))
            .collect();
        
        self.lower_jac.set_mass_influence_from_dij(
            &self.dij,
            &self.lower_bl_to_panel,
            &lower_ue,
            &lower_theta,
            &lower_h,
            self.config.reynolds,
        );
    }
    
    /// Solve the Newton system for both surfaces.
    /// 
    /// Returns true if both surfaces converged.
    pub fn solve(&mut self) -> bool {
        let upper_ok = self.upper_jac.solve();
        let lower_ok = self.lower_jac.solve();
        
        // Compute combined residual norm
        let upper_norm = self.upper_jac.rhs_norm();
        let lower_norm = self.lower_jac.rhs_norm();
        self.residual_norm = (upper_norm.powi(2) + lower_norm.powi(2)).sqrt();
        
        upper_ok && lower_ok
    }
    
    /// Apply the solution update to BL states.
    /// 
    /// Returns the updated mass arrays for upper and lower surfaces.
    pub fn apply_update(
        &self,
        upper_states: &mut [StationState],
        lower_states: &mut [StationState],
        damping: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        // Apply upper surface updates
        for (j, sol) in self.upper_jac.solution.iter().enumerate() {
            if j < upper_states.len() {
                upper_states[j].ctau_or_ampl = (upper_states[j].ctau_or_ampl + damping * sol[0]).max(0.0);
                upper_states[j].theta = (upper_states[j].theta + damping * sol[1]).max(1e-12);
                upper_states[j].mass = (upper_states[j].mass + damping * sol[2]).max(1e-12);
                upper_states[j].update_h();
            }
        }
        
        // Apply lower surface updates
        for (j, sol) in self.lower_jac.solution.iter().enumerate() {
            if j < lower_states.len() {
                lower_states[j].ctau_or_ampl = (lower_states[j].ctau_or_ampl + damping * sol[0]).max(0.0);
                lower_states[j].theta = (lower_states[j].theta + damping * sol[1]).max(1e-12);
                lower_states[j].mass = (lower_states[j].mass + damping * sol[2]).max(1e-12);
                lower_states[j].update_h();
            }
        }
        
        // Return updated mass arrays
        let mass_upper: Vec<f64> = upper_states.iter().map(|s| s.mass).collect();
        let mass_lower: Vec<f64> = lower_states.iter().map(|s| s.mass).collect();
        (mass_upper, mass_lower)
    }
    
    /// Detect and enable inverse mode for separated stations.
    /// 
    /// Checks each station's shape factor against separation thresholds:
    /// - Laminar: Hk > HK_MAX_LAMINAR (3.8)
    /// - Turbulent: Hk > HK_MAX_TURBULENT (2.5)
    /// 
    /// When separation is detected, enables inverse mode and computes
    /// target Hk using XFOIL's approach.
    pub fn detect_and_set_inverse_mode(
        &self,
        upper_states: &mut [StationState],
        lower_states: &mut [StationState],
    ) {
        use crate::boundary_layer::{HK_MAX_LAMINAR, HK_MAX_TURBULENT, compute_target_hk};
        
        // Process upper surface
        let mut in_inverse_upper = false;
        for j in 1..upper_states.len() {
            let hk_max = if upper_states[j].is_turbulent { HK_MAX_TURBULENT } else { HK_MAX_LAMINAR };
            let hk = upper_states[j].h;
            
            // Check if we should enter inverse mode
            if hk >= hk_max || in_inverse_upper {
                in_inverse_upper = true;
                
                // Compute target Hk (gradually relax toward equilibrium)
                let ds = (upper_states[j].s - upper_states[j - 1].s).abs().max(1e-10);
                let hk_target = compute_target_hk(
                    upper_states[j - 1].h,
                    ds,
                    upper_states[j - 1].theta,
                    upper_states[j].is_turbulent,
                    false, // not wake
                );
                
                upper_states[j].set_inverse_mode(hk_target);
                
                // Check for reattachment
                if hk < hk_max * 0.9 {
                    in_inverse_upper = false;
                    upper_states[j].clear_inverse_mode();
                }
            } else {
                upper_states[j].clear_inverse_mode();
            }
        }
        
        // Process lower surface
        let mut in_inverse_lower = false;
        for j in 1..lower_states.len() {
            let hk_max = if lower_states[j].is_turbulent { HK_MAX_TURBULENT } else { HK_MAX_LAMINAR };
            let hk = lower_states[j].h;
            
            if hk >= hk_max || in_inverse_lower {
                in_inverse_lower = true;
                
                let ds = (lower_states[j].s - lower_states[j - 1].s).abs().max(1e-10);
                let hk_target = compute_target_hk(
                    lower_states[j - 1].h,
                    ds,
                    lower_states[j - 1].theta,
                    lower_states[j].is_turbulent,
                    false,
                );
                
                lower_states[j].set_inverse_mode(hk_target);
                
                if hk < hk_max * 0.9 {
                    in_inverse_lower = false;
                    lower_states[j].clear_inverse_mode();
                }
            } else {
                lower_states[j].clear_inverse_mode();
            }
        }
    }
    
    /// Run the full Newton iteration loop.
    /// 
    /// This performs XFOIL-style coupled Newton iteration:
    /// 1. Detect separation and enable inverse mode
    /// 2. Compute BL residuals and Jacobian
    /// 3. Set up mass influence from DIJ
    /// 4. Solve block-tridiagonal system
    /// 5. Update BL states and Ue
    /// 6. Check convergence
    pub fn iterate(
        &mut self,
        upper_states: &mut [StationState],
        lower_states: &mut [StationState],
    ) -> bool {
        let mut damping = self.config.initial_damping;
        let mut prev_norm = f64::MAX;
        
        for iter in 0..self.config.max_iterations {
            self.iterations = iter + 1;
            
            // Update Ue from current mass defect
            let mass_upper: Vec<f64> = upper_states.iter().map(|s| s.mass).collect();
            let mass_lower: Vec<f64> = lower_states.iter().map(|s| s.mass).collect();
            self.update_ue_from_mass(&mass_upper, &mass_lower);
            
            // Update Ue in station states
            for (j, &panel_idx) in self.upper_bl_to_panel.iter().enumerate() {
                if j < upper_states.len() && panel_idx < self.ue.len() {
                    upper_states[j].ue = self.ue[panel_idx];
                    upper_states[j].update_h();
                }
            }
            for (j, &panel_idx) in self.lower_bl_to_panel.iter().enumerate() {
                if j < lower_states.len() && panel_idx < self.ue.len() {
                    lower_states[j].ue = self.ue[panel_idx];
                    lower_states[j].update_h();
                }
            }
            
            // Detect and enable inverse mode for separated stations
            self.detect_and_set_inverse_mode(upper_states, lower_states);
            
            // Set up mass influence from DIJ
            self.setup_mass_influence(upper_states, lower_states);
            
            // Solve block systems
            if !self.solve() {
                // Reduce damping and retry
                damping = (damping * 0.5).max(self.config.min_damping);
                continue;
            }
            
            // Check convergence
            if self.residual_norm < self.config.tolerance {
                self.converged = true;
                return true;
            }
            
            // Adaptive damping
            if self.residual_norm > prev_norm {
                damping = (damping * 0.7).max(self.config.min_damping);
            } else if self.residual_norm < 0.5 * prev_norm {
                damping = (damping * 1.2).min(1.0);
            }
            prev_norm = self.residual_norm;
            
            // Apply updates
            let (new_mass_upper, new_mass_lower) = self.apply_update(upper_states, lower_states, damping);
            
            // Update Ue for next iteration
            self.update_ue_from_mass(&new_mass_upper, &new_mass_lower);
        }
        
        self.converged = false;
        false
    }
}

/// State at a single BL station for local Newton iteration.
#[derive(Debug, Clone, Copy)]
pub struct StationState {
    /// First variable: Ctau (turbulent) or Ampl (laminar)
    pub ctau_or_ampl: f64,
    /// Momentum thickness theta
    pub theta: f64,
    /// Mass defect m = Ue * delta* (XFOIL uses this instead of delta*)
    pub mass: f64,
    /// Edge velocity at this station
    pub ue: f64,
    /// Arc-length coordinate
    pub s: f64,
    /// Whether turbulent
    pub is_turbulent: bool,
    /// Shape factor H (derived from mass/theta/Ue)
    pub h: f64,
    /// Whether in inverse mode (Hk exceeds separation threshold)
    /// In inverse mode, Hk is prescribed and Ue is solved for.
    pub inverse_mode: bool,
    /// Target Hk when in inverse mode (from compute_target_hk)
    pub hk_target: f64,
}

impl Default for StationState {
    fn default() -> Self {
        Self {
            ctau_or_ampl: 0.0,
            theta: 1e-6,
            mass: 0.0,
            ue: 1.0,
            s: 0.0,
            is_turbulent: false,
            h: 2.6,
            inverse_mode: false,
            hk_target: 0.0,
        }
    }
}

impl StationState {
    /// Create from BL variables.
    pub fn new(ctau_or_ampl: f64, theta: f64, h: f64, ue: f64, s: f64, is_turbulent: bool) -> Self {
        let delta_star = theta * h;
        let mass = ue * delta_star;
        Self {
            ctau_or_ampl,
            theta,
            mass,
            ue,
            s,
            is_turbulent,
            h,
            inverse_mode: false,
            hk_target: 0.0,
        }
    }
    
    /// Create with inverse mode enabled.
    pub fn new_inverse(
        ctau_or_ampl: f64,
        theta: f64,
        h: f64,
        ue: f64,
        s: f64,
        is_turbulent: bool,
        hk_target: f64,
    ) -> Self {
        let delta_star = theta * h;
        let mass = ue * delta_star;
        Self {
            ctau_or_ampl,
            theta,
            mass,
            ue,
            s,
            is_turbulent,
            h,
            inverse_mode: true,
            hk_target,
        }
    }
    
    /// Check if inverse mode should be enabled based on current Hk.
    pub fn should_use_inverse(&self) -> bool {
        use crate::boundary_layer::{HK_MAX_LAMINAR, HK_MAX_TURBULENT};
        let hk_max = if self.is_turbulent { HK_MAX_TURBULENT } else { HK_MAX_LAMINAR };
        self.h >= hk_max
    }
    
    /// Enable inverse mode with a target Hk.
    pub fn set_inverse_mode(&mut self, hk_target: f64) {
        self.inverse_mode = true;
        self.hk_target = hk_target;
    }
    
    /// Disable inverse mode.
    pub fn clear_inverse_mode(&mut self) {
        self.inverse_mode = false;
        self.hk_target = 0.0;
    }
    
    /// Get displacement thickness delta*.
    pub fn delta_star(&self) -> f64 {
        if self.ue.abs() > 1e-10 {
            self.mass / self.ue
        } else {
            self.theta * self.h
        }
    }
    
    /// Update H from mass, theta, Ue.
    pub fn update_h(&mut self) {
        let delta_star = self.delta_star();
        if self.theta > 1e-12 {
            self.h = delta_star / self.theta;
        }
    }
    
    /// Pack into array [ctau, theta, mass].
    pub fn to_array(&self) -> [f64; 3] {
        [self.ctau_or_ampl, self.theta, self.mass]
    }
    
    /// Unpack from array.
    pub fn from_array(&mut self, arr: &[f64; 3]) {
        self.ctau_or_ampl = arr[0];
        self.theta = arr[1];
        self.mass = arr[2];
        self.update_h();
    }
}

/// Newton VII solver.
pub struct NewtonVIISolver {
    config: NewtonConfig,
}

impl NewtonVIISolver {
    /// Create a new Newton solver.
    pub fn new(config: NewtonConfig) -> Self {
        Self { config }
    }
    
    /// Solve the coupled VII system.
    pub fn solve(
        &self,
        initial_state: NewtonState,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
    ) -> NewtonResult {
        let mut state = initial_state;
        let mut residual_history = Vec::new();
        let mut converged = false;
        let mut iteration = 0;
        
        // Initial residual
        let mut residuals = self.compute_residuals(&state, factorized, flow, geometry);
        let mut r_norm = residuals.norm();
        residual_history.push(r_norm);
        
        let r_norm_0 = r_norm.max(1e-10);
        
        for iter in 0..self.config.max_iterations {
            iteration = iter + 1;
            
            // Check convergence
            if r_norm / r_norm_0 < self.config.tolerance {
                converged = true;
                break;
            }
            
            // Build Jacobian using finite differences
            let jacobian = self.build_jacobian_fd(&state, &residuals, factorized, flow, geometry);
            
            // Solve Newton step: J * dx = -R
            let step = match self.solve_newton_step(&jacobian, &residuals) {
                Some(s) => s,
                None => {
                    // Singular Jacobian - try smaller step
                    #[cfg(debug_assertions)]
                    eprintln!("Newton iter {}: Singular Jacobian", iter);
                    break;
                }
            };
            
            // Line search
            let alpha = self.line_search(&state, &step, &residuals, factorized, flow, geometry);
            
            // Update state
            state = state.apply_step(&step, alpha);
            
            // Recompute residuals
            residuals = self.compute_residuals(&state, factorized, flow, geometry);
            r_norm = residuals.norm();
            residual_history.push(r_norm);
            
            #[cfg(debug_assertions)]
            if iter % 5 == 0 {
                eprintln!(
                    "Newton iter {}: ||R||={:.2e}, α={:.3}, ||dx||={:.2e}",
                    iter, r_norm, alpha, step.norm()
                );
            }
        }
        
        NewtonResult {
            state,
            residual_norm: r_norm,
            iterations: iteration,
            converged,
            residual_history,
        }
    }
    
    /// Solve the coupled VII system with DIJ mass influence matrix.
    /// 
    /// This is the full XFOIL-style coupling where:
    /// - DIJ[i][j] = dGamma(i)/dSigma(j) (source influence on surface velocity)
    /// - Mass defect changes dm affect edge velocity: dUe = DIJ * dm
    /// - Newton iteration includes this coupling in the Jacobian
    /// - Uses adaptive damping: start conservative, increase on success
    pub fn solve_with_dij(
        &self,
        initial_state: NewtonState,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
        dij: &[Vec<f64>],
    ) -> NewtonResult {
        let mut state = initial_state;
        let mut residual_history = Vec::new();
        let mut converged = false;
        let mut iteration = 0;
        
        // Adaptive damping: start conservative, increase on success
        let mut damping = self.config.initial_damping;
        let mut prev_r_norm = f64::MAX;
        let mut stall_count = 0;
        
        // Initial residual (using DIJ for mass influence)
        let mut residuals = self.compute_residuals_with_dij(&state, factorized, flow, geometry, dij);
        let mut r_norm = residuals.norm();
        residual_history.push(r_norm);
        
        let r_norm_0 = r_norm.max(1e-10);
        
        for iter in 0..self.config.max_iterations {
            iteration = iter + 1;
            
            // Check convergence
            if r_norm / r_norm_0 < self.config.tolerance {
                converged = true;
                break;
            }
            
            // Build Jacobian with DIJ mass influence
            let jacobian = self.build_jacobian_with_dij(&state, &residuals, factorized, flow, geometry, dij);
            
            // Solve Newton step: J * dx = -R
            let step = match self.solve_newton_step(&jacobian, &residuals) {
                Some(s) => s,
                None => {
                    #[cfg(debug_assertions)]
                    eprintln!("Newton iter {}: Singular Jacobian (DIJ)", iter);
                    // Reduce damping and try again
                    damping = (damping * 0.5).max(self.config.min_damping);
                    stall_count += 1;
                    if stall_count > 5 {
                        break;
                    }
                    continue;
                }
            };
            
            // Line search with DIJ residuals, scaled by damping
            let alpha_ls = self.line_search_with_dij(&state, &step, &residuals, factorized, flow, geometry, dij);
            let alpha = alpha_ls * damping;
            
            // Update state
            let state_new = state.apply_step(&step, alpha);
            
            // Apply inverse mode gamma update (XFOIL-style: solve for Ue from prescribed Hk)
            // Note: DIJ coupling is now fully in the Jacobian, so no separate post-step correction
            let mut state_corrected = state_new;
            self.apply_inverse_mode_gamma_update(&mut state_corrected, geometry);
            
            // Recompute residuals
            let residuals_new = self.compute_residuals_with_dij(&state_corrected, factorized, flow, geometry, dij);
            let r_norm_new = residuals_new.norm();
            
            // Adaptive damping update
            if r_norm_new < r_norm * 0.99 {
                // Making progress: increase damping
                damping = (damping * 1.2).min(1.0);
                stall_count = 0;
            } else if r_norm_new > prev_r_norm * 1.1 {
                // Diverging: decrease damping
                damping = (damping * 0.5).max(self.config.min_damping);
                stall_count += 1;
            }
            
            // Accept step
            state = state_corrected;
            residuals = residuals_new;
            prev_r_norm = r_norm;
            r_norm = r_norm_new;
            residual_history.push(r_norm);
            
            #[cfg(debug_assertions)]
            if iter % 5 == 0 {
                eprintln!(
                    "Newton (DIJ) iter {}: ||R||={:.2e}, α={:.3}, damp={:.2}, ||dx||={:.2e}",
                    iter, r_norm, alpha, damping, step.norm()
                );
            }
            
            // Check for stalled convergence
            if stall_count > 10 {
                #[cfg(debug_assertions)]
                eprintln!("Newton (DIJ): Stalled after {} iterations", iter);
                break;
            }
        }
        
        NewtonResult {
            state,
            residual_norm: r_norm,
            iterations: iteration,
            converged,
            residual_history,
        }
    }
    
    /// Compute residuals including DIJ mass influence on edge velocity.
    fn compute_residuals_with_dij(
        &self,
        state: &NewtonState,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
        dij: &[Vec<f64>],
    ) -> Residuals {
        use crate::viscous::compute_transpiration;
        
        let n_panels = geometry.n_panels;
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        let n_bl = n_upper + n_lower;
        
        let mut residuals = Residuals::new(n_panels, n_bl);
        
        // 1. Compute base edge velocity from gamma
        let ue_base: Vec<f64> = state.gamma.iter().map(|g| g.abs()).collect();
        
        // 2. Compute mass defect m = Ue * delta_star at each BL station
        let delta_star = state.delta_star();
        let mass_defect: Vec<f64> = (0..n_bl).map(|j| {
            let panel_idx = if j < n_upper {
                geometry.upper_indices[j]
            } else {
                geometry.lower_indices[j - n_upper]
            };
            ue_base[panel_idx] * delta_star[j]
        }).collect();
        
        // 3. Compute edge velocity correction from DIJ coupling
        // Map mass defect back to panel nodes
        let mut dm_panels = vec![0.0; n_panels];
        for j in 0..n_upper {
            let panel_idx = geometry.upper_indices[j];
            dm_panels[panel_idx] = mass_defect[j];
        }
        for j in 0..n_lower {
            let panel_idx = geometry.lower_indices[j];
            dm_panels[panel_idx] = mass_defect[n_upper + j];
        }
        
        // dUe from DIJ: dUe[i] = sum_j DIJ[i][j] * m[j]
        let due_from_dij: Vec<f64> = crate::inviscid::compute_ue_change_from_mass(dij, &dm_panels);
        
        // 4. Effective edge velocity includes DIJ correction (XFOIL-style: Ue = Uinv + DIJ*m)
        let ue_effective: Vec<f64> = ue_base.iter()
            .zip(due_from_dij.iter())
            .map(|(&ue, &due)| (ue + due).max(1e-10))
            .collect();
        
        // 5. Map delta_star to panel nodes for transpiration
        let delta_star_nodes = self.map_bl_to_nodes(&delta_star, geometry);
        
        // 6. Compute transpiration Vn = d(Ue*δ*)/ds using DIJ-corrected Ue
        let vn = compute_transpiration(&ue_effective, &delta_star_nodes, &geometry.s_coords);
        
        // 7. Inviscid residuals with transpiration
        let inviscid_with_vn = factorized.solve_with_transpiration(flow, &vn, &geometry.s_coords);
        for i in 0..n_panels {
            residuals.r_inviscid[i] = state.gamma[i] - inviscid_with_vn.gamma[i];
        }
        
        // 8. BL residuals (momentum and shape equations) using DIJ-corrected Ue
        self.compute_bl_residuals(
            state,
            &ue_effective,
            geometry,
            &mut residuals.r_momentum,
            &mut residuals.r_shape,
        );
        
        residuals
    }
    
    /// Build Jacobian with DIJ mass influence included.
    fn build_jacobian_with_dij(
        &self,
        state: &NewtonState,
        residuals: &Residuals,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
        dij: &[Vec<f64>],
    ) -> DMatrix<f64> {
        // Start with base Jacobian from finite differences
        let mut jacobian = self.build_jacobian_fd(state, residuals, factorized, flow, geometry);
        
        let n_panels = geometry.n_panels;
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        let n_bl = n_upper + n_lower;
        
        // Add mass influence coupling: dR_inviscid/dm through DIJ
        // For each inviscid equation i, add sensitivity to mass defect changes
        // 
        // R_inviscid[i] = gamma[i] - gamma_target[i](transpiration)
        // dR_inviscid[i]/dm[j] = -d(gamma_target[i])/dm[j] = -DIJ[i][panel_j]
        
        for i in 0..n_panels {
            for j in 0..n_bl {
                // Map BL station j to panel index
                let panel_j = if j < n_upper {
                    geometry.upper_indices[j]
                } else {
                    geometry.lower_indices[j - n_upper]
                };
                
                // Get DIJ influence
                let dij_val = dij.get(i)
                    .and_then(|row| row.get(panel_j))
                    .copied()
                    .unwrap_or(0.0);
                
                // Add to Jacobian: column for shape factor H relates to delta_star
                // h[j] is at column n_panels + n_bl + j
                let col_h = n_panels + n_bl + j;
                if col_h < jacobian.ncols() {
                    // dR_i/dH_j = -DIJ[i][panel_j] * dδ*/dH * Ue
                    // dδ*/dH = θ (since δ* = θ*H)
                    // Use moderate coupling (0.5) for stability with line search
                    let theta_j = state.theta[j].max(1e-10);
                    let ue_j = state.gamma[panel_j].abs().max(1e-10);
                    jacobian[(i, col_h)] += -dij_val * theta_j * ue_j * 0.5;
                }
            }
        }
        
        jacobian
    }
    
    /// Line search with DIJ residuals.
    fn line_search_with_dij(
        &self,
        state: &NewtonState,
        step: &NewtonStep,
        residuals: &Residuals,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
        dij: &[Vec<f64>],
    ) -> f64 {
        let mut alpha = 1.0;
        let r_norm_0 = residuals.norm();
        
        for _ in 0..self.config.line_search_max_iter {
            let state_new = state.apply_step(step, alpha);
            let r_new = self.compute_residuals_with_dij(&state_new, factorized, flow, geometry, dij);
            let r_norm = r_new.norm();
            
            // Armijo condition
            if r_norm < r_norm_0 * (1.0 - self.config.armijo_c * alpha) {
                return alpha;
            }
            
            alpha *= 0.5;
        }
        
        alpha
    }
    
    /// Apply DIJ-based edge velocity correction after Newton step.
    /// 
    /// The mass defect changes dm propagate to edge velocity through:
    /// dUe = DIJ * dm
    fn apply_dij_velocity_correction(
        &self,
        state: &mut NewtonState,
        step: &NewtonStep,
        alpha: f64,
        geometry: &NewtonGeometry,
        dij: &[Vec<f64>],
    ) {
        let n_panels = geometry.n_panels;
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        let n_bl = n_upper + n_lower;
        
        // Compute mass defect changes from the step
        // dm = d(Ue*δ*) ≈ Ue * d(θ*H) + θ*H * dUe
        // For simplicity, use dm ≈ Ue * θ * dH (primary contribution)
        let mut dm_panels = vec![0.0; n_panels];
        
        for j in 0..n_bl {
            let panel_idx = if j < n_upper {
                geometry.upper_indices[j]
            } else {
                geometry.lower_indices[j - n_upper]
            };
            
            // d_mass is now the direct change in mass defect
            dm_panels[panel_idx] = step.d_mass[j] * alpha;
        }
        
        // Compute Ue correction from DIJ
        let due = crate::inviscid::compute_ue_change_from_mass(dij, &dm_panels);
        
        // Apply full correction to gamma (gamma ≈ Ue for attached flow)
        // XFOIL-style: full DIJ coupling, stability from Newton damping
        let n_apply = n_panels.min(due.len()).min(state.gamma.len());
        for i in 0..n_apply {
            let sign = state.gamma[i].signum();
            state.gamma[i] += sign * due[i];
        }
    }
    
    /// Update gamma for inverse mode stations (XFOIL-style MRCHUE/QVFUE coupling).
    /// 
    /// When in inverse mode, the shape factor Hk is prescribed. We compute the
    /// edge velocity Ue that gives this Hk and update gamma accordingly.
    /// 
    /// From H = mass / (Ue * theta) and H = Hk_target:
    ///   Ue = mass / (theta * Hk_target)
    fn apply_inverse_mode_gamma_update(
        &self,
        state: &mut NewtonState,
        geometry: &NewtonGeometry,
    ) {
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        let n_bl = n_upper + n_lower;
        
        for j in 0..n_bl {
            if !state.inverse_mode[j] {
                continue;
            }
            
            // Get panel index for this BL station
            let panel_idx = if j < n_upper {
                geometry.upper_indices[j]
            } else {
                geometry.lower_indices[j - n_upper]
            };
            
            if panel_idx >= state.gamma.len() {
                continue;
            }
            
            // Compute Ue that gives the target Hk
            // H = mass / (Ue * theta), so Ue = mass / (theta * Hk_target)
            let theta = state.theta[j].max(1e-10);
            let mass = state.mass[j].max(1e-12);
            let hk_target = state.hk_target[j].max(1.05);
            
            let ue_from_inverse = mass / (theta * hk_target);
            let ue_bounded = ue_from_inverse.clamp(0.01, 3.0); // Physical bounds
            
            // Update gamma to reflect the inverse mode Ue
            // Preserve sign of gamma
            let sign = state.gamma[panel_idx].signum();
            let current_ue = state.gamma[panel_idx].abs();
            
            // Blend toward inverse mode Ue (relaxed update for stability)
            let relaxation = 0.5;
            let new_ue = current_ue * (1.0 - relaxation) + ue_bounded * relaxation;
            
            state.gamma[panel_idx] = sign * new_ue;
        }
    }
    
    /// Compute all residuals at current state.
    pub fn compute_residuals(
        &self,
        state: &NewtonState,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
    ) -> Residuals {
        let n_panels = geometry.n_panels;
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        let n_bl = n_upper + n_lower;
        
        let mut residuals = Residuals::new(n_panels, n_bl);
        
        // 1. Compute edge velocity from gamma
        let ue: Vec<f64> = state.gamma.iter().map(|g| g.abs()).collect();
        
        // 2. Compute displacement thickness
        let delta_star = state.delta_star();
        
        // 3. Map delta_star to panel nodes for transpiration
        let delta_star_nodes = self.map_bl_to_nodes(&delta_star, geometry);
        let ue_nodes = &ue;
        
        // 4. Compute transpiration
        let vn = compute_transpiration(ue_nodes, &delta_star_nodes, &geometry.s_coords);
        
        // 5. Inviscid residuals: compare gamma to what factorized solution gives with transpiration
        let inviscid_with_vn = factorized.solve_with_transpiration(flow, &vn, &geometry.s_coords);
        for i in 0..n_panels {
            residuals.r_inviscid[i] = state.gamma[i] - inviscid_with_vn.gamma[i];
        }
        
        // 6. BL residuals (momentum and shape equations)
        self.compute_bl_residuals(
            state,
            &ue,
            geometry,
            &mut residuals.r_momentum,
            &mut residuals.r_shape,
        );
        
        residuals
    }
    
    /// Compute BL residuals for momentum and shape equations.
    fn compute_bl_residuals(
        &self,
        state: &NewtonState,
        ue: &[f64],
        geometry: &NewtonGeometry,
        r_momentum: &mut [f64],
        r_shape: &mut [f64],
    ) {
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        
        // Process upper surface
        self.compute_surface_bl_residuals(
            state,
            ue,
            geometry,
            &geometry.upper_indices,
            0, // offset in BL arrays
            r_momentum,
            r_shape,
        );
        
        // Process lower surface  
        self.compute_surface_bl_residuals(
            state,
            ue,
            geometry,
            &geometry.lower_indices,
            n_upper, // offset in BL arrays
            r_momentum,
            r_shape,
        );
    }
    
    /// Compute BL residuals for one surface.
    fn compute_surface_bl_residuals(
        &self,
        state: &NewtonState,
        ue: &[f64],
        geometry: &NewtonGeometry,
        surface_indices: &[usize],
        bl_offset: usize,
        r_momentum: &mut [f64],
        r_shape: &mut [f64],
    ) {
        let n_stations = surface_indices.len();
        if n_stations < 2 {
            return;
        }
        
        // First station: initial condition residual
        // theta should be small (stagnation point)
        let theta_init = 1e-4 / self.config.reynolds.sqrt();
        let panel_0 = surface_indices[0];
        let ue_0 = ue[panel_0].max(1e-10);
        
        r_momentum[bl_offset] = state.theta[bl_offset] - theta_init;
        
        // For mass, initial condition: mass = Ue * θ * H_blasius
        let h_blasius = 2.59;
        let mass_init = ue_0 * theta_init * h_blasius;
        r_shape[bl_offset] = state.mass[bl_offset] - mass_init;
        
        // March along surface
        for j in 1..n_stations {
            let bl_idx = bl_offset + j;
            let bl_idx_prev = bl_offset + j - 1;
            
            let panel_idx = surface_indices[j];
            let panel_idx_prev = surface_indices[j - 1];
            
            // Bounds check
            if panel_idx >= ue.len() || panel_idx_prev >= ue.len() ||
               panel_idx >= geometry.s_coords.len() || panel_idx_prev >= geometry.s_coords.len() {
                continue;
            }
            
            // Get local values
            let s = geometry.s_coords[panel_idx];
            let s_prev = geometry.s_coords[panel_idx_prev];
            let ds = (s - s_prev).abs().max(1e-10);
            
            let ue_local = ue[panel_idx].max(1e-10);
            let ue_prev = ue[panel_idx_prev].max(1e-10);
            let ue_avg = 0.5 * (ue_local + ue_prev);
            
            let theta = state.theta[bl_idx];
            let theta_prev = state.theta[bl_idx_prev];
            
            // Compute H from mass: H = mass / (Ue * θ)
            let h = (state.mass[bl_idx] / (ue_local * theta.max(1e-10))).clamp(1.05, 10.0);
            let h_prev = (state.mass[bl_idx_prev] / (ue_prev * theta_prev.max(1e-10))).clamp(1.05, 10.0);
            
            // Velocity gradient
            let due_ds = (ue_local - ue_prev) / ds;
            
            // Skin friction
            let re_theta = ue_avg * theta_prev * self.config.reynolds;
            let cf = if state.is_turbulent[bl_idx] {
                ludwieg_tillmann_cf(re_theta, h_prev)
            } else {
                // Laminar Cf
                0.664 / re_theta.sqrt().max(1.0)
            };
            
            // Momentum equation residual:
            // dθ/ds = Cf/2 - (H+2)(θ/Ue)(dUe/ds)
            let dtheta_ds_computed = cf / 2.0 - (h_prev + 2.0) * theta_prev * due_ds / ue_avg;
            let dtheta_ds_actual = (theta - theta_prev) / ds;
            r_momentum[bl_idx] = dtheta_ds_actual - dtheta_ds_computed;
            
            // Shape equation residual
            // In inverse mode: Hk - Hk_target = 0
            // In direct mode: d(H1*θ)/ds = Ce
            if state.inverse_mode[bl_idx] {
                let hk_target = state.hk_target[bl_idx];
                r_shape[bl_idx] = h - hk_target;
            } else {
                // Direct mode: Head's entrainment
                let h1 = head_h1(h);
                let h1_prev = head_h1(h_prev);
                let ce = head_entrainment(h1_prev);
                
                let d_h1_theta_ds_actual = (h1 * theta - h1_prev * theta_prev) / ds;
                r_shape[bl_idx] = d_h1_theta_ds_actual - ce;
            }
        }
    }
    
    /// Map BL values to panel node indices.
    fn map_bl_to_nodes(&self, bl_values: &[f64], geometry: &NewtonGeometry) -> Vec<f64> {
        let mut node_values = vec![0.0; geometry.n_panels];
        
        let n_upper = geometry.upper_indices.len();
        
        // Upper surface
        for (j, &panel_idx) in geometry.upper_indices.iter().enumerate() {
            if panel_idx < node_values.len() && j < bl_values.len() {
                node_values[panel_idx] = bl_values[j];
            }
        }
        
        // Lower surface
        for (j, &panel_idx) in geometry.lower_indices.iter().enumerate() {
            let bl_idx = n_upper + j;
            if panel_idx < node_values.len() && bl_idx < bl_values.len() {
                node_values[panel_idx] = bl_values[bl_idx];
            }
        }
        
        node_values
    }
    
    /// Build Jacobian matrix using finite differences.
    pub fn build_jacobian_fd(
        &self,
        state: &NewtonState,
        residuals: &Residuals,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
    ) -> DMatrix<f64> {
        let n_vars_state = state.total_size();
        let n_vars_res = residuals.total_size();
        let eps = self.config.fd_epsilon;
        
        // State and residual sizes must match for square Jacobian
        debug_assert_eq!(
            n_vars_state, n_vars_res,
            "State size {} != Residual size {}. State: gamma={}, theta={}, mass={}. \
             Residuals: inv={}, mom={}, shape={}",
            n_vars_state, n_vars_res,
            state.gamma.len(), state.theta.len(), state.mass.len(),
            residuals.r_inviscid.len(), residuals.r_momentum.len(), residuals.r_shape.len()
        );
        
        // Use the minimum to avoid out-of-bounds in non-debug builds
        let n_vars = n_vars_state.min(n_vars_res);
        
        let mut jacobian = DMatrix::zeros(n_vars, n_vars);
        
        for i in 0..n_vars {
            // Perturb state[i]
            let mut state_plus = state.clone();
            state_plus.perturb(i, eps);
            
            let r_plus = self.compute_residuals(&state_plus, factorized, flow, geometry);
            
            // J[:,i] = (R_plus - R) / eps
            for j in 0..n_vars.min(r_plus.total_size()) {
                let deriv = (r_plus.get(j) - residuals.get(j)) / eps;
                if deriv.abs() > 1e-15 {
                    jacobian[(j, i)] = deriv;
                }
            }
        }
        
        jacobian
    }
    
    /// Solve Newton step: J * dx = -R
    /// 
    /// Uses Tikhonov regularization to improve conditioning when the
    /// Jacobian is ill-conditioned or near-singular.
    fn solve_newton_step(
        &self,
        jacobian: &DMatrix<f64>,
        residuals: &Residuals,
    ) -> Option<NewtonStep> {
        let r_vec = DVector::from_vec(residuals.to_vec());
        let neg_r = -&r_vec;
        
        // Tikhonov regularization: J_reg = J + λI
        // Start with small λ and increase if needed
        let n = jacobian.nrows();
        let mut jacobian_reg = jacobian.clone();
        
        // Compute adaptive regularization parameter based on diagonal magnitude
        let diag_mean = (0..n)
            .map(|i| jacobian[(i, i)].abs())
            .sum::<f64>() / n as f64;
        let mut lambda = 1e-8 * diag_mean.max(1e-10);
        
        // Try increasing regularization until solve succeeds
        for attempt in 0..5 {
            // Add regularization to diagonal
            if attempt > 0 {
                for i in 0..n {
                    jacobian_reg[(i, i)] = jacobian[(i, i)] + lambda;
                }
            }
            
            // Use LU decomposition
            let lu = jacobian_reg.clone().lu();
            if let Some(dx) = lu.solve(&neg_r) {
                let n_gamma = residuals.r_inviscid.len();
                let n_bl = residuals.r_momentum.len();
                return Some(NewtonStep::from_vec(dx.as_slice(), n_gamma, n_bl));
            }
            
            // Increase regularization for next attempt
            lambda *= 10.0;
        }
        
        None // Failed even with regularization
    }
    
    /// Backtracking line search to ensure residual reduction.
    fn line_search(
        &self,
        state: &NewtonState,
        step: &NewtonStep,
        residuals: &Residuals,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
    ) -> f64 {
        let mut alpha = 1.0;
        let r_norm_0 = residuals.norm();
        
        for _ in 0..self.config.line_search_max_iter {
            let state_new = state.apply_step(step, alpha);
            let r_new = self.compute_residuals(&state_new, factorized, flow, geometry);
            let r_norm = r_new.norm();
            
            // Armijo condition
            if r_norm < r_norm_0 * (1.0 - self.config.armijo_c * alpha) {
                return alpha;
            }
            
            alpha *= 0.5;
        }
        
        alpha // Return best found
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_newton_state_indexing() {
        let mut state = NewtonState::new(10, 5);
        
        // Set some values
        state.gamma[0] = 1.0;
        state.theta[0] = 2.0;
        state.mass[0] = 3.0;
        
        // Test get
        assert!((state.get(0) - 1.0).abs() < 1e-10); // gamma[0]
        assert!((state.get(10) - 2.0).abs() < 1e-10); // theta[0]
        assert!((state.get(15) - 3.0).abs() < 1e-10); // mass[0]
        
        // Test set
        state.set(0, 10.0);
        assert!((state.gamma[0] - 10.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_newton_state_total_size() {
        let state = NewtonState::new(100, 50);
        assert_eq!(state.total_size(), 100 + 50 + 50); // gamma + theta + mass
    }
    
    #[test]
    fn test_residuals_norm() {
        let mut residuals = Residuals::new(10, 5);
        residuals.r_inviscid[0] = 3.0;
        residuals.r_momentum[0] = 4.0;
        
        let norm = residuals.norm();
        assert!((norm - 5.0).abs() < 1e-10); // sqrt(9 + 16) = 5
    }
    
    #[test]
    fn test_newton_config_default() {
        let config = NewtonConfig::default();
        assert!((config.reynolds - 1e6).abs() < 1.0);
        assert_eq!(config.max_iterations, 50);
    }
    
    #[test]
    fn test_bl_block_solve() {
        // Test 3x3 block solve
        let block = BLBlock {
            data: [
                [2.0, 1.0, 0.0],
                [1.0, 3.0, 1.0],
                [0.0, 1.0, 2.0],
            ],
        };
        
        let rhs = [1.0, 2.0, 3.0];
        let solution = block.solve(&rhs).expect("Block should be non-singular");
        
        // Verify: A * x should equal b
        let ax = block.mul_vec(&solution);
        for i in 0..3 {
            assert!((ax[i] - rhs[i]).abs() < 1e-10, "ax[{}] = {}, expected {}", i, ax[i], rhs[i]);
        }
    }
    
    #[test]
    fn test_bl_block_singular() {
        // Test singular block detection
        let singular = BLBlock {
            data: [
                [1.0, 2.0, 3.0],
                [2.0, 4.0, 6.0], // Row 2 = 2 * Row 1
                [1.0, 1.0, 1.0],
            ],
        };
        
        let rhs = [1.0, 2.0, 3.0];
        assert!(singular.solve(&rhs).is_none(), "Should detect singular block");
    }
    
    #[test]
    fn test_block_tridiag_jacobian_simple() {
        // Test block-tridiagonal solve for a simple 3-station case
        let mut jacobian = BlockTridiagJacobian::new(3, 3);
        
        // Set up identity-ish diagonal blocks
        for iv in 0..3 {
            let mut diag = BLBlock::identity();
            // Make it slightly more interesting
            diag.data[0][0] = 2.0;
            diag.data[1][1] = 2.0;
            diag.data[2][2] = 2.0;
            jacobian.set_diag(iv, diag);
        }
        
        // Set up sub-diagonal blocks (coupling to previous)
        for iv in 1..3 {
            let mut sub = BLBlock::zero();
            sub.data[0][0] = -1.0;
            sub.data[1][1] = -1.0;
            sub.data[2][2] = -1.0;
            jacobian.set_subdiag(iv, sub);
        }
        
        // Set RHS
        jacobian.set_rhs(0, [1.0, 1.0, 1.0]);
        jacobian.set_rhs(1, [0.0, 0.0, 0.0]);
        jacobian.set_rhs(2, [0.0, 0.0, 0.0]);
        
        // Solve
        let success = jacobian.solve();
        assert!(success, "Block solve should succeed");
        
        // Solution should be finite
        for iv in 0..3 {
            let sol = jacobian.get_solution(iv).unwrap();
            for k in 0..3 {
                assert!(sol[k].is_finite(), "Solution[{}][{}] = {}", iv, k, sol[k]);
            }
        }
    }
    
    #[test]
    fn test_vaccel_params() {
        let params = VaccelParams::default();
        assert!((params.vaccel - 0.01).abs() < 1e-10);
        
        let scaled = VaccelParams::from_arc_length(0.01, 2.0);
        assert!((scaled.vacc2 - 0.01).abs() < 1e-10);
    }
    
    #[test]
    fn test_newton_scaling() {
        let scaling = NewtonScaling::from_flow(1e6, 1.0, 1.5);
        
        // Theta scale should be ~ chord/sqrt(Re) = 1/1000 = 0.001
        assert!(scaling.s_theta > 0.0 && scaling.s_theta < 0.01);
        
        // Mass scale should include Ue
        assert!(scaling.s_mass > scaling.s_theta);
        
        // Test scale/unscale roundtrip
        let mut v = [0.01, 0.001, 0.002];
        let v_orig = v;
        scaling.scale(&mut v);
        scaling.unscale(&mut v);
        for i in 0..3 {
            assert!((v[i] - v_orig[i]).abs() < 1e-12);
        }
    }
    
    #[test]
    fn test_station_state() {
        let state = StationState::new(0.03, 0.001, 2.5, 1.0, 0.5, true);
        
        assert!((state.theta - 0.001).abs() < 1e-10);
        assert!((state.h - 2.5).abs() < 1e-10);
        assert!((state.ctau_or_ampl - 0.03).abs() < 1e-10);
        
        // Delta* should be theta * H
        let expected_delta_star = 0.001 * 2.5;
        assert!((state.delta_star() - expected_delta_star).abs() < 1e-10);
        
        // Mass should be Ue * delta*
        let expected_mass = 1.0 * expected_delta_star;
        assert!((state.mass - expected_mass).abs() < 1e-10);
    }
    
    #[test]
    fn test_mass_influence_from_dij() {
        // Create a simple 4-station system
        let n_stations = 4;
        let mut jacobian = BlockTridiagJacobian::new(n_stations, n_stations);
        
        // Simple DIJ matrix (4x4)
        let dij = vec![
            vec![0.0, 0.1, 0.05, 0.02],
            vec![0.1, 0.0, 0.1, 0.05],
            vec![0.05, 0.1, 0.0, 0.1],
            vec![0.02, 0.05, 0.1, 0.0],
        ];
        
        // BL-to-panel mapping (identity in this case)
        let bl_to_panel: Vec<usize> = (0..n_stations).collect();
        
        // Mock BL data
        let ue = vec![1.0, 1.1, 1.05, 0.95];
        let theta = vec![0.001, 0.0012, 0.0015, 0.002];
        let h = vec![2.5, 2.3, 2.2, 2.0];
        let reynolds = 1e6;
        
        // Set mass influence from DIJ
        jacobian.set_mass_influence_from_dij(&dij, &bl_to_panel, &ue, &theta, &h, reynolds);
        
        // Check that mass influence was set
        // The off-diagonal entries should be non-zero where DIJ is non-zero
        for iv in 0..n_stations {
            for jv in 0..n_stations {
                if iv != jv && dij[iv][jv].abs() > 1e-10 {
                    // At least equation 1 (momentum) should have influence
                    let vm = jacobian.mass_influence[iv][jv];
                    // Equation 1 has strongest Ue dependence
                    assert!(vm[1].abs() > 0.0 || vm[2].abs() > 0.0,
                        "Mass influence VM[{}][{}] should be non-zero", iv, jv);
                }
            }
        }
    }
    
    #[test]
    fn test_compute_ue_update() {
        let n_stations = 4;
        let mut jacobian = BlockTridiagJacobian::new(n_stations, n_stations);
        
        // Set solution (third component is mass defect change)
        jacobian.solution[0] = [0.0, 0.0, 0.001];
        jacobian.solution[1] = [0.0, 0.0, 0.002];
        jacobian.solution[2] = [0.0, 0.0, 0.001];
        jacobian.solution[3] = [0.0, 0.0, 0.0005];
        
        // Simple DIJ matrix
        let dij = vec![
            vec![0.0, 0.1, 0.05, 0.02],
            vec![0.1, 0.0, 0.1, 0.05],
            vec![0.05, 0.1, 0.0, 0.1],
            vec![0.02, 0.05, 0.1, 0.0],
        ];
        
        let bl_to_panel: Vec<usize> = (0..n_stations).collect();
        
        let due = jacobian.compute_ue_update(&dij, &bl_to_panel);
        
        assert_eq!(due.len(), n_stations);
        
        // dUe[0] = 0.1*0.002 + 0.05*0.001 + 0.02*0.0005 = 0.00026
        let expected_due_0 = 0.1 * 0.002 + 0.05 * 0.001 + 0.02 * 0.0005;
        assert!((due[0] - expected_due_0).abs() < 1e-10,
            "dUe[0] = {} expected {}", due[0], expected_due_0);
    }
}
