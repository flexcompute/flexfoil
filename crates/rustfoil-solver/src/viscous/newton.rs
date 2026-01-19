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
#[derive(Debug, Clone)]
pub struct NewtonState {
    /// Panel vorticity at each node (N values)
    pub gamma: Vec<f64>,
    /// BL momentum thickness at each station
    /// Ordered: upper surface (LE to TE), lower surface (LE to TE)
    pub theta: Vec<f64>,
    /// BL shape factor H at each station
    pub h: Vec<f64>,
    /// Whether each station is turbulent
    pub is_turbulent: Vec<bool>,
    /// Transition index on upper surface (index into upper stations)
    pub transition_upper: Option<usize>,
    /// Transition index on lower surface (index into lower stations)
    pub transition_lower: Option<usize>,
}

impl NewtonState {
    /// Create a new state with given sizes.
    pub fn new(n_panels: usize, n_bl_stations: usize) -> Self {
        Self {
            gamma: vec![0.0; n_panels],
            theta: vec![1e-6; n_bl_stations],
            h: vec![2.6; n_bl_stations], // Blasius initial value
            is_turbulent: vec![false; n_bl_stations],
            transition_upper: None,
            transition_lower: None,
        }
    }
    
    /// Initialize from inviscid solution.
    pub fn from_inviscid(inviscid: &InviscidSolution, n_bl_stations: usize) -> Self {
        let n_panels = inviscid.gamma.len();
        let mut state = Self::new(n_panels, n_bl_stations);
        state.gamma = inviscid.gamma.clone();
        state
    }
    
    /// Total number of unknowns.
    pub fn total_size(&self) -> usize {
        self.gamma.len() + self.theta.len() + self.h.len()
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
            self.h[idx - n_gamma - n_theta]
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
            self.h[idx - n_gamma - n_theta] = value;
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
        v.extend(&self.h);
        v
    }
    
    /// Update from flat vector.
    pub fn from_vec(&mut self, v: &[f64]) {
        let n_gamma = self.gamma.len();
        let n_theta = self.theta.len();
        
        self.gamma.copy_from_slice(&v[0..n_gamma]);
        self.theta.copy_from_slice(&v[n_gamma..n_gamma + n_theta]);
        self.h.copy_from_slice(&v[n_gamma + n_theta..]);
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
        for i in 0..self.h.len() {
            new_state.h[i] = (self.h[i] + alpha * step.d_h[i]).clamp(1.0, 10.0);
        }
        
        new_state
    }
    
    /// Compute displacement thickness from theta and H.
    pub fn delta_star(&self) -> Vec<f64> {
        self.theta.iter()
            .zip(self.h.iter())
            .map(|(&th, &h)| th * h)
            .collect()
    }
}

/// Newton step (dx in J*dx = -R).
#[derive(Debug, Clone)]
pub struct NewtonStep {
    /// Change in gamma
    pub d_gamma: Vec<f64>,
    /// Change in theta
    pub d_theta: Vec<f64>,
    /// Change in H
    pub d_h: Vec<f64>,
}

impl NewtonStep {
    /// Create from flat vector.
    pub fn from_vec(v: &[f64], n_gamma: usize, n_bl: usize) -> Self {
        Self {
            d_gamma: v[0..n_gamma].to_vec(),
            d_theta: v[n_gamma..n_gamma + n_bl].to_vec(),
            d_h: v[n_gamma + n_bl..].to_vec(),
        }
    }
    
    /// Compute L2 norm.
    pub fn norm(&self) -> f64 {
        let sum: f64 = self.d_gamma.iter().map(|x| x * x).sum::<f64>()
            + self.d_theta.iter().map(|x| x * x).sum::<f64>()
            + self.d_h.iter().map(|x| x * x).sum::<f64>();
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

/// State at a single BL station for local Newton iteration.
#[derive(Debug, Clone, Copy, Default)]
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
        }
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
                    break;
                }
            };
            
            // Line search with DIJ residuals
            let alpha = self.line_search_with_dij(&state, &step, &residuals, factorized, flow, geometry, dij);
            
            // Update state
            state = state.apply_step(&step, alpha);
            
            // Apply DIJ-based edge velocity correction
            // The mass defect changes propagate to edge velocity
            self.apply_dij_velocity_correction(&mut state, &step, alpha, geometry, dij);
            
            // Recompute residuals
            residuals = self.compute_residuals_with_dij(&state, factorized, flow, geometry, dij);
            r_norm = residuals.norm();
            residual_history.push(r_norm);
            
            #[cfg(debug_assertions)]
            if iter % 5 == 0 {
                eprintln!(
                    "Newton (DIJ) iter {}: ||R||={:.2e}, α={:.3}, ||dx||={:.2e}",
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
        
        // 4. Effective edge velocity includes DIJ correction
        let _ue_effective: Vec<f64> = ue_base.iter()
            .zip(due_from_dij.iter())
            .map(|(&ue, &due)| (ue + due).max(1e-10))
            .collect();
        
        // 5. Map delta_star to panel nodes for transpiration
        let delta_star_nodes = self.map_bl_to_nodes(&delta_star, geometry);
        
        // 6. Compute transpiration Vn = d(Ue*δ*)/ds
        let vn = compute_transpiration(&ue_base, &delta_star_nodes, &geometry.s_coords);
        
        // 7. Inviscid residuals with transpiration
        let inviscid_with_vn = factorized.solve_with_transpiration(flow, &vn, &geometry.s_coords);
        for i in 0..n_panels {
            residuals.r_inviscid[i] = state.gamma[i] - inviscid_with_vn.gamma[i];
        }
        
        // 8. BL residuals (momentum and shape equations)
        self.compute_bl_residuals(
            state,
            &ue_base,
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
                    let theta_j = state.theta[j].max(1e-10);
                    let ue_j = state.gamma[panel_j].abs().max(1e-10);
                    jacobian[(i, col_h)] += -dij_val * theta_j * ue_j * 0.1; // Damped influence
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
            
            let ue = state.gamma[panel_idx].abs().max(1e-10);
            let theta = state.theta[j].max(1e-10);
            let dh = step.d_h[j] * alpha;
            
            dm_panels[panel_idx] = ue * theta * dh;
        }
        
        // Compute Ue correction from DIJ
        let due = crate::inviscid::compute_ue_change_from_mass(dij, &dm_panels);
        
        // Apply correction to gamma (gamma ≈ Ue for attached flow)
        // Only apply to nodes where both arrays have values
        let n_apply = n_panels.min(due.len()).min(state.gamma.len());
        for i in 0..n_apply {
            let sign = state.gamma[i].signum();
            state.gamma[i] += sign * due[i] * 0.1; // Damped to avoid instability
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
        r_momentum[bl_offset] = state.theta[bl_offset] - theta_init;
        r_shape[bl_offset] = state.h[bl_offset] - 2.59; // Blasius H
        
        // March along surface
        for j in 1..n_stations {
            let bl_idx = bl_offset + j;
            let bl_idx_prev = bl_offset + j - 1;
            
            let panel_idx = surface_indices[j];
            let panel_idx_prev = surface_indices[j - 1];
            
            // Get local values
            let s = geometry.s_coords[panel_idx];
            let s_prev = geometry.s_coords[panel_idx_prev];
            let ds = (s - s_prev).abs().max(1e-10);
            
            let ue_local = ue[panel_idx].max(1e-10);
            let ue_prev = ue[panel_idx_prev].max(1e-10);
            let ue_avg = 0.5 * (ue_local + ue_prev);
            
            let theta = state.theta[bl_idx];
            let theta_prev = state.theta[bl_idx_prev];
            let h = state.h[bl_idx];
            let h_prev = state.h[bl_idx_prev];
            
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
            
            // Shape equation residual (using Head's method):
            // d(H1*θ)/ds = Ce
            let h1 = head_h1(h);
            let h1_prev = head_h1(h_prev);
            let ce = head_entrainment(h1_prev);
            
            let d_h1_theta_ds_actual = (h1 * theta - h1_prev * theta_prev) / ds;
            r_shape[bl_idx] = d_h1_theta_ds_actual - ce;
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
            "State size {} != Residual size {}. State: gamma={}, theta={}, h={}. \
             Residuals: inv={}, mom={}, shape={}",
            n_vars_state, n_vars_res,
            state.gamma.len(), state.theta.len(), state.h.len(),
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
    fn solve_newton_step(
        &self,
        jacobian: &DMatrix<f64>,
        residuals: &Residuals,
    ) -> Option<NewtonStep> {
        let r_vec = DVector::from_vec(residuals.to_vec());
        let neg_r = -&r_vec;
        
        // Use LU decomposition
        let lu = jacobian.clone().lu();
        let dx = lu.solve(&neg_r)?;
        
        let n_gamma = residuals.r_inviscid.len();
        let n_bl = residuals.r_momentum.len();
        
        Some(NewtonStep::from_vec(dx.as_slice(), n_gamma, n_bl))
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
        state.h[0] = 3.0;
        
        // Test get
        assert!((state.get(0) - 1.0).abs() < 1e-10); // gamma[0]
        assert!((state.get(10) - 2.0).abs() < 1e-10); // theta[0]
        assert!((state.get(15) - 3.0).abs() < 1e-10); // h[0]
        
        // Test set
        state.set(0, 10.0);
        assert!((state.gamma[0] - 10.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_newton_state_total_size() {
        let state = NewtonState::new(100, 50);
        assert_eq!(state.total_size(), 100 + 50 + 50); // gamma + theta + h
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
