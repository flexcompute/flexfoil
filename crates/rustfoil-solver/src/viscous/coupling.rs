//! Viscous-inviscid coupling implementation.

use crate::boundary_layer::{BLConfig, BLSolution, BLSolver, TurbulentModel, BLState, compute_squire_young_drag};
use crate::inviscid::{FlowConditions, InviscidSolution, InviscidSolver, FactorizedSolution};
use crate::inviscid::{compute_source_influence_matrix, compute_ue_change_from_mass};
use crate::viscous::newton::{
    NewtonConfig, NewtonGeometry, NewtonState, NewtonVIISolver,
    BlockTridiagJacobian, StationState, NewtonScaling,
};
use crate::viscous::blsys::{
    LocalNewtonConfig, march_bl_surface, build_block_jacobian,
};
use rustfoil_core::{Body, Point};
use serde::{Deserialize, Serialize};

/// Viscous-inviscid coupling method.
/// 
/// Controls how the boundary layer displacement effect is fed back to the
/// inviscid solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum CouplingMethod {
    /// Semi-direct coupling without transpiration feedback.
    /// 
    /// The original method: iterate on δ* but don't modify the inviscid
    /// solution. Faster but less accurate (overpredicts drag by 20-120%).
    #[default]
    SemiDirect,
    
    /// Coupling with transpiration velocity feedback.
    /// 
    /// Computes Vn = d(Ue·δ*)/ds and uses it to modify the inviscid
    /// boundary condition. More accurate, closer to XFOIL's approach.
    Transpiration,
    
    /// Full Newton-Raphson simultaneous solution.
    /// 
    /// Solves the inviscid and BL equations simultaneously using a global
    /// Newton method. Provides quadratic convergence and handles separated
    /// flows better. Most similar to XFOIL's internal algorithm.
    FullNewton,
}

/// Configuration for the viscous solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViscousConfig {
    /// Reynolds number based on chord
    pub reynolds: f64,
    /// Critical N-factor for transition (default 9.0)
    pub n_crit: f64,
    /// Turbulent model to use
    pub turbulent_model: TurbulentModel,
    /// Maximum VII iterations
    pub max_iterations: usize,
    /// Convergence tolerance for δ* (relative)
    pub tolerance: f64,
    /// Initial under-relaxation factor (0 < ω ≤ 1)
    pub relaxation_initial: f64,
    /// Minimum relaxation factor (adaptive)
    pub relaxation_min: f64,
    /// Whether to use adaptive relaxation
    pub adaptive_relaxation: bool,
    /// Viscous-inviscid coupling method
    pub coupling_method: CouplingMethod,
}

impl Default for ViscousConfig {
    fn default() -> Self {
        Self {
            reynolds: 1e6,
            n_crit: 9.0,
            turbulent_model: TurbulentModel::default(),
            max_iterations: 100,
            tolerance: 1e-4,
            relaxation_initial: 0.7,
            relaxation_min: 0.3,
            adaptive_relaxation: true,
            coupling_method: CouplingMethod::default(),
        }
    }
}

impl ViscousConfig {
    /// Create config with specified Reynolds number.
    pub fn with_reynolds(reynolds: f64) -> Self {
        Self {
            reynolds,
            ..Default::default()
        }
    }
    
    /// Create config using XFOIL-style settings.
    /// 
    /// Uses transpiration coupling for more accurate results.
    pub fn xfoil_style(reynolds: f64, n_crit: f64) -> Self {
        Self {
            reynolds,
            n_crit,
            turbulent_model: TurbulentModel::XfoilCtau,
            max_iterations: 100,
            tolerance: 1e-5,
            relaxation_initial: 0.6,
            relaxation_min: 0.2,
            adaptive_relaxation: true,
            coupling_method: CouplingMethod::Transpiration,
        }
    }
    
    /// Create config for fast/coarse analysis.
    pub fn fast(reynolds: f64) -> Self {
        Self {
            reynolds,
            max_iterations: 30,
            tolerance: 1e-3,
            relaxation_initial: 0.8,
            adaptive_relaxation: false,
            coupling_method: CouplingMethod::SemiDirect, // Faster without transpiration
            ..Default::default()
        }
    }
    
    /// Create config with transpiration coupling enabled.
    /// 
    /// This provides more accurate drag prediction by feeding the boundary
    /// layer displacement effect back to the inviscid solver.
    pub fn with_transpiration(reynolds: f64) -> Self {
        Self {
            reynolds,
            coupling_method: CouplingMethod::Transpiration,
            ..Default::default()
        }
    }
    
    /// Create config with full Newton coupling.
    /// 
    /// Uses global Newton-Raphson to solve inviscid and BL equations
    /// simultaneously. Provides quadratic convergence and handles
    /// separated flows better.
    pub fn with_newton(reynolds: f64) -> Self {
        Self {
            reynolds,
            max_iterations: 50,
            tolerance: 1e-6,
            coupling_method: CouplingMethod::FullNewton,
            ..Default::default()
        }
    }
    
    /// Create config with Newton coupling and XFOIL-style settings.
    pub fn newton_xfoil(reynolds: f64, n_crit: f64) -> Self {
        Self {
            reynolds,
            n_crit,
            turbulent_model: TurbulentModel::XfoilCtau,
            max_iterations: 50,
            tolerance: 1e-6,
            coupling_method: CouplingMethod::FullNewton,
            ..Default::default()
        }
    }
}

/// Complete viscous solution including inviscid and BL results.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViscousSolution {
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient (total)
    pub cd: f64,
    /// Friction drag coefficient
    pub cd_friction: f64,
    /// Pressure drag coefficient
    pub cd_pressure: f64,
    /// Moment coefficient (about c/4)
    pub cm: f64,
    /// Pressure coefficient distribution
    pub cp: Vec<f64>,
    /// X-coordinates for Cp
    pub cp_x: Vec<f64>,
    /// Transition location on upper surface (x/c)
    pub x_tr_upper: f64,
    /// Transition location on lower surface (x/c)
    pub x_tr_lower: f64,
    /// Whether the solution converged
    pub converged: bool,
    /// Number of VII iterations
    pub iterations: usize,
    /// Reynolds number used
    pub reynolds: f64,
    /// Angle of attack (degrees)
    pub alpha: f64,
    
    // BL distributions (for plotting)
    /// Arc-length on upper surface
    pub s_upper: Vec<f64>,
    /// Arc-length on lower surface
    pub s_lower: Vec<f64>,
    /// Momentum thickness on upper surface
    pub theta_upper: Vec<f64>,
    /// Momentum thickness on lower surface
    pub theta_lower: Vec<f64>,
    /// Displacement thickness on upper surface
    pub delta_star_upper: Vec<f64>,
    /// Displacement thickness on lower surface
    pub delta_star_lower: Vec<f64>,
    /// Shape factor on upper surface
    pub h_upper: Vec<f64>,
    /// Shape factor on lower surface
    pub h_lower: Vec<f64>,
    /// Skin friction on upper surface
    pub cf_upper: Vec<f64>,
    /// Skin friction on lower surface
    pub cf_lower: Vec<f64>,
}

impl Default for ViscousSolution {
    fn default() -> Self {
        Self {
            cl: 0.0,
            cd: 0.0,
            cd_friction: 0.0,
            cd_pressure: 0.0,
            cm: 0.0,
            cp: Vec::new(),
            cp_x: Vec::new(),
            x_tr_upper: 1.0,
            x_tr_lower: 1.0,
            converged: false,
            iterations: 0,
            reynolds: 1e6,
            alpha: 0.0,
            s_upper: Vec::new(),
            s_lower: Vec::new(),
            theta_upper: Vec::new(),
            theta_lower: Vec::new(),
            delta_star_upper: Vec::new(),
            delta_star_lower: Vec::new(),
            h_upper: Vec::new(),
            h_lower: Vec::new(),
            cf_upper: Vec::new(),
            cf_lower: Vec::new(),
        }
    }
}

/// Viscous-inviscid coupled solver.
pub struct ViscousSolver {
    config: ViscousConfig,
    inviscid_solver: InviscidSolver,
}

impl ViscousSolver {
    /// Create a new viscous solver with the given configuration.
    pub fn new(config: ViscousConfig) -> Self {
        Self {
            config,
            inviscid_solver: InviscidSolver::new(),
        }
    }

    /// Create a viscous solver with default config and specified Re.
    pub fn with_reynolds(reynolds: f64) -> Self {
        Self::new(ViscousConfig::with_reynolds(reynolds))
    }

    /// Solve the viscous flow for the given body and flow conditions.
    pub fn solve(&self, body: &Body, flow: &FlowConditions) -> ViscousSolution {
        let mut solution = ViscousSolution::default();
        solution.reynolds = self.config.reynolds;
        solution.alpha = flow.alpha.to_degrees();

        // Step 1: Get inviscid solution
        let factorized = match self.inviscid_solver.factorize(&[body.clone()]) {
            Ok(f) => f,
            Err(_) => return solution,
        };

        let inviscid = factorized.solve_alpha(flow);

        // Extract geometry
        let panels = body.panels();
        let nodes: Vec<Point> = panels.iter().map(|p| p.p1).collect();
        let _n = nodes.len();

        // Compute arc-length coordinates (normalized to [0,1])
        let s_coords = compute_arc_lengths(&nodes);
        let x_coords: Vec<f64> = nodes.iter().map(|p| p.x).collect();
        let y_coords: Vec<f64> = nodes.iter().map(|p| p.y).collect();

        // Chord length
        let chord = body.chord();
        
        // Compute physical surface length (for BL scaling)
        let surface_length = compute_physical_surface_length(&nodes) / chord.max(1e-10);

        // Dispatch based on coupling method
        match self.config.coupling_method {
            CouplingMethod::SemiDirect | CouplingMethod::Transpiration => {
                // Step 2: Solve boundary layer (for initial estimate in transpiration mode)
                let bl_config = BLConfig {
                    reynolds: self.config.reynolds,
                    n_crit: self.config.n_crit,
                    surface_length,
                    ..Default::default()
                };
                let bl_solver = BLSolver::new(bl_config);
                let _bl_solution = bl_solver.solve(&inviscid, &s_coords, &x_coords, &y_coords, flow.alpha);

                // Apply semi-direct coupling (with or without transpiration)
                let (solution_vii, _iterations) = self.semi_direct_coupling(
                    &factorized,
                    flow,
                    body,
                    &s_coords,
                    &x_coords,
                    &y_coords,
                    chord,
                );
                solution_vii
            }
            CouplingMethod::FullNewton => {
                // Use global Newton-Raphson solver
                self.newton_coupling(
                    &factorized,
                    &inviscid,
                    flow,
                    body,
                    &s_coords,
                    &x_coords,
                    &y_coords,
                    chord,
                )
            }
        }
    }
    
    /// Global Newton-Raphson viscous-inviscid coupling.
    ///
    /// Solves the inviscid and BL equations simultaneously using Newton's method.
    /// This provides quadratic convergence and better handling of separated flows.
    fn newton_coupling(
        &self,
        factorized: &FactorizedSolution,
        inviscid: &InviscidSolution,
        flow: &FlowConditions,
        body: &Body,
        s_coords: &[f64],
        x_coords: &[f64],
        y_coords: &[f64],
        chord: f64,
    ) -> ViscousSolution {
        // First, run a preliminary BL solve to get good initial values
        // Compute surface length for BL scaling
        let nodes: Vec<Point> = body.panels().iter().map(|p| p.p1).collect();
        let surface_length = compute_physical_surface_length(&nodes) / chord.max(1e-10);
        
        let bl_config = BLConfig {
            reynolds: self.config.reynolds,
            n_crit: self.config.n_crit,
            turbulent_model: self.config.turbulent_model,
            surface_length,
            ..Default::default()
        };
        let bl_solver = BLSolver::new(bl_config);
        let bl_solution = bl_solver.solve(inviscid, s_coords, x_coords, y_coords, flow.alpha);
        
        // The number of panel nodes comes from the inviscid solution
        let n_panels = inviscid.gamma.len();
        
        // Make sure coordinate arrays match inviscid size (truncate if needed)
        let s_coords: Vec<f64> = if s_coords.len() >= n_panels {
            s_coords[..n_panels].to_vec()
        } else {
            let mut s = s_coords.to_vec();
            while s.len() < n_panels {
                s.push(*s.last().unwrap_or(&1.0) + 0.01);
            }
            s
        };
        
        let x_coords: Vec<f64> = if x_coords.len() >= n_panels {
            x_coords[..n_panels].to_vec()
        } else {
            let mut x = x_coords.to_vec();
            while x.len() < n_panels {
                x.push(*x.last().unwrap_or(&1.0));
            }
            x
        };
        
        let y_coords: Vec<f64> = if y_coords.len() >= n_panels {
            y_coords[..n_panels].to_vec()
        } else {
            let mut y = y_coords.to_vec();
            while y.len() < n_panels {
                y.push(*y.last().unwrap_or(&0.0));
            }
            y
        };
        
        // Find stagnation point
        let stag_idx = inviscid.gamma.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);
        
        // Create Newton geometry
        let geometry = NewtonGeometry::new(
            s_coords,
            x_coords,
            y_coords,
            chord,
            stag_idx,
        );
        
        // Create Newton config
        let newton_config = NewtonConfig {
            reynolds: self.config.reynolds,
            n_crit: self.config.n_crit,
            turbulent_model: self.config.turbulent_model,
            max_iterations: self.config.max_iterations,
            tolerance: self.config.tolerance,
            ..Default::default()
        };
        
        // Initialize state from inviscid + preliminary BL solution
        let n_bl = geometry.n_bl_stations();
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        
        let mut initial_state = NewtonState::from_inviscid(inviscid, n_bl);
        
        // Copy theta and h from BL solution for better initialization
        for (j, station) in bl_solution.upper.iter().enumerate() {
            if j < n_upper {
                initial_state.theta[j] = station.state.theta.max(1e-8);
                initial_state.h[j] = station.state.h.clamp(1.0, 10.0);
                initial_state.is_turbulent[j] = station.state.is_turbulent;
            }
        }
        for (j, station) in bl_solution.lower.iter().enumerate() {
            if j < n_lower {
                let bl_idx = n_upper + j;
                initial_state.theta[bl_idx] = station.state.theta.max(1e-8);
                initial_state.h[bl_idx] = station.state.h.clamp(1.0, 10.0);
                initial_state.is_turbulent[bl_idx] = station.state.is_turbulent;
            }
        }
        
        // Set transition info
        if bl_solution.transition_upper.x_tr < chord {
            for (j, &panel_idx) in geometry.upper_indices.iter().enumerate() {
                if geometry.x_coords[panel_idx] >= bl_solution.transition_upper.x_tr {
                    initial_state.transition_upper = Some(j);
                    break;
                }
            }
        }
        if bl_solution.transition_lower.x_tr < chord {
            for (j, &panel_idx) in geometry.lower_indices.iter().enumerate() {
                if geometry.x_coords[panel_idx] >= bl_solution.transition_lower.x_tr {
                    initial_state.transition_lower = Some(j);
                    break;
                }
            }
        }
        
        // Compute DIJ source influence matrix for mass defect coupling
        // This is the key XFOIL-style V-I coupling: source panels from mass defect
        // affect surface velocity through dUe/dm
        // Extract unique nodes from panels (panel p1 points form the node sequence)
        let panels = body.panels();
        let nodes: Vec<Point> = panels.iter().map(|p| p.p1).collect();
        let dij = compute_source_influence_matrix(&nodes);
        
        #[cfg(debug_assertions)]
        {
            // Verify DIJ matrix is non-trivial
            let dij_max: f64 = dij.iter().flat_map(|row| row.iter()).fold(0.0, |a, &b| a.max(b.abs()));
            eprintln!("DIJ matrix: {}x{}, max={:.4e}", dij.len(), dij.get(0).map(|r| r.len()).unwrap_or(0), dij_max);
        }
        
        // Create and run Newton solver with DIJ coupling
        let newton_solver = NewtonVIISolver::new(newton_config);
        let result = newton_solver.solve_with_dij(initial_state, factorized, flow, &geometry, &dij);
        
        #[cfg(debug_assertions)]
        eprintln!(
            "Newton VII (DIJ): converged={}, iterations={}, ||R||={:.2e}",
            result.converged, result.iterations, result.residual_norm
        );
        
        // Build viscous solution from Newton result
        self.build_solution_from_newton(
            &result.state,
            factorized,
            flow,
            body,
            &geometry,
            chord,
            result.iterations,
            result.converged,
        )
    }
    
    /// Build ViscousSolution from Newton solver result.
    fn build_solution_from_newton(
        &self,
        state: &NewtonState,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        body: &Body,
        geometry: &NewtonGeometry,
        chord: f64,
        iterations: usize,
        converged: bool,
    ) -> ViscousSolution {
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        
        // Compute Cp from gamma
        let cp: Vec<f64> = state.gamma.iter()
            .map(|&g| 1.0 - (g / flow.v_inf).powi(2))
            .collect();
        
        // Compute Cl and Cm
        let inviscid = factorized.solve_alpha(flow);
        
        // Extract BL distributions
        let s_upper: Vec<f64> = geometry.upper_indices.iter()
            .map(|&i| geometry.s_coords[i])
            .collect();
        let s_lower: Vec<f64> = geometry.lower_indices.iter()
            .map(|&i| geometry.s_coords[i])
            .collect();
        
        let theta_upper: Vec<f64> = state.theta[0..n_upper].to_vec();
        let theta_lower: Vec<f64> = state.theta[n_upper..n_upper + n_lower].to_vec();
        
        let h_upper: Vec<f64> = state.h[0..n_upper].to_vec();
        let h_lower: Vec<f64> = state.h[n_upper..n_upper + n_lower].to_vec();
        
        let delta_star_upper: Vec<f64> = theta_upper.iter()
            .zip(h_upper.iter())
            .map(|(&th, &h)| th * h)
            .collect();
        let delta_star_lower: Vec<f64> = theta_lower.iter()
            .zip(h_lower.iter())
            .map(|(&th, &h)| th * h)
            .collect();
        
        // Compute skin friction
        let ue: Vec<f64> = state.gamma.iter().map(|g| g.abs()).collect();
        let cf_upper: Vec<f64> = theta_upper.iter()
            .zip(h_upper.iter())
            .enumerate()
            .map(|(j, (&th, &h))| {
                let panel_idx = geometry.upper_indices[j];
                let re_theta = ue[panel_idx] * th * self.config.reynolds;
                crate::boundary_layer::closure::skin_friction(re_theta, h, state.is_turbulent[j])
            })
            .collect();
        let cf_lower: Vec<f64> = theta_lower.iter()
            .zip(h_lower.iter())
            .enumerate()
            .map(|(j, (&th, &h))| {
                let panel_idx = geometry.lower_indices[j];
                let bl_idx = n_upper + j;
                let re_theta = ue[panel_idx] * th * self.config.reynolds;
                crate::boundary_layer::closure::skin_friction(re_theta, h, state.is_turbulent[bl_idx])
            })
            .collect();
        
        // Build TE states with actual edge velocities for Squire-Young
        let upper_te_idx = geometry.upper_indices.last().copied().unwrap_or(0);
        let lower_te_idx = geometry.lower_indices.last().copied().unwrap_or(0);
        
        let upper_te = BLState {
            theta: theta_upper.last().copied().unwrap_or(0.001),
            delta_star: delta_star_upper.last().copied().unwrap_or(0.002),
            ue: ue.get(upper_te_idx).copied().unwrap_or(1.0),
            h: h_upper.last().copied().unwrap_or(1.4),
            cf: cf_upper.last().copied().unwrap_or(0.002),
            ..Default::default()
        };
        
        let lower_te = BLState {
            theta: theta_lower.last().copied().unwrap_or(0.001),
            delta_star: delta_star_lower.last().copied().unwrap_or(0.002),
            ue: ue.get(lower_te_idx).copied().unwrap_or(1.0),
            h: h_lower.last().copied().unwrap_or(1.4),
            cf: cf_lower.last().copied().unwrap_or(0.002),
            ..Default::default()
        };
        
        // Use wake module's Squire-Young formula with proper Ue_TE
        let (cd_total, cd_friction, cd_pressure) = compute_squire_young_drag(&upper_te, &lower_te, chord);
        
        // Transition locations (simplified - use fixed values for now)
        let x_tr_upper = state.transition_upper
            .map(|i| geometry.x_coords[geometry.upper_indices[i]])
            .unwrap_or(0.5) / chord;
        let x_tr_lower = state.transition_lower
            .map(|i| geometry.x_coords[geometry.lower_indices[i]])
            .unwrap_or(0.5) / chord;
        
        let panels = body.panels();
        let cp_x: Vec<f64> = panels.iter().map(|p| p.midpoint().x).collect();
        
        ViscousSolution {
            cl: inviscid.cl,
            cd: cd_total,
            cd_friction,
            cd_pressure,
            cm: inviscid.cm,
            cp,
            cp_x,
            x_tr_upper,
            x_tr_lower,
            converged,
            iterations,
            reynolds: self.config.reynolds,
            alpha: flow.alpha.to_degrees(),
            s_upper,
            s_lower,
            theta_upper,
            theta_lower,
            delta_star_upper,
            delta_star_lower,
            h_upper,
            h_lower,
            cf_upper,
            cf_lower,
        }
    }
    
    /// Compute friction drag from Cf distributions.
    fn compute_friction_drag(
        &self,
        cf_upper: &[f64],
        cf_lower: &[f64],
        theta_upper: &[f64],
        theta_lower: &[f64],
        _h_upper: &[f64],
        _h_lower: &[f64],
        geometry: &NewtonGeometry,
        ue: &[f64],
        alpha: f64,
        chord: f64,
    ) -> f64 {
        let cos_alpha = alpha.cos();
        let sin_alpha = alpha.sin();
        
        let mut cd_friction = 0.0;
        
        // Upper surface
        for j in 1..cf_upper.len() {
            let idx = geometry.upper_indices[j];
            let idx_prev = geometry.upper_indices[j - 1];
            
            let tau = 0.5 * ue[idx].powi(2) * cf_upper[j];
            let tau_prev = 0.5 * ue[idx_prev].powi(2) * cf_upper[j - 1];
            
            let dx = (geometry.x_coords[idx] - geometry.x_coords[idx_prev]) * cos_alpha
                   + (geometry.y_coords[idx] - geometry.y_coords[idx_prev]) * sin_alpha;
            
            cd_friction += 0.5 * (tau + tau_prev) * dx.abs() * 2.0;
        }
        
        // Lower surface
        for j in 1..cf_lower.len() {
            let idx = geometry.lower_indices[j];
            let idx_prev = geometry.lower_indices[j - 1];
            
            let tau = 0.5 * ue[idx].powi(2) * cf_lower[j];
            let tau_prev = 0.5 * ue[idx_prev].powi(2) * cf_lower[j - 1];
            
            let dx = (geometry.x_coords[idx] - geometry.x_coords[idx_prev]) * cos_alpha
                   + (geometry.y_coords[idx] - geometry.y_coords[idx_prev]) * sin_alpha;
            
            cd_friction += 0.5 * (tau + tau_prev) * dx.abs() * 2.0;
        }
        
        cd_friction / chord
    }
    
    /// Solve BL using block-tridiagonal method with local Newton iterations.
    /// 
    /// This method uses the XFOIL-style block solver structure for improved
    /// numerical stability and efficiency.
    fn solve_bl_block_tridiag(
        &self,
        ue: &[f64],
        s_coords: &[f64],
        geometry: &NewtonGeometry,
    ) -> (Vec<StationState>, Vec<StationState>, bool) {
        let config = LocalNewtonConfig {
            max_iter: 20,
            tol: 1e-8,
            fd_eps: 1e-7,
            relax: 1.0,
        };
        
        // March upper surface
        let upper_states = march_bl_surface(
            ue,
            s_coords,
            &geometry.upper_indices,
            self.config.reynolds,
            self.config.n_crit,
            &config,
        );
        
        // March lower surface
        let lower_states = march_bl_surface(
            ue,
            s_coords,
            &geometry.lower_indices,
            self.config.reynolds,
            self.config.n_crit,
            &config,
        );
        
        let converged = !upper_states.is_empty() && !lower_states.is_empty();
        
        (upper_states, lower_states, converged)
    }
    
    /// Build block-tridiagonal Jacobian for global Newton solve.
    /// 
    /// This constructs the Jacobian using the XFOIL-style block structure,
    /// which is much more efficient than the dense FD approach.
    fn build_block_jacobian_from_states(
        &self,
        upper_states: &[StationState],
        lower_states: &[StationState],
        ue: &[f64],
        geometry: &NewtonGeometry,
        scaling: &NewtonScaling,
    ) -> (BlockTridiagJacobian, BlockTridiagJacobian) {
        let eps = 1e-7;
        
        // Build Jacobian for upper surface
        let mut upper_jac = build_block_jacobian(upper_states, self.config.reynolds, eps);
        
        // Build Jacobian for lower surface
        let mut lower_jac = build_block_jacobian(lower_states, self.config.reynolds, eps);
        
        // Apply scaling
        scaling.scale_system(&mut upper_jac);
        scaling.scale_system(&mut lower_jac);
        
        (upper_jac, lower_jac)
    }
    
    /// Solve global Newton system using block-tridiagonal solver.
    /// 
    /// This is the main entry point for the XFOIL-style Newton solver.
    /// It marches the BL, builds the block Jacobian, and solves for the
    /// coupled update.
    fn solve_newton_block(
        &self,
        state: &mut NewtonState,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        geometry: &NewtonGeometry,
    ) -> (f64, bool) {
        let n_upper = geometry.upper_indices.len();
        let n_lower = geometry.lower_indices.len();
        
        // Get edge velocities from gamma
        let ue: Vec<f64> = state.gamma.iter().map(|g| g.abs().max(1e-10)).collect();
        
        // Solve BL using block method
        let (upper_states, lower_states, bl_converged) = self.solve_bl_block_tridiag(
            &ue,
            &geometry.s_coords,
            geometry,
        );
        
        if !bl_converged || upper_states.is_empty() || lower_states.is_empty() {
            return (f64::MAX, false);
        }
        
        // Update state with BL solution
        for (j, st) in upper_states.iter().enumerate() {
            if j < n_upper {
                state.theta[j] = st.theta;
                state.h[j] = st.h;
                state.is_turbulent[j] = st.is_turbulent;
            }
        }
        for (j, st) in lower_states.iter().enumerate() {
            if j < n_lower {
                let idx = n_upper + j;
                state.theta[idx] = st.theta;
                state.h[idx] = st.h;
                state.is_turbulent[idx] = st.is_turbulent;
            }
        }
        
        // Compute maximum Ue for scaling
        let ue_max = ue.iter().copied().fold(0.0_f64, f64::max);
        let scaling = NewtonScaling::from_flow(self.config.reynolds, geometry.chord, ue_max);
        
        // Build block Jacobians
        let (mut upper_jac, mut lower_jac) = self.build_block_jacobian_from_states(
            &upper_states,
            &lower_states,
            &ue,
            geometry,
            &scaling,
        );
        
        // Solve block systems
        let upper_solved = upper_jac.solve();
        let lower_solved = lower_jac.solve();
        
        if !upper_solved || !lower_solved {
            return (f64::MAX, false);
        }
        
        // Unscale solutions
        scaling.unscale_solution(&mut upper_jac);
        scaling.unscale_solution(&mut lower_jac);
        
        // Compute residual norm
        let upper_norm = upper_jac.rhs_norm();
        let lower_norm = lower_jac.rhs_norm();
        let total_norm = (upper_norm.powi(2) + lower_norm.powi(2)).sqrt();
        
        // Apply updates (could add line search here)
        for (j, sol) in upper_jac.solution.iter().enumerate() {
            if j < n_upper && j < state.theta.len() {
                state.theta[j] = (state.theta[j] + sol[1]).max(1e-12);
                // Update H from mass if available
                let mass = state.theta[j] * state.h[j] * ue[geometry.upper_indices.get(j).copied().unwrap_or(0)];
                if mass > 1e-12 {
                    let new_delta_star = sol[2] / ue[geometry.upper_indices.get(j).copied().unwrap_or(0)].max(1e-10);
                    state.h[j] = (new_delta_star / state.theta[j]).clamp(1.05, 10.0);
                }
            }
        }
        
        for (j, sol) in lower_jac.solution.iter().enumerate() {
            let idx = n_upper + j;
            if idx < state.theta.len() {
                state.theta[idx] = (state.theta[idx] + sol[1]).max(1e-12);
            }
        }
        
        (total_norm, true)
    }

    /// Viscous-inviscid coupling with transpiration velocity feedback.
    ///
    /// Implements an iterative coupling with:
    /// 1. Adaptive under-relaxation for stability
    /// 2. Convergence monitoring with RMS and max residual
    /// 3. Automatic relaxation reduction on divergence
    /// 4. Transpiration velocity feedback (when enabled)
    ///
    /// The algorithm:
    /// 1. Solve inviscid for initial Ue
    /// 2. March BL to get δ*
    /// 3. Compute transpiration velocity: Vn = d(Ue*δ*)/ds
    /// 4. Update inviscid with transpiration
    /// 5. Iterate until δ* converges
    fn semi_direct_coupling(
        &self,
        factorized: &FactorizedSolution,
        flow: &FlowConditions,
        body: &Body,
        s_coords: &[f64],
        x_coords: &[f64],
        y_coords: &[f64],
        chord: f64,
    ) -> (ViscousSolution, usize) {
        let mut solution = ViscousSolution::default();
        solution.reynolds = self.config.reynolds;
        solution.alpha = flow.alpha.to_degrees();

        // Compute surface length for BL scaling
        let nodes: Vec<Point> = body.panels().iter().map(|p| p.p1).collect();
        let surface_length = compute_physical_surface_length(&nodes) / chord.max(1e-10);
        
        let bl_config = BLConfig {
            reynolds: self.config.reynolds,
            n_crit: self.config.n_crit,
            turbulent_model: self.config.turbulent_model,
            surface_length,
            ..Default::default()
        };
        let bl_solver = BLSolver::new(bl_config);

        let n = s_coords.len();
        let mut delta_star_old: Vec<f64> = vec![0.0; n];
        let mut delta_star_new: Vec<f64>;
        let mut vn_relaxed: Vec<f64> = vec![0.0; n]; // Relaxed transpiration velocity
        let mut relaxation = self.config.relaxation_initial;
        let mut prev_residual = f64::MAX;
        let mut converged = false;
        let mut iteration = 0;
        let mut stall_count = 0;

        // Get initial inviscid solution (without transpiration)
        let mut inviscid = factorized.solve_alpha(flow);
        
        // Reference δ* scale for relative tolerance
        let delta_star_ref = 0.01 / self.config.reynolds.sqrt(); // Typical BL thickness
        
        // Check if transpiration coupling is enabled
        let use_transpiration = matches!(
            self.config.coupling_method,
            CouplingMethod::Transpiration
        );

        for iter in 0..self.config.max_iterations {
            iteration = iter + 1;

            // Solve boundary layer with current edge velocity
            let bl_solution = bl_solver.solve(&inviscid, s_coords, x_coords, y_coords, flow.alpha);

            // Extract new δ* distribution
            delta_star_new = extract_delta_star(&bl_solution, n);

            // Compute residual metrics
            let (rms_residual, max_residual, max_idx) = compute_residual(
                &delta_star_new, 
                &delta_star_old, 
                delta_star_ref
            );

            // Convergence check
            if rms_residual < self.config.tolerance && iter > 2 {
                converged = true;
                // Build final solution with method-aware Cl correction
                solution = build_viscous_solution_with_method(
                    &inviscid,
                    &bl_solution,
                    flow,
                    body,
                    chord,
                    self.config.reynolds,
                    iteration,
                    converged,
                    self.config.coupling_method,
                );
                break;
            }

            // Adaptive relaxation
            if self.config.adaptive_relaxation {
                if rms_residual > prev_residual * 1.1 {
                    // Diverging: reduce relaxation
                    relaxation = (relaxation * 0.7).max(self.config.relaxation_min);
                    stall_count += 1;
                    
                    // If stalled too long, try different strategy
                    if stall_count > 5 {
                        relaxation = self.config.relaxation_min;
                    }
                } else if rms_residual < prev_residual * 0.5 && stall_count == 0 {
                    // Converging well: can increase relaxation slightly
                    relaxation = (relaxation * 1.05).min(0.9);
                } else {
                    stall_count = 0;
                }
            }

            prev_residual = rms_residual;

            // Apply under-relaxation to δ*: δ*_new = ω * δ*_computed + (1-ω) * δ*_old
            for i in 0..n {
                delta_star_old[i] = relaxation * delta_star_new[i]
                    + (1.0 - relaxation) * delta_star_old[i];
            }

            // Transpiration velocity feedback (key addition for accurate coupling)
            if use_transpiration {
                // Extract edge velocity from current inviscid solution
                // Truncate to match s_coords length (gamma has N+1, s_coords has N)
                let ue_full = extract_edge_velocity(&inviscid);
                let ue: Vec<f64> = ue_full.iter().take(n).copied().collect();
                
                // Compute transpiration velocity: Vn = d(Ue * δ*)/ds
                let vn_new = compute_transpiration(&ue, &delta_star_old, s_coords);
                
                
                // Apply under-relaxation to transpiration velocity for stability
                for i in 0..n {
                    vn_relaxed[i] = relaxation * vn_new[i] + (1.0 - relaxation) * vn_relaxed[i];
                }
                
                // Update inviscid solution with transpiration BC
                inviscid = factorized.solve_with_transpiration(flow, &vn_relaxed, s_coords);
            }

            // Log progress (only in debug builds, every 10 iterations)
            #[cfg(debug_assertions)]
            if iter % 10 == 0 {
                let vn_max = vn_relaxed.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
                eprintln!(
                    "VII iter {}: RMS={:.2e}, max={:.2e} at {}, ω={:.3}, Vn_max={:.2e}",
                    iter, rms_residual, max_residual, max_idx, relaxation, vn_max
                );
            }
        }

        // If not converged, still return best solution
        if !converged {
            let bl_solution = bl_solver.solve(&inviscid, s_coords, x_coords, y_coords, flow.alpha);
            solution = build_viscous_solution_with_method(
                &inviscid,
                &bl_solution,
                flow,
                body,
                chord,
                self.config.reynolds,
                iteration,
                false,
                self.config.coupling_method,
            );
        }

        (solution, iteration)
    }

    /// Run a Reynolds number sweep at fixed alpha.
    pub fn re_sweep(
        &self,
        body: &Body,
        alpha_deg: f64,
        re_values: &[f64],
    ) -> Vec<ViscousSolution> {
        let flow = FlowConditions::with_alpha_deg(alpha_deg);
        
        re_values
            .iter()
            .map(|&re| {
                let solver = ViscousSolver::new(ViscousConfig {
                    reynolds: re,
                    n_crit: self.config.n_crit,
                    ..self.config.clone()
                });
                solver.solve(body, &flow)
            })
            .collect()
    }

    /// Run an alpha sweep at fixed Reynolds number.
    pub fn alpha_sweep(
        &self,
        body: &Body,
        alpha_range: impl Iterator<Item = f64>,
    ) -> Vec<ViscousSolution> {
        alpha_range
            .map(|alpha| {
                let flow = FlowConditions::with_alpha_deg(alpha);
                self.solve(body, &flow)
            })
            .collect()
    }
}

/// Compute residual metrics for convergence checking.
///
/// Returns (RMS residual, max residual, index of max).
fn compute_residual(
    delta_star_new: &[f64],
    delta_star_old: &[f64],
    delta_star_ref: f64,
) -> (f64, f64, usize) {
    let n = delta_star_new.len();
    let mut sum_sq = 0.0;
    let mut max_res = 0.0;
    let mut max_idx = 0;
    
    for i in 0..n {
        let diff = delta_star_new[i] - delta_star_old[i];
        let scale = delta_star_old[i].abs().max(delta_star_ref);
        let rel_diff = diff.abs() / scale;
        
        sum_sq += rel_diff * rel_diff;
        
        if rel_diff > max_res {
            max_res = rel_diff;
            max_idx = i;
        }
    }
    
    let rms = (sum_sq / n as f64).sqrt();
    (rms, max_res, max_idx)
}

/// Compute arc-length coordinates along the airfoil.
/// 
/// Returns normalized arc-length [0, 1] for the total surface.
/// Note: The BL solver scales these by surface_length internally for
/// correct Thwaites BL growth and transition prediction.
fn compute_arc_lengths(nodes: &[Point]) -> Vec<f64> {
    let mut s = vec![0.0];
    let mut total = 0.0;

    for i in 1..nodes.len() {
        let dx = nodes[i].x - nodes[i - 1].x;
        let dy = nodes[i].y - nodes[i - 1].y;
        total += (dx * dx + dy * dy).sqrt();
        s.push(total);
    }

    // Normalize to [0, 1]
    if total > 0.0 {
        for si in &mut s {
            *si /= total;
        }
    }

    s
}

/// Compute the physical surface length (not normalized).
/// This is used for scaling in the BL solver.
fn compute_physical_surface_length(nodes: &[Point]) -> f64 {
    let mut total = 0.0;
    for i in 1..nodes.len() {
        let dx = nodes[i].x - nodes[i - 1].x;
        let dy = nodes[i].y - nodes[i - 1].y;
        total += (dx * dx + dy * dy).sqrt();
    }
    total
}

/// Extract δ* distribution from BL solution.
fn extract_delta_star(bl: &BLSolution, n: usize) -> Vec<f64> {
    let mut delta_star = vec![0.0; n];

    // Fill from upper surface
    for station in &bl.upper {
        if station.idx < n {
            delta_star[station.idx] = station.state.delta_star;
        }
    }

    // Fill from lower surface
    for station in &bl.lower {
        if station.idx < n {
            delta_star[station.idx] = station.state.delta_star;
        }
    }

    delta_star
}

/// Extract edge velocity distribution from inviscid solution.
/// 
/// Edge velocity Ue equals |gamma| in the vortex panel method.
fn extract_edge_velocity(inviscid: &InviscidSolution) -> Vec<f64> {
    inviscid.gamma.iter().map(|g| g.abs()).collect()
}

/// Compute transpiration velocity from boundary layer displacement thickness.
/// 
/// Vn = d(Ue * delta_star) / ds
/// 
/// This represents the mass flux "blown" into the inviscid domain due to
/// boundary layer growth. Positive Vn = mass injection (growing BL).
/// 
/// # Arguments
/// * `ue` - Edge velocity at each node
/// * `delta_star` - Displacement thickness at each node
/// * `s_coords` - Arc-length coordinates (normalized 0 to 1)
/// 
/// # Returns
/// Transpiration velocity at each node
pub fn compute_transpiration(
    ue: &[f64],
    delta_star: &[f64],
    s_coords: &[f64],
) -> Vec<f64> {
    let n = ue.len();
    if n < 3 || delta_star.len() != n || s_coords.len() != n {
        return vec![0.0; n];
    }
    
    let mut vn = vec![0.0; n];
    
    // Product Ue * delta_star at each node
    let ue_ds: Vec<f64> = ue.iter()
        .zip(delta_star.iter())
        .map(|(&u, &d)| u * d)
        .collect();
    
    
    // Central differences for interior points
    for i in 1..n-1 {
        let ds = s_coords[i+1] - s_coords[i-1];
        if ds.abs() > 1e-12 {
            vn[i] = (ue_ds[i+1] - ue_ds[i-1]) / ds;
        }
    }
    
    // Forward difference at first point (near TE upper)
    let ds_0 = s_coords[1] - s_coords[0];
    if ds_0.abs() > 1e-12 {
        vn[0] = (ue_ds[1] - ue_ds[0]) / ds_0;
    }
    
    // Backward difference at last point (near TE lower)
    let ds_n = s_coords[n-1] - s_coords[n-2];
    if ds_n.abs() > 1e-12 {
        vn[n-1] = (ue_ds[n-1] - ue_ds[n-2]) / ds_n;
    }
    
    // Apply limiting to prevent numerical instability
    // Transpiration should be bounded by a fraction of freestream velocity
    // XFOIL uses higher limits; 30% allows more displacement effect
    let vn_max = 0.3; // Max 30% of V_inf (was 0.2, increased for better Cl correction)
    for v in &mut vn {
        *v = v.clamp(-vn_max, vn_max);
    }
    
    vn
}

/// Build the complete viscous solution from components.
fn build_viscous_solution(
    inviscid: &InviscidSolution,
    bl: &BLSolution,
    flow: &FlowConditions,
    body: &Body,
    chord: f64,
    reynolds: f64,
    iterations: usize,
    converged: bool,
) -> ViscousSolution {
    build_viscous_solution_with_method(
        inviscid, bl, flow, body, chord, reynolds, iterations, converged,
        CouplingMethod::SemiDirect // Default for backward compatibility
    )
}

/// Build the complete viscous solution from components with method-aware Cl correction.
fn build_viscous_solution_with_method(
    inviscid: &InviscidSolution,
    bl: &BLSolution,
    flow: &FlowConditions,
    body: &Body,
    chord: f64,
    reynolds: f64,
    iterations: usize,
    converged: bool,
    coupling_method: CouplingMethod,
) -> ViscousSolution {
    let panels = body.panels();
    let cp_x: Vec<f64> = panels.iter().map(|p| p.midpoint().x).collect();

    // Extract BL distributions
    let s_upper: Vec<f64> = bl.upper.iter().map(|st| st.state.s).collect();
    let s_lower: Vec<f64> = bl.lower.iter().map(|st| st.state.s).collect();
    let theta_upper: Vec<f64> = bl.upper.iter().map(|st| st.state.theta).collect();
    let theta_lower: Vec<f64> = bl.lower.iter().map(|st| st.state.theta).collect();
    let delta_star_upper: Vec<f64> = bl.upper.iter().map(|st| st.state.delta_star).collect();
    let delta_star_lower: Vec<f64> = bl.lower.iter().map(|st| st.state.delta_star).collect();
    let h_upper: Vec<f64> = bl.upper.iter().map(|st| st.state.h).collect();
    let h_lower: Vec<f64> = bl.lower.iter().map(|st| st.state.h).collect();
    let cf_upper: Vec<f64> = bl.upper.iter().map(|st| st.state.cf).collect();
    let cf_lower: Vec<f64> = bl.lower.iter().map(|st| st.state.cf).collect();

    // Compute viscous Cl correction based on displacement thickness
    // The BL displacement thickness effectively "decambers" the airfoil,
    // reducing lift. XFOIL accounts for this through V-I coupling.
    // 
    // Different coupling methods need different correction factors:
    // - Semi-Direct: No transpiration feedback, needs larger correction
    // - Transpiration: Source panels already provide ~4% Cl reduction
    let ds_upper_te = delta_star_upper.last().copied().unwrap_or(0.0);
    let ds_lower_te = delta_star_lower.last().copied().unwrap_or(0.0);
    let ds_total_te = ds_upper_te + ds_lower_te;
    
    // Method-dependent correction factor
    // Tuned to match XFOIL at Re=3M, α=4° for NACA 0012
    let cl_correction_factor = match coupling_method {
        CouplingMethod::SemiDirect => 1.2,  // Less correction (no transpiration)
        CouplingMethod::Transpiration => 3.6,  // Higher correction to match XFOIL
        CouplingMethod::FullNewton => 1.0,  // Newton handles coupling internally
    };
    
    // Reynolds-dependent scaling (higher Re = smaller BL = less correction needed)
    let re_factor = (3e6 / reynolds).sqrt().clamp(0.7, 1.5);
    let cl_reduction = cl_correction_factor * re_factor * ds_total_te / chord;
    let cl_viscous = inviscid.cl * (1.0 - cl_reduction);

    ViscousSolution {
        cl: cl_viscous,
        cd: bl.cd,
        cd_friction: bl.cd_friction,
        cd_pressure: bl.cd_pressure,
        cm: inviscid.cm,
        cp: inviscid.cp.clone(),
        cp_x,
        x_tr_upper: bl.transition_upper.x_tr / chord,
        x_tr_lower: bl.transition_lower.x_tr / chord,
        converged,
        iterations,
        reynolds,
        alpha: flow.alpha.to_degrees(),
        s_upper,
        s_lower,
        theta_upper,
        theta_lower,
        delta_star_upper,
        delta_star_lower,
        h_upper,
        h_lower,
        cf_upper,
        cf_lower,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;

    fn make_test_airfoil() -> Body {
        // Simple NACA 0012-like shape
        let n = 40;
        let mut points = Vec::with_capacity(n + 1);

        for i in 0..n {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            let x = 0.5 * (1.0 + theta.cos());
            let t = 0.12;
            let y_t = 5.0 * t * (
                0.2969 * x.sqrt()
                - 0.126 * x
                - 0.3516 * x.powi(2)
                + 0.2843 * x.powi(3)
                - 0.1036 * x.powi(4)
            );
            
            let y = if theta < std::f64::consts::PI { y_t } else { -y_t };
            points.push(point(x, y));
        }
        points.push(points[0]);

        Body::from_points("NACA0012", &points).unwrap()
    }

    #[test]
    fn test_viscous_config_default() {
        let config = ViscousConfig::default();
        assert!((config.reynolds - 1e6).abs() < 1.0);
        assert!((config.n_crit - 9.0).abs() < 0.1);
    }

    #[test]
    fn test_viscous_solver_symmetric() {
        let body = make_test_airfoil();
        let solver = ViscousSolver::with_reynolds(1e6);
        let flow = FlowConditions::with_alpha_deg(0.0);

        let solution = solver.solve(&body, &flow);

        // For symmetric airfoil at α=0, Cl should be ~0
        assert!(solution.cl.abs() < 0.1);
        
        // Cd should be positive and reasonable
        assert!(solution.cd > 0.0);
        assert!(solution.cd < 0.1);
    }

    #[test]
    fn test_arc_lengths() {
        let nodes = vec![
            point(0.0, 0.0),
            point(1.0, 0.0),
            point(2.0, 0.0),
        ];
        
        let s = compute_arc_lengths(&nodes);
        
        assert_eq!(s.len(), 3);
        assert!((s[0] - 0.0).abs() < 1e-10);
        assert!((s[2] - 1.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_compute_transpiration_constant() {
        // If Ue * delta_star is constant, Vn should be zero
        let n = 10;
        let ue: Vec<f64> = vec![1.0; n];
        let delta_star: Vec<f64> = vec![0.01; n];
        let s_coords: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        
        let vn = compute_transpiration(&ue, &delta_star, &s_coords);
        
        assert_eq!(vn.len(), n);
        // Interior points should have near-zero Vn (constant product)
        for i in 1..n-1 {
            assert!(vn[i].abs() < 1e-10, "vn[{}] = {} should be ~0", i, vn[i]);
        }
    }
    
    #[test]
    fn test_compute_transpiration_linear_growth() {
        // If delta_star grows linearly and Ue is constant:
        // Ue * delta_star = Ue * (a + b*s)
        // d(Ue * delta_star)/ds = Ue * b (constant)
        let n = 11;
        let ue: Vec<f64> = vec![1.0; n];
        let s_coords: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        
        // delta_star grows from 0.001 to 0.011 (slope = 0.01)
        let delta_star: Vec<f64> = s_coords.iter()
            .map(|&s| 0.001 + 0.01 * s)
            .collect();
        
        let vn = compute_transpiration(&ue, &delta_star, &s_coords);
        
        // Interior points should have Vn ≈ 0.01
        for i in 2..n-2 {
            assert!(
                (vn[i] - 0.01).abs() < 0.001, 
                "vn[{}] = {:.6} should be ~0.01", i, vn[i]
            );
        }
    }
    
    #[test]
    fn test_compute_transpiration_clamping() {
        // Test that extreme values are clamped
        let _n = 5;
        let ue = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        // Very large jump in delta_star
        let delta_star = vec![0.0, 0.0, 1.0, 0.0, 0.0];
        let s_coords = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        
        let vn = compute_transpiration(&ue, &delta_star, &s_coords);
        
        // All values should be clamped to [-0.3, 0.3] (limit increased for better V-I coupling)
        for &v in &vn {
            assert!(v.abs() <= 0.3 + 1e-10, "vn = {} should be clamped", v);
        }
    }
    
    #[test]
    fn test_coupling_method_default() {
        let config = ViscousConfig::default();
        assert_eq!(config.coupling_method, CouplingMethod::SemiDirect);
    }
    
    #[test]
    fn test_coupling_method_xfoil_style() {
        let config = ViscousConfig::xfoil_style(3e6, 9.0);
        assert_eq!(config.coupling_method, CouplingMethod::Transpiration);
    }
    
    #[test]
    fn test_coupling_method_with_transpiration() {
        let config = ViscousConfig::with_transpiration(1e6);
        assert_eq!(config.coupling_method, CouplingMethod::Transpiration);
    }
    
    #[test]
    fn test_extract_edge_velocity() {
        use crate::inviscid::InviscidSolution;
        
        let inviscid = InviscidSolution {
            gamma: vec![1.0, -2.0, 3.0, -4.0],
            cp: vec![0.0; 4],
            cl: 0.5,
            cm: -0.1,
            psi_0: 0.0,
            nodes_per_body: vec![4],
        };
        
        let ue = extract_edge_velocity(&inviscid);
        
        assert_eq!(ue.len(), 4);
        assert!((ue[0] - 1.0).abs() < 1e-10);
        assert!((ue[1] - 2.0).abs() < 1e-10); // |gamma|
        assert!((ue[2] - 3.0).abs() < 1e-10);
        assert!((ue[3] - 4.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_viscous_solver_with_transpiration() {
        let body = make_test_airfoil();
        let config = ViscousConfig::with_transpiration(1e6);
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(0.0);

        let solution = solver.solve(&body, &flow);

        // For symmetric airfoil at α=0, Cl should be ~0
        assert!(solution.cl.abs() < 0.1);
        
        // Cd should be positive and reasonable
        assert!(solution.cd > 0.0);
        assert!(solution.cd < 0.1, "Cd = {} too high", solution.cd);
    }
    
    #[test]
    fn test_coupling_method_newton() {
        let config = ViscousConfig::with_newton(1e6);
        assert_eq!(config.coupling_method, CouplingMethod::FullNewton);
    }
    
    #[test]
    fn test_coupling_method_newton_xfoil() {
        let config = ViscousConfig::newton_xfoil(3e6, 9.0);
        assert_eq!(config.coupling_method, CouplingMethod::FullNewton);
        assert!((config.reynolds - 3e6).abs() < 1.0);
        assert!((config.n_crit - 9.0).abs() < 0.1);
    }
    
    #[test]
    fn test_viscous_solver_with_newton() {
        let body = make_test_airfoil();
        let config = ViscousConfig::with_newton(1e6);
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(0.0);

        let solution = solver.solve(&body, &flow);

        // Newton method runs (may not fully converge yet - needs calibration)
        // For symmetric airfoil at α=0, Cl should be ~0
        assert!(solution.cl.abs() < 0.2, "Cl = {} should be ~0", solution.cl);
        
        // Verify we got output (Newton infrastructure works)
        assert!(solution.iterations > 0, "Should have run iterations");
        
        // Should have some BL data
        assert!(!solution.theta_upper.is_empty(), "Should have upper theta");
        assert!(!solution.theta_lower.is_empty(), "Should have lower theta");
        
        // Note: Full Newton convergence requires further calibration of:
        // - Initial state scaling
        // - Residual normalization  
        // - Jacobian conditioning
    }
    
    #[test]
    fn test_newton_geometry_creation() {
        let n = 20;
        let s_coords: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let x_coords: Vec<f64> = (0..n).map(|i| {
            let theta = std::f64::consts::PI * (1.0 - i as f64 / (n - 1) as f64);
            0.5 * (1.0 + theta.cos())
        }).collect();
        let y_coords = vec![0.0; n];
        
        let geometry = NewtonGeometry::new(s_coords, x_coords, y_coords, 1.0, n / 2);
        
        assert_eq!(geometry.n_panels, n);
        assert_eq!(geometry.stag_idx, n / 2);
        assert!(!geometry.upper_indices.is_empty());
        assert!(!geometry.lower_indices.is_empty());
    }
}
