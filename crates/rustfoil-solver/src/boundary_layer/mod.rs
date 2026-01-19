//! Boundary Layer Module
//!
//! Implements integral boundary layer equations for viscous analysis:
//! - Thwaites method for laminar boundary layers
//! - Multiple turbulent models: Head's entrainment, XFOIL Cτ lag
//! - eN transition prediction
//!
//! See `docs/VISCOUS_MODELS.md` for full mathematical documentation.

mod state;
mod laminar;
mod turbulent;
mod transition;
/// Closure relations for boundary layer equations.
pub mod closure;
/// XFOIL-style Cτ lag-dissipation turbulent model.
pub mod xfoil_turb;
/// Wake boundary layer marching.
pub mod wake;

pub use state::{BLState, BLStation, Surface};
pub use laminar::thwaites_solve;
pub use turbulent::head_solve;
pub use transition::{TransitionInfo, compute_amplification};
pub use closure::{shape_factor_correlations, skin_friction, thwaites_h};
pub use xfoil_turb::{XfoilConstants, xfoil_turb_solve};
pub use wake::{WakeConfig, WakeStation, initialize_wake, march_wake, squire_young_drag, compute_squire_young_drag, compute_wake_drag};

use serde::{Deserialize, Serialize};

/// Available turbulent boundary layer models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum TurbulentModel {
    /// Head's entrainment method (1958).
    /// Two-equation: θ and H₁ (entrainment shape factor).
    /// Fast and stable, good for attached flows.
    #[default]
    Head,
    
    /// XFOIL's Cτ lag-dissipation method (Drela, 1987).
    /// Two-equation: θ and Cτ (shear lag coefficient).
    /// More accurate for separating flows.
    XfoilCtau,
    
    /// Green's lag-entrainment method (1977).
    /// Three-equation: θ, H₁, and Cτ.
    /// Best history effects, more complex.
    GreenLag,
}

use crate::inviscid::InviscidSolution;

/// Complete boundary layer solution for both surfaces.
#[derive(Debug, Clone)]
pub struct BLSolution {
    /// Upper surface BL stations (from stagnation point to TE)
    pub upper: Vec<BLStation>,
    /// Lower surface BL stations (from stagnation point to TE)
    pub lower: Vec<BLStation>,
    /// Transition info for upper surface
    pub transition_upper: TransitionInfo,
    /// Transition info for lower surface  
    pub transition_lower: TransitionInfo,
    /// Critical amplification factor used
    pub n_crit: f64,
    /// Total skin friction drag coefficient
    pub cd_friction: f64,
    /// Total pressure drag coefficient (from momentum deficit)
    pub cd_pressure: f64,
    /// Total drag coefficient
    pub cd: f64,
    /// Whether the BL solution converged
    pub converged: bool,
}

impl BLSolution {
    /// Create a new empty BL solution.
    pub fn new(n_crit: f64) -> Self {
        Self {
            upper: Vec::new(),
            lower: Vec::new(),
            transition_upper: TransitionInfo::default(),
            transition_lower: TransitionInfo::default(),
            n_crit,
            cd_friction: 0.0,
            cd_pressure: 0.0,
            cd: 0.0,
            converged: false,
        }
    }
}

/// Boundary layer solver configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLConfig {
    /// Critical N-factor for transition (default 9.0)
    /// - 9.0: Clean wind tunnel / flight
    /// - 11-14: Very low turbulence
    /// - 4-8: High turbulence / rough surface
    pub n_crit: f64,
    
    /// Reynolds number based on chord
    pub reynolds: f64,
    
    /// Turbulent model selection
    pub turbulent_model: TurbulentModel,
    
    /// XFOIL constants (used if turbulent_model is XfoilCtau)
    #[serde(skip)]
    pub xfoil_constants: XfoilConstants,
    
    /// Maximum iterations for turbulent closure
    pub max_iterations: usize,
    
    /// Convergence tolerance for iterative procedures
    pub tolerance: f64,
    
    /// Force transition at this x/c location (None = natural)
    pub forced_transition_upper: Option<f64>,
    
    /// Force transition at this x/c location (None = natural)
    pub forced_transition_lower: Option<f64>,
}

impl Default for BLConfig {
    fn default() -> Self {
        Self {
            n_crit: 9.0,
            reynolds: 1e6,
            turbulent_model: TurbulentModel::default(),
            xfoil_constants: XfoilConstants::default(),
            max_iterations: 50,
            tolerance: 1e-6,
            forced_transition_upper: None,
            forced_transition_lower: None,
        }
    }
}

impl BLConfig {
    /// Create a config for low turbulence conditions.
    pub fn low_turbulence(reynolds: f64) -> Self {
        Self {
            n_crit: 11.0,
            reynolds,
            ..Default::default()
        }
    }
    
    /// Create a config for high turbulence conditions.
    pub fn high_turbulence(reynolds: f64) -> Self {
        Self {
            n_crit: 5.0,
            reynolds,
            ..Default::default()
        }
    }
    
    /// Create a config using the XFOIL Cτ model.
    pub fn xfoil_model(reynolds: f64) -> Self {
        Self {
            reynolds,
            turbulent_model: TurbulentModel::XfoilCtau,
            ..Default::default()
        }
    }
}

/// Boundary layer solver.
pub struct BLSolver {
    config: BLConfig,
}

impl BLSolver {
    /// Create a new boundary layer solver.
    pub fn new(config: BLConfig) -> Self {
        Self { config }
    }

    /// Solve the boundary layer given an inviscid solution.
    ///
    /// # Arguments
    /// * `inviscid` - The inviscid flow solution (provides edge velocity Ue)
    /// * `s_coords` - Arc-length coordinates along the surface
    /// * `x_coords` - X-coordinates of surface points
    /// * `y_coords` - Y-coordinates of surface points
    /// * `alpha` - Angle of attack in radians
    ///
    /// # Returns
    /// BL solution with theta, delta*, Cf distributions on both surfaces.
    pub fn solve(
        &self,
        inviscid: &InviscidSolution,
        s_coords: &[f64],
        x_coords: &[f64],
        y_coords: &[f64],
        alpha: f64,
    ) -> BLSolution {
        let _n = inviscid.gamma.len();
        let mut solution = BLSolution::new(self.config.n_crit);

        // Find stagnation point (minimum |gamma|, which equals Ue)
        let stag_idx = self.find_stagnation_point(&inviscid.gamma);

        // Edge velocity is |gamma| (from vortex panel method)
        let ue: Vec<f64> = inviscid.gamma.iter().map(|g| g.abs()).collect();

        // March upper surface: stagnation -> TE (indices decreasing in typical ordering)
        solution.upper = self.march_surface(
            Surface::Upper,
            stag_idx,
            &ue,
            s_coords,
            x_coords,
            &mut solution.transition_upper,
        );

        // March lower surface: stagnation -> TE (indices increasing)
        solution.lower = self.march_surface(
            Surface::Lower,
            stag_idx,
            &ue,
            s_coords,
            x_coords,
            &mut solution.transition_lower,
        );

        // Compute drag coefficients using XFOIL method
        self.compute_drag(&mut solution, x_coords, y_coords, alpha);

        solution.converged = true;
        solution
    }

    /// Find the stagnation point index (where Ue is minimum).
    fn find_stagnation_point(&self, gamma: &[f64]) -> usize {
        gamma
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// March the boundary layer along one surface.
    fn march_surface(
        &self,
        surface: Surface,
        stag_idx: usize,
        ue: &[f64],
        s_coords: &[f64],
        x_coords: &[f64],
        transition: &mut TransitionInfo,
    ) -> Vec<BLStation> {
        let n = ue.len();
        let mut stations = Vec::new();

        // Determine marching direction and range
        let (indices, _): (Vec<usize>, i32) = match surface {
            Surface::Upper => {
                // March from stag_idx down to 0 (toward upper TE)
                ((0..=stag_idx).rev().collect(), -1)
            }
            Surface::Lower => {
                // March from stag_idx up to n-1 (toward lower TE)
                ((stag_idx..n).collect(), 1)
            }
        };

        if indices.len() < 2 {
            return stations;
        }

        // Initialize at stagnation point with Thwaites similarity solution
        let mut state = BLState::stagnation_point(self.config.reynolds);
        let mut is_turbulent = false;
        let mut amplification = 0.0;

        // Compute initial s from stagnation
        let s_stag = s_coords[stag_idx];

        // Get consistent bounds
        let n_s = s_coords.len();
        let n_x = x_coords.len();
        let n_ue = ue.len();
        let n_max = n_s.min(n_x).min(n_ue);
        
        for (step, &idx) in indices.iter().enumerate() {
            if idx >= n_max {
                continue;
            }

            let s = (s_coords[idx] - s_stag).abs();
            let x = x_coords[idx];
            let ue_local = ue[idx].max(1e-10); // Prevent division by zero

            // Compute velocity gradient (dUe/ds)
            let n_pts = ue.len();
            let due_ds = if step > 0 && step < indices.len() - 1 {
                let idx_prev = indices[step - 1];
                let idx_next = indices[step + 1];
                // Bounds check
                if idx_prev < n_pts && idx_next < n_pts && idx_prev < s_coords.len() && idx_next < s_coords.len() {
                    let ds = (s_coords[idx_next] - s_coords[idx_prev]).abs();
                    if ds > 1e-12 {
                        (ue[idx_next] - ue[idx_prev]) / ds
                    } else {
                        0.0
                    }
                } else {
                    0.0
                }
            } else {
                0.0
            };

            // Check for forced transition
            let force_tr = match surface {
                Surface::Upper => self.config.forced_transition_upper,
                Surface::Lower => self.config.forced_transition_lower,
            };
            let forced = force_tr.map(|x_tr| x >= x_tr).unwrap_or(false);
            
            if !is_turbulent && !forced {
                // Laminar: use Thwaites method
                let theta_sq_times_6 = thwaites_solve(
                    s,
                    ue_local,
                    self.config.reynolds,
                );
                
                state.theta = (theta_sq_times_6 / 6.0).sqrt().max(1e-10);
                
                // Thwaites pressure gradient parameter
                let lambda = state.theta.powi(2) * due_ds * self.config.reynolds / ue_local;
                
                // Shape factor from Thwaites correlation
                state.h = closure::thwaites_h(lambda);
                state.delta_star = state.h * state.theta;
                
                // Skin friction from Thwaites
                state.cf = closure::thwaites_cf(lambda, state.theta, self.config.reynolds, ue_local);

                // Check for transition using eN method
                let re_theta = ue_local * state.theta * self.config.reynolds;
                
                // Compute step size (ds) for amplification calculation
                let ds_step = if step > 0 {
                    let idx_prev = indices[step - 1];
                    if idx_prev < s_coords.len() {
                        (s_coords[idx] - s_coords[idx_prev]).abs()
                    } else {
                        0.01
                    }
                } else {
                    0.01
                };
                
                amplification = compute_amplification(
                    amplification, 
                    state.h, 
                    state.theta, 
                    re_theta, 
                    ds_step
                );

                if amplification >= self.config.n_crit || forced {
                    is_turbulent = true;
                    transition.x_tr = x;
                    transition.s_tr = s;
                    transition.n_factor = amplification;
                    transition.re_theta_tr = re_theta;
                    
                    // Initialize turbulent state for XFOIL model
                    if self.config.turbulent_model == TurbulentModel::XfoilCtau {
                        state = xfoil_turb::initialize_turbulent(&state, re_theta);
                    }
                }
            } else {
                // Turbulent: use selected model
                state = match self.config.turbulent_model {
                    TurbulentModel::Head => {
                        head_solve(&state, s, ue_local, due_ds, self.config.reynolds)
                    }
                    TurbulentModel::XfoilCtau => {
                        xfoil_turb_solve(
                            &state, s, ue_local, due_ds, 
                            self.config.reynolds, &self.config.xfoil_constants
                        )
                    }
                    TurbulentModel::GreenLag => {
                        // Fall back to Head for now (Green not fully implemented)
                        turbulent::green_lag_solve(
                            &state, s, ue_local, due_ds, 
                            self.config.reynolds, state.ctau
                        )
                    }
                };
            }

            state.s = s;
            state.x = x;
            state.ue = ue_local;
            state.is_turbulent = is_turbulent;

            stations.push(BLStation {
                surface,
                idx,
                state: state.clone(),
            });
        }

        stations
    }

    /// Compute drag coefficients from the BL solution using XFOIL method.
    ///
    /// Friction drag is computed by integrating wall shear stress:
    ///   Cd_f = ∫ τ·dx / (0.5·ρ·U∞²·c)
    /// where τ = 0.5·ρ·Ue²·Cf (wall shear stress)
    ///
    /// Total drag is from wake momentum deficit (Squire-Young):
    ///   Cd = 2·θ_wake · (Ue_wake/U∞)^((5+H_wake)/2)
    ///
    /// Pressure drag = total - friction
    fn compute_drag(
        &self,
        solution: &mut BLSolution,
        x_coords: &[f64],
        y_coords: &[f64],
        alpha: f64,
    ) {
        let cos_alpha = alpha.cos();
        let sin_alpha = alpha.sin();
        
        // Chord length (for normalization)
        let chord = x_coords.iter().cloned().fold(0.0_f64, f64::max)
            - x_coords.iter().cloned().fold(f64::INFINITY, f64::min);
        let chord = chord.max(1e-10);
        
        // U_inf = 1 for normalized flow
        let u_inf_sq = 1.0;
        
        // === Friction Drag (XFOIL method) ===
        // CDF = Σ 0.5*(τ[i] + τ[i-1]) * dx * 2/Q_inf²
        // where τ = 0.5 * Ue² * Cf (wall shear stress, ρ=1)
        // dx = projection onto freestream direction
        
        let mut cd_friction = 0.0;
        
        // Integrate upper surface
        for i in 1..solution.upper.len() {
            let st = &solution.upper[i];
            let st_prev = &solution.upper[i - 1];
            
            // Wall shear stress: τ = 0.5 * ρ * Ue² * Cf
            let tau_i = 0.5 * st.state.ue.powi(2) * st.state.cf;
            let tau_prev = 0.5 * st_prev.state.ue.powi(2) * st_prev.state.cf;
            
            // Get coordinates from indices
            let idx_i = st.idx.min(x_coords.len() - 1);
            let idx_prev = st_prev.idx.min(x_coords.len() - 1);
            
            // Project dx into freestream direction
            let dx = (x_coords[idx_i] - x_coords[idx_prev]) * cos_alpha 
                   + (y_coords[idx_i] - y_coords[idx_prev]) * sin_alpha;
            
            // Trapezoidal integration
            cd_friction += 0.5 * (tau_i + tau_prev) * dx.abs() * 2.0 / u_inf_sq;
        }
        
        // Integrate lower surface
        for i in 1..solution.lower.len() {
            let st = &solution.lower[i];
            let st_prev = &solution.lower[i - 1];
            
            let tau_i = 0.5 * st.state.ue.powi(2) * st.state.cf;
            let tau_prev = 0.5 * st_prev.state.ue.powi(2) * st_prev.state.cf;
            
            let idx_i = st.idx.min(x_coords.len() - 1);
            let idx_prev = st_prev.idx.min(x_coords.len() - 1);
            
            let dx = (x_coords[idx_i] - x_coords[idx_prev]) * cos_alpha 
                   + (y_coords[idx_i] - y_coords[idx_prev]) * sin_alpha;
            
            cd_friction += 0.5 * (tau_i + tau_prev) * dx.abs() * 2.0 / u_inf_sq;
        }
        
        // Normalize by chord
        solution.cd_friction = cd_friction / chord;
        
        // === Total Drag via Wake Marching and Squire-Young ===
        // March the merged BL into the wake, then apply Squire-Young
        // at the wake end for accurate drag prediction
        
        // Get TE states
        let upper_te = solution.upper.last()
            .map(|st| st.state.clone())
            .unwrap_or_default();
        let lower_te = solution.lower.last()
            .map(|st| st.state.clone())
            .unwrap_or_default();
        
        // Get TE x-coordinate
        let x_te = solution.upper.last()
            .map(|st| st.state.x)
            .unwrap_or(1.0);
        
        // Initialize wake at TE (merge upper + lower BLs)
        let wake_te = wake::initialize_wake(&upper_te, &lower_te, x_te);
        
        // March wake downstream
        let wake_config = wake::WakeConfig::default();
        let wake_stations = wake::march_wake(&wake_te, chord, &wake_config);
        
        // Get wake-end values for Squire-Young
        let wake_end = wake_stations.last().unwrap_or(&wake_te);
        
        // Total drag from Squire-Young at wake end
        let cd_total = wake::squire_young_drag(wake_end, chord);
        
        // Pressure drag is the remainder
        solution.cd_pressure = (cd_total - solution.cd_friction).max(0.0);
        solution.cd = cd_total;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bl_config_default() {
        let config = BLConfig::default();
        assert!((config.n_crit - 9.0).abs() < 1e-10);
        assert!((config.reynolds - 1e6).abs() < 1e-10);
    }

    #[test]
    fn test_bl_solution_new() {
        let solution = BLSolution::new(9.0);
        assert!(solution.upper.is_empty());
        assert!(solution.lower.is_empty());
        assert!(!solution.converged);
    }
}
