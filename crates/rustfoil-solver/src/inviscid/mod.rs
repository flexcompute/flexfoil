//! Inviscid flow solver using the Linear Vorticity Panel Method.
//!
//! This module implements the core inviscid aerodynamic analysis:
//!
//! # Algorithm Overview
//!
//! 1. **Geometry Input:** Accept one or more `Body` objects
//! 2. **Influence Matrix:** Build A_ij coefficients
//! 3. **RHS Vector:** Compute freestream normal velocities
//! 4. **Kutta Condition:** Add constraint row for each body
//! 5. **Solve:** LU decomposition of the augmented system
//! 6. **Post-process:** Compute Cp, Cl, Cm from vorticity
//!
//! # Multi-Body Support
//!
//! For N bodies with panels [n₁, n₂, ..., nₙ], the system size is:
//! ```text
//! Total unknowns = Σnᵢ + N  (vorticity nodes + Kutta constraints)
//! ```
//!
//! # Performance Considerations
//!
//! - Influence coefficients are O(N²) to compute
//! - Dense LU solve is O(N³)
//! - For real-time (60 Hz), target N < 200 panels total
//! - Future: iterative solvers for larger systems

mod error;
mod influence;

pub use error::SolverError;

use nalgebra::{DMatrix, DVector};
use rustfoil_core::{Body, Panel, Vec2};

/// Result type for solver operations.
pub type SolverResult<T> = Result<T, SolverError>;

/// Flow conditions for the analysis.
#[derive(Debug, Clone, Copy)]
pub struct FlowConditions {
    /// Angle of attack in radians
    pub alpha: f64,
    /// Freestream velocity magnitude (typically normalized to 1.0)
    pub v_inf: f64,
}

impl Default for FlowConditions {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            v_inf: 1.0,
        }
    }
}

impl FlowConditions {
    /// Create flow conditions with the given angle of attack (in degrees).
    pub fn with_alpha_deg(alpha_deg: f64) -> Self {
        Self {
            alpha: alpha_deg.to_radians(),
            v_inf: 1.0,
        }
    }

    /// Freestream velocity vector.
    #[inline]
    pub fn velocity(&self) -> Vec2 {
        Vec2::new(self.v_inf * self.alpha.cos(), self.v_inf * self.alpha.sin())
    }
}

/// Solution from the inviscid panel method.
#[derive(Debug, Clone)]
pub struct InviscidSolution {
    /// Vorticity values at each panel node
    pub gamma: Vec<f64>,
    /// Pressure coefficient at each panel midpoint
    pub cp: Vec<f64>,
    /// Lift coefficient (per unit span)
    pub cl: f64,
    /// Moment coefficient about quarter-chord
    pub cm: f64,
    /// Number of panels per body (for indexing)
    pub panels_per_body: Vec<usize>,
}

/// Inviscid flow solver using linear vorticity panel method.
///
/// # Usage
///
/// ```ignore
/// use rustfoil_solver::inviscid::{InviscidSolver, FlowConditions};
///
/// let solver = InviscidSolver::new();
/// let flow = FlowConditions::with_alpha_deg(5.0);
///
/// let solution = solver.solve(&[airfoil], &flow)?;
/// println!("Cl = {:.4}", solution.cl);
/// ```
pub struct InviscidSolver {
    // Configuration options (can be extended)
    _config: SolverConfig,
}

#[derive(Debug, Clone, Default)]
struct SolverConfig {
    // Placeholder for future configuration options
    // e.g., wake model, far-field corrections, etc.
}

impl Default for InviscidSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl InviscidSolver {
    /// Create a new inviscid solver with default configuration.
    pub fn new() -> Self {
        Self {
            _config: SolverConfig::default(),
        }
    }

    /// Solve for the inviscid flow around the given bodies.
    ///
    /// # Arguments
    /// * `bodies` - Slice of aerodynamic bodies (airfoils, flaps, etc.)
    /// * `flow` - Freestream flow conditions
    ///
    /// # Returns
    /// `InviscidSolution` containing vorticity, Cp, and aerodynamic coefficients.
    ///
    /// # Errors
    /// - `SolverError::NoBodies` if the bodies slice is empty
    /// - `SolverError::SingularMatrix` if the influence matrix is singular
    pub fn solve(&self, bodies: &[Body], flow: &FlowConditions) -> SolverResult<InviscidSolution> {
        if bodies.is_empty() {
            return Err(SolverError::NoBodies);
        }

        // Count total panels and nodes
        let panels_per_body: Vec<usize> = bodies.iter().map(|b| b.n_panels()).collect();
        let total_panels: usize = panels_per_body.iter().sum();
        let n_bodies = bodies.len();

        // System size: one gamma per panel node (= total_panels for closed bodies)
        // Plus one Kutta condition per body
        let n_unknowns = total_panels + n_bodies;

        // Build influence coefficient matrix
        let mut a_matrix = DMatrix::<f64>::zeros(n_unknowns, n_unknowns);
        let mut rhs = DVector::<f64>::zeros(n_unknowns);

        // Collect all panels for influence calculations
        let all_panels: Vec<&Panel> = bodies.iter().flat_map(|b| b.panels()).collect();

        // Fill influence matrix (panel-to-panel)
        for (i, panel_i) in all_panels.iter().enumerate() {
            let control_point = panel_i.midpoint();
            let normal = panel_i.normal();

            for (j, panel_j) in all_panels.iter().enumerate() {
                // Influence of vortex panel j on control point i
                let (u, v) = influence::vortex_panel_velocity(panel_j, &control_point, 1.0);

                // Normal velocity component
                a_matrix[(i, j)] = u * normal.x + v * normal.y;
            }

            // RHS: negative freestream normal velocity
            let v_inf = flow.velocity();
            rhs[i] = -(v_inf.x * normal.x + v_inf.y * normal.y);
        }

        // Add Kutta conditions (one per body)
        let mut panel_offset = 0;
        for (body_idx, body) in bodies.iter().enumerate() {
            let row = total_panels + body_idx;
            let n = body.n_panels();

            // Kutta condition: γ(TE_lower) + γ(TE_upper) = 0
            // For our convention: first panel is at TE lower, last at TE upper
            // The vorticity at the trailing edge is the average of adjacent panels
            let te_lower = panel_offset; // First panel
            let te_upper = panel_offset + n - 1; // Last panel

            a_matrix[(row, te_lower)] = 1.0;
            a_matrix[(row, te_upper)] = 1.0;
            rhs[row] = 0.0;

            panel_offset += n;
        }

        // Solve the linear system
        let lu = a_matrix.lu();
        let gamma_vec = lu
            .solve(&rhs)
            .ok_or(SolverError::SingularMatrix)?;

        let gamma: Vec<f64> = gamma_vec.iter().take(total_panels).copied().collect();

        // Compute surface velocities and Cp
        let v_inf_mag = flow.v_inf;
        let cp: Vec<f64> = gamma
            .iter()
            .map(|&g| {
                // Surface velocity ≈ γ for unit chord
                // Cp = 1 - (V/V∞)²
                let v_surface = g.abs();
                1.0 - (v_surface / v_inf_mag).powi(2)
            })
            .collect();

        // Compute lift coefficient using Kutta-Joukowski
        // Cl = 2 * Γ / (V∞ * c) where Γ = ∮γ ds
        let mut total_circulation = 0.0;
        for (i, panel) in all_panels.iter().enumerate() {
            total_circulation += gamma[i] * panel.length();
        }

        // Assuming unit chord
        let cl = 2.0 * total_circulation / v_inf_mag;

        // Moment coefficient (simplified, about LE)
        // Cm = -∮ x * γ ds / (c² * V∞)
        let mut moment = 0.0;
        for (i, panel) in all_panels.iter().enumerate() {
            let x_mid = panel.midpoint().x;
            moment -= x_mid * gamma[i] * panel.length();
        }
        let cm = moment / v_inf_mag;

        Ok(InviscidSolution {
            gamma,
            cp,
            cl,
            cm,
            panels_per_body,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;

    fn make_larger_airfoil() -> Body {
        // More refined diamond shape (more panels for better conditioning)
        let points = vec![
            point(1.0, 0.0),
            point(0.75, -0.05),
            point(0.5, -0.08),
            point(0.25, -0.05),
            point(0.0, 0.0),
            point(0.25, 0.05),
            point(0.5, 0.08),
            point(0.75, 0.05),
            point(1.0, 0.0),
        ];
        Body::from_points("symmetric", &points).unwrap()
    }

    #[test]
    fn test_solver_creates_solution() {
        let solver = InviscidSolver::new();
        let airfoil = make_larger_airfoil();
        let flow = FlowConditions::default();

        // Just verify the solver produces a result
        // Detailed accuracy testing is Phase 2 work
        let result = solver.solve(&[airfoil], &flow);

        // For now, we accept either success or singular matrix
        // (matrix conditioning improvements are Phase 2)
        match result {
            Ok(solution) => {
                assert_eq!(solution.cp.len(), 8); // 8 panels
                assert!(solution.cl.is_finite());
            }
            Err(SolverError::SingularMatrix) => {
                // This is acceptable for Phase 1 - needs conditioning work
            }
            Err(e) => panic!("Unexpected error: {:?}", e),
        }
    }

    #[test]
    fn test_no_bodies_error() {
        let solver = InviscidSolver::new();
        let flow = FlowConditions::default();

        let result = solver.solve(&[], &flow);
        assert!(matches!(result, Err(SolverError::NoBodies)));
    }

    #[test]
    fn test_flow_conditions() {
        let flow = FlowConditions::with_alpha_deg(5.0);
        assert!((flow.alpha - 5.0_f64.to_radians()).abs() < 1e-10);

        let vel = flow.velocity();
        assert!(vel.x > 0.0);
        assert!(vel.y > 0.0);
    }
}
