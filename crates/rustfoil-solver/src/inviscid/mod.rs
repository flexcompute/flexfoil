//! Stub inviscid module - will be replaced by flexfoil's implementation.
//!
//! MERGE NOTE: Delete this entire module after merging with flexfoil.
//! The real implementation lives in flexfoil/crates/rustfoil-solver/src/inviscid/.
//!
//! This stub provides type-compatible interfaces for compilation and testing
//! of the viscous module before the repos are merged.

use crate::error::{SolverError, SolverResult};

/// Flow conditions for inviscid analysis.
///
/// MERGE NOTE: Replace with flexfoil's FlowConditions which has additional
/// fields for Mach number, reference quantities, etc.
#[derive(Debug, Clone, Copy)]
pub struct FlowConditions {
    /// Angle of attack in radians
    pub alpha: f64,
    /// Freestream velocity magnitude
    pub v_inf: f64,
}

impl FlowConditions {
    /// Create flow conditions with angle of attack in degrees.
    pub fn with_alpha_deg(alpha_deg: f64) -> Self {
        Self {
            alpha: alpha_deg.to_radians(),
            v_inf: 1.0,
        }
    }

    /// Create flow conditions with angle of attack in radians.
    pub fn with_alpha_rad(alpha_rad: f64) -> Self {
        Self {
            alpha: alpha_rad,
            v_inf: 1.0,
        }
    }
}

impl Default for FlowConditions {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            v_inf: 1.0,
        }
    }
}

/// Result of inviscid panel method solution.
///
/// MERGE NOTE: Replace with flexfoil's InviscidSolution which includes
/// additional data like pressure coefficient at control points, etc.
#[derive(Debug, Clone)]
pub struct InviscidSolution {
    /// Circulation/surface velocity at each panel (gamma)
    /// This IS the edge velocity Ue that the boundary layer needs.
    pub gamma: Vec<f64>,
    /// Pressure coefficient at each panel
    pub cp: Vec<f64>,
    /// Lift coefficient
    pub cl: f64,
    /// Moment coefficient (about quarter-chord)
    pub cm: f64,
}

impl InviscidSolution {
    /// Create a stub solution for testing.
    ///
    /// This generates a simplified velocity distribution for BL testing.
    /// The real solution comes from the panel method in flexfoil.
    #[allow(dead_code)]
    pub fn stub(n_panels: usize, alpha: f64) -> Self {
        // Generate a simplified velocity distribution for a flat plate
        // Real implementation uses linear vorticity panel method
        let mut gamma = Vec::with_capacity(n_panels);
        let mut cp = Vec::with_capacity(n_panels);

        for i in 0..n_panels {
            // Simple model: velocity varies around the airfoil
            let theta = std::f64::consts::PI * (i as f64) / (n_panels as f64 - 1.0);
            let ue = 1.0 + 0.5 * alpha.sin() * theta.sin();
            gamma.push(ue);
            cp.push(1.0 - ue * ue);
        }

        // Approximate lift from thin airfoil theory
        let cl = 2.0 * std::f64::consts::PI * alpha;

        Self {
            gamma,
            cp,
            cl,
            cm: -cl / 4.0, // Approximate moment
        }
    }
}

/// Inviscid panel method solver.
///
/// MERGE NOTE: Replace with flexfoil's InviscidSolver which implements
/// the Linear Vorticity Panel Method with Kutta condition.
pub struct InviscidSolver;

impl InviscidSolver {
    /// Create a new inviscid solver.
    pub fn new() -> Self {
        Self
    }

    /// Solve for a given angle of attack.
    ///
    /// STUB: Returns a simplified solution for BL testing.
    /// Real implementation factorizes the panel influence matrix once,
    /// then solves for any alpha in O(N) time.
    #[allow(dead_code)]
    pub fn solve(&self, n_panels: usize, flow: &FlowConditions) -> SolverResult<InviscidSolution> {
        if n_panels < 10 {
            return Err(SolverError::InsufficientPanels);
        }
        Ok(InviscidSolution::stub(n_panels, flow.alpha))
    }
}

impl Default for InviscidSolver {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flow_conditions_default() {
        let flow = FlowConditions::default();
        assert_eq!(flow.alpha, 0.0);
        assert_eq!(flow.v_inf, 1.0);
    }

    #[test]
    fn test_flow_conditions_with_alpha() {
        let flow = FlowConditions::with_alpha_deg(5.0);
        assert!((flow.alpha - 5.0_f64.to_radians()).abs() < 1e-10);
    }

    #[test]
    fn test_stub_solution() {
        let sol = InviscidSolution::stub(100, 0.1);
        assert_eq!(sol.gamma.len(), 100);
        assert_eq!(sol.cp.len(), 100);
        // CL should be approximately 2*pi*alpha for small angles
        assert!((sol.cl - 2.0 * std::f64::consts::PI * 0.1).abs() < 0.01);
    }
}
