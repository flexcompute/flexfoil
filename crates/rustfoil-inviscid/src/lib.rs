//! XFOIL-matching inviscid panel method solver.
//!
//! This crate implements Mark Drela's linear vorticity stream function panel method
//! exactly as in XFOIL, producing identical results for the inviscid flow solution.
//!
//! # Architecture
//!
//! The solver is organized into modules matching XFOIL's structure:
//!
//! - [`geometry`] - Panel geometry processing (NCALC, APCALC, TECALC)
//! - [`influence`] - Influence coefficient calculation (PSILIN)
//! - [`system`] - Matrix assembly and solution (GGCALC)
//! - [`solution`] - Solution computation and forces (SPECAL, CLCALC)
//! - [`stagnation`] - Stagnation point detection (STFIND)
//!
//! # Key Insight: γ = velocity
//!
//! In XFOIL's formulation, the vortex strength γ at each node equals the surface
//! tangential velocity. This allows direct use of γ as the boundary layer edge velocity.
//!
//! # Usage
//!
//! ```rust,ignore
//! use rustfoil_inviscid::{InviscidSolver, FlowConditions};
//! use rustfoil_core::Body;
//!
//! // Create airfoil geometry
//! let body = Body::from_naca("0012", 160)?;
//!
//! // Solve inviscid flow
//! let solver = InviscidSolver::new();
//! let factorized = solver.factorize(&body)?;
//!
//! // Get solution at any angle of attack
//! let flow = FlowConditions::with_alpha_deg(4.0);
//! let solution = factorized.solve_alpha(&flow);
//!
//! println!("CL = {:.4}", solution.cl);
//! ```
//!
//! # References
//!
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - XFOIL 6.99 source: `xpanel.f`, `xoper.f`, `xsolve.f`
//! - Katz & Plotkin, "Low-Speed Aerodynamics", Chapter 11

pub mod geometry;
pub mod influence;
pub mod system;
pub mod solution;
pub mod stagnation;

mod error;

pub use error::InviscidError;
pub use geometry::AirfoilGeometry;
pub use system::FactorizedSystem;
pub use solution::{InviscidSolution, FlowConditions};

use std::f64::consts::PI;

/// 1/(4π) - used for stream function influence coefficients
pub const QOPI: f64 = 0.25 / PI;

/// 1/(2π) - used for TE panel source/vortex influence
pub const HOPI: f64 = 0.5 / PI;

/// Result type for inviscid solver operations.
pub type Result<T> = std::result::Result<T, InviscidError>;

/// Main inviscid solver implementing XFOIL's panel method.
pub struct InviscidSolver {
    // Configuration options (reserved for future use)
}

impl Default for InviscidSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl InviscidSolver {
    /// Create a new inviscid solver.
    pub fn new() -> Self {
        Self {}
    }

    /// Factorize the geometry and solve for the two base solutions (α=0° and α=90°).
    ///
    /// This is XFOIL's key optimization: the expensive O(N³) factorization is done once,
    /// then any angle of attack can be computed in O(N) time by linear combination.
    ///
    /// # Arguments
    ///
    /// * `points` - Airfoil coordinates as (x, y) pairs, ordered counter-clockwise from upper TE
    ///
    /// # Returns
    ///
    /// A `FactorizedSystem` that can efficiently compute solutions at any angle of attack.
    pub fn factorize(&self, points: &[(f64, f64)]) -> Result<FactorizedSystem> {
        // Build geometry from points
        let geom = AirfoilGeometry::from_points(points)?;
        
        // Build and solve the influence coefficient system
        system::build_and_factorize(&geom)
    }

    /// Solve for the inviscid flow at a given angle of attack.
    ///
    /// This is a convenience method that factorizes and solves in one step.
    /// For multiple angles of attack, use `factorize()` followed by `solve_alpha()`.
    pub fn solve(&self, points: &[(f64, f64)], flow: &FlowConditions) -> Result<InviscidSolution> {
        let factorized = self.factorize(points)?;
        Ok(factorized.solve_alpha(flow))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants() {
        assert!((QOPI - 0.25 / PI).abs() < 1e-15);
        assert!((HOPI - 0.5 / PI).abs() < 1e-15);
    }
}
