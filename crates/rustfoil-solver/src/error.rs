//! Solver error types for inviscid and viscous solvers.
//!
//! MERGE NOTE: When merging with flexfoil, combine this with the existing
//! SolverError enum in flexfoil/crates/rustfoil-solver/src/error.rs.
//! The inviscid error variants already exist there; add the viscous variants.

use std::fmt;

/// Error types for solver operations.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    // ========== Inviscid errors (from flexfoil) ==========
    /// No bodies provided to solver
    NoBodies,
    /// Insufficient panels for accurate solution
    InsufficientPanels,
    /// Matrix factorization failed (singular system)
    SingularMatrix,
    /// Solver did not converge within iteration limit
    NonConvergence {
        iterations: usize,
        residual: f64,
    },
    /// Invalid flow conditions specified
    InvalidFlowConditions {
        reason: &'static str,
    },

    // ========== Viscous-specific errors (new) ==========
    /// Boundary layer separated before trailing edge
    BoundaryLayerSeparation {
        /// x/c location where separation occurred
        x_sep: f64,
    },
    /// Transition prediction failed
    TransitionFailure {
        reason: String,
    },
    /// Invalid Reynolds number (must be positive)
    InvalidReynolds,
    /// BL marching failed to converge
    MarchingFailure {
        station: usize,
        reason: String,
    },
    /// Newton system solve failed
    NewtonSolveFailure {
        iteration: usize,
        residual: f64,
    },
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            // Inviscid errors
            SolverError::NoBodies => write!(f, "No bodies provided to solver"),
            SolverError::InsufficientPanels => {
                write!(f, "Insufficient panels for accurate solution")
            }
            SolverError::SingularMatrix => write!(f, "Matrix factorization failed (singular system)"),
            SolverError::NonConvergence {
                iterations,
                residual,
            } => write!(
                f,
                "Solver did not converge after {} iterations (residual: {:.2e})",
                iterations, residual
            ),
            SolverError::InvalidFlowConditions { reason } => {
                write!(f, "Invalid flow conditions: {}", reason)
            }

            // Viscous errors
            SolverError::BoundaryLayerSeparation { x_sep } => {
                write!(f, "Boundary layer separation at x/c = {:.4}", x_sep)
            }
            SolverError::TransitionFailure { reason } => {
                write!(f, "Transition prediction failed: {}", reason)
            }
            SolverError::InvalidReynolds => {
                write!(f, "Invalid Reynolds number (must be positive)")
            }
            SolverError::MarchingFailure { station, reason } => {
                write!(f, "BL marching failed at station {}: {}", station, reason)
            }
            SolverError::NewtonSolveFailure { iteration, residual } => {
                write!(
                    f,
                    "Newton solve failed at iteration {} (residual: {:.2e})",
                    iteration, residual
                )
            }
        }
    }
}

impl std::error::Error for SolverError {}

/// Result type alias for solver operations.
pub type SolverResult<T> = Result<T, SolverError>;
