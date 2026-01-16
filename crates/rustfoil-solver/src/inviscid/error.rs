//! Error types for the inviscid solver.

use core::fmt;

/// Errors that can occur during inviscid flow solution.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    /// No bodies were provided to the solver.
    NoBodies,

    /// Insufficient panels for a valid solution (need at least 3).
    InsufficientPanels,

    /// The influence coefficient matrix is singular.
    ///
    /// This typically indicates:
    /// - Degenerate geometry (panels overlapping)
    /// - Numerical issues with very small panels
    SingularMatrix,

    /// The solver failed to converge (for iterative methods).
    NonConvergence {
        /// Number of iterations performed
        iterations: usize,
        /// Final residual norm
        residual: f64,
    },

    /// Invalid flow conditions were specified.
    InvalidFlowConditions {
        /// Description of the problem
        reason: &'static str,
    },
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverError::NoBodies => {
                write!(f, "No bodies provided to solver")
            }
            SolverError::InsufficientPanels => {
                write!(f, "Insufficient panels (need at least 3)")
            }
            SolverError::SingularMatrix => {
                write!(f, "Influence matrix is singular (check geometry)")
            }
            SolverError::NonConvergence {
                iterations,
                residual,
            } => {
                write!(
                    f,
                    "Solver did not converge after {} iterations (residual: {:.2e})",
                    iterations, residual
                )
            }
            SolverError::InvalidFlowConditions { reason } => {
                write!(f, "Invalid flow conditions: {}", reason)
            }
        }
    }
}

impl std::error::Error for SolverError {}
