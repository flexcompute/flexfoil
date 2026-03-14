//! Error types for inviscid and viscous solvers.

use core::fmt;

/// Errors that can occur during inviscid or viscous flow solution.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    // ========== Inviscid errors ==========
    
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

    // ========== Viscous errors ==========
    
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
