//! Error types for the inviscid solver.

use thiserror::Error;

/// Errors that can occur during inviscid flow analysis.
#[derive(Error, Debug)]
pub enum InviscidError {
    /// Not enough points to form panels.
    #[error("Insufficient points: need at least 10, got {0}")]
    InsufficientPoints(usize),

    /// Zero or negative chord length.
    #[error("Invalid chord length: {0}")]
    InvalidChord(f64),

    /// Duplicate consecutive points detected.
    #[error("Duplicate points at index {0}")]
    DuplicatePoints(usize),

    /// Matrix is singular and cannot be factorized.
    #[error("Singular influence matrix - check geometry")]
    SingularMatrix,

    /// Panel has zero length.
    #[error("Zero-length panel at index {0}")]
    ZeroLengthPanel(usize),

    /// Spline computation failed.
    #[error("Spline computation failed: {0}")]
    SplineError(String),

    /// Invalid angle of attack.
    #[error("Invalid angle of attack: {0} radians")]
    InvalidAlpha(f64),
}
