//! Custom error types for RustFoil geometry operations.
//!
//! This module provides a dependency-free error enum for all geometry-related
//! failures. The core crate intentionally avoids external error libraries
//! (like `thiserror`) to keep the physics engine lean for WASM compilation.

use core::fmt;

/// Errors that can occur during geometry construction and manipulation.
///
/// # Design Philosophy
/// These errors represent *programmer errors* or *invalid input data*, not
/// runtime physics failures. For physics failures (e.g., non-convergence),
/// see the solver crate's error types.
#[derive(Debug, Clone, PartialEq)]
pub enum GeometryError {
    /// Not enough points to construct the requested geometry.
    ///
    /// For example, constructing a `Body` requires at least 3 points
    /// to form a closed contour with at least one panel.
    InsufficientPoints {
        /// Minimum number of points required
        required: usize,
        /// Number of points actually provided
        provided: usize,
    },

    /// A panel has zero or near-zero length.
    ///
    /// Degenerate panels cannot have well-defined normals or tangents,
    /// which breaks the panel method influence coefficient calculations.
    DegeneratePanel {
        /// Index of the degenerate panel
        index: usize,
    },

    /// Spline interpolation failed.
    ///
    /// This typically occurs when input points are collinear, coincident,
    /// or have non-monotonic parameter values.
    SplineInterpolationFailed {
        /// Description of what went wrong
        reason: &'static str,
    },

    /// A numerical parameter is out of valid range.
    ///
    /// Examples: negative panel count, angle of attack outside [-90°, 90°].
    InvalidParameter {
        /// Name of the invalid parameter
        name: &'static str,
        /// The invalid value
        value: f64,
    },

    /// Airfoil coordinate file parsing failed.
    ParseError {
        /// Line number where parsing failed
        line: usize,
        /// Description of the parse error
        message: &'static str,
    },

    /// The trailing edge could not be identified.
    ///
    /// This occurs when the geometry doesn't have a clear trailing edge
    /// (e.g., a circle or unclosed contour).
    TrailingEdgeNotFound,
}

impl fmt::Display for GeometryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeometryError::InsufficientPoints { required, provided } => {
                write!(
                    f,
                    "Insufficient points: required {}, provided {}",
                    required, provided
                )
            }
            GeometryError::DegeneratePanel { index } => {
                write!(f, "Degenerate panel at index {} (zero length)", index)
            }
            GeometryError::SplineInterpolationFailed { reason } => {
                write!(f, "Spline interpolation failed: {}", reason)
            }
            GeometryError::InvalidParameter { name, value } => {
                write!(f, "Invalid parameter '{}': {}", name, value)
            }
            GeometryError::ParseError { line, message } => {
                write!(f, "Parse error at line {}: {}", line, message)
            }
            GeometryError::TrailingEdgeNotFound => {
                write!(f, "Could not identify trailing edge in geometry")
            }
        }
    }
}

// Implement std::error::Error
impl std::error::Error for GeometryError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = GeometryError::InsufficientPoints {
            required: 3,
            provided: 2,
        };
        assert_eq!(
            err.to_string(),
            "Insufficient points: required 3, provided 2"
        );

        let err = GeometryError::DegeneratePanel { index: 5 };
        assert_eq!(err.to_string(), "Degenerate panel at index 5 (zero length)");
    }
}
