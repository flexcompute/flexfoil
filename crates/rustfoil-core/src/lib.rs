//! RustFoil Core - Geometry and math primitives for airfoil analysis.
//!
//! This crate provides the foundational types for the RustFoil airfoil
//! analysis engine. It is intentionally dependency-light (only `nalgebra`)
//! to ensure fast compilation and small WASM bundle size.
//!
//! # Modules
//!
//! - [`point`] - 2D points and vectors with convenience functions
//! - [`panel`] - Panel discretization for panel methods
//! - [`body`] - Aerodynamic body representation (airfoils, flaps, etc.)
//! - [`spline`] - Cubic spline interpolation for geometry smoothing
//! - [`error`] - Custom error types for geometry operations
//!
//! # Design Philosophy
//!
//! ## Multi-Body First
//! Unlike XFOIL's single-body assumption, RustFoil is designed from the
//! ground up to support multi-element airfoil configurations:
//!
//! ```rust
//! use rustfoil_core::body::Body;
//! use rustfoil_core::point::point;
//!
//! // Define coordinates for each element
//! let slat_coords = vec![
//!     point(0.0, 0.0), point(0.1, -0.01), point(0.05, 0.02), point(0.0, 0.0)
//! ];
//! let main_coords = vec![
//!     point(1.0, 0.0), point(0.5, -0.05), point(0.0, 0.0), point(0.5, 0.05), point(1.0, 0.0)
//! ];
//! let flap_coords = vec![
//!     point(1.2, -0.05), point(1.1, -0.06), point(1.0, -0.03), point(1.1, -0.02), point(1.2, -0.05)
//! ];
//!
//! // Multi-element configuration
//! let slat = Body::from_points("slat", &slat_coords).unwrap();
//! let main = Body::from_points("main", &main_coords).unwrap();
//! let flap = Body::from_points("flap", &flap_coords).unwrap();
//!
//! let configuration = vec![slat, main, flap];
//! assert_eq!(configuration.len(), 3);
//! ```
//!
//! ## Cached Geometry
//! Panel normals, tangents, and midpoints are computed once at construction
//! time, not in hot loops. This is critical for 60 Hz real-time feedback.
//!
//! ## Coordinate Convention
//! - **X-axis:** Downstream (freestream direction). LE at x≈0, TE at x≈1.
//! - **Y-axis:** Upward. Upper surface has positive y.
//! - **Panel ordering:** Counter-clockwise from TE lower to TE upper.
//!
//! # Example: Creating an Airfoil
//!
//! ```rust
//! use rustfoil_core::body::Body;
//! use rustfoil_core::point::point;
//!
//! // Simple diamond airfoil (for testing)
//! let points = vec![
//!     point(1.0, 0.0),   // Trailing edge
//!     point(0.5, -0.05), // Lower surface
//!     point(0.0, 0.0),   // Leading edge
//!     point(0.5, 0.05),  // Upper surface
//!     point(1.0, 0.0),   // Back to TE (closed contour)
//! ];
//!
//! let airfoil = Body::from_points("diamond", &points).unwrap();
//!
//! println!("Panels: {}", airfoil.n_panels());
//! println!("Chord: {:.3}", airfoil.chord());
//! println!("Arc length: {:.3}", airfoil.arc_length());
//! ```

#![warn(missing_docs)]
#![warn(clippy::all)]

pub mod body;
pub mod error;
pub mod flap;
pub mod naca;
pub mod panel;
pub mod point;
pub mod spline;
pub mod xfoil_spline;

#[cfg(test)]
mod xfoil_spline_test;

// Re-export commonly used types at the crate root
pub use body::Body;
pub use error::GeometryError;
pub use panel::Panel;
pub use point::{point, vec2, Point, Vec2};
pub use spline::{CubicSpline, PanelingParams};
pub use xfoil_spline::XfoilSpline;
