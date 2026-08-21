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
//! ## Multi-Element Geometry (Design Direction)
//! XFOIL assumes a single body throughout. This crate instead treats a
//! [`Body`] as a self-contained element — its own contour, panels, and cached
//! geometry — so that a multi-element configuration can be expressed as a
//! collection of bodies rather than as a special case. That representation is
//! the direction the geometry layer is built for; the solver work it implies
//! (paneling a configuration as a whole, inviscid interaction between
//! elements, and viscous treatment of the resulting wakes and gaps) is not
//! implemented, and the current solve path operates on a single body. See the
//! development phases in `README.md` for status.
//!
//! The intended shape of a configuration:
//!
//! ```text
//! // Each element is an independent Body.
//! let slat = Body::from_points("slat", &slat_coords)?;
//! let main = Body::from_points("main", &main_coords)?;
//! let flap = Body::from_points("flap", &flap_coords)?;
//!
//! // Building the collection is supported today; solving the aerodynamic
//! // interaction between its members is not yet.
//! let configuration = vec![slat, main, flap];
//! ```
//!
//! Flap deflection in [`flap`] is XFOIL's plain flap: it rotates a region of
//! one contour about a hinge point, so the result is still a single body, not
//! a separate element.
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
