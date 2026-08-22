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
//! - [`placement`] - Rigid placement of one element within a configuration
//! - [`configuration`] - Multi-element configurations and reference quantities
//! - [`layout`] - Node numbering across a configuration's elements
//! - [`paneling`] - Per-element paneling of a configuration
//! - [`clearance`] - Clearance, gap and overlap between elements
//! - [`config_io`] - Configuration import, export and serialisation
//!
//! # Design Philosophy
//!
//! ## Multi-Element Geometry (Design Direction)
//! XFOIL assumes a single body throughout. This crate instead treats a
//! [`Body`] as a self-contained element — its own contour, panels, and cached
//! geometry — and a [`Configuration`] as an ordered set of placed elements.
//! That representation is the direction the geometry layer is built for. The
//! geometry side of it exists: elements can be placed, paneled per element,
//! checked for clearance, and serialised. The *solver* work it implies —
//! inviscid interaction between elements, and viscous treatment of the
//! resulting wakes — is not implemented, and the current solve path operates on
//! a single body. See the development phases in `README.md` for status.
//!
//! The shape of a configuration:
//!
//! ```rust
//! use rustfoil_core::{Body, Configuration, Element, Layout, Placement};
//! use rustfoil_core::point::point;
//!
//! // Placeholder contour, one per element.
//! let coords = |chord: f64| vec![
//!     point(chord, 0.0),
//!     point(0.5 * chord, -0.05 * chord),
//!     point(0.0, 0.0),
//!     point(0.5 * chord, 0.05 * chord),
//!     point(chord, 0.0),
//! ];
//!
//! // Each element is an independent Body, positioned by a Placement rather
//! // than by baking the position into its coordinates.
//! let mut slat = Element::from_body(Body::from_points("slat", &coords(0.15)).unwrap());
//! slat.placement = Placement::from_translation(-0.12, 0.03);
//! let main = Element::from_body(Body::from_points("main", &coords(1.0)).unwrap());
//! let mut flap = Element::from_body(Body::from_points("flap", &coords(0.3)).unwrap());
//! flap.placement = Placement::rotation_about(point(0.0, 0.0), -30.0);
//!
//! let configuration = Configuration::new(vec![slat, main, flap]);
//! assert_eq!(configuration.len(), 3);
//!
//! // The node numbering the elements share, with per-element panel closure.
//! let layout = Layout::from_configuration(&configuration).unwrap();
//! assert_eq!(layout.n_elements(), 3);
//!
//! // Building and numbering a configuration is supported today; solving the
//! // aerodynamic interaction between its members is not yet.
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
pub mod clearance;
pub mod config_io;
pub mod configuration;
pub mod error;
pub mod flap;
pub mod layout;
pub mod naca;
pub mod panel;
pub mod paneling;
pub mod placement;
pub mod point;
pub mod spline;
pub mod xfoil_spline;

#[cfg(test)]
mod xfoil_spline_test;

// Re-export commonly used types at the crate root
pub use body::{contour_is_closed, contour_is_closed_within, Body, CONTOUR_CLOSURE_TOLERANCE};
pub use configuration::{Configuration, Element};
pub use error::GeometryError;
pub use layout::{ElementSpan, Layout};
pub use panel::Panel;
pub use placement::Placement;
pub use point::{point, vec2, Point, Vec2};
pub use spline::{CubicSpline, PanelingParams};
pub use xfoil_spline::XfoilSpline;
