//! Aerodynamic body representation.
//!
//! A `Body` represents a single aerodynamic element (airfoil, slat, flap, etc.)
//! as a closed contour discretized into panels.
//!
//! # Multi-Body Support
//! Unlike XFOIL's single-body limitation, RustFoil is designed from the ground
//! up to support multiple bodies. A typical multi-element configuration might be:
//! ```text
//! bodies: Vec<Body> = vec![slat, main_wing, flap]
//! ```
//!
//! # Panel Ordering Convention
//! Panels are ordered **counter-clockwise** starting from the trailing edge:
//! 1. Lower surface (TE → LE)
//! 2. Upper surface (LE → TE)
//!
//! This convention ensures:
//! - Normal vectors point outward (into the flow)
//! - The trailing edge is easily identified as the first/last point
//! - Boundary layer marching proceeds naturally from stagnation point

use crate::error::GeometryError;
use crate::panel::Panel;
use crate::point::{points_coincident, Point};

/// A single aerodynamic body discretized into panels.
///
/// # Invariants
/// - The body forms a closed contour (first and last points are coincident,
///   or close enough that they form the trailing edge)
/// - All panels have non-zero length
/// - Panels are ordered counter-clockwise
#[derive(Debug, Clone)]
pub struct Body {
    /// Human-readable identifier (e.g., "main", "slat", "flap")
    pub name: String,

    /// Ordered panels forming the closed contour
    panels: Vec<Panel>,

    /// Index of the trailing edge panel (upper surface side).
    ///
    /// For the Kutta condition, we need to know which panel is at the TE.
    /// This is typically the last panel (index = n_panels - 1), but may
    /// differ for blunt trailing edges or special geometries.
    te_panel_upper: usize,

    /// Index of the trailing edge panel (lower surface side).
    ///
    /// For sharp trailing edges, this is panel 0. For blunt TE, there may
    /// be a "base" panel connecting the two TE points.
    te_panel_lower: usize,

    /// Index of the leading edge point (approximate, for reference).
    ///
    /// The LE is identified as the point with minimum x-coordinate.
    le_point_idx: usize,
}

impl Body {
    /// Construct a body from raw coordinate points.
    ///
    /// # Arguments
    /// * `name` - Identifier for the body (e.g., "main", "slat")
    /// * `points` - Ordered points forming the airfoil contour. Should be
    ///   ordered CCW from the trailing edge lower surface, around the leading
    ///   edge, back to the trailing edge upper surface.
    ///
    /// # Point Closure
    /// The points should form a closed contour:
    /// - Either the first and last points are coincident (sharp TE)
    /// - Or they are distinct (blunt TE), in which case no closing panel is added
    ///
    /// # Errors
    /// - `InsufficientPoints` if fewer than 3 points are provided
    /// - `DegeneratePanel` if any resulting panel has zero length
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::body::Body;
    /// use rustfoil_core::point::point;
    ///
    /// // Simple diamond airfoil
    /// let points = vec![
    ///     point(1.0, 0.0),   // TE
    ///     point(0.5, -0.1),  // Lower
    ///     point(0.0, 0.0),   // LE
    ///     point(0.5, 0.1),   // Upper
    ///     point(1.0, 0.0),   // Back to TE (closed)
    /// ];
    /// let body = Body::from_points("diamond", &points).unwrap();
    /// assert_eq!(body.n_panels(), 4);
    /// ```
    pub fn from_points(name: &str, points: &[Point]) -> Result<Self, GeometryError> {
        const MIN_POINTS: usize = 3;

        if points.len() < MIN_POINTS {
            return Err(GeometryError::InsufficientPoints {
                required: MIN_POINTS,
                provided: points.len(),
            });
        }

        // Check if the contour is closed (first ≈ last point)
        let _is_closed = points_coincident(&points[0], &points[points.len() - 1]);

        // Number of panels:
        // - Closed contour: n_points - 1 (last point is duplicate)
        // - Open contour: n_points - 1 (no closing panel added)
        let n_panels = points.len() - 1;

        // Build panels
        let mut panels = Vec::with_capacity(n_panels);
        for i in 0..n_panels {
            let panel = Panel::new_with_index(points[i], points[i + 1], i)?;
            panels.push(panel);
        }

        // If contour is not closed and we want to close it, we'd add a panel here.
        // For now, we assume input is properly formatted (closed or intentionally open).

        // Find trailing edge panels (first and last by convention)
        let te_panel_lower = 0;
        let te_panel_upper = panels.len() - 1;

        // Find leading edge (minimum x-coordinate)
        let le_point_idx = points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);

        Ok(Self {
            name: name.to_string(),
            panels,
            te_panel_upper,
            te_panel_lower,
            le_point_idx,
        })
    }

    /// Number of panels in the body.
    #[inline]
    pub fn n_panels(&self) -> usize {
        self.panels.len()
    }

    /// Access the panels as a slice.
    #[inline]
    pub fn panels(&self) -> &[Panel] {
        &self.panels
    }

    /// Mutable access to panels (for repaneling operations).
    #[inline]
    pub fn panels_mut(&mut self) -> &mut [Panel] {
        &mut self.panels
    }

    /// Index of the upper-surface trailing edge panel.
    ///
    /// This is used for enforcing the Kutta condition, which requires
    /// the vorticity at the trailing edge to satisfy γ_upper + γ_lower = 0
    /// (for a sharp TE) to ensure finite velocity.
    #[inline]
    pub fn te_upper_index(&self) -> usize {
        self.te_panel_upper
    }

    /// Index of the lower-surface trailing edge panel.
    #[inline]
    pub fn te_lower_index(&self) -> usize {
        self.te_panel_lower
    }

    /// Index of the leading edge point.
    #[inline]
    pub fn le_index(&self) -> usize {
        self.le_point_idx
    }

    /// Total surface arc length of the body.
    pub fn arc_length(&self) -> f64 {
        self.panels.iter().map(|p| p.length()).sum()
    }

    /// Get all panel midpoints (control points) as a vector.
    ///
    /// Useful for setting up the influence coefficient matrix.
    pub fn control_points(&self) -> Vec<Point> {
        self.panels.iter().map(|p| p.midpoint()).collect()
    }

    /// Compute the approximate chord length.
    ///
    /// Defined as the distance from leading edge to trailing edge.
    pub fn chord(&self) -> f64 {
        if self.panels.is_empty() {
            return 0.0;
        }

        let le = &self.panels[self.le_point_idx.min(self.panels.len() - 1)].p1;
        let te = &self.panels[self.te_panel_lower].p1;

        (te - le).norm()
    }

    /// Iterator over panel indices for the lower surface (TE → LE).
    pub fn lower_surface_indices(&self) -> impl Iterator<Item = usize> {
        0..=self.le_point_idx.min(self.panels.len() - 1)
    }

    /// Iterator over panel indices for the upper surface (LE → TE).
    pub fn upper_surface_indices(&self) -> impl Iterator<Item = usize> {
        self.le_point_idx.min(self.panels.len() - 1)..self.panels.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::point;
    use approx::assert_relative_eq;

    fn make_diamond_airfoil() -> Body {
        // Simple closed diamond shape
        let points = vec![
            point(1.0, 0.0),  // TE
            point(0.5, -0.1), // Lower
            point(0.0, 0.0),  // LE
            point(0.5, 0.1),  // Upper
            point(1.0, 0.0),  // Back to TE
        ];
        Body::from_points("diamond", &points).unwrap()
    }

    #[test]
    fn test_body_construction() {
        let body = make_diamond_airfoil();

        assert_eq!(body.n_panels(), 4);
        assert_eq!(body.name, "diamond");
        assert_eq!(body.te_lower_index(), 0);
        assert_eq!(body.te_upper_index(), 3);
    }

    #[test]
    fn test_leading_edge_detection() {
        let body = make_diamond_airfoil();

        // LE should be at index 2 (the point at x=0)
        assert_eq!(body.le_index(), 2);
    }

    #[test]
    fn test_chord_length() {
        let body = make_diamond_airfoil();

        // Chord from LE (0,0) to TE (1,0) should be 1.0
        assert_relative_eq!(body.chord(), 1.0, epsilon = 0.1);
    }

    #[test]
    fn test_insufficient_points() {
        let points = vec![point(0.0, 0.0), point(1.0, 0.0)];
        let result = Body::from_points("test", &points);

        assert!(matches!(
            result,
            Err(GeometryError::InsufficientPoints {
                required: 3,
                provided: 2
            })
        ));
    }

    #[test]
    fn test_arc_length() {
        // Square: 4 panels of length 1
        let points = vec![
            point(1.0, 0.0),
            point(0.0, 0.0),
            point(0.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 0.0),
        ];
        let body = Body::from_points("square", &points).unwrap();

        assert_relative_eq!(body.arc_length(), 4.0);
    }

    #[test]
    fn test_control_points() {
        let body = make_diamond_airfoil();
        let cps = body.control_points();

        assert_eq!(cps.len(), 4);
        // First panel midpoint should be between (1,0) and (0.5,-0.1)
        assert_relative_eq!(cps[0].x, 0.75);
        assert_relative_eq!(cps[0].y, -0.05);
    }
}
