//! Panel representation for the panel method.
//!
//! A panel is the fundamental discretization unit in panel methods. Each panel
//! represents a linear segment of the airfoil surface with associated geometric
//! properties (normal, tangent, midpoint) that are used in the influence
//! coefficient calculations.
//!
//! # Aerodynamic Convention
//! - Panels are ordered **counter-clockwise** around the body.
//! - Starting from the trailing edge lower surface, going around the leading edge,
//!   and back to the trailing edge upper surface.
//! - The **normal vector** points **outward** from the body (into the flow).
//! - The **tangent vector** points in the direction of increasing panel index (CCW).

use crate::error::GeometryError;
use crate::point::{cross_2d, Point, Vec2, GEOMETRY_TOLERANCE};

/// A single panel representing a linear segment of the airfoil surface.
///
/// # Cached Properties
/// Geometric properties (midpoint, normal, tangent, length) are computed once
/// at construction time and cached. This avoids redundant calculations in
/// the hot loops of the influence coefficient matrix assembly.
///
/// # Panel Method Context
/// In the Linear Vorticity Panel Method, each panel carries:
/// - A linearly-varying vorticity distribution (γ varies from γ₁ to γ₂)
/// - A control point at the midpoint where boundary conditions are enforced
#[derive(Debug, Clone)]
pub struct Panel {
    /// Start point of the panel (in CCW ordering)
    pub p1: Point,
    /// End point of the panel
    pub p2: Point,

    // ---- Cached geometric properties (computed at construction) ----
    /// Midpoint of the panel (control point for boundary conditions)
    midpoint: Point,
    /// Outward-pointing unit normal vector
    normal: Vec2,
    /// Unit tangent vector (p1 → p2 direction)
    tangent: Vec2,
    /// Panel length
    length: f64,
}

impl Panel {
    /// Construct a new panel from two endpoints.
    ///
    /// # Arguments
    /// * `p1` - Start point (in counter-clockwise panel ordering)
    /// * `p2` - End point
    ///
    /// # Errors
    /// Returns `GeometryError::DegeneratePanel` if the panel length is below
    /// the geometric tolerance (effectively zero-length).
    ///
    /// # Normal Direction Convention
    /// The normal is computed by rotating the tangent 90° counter-clockwise:
    /// ```text
    /// tangent = (p2 - p1) / |p2 - p1|
    /// normal  = (-tangent.y, tangent.x)
    /// ```
    /// For CCW-ordered panels, this points outward from the body.
    pub fn new(p1: Point, p2: Point) -> Result<Self, GeometryError> {
        let delta = p2 - p1;
        let length = delta.norm();

        if length < GEOMETRY_TOLERANCE {
            return Err(GeometryError::DegeneratePanel { index: 0 });
        }

        // Unit tangent: direction from p1 to p2
        let tangent = delta / length;

        // Unit normal: rotate tangent 90° CCW (points outward for CCW panels)
        let normal = Vec2::new(-tangent.y, tangent.x);

        // Midpoint: control point for boundary condition
        let midpoint = Point::from((p1.coords + p2.coords) * 0.5);

        Ok(Self {
            p1,
            p2,
            midpoint,
            normal,
            tangent,
            length,
        })
    }

    /// Create a panel with an explicit index for error reporting.
    ///
    /// This is used when constructing panels from a point array, where the
    /// panel index is meaningful for debugging.
    pub fn new_with_index(p1: Point, p2: Point, index: usize) -> Result<Self, GeometryError> {
        Self::new(p1, p2).map_err(|_| GeometryError::DegeneratePanel { index })
    }

    /// Panel midpoint (control point for panel method).
    ///
    /// In the panel method, boundary conditions (no flow through surface)
    /// are enforced at panel midpoints, not endpoints. This gives better
    /// accuracy than endpoint collocation.
    #[inline]
    pub fn midpoint(&self) -> Point {
        self.midpoint
    }

    /// Outward-pointing unit normal vector.
    ///
    /// # Panel Method Usage
    /// The normal is used to compute the "normal velocity" component:
    /// ```text
    /// V_n = V · n = 0  (no-penetration condition)
    /// ```
    #[inline]
    pub fn normal(&self) -> Vec2 {
        self.normal
    }

    /// Unit tangent vector (p1 → p2 direction).
    ///
    /// # Panel Method Usage
    /// The tangent is used to compute surface velocity:
    /// ```text
    /// V_t = V · t  (tangential velocity for Cp calculation)
    /// ```
    #[inline]
    pub fn tangent(&self) -> Vec2 {
        self.tangent
    }

    /// Panel length.
    ///
    /// Used for:
    /// - Influence coefficient integration
    /// - Surface arc length calculations
    /// - Boundary layer marching (ds = panel length)
    #[inline]
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Panel angle relative to horizontal (in radians).
    ///
    /// Measured counter-clockwise from the positive x-axis.
    /// Range: (-π, π]
    #[inline]
    pub fn angle(&self) -> f64 {
        self.tangent.y.atan2(self.tangent.x)
    }

    /// Compute the signed perpendicular distance from a point to the panel line.
    ///
    /// Positive distance means the point is on the normal side (outside the body).
    /// Negative means inside.
    ///
    /// # Mathematical Background
    /// ```text
    /// distance = (point - midpoint) · normal
    /// ```
    pub fn signed_distance(&self, point: &Point) -> f64 {
        (point - self.midpoint).dot(&self.normal)
    }

    /// Project a point onto the panel line.
    ///
    /// Returns the parameter `t` in [0, 1] where the projection lies:
    /// - t = 0: projection is at p1
    /// - t = 1: projection is at p2
    /// - t < 0 or t > 1: projection is outside the panel segment
    pub fn project_parameter(&self, point: &Point) -> f64 {
        let v = point - self.p1;
        v.dot(&self.tangent) / self.length
    }

    /// Evaluate a point on the panel at parameter t ∈ [0, 1].
    ///
    /// ```text
    /// P(t) = p1 + t * (p2 - p1)
    /// ```
    pub fn point_at(&self, t: f64) -> Point {
        Point::from(self.p1.coords + t * (self.p2 - self.p1))
    }
}

/// Compute the turning angle between two adjacent panels.
///
/// This is the exterior angle at the shared vertex. Positive values indicate
/// the surface curves "outward" (convex), negative indicates "inward" (concave).
///
/// # Panel Method Usage
/// Large turning angles indicate high curvature regions where:
/// - More panels may be needed for accuracy
/// - Flow acceleration is significant
pub fn turning_angle(panel1: &Panel, panel2: &Panel) -> f64 {
    let cross = cross_2d(&panel1.tangent, &panel2.tangent);
    let dot = panel1.tangent.dot(&panel2.tangent);
    cross.atan2(dot)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::point;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_horizontal_panel() {
        let p = Panel::new(point(0.0, 0.0), point(1.0, 0.0)).unwrap();

        assert_relative_eq!(p.length(), 1.0);
        assert_relative_eq!(p.midpoint().x, 0.5);
        assert_relative_eq!(p.midpoint().y, 0.0);
        assert_relative_eq!(p.tangent().x, 1.0);
        assert_relative_eq!(p.tangent().y, 0.0);
        // Normal points up (CCW rotation of tangent)
        assert_relative_eq!(p.normal().x, 0.0);
        assert_relative_eq!(p.normal().y, 1.0);
        assert_relative_eq!(p.angle(), 0.0);
    }

    #[test]
    fn test_vertical_panel() {
        let p = Panel::new(point(0.0, 0.0), point(0.0, 1.0)).unwrap();

        assert_relative_eq!(p.tangent().x, 0.0);
        assert_relative_eq!(p.tangent().y, 1.0);
        // Normal points left (CCW rotation)
        assert_relative_eq!(p.normal().x, -1.0);
        assert_relative_eq!(p.normal().y, 0.0);
        assert_relative_eq!(p.angle(), PI / 2.0);
    }

    #[test]
    fn test_degenerate_panel() {
        let result = Panel::new(point(1.0, 2.0), point(1.0, 2.0));
        assert!(matches!(result, Err(GeometryError::DegeneratePanel { .. })));
    }

    #[test]
    fn test_signed_distance() {
        let p = Panel::new(point(0.0, 0.0), point(1.0, 0.0)).unwrap();

        // Point above the panel (on normal side)
        assert_relative_eq!(p.signed_distance(&point(0.5, 1.0)), 1.0);
        // Point below the panel
        assert_relative_eq!(p.signed_distance(&point(0.5, -1.0)), -1.0);
        // Point on the panel
        assert_relative_eq!(p.signed_distance(&point(0.5, 0.0)), 0.0);
    }

    #[test]
    fn test_project_parameter() {
        let p = Panel::new(point(0.0, 0.0), point(2.0, 0.0)).unwrap();

        assert_relative_eq!(p.project_parameter(&point(0.0, 5.0)), 0.0);
        assert_relative_eq!(p.project_parameter(&point(1.0, 3.0)), 0.5);
        assert_relative_eq!(p.project_parameter(&point(2.0, -1.0)), 1.0);
        // Outside the panel
        assert_relative_eq!(p.project_parameter(&point(-1.0, 0.0)), -0.5);
    }

    #[test]
    fn test_turning_angle() {
        // Two panels forming a 90° corner
        let p1 = Panel::new(point(0.0, 0.0), point(1.0, 0.0)).unwrap();
        let p2 = Panel::new(point(1.0, 0.0), point(1.0, 1.0)).unwrap();

        assert_relative_eq!(turning_angle(&p1, &p2), PI / 2.0, epsilon = 1e-10);
    }
}
