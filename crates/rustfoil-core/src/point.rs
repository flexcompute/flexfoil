//! Point and vector types for 2D aerodynamic geometry.
//!
//! This module provides type aliases and convenience functions for working
//! with 2D points and vectors in the airfoil coordinate system.
//!
//! # Coordinate Convention
//! - **X-axis:** Increases in the downstream (freestream) direction.
//!   The leading edge is typically at x ≈ 0, trailing edge at x ≈ 1.
//! - **Y-axis:** Increases upward. Upper surface has positive y, lower surface negative.
//! - **Origin:** Typically at the leading edge for normalized coordinates.

use nalgebra::{Point2, Vector2};

/// A 2D point in the airfoil coordinate system.
///
/// This is a type alias for `nalgebra::Point2<f64>`, which provides:
/// - SIMD-optimized operations (when available)
/// - Proper affine point semantics (point - point = vector)
/// - Integration with nalgebra's matrix operations
pub type Point = Point2<f64>;

/// A 2D vector in the airfoil coordinate system.
///
/// Used for directions, normals, tangents, and velocities.
pub type Vec2 = Vector2<f64>;

/// Convenience constructor for creating a 2D point.
///
/// # Example
/// ```
/// use rustfoil_core::point::point;
///
/// let leading_edge = point(0.0, 0.0);
/// let trailing_edge = point(1.0, 0.0);
/// ```
#[inline]
pub fn point(x: f64, y: f64) -> Point {
    Point2::new(x, y)
}

/// Convenience constructor for creating a 2D vector.
///
/// # Example
/// ```
/// use rustfoil_core::point::vec2;
///
/// let freestream = vec2(1.0, 0.0);  // Unit velocity in x-direction
/// ```
#[inline]
pub fn vec2(x: f64, y: f64) -> Vec2 {
    Vector2::new(x, y)
}

/// Tolerance for geometric comparisons.
///
/// Two points closer than this distance are considered coincident.
/// This value is chosen to be:
/// - Small enough for accurate panel methods (panels ~ 0.01 chord)
/// - Large enough to avoid floating-point comparison issues
pub const GEOMETRY_TOLERANCE: f64 = 1e-12;

/// Check if two points are approximately equal within geometric tolerance.
#[inline]
pub fn points_coincident(a: &Point, b: &Point) -> bool {
    (a - b).norm_squared() < GEOMETRY_TOLERANCE * GEOMETRY_TOLERANCE
}

/// Compute the signed area of the triangle formed by three points.
///
/// Returns positive if the points are in counter-clockwise order,
/// negative if clockwise, and zero if collinear.
///
/// # Mathematical Background
/// This is the 2D cross product (z-component of 3D cross product):
/// ```text
/// signed_area = 0.5 * ((p2 - p1) × (p3 - p1))
/// ```
#[inline]
pub fn signed_triangle_area(p1: &Point, p2: &Point, p3: &Point) -> f64 {
    let v1 = p2 - p1;
    let v2 = p3 - p1;
    0.5 * (v1.x * v2.y - v1.y * v2.x)
}

/// Compute the 2D cross product (perpendicular dot product).
///
/// Returns the z-component of the 3D cross product of two 2D vectors
/// embedded in the xy-plane. Useful for determining handedness.
#[inline]
pub fn cross_2d(a: &Vec2, b: &Vec2) -> f64 {
    a.x * b.y - a.y * b.x
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_point_construction() {
        let p = point(1.5, -0.3);
        assert_relative_eq!(p.x, 1.5);
        assert_relative_eq!(p.y, -0.3);
    }

    #[test]
    fn test_points_coincident() {
        let p1 = point(1.0, 2.0);
        let p2 = point(1.0 + 1e-14, 2.0 - 1e-14);
        let p3 = point(1.1, 2.0);

        assert!(points_coincident(&p1, &p2));
        assert!(!points_coincident(&p1, &p3));
    }

    #[test]
    fn test_signed_triangle_area() {
        // Counter-clockwise triangle
        let p1 = point(0.0, 0.0);
        let p2 = point(1.0, 0.0);
        let p3 = point(0.0, 1.0);
        assert_relative_eq!(signed_triangle_area(&p1, &p2, &p3), 0.5);

        // Clockwise (reversed)
        assert_relative_eq!(signed_triangle_area(&p1, &p3, &p2), -0.5);

        // Collinear
        let p4 = point(0.5, 0.0);
        assert_relative_eq!(signed_triangle_area(&p1, &p2, &p4), 0.0);
    }

    #[test]
    fn test_cross_2d() {
        let a = vec2(1.0, 0.0);
        let b = vec2(0.0, 1.0);
        assert_relative_eq!(cross_2d(&a, &b), 1.0);
        assert_relative_eq!(cross_2d(&b, &a), -1.0);
    }
}
