//! Cubic spline interpolation for airfoil geometry.
//!
//! This module provides spline interpolation for smoothing and repaneling
//! airfoil coordinates. The primary use cases are:
//!
//! 1. **Smoothing**: Raw digitized airfoil data often has noise. Spline
//!    interpolation produces smooth derivatives needed for accurate
//!    boundary layer calculations.
//!
//! 2. **Repaneling**: XFOIL's "cosine spacing" clusters panels near the
//!    leading and trailing edges where gradients are steep. Splines allow
//!    evaluating the airfoil shape at arbitrary arc-length positions.
//!
//! # Algorithm
//! We use natural cubic splines with arc-length parameterization:
//! - Parameter `s` = cumulative arc length along the airfoil
//! - `x(s)` and `y(s)` are each represented by cubic splines
//!
//! # Future Work
//! - Implement parametric cubic spline interpolation
//! - Add cosine spacing repaneling
//! - Support tension splines for sharp features (optional)

use crate::error::GeometryError;
use crate::point::Point;

/// A parametric cubic spline representing an airfoil contour.
///
/// The spline is parameterized by arc length `s ∈ [0, s_max]`, where
/// `s_max` is the total arc length of the original point set.
///
/// # Construction
/// Use `CubicSpline::from_points()` to create a spline from discrete points.
/// The spline will pass exactly through all input points.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    /// Arc-length parameter values at each input point
    s_values: Vec<f64>,
    /// Spline coefficients for x(s)
    x_coeffs: SplineCoeffs,
    /// Spline coefficients for y(s)
    y_coeffs: SplineCoeffs,
}

/// Coefficients for a 1D cubic spline.
///
/// For segment `i`, the spline is:
/// ```text
/// f(t) = a[i] + b[i]*t + c[i]*t² + d[i]*t³
/// ```
/// where `t = s - s[i]` is the local parameter within the segment.
#[derive(Debug, Clone)]
struct SplineCoeffs {
    a: Vec<f64>, // Function values at knots
    b: Vec<f64>, // First derivative coefficients
    c: Vec<f64>, // Second derivative coefficients  
    d: Vec<f64>, // Third derivative coefficients
}

impl CubicSpline {
    /// Construct a cubic spline from a sequence of points.
    ///
    /// # Arguments
    /// * `points` - Ordered points along the airfoil contour
    ///
    /// # Errors
    /// - `InsufficientPoints` if fewer than 2 points provided
    /// - `SplineInterpolationFailed` if points are degenerate
    ///
    /// # Boundary Conditions
    /// Uses "natural" boundary conditions (zero second derivative at endpoints).
    /// This is appropriate for closed airfoil contours where the actual
    /// boundary condition is periodic.
    pub fn from_points(points: &[Point]) -> Result<Self, GeometryError> {
        if points.len() < 2 {
            return Err(GeometryError::InsufficientPoints {
                required: 2,
                provided: points.len(),
            });
        }

        // Compute arc-length parameterization
        let s_values = compute_arc_lengths(points);

        // Extract x and y coordinates
        let x_vals: Vec<f64> = points.iter().map(|p| p.x).collect();
        let y_vals: Vec<f64> = points.iter().map(|p| p.y).collect();

        // Build spline coefficients for x(s) and y(s)
        let x_coeffs = build_natural_spline(&s_values, &x_vals)?;
        let y_coeffs = build_natural_spline(&s_values, &y_vals)?;

        Ok(Self {
            s_values,
            x_coeffs,
            y_coeffs,
        })
    }

    /// Total arc length of the spline.
    pub fn total_arc_length(&self) -> f64 {
        *self.s_values.last().unwrap_or(&0.0)
    }

    /// Evaluate the spline at arc-length parameter `s`.
    ///
    /// # Arguments
    /// * `s` - Arc-length parameter in `[0, total_arc_length()]`
    ///
    /// # Returns
    /// The interpolated point `(x(s), y(s))`.
    ///
    /// # Panics
    /// May panic if `s` is significantly outside the valid range.
    pub fn evaluate(&self, s: f64) -> Point {
        let x = evaluate_spline(&self.s_values, &self.x_coeffs, s);
        let y = evaluate_spline(&self.s_values, &self.y_coeffs, s);
        Point::new(x, y)
    }

    /// Evaluate the first derivative (tangent) at arc-length `s`.
    ///
    /// Returns `(dx/ds, dy/ds)`, which is the unit tangent vector
    /// (since `s` is arc-length parameterized).
    pub fn derivative(&self, s: f64) -> (f64, f64) {
        let dx = evaluate_spline_derivative(&self.s_values, &self.x_coeffs, s);
        let dy = evaluate_spline_derivative(&self.s_values, &self.y_coeffs, s);
        (dx, dy)
    }

    /// Evaluate the second derivative (curvature-related) at arc-length `s`.
    pub fn second_derivative(&self, s: f64) -> (f64, f64) {
        let d2x = evaluate_spline_second_derivative(&self.s_values, &self.x_coeffs, s);
        let d2y = evaluate_spline_second_derivative(&self.s_values, &self.y_coeffs, s);
        (d2x, d2y)
    }

    /// Compute curvature at arc-length `s`.
    ///
    /// # Mathematical Background
    /// For a parametric curve `(x(s), y(s))`:
    /// ```text
    /// κ = (x'*y'' - y'*x'') / (x'² + y'²)^(3/2)
    /// ```
    /// Since `s` is arc-length, `x'² + y'² = 1`, so:
    /// ```text
    /// κ = x'*y'' - y'*x''
    /// ```
    pub fn curvature(&self, s: f64) -> f64 {
        let (dx, dy) = self.derivative(s);
        let (d2x, d2y) = self.second_derivative(s);

        // For arc-length parameterization, denominator ≈ 1
        let denom = (dx * dx + dy * dy).powf(1.5);
        if denom < 1e-12 {
            return 0.0;
        }

        (dx * d2y - dy * d2x) / denom
    }

    /// Resample the spline at `n` uniformly-spaced arc-length positions.
    ///
    /// This is useful for creating evenly-spaced panels (though cosine
    /// spacing is usually preferred for airfoils).
    pub fn resample_uniform(&self, n: usize) -> Vec<Point> {
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![self.evaluate(0.0)];
        }

        let s_max = self.total_arc_length();
        let ds = s_max / (n - 1) as f64;

        (0..n).map(|i| self.evaluate(i as f64 * ds)).collect()
    }

    /// Resample using cosine spacing (clusters points at LE and TE).
    ///
    /// # Cosine Spacing
    /// Instead of uniform `s` values, we use:
    /// ```text
    /// s_i = s_max * (1 - cos(π * i / (n-1))) / 2
    /// ```
    /// This produces dense spacing near `s = 0` and `s = s_max` (the trailing
    /// edge) and at the leading edge (if properly set up).
    ///
    /// # Arguments
    /// * `n` - Number of output points
    pub fn resample_cosine(&self, n: usize) -> Vec<Point> {
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![self.evaluate(0.0)];
        }

        let s_max = self.total_arc_length();
        let pi = std::f64::consts::PI;

        (0..n)
            .map(|i| {
                let theta = pi * i as f64 / (n - 1) as f64;
                let s = s_max * (1.0 - theta.cos()) / 2.0;
                self.evaluate(s)
            })
            .collect()
    }
}

/// Compute cumulative arc lengths for a sequence of points.
fn compute_arc_lengths(points: &[Point]) -> Vec<f64> {
    let mut s = Vec::with_capacity(points.len());
    s.push(0.0);

    for i in 1..points.len() {
        let ds = (points[i] - points[i - 1]).norm();
        s.push(s[i - 1] + ds);
    }

    s
}

/// Build natural cubic spline coefficients using the Thomas algorithm.
///
/// Natural spline: second derivative is zero at both endpoints.
fn build_natural_spline(s: &[f64], f: &[f64]) -> Result<SplineCoeffs, GeometryError> {
    let n = s.len();
    if n < 2 {
        return Err(GeometryError::SplineInterpolationFailed {
            reason: "need at least 2 points",
        });
    }

    if n == 2 {
        // Linear interpolation
        let h = s[1] - s[0];
        if h.abs() < 1e-14 {
            return Err(GeometryError::SplineInterpolationFailed {
                reason: "coincident points",
            });
        }
        return Ok(SplineCoeffs {
            a: vec![f[0]],
            b: vec![(f[1] - f[0]) / h],
            c: vec![0.0],
            d: vec![0.0],
        });
    }

    // Number of segments
    let n_seg = n - 1;

    // Interval sizes
    let h: Vec<f64> = (0..n_seg).map(|i| s[i + 1] - s[i]).collect();

    // Check for zero intervals
    for &hi in h.iter() {
        if hi.abs() < 1e-14 {
            return Err(GeometryError::SplineInterpolationFailed {
                reason: "coincident knots in spline",
            });
        }
    }

    // Build tridiagonal system for second derivatives (c values)
    // Natural BCs: c[0] = 0, c[n-1] = 0
    let mut c = vec![0.0; n];

    if n > 2 {
        // Interior equations
        let n_interior = n - 2;
        let mut lower = vec![0.0; n_interior];
        let mut diag = vec![0.0; n_interior];
        let mut upper = vec![0.0; n_interior];
        let mut rhs = vec![0.0; n_interior];

        for i in 0..n_interior {
            let j = i + 1; // Index in original arrays
            lower[i] = h[j - 1];
            diag[i] = 2.0 * (h[j - 1] + h[j]);
            upper[i] = h[j];
            rhs[i] = 3.0 * ((f[j + 1] - f[j]) / h[j] - (f[j] - f[j - 1]) / h[j - 1]);
        }

        // Solve tridiagonal system (Thomas algorithm)
        let c_interior = solve_tridiagonal(&lower, &diag, &upper, &rhs);

        for (i, &ci) in c_interior.iter().enumerate() {
            c[i + 1] = ci;
        }
    }

    // Compute remaining coefficients
    let a: Vec<f64> = f[..n_seg].to_vec();

    let b: Vec<f64> = (0..n_seg)
        .map(|i| (f[i + 1] - f[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0)
        .collect();

    let d: Vec<f64> = (0..n_seg)
        .map(|i| (c[i + 1] - c[i]) / (3.0 * h[i]))
        .collect();

    let c_seg: Vec<f64> = c[..n_seg].to_vec();

    Ok(SplineCoeffs {
        a,
        b,
        c: c_seg,
        d,
    })
}

/// Solve a tridiagonal system using the Thomas algorithm.
fn solve_tridiagonal(lower: &[f64], diag: &[f64], upper: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    if n == 0 {
        return vec![];
    }

    // Forward elimination
    let mut c_prime = vec![0.0; n];
    let mut d_prime = vec![0.0; n];

    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];

    for i in 1..n {
        let denom = diag[i] - lower[i] * c_prime[i - 1];
        c_prime[i] = if i < n - 1 {
            upper[i] / denom
        } else {
            0.0
        };
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    let mut x = vec![0.0; n];
    x[n - 1] = d_prime[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    x
}

/// Find the segment index for parameter value `s`.
fn find_segment(s_values: &[f64], s: f64) -> usize {
    // Binary search for the segment containing s
    match s_values.binary_search_by(|&si| si.partial_cmp(&s).unwrap()) {
        Ok(i) => i.min(s_values.len() - 2),
        Err(i) => (i.saturating_sub(1)).min(s_values.len() - 2),
    }
}

/// Evaluate spline at parameter `s`.
fn evaluate_spline(s_values: &[f64], coeffs: &SplineCoeffs, s: f64) -> f64 {
    let i = find_segment(s_values, s);
    let t = s - s_values[i];
    coeffs.a[i] + coeffs.b[i] * t + coeffs.c[i] * t * t + coeffs.d[i] * t * t * t
}

/// Evaluate first derivative at parameter `s`.
fn evaluate_spline_derivative(s_values: &[f64], coeffs: &SplineCoeffs, s: f64) -> f64 {
    let i = find_segment(s_values, s);
    let t = s - s_values[i];
    coeffs.b[i] + 2.0 * coeffs.c[i] * t + 3.0 * coeffs.d[i] * t * t
}

/// Evaluate second derivative at parameter `s`.
fn evaluate_spline_second_derivative(s_values: &[f64], coeffs: &SplineCoeffs, s: f64) -> f64 {
    let i = find_segment(s_values, s);
    let t = s - s_values[i];
    2.0 * coeffs.c[i] + 6.0 * coeffs.d[i] * t
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::point;
    use approx::assert_relative_eq;

    #[test]
    fn test_linear_interpolation() {
        let points = vec![point(0.0, 0.0), point(1.0, 1.0)];
        let spline = CubicSpline::from_points(&points).unwrap();

        let mid = spline.evaluate(spline.total_arc_length() / 2.0);
        assert_relative_eq!(mid.x, 0.5, epsilon = 0.01);
        assert_relative_eq!(mid.y, 0.5, epsilon = 0.01);
    }

    #[test]
    fn test_spline_passes_through_points() {
        let points = vec![
            point(0.0, 0.0),
            point(1.0, 1.0),
            point(2.0, 0.0),
            point(3.0, 1.0),
        ];
        let spline = CubicSpline::from_points(&points).unwrap();

        // Should pass through first and last points
        let p0 = spline.evaluate(0.0);
        assert_relative_eq!(p0.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(p0.y, 0.0, epsilon = 1e-10);

        let p_end = spline.evaluate(spline.total_arc_length());
        assert_relative_eq!(p_end.x, 3.0, epsilon = 1e-10);
        assert_relative_eq!(p_end.y, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cosine_spacing() {
        let points = vec![point(0.0, 0.0), point(1.0, 0.0), point(2.0, 0.0)];
        let spline = CubicSpline::from_points(&points).unwrap();

        let resampled = spline.resample_cosine(5);
        assert_eq!(resampled.len(), 5);

        // First and last should match original endpoints
        assert_relative_eq!(resampled[0].x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(resampled[4].x, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_insufficient_points() {
        let points = vec![point(0.0, 0.0)];
        let result = CubicSpline::from_points(&points);
        assert!(matches!(
            result,
            Err(GeometryError::InsufficientPoints { .. })
        ));
    }
}
