//! Float comparison utilities for testing numerical code.

/// Compute the relative error between two values.
///
/// Returns `|a - b| / max(|a|, |b|, 1e-10)` to handle near-zero values gracefully.
pub fn relative_error(a: f64, b: f64) -> f64 {
    let diff = (a - b).abs();
    let scale = a.abs().max(b.abs()).max(1e-10);
    diff / scale
}

/// Assert that two values are close within a relative tolerance.
///
/// # Panics
/// Panics with a descriptive message if the relative error exceeds `tol`.
pub fn assert_close(a: f64, b: f64, tol: f64, name: &str) {
    let rel_err = relative_error(a, b);
    if rel_err > tol {
        panic!(
            "{}: values differ by {:.2e} (rel), expected {:.8}, got {:.8}",
            name, rel_err, a, b
        );
    }
}

/// Assert that two values are close within an absolute tolerance.
///
/// # Panics
/// Panics with a descriptive message if `|a - b| > tol`.
pub fn assert_close_abs(a: f64, b: f64, tol: f64, name: &str) {
    let diff = (a - b).abs();
    if diff > tol {
        panic!(
            "{}: values differ by {:.2e} (abs), expected {:.8}, got {:.8}",
            name, diff, a, b
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_relative_error_identical() {
        assert_eq!(relative_error(1.0, 1.0), 0.0);
    }

    #[test]
    fn test_relative_error_small_diff() {
        let err = relative_error(1.0, 1.001);
        // err = |1.0 - 1.001| / max(1.0, 1.001) = 0.001 / 1.001 ≈ 0.000999
        let expected = 0.001 / 1.001;
        assert!((err - expected).abs() < 1e-10);
    }

    #[test]
    fn test_relative_error_near_zero() {
        // Should not divide by zero
        let err = relative_error(1e-15, 2e-15);
        assert!(err.is_finite());
    }

    #[test]
    fn test_assert_close_passes() {
        assert_close(1.0, 1.0001, 1e-3, "test");
    }

    #[test]
    #[should_panic(expected = "values differ")]
    fn test_assert_close_fails() {
        assert_close(1.0, 1.1, 1e-3, "test");
    }
}
