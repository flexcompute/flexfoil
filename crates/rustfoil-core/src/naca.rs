//! NACA airfoil generators.
//!
//! This module provides airfoil coordinate generation that matches XFOIL exactly.
//!
//! Reference: XFOIL source `naca.f` (Mark Drela, 2000)

use crate::point::{point, Point};

/// Default number of points per side (XFOIL's IQX/3 = 370/3 = 123).
pub const XFOIL_NSIDE: usize = 123;

/// Generate a NACA 4-digit airfoil using XFOIL's exact algorithm.
///
/// # Arguments
/// * `designation` - 4-digit NACA designation (e.g., 12 for NACA 0012, 2412 for NACA 2412)
/// * `nside` - Number of points per side (default: [`XFOIL_NSIDE`] = 123)
///
/// # Returns
/// Coordinates in XFOIL order: upper surface (TE→LE), then lower surface (LE→TE, no duplicate LE).
/// Total points = 2*nside - 1.
///
/// # Reference
/// XFOIL `naca.f` SUBROUTINE NACA4 (lines 21-82)
///
/// # Example
/// ```rust
/// use rustfoil_core::naca::naca4;
///
/// // Generate NACA 0012 with default resolution
/// let coords = naca4(12, None);
/// assert_eq!(coords.len(), 2 * 123 - 1);  // 245 points
///
/// // Generate NACA 2412 with custom resolution
/// let coords = naca4(2412, Some(50));
/// assert_eq!(coords.len(), 2 * 50 - 1);  // 99 points
/// ```
pub fn naca4(designation: u32, nside: Option<usize>) -> Vec<Point> {
    let nside = nside.unwrap_or(XFOIL_NSIDE);

    // Parse 4-digit designation: e.g., 2412 -> m=0.02, p=0.4, t=0.12
    let n4 = designation / 1000;
    let n3 = (designation - n4 * 1000) / 100;
    let n2 = (designation - n4 * 1000 - n3 * 100) / 10;
    let n1 = designation - n4 * 1000 - n3 * 100 - n2 * 10;

    let m = n4 as f64 / 100.0; // max camber
    let p = n3 as f64 / 10.0; // camber position
    let t = (n2 * 10 + n1) as f64 / 100.0; // thickness

    // TE point bunching parameter (XFOIL default)
    let an = 1.5_f64;
    let anp = an + 1.0;

    // Generate x and thickness/camber distributions
    let mut xx = vec![0.0; nside];
    let mut yt = vec![0.0; nside];
    let mut yc = vec![0.0; nside];

    for i in 0..nside {
        let frac = i as f64 / (nside - 1) as f64;

        // XFOIL's x-distribution with TE bunching (line 48)
        if i == nside - 1 {
            xx[i] = 1.0;
        } else {
            xx[i] = 1.0 - anp * frac * (1.0 - frac).powf(an) - (1.0 - frac).powf(anp);
        }

        // NACA 4-digit thickness formula (lines 50-54)
        // Note: XFOIL uses the original formula (blunt TE), not the modified one
        yt[i] = (0.29690 * xx[i].sqrt()
            - 0.12600 * xx[i]
            - 0.35160 * xx[i].powi(2)
            + 0.28430 * xx[i].powi(3)
            - 0.10150 * xx[i].powi(4))
            * t
            / 0.20;

        // Camber line (lines 55-59)
        if p > 0.0 {
            if xx[i] < p {
                yc[i] = m / (p * p) * (2.0 * p * xx[i] - xx[i] * xx[i]);
            } else {
                yc[i] = m / ((1.0 - p) * (1.0 - p))
                    * ((1.0 - 2.0 * p) + 2.0 * p * xx[i] - xx[i] * xx[i]);
            }
        }
    }

    // Build output coordinates in XFOIL order
    let mut coords = Vec::with_capacity(2 * nside - 1);

    // Upper surface: from TE (i=nside-1) to LE (i=0) - reversed
    for i in (0..nside).rev() {
        coords.push(point(xx[i], yc[i] + yt[i]));
    }

    // Lower surface: from LE (i=1, skip 0) to TE (i=nside-1)
    for i in 1..nside {
        coords.push(point(xx[i], yc[i] - yt[i]));
    }

    coords
}

/// Check if coordinates are symmetric (for symmetric airfoils like NACA 00xx).
pub fn is_symmetric(coords: &[Point], tol: f64) -> bool {
    let n = coords.len();
    let nside = (n + 1) / 2;

    for i in 0..nside {
        let j = n - 1 - i;
        if i >= j {
            break;
        }

        let x_diff = (coords[i].x - coords[j].x).abs();
        let y_sum = (coords[i].y + coords[j].y).abs();

        if x_diff > tol || y_sum > tol {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naca0012() {
        let coords = naca4(12, None);

        // Should have 2*123 - 1 = 245 points
        assert_eq!(coords.len(), 2 * XFOIL_NSIDE - 1);

        // First and last points should be at TE
        assert!((coords[0].x - 1.0).abs() < 1e-10);
        assert!((coords.last().unwrap().x - 1.0).abs() < 1e-10);

        // LE should be at x=0
        let le_idx = (coords.len() + 1) / 2 - 1;
        assert!(coords[le_idx].x < 1e-10);

        // Should be symmetric
        assert!(is_symmetric(&coords, 1e-12));
    }

    #[test]
    fn test_naca2412() {
        let coords = naca4(2412, Some(100));

        // Should have 2*100 - 1 = 199 points
        assert_eq!(coords.len(), 199);

        // First and last points should be at TE
        assert!((coords[0].x - 1.0).abs() < 1e-10);

        // Should NOT be symmetric (cambered airfoil)
        assert!(!is_symmetric(&coords, 1e-6));

        // Upper surface should have more y than lower surface at max camber location
        // Max camber is at x = 0.4 for NACA 2412
    }

    #[test]
    fn test_te_thickness() {
        let coords = naca4(12, Some(100));

        // NACA 0012 TE thickness should be 2 * y_te
        // Original formula: y = 0.12 * (0.29690*sqrt(1) - 0.12600*1 - ...) = 0.12 * 0.0021 = 0.00252
        // But XFOIL formula is t/0.20 * (...), so y = 0.12/0.20 * 0.0021 = 0.00126
        let y_te_upper = coords[0].y;
        let y_te_lower = coords.last().unwrap().y;

        assert!((y_te_upper - 0.00126).abs() < 1e-5);
        assert!((y_te_lower + 0.00126).abs() < 1e-5);
    }
}
