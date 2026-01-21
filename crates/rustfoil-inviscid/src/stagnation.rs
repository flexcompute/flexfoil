//! Stagnation point detection (XFOIL's STFIND).
//!
//! This module finds the stagnation point on the airfoil surface,
//! which is critical for boundary layer initialization.
//!
//! # Algorithm
//!
//! The stagnation point is where the surface velocity (gamma) crosses zero.
//! XFOIL uses linear interpolation to find the exact arc-length position.
//!
//! # XFOIL Reference
//!
//! - `xpanel.f`: STFIND subroutine
//! - `xbl.f`: IBLPAN, XICALC

use crate::geometry::AirfoilGeometry;

/// Stagnation point information.
#[derive(Debug, Clone, Copy)]
pub struct StagnationPoint {
    /// Panel index (stagnation is between IST and IST+1)
    pub ist: usize,
    /// Exact arc length at stagnation point
    pub sst: f64,
    /// X coordinate at stagnation
    pub xst: f64,
    /// Y coordinate at stagnation
    pub yst: f64,
}

/// Find the stagnation point on the airfoil surface.
///
/// The stagnation point is where the surface tangential velocity (gamma)
/// crosses from positive to negative when traversing counter-clockwise.
///
/// # Arguments
///
/// * `gamma` - Surface vorticity distribution (= tangential velocity)
/// * `geom` - Airfoil geometry
///
/// # Returns
///
/// `StagnationPoint` with the exact location, or `None` if not found.
///
/// # XFOIL Algorithm
///
/// ```fortran
/// DO I=1, N-1
///     IF(GAM(I).GE.0.0 .AND. GAM(I+1).LT.0.0) GO TO 11
/// ENDDO
/// 11 CONTINUE
/// IST = I
///
/// DGAM = GAM(I+1) - GAM(I)
/// DS = S(I+1) - S(I)
///
/// C---- minimize roundoff for very small GAM
/// IF(GAM(I) .LT. -GAM(I+1)) THEN
///     SST = S(I)   - DS*(GAM(I)  /DGAM)
/// ELSE
///     SST = S(I+1) - DS*(GAM(I+1)/DGAM)
/// ENDIF
/// ```
pub fn find_stagnation(gamma: &[f64], geom: &AirfoilGeometry) -> Option<StagnationPoint> {
    let n = gamma.len();
    
    if n < 2 || n != geom.n {
        return None;
    }

    // Find where gamma crosses from positive to negative
    let mut ist = None;
    for i in 0..n - 1 {
        if gamma[i] >= 0.0 && gamma[i + 1] < 0.0 {
            ist = Some(i);
            break;
        }
    }

    let ist = ist?;
    
    let dgam = gamma[ist + 1] - gamma[ist];
    let ds = geom.s[ist + 1] - geom.s[ist];

    // XFOIL's roundoff minimization logic
    let sst = if dgam.abs() < 1e-20 {
        // Gamma is essentially constant - use midpoint
        0.5 * (geom.s[ist] + geom.s[ist + 1])
    } else if gamma[ist] < -gamma[ist + 1] {
        // GAM(I) is smaller in magnitude - interpolate from start
        geom.s[ist] - ds * (gamma[ist] / dgam)
    } else {
        // GAM(I+1) is smaller in magnitude - interpolate from end
        geom.s[ist + 1] - ds * (gamma[ist + 1] / dgam)
    };

    // Clamp to panel range with small offset (XFOIL's tweak)
    let sst = if sst <= geom.s[ist] {
        geom.s[ist] + 1e-7
    } else if sst >= geom.s[ist + 1] {
        geom.s[ist + 1] - 1e-7
    } else {
        sst
    };

    // Interpolate position
    let (xst, yst) = interpolate_position(sst, geom);

    Some(StagnationPoint { ist, sst, xst, yst })
}

/// Interpolate position on airfoil surface at arc-length s.
fn interpolate_position(s: f64, geom: &AirfoilGeometry) -> (f64, f64) {
    let n = geom.n;
    
    // Find segment
    let mut i = 0;
    for j in 0..n - 1 {
        if s >= geom.s[j] && s <= geom.s[j + 1] {
            i = j;
            break;
        }
    }

    let ds = geom.s[i + 1] - geom.s[i];
    if ds < 1e-20 {
        return (geom.x[i], geom.y[i]);
    }

    let t = (s - geom.s[i]) / ds;
    let x = geom.x[i] + t * (geom.x[i + 1] - geom.x[i]);
    let y = geom.y[i] + t * (geom.y[i + 1] - geom.y[i]);

    (x, y)
}

/// Extract upper and lower surface indices from stagnation point.
///
/// This matches XFOIL's IBLPAN logic.
///
/// # Returns
///
/// (upper_indices, lower_indices) where:
/// - upper_indices: panels from stagnation toward upper TE (decreasing i)
/// - lower_indices: panels from stagnation toward lower TE (increasing i)
pub fn surface_panels(stag: &StagnationPoint, n: usize) -> (Vec<usize>, Vec<usize>) {
    let ist = stag.ist;

    // Upper surface: from stagnation toward node 0 (upper TE)
    // Panels: IST, IST-1, IST-2, ..., 0
    let upper: Vec<usize> = (0..=ist).rev().collect();

    // Lower surface: from stagnation toward node N-1 (lower TE)
    // Panels: IST+1, IST+2, ..., N-1
    let lower: Vec<usize> = (ist + 1..n).collect();

    (upper, lower)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_test_geometry() -> AirfoilGeometry {
        let n_panels = 40;
        let n_half = n_panels / 2;
        let t = 0.12;

        let x_coords: Vec<f64> = (0..=n_half)
            .map(|i| {
                let beta = PI * (i as f64) / (n_half as f64);
                0.5 * (1.0 - beta.cos())
            })
            .collect();

        let thickness = |x: f64| -> f64 {
            5.0 * t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) 
                + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
        };

        let mut points = Vec::with_capacity(2 * n_half);

        for i in (0..=n_half).rev() {
            let x = x_coords[i];
            let y = thickness(x);
            points.push((x, y));
        }

        for i in 1..=n_half {
            let x = x_coords[i];
            let y = -thickness(x);
            points.push((x, y));
        }

        AirfoilGeometry::from_points(&points).unwrap()
    }

    #[test]
    fn test_stagnation_at_zero_alpha() {
        let geom = make_test_geometry();
        let n = geom.n;

        // At α=0 for symmetric airfoil, gamma should cross zero near the LE
        // Create a simple gamma distribution: positive on upper, negative on lower
        let mut gamma = vec![0.0; n];
        let n_half = n / 2;
        
        for i in 0..n_half {
            gamma[i] = 1.0 - (i as f64) / (n_half as f64); // Decreasing on upper
        }
        for i in n_half..n {
            gamma[i] = -1.0 + ((i - n_half) as f64) / ((n - n_half) as f64); // Increasing on lower
        }

        let stag = find_stagnation(&gamma, &geom);
        assert!(stag.is_some(), "Should find stagnation point");

        let stag = stag.unwrap();
        assert!(stag.ist < n - 1, "IST should be valid");
        assert!(stag.sst > 0.0, "SST should be positive");
    }

    #[test]
    fn test_surface_panels() {
        let stag = StagnationPoint {
            ist: 20,
            sst: 0.5,
            xst: 0.0,
            yst: 0.0,
        };

        let (upper, lower) = surface_panels(&stag, 41);

        // Upper should be 20, 19, 18, ..., 0
        assert_eq!(upper.len(), 21);
        assert_eq!(upper[0], 20);
        assert_eq!(upper[20], 0);

        // Lower should be 21, 22, ..., 40
        assert_eq!(lower.len(), 20);
        assert_eq!(lower[0], 21);
        assert_eq!(lower[19], 40);
    }
}
