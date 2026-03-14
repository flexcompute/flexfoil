//! Airfoil geometry processing matching XFOIL's NCALC, APCALC, and TECALC.
//!
//! This module computes:
//! - Node positions and arc-length parameterization
//! - Spline derivatives (dX/dS, dY/dS)
//! - Outward unit normal vectors (NCALC)
//! - Panel angles (APCALC)
//! - Trailing edge geometry (TECALC)
//! - Leading edge location (LEFIND)
//!
//! # XFOIL Reference
//!
//! - `xpanel.f`: NCALC (lines 51-96), APCALC (lines 22-48)
//! - `xgeom.f`: LEFIND, TECALC
//! - `XFOIL.INC`: Common block definitions

use crate::{InviscidError, Result};
use std::f64::consts::PI;

/// Complete airfoil geometry with all derived quantities needed for panel method.
///
/// This structure matches XFOIL's internal representation in `XFOIL.INC`.
#[derive(Debug, Clone)]
pub struct AirfoilGeometry {
    // === Node arrays (indexed 0..n-1) ===
    
    /// X coordinates at each node
    pub x: Vec<f64>,
    /// Y coordinates at each node
    pub y: Vec<f64>,
    /// Arc length from first point to each node
    pub s: Vec<f64>,
    /// dX/dS at each node (spline derivative)
    pub xp: Vec<f64>,
    /// dY/dS at each node (spline derivative)
    pub yp: Vec<f64>,
    /// Outward unit normal X component at each node
    pub nx: Vec<f64>,
    /// Outward unit normal Y component at each node
    pub ny: Vec<f64>,
    /// Panel angle for each panel (atan2 formulation)
    pub apanel: Vec<f64>,
    
    // === Trailing edge geometry ===
    
    /// TE midpoint X coordinate
    pub xte: f64,
    /// TE midpoint Y coordinate  
    pub yte: f64,
    /// TE gap length (distance from node 0 to node n-1)
    pub dste: f64,
    /// TE normal-projected gap (ANTE in XFOIL)
    pub ante: f64,
    /// TE tangent-projected gap (ASTE in XFOIL)
    pub aste: f64,
    /// True if trailing edge is sharp (DSTE < 0.0001 * chord)
    pub sharp: bool,
    
    // === Leading edge geometry ===
    
    /// LE X coordinate
    pub xle: f64,
    /// LE Y coordinate
    pub yle: f64,
    /// Arc length at leading edge
    pub sle: f64,
    
    // === Dimensions ===
    
    /// Number of nodes
    pub n: usize,
    /// Chord length
    pub chord: f64,
}

impl AirfoilGeometry {
    /// Build geometry from airfoil coordinate points.
    ///
    /// # Arguments
    ///
    /// * `points` - (x, y) coordinates ordered counter-clockwise from upper trailing edge
    ///
    /// # XFOIL Convention
    ///
    /// - Node 0 is upper trailing edge
    /// - Node n-1 is lower trailing edge
    /// - For sharp TE, nodes 0 and n-1 are very close but not identical
    /// - Panel i goes from node i to node i+1
    /// - Panel n-1 (the "TE panel") goes from node n-1 to node 0
    pub fn from_points(points: &[(f64, f64)]) -> Result<Self> {
        let n = points.len();
        
        if n < 10 {
            return Err(InviscidError::InsufficientPoints(n));
        }

        // Extract x, y coordinates
        let x: Vec<f64> = points.iter().map(|p| p.0).collect();
        let y: Vec<f64> = points.iter().map(|p| p.1).collect();

        // Compute arc-length parameterization (SCALC in XFOIL)
        let s = Self::compute_arc_length(&x, &y);

        // Check for duplicate points
        for i in 0..n - 1 {
            if (s[i + 1] - s[i]).abs() < 1e-12 {
                return Err(InviscidError::DuplicatePoints(i));
            }
        }

        // Compute spline derivatives using XFOIL's SEGSPL (zero third derivative end conditions)
        let xp = Self::segspl(&s, &x);
        let yp = Self::segspl(&s, &y);

        // Compute normals (NCALC)
        let (nx, ny) = Self::compute_normals(&xp, &yp, &s);

        // Compute panel angles (APCALC)
        let apanel = Self::compute_panel_angles(&x, &y, &nx, &ny);

        // Compute trailing edge geometry (TECALC)
        let (xte, yte, dste, ante, aste) = Self::compute_te_geometry(&x, &y, &apanel);

        // Find leading edge (LEFIND)
        let (xle, yle, sle) = Self::find_leading_edge(&x, &y, &s, &xp, &yp, xte, yte);

        // XFOIL defines CHORD from the true LE/TE geometry, not the raw
        // coordinate extents. This matters for paneled files whose minimum x
        // is slightly ahead of the exact spline-leading-edge location.
        let chord = Self::compute_chord(xte, yte, xle, yle);

        if chord <= 0.0 {
            return Err(InviscidError::InvalidChord(chord));
        }

        // Determine if sharp TE
        let sharp = dste < 0.0001 * chord;

        Ok(Self {
            x,
            y,
            s,
            xp,
            yp,
            nx,
            ny,
            apanel,
            xte,
            yte,
            dste,
            ante,
            aste,
            sharp,
            xle,
            yle,
            sle,
            n,
            chord,
        })
    }

    /// Compute arc-length parameterization (XFOIL's SCALC).
    fn compute_arc_length(x: &[f64], y: &[f64]) -> Vec<f64> {
        let n = x.len();
        let mut s = vec![0.0; n];
        
        for i in 1..n {
            let dx = x[i] - x[i - 1];
            let dy = y[i] - y[i - 1];
            s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
        }
        
        s
    }

    /// XFOIL's SEGSPL - compute spline derivatives with zero third derivative end conditions.
    ///
    /// This matches XFOIL's `spline.f` SEGSPL subroutine exactly.
    fn segspl(s: &[f64], f: &[f64]) -> Vec<f64> {
        let n = s.len();
        if n < 2 {
            return vec![0.0; n];
        }
        if n == 2 {
            let df = (f[1] - f[0]) / (s[1] - s[0]);
            return vec![df, df];
        }

        let mut a = vec![0.0; n];  // Main diagonal
        let mut b = vec![0.0; n];  // Upper diagonal
        let mut c = vec![0.0; n];  // Lower diagonal
        let mut fs = vec![0.0; n]; // RHS, then solution

        // Interior points (XFOIL lines 87-94 in spline.f)
        for i in 1..n - 1 {
            let dsm = s[i] - s[i - 1];
            let dsp = s[i + 1] - s[i];
            b[i] = dsp;
            a[i] = 2.0 * (dsm + dsp);
            c[i] = dsm;
            fs[i] = 3.0 * ((f[i + 1] - f[i]) * dsm / dsp + (f[i] - f[i - 1]) * dsp / dsm);
        }

        // Zero THIRD derivative end conditions (XFOIL lines 101-105, 117-120)
        // This is SEGSPL with XS1 = XS2 = -999.0
        a[0] = 1.0;
        c[0] = 1.0;
        fs[0] = 2.0 * (f[1] - f[0]) / (s[1] - s[0]);

        b[n - 1] = 1.0;
        a[n - 1] = 1.0;
        fs[n - 1] = 2.0 * (f[n - 1] - f[n - 2]) / (s[n - 1] - s[n - 2]);

        // Solve tridiagonal system (TRISOL)
        Self::trisol(&a, &b, &c, &mut fs);

        fs
    }

    /// XFOIL's TRISOL - tridiagonal matrix solver.
    fn trisol(a: &[f64], b: &[f64], c: &[f64], x: &mut [f64]) {
        let n = a.len();
        if n < 2 {
            return;
        }

        let mut aa = a.to_vec();
        
        // Forward elimination
        for i in 1..n {
            let piv = b[i] / aa[i - 1];
            aa[i] = aa[i] - c[i - 1] * piv;
            x[i] = x[i] - x[i - 1] * piv;
        }

        // Back substitution
        x[n - 1] = x[n - 1] / aa[n - 1];
        for i in (0..n - 1).rev() {
            x[i] = (x[i] - c[i] * x[i + 1]) / aa[i];
        }
    }

    /// Compute outward unit normal vectors at each node (XFOIL's NCALC).
    ///
    /// The normal is computed by rotating the tangent 90° counter-clockwise:
    /// - tangent = (dX/dS, dY/dS)
    /// - normal = (dY/dS, -dX/dS) / |tangent|
    ///
    /// For counter-clockwise traversal, this gives outward-pointing normals.
    fn compute_normals(xp: &[f64], yp: &[f64], s: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = xp.len();
        let mut nx = vec![0.0; n];
        let mut ny = vec![0.0; n];

        for i in 0..n {
            // Rotate tangent 90° CCW to get outward normal
            let sx = yp[i];
            let sy = -xp[i];
            let smod = (sx * sx + sy * sy).sqrt();

            if smod < 1e-12 {
                // Degenerate case: use -x direction (XFOIL default)
                nx[i] = -1.0;
                ny[i] = 0.0;
            } else {
                nx[i] = sx / smod;
                ny[i] = sy / smod;
            }
        }

        // Average normal vectors at corner points (where s[i] == s[i+1])
        // This handles sharp corners in the geometry
        for i in 0..n - 1 {
            if (s[i] - s[i + 1]).abs() < 1e-12 {
                let sx = 0.5 * (nx[i] + nx[i + 1]);
                let sy = 0.5 * (ny[i] + ny[i + 1]);
                let smod = (sx * sx + sy * sy).sqrt();

                if smod < 1e-12 {
                    nx[i] = -1.0;
                    ny[i] = 0.0;
                    nx[i + 1] = -1.0;
                    ny[i + 1] = 0.0;
                } else {
                    nx[i] = sx / smod;
                    ny[i] = sy / smod;
                    nx[i + 1] = sx / smod;
                    ny[i + 1] = sy / smod;
                }
            }
        }

        (nx, ny)
    }

    /// Compute panel angles (XFOIL's APCALC).
    ///
    /// Panel angle is the angle of the panel's outward normal.
    /// For panel i (from node i to node i+1):
    ///   APANEL[i] = atan2(SX, -SY)
    /// where SX = X[i+1] - X[i], SY = Y[i+1] - Y[i]
    fn compute_panel_angles(x: &[f64], y: &[f64], nx: &[f64], ny: &[f64]) -> Vec<f64> {
        let n = x.len();
        let mut apanel = vec![0.0; n];

        // Regular panels (i = 0 to n-2)
        for i in 0..n - 1 {
            let sx = x[i + 1] - x[i];
            let sy = y[i + 1] - y[i];

            if sx.abs() < 1e-12 && sy.abs() < 1e-12 {
                // Zero-length panel: use node normal
                apanel[i] = (-ny[i]).atan2(-nx[i]);
            } else {
                // Normal formula: atan2(SX, -SY)
                apanel[i] = sx.atan2(-sy);
            }
        }

        // TE panel (from node n-1 to node 0)
        let sx = x[0] - x[n - 1];
        let sy = y[0] - y[n - 1];
        
        // XFOIL: APANEL(N) = ATAN2(-SX, SY) + PI
        apanel[n - 1] = (-sx).atan2(sy) + PI;

        apanel
    }

    /// Compute trailing edge geometry (XFOIL's TECALC).
    ///
    /// Returns (xte, yte, dste, ante, aste) where:
    /// - xte, yte: TE midpoint coordinates
    /// - dste: TE gap length
    /// - ante: normal-projected gap (ANTE)
    /// - aste: tangent-projected gap (ASTE)
    fn compute_te_geometry(x: &[f64], y: &[f64], apanel: &[f64]) -> (f64, f64, f64, f64, f64) {
        let n = x.len();

        // TE gap vector (from upper TE to lower TE)
        let dxte = x[n - 1] - x[0];
        let dyte = y[n - 1] - y[0];
        let dste = (dxte * dxte + dyte * dyte).sqrt();

        // TE midpoint
        let xte = 0.5 * (x[0] + x[n - 1]);
        let yte = 0.5 * (y[0] + y[n - 1]);

        // Mean TE tangent direction (average of panel 0 and panel n-2 tangents)
        // Panel 0 tangent: from node 0 to node 1
        // Panel n-2 tangent: from node n-2 to node n-1 (but we want outward, so negate)
        let ap0 = apanel[0];
        let apn = if n >= 2 { apanel[n - 2] } else { ap0 };

        // XFOIL computes the mean tangent as average of adjacent panel tangents
        // tangent direction is perpendicular to panel normal
        let t0x = ap0.sin();  // tangent of panel 0
        let t0y = -ap0.cos();
        let tnx = -apn.sin(); // tangent of last panel (reversed for consistency)
        let tny = apn.cos();

        let dxs = 0.5 * (t0x + tnx);
        let dys = 0.5 * (t0y + tny);

        // ANTE: normal-projected gap (gap · normal direction)
        // ASTE: tangent-projected gap (gap · tangent direction)
        let ante = dxs * dyte - dys * dxte;
        let aste = dxs * dxte + dys * dyte;

        (xte, yte, dste, ante, aste)
    }

    /// Compute chord length from the true leading/trailing-edge geometry.
    ///
    /// XFOIL sets:
    ///   CHORD = SQRT((XTE-XLE)^2 + (YTE-YLE)^2)
    fn compute_chord(xte: f64, yte: f64, xle: f64, yle: f64) -> f64 {
        let dx = xte - xle;
        let dy = yte - yle;
        (dx * dx + dy * dy).sqrt()
    }

    /// Find leading edge location (XFOIL's LEFIND).
    ///
    /// The LE is defined as the point where the surface tangent is
    /// perpendicular to the chord line (TE to LE).
    fn find_leading_edge(
        x: &[f64],
        y: &[f64],
        s: &[f64],
        xp: &[f64],
        yp: &[f64],
        xte: f64,
        yte: f64,
    ) -> (f64, f64, f64) {
        let n = x.len();

        // Initial guess: find where dot product with TE changes sign
        let mut i_le = n / 2;
        for i in 2..n - 2 {
            let dx_te = x[i] - xte;
            let dy_te = y[i] - yte;
            let dx = x[i + 1] - x[i];
            let dy = y[i + 1] - y[i];
            let dotp = dx_te * dx + dy_te * dy;
            if dotp < 0.0 {
                i_le = i;
                break;
            }
        }

        let mut s_le = s[i_le];
        let dseps = (s[n - 1] - s[0]) * 1e-5;

        // Newton iteration for exact SLE
        for _iter in 0..50 {
            // Evaluate spline at s_le
            let (x_le, y_le) = Self::seval_point(s_le, s, x, y, xp, yp);
            let (dxds, dyds) = Self::deval_point(s_le, s, x, y, xp, yp);
            let (dxdd, dydd) = Self::d2val_point(s_le, s, x, xp, y, yp);

            let x_chord = x_le - xte;
            let y_chord = y_le - yte;

            // Drive dot product between chord line and LE tangent to zero
            let res = x_chord * dxds + y_chord * dyds;
            let ress = dxds * dxds + dyds * dyds + x_chord * dxdd + y_chord * dydd;

            if ress.abs() < 1e-20 {
                break;
            }

            let mut ds_le = -res / ress;

            // Match XFOIL LEFIND exactly: limit the Newton step using
            // ABS(XCHORD + YCHORD) without an additional floor term.
            let chord_scale = (x_chord + y_chord).abs();
            ds_le = ds_le.max(-0.02 * chord_scale).min(0.02 * chord_scale);
            s_le += ds_le;

            if ds_le.abs() < dseps {
                break;
            }
        }

        let (xle, yle) = Self::seval_point(s_le, s, x, y, xp, yp);
        (xle, yle, s_le)
    }

    /// Evaluate spline at parameter ss (1D).
    fn seval_1d(ss: f64, s: &[f64], f: &[f64], fs: &[f64]) -> f64 {
        let n = s.len();
        if n < 2 {
            return f.first().copied().unwrap_or(0.0);
        }

        let i = Self::find_segment(ss, s);

        let ds = s[i] - s[i - 1];
        let t = (ss - s[i - 1]) / ds;
        let cx1 = ds * fs[i - 1] - f[i] + f[i - 1];
        let cx2 = ds * fs[i] - f[i] + f[i - 1];

        t * f[i] + (1.0 - t) * f[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
    }

    /// Evaluate spline derivative at parameter ss (1D).
    fn deval_1d(ss: f64, s: &[f64], f: &[f64], fs: &[f64]) -> f64 {
        let n = s.len();
        if n < 2 {
            return 0.0;
        }

        let i = Self::find_segment(ss, s);

        let ds = s[i] - s[i - 1];
        let t = (ss - s[i - 1]) / ds;
        let cx1 = ds * fs[i - 1] - f[i] + f[i - 1];
        let cx2 = ds * fs[i] - f[i] + f[i - 1];

        let deval = f[i] - f[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
        deval / ds
    }

    /// Evaluate spline second derivative at parameter ss (1D).
    fn d2val_1d(ss: f64, s: &[f64], f: &[f64], fs: &[f64]) -> f64 {
        let n = s.len();
        if n < 2 {
            return 0.0;
        }

        let i = Self::find_segment(ss, s);

        let ds = s[i] - s[i - 1];
        let t = (ss - s[i - 1]) / ds;
        let cx1 = ds * fs[i - 1] - f[i] + f[i - 1];
        let cx2 = ds * fs[i] - f[i] + f[i - 1];

        let d2val = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
        d2val / (ds * ds)
    }

    /// Evaluate point on spline.
    fn seval_point(ss: f64, s: &[f64], x: &[f64], y: &[f64], xp: &[f64], yp: &[f64]) -> (f64, f64) {
        let x_val = Self::seval_1d(ss, s, x, xp);
        let y_val = Self::seval_1d(ss, s, y, yp);
        (x_val, y_val)
    }

    /// Evaluate derivative on spline.
    fn deval_point(ss: f64, s: &[f64], x: &[f64], y: &[f64], xp: &[f64], yp: &[f64]) -> (f64, f64) {
        let dx = Self::deval_1d(ss, s, x, xp);
        let dy = Self::deval_1d(ss, s, y, yp);
        (dx, dy)
    }

    /// Evaluate second derivative on spline.
    fn d2val_point(ss: f64, s: &[f64], x: &[f64], xp: &[f64], y: &[f64], yp: &[f64]) -> (f64, f64) {
        let d2x = Self::d2val_1d(ss, s, x, xp);
        let d2y = Self::d2val_1d(ss, s, y, yp);
        (d2x, d2y)
    }

    /// Find segment index for parameter ss (binary search).
    fn find_segment(ss: f64, s: &[f64]) -> usize {
        let n = s.len();
        if n < 2 {
            return 1;
        }

        if ss <= s[0] {
            return 1;
        }
        if ss >= s[n - 1] {
            return n - 1;
        }

        let mut ilow = 0;
        let mut i = n - 1;

        while i - ilow > 1 {
            let imid = (i + ilow) / 2;
            if ss < s[imid] {
                i = imid;
            } else {
                ilow = imid;
            }
        }

        i
    }

    /// Get TE panel source/vortex coefficients (SCS, SDS).
    ///
    /// For blunt TE: SCS = ANTE/DSTE, SDS = ASTE/DSTE
    /// For sharp TE: SCS = 1.0, SDS = 0.0
    pub fn te_coefficients(&self) -> (f64, f64) {
        if self.sharp {
            (1.0, 0.0)
        } else {
            (self.ante / self.dste, self.aste / self.dste)
        }
    }

    /// Return XFOIL's sharp trailing-edge bisector control point and normal.
    ///
    /// The control point sits slightly inside the TE corner along the bisector,
    /// and the returned normal corresponds to the tangential-velocity probe used
    /// by XFOIL's sharp-TE row in `GGCALC/QDCALC`.
    pub fn sharp_te_bisector_control(&self) -> Option<(f64, f64, f64, f64)> {
        if !self.sharp || self.n < 3 {
            return None;
        }

        let upper_tx = -self.xp[0];
        let upper_ty = -self.yp[0];
        let lower_tx = self.xp[self.n - 1];
        let lower_ty = self.yp[self.n - 1];

        let bis_x = upper_tx + lower_tx;
        let bis_y = upper_ty + lower_ty;
        let bis_norm = (bis_x * bis_x + bis_y * bis_y).sqrt().max(1.0e-12);
        let cbis = bis_x / bis_norm;
        let sbis = bis_y / bis_norm;

        let ds1 = ((self.x[0] - self.x[1]).powi(2) + (self.y[0] - self.y[1]).powi(2)).sqrt();
        let ds2 = ((self.x[self.n - 1] - self.x[self.n - 2]).powi(2)
            + (self.y[self.n - 1] - self.y[self.n - 2]).powi(2))
        .sqrt();
        let dsmin = ds1.min(ds2);
        let bwt = 0.1;

        let xbis = self.xte - bwt * dsmin * cbis;
        let ybis = self.yte - bwt * dsmin * sbis;

        Some((xbis, ybis, -sbis, cbis))
    }

    /// Get total arc length.
    pub fn total_arc_length(&self) -> f64 {
        self.s.last().copied().unwrap_or(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Generate NACA 0012 coordinates with cosine spacing.
    fn make_naca0012(n_panels: usize) -> Vec<(f64, f64)> {
        let n_half = n_panels / 2;
        let t = 0.12; // thickness ratio
        
        // Generate x coordinates with cosine spacing
        let x_coords: Vec<f64> = (0..=n_half)
            .map(|i| {
                let beta = PI * (i as f64) / (n_half as f64);
                0.5 * (1.0 - beta.cos())
            })
            .collect();

        // NACA 0012 thickness distribution (closed TE)
        let thickness = |x: f64| -> f64 {
            5.0 * t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) 
                + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
        };

        let mut points = Vec::with_capacity(2 * n_half);

        // Upper surface: TE to LE (x decreasing)
        for i in (0..=n_half).rev() {
            let x = x_coords[i];
            let y = thickness(x);
            points.push((x, y));
        }

        // Lower surface: LE+1 to TE (x increasing, skip LE to avoid duplicate)
        for i in 1..=n_half {
            let x = x_coords[i];
            let y = -thickness(x);
            points.push((x, y));
        }

        points
    }

    #[test]
    fn test_geometry_creation() {
        let points = make_naca0012(80);
        let geom = AirfoilGeometry::from_points(&points).unwrap();

        // Basic sanity checks
        assert_eq!(geom.n, points.len());
        assert!(geom.chord > 0.9 && geom.chord < 1.1, "Chord should be ~1.0");
        assert!(geom.dste < 0.01, "TE gap should be small for closed airfoil");
    }

    #[test]
    fn test_normals_outward() {
        let points = make_naca0012(80);
        let geom = AirfoilGeometry::from_points(&points).unwrap();

        // Upper surface normals should point upward (positive ny)
        // Lower surface normals should point downward (negative ny)
        let n_upper = geom.n / 2;
        
        // Check a few upper surface points (near middle, not at LE/TE)
        for i in 5..n_upper - 5 {
            assert!(geom.ny[i] > 0.0, "Upper surface normal at {} should point up", i);
        }
        
        // Check a few lower surface points
        for i in n_upper + 5..geom.n - 5 {
            assert!(geom.ny[i] < 0.0, "Lower surface normal at {} should point down", i);
        }
    }

    #[test]
    fn test_arc_length_monotonic() {
        let points = make_naca0012(80);
        let geom = AirfoilGeometry::from_points(&points).unwrap();

        // Arc length should be strictly increasing
        for i in 1..geom.n {
            assert!(geom.s[i] > geom.s[i - 1], 
                "Arc length not monotonic at {}: {} <= {}", 
                i, geom.s[i], geom.s[i - 1]);
        }
    }

    #[test]
    fn test_leading_edge_location() {
        let points = make_naca0012(80);
        let geom = AirfoilGeometry::from_points(&points).unwrap();

        // LE should be near x=0
        assert!(geom.xle.abs() < 0.01, "LE x should be near 0, got {}", geom.xle);
        assert!(geom.yle.abs() < 0.01, "LE y should be near 0, got {}", geom.yle);
    }

    #[test]
    fn test_te_coefficients() {
        let points = make_naca0012(80);
        let geom = AirfoilGeometry::from_points(&points).unwrap();

        let (scs, sds) = geom.te_coefficients();
        
        // For closed airfoil with sharp TE
        if geom.sharp {
            assert!((scs - 1.0).abs() < 1e-10);
            assert!(sds.abs() < 1e-10);
        } else {
            // SCS and SDS should be finite
            assert!(scs.is_finite());
            assert!(sds.is_finite());
        }
    }
}
