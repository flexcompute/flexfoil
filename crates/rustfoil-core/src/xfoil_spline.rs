//! XFOIL-compatible Hermite cubic spline implementation.
//!
//! This module provides spline interpolation that matches XFOIL's `spline.f` exactly.
//! XFOIL uses Hermite cubic splines where the spline coefficients are first derivatives
//! (dX/dS) at each knot, solved via a tridiagonal system.
//!
//! Reference: XFOIL source `spline.f` (Mark Drela, 2000)

use crate::point::Point;

// ============================================================================
// 1D Hermite Spline (for scalar values like curvature)
// ============================================================================

/// XFOIL-compatible 1D Hermite cubic spline for scalar values.
///
/// Used for interpolating curvature values during paneling (W5/W6 arrays in XFOIL).
#[derive(Debug, Clone)]
pub struct Spline1D {
    /// Parameter values at each knot
    s: Vec<f64>,
    /// Function values at each knot
    f: Vec<f64>,
    /// dF/dS at each knot (computed by SPLINE)
    fs: Vec<f64>,
}

impl Spline1D {
    /// Build 1D spline from parameter-value pairs using XFOIL's SEGSPL routine.
    /// Uses zero third derivative end conditions (like XFOIL's SEGSPL with XS=-999).
    pub fn new(s: &[f64], f: &[f64]) -> Option<Self> {
        let n = s.len();
        if n < 2 || f.len() != n {
            return None;
        }

        let fs = Self::spline_coeffs_segspl(s, f);

        Some(Self {
            s: s.to_vec(),
            f: f.to_vec(),
            fs,
        })
    }

    /// XFOIL's SEGSPL routine - uses zero THIRD derivative end conditions.
    /// This matches XFOIL's behavior for curvature splines (SEGSPL with XS=-999).
    fn spline_coeffs_segspl(s: &[f64], f: &[f64]) -> Vec<f64> {
        let n = s.len();
        if n < 2 {
            return vec![0.0; n];
        }
        if n == 2 {
            let df = (f[1] - f[0]) / (s[1] - s[0]);
            return vec![df, df];
        }

        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n];
        let mut fs = vec![0.0; n];

        // Interior points (XFOIL lines 39-46 / 87-94)
        for i in 1..n - 1 {
            let dsm = s[i] - s[i - 1];
            let dsp = s[i + 1] - s[i];
            b[i] = dsp;
            a[i] = 2.0 * (dsm + dsp);
            c[i] = dsm;
            fs[i] = 3.0 * ((f[i + 1] - f[i]) * dsm / dsp + (f[i] - f[i - 1]) * dsp / dsm);
        }

        // Zero THIRD derivative end conditions (XFOIL lines 101-105, 117-120)
        // This is what SEGSPL uses with XS1 = XS2 = -999.0
        a[0] = 1.0;
        c[0] = 1.0;
        fs[0] = 2.0 * (f[1] - f[0]) / (s[1] - s[0]);

        b[n - 1] = 1.0;
        a[n - 1] = 1.0;
        fs[n - 1] = 2.0 * (f[n - 1] - f[n - 2]) / (s[n - 1] - s[n - 2]);

        // Solve tridiagonal system
        Self::trisol(&a, &b, &c, &mut fs);

        fs
    }

    /// XFOIL's TRISOL routine.
    fn trisol(a: &[f64], b: &[f64], c: &[f64], x: &mut [f64]) {
        let n = a.len();
        if n < 2 {
            return;
        }

        let mut aa = a.to_vec();
        for i in 1..n {
            let piv = b[i] / aa[i - 1];
            aa[i] = aa[i] - c[i - 1] * piv;
            x[i] = x[i] - x[i - 1] * piv;
        }

        x[n - 1] = x[n - 1] / aa[n - 1];
        for i in (0..n - 1).rev() {
            x[i] = (x[i] - c[i] * x[i + 1]) / aa[i];
        }
    }

    /// XFOIL's SEVAL - evaluate spline at parameter ss.
    pub fn eval(&self, ss: f64) -> f64 {
        let n = self.s.len();
        if n < 2 {
            return self.f.first().copied().unwrap_or(0.0);
        }

        let i = self.find_segment(ss);

        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;
        let cx1 = ds * self.fs[i - 1] - self.f[i] + self.f[i - 1];
        let cx2 = ds * self.fs[i] - self.f[i] + self.f[i - 1];

        t * self.f[i] + (1.0 - t) * self.f[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
    }

    /// XFOIL's DEVAL - evaluate first derivative at parameter ss.
    pub fn deriv(&self, ss: f64) -> f64 {
        let n = self.s.len();
        if n < 2 {
            return 0.0;
        }

        let i = self.find_segment(ss);

        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;
        let cx1 = ds * self.fs[i - 1] - self.f[i] + self.f[i - 1];
        let cx2 = ds * self.fs[i] - self.f[i] + self.f[i - 1];

        let deval = self.f[i] - self.f[i - 1] 
            + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 
            + t * (3.0 * t - 2.0) * cx2;
        deval / ds
    }

    /// Find segment index for parameter ss (binary search).
    fn find_segment(&self, ss: f64) -> usize {
        let n = self.s.len();
        if n < 2 {
            return 1;
        }

        if ss <= self.s[0] {
            return 1;
        }
        if ss >= self.s[n - 1] {
            return n - 1;
        }

        let mut ilow = 0;
        let mut i = n - 1;

        while i - ilow > 1 {
            let imid = (i + ilow) / 2;
            if ss < self.s[imid] {
                i = imid;
            } else {
                ilow = imid;
            }
        }

        i
    }
}

// ============================================================================
// 2D Parametric Hermite Spline (for airfoil geometry)
// ============================================================================

/// XFOIL-compatible parametric Hermite cubic spline.
///
/// Stores x(s), y(s) and their first derivatives xs(s), ys(s) at each knot.
#[derive(Debug, Clone)]
pub struct XfoilSpline {
    /// Arc-length parameter values at each knot
    s: Vec<f64>,
    /// X coordinates at each knot
    x: Vec<f64>,
    /// Y coordinates at each knot  
    y: Vec<f64>,
    /// dX/dS at each knot (computed by SPLINE)
    xs: Vec<f64>,
    /// dY/dS at each knot (computed by SPLINE)
    ys: Vec<f64>,
}

impl XfoilSpline {
    /// Build spline from points using XFOIL's SPLINE routine.
    ///
    /// This matches XFOIL `spline.f` lines 21-60 exactly.
    /// Uses zero second derivative end conditions.
    pub fn from_points(points: &[Point]) -> Option<Self> {
        let n = points.len();
        if n < 2 {
            return None;
        }

        // Compute arc-length parameterization (XFOIL's SCALC)
        let mut s = vec![0.0; n];
        for i in 1..n {
            let dx = points[i].x - points[i - 1].x;
            let dy = points[i].y - points[i - 1].y;
            s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
        }

        let x: Vec<f64> = points.iter().map(|p| p.x).collect();
        let y: Vec<f64> = points.iter().map(|p| p.y).collect();

        // Compute spline derivatives using XFOIL's SPLINE routine
        let xs = Self::spline_coeffs(&s, &x);
        let ys = Self::spline_coeffs(&s, &y);

        Some(Self { s, x, y, xs, ys })
    }

    /// XFOIL's SPLINE routine - computes first derivative coefficients.
    /// 
    /// Reference: spline.f lines 21-60
    /// Zero second derivative end conditions.
    fn spline_coeffs(s: &[f64], x: &[f64]) -> Vec<f64> {
        let n = s.len();
        if n < 2 {
            return vec![0.0; n];
        }
        if n == 2 {
            let dx = (x[1] - x[0]) / (s[1] - s[0]);
            return vec![dx, dx];
        }

        // Build tridiagonal system (XFOIL lines 39-46)
        // A[i]*XS[i-1] + B[i]*XS[i] + C[i]*XS[i+1] = RHS[i]
        // But XFOIL uses: B[i]=dsp, A[i]=2*(dsm+dsp), C[i]=dsm
        let mut a = vec![0.0; n];  // Main diagonal
        let mut b = vec![0.0; n];  // Upper diagonal (offset)
        let mut c = vec![0.0; n];  // Lower diagonal (offset)
        let mut xs = vec![0.0; n]; // RHS, then solution

        for i in 1..n - 1 {
            let dsm = s[i] - s[i - 1];
            let dsp = s[i + 1] - s[i];
            b[i] = dsp;                    // XFOIL: B(I) = DSP
            a[i] = 2.0 * (dsm + dsp);      // XFOIL: A(I) = 2.0*(DSM+DSP)
            c[i] = dsm;                    // XFOIL: C(I) = DSM
            xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
        }

        // Zero second derivative end conditions (XFOIL lines 48-54)
        a[0] = 2.0;
        c[0] = 1.0;
        xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);

        b[n - 1] = 1.0;
        a[n - 1] = 2.0;
        xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);

        // Solve tridiagonal system (XFOIL's TRISOL)
        Self::trisol(&a, &b, &c, &mut xs);

        xs
    }

    /// XFOIL's TRISOL routine - tridiagonal matrix solver.
    /// 
    /// Reference: spline.f lines 182-215
    /// Solves: B[i]*X[i-1] + A[i]*X[i] + C[i]*X[i+1] = X[i] (RHS)
    fn trisol(a: &[f64], b: &[f64], c: &[f64], x: &mut [f64]) {
        let n = a.len();
        if n < 2 {
            return;
        }

        // XFOIL's TRISOL forward sweep (lines 195-201)
        let mut aa = a.to_vec();
        for i in 1..n {
            let piv = b[i] / aa[i - 1];
            aa[i] = aa[i] - c[i - 1] * piv;
            x[i] = x[i] - x[i - 1] * piv;
        }

        // Back substitution (XFOIL lines 203-206)
        x[n - 1] = x[n - 1] / aa[n - 1];
        for i in (0..n - 1).rev() {
            x[i] = (x[i] - c[i] * x[i + 1]) / aa[i];
        }
    }

    /// Total arc length of the spline.
    pub fn total_arc_length(&self) -> f64 {
        *self.s.last().unwrap_or(&0.0)
    }

    /// Number of knots.
    pub fn len(&self) -> usize {
        self.s.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.s.is_empty()
    }

    /// XFOIL's SEVAL - evaluate spline at parameter ss.
    ///
    /// Reference: spline.f lines 218-243
    pub fn seval(&self, ss: f64) -> Point {
        let x = self.seval_1d(ss, &self.x, &self.xs);
        let y = self.seval_1d(ss, &self.y, &self.ys);
        Point::new(x, y)
    }

    /// 1D spline evaluation (XFOIL's SEVAL)
    fn seval_1d(&self, ss: f64, x: &[f64], xs: &[f64]) -> f64 {
        let n = self.s.len();
        if n < 2 {
            return x.first().copied().unwrap_or(0.0);
        }

        // Binary search for segment (XFOIL lines 224-235)
        let i = self.find_segment(ss);

        // Evaluate (XFOIL lines 237-241)
        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;
        let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
        let cx2 = ds * xs[i] - x[i] + x[i - 1];
        
        // SEVAL = T*X(I) + (1.0-T)*X(I-1) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
        t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
    }

    /// XFOIL's DEVAL - evaluate first derivative at parameter ss.
    ///
    /// Reference: spline.f lines 245-271
    pub fn deval(&self, ss: f64) -> (f64, f64) {
        let dx = self.deval_1d(ss, &self.x, &self.xs);
        let dy = self.deval_1d(ss, &self.y, &self.ys);
        (dx, dy)
    }

    /// 1D derivative evaluation (XFOIL's DEVAL)
    fn deval_1d(&self, ss: f64, x: &[f64], xs: &[f64]) -> f64 {
        let n = self.s.len();
        if n < 2 {
            return 0.0;
        }

        let i = self.find_segment(ss);

        // XFOIL lines 264-269
        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;
        let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
        let cx2 = ds * xs[i] - x[i] + x[i - 1];

        // DEVAL = (X(I) - X(I-1) + (1.-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.)*CX2) / DS
        let deval = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
        deval / ds
    }

    /// XFOIL's D2VAL - evaluate second derivative at parameter ss.
    ///
    /// Reference: spline.f lines 273-299
    pub fn d2val(&self, ss: f64) -> (f64, f64) {
        let d2x = self.d2val_1d(ss, &self.x, &self.xs);
        let d2y = self.d2val_1d(ss, &self.y, &self.ys);
        (d2x, d2y)
    }

    /// 1D second derivative evaluation (XFOIL's D2VAL)
    fn d2val_1d(&self, ss: f64, x: &[f64], xs: &[f64]) -> f64 {
        let n = self.s.len();
        if n < 2 {
            return 0.0;
        }

        let i = self.find_segment(ss);

        // XFOIL lines 292-297
        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;
        let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
        let cx2 = ds * xs[i] - x[i] + x[i - 1];

        // D2VAL = ((6.*T-4.)*CX1 + (6.*T-2.0)*CX2) / DS**2
        let d2val = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
        d2val / (ds * ds)
    }

    /// XFOIL's CURV - evaluate curvature at parameter ss.
    ///
    /// Reference: spline.f lines 302-344
    pub fn curvature(&self, ss: f64) -> f64 {
        let n = self.s.len();
        if n < 2 {
            return 0.0;
        }

        let i = self.find_segment(ss);

        // XFOIL lines 327-343
        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;

        let cx1 = ds * self.xs[i - 1] - self.x[i] + self.x[i - 1];
        let cx2 = ds * self.xs[i] - self.x[i] + self.x[i - 1];
        let xd = self.x[i] - self.x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
        let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

        let cy1 = ds * self.ys[i - 1] - self.y[i] + self.y[i - 1];
        let cy2 = ds * self.ys[i] - self.y[i] + self.y[i - 1];
        let yd = self.y[i] - self.y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
        let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;

        // SD = SQRT(XD*XD + YD*YD)
        // SD = MAX(SD, 0.001*DS)
        // CURV = (XD*YDD - YD*XDD) / SD**3
        let sd = (xd * xd + yd * yd).sqrt().max(0.001 * ds);
        (xd * ydd - yd * xdd) / (sd * sd * sd)
    }

    /// XFOIL's LEFIND - find leading edge arc-length position.
    ///
    /// The LE is defined as the point where the surface tangent is
    /// perpendicular to the chord line (TE to LE).
    ///
    /// Reference: xgeom.f lines 21-87
    pub fn lefind(&self) -> f64 {
        let n = self.s.len();
        if n < 5 {
            return self.s[n / 2];
        }

        let dseps = (self.s[n - 1] - self.s[0]) * 1.0e-5;

        // Trailing edge point
        let x_te = 0.5 * (self.x[0] + self.x[n - 1]);
        let y_te = 0.5 * (self.y[0] + self.y[n - 1]);

        // Get first guess: find where dot product with TE changes sign
        let mut i_le = n / 2;
        for i in 2..n - 2 {
            let dx_te = self.x[i] - x_te;
            let dy_te = self.y[i] - y_te;
            let dx = self.x[i + 1] - self.x[i];
            let dy = self.y[i + 1] - self.y[i];
            let dotp = dx_te * dx + dy_te * dy;
            if dotp < 0.0 {
                i_le = i;
                break;
            }
        }

        let mut s_le = self.s[i_le];

        // Check for sharp LE (doubled point)
        if i_le > 0 && (self.s[i_le] - self.s[i_le - 1]).abs() < 1e-12 {
            return s_le;
        }

        // Newton iteration for exact SLE
        for _iter in 0..50 {
            let pt = self.seval(s_le);
            let (dxds, dyds) = self.deval(s_le);
            let (dxdd, dydd) = self.d2val(s_le);

            let x_chord = pt.x - x_te;
            let y_chord = pt.y - y_te;

            // Drive dot product between chord line and LE tangent to zero
            let res = x_chord * dxds + y_chord * dyds;
            let ress = dxds * dxds + dyds * dyds + x_chord * dxdd + y_chord * dydd;

            if ress.abs() < 1e-20 {
                break;
            }

            let mut ds_le = -res / ress;

            // Limit step size
            let chord_scale = (x_chord.abs() + y_chord.abs()).max(0.01);
            ds_le = ds_le.max(-0.02 * chord_scale).min(0.02 * chord_scale);
            s_le += ds_le;

            if ds_le.abs() < dseps {
                break;
            }
        }

        s_le
    }

    /// Find segment index for parameter ss (binary search).
    fn find_segment(&self, ss: f64) -> usize {
        let n = self.s.len();
        if n < 2 {
            return 1;
        }

        // Clamp to valid range
        if ss <= self.s[0] {
            return 1;
        }
        if ss >= self.s[n - 1] {
            return n - 1;
        }

        // Binary search (XFOIL lines 224-235)
        let mut ilow = 0;
        let mut i = n - 1;

        while i - ilow > 1 {
            let imid = (i + ilow) / 2;
            if ss < self.s[imid] {
                i = imid;
            } else {
                ilow = imid;
            }
        }

        i
    }

    /// Resample spline at n uniformly-spaced arc-length positions.
    pub fn resample_uniform(&self, n: usize) -> Vec<Point> {
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![self.seval(0.0)];
        }

        let s_max = self.total_arc_length();
        let ds = s_max / (n - 1) as f64;

        (0..n).map(|i| self.seval(i as f64 * ds)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::point;

    #[test]
    fn test_xfoil_spline_endpoints() {
        let points = vec![point(0.0, 0.0), point(1.0, 1.0), point(2.0, 0.0)];
        let spline = XfoilSpline::from_points(&points).unwrap();

        // Should pass through endpoints
        let p0 = spline.seval(0.0);
        assert!((p0.x - 0.0).abs() < 1e-10);
        assert!((p0.y - 0.0).abs() < 1e-10);

        let p_end = spline.seval(spline.total_arc_length());
        assert!((p_end.x - 2.0).abs() < 1e-10);
        assert!((p_end.y - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_xfoil_spline_curvature() {
        // Circle segment - constant curvature
        let n = 50;
        let pi = std::f64::consts::PI;
        let r = 1.0;
        
        let points: Vec<_> = (0..n)
            .map(|i| {
                let theta = pi * i as f64 / (n - 1) as f64;
                point(r * theta.cos(), r * theta.sin())
            })
            .collect();

        let spline = XfoilSpline::from_points(&points).unwrap();
        
        // Curvature of circle with radius 1 should be 1
        let mid_s = spline.total_arc_length() / 2.0;
        let k = spline.curvature(mid_s);
        assert!((k.abs() - 1.0).abs() < 0.05, "Expected curvature ~1.0, got {}", k);
    }
}
