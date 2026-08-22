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

    /// Total arc-length range of the spline.
    pub fn s_range(&self) -> f64 {
        if self.s.len() < 2 { return 0.0; }
        (self.s.last().unwrap() - self.s.first().unwrap()).abs()
    }

    /// Second derivative (D2VAL equivalent for 1D spline).
    pub fn d2val(&self, ss: f64) -> f64 {
        let n = self.s.len();
        if n < 2 { return 0.0; }
        let i = self.find_segment(ss);
        let ds = self.s[i] - self.s[i - 1];
        let t = (ss - self.s[i - 1]) / ds;
        let cx1 = ds * self.fs[i - 1] - self.f[i] + self.f[i - 1];
        let cx2 = ds * self.fs[i] - self.f[i] + self.f[i - 1];
        let d2f = 2.0 * cx1 * (1.0 - 3.0 * t) - 2.0 * cx2 * (2.0 - 3.0 * t);
        d2f / (ds * ds)
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

    /// XFOIL's SEGSPL-style spline - computes first derivative coefficients.
    /// 
    /// Reference: spline.f SEGSPL which calls SPLIND with XS1=XS2=-999.0
    /// This uses zero THIRD derivative end conditions, which is what XFOIL's
    /// PANGEN uses for geometry splines (xfoil.f lines 1662-1663).
    /// 
    /// NOTE: This differs from XFOIL's plain SPLINE routine which uses zero
    /// SECOND derivative end conditions. Using SEGSPL-style end conditions
    /// is critical for matching XFOIL's curvature computation in PANGEN.
    fn spline_coeffs(s: &[f64], x: &[f64]) -> Vec<f64> {
        let n = s.len();
        if n < 2 {
            return vec![0.0; n];
        }
        if n == 2 {
            let dx = (x[1] - x[0]) / (s[1] - s[0]);
            return vec![dx, dx];
        }

        // Build tridiagonal system (XFOIL lines 87-94 in SPLIND)
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

        // Zero THIRD derivative end conditions (XFOIL SPLIND lines 101-105, 117-120)
        // This is what SEGSPL uses with XS1 = XS2 = -999.0
        // This is CRITICAL for matching XFOIL's PANGEN curvature computation
        a[0] = 1.0;
        c[0] = 1.0;
        xs[0] = 2.0 * (x[1] - x[0]) / (s[1] - s[0]);

        b[n - 1] = 1.0;
        a[n - 1] = 1.0;
        xs[n - 1] = 2.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);

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

    /// Node furthest from the trailing-edge midpoint.
    ///
    /// The leading edge is the point of the contour furthest from the
    /// trailing-edge midpoint, so the node where that distance is largest is
    /// within one panel of it. One O(n) pass, and it reads the whole contour, so
    /// a straight or vertical stretch elsewhere on the surface cannot stand in
    /// for the leading edge. Ties go to the lower index.
    fn node_furthest_from(&self, x_te: f64, y_te: f64) -> usize {
        let mut i_far = 0usize;
        let mut d_far = f64::NEG_INFINITY;
        for i in 0..self.s.len() {
            let dx = self.x[i] - x_te;
            let dy = self.y[i] - y_te;
            let d = dx * dx + dy * dy;
            if d > d_far {
                d_far = d;
                i_far = i;
            }
        }
        i_far
    }

    /// XFOIL's LEFIND - find leading edge arc-length position.
    ///
    /// The LE is defined as the point where the surface tangent is
    /// perpendicular to the chord line (TE to LE).
    ///
    /// Reference: xgeom.f lines 21-87
    ///
    /// # Where the Newton iteration starts, and what it is allowed to return
    /// XFOIL's initial guess is a forward scan for the first node whose step
    /// away from the trailing edge turns back towards it. That is a *local*
    /// test, and it stops at the first place the contour turns back — on a
    /// conventional airfoil the leading edge, but on an element with a deep cove
    /// the cove face, which for the main element of a slat/main/flap section is
    /// most of a chord downstream of the leading edge. The Newton iteration then
    /// converges to a stationary point of the same residual on that face and
    /// reports it as the leading edge.
    ///
    /// So the scan is paired with [`node_furthest_from`](Self::node_furthest_from),
    /// which is global:
    ///
    /// - where the scan node and the furthest-from-trailing-edge node agree to
    ///   within one node — every contour whose first turn back *is* the leading
    ///   edge — the iteration starts from the scan node and nothing else
    ///   changes, so the result is the same value bit for bit;
    /// - where they disagree, at most one of them is the leading edge, and
    ///   neither rule settles it in general: the scan stops at the first turn
    ///   back, which a cove face aft of the leading edge also satisfies, while
    ///   the furthest node is not the leading edge on every contour either. So
    ///   the iteration is run from both and the candidate genuinely furthest
    ///   from the trailing edge wins. Choosing on the quantity itself means
    ///   neither heuristic has to be right, and a correct scan result is never
    ///   discarded in favour of a worse one.
    ///
    /// `AirfoilGeometry::find_leading_edge` in `rustfoil-inviscid` seeds the same
    /// iteration the same way, and agrees with this one bit for bit on the
    /// sections both pin. `lefind_is_bit_identical_on_conventional_sections`
    /// pins the unchanged case, and
    /// `lefind_reports_the_furthest_point_from_the_trailing_edge` holds the
    /// property over the whole bundled coordinate library.
    pub fn lefind(&self) -> f64 {
        let n = self.s.len();
        if n < 5 {
            return self.s[n / 2];
        }

        let dseps = (self.s[n - 1] - self.s[0]) * 1.0e-5;

        // Trailing edge point
        let x_te = 0.5 * (self.x[0] + self.x[n - 1]);
        let y_te = 0.5 * (self.y[0] + self.y[n - 1]);

        // Node furthest from the trailing edge: a global candidate, independent
        // of where the contour first turns back.
        let i_far = self.node_furthest_from(x_te, y_te);

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

        // Newton iteration for exact SLE, unchanged, run from a given start.
        let refine = |mut s_le: f64| -> f64 {
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

                // Limit step size - XFOIL uses ABS(XCHORD+YCHORD), not sum of absolutes
                // This matches XFOIL xgeom.f lines 79-80 exactly
                let chord_scale = (x_chord + y_chord).abs();
                ds_le = ds_le.max(-0.02 * chord_scale).min(0.02 * chord_scale);
                s_le += ds_le;

                if ds_le.abs() < dseps {
                    break;
                }
            }
            s_le
        };

        // Refine from a node, keeping the sharp-LE check: a zero-length segment
        // ending at the start node is a leading edge written as a doubled point,
        // and there is nothing there for the iteration to refine.
        let refine_from_node = |i: usize| -> f64 {
            if i > 0 && (self.s[i] - self.s[i - 1]).abs() < 1e-12 {
                return self.s[i];
            }
            refine(self.s[i])
        };

        // Distance from the trailing edge, which is the quantity "leading edge"
        // actually names.
        let te_distance = |s_at: f64| -> f64 {
            let pt = self.seval(s_at);
            (pt.x - x_te).powi(2) + (pt.y - y_te).powi(2)
        };

        if i_le.max(i_far) - i_le.min(i_far) <= 1 {
            // Every contour whose first turn back *is* the leading edge. Same
            // start, same steps, same result as before, bit for bit.
            return refine_from_node(i_le);
        }

        // The scan and the furthest node disagree, so at most one of them is the
        // leading edge and neither rule decides it in general. Refine from both
        // and keep whichever is genuinely further from the trailing edge, so the
        // answer is chosen on the quantity itself rather than on either
        // heuristic being right. The unrefined nodes stay eligible too, so a
        // start that the iteration walks away from cannot lose to a worse point.
        let from_scan = refine_from_node(i_le);
        let from_far = refine_from_node(i_far);
        let mut best = from_far;
        let mut best_d = te_distance(from_far);
        for candidate in [from_scan, self.s[i_far], self.s[i_le]] {
            let d = te_distance(candidate);
            if d > best_d {
                best = candidate;
                best_d = d;
            }
        }
        best
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

    fn repo_file(relative: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(relative)
    }

    /// Read a Selig/XFOIL `.dat` file into one coordinate block per element.
    ///
    /// Blank and comment lines separate blocks. Deliberately minimal: the
    /// production import lives in the CLI and the UI, and this only has to open
    /// the fixtures below.
    fn dat_blocks(path: &std::path::Path) -> Vec<Vec<Point>> {
        let text = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let mut blocks: Vec<Vec<Point>> = Vec::new();
        let mut current: Vec<Point> = Vec::new();
        for line in text.lines() {
            let mut parts = line.trim().split_whitespace();
            match (
                parts.next().and_then(|s| s.parse::<f64>().ok()),
                parts.next().and_then(|s| s.parse::<f64>().ok()),
            ) {
                (Some(x), Some(y)) => current.push(point(x, y)),
                _ if !current.is_empty() => blocks.push(std::mem::take(&mut current)),
                _ => {}
            }
        }
        if !current.is_empty() {
            blocks.push(current);
        }
        blocks
    }

    /// Trailing-edge midpoint and the largest node distance from it.
    fn te_and_reach(pts: &[Point]) -> (f64, f64, f64) {
        let n = pts.len();
        let x_te = 0.5 * (pts[0].x + pts[n - 1].x);
        let y_te = 0.5 * (pts[0].y + pts[n - 1].y);
        let reach = pts
            .iter()
            .map(|p| ((p.x - x_te).powi(2) + (p.y - y_te).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);
        (x_te, y_te, reach)
    }

    /// Pairing the leading-edge scan with the furthest node leaves a
    /// conventional section's leading edge exactly where it was.
    ///
    /// The values below were measured before the pairing was introduced. A
    /// conventional airfoil's first turn back towards the trailing edge *is* its
    /// leading edge, so the two candidates agree, the iteration starts in the
    /// same place and lands in the same place, and these are raw bit patterns
    /// rather than tolerances because nothing about the arithmetic changed.
    ///
    /// The same numbers appear in `rustfoil-inviscid`'s
    /// `improved_le_seeding_leaves_single_element_landmarks_bit_identical`, which
    /// measures the same quantity through `AirfoilGeometry`.
    #[test]
    fn lefind_is_bit_identical_on_conventional_sections() {
        // (file, sle, xle, yle)
        let pinned: [(&str, u64, u64, u64); 5] = [
            (
                "naca0012.dat",
                0x3ff050471ff1e073,
                0x0000000000000000,
                0x0000000000000000,
            ),
            (
                "naca2412.dat",
                0x3ff06c2ef11e1d09,
                0xbf1444200ec3adf8,
                0x3f59eb462913dccc,
            ),
            (
                "naca0012_xfoil_paneled.dat",
                0x3ff0505e655529ac,
                0xbe6bad1432c56e00,
                0x0000000000000000,
            ),
            (
                "naca0012_repaneled.dat",
                0x3ff01c98f7f3191b,
                0x0000000000000000,
                0x0000000000000000,
            ),
            (
                "naca0012_buffer_real.dat",
                0x3ff05061b1e2a4fd,
                0x0000000000000000,
                0x0000000000000000,
            ),
        ];

        for (file, sle, xle, yle) in pinned {
            let blocks = dat_blocks(&repo_file(&format!("testdata/{file}")));
            assert_eq!(blocks.len(), 1, "{file}: expected a single element");
            let spline = XfoilSpline::from_points(&blocks[0]).unwrap();
            let s_le = spline.lefind();
            let le = spline.seval(s_le);
            for (name, actual, expected) in [
                ("sle", s_le.to_bits(), sle),
                ("xle", le.x.to_bits(), xle),
                ("yle", le.y.to_bits(), yle),
            ] {
                assert_eq!(
                    actual,
                    expected,
                    "{file}.{name}: {actual:#018x} vs {expected:#018x} ({})",
                    f64::from_bits(actual)
                );
            }
        }
    }

    /// The reported leading edge is the furthest point of the contour from the
    /// trailing edge, on every contour in the bundled coordinate library.
    ///
    /// This is the property the local scan alone does not have: it stops at the
    /// first place the contour turns back towards the trailing edge, which a
    /// cove face aft of the leading edge also satisfies. Measured over the
    /// library, six contours failed this before the furthest node was brought
    /// into the search — the slat and the main element of the 30P-30N section
    /// (in both the combined file and the per-element ones), the truncated main
    /// element `ua79sfm.dat`, and the open cowl curve `naca1.dat` — the worst at
    /// 21% of the true distance.
    ///
    /// The bound is one-sided above: the spline can place the true maximum a
    /// little beyond the furthest *node*, and does, by up to 0.5%.
    ///
    /// Contours carrying a repeated node are left out, because the spline itself
    /// is undefined on them rather than the leading edge being in the wrong
    /// place: a zero-length segment divides by zero in
    /// [`XfoilSpline::from_points`], and XFOIL splits its spline at a doubled
    /// point (SEGSPL) where this port does not. Three contours in the library are
    /// like that — `e337.dat`, `e817.dat` and `fxlv152.dat` — and they behave the
    /// same either side of this search. Only a repeated node earns the skip, so a
    /// contour that comes out unusable for any other reason still fails here.
    #[test]
    fn lefind_reports_the_furthest_point_from_the_trailing_edge() {
        let mut checked = 0usize;
        for directory in ["flexfoil-ui/public/airfoils", "testdata"] {
            // The corpus is not vendored in every checkout; the targeted tests
            // above still cover the behaviour.
            let Ok(entries) = std::fs::read_dir(repo_file(directory)) else {
                continue;
            };
            for entry in entries.flatten() {
                let path = entry.path();
                if path.extension().and_then(|e| e.to_str()) != Some("dat") {
                    continue;
                }
                for (block, pts) in dat_blocks(&path).iter().enumerate() {
                    if pts.len() < 5 {
                        continue;
                    }
                    let repeated_node =
                        pts.windows(2).any(|w| w[0].x == w[1].x && w[0].y == w[1].y);
                    if repeated_node {
                        continue;
                    }
                    let Some(spline) = XfoilSpline::from_points(pts) else {
                        continue;
                    };
                    let (x_te, y_te, reach) = te_and_reach(pts);
                    if reach == 0.0 {
                        continue;
                    }
                    let le = spline.seval(spline.lefind());
                    let distance = ((le.x - x_te).powi(2) + (le.y - y_te).powi(2)).sqrt();
                    assert!(
                        distance >= reach * (1.0 - 1e-6),
                        "{} block {block}: leading edge ({:.6}, {:.6}) is {:.6} from the \
                         trailing edge, against a node {:.6} away",
                        path.display(),
                        le.x,
                        le.y,
                        distance,
                        reach
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no .dat contours were checked");
    }

    /// Every element of a real slat/main/flap section resolves to its own
    /// leading edge.
    ///
    /// `30p-30n.dat` is a McDonnell Douglas 30P-30N section in configuration
    /// coordinates: blocks of 201, 221 and 242 points, the slat and flap
    /// deflected. The slat and the main element both have a cove, and the cove
    /// face is where the local scan stops — the main element's a full 0.66
    /// downstream of its leading edge.
    #[test]
    fn lefind_finds_the_leading_edge_of_every_element_of_a_real_section() {
        let blocks = dat_blocks(&repo_file("flexfoil-ui/public/airfoils/30p-30n.dat"));
        assert_eq!(blocks.len(), 3, "expected slat, main and flap");

        // (role, leading-edge x, leading-edge y)
        let expected = [
            ("slat", -0.0806, -0.1108),
            ("main", 0.0438, -0.0171),
            ("flap", 0.8735, 0.0137),
        ];

        for ((role, xle, yle), pts) in expected.iter().zip(&blocks) {
            let spline = XfoilSpline::from_points(pts).unwrap();
            let le = spline.seval(spline.lefind());
            assert!(
                (le.x - xle).abs() < 5e-4 && (le.y - yle).abs() < 5e-4,
                "{role}: leading edge ({:.6}, {:.6}) against ({xle}, {yle})",
                le.x,
                le.y
            );
        }
    }
}
