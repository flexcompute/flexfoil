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
//!
//! # One body or several
//! [`AirfoilGeometry`] is one body: one trailing edge, one leading edge, one
//! chord. [`ConfigGeometry`] is the multi-element form of the same thing — one
//! flat set of node arrays for the whole configuration plus one
//! [`ElementGeometry`] per element, indexed through a
//! [`Layout`]. It is built out of
//! [`AirfoilGeometry`] rather than beside it, so NCALC, APCALC, TECALC and
//! LEFIND exist once and a single-element `ConfigGeometry` reproduces
//! `AirfoilGeometry` exactly.

use crate::{InviscidError, Result};
use core::ops::Range;
use rustfoil_core::layout::Layout;
use rustfoil_core::paneling::{PanelCounts, PaneledConfiguration};
use rustfoil_core::Configuration;
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
    /// - Panel i goes from node i to node i+1
    /// - Panel n-1 (the "TE panel") goes from node n-1 to node 0
    ///
    /// ## Both trailing-edge nodes are present
    /// `dste` is measured between node 0 and node n-1, so those two have to be
    /// the two trailing-edge nodes for it to be the trailing-edge gap and for
    /// `sharp` to mean what it says. A **sharp** trailing edge therefore arrives
    /// as a closed contour: node n-1 is the closing point, coincident with node 0
    /// (or within rounding of it), `dste` is zero and the TE panel is the
    /// degenerate one XFOIL gives `APANEL = PI`. A **blunt** trailing edge
    /// arrives as an open contour whose two ends are the two trailing-edge
    /// corners.
    ///
    /// A contour with the closing point dropped — node n-1 being the last
    /// lower-surface node *before* the trailing edge — does not satisfy this.
    /// `dste` would then be the length of the final surface panel, which for a
    /// sharp element is the wrong quantity and reports it as blunt. Only
    /// *consecutive* duplicate nodes are rejected, because node 0 and node n-1
    /// coinciding is the sharp case rather than an error.
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

    /// Node furthest from the trailing-edge midpoint, and the arc-length
    /// interval that brackets the leading edge.
    ///
    /// The leading edge is the point of the contour furthest from the
    /// trailing-edge midpoint, so the node where that distance is largest is
    /// within one panel of it and `(s[i-1], s[i+1])` brackets it. One O(n) pass,
    /// and it reads the whole contour, so a straight or vertical stretch
    /// elsewhere on the surface cannot stand in for the leading edge.
    fn leading_edge_bracket(
        x: &[f64],
        y: &[f64],
        s: &[f64],
        xte: f64,
        yte: f64,
    ) -> (usize, f64, f64) {
        let n = x.len();
        let mut i_far = 0usize;
        let mut d_far = f64::NEG_INFINITY;
        for i in 0..n {
            let dx = x[i] - xte;
            let dy = y[i] - yte;
            let d = dx * dx + dy * dy;
            if d > d_far {
                d_far = d;
                i_far = i;
            }
        }
        let lo = s[i_far.saturating_sub(1)];
        let hi = s[(i_far + 1).min(n - 1)];
        (i_far, lo, hi)
    }

    /// Find leading edge location (XFOIL's LEFIND).
    ///
    /// The LE is defined as the point where the surface tangent is
    /// perpendicular to the chord line (TE to LE).
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
    /// So the scan is bracketed by [`leading_edge_bracket`](Self::leading_edge_bracket),
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
    /// `improved_le_seeding_leaves_single_element_landmarks_bit_identical` pins
    /// the unchanged case on the coordinate corpus.
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

        // Node furthest from the trailing edge: a global candidate, independent
        // of where the contour first turns back.
        let (i_far, _s_lo, _s_hi) = Self::leading_edge_bracket(x, y, s, xte, yte);

        // XFOIL's initial guess: find where dot product with TE changes sign
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

        let dseps = (s[n - 1] - s[0]) * 1e-5;

        // XFOIL's Newton iteration for the exact leading-edge arc length,
        // unchanged, run from a given start.
        let refine = |mut s_le: f64| -> f64 {
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
            s_le
        };

        // Distance from the trailing edge, which is the quantity "leading edge"
        // actually names.
        let te_distance = |s_at: f64| -> f64 {
            let (px, py) = Self::seval_point(s_at, s, x, y, xp, yp);
            (px - xte).powi(2) + (py - yte).powi(2)
        };

        let scan_agrees_with_bracket = i_le.max(i_far) - i_le.min(i_far) <= 1;

        let s_le = if scan_agrees_with_bracket {
            // Every contour whose first turn back *is* the leading edge. Same
            // start, same steps, same result as before, bit for bit.
            refine(s[i_le])
        } else {
            // The scan and the bracket disagree, so at most one of them is the
            // leading edge and neither rule decides it in general: the scan
            // stops at the first turn back, which a cove face aft of the
            // leading edge satisfies, while the furthest node from the
            // trailing edge is not the leading edge for every contour either.
            // Refine from both and keep whichever is genuinely further from the
            // trailing edge, so the answer is chosen on the quantity itself
            // rather than on either heuristic being right. Candidates outside
            // the bracket stay eligible, so a correct scan result is never
            // discarded in favour of a worse bracketed one.
            let from_scan = refine(s[i_le]);
            let from_bracket = refine(s[i_far]);
            let mut best = from_bracket;
            let mut best_d = te_distance(from_bracket);
            for candidate in [from_scan, s[i_far], s[i_le]] {
                let d = te_distance(candidate);
                if d > best_d {
                    best = candidate;
                    best_d = d;
                }
            }
            best
        };

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

// ===========================================================================
// Multi-element geometry
// ===========================================================================

/// One element's trailing-edge and leading-edge geometry, and the stretch of a
/// [`ConfigGeometry`]'s node arrays it owns.
///
/// These are exactly the scalars [`AirfoilGeometry`] holds one copy of, plus
/// `start` and `n`. A configuration has one of these per element; nothing here
/// is shared between elements, because a slat, a main element and a flap each
/// have their own trailing edge, leading edge and chord.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementGeometry {
    // === Node range in the global arrays ===
    /// Global index of this element's first node.
    pub start: usize,
    /// Number of nodes belonging to this element.
    pub n: usize,

    // === Trailing edge geometry ===
    /// TE midpoint X coordinate
    pub xte: f64,
    /// TE midpoint Y coordinate
    pub yte: f64,
    /// TE gap length (distance from this element's first node to its last)
    pub dste: f64,
    /// TE normal-projected gap (ANTE in XFOIL)
    pub ante: f64,
    /// TE tangent-projected gap (ASTE in XFOIL)
    pub aste: f64,
    /// True if this element's trailing edge is sharp: `dste < 0.0001 * chord`.
    ///
    /// # This is the authoritative sharpness test
    /// It is the *relative* test the numerics actually gate on — the same
    /// `dste < 0.0001 * chord` as [`AirfoilGeometry::sharp`], consumed by the
    /// trailing-edge branches in `influence.rs` and `system.rs`. Anything that
    /// decides how a trailing edge is treated aerodynamically reads this.
    ///
    /// [`ElementSpan::has_blunt_te`](rustfoil_core::layout::ElementSpan::has_blunt_te)
    /// is not a substitute. That is a *connectivity* question — are this
    /// element's upper and lower trailing-edge nodes two distinct nodes? —
    /// answered by an absolute componentwise `1e-10` closure test in
    /// rustfoil-core. The two disagree in both directions: a contour whose ends
    /// sit `1e-6` apart on a unit chord has distinct nodes and is sharp here,
    /// and a contour on a chord below the closure tolerance can be one node and
    /// blunt here. Substituting one for the other changes which trailing-edge
    /// treatment the solver applies.
    pub sharp: bool,

    // === Leading edge geometry ===
    /// LE X coordinate
    pub xle: f64,
    /// LE Y coordinate
    pub yle: f64,
    /// Arc length at the leading edge, measured from this element's first node.
    pub sle: f64,

    // === Dimensions ===
    /// Chord length, from this element's leading edge to its TE midpoint.
    ///
    /// # This is the authoritative chord
    /// Per decision D4 this element's own coefficients normalise by this, and a
    /// configuration total normalises by
    /// [`default_ref_chord`](ConfigGeometry::default_ref_chord), which is the
    /// largest of these. Both ends of the measurement are the solver's own:
    /// [`xle`](Self::xle)/[`yle`](Self::yle) from LEFIND on this element's
    /// spline, and [`xte`](Self::xte)/[`yte`](Self::yte) from TECALC.
    ///
    /// [`Element::chord`](rustfoil_core::Element::chord) is a different length
    /// and is not a substitute. It measures from the contour's minimum-x node,
    /// which is a coarser stand-in for the leading edge on a paneled contour and
    /// not the leading edge at all for an element whose deflection is baked into
    /// its coordinates rather than carried in its
    /// [`Placement`](rustfoil_core::Placement). Measured on this repository's
    /// fixtures the two differ by 7.9e-5 relative on `naca2412.dat` and 2.6e-5
    /// on `naca0012_xfoil_paneled.dat`, and by 3.3% on the deflected slat and
    /// 0.8% on the deflected flap of the 30P-30N section.
    /// `the_two_chord_definitions_disagree_by_a_pinned_amount` measures it, and
    /// [`Element::chord`](rustfoil_core::Element::chord) documents what it is
    /// for.
    pub chord: f64,
}

impl ElementGeometry {
    /// Global index one past this element's last node.
    #[inline]
    pub fn end(&self) -> usize {
        self.start + self.n
    }

    /// This element's node range in the global arrays.
    #[inline]
    pub fn range(&self) -> Range<usize> {
        self.start..self.end()
    }

    /// Whether `global` is one of this element's nodes.
    #[inline]
    pub fn contains(&self, global: usize) -> bool {
        global >= self.start && global < self.end()
    }

    /// This element's TE panel source/vortex coefficients (SCS, SDS).
    ///
    /// The single-element [`AirfoilGeometry::te_coefficients`] for this
    /// element's own trailing edge: `(1.0, 0.0)` when [`sharp`](Self::sharp),
    /// otherwise `(ante / dste, aste / dste)`.
    pub fn te_coefficients(&self) -> (f64, f64) {
        if self.sharp {
            (1.0, 0.0)
        } else {
            (self.ante / self.dste, self.aste / self.dste)
        }
    }
}

/// A configuration's inviscid geometry: one flat set of node arrays for every
/// element, plus one [`ElementGeometry`] per element.
///
/// # Layout of the node arrays
/// `x`, `y`, `s`, `xp`, `yp`, `nx`, `ny` and `apanel` are each
/// [`total_nodes`](Self::total_nodes) long and hold every element's nodes
/// concatenated in configuration order — one contiguous array plus a segment
/// table, not a per-element fragmentation, because that is what the influence
/// kernels want. Element `k` owns `elements()[k].range()`, and
/// [`layout`](Self::layout) is the same table in
/// [`Layout`] form.
///
/// # Per element, not across elements
/// Everything derived is derived per element and never across a boundary:
///
/// - `s` **restarts at 0.0 at each element's first node**. It is that element's
///   own arc length, not a running total over the configuration. A cumulative
///   `s` would put a spurious segment between one element's last node and the
///   next element's first.
/// - `xp`, `yp` come from one SEGSPL over that element's nodes alone, so an
///   element's spline end conditions are its own trailing edge and not its
///   neighbour's leading edge.
/// - `nx`, `ny` corner-average only within an element.
/// - `apanel[i]` is the angle of the panel from node `i` to
///   [`next_node(i)`](Self::next_node), which wraps inside the owning element.
///   At an element's last node that is the element's TE panel, closing onto the
///   element's own first node — the same slot `apanel[n-1]` holds for a single
///   airfoil.
///
/// # What counts as a node
/// The same convention as the single-element path, node for node: an element's
/// node 0 is its upper-surface trailing edge and its last node its
/// lower-surface trailing edge, and **both** trailing-edge nodes are present.
/// For a blunt trailing edge those are the two trailing-edge corners; for a
/// sharp one they coincide, the last node being the contour's closing point.
/// See [`AirfoilGeometry::from_points`] for why: `dste`, and so
/// [`sharp`](ElementGeometry::sharp) and the trailing-edge treatment that gates
/// on it, are measured between an element's first and last node.
///
/// [`PaneledConfiguration`] keeps only an element's *distinct* nodes and omits a
/// closed contour's closing point, so it is one node short of this convention
/// per sharp-trailing-edge element.
/// [`from_paneled`](Self::from_paneled) restores it. The two numberings are
/// therefore not interchangeable for a configuration with a sharp trailing edge:
/// index this geometry with [`layout`](Self::layout) or with its
/// [`ElementGeometry`] ranges, not with a [`Layout`] built from the same
/// [`Configuration`] elsewhere.
///
/// # Landmarks live in `ElementGeometry`
/// The [`Layout`] here carries **connectivity
/// only** — node counts and the closure rule. Its spans deliberately make no
/// landmark claim, so `layout().span(k).has_le()` is `false` and
/// `has_blunt_te()` is `false` for every element however the geometry was
/// built. Ask [`ElementGeometry`] for a leading edge (`xle`, `yle`, `sle`), a
/// trailing edge (`xte`, `yte`, `dste`) or a sharpness verdict (`sharp`); those
/// are computed from the geometry and are the only answers.
///
/// # Fields are private
/// The arrays and the element table have to agree with each other — each array
/// is `layout.total_nodes()` long, and `elements[k]` matches `layout.span(k)`.
/// Construction is what establishes that, so it is also the only way in. Slices
/// read out through the accessors are as cheap as a field read.
///
/// # Example
/// ```
/// use rustfoil_core::naca::naca4;
/// use rustfoil_inviscid::geometry::{AirfoilGeometry, ConfigGeometry};
///
/// let points: Vec<(f64, f64)> = naca4(2412, Some(80)).iter().map(|p| (p.x, p.y)).collect();
///
/// let single = AirfoilGeometry::from_points(&points).unwrap();
/// let config = ConfigGeometry::from_points(&points).unwrap();
///
/// // One element, and the same derived geometry.
/// assert_eq!(config.n_elements(), 1);
/// assert_eq!(config.total_nodes(), single.n);
/// assert_eq!(config.element(0).chord.to_bits(), single.chord.to_bits());
/// assert_eq!(config.element(0).sharp, single.sharp);
///
/// // A single element still closes onto itself.
/// assert_eq!(config.next_node(config.total_nodes() - 1), 0);
/// ```
#[derive(Debug, Clone)]
pub struct ConfigGeometry {
    x: Vec<f64>,
    y: Vec<f64>,
    s: Vec<f64>,
    xp: Vec<f64>,
    yp: Vec<f64>,
    nx: Vec<f64>,
    ny: Vec<f64>,
    apanel: Vec<f64>,
    elements: Vec<ElementGeometry>,
    layout: Layout,
}

impl ConfigGeometry {
    /// Build a one-element configuration geometry from one body's coordinates.
    ///
    /// Equivalent to [`AirfoilGeometry::from_points`] on the same points: the
    /// node arrays and the element's scalars are produced by that call and are
    /// bit-for-bit the same values.
    ///
    /// # Errors
    /// Whatever [`AirfoilGeometry::from_points`] reports.
    pub fn from_points(points: &[(f64, f64)]) -> Result<Self> {
        Self::from_element_points(&[points])
    }

    /// Build from one coordinate list per element, in configuration order.
    ///
    /// The coordinates are taken as given, in configuration coordinates — any
    /// [`Placement`](rustfoil_core::Placement) has to be applied before this
    /// point. Each element must satisfy [`AirfoilGeometry::from_points`] on its
    /// own, so at least ten nodes each and no duplicate consecutive nodes.
    ///
    /// Accepts anything that borrows as a coordinate slice, so `&[Vec<_>]` and
    /// `&[&[_]]` both work.
    ///
    /// # A closed element keeps its closing point
    /// An element whose first and last coordinates coincide is a sharp
    /// trailing edge, and that closing point is its last node: it is neither
    /// stripped nor rejected here. That is the node convention on this type, and
    /// it is what makes the element's `dste` its trailing-edge gap rather than
    /// the length of its final surface panel. So the blocks of a multi-element
    /// coordinate file transfer straight in, closed blocks and open blocks alike,
    /// and a closed one comes out `sharp` with `dste == 0.0` and a degenerate
    /// trailing-edge panel at `APANEL = PI`, exactly as the same block would
    /// through [`AirfoilGeometry::from_points`].
    ///
    /// Consecutive duplicates elsewhere on a contour remain an error; only the
    /// first-to-last coincidence is meaningful.
    ///
    /// # Errors
    /// Whatever [`AirfoilGeometry::from_points`] reports for the first element
    /// that fails. Indices inside that error are local to the failing element,
    /// and construction stops there.
    pub fn from_element_points<P: AsRef<[(f64, f64)]>>(elements: &[P]) -> Result<Self> {
        let mut per_element = Vec::with_capacity(elements.len());
        for points in elements {
            per_element.push(AirfoilGeometry::from_points(points.as_ref())?);
        }
        Ok(Self::from_element_geometries(&per_element))
    }

    /// Build from an already-paneled configuration.
    ///
    /// [`PaneledConfiguration`] holds each element's nodes in configuration
    /// coordinates, which transfer unchanged, with one adjustment: it keeps only
    /// an element's *distinct* nodes, so a sharp-trailing-edge element's closing
    /// point — coincident with its first node — is not among them, and this
    /// type's node convention needs it back. It is restored here, for the
    /// elements
    /// [`element_is_closed`](rustfoil_core::paneling::PaneledConfiguration::element_is_closed)
    /// reports closed, by repeating that element's first node.
    ///
    /// Without it an element's `dste` would be the distance from its first node
    /// to the last lower-surface node before its trailing edge — the length of
    /// its final surface panel — and a sharp element would come out
    /// [`sharp`](ElementGeometry::sharp)`== false`, taking the blunt
    /// trailing-edge treatment. So a closed element here has one more node than
    /// the same element in the `PaneledConfiguration` it came from, and an open
    /// one has the same number.
    ///
    /// # Errors
    /// As [`from_element_points`](Self::from_element_points). An element paneled
    /// to fewer than ten nodes is the likely one.
    pub fn from_paneled(paneled: &PaneledConfiguration) -> Result<Self> {
        let mut contours: Vec<Vec<(f64, f64)>> = Vec::with_capacity(paneled.n_elements());
        for element in 0..paneled.n_elements() {
            let mut nodes: Vec<(f64, f64)> = paneled
                .element_nodes(element)
                .iter()
                .map(|p| (p.x, p.y))
                .collect();
            if paneled.element_is_closed(element) {
                if let Some(&first) = nodes.first() {
                    nodes.push(first);
                }
            }
            contours.push(nodes);
        }
        Self::from_element_points(&contours)
    }

    /// Panel a configuration and build its geometry in one step.
    ///
    /// [`Configuration::panel_all`](rustfoil_core::Configuration::panel_all)
    /// followed by [`from_paneled`](Self::from_paneled).
    ///
    /// # Errors
    /// - `SplineError` wrapping the paneling failure if the configuration
    ///   cannot be paneled at the requested counts.
    /// - Otherwise as [`from_element_points`](Self::from_element_points).
    pub fn from_configuration(config: &Configuration, counts: &PanelCounts) -> Result<Self> {
        let paneled = config.panel_all(counts).map_err(|error| {
            InviscidError::SplineError(format!("paneling the configuration: {error}"))
        })?;
        Self::from_paneled(&paneled)
    }

    /// Concatenate per-element geometries into the flat form.
    fn from_element_geometries(per_element: &[AirfoilGeometry]) -> Self {
        let counts: Vec<usize> = per_element.iter().map(|geom| geom.n).collect();
        // Connectivity only, so no landmark claim is made; see the type's
        // documentation. `from_node_counts` rejects only a zero-node element,
        // and `AirfoilGeometry::from_points` has already refused anything under
        // ten nodes, so this cannot fail.
        let layout = Layout::from_node_counts(&counts)
            .expect("every element geometry has at least ten nodes");
        let total = layout.total_nodes();

        let mut x = Vec::with_capacity(total);
        let mut y = Vec::with_capacity(total);
        let mut s = Vec::with_capacity(total);
        let mut xp = Vec::with_capacity(total);
        let mut yp = Vec::with_capacity(total);
        let mut nx = Vec::with_capacity(total);
        let mut ny = Vec::with_capacity(total);
        let mut apanel = Vec::with_capacity(total);
        let mut elements = Vec::with_capacity(per_element.len());

        let mut start = 0usize;
        for geom in per_element {
            x.extend_from_slice(&geom.x);
            y.extend_from_slice(&geom.y);
            s.extend_from_slice(&geom.s);
            xp.extend_from_slice(&geom.xp);
            yp.extend_from_slice(&geom.yp);
            nx.extend_from_slice(&geom.nx);
            ny.extend_from_slice(&geom.ny);
            apanel.extend_from_slice(&geom.apanel);

            elements.push(ElementGeometry {
                start,
                n: geom.n,
                xte: geom.xte,
                yte: geom.yte,
                dste: geom.dste,
                ante: geom.ante,
                aste: geom.aste,
                sharp: geom.sharp,
                xle: geom.xle,
                yle: geom.yle,
                sle: geom.sle,
                chord: geom.chord,
            });
            start += geom.n;
        }

        Self {
            x,
            y,
            s,
            xp,
            yp,
            nx,
            ny,
            apanel,
            elements,
            layout,
        }
    }

    // --- dimensions and tables -------------------------------------------

    /// Number of elements.
    #[inline]
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Total number of nodes across all elements — the length of every node
    /// array.
    #[inline]
    pub fn total_nodes(&self) -> usize {
        self.layout.total_nodes()
    }

    /// True if there are no elements, and therefore no nodes.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// The node numbering, one span per element.
    ///
    /// Connectivity only — see the type's documentation before reading a
    /// landmark off a span.
    #[inline]
    pub fn layout(&self) -> &Layout {
        &self.layout
    }

    /// Every element's geometry, in configuration order.
    #[inline]
    pub fn elements(&self) -> &[ElementGeometry] {
        &self.elements
    }

    /// One element's geometry.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element(&self, element: usize) -> &ElementGeometry {
        &self.elements[element]
    }

    /// One element's geometry, or `None` if out of range.
    #[inline]
    pub fn try_element(&self, element: usize) -> Option<&ElementGeometry> {
        self.elements.get(element)
    }

    /// Index of the largest-chord element — the main element, by the D4 rule.
    ///
    /// Compared on [`ElementGeometry::chord`], so on the same chords the
    /// per-element coefficients normalise by. Ties go to the lowest index;
    /// `None` for an empty geometry.
    ///
    /// [`Configuration::main_element_index`](rustfoil_core::Configuration::main_element_index)
    /// answers the same question before there is any inviscid geometry, from
    /// [`Element::chord`](rustfoil_core::Element::chord). The two can only
    /// disagree over elements whose chords are closer together than the
    /// difference between the two definitions — a few parts in 1e5 for elements
    /// in their own coordinates, a few percent for one whose deflection is baked
    /// into its coordinates. See [`ElementGeometry::chord`].
    pub fn main_element_index(&self) -> Option<usize> {
        self.elements
            .iter()
            .enumerate()
            .fold(None, |best, (i, element)| match best {
                // Strictly greater, so the first of equal chords wins.
                Some((_, best_chord)) if element.chord <= best_chord => best,
                _ => Some((i, element.chord)),
            })
            .map(|(i, _)| i)
    }

    /// The D4 default reference chord for configuration totals: the chord of the
    /// largest-chord element.
    ///
    /// # Which reference chord a total should use
    /// This is the authoritative form of the D4 default, because it is measured
    /// between the landmarks the solver uses for its own geometry — see
    /// [`ElementGeometry::chord`].
    /// [`Configuration::resolved_ref_chord`](rustfoil_core::Configuration::resolved_ref_chord)
    /// resolves the same default at geometry time, from
    /// [`Element::chord`](rustfoil_core::Element::chord), and is what the
    /// geometric clearance diagnostics scale by; it differs from this by the
    /// amount recorded on [`ElementGeometry::chord`].
    ///
    /// This is only the *default*. An explicit
    /// [`Configuration::ref_chord`](rustfoil_core::Configuration::ref_chord)
    /// overrides it and is not visible from here, so a caller that sets one has
    /// to carry it through itself. A high-lift rigging table normally quotes
    /// against the retracted chord of the whole section rather than the main
    /// element alone, which is exactly that case.
    ///
    /// Returns `1.0` for an empty geometry, or one whose main element has no
    /// measurable chord, so a caller dividing by this never divides by zero.
    pub fn default_ref_chord(&self) -> f64 {
        let chord = self
            .main_element_index()
            .map(|i| self.elements[i].chord)
            .unwrap_or(0.0);
        if chord > 0.0 {
            chord
        } else {
            1.0
        }
    }

    // --- node arrays -----------------------------------------------------

    /// X coordinates at every node.
    #[inline]
    pub fn x(&self) -> &[f64] {
        &self.x
    }

    /// Y coordinates at every node.
    #[inline]
    pub fn y(&self) -> &[f64] {
        &self.y
    }

    /// Arc length at every node, restarting at 0.0 at each element.
    #[inline]
    pub fn s(&self) -> &[f64] {
        &self.s
    }

    /// dX/dS at every node, from that element's own spline.
    #[inline]
    pub fn xp(&self) -> &[f64] {
        &self.xp
    }

    /// dY/dS at every node, from that element's own spline.
    #[inline]
    pub fn yp(&self) -> &[f64] {
        &self.yp
    }

    /// Outward unit normal X component at every node.
    #[inline]
    pub fn nx(&self) -> &[f64] {
        &self.nx
    }

    /// Outward unit normal Y component at every node.
    #[inline]
    pub fn ny(&self) -> &[f64] {
        &self.ny
    }

    /// Panel angle for the panel starting at every node — see the type's
    /// documentation for what the panel at an element's last node is.
    #[inline]
    pub fn apanel(&self) -> &[f64] {
        &self.apanel
    }

    // --- per-element views -----------------------------------------------

    /// One element's node range in the global arrays.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_range(&self, element: usize) -> Range<usize> {
        self.elements[element].range()
    }

    /// One element's X coordinates.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_x(&self, element: usize) -> &[f64] {
        &self.x[self.element_range(element)]
    }

    /// One element's Y coordinates.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_y(&self, element: usize) -> &[f64] {
        &self.y[self.element_range(element)]
    }

    /// One element's arc lengths, starting at 0.0.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_s(&self, element: usize) -> &[f64] {
        &self.s[self.element_range(element)]
    }

    /// One element's dX/dS.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_xp(&self, element: usize) -> &[f64] {
        &self.xp[self.element_range(element)]
    }

    /// One element's dY/dS.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_yp(&self, element: usize) -> &[f64] {
        &self.yp[self.element_range(element)]
    }

    /// One element's outward normal X components.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_nx(&self, element: usize) -> &[f64] {
        &self.nx[self.element_range(element)]
    }

    /// One element's outward normal Y components.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_ny(&self, element: usize) -> &[f64] {
        &self.ny[self.element_range(element)]
    }

    /// One element's panel angles, its TE panel last.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    #[inline]
    pub fn element_apanel(&self, element: usize) -> &[f64] {
        &self.apanel[self.element_range(element)]
    }

    /// One element's total arc length.
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    pub fn element_arc_length(&self, element: usize) -> f64 {
        self.element_s(element).last().copied().unwrap_or(0.0)
    }

    // --- connectivity ----------------------------------------------------

    /// Which element owns a global node index.
    ///
    /// # Panics
    /// If `global >= self.total_nodes()`.
    #[inline]
    pub fn element_of(&self, global: usize) -> usize {
        self.layout.element_of(global)
    }

    /// The node that closes the panel starting at `global` — the next node
    /// within the same element.
    ///
    /// The sanctioned replacement for `(global + 1) % n`. For a single-element
    /// geometry the two are the same value.
    ///
    /// # Panics
    /// If `global >= self.total_nodes()`.
    #[inline]
    pub fn next_node(&self, global: usize) -> usize {
        self.layout.next_node(global)
    }

    /// The node before `global` within the same element.
    ///
    /// # Panics
    /// If `global >= self.total_nodes()`.
    #[inline]
    pub fn prev_node(&self, global: usize) -> usize {
        self.layout.prev_node(global)
    }

    // --- trailing-edge quantities ----------------------------------------

    /// One element's TE panel source/vortex coefficients (SCS, SDS).
    ///
    /// # Panics
    /// If `element >= self.n_elements()`.
    pub fn te_coefficients(&self, element: usize) -> (f64, f64) {
        self.elements[element].te_coefficients()
    }

    /// Every element's TE panel coefficients, in configuration order.
    pub fn te_coefficients_all(&self) -> Vec<(f64, f64)> {
        self.elements
            .iter()
            .map(ElementGeometry::te_coefficients)
            .collect()
    }

    /// One element's sharp trailing-edge bisector control point and normal.
    ///
    /// The per-element form of [`AirfoilGeometry::sharp_te_bisector_control`],
    /// using this element's own first and last nodes and its own TE midpoint.
    /// `None` when the element's trailing edge is not
    /// [`sharp`](ElementGeometry::sharp), when it has fewer than three nodes, or
    /// when `element` is out of range.
    pub fn sharp_te_bisector_control(&self, element: usize) -> Option<(f64, f64, f64, f64)> {
        let geom = self.try_element(element)?;
        if !geom.sharp || geom.n < 3 {
            return None;
        }
        // The element's own first and last node, standing in for node 0 and
        // node n-1 of a single airfoil.
        let first = geom.start;
        let last = geom.end() - 1;

        let upper_tx = -self.xp[first];
        let upper_ty = -self.yp[first];
        let lower_tx = self.xp[last];
        let lower_ty = self.yp[last];

        let bis_x = upper_tx + lower_tx;
        let bis_y = upper_ty + lower_ty;
        let bis_norm = (bis_x * bis_x + bis_y * bis_y).sqrt().max(1.0e-12);
        let cbis = bis_x / bis_norm;
        let sbis = bis_y / bis_norm;

        let ds1 = ((self.x[first] - self.x[first + 1]).powi(2)
            + (self.y[first] - self.y[first + 1]).powi(2))
        .sqrt();
        let ds2 = ((self.x[last] - self.x[last - 1]).powi(2)
            + (self.y[last] - self.y[last - 1]).powi(2))
        .sqrt();
        let dsmin = ds1.min(ds2);
        let bwt = 0.1;

        let xbis = geom.xte - bwt * dsmin * cbis;
        let ybis = geom.yte - bwt * dsmin * sbis;

        Some((xbis, ybis, -sbis, cbis))
    }

    // --- compatibility bridge --------------------------------------------

    /// One element as a standalone [`AirfoilGeometry`].
    ///
    /// A bridge for reusing the single-element kernels one element at a time:
    /// the returned value is exactly what
    /// [`AirfoilGeometry::from_points`] would have produced from that element's
    /// coordinates alone. Its node indices are element-local, so a global index
    /// from this configuration does not index it — subtract
    /// `element(k).start`, or use [`Layout::from_global`].
    ///
    /// Copies the element's arrays, so this is a setup-time convenience rather
    /// than something to call from a loop. `None` if `element` is out of range.
    pub fn element_as_airfoil_geometry(&self, element: usize) -> Option<AirfoilGeometry> {
        let geom = self.try_element(element)?;
        let range = geom.range();
        Some(AirfoilGeometry {
            x: self.x[range.clone()].to_vec(),
            y: self.y[range.clone()].to_vec(),
            s: self.s[range.clone()].to_vec(),
            xp: self.xp[range.clone()].to_vec(),
            yp: self.yp[range.clone()].to_vec(),
            nx: self.nx[range.clone()].to_vec(),
            ny: self.ny[range.clone()].to_vec(),
            apanel: self.apanel[range].to_vec(),
            xte: geom.xte,
            yte: geom.yte,
            dste: geom.dste,
            ante: geom.ante,
            aste: geom.aste,
            sharp: geom.sharp,
            xle: geom.xle,
            yle: geom.yle,
            sle: geom.sle,
            n: geom.n,
            chord: geom.chord,
        })
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

    // =====================================================================
    // ConfigGeometry
    // =====================================================================

    use std::path::PathBuf;

    /// A NACA 0012 whose trailing edge has been opened to `gap`, split evenly
    /// about `y = 0`.
    fn make_naca0012_with_te_gap(n_panels: usize, gap: f64) -> Vec<(f64, f64)> {
        let mut points = make_naca0012(n_panels);
        let last = points.len() - 1;
        points[0].1 += 0.5 * gap;
        points[last].1 -= 0.5 * gap;
        points
    }

    fn testdata(file: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata")
            .join(file)
    }

    /// A path relative to the repository root.
    fn repo_file(relative: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(relative)
    }

    /// Read a Selig/XFOIL `.dat` file into one coordinate block per element.
    ///
    /// Blank lines and comment lines separate blocks. Deliberately minimal: the
    /// production import lives in the CLI and the UI, and this only has to open
    /// the fixtures below.
    fn read_dat_blocks(file: &str) -> Vec<Vec<(f64, f64)>> {
        let path = testdata(file);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));

        let mut blocks: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut current: Vec<(f64, f64)> = Vec::new();
        for line in text.lines() {
            let trimmed = line.trim();
            let mut parts = trimmed.split_whitespace();
            let parsed = match (
                parts.next().and_then(|s| s.parse::<f64>().ok()),
                parts.next().and_then(|s| s.parse::<f64>().ok()),
            ) {
                (Some(x), Some(y)) => Some((x, y)),
                _ => None,
            };
            match parsed {
                Some(pair) => current.push(pair),
                None if !current.is_empty() => blocks.push(std::mem::take(&mut current)),
                None => {}
            }
        }
        if !current.is_empty() {
            blocks.push(current);
        }
        blocks
    }

    /// Every element of a real slat/main/flap fixture.
    fn mda_blocks() -> Vec<Vec<(f64, f64)>> {
        let blocks = read_dat_blocks("mda_30p_30n_trimmed.dat");
        assert_eq!(blocks.len(), 3, "expected slat, main and flap");
        blocks
    }

    /// Compare two f64 arrays bit for bit.
    fn assert_bits_eq(actual: &[f64], expected: &[f64], what: &str) {
        assert_eq!(actual.len(), expected.len(), "{what}: length");
        for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "{what}[{i}]: {a:?} ({:#x}) vs {b:?} ({:#x})",
                a.to_bits(),
                b.to_bits()
            );
        }
    }

    /// Every derived quantity of a one-element `ConfigGeometry` against the
    /// `AirfoilGeometry` for the same points, compared on the raw bits.
    fn assert_single_element_bits_identical(
        config: &ConfigGeometry,
        single: &AirfoilGeometry,
        what: &str,
    ) {
        assert_eq!(config.n_elements(), 1, "{what}: element count");
        assert_eq!(config.total_nodes(), single.n, "{what}: total nodes");

        let element = config.element(0);
        assert_eq!(element.start, 0, "{what}: start");
        assert_eq!(element.n, single.n, "{what}: n");
        assert_eq!(element.range(), 0..single.n, "{what}: range");

        // Node arrays.
        assert_bits_eq(config.x(), &single.x, &format!("{what}.x"));
        assert_bits_eq(config.y(), &single.y, &format!("{what}.y"));
        assert_bits_eq(config.s(), &single.s, &format!("{what}.s"));
        assert_bits_eq(config.xp(), &single.xp, &format!("{what}.xp"));
        assert_bits_eq(config.yp(), &single.yp, &format!("{what}.yp"));
        assert_bits_eq(config.nx(), &single.nx, &format!("{what}.nx"));
        assert_bits_eq(config.ny(), &single.ny, &format!("{what}.ny"));
        assert_bits_eq(config.apanel(), &single.apanel, &format!("{what}.apanel"));

        // The per-element slices are the same arrays for one element.
        assert_bits_eq(config.element_x(0), &single.x, &format!("{what}.element_x"));
        assert_bits_eq(config.element_s(0), &single.s, &format!("{what}.element_s"));
        assert_eq!(
            config.element_arc_length(0).to_bits(),
            single.total_arc_length().to_bits(),
            "{what}: arc length"
        );

        // Trailing-edge, leading-edge and chord scalars.
        for (name, actual, expected) in [
            ("xte", element.xte, single.xte),
            ("yte", element.yte, single.yte),
            ("dste", element.dste, single.dste),
            ("ante", element.ante, single.ante),
            ("aste", element.aste, single.aste),
            ("xle", element.xle, single.xle),
            ("yle", element.yle, single.yle),
            ("sle", element.sle, single.sle),
            ("chord", element.chord, single.chord),
        ] {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "{what}.{name}: {actual:?} vs {expected:?}"
            );
        }
        assert_eq!(element.sharp, single.sharp, "{what}.sharp");

        // Derived trailing-edge quantities.
        let (scs, sds) = config.te_coefficients(0);
        let (scs_single, sds_single) = single.te_coefficients();
        assert_eq!(scs.to_bits(), scs_single.to_bits(), "{what}: scs");
        assert_eq!(sds.to_bits(), sds_single.to_bits(), "{what}: sds");
        assert_eq!(
            config.te_coefficients_all(),
            vec![(scs_single, sds_single)],
            "{what}: te_coefficients_all"
        );

        match (
            config.sharp_te_bisector_control(0),
            single.sharp_te_bisector_control(),
        ) {
            (None, None) => {}
            (Some(a), Some(b)) => {
                assert_eq!(a.0.to_bits(), b.0.to_bits(), "{what}: bisector x");
                assert_eq!(a.1.to_bits(), b.1.to_bits(), "{what}: bisector y");
                assert_eq!(a.2.to_bits(), b.2.to_bits(), "{what}: bisector nx");
                assert_eq!(a.3.to_bits(), b.3.to_bits(), "{what}: bisector ny");
            }
            (a, b) => panic!("{what}: bisector control disagrees: {a:?} vs {b:?}"),
        }

        // And the bridge back to the single-element type.
        let bridged = config
            .element_as_airfoil_geometry(0)
            .expect("element 0 exists");
        assert_bits_eq(&bridged.x, &single.x, &format!("{what}.bridged.x"));
        assert_bits_eq(&bridged.y, &single.y, &format!("{what}.bridged.y"));
        assert_bits_eq(&bridged.s, &single.s, &format!("{what}.bridged.s"));
        assert_bits_eq(&bridged.xp, &single.xp, &format!("{what}.bridged.xp"));
        assert_bits_eq(&bridged.yp, &single.yp, &format!("{what}.bridged.yp"));
        assert_bits_eq(&bridged.nx, &single.nx, &format!("{what}.bridged.nx"));
        assert_bits_eq(&bridged.ny, &single.ny, &format!("{what}.bridged.ny"));
        assert_bits_eq(
            &bridged.apanel,
            &single.apanel,
            &format!("{what}.bridged.apanel"),
        );
        assert_eq!(bridged.n, single.n, "{what}.bridged.n");
        assert_eq!(bridged.sharp, single.sharp, "{what}.bridged.sharp");
        assert_eq!(
            bridged.chord.to_bits(),
            single.chord.to_bits(),
            "{what}.bridged.chord"
        );
        assert_eq!(
            bridged.sle.to_bits(),
            single.sle.to_bits(),
            "{what}.bridged.sle"
        );
    }

    /// The constraint the whole refactor rests on: with one element, every
    /// number a `ConfigGeometry` derives is the number `AirfoilGeometry` derives
    /// today, to the last bit. Checked on two real coordinate files and on the
    /// analytic section, so it covers real paneling, a blunt trailing edge and a
    /// sharp one.
    #[test]
    fn one_element_config_geometry_is_bit_identical_to_airfoil_geometry() {
        let mut cases: Vec<(String, Vec<(f64, f64)>)> = vec![
            ("analytic naca0012".to_string(), make_naca0012(80)),
            (
                "analytic naca0012, 1e-6 te gap".to_string(),
                make_naca0012_with_te_gap(80, 1e-6),
            ),
        ];
        for file in ["naca0012.dat", "naca2412.dat", "naca0012_xfoil_paneled.dat"] {
            let blocks = read_dat_blocks(file);
            assert_eq!(blocks.len(), 1, "{file}: expected a single element");
            cases.push((file.to_string(), blocks.into_iter().next().unwrap()));
        }

        for (what, points) in &cases {
            let single = AirfoilGeometry::from_points(points)
                .unwrap_or_else(|e| panic!("{what}: AirfoilGeometry::from_points: {e}"));
            let config = ConfigGeometry::from_points(points)
                .unwrap_or_else(|e| panic!("{what}: ConfigGeometry::from_points: {e}"));
            assert_single_element_bits_identical(&config, &single, what);

            // The one-element form of the list constructor has to agree too.
            let from_list = ConfigGeometry::from_element_points(&[points.clone()]).unwrap();
            assert_single_element_bits_identical(&from_list, &single, what);
        }

        // Both trailing-edge branches were actually exercised.
        let sharp = ConfigGeometry::from_points(&make_naca0012(80)).unwrap();
        let blunt = ConfigGeometry::from_points(&read_dat_blocks("naca0012.dat")[0]).unwrap();
        assert!(sharp.element(0).sharp, "analytic section should be sharp");
        assert!(!blunt.element(0).sharp, "naca0012.dat should be blunt");
        assert!(sharp.sharp_te_bisector_control(0).is_some());
        assert!(blunt.sharp_te_bisector_control(0).is_none());
    }

    /// Sharpness is the solver's relative test, and it can disagree with the
    /// node-identity question rustfoil-core's closure test answers.
    #[test]
    fn sharpness_is_the_relative_test_not_node_identity() {
        // A unit chord with a 1e-6 trailing-edge gap. The gap is 10^4 times the
        // 1e-10 closure tolerance, so the two trailing-edge nodes are distinct
        // and `ElementSpan::has_blunt_te` would report a blunt trailing edge.
        // The relative test is what the numerics use, and 1e-6 < 1e-4 * chord.
        let gap = 1e-6;
        let points = make_naca0012_with_te_gap(80, gap);
        let config = ConfigGeometry::from_points(&points).unwrap();
        let element = config.element(0);

        assert!(element.dste > rustfoil_core::CONTOUR_CLOSURE_TOLERANCE);
        assert!(element.dste < 1e-4 * element.chord);
        assert!(
            element.sharp,
            "dste = {} on chord {} is sharp by the relative test",
            element.dste, element.chord
        );

        // rustfoil-core's connectivity test disagrees, which is the point.
        let contour: Vec<rustfoil_core::Point> = points
            .iter()
            .map(|&(x, y)| rustfoil_core::point(x, y))
            .collect();
        assert!(
            !rustfoil_core::contour_is_closed(&contour),
            "the contour's ends are further apart than the closure tolerance"
        );

        // And the sharp branch is the one taken.
        assert_eq!(config.te_coefficients(0), (1.0, 0.0));
        assert!(config.sharp_te_bisector_control(0).is_some());
    }

    /// A `ConfigGeometry`'s layout answers connectivity questions only. Nobody
    /// should be reading a leading edge or a sharpness verdict off it.
    #[test]
    fn config_geometry_layout_claims_no_landmarks() {
        let config = ConfigGeometry::from_element_points(&mda_blocks()).unwrap();
        assert_eq!(config.n_elements(), 3);

        for element in 0..config.n_elements() {
            let span = config.layout().span(element);
            // Read from another crate, through the accessors and through the
            // fields, so both forms stay usable outside rustfoil-core.
            assert_eq!(span.start(), config.element(element).start);
            assert_eq!(span.len(), config.element(element).n);
            assert_eq!(span.len, span.len());
            assert!(!span.has_le(), "element {element} claims a leading edge");
            assert!(
                !span.has_blunt_te(),
                "element {element} claims a trailing-edge node count"
            );
            // The geometry does carry a leading edge, and it is here.
            let geom = config.element(element);
            assert!(geom.sle > 0.0 && geom.sle < config.element_arc_length(element));
        }
    }

    /// Three real elements: each one's arrays are exactly what that element
    /// alone produces, so nothing is splined, normalised or arc-lengthed across
    /// an element boundary.
    #[test]
    fn each_element_is_derived_from_its_own_nodes_alone() {
        let blocks = mda_blocks();
        let config = ConfigGeometry::from_element_points(&blocks).unwrap();

        assert_eq!(config.n_elements(), 3);
        let expected_total: usize = blocks.iter().map(|b| b.len()).sum();
        assert_eq!(config.total_nodes(), expected_total);
        assert_eq!(config.x().len(), expected_total);

        let mut start = 0usize;
        for (element, block) in blocks.iter().enumerate() {
            let alone = AirfoilGeometry::from_points(block).unwrap();
            let what = format!("element {element}");

            assert_eq!(config.element(element).start, start, "{what}: start");
            assert_eq!(config.element(element).n, block.len(), "{what}: n");

            assert_bits_eq(config.element_x(element), &alone.x, &format!("{what}.x"));
            assert_bits_eq(config.element_y(element), &alone.y, &format!("{what}.y"));
            assert_bits_eq(config.element_s(element), &alone.s, &format!("{what}.s"));
            assert_bits_eq(config.element_xp(element), &alone.xp, &format!("{what}.xp"));
            assert_bits_eq(config.element_yp(element), &alone.yp, &format!("{what}.yp"));
            assert_bits_eq(config.element_nx(element), &alone.nx, &format!("{what}.nx"));
            assert_bits_eq(config.element_ny(element), &alone.ny, &format!("{what}.ny"));
            assert_bits_eq(
                config.element_apanel(element),
                &alone.apanel,
                &format!("{what}.apanel"),
            );

            let geom = config.element(element);
            assert_eq!(geom.xte.to_bits(), alone.xte.to_bits(), "{what}: xte");
            assert_eq!(geom.yte.to_bits(), alone.yte.to_bits(), "{what}: yte");
            assert_eq!(geom.dste.to_bits(), alone.dste.to_bits(), "{what}: dste");
            assert_eq!(geom.ante.to_bits(), alone.ante.to_bits(), "{what}: ante");
            assert_eq!(geom.aste.to_bits(), alone.aste.to_bits(), "{what}: aste");
            assert_eq!(geom.xle.to_bits(), alone.xle.to_bits(), "{what}: xle");
            assert_eq!(geom.yle.to_bits(), alone.yle.to_bits(), "{what}: yle");
            assert_eq!(geom.sle.to_bits(), alone.sle.to_bits(), "{what}: sle");
            assert_eq!(geom.chord.to_bits(), alone.chord.to_bits(), "{what}: chord");
            assert_eq!(geom.sharp, alone.sharp, "{what}: sharp");

            // Arc length restarts at zero, so it is this element's own.
            assert_eq!(config.s()[start], 0.0, "{what}: s restarts");
            start += block.len();
        }
        assert_eq!(start, expected_total);

        // The three elements have genuinely different chords, so the per-element
        // normalisation is not a single shared number.
        let chords: Vec<f64> = config.elements().iter().map(|e| e.chord).collect();
        assert!(chords[0] < chords[1], "slat chord below main chord");
        assert!(chords[2] < chords[1], "flap chord below main chord");
    }

    /// Panel closure stays inside each element, and the last slot of an
    /// element's `apanel` is that element's own TE panel.
    #[test]
    fn panel_closure_stays_within_each_element() {
        let config = ConfigGeometry::from_element_points(&mda_blocks()).unwrap();

        for global in 0..config.total_nodes() {
            assert_eq!(
                config.element_of(config.next_node(global)),
                config.element_of(global),
                "next_node({global}) left its element"
            );
            assert_eq!(config.prev_node(config.next_node(global)), global);
        }

        for element in 0..config.n_elements() {
            let geom = config.element(element);
            assert_eq!(config.next_node(geom.end() - 1), geom.start);

            // The last panel angle of the element is XFOIL's TE-panel formula
            // applied to that element's own first and last nodes.
            let first = geom.start;
            let last = geom.end() - 1;
            let sx = config.x()[first] - config.x()[last];
            let sy = config.y()[first] - config.y()[last];
            let expected = (-sx).atan2(sy) + PI;
            assert_eq!(
                config.apanel()[last].to_bits(),
                expected.to_bits(),
                "element {element}: TE panel angle"
            );
        }
    }

    /// The `Configuration` route: paneling a real three-element configuration
    /// and building its geometry gives one element per configuration element,
    /// with the node counts the paneling produced plus the closing node a sharp
    /// trailing edge needs.
    #[test]
    fn config_geometry_from_a_paneled_configuration() {
        use rustfoil_core::paneling::PanelCounts;
        use rustfoil_core::{point, Body, Configuration, Element};

        let elements: Vec<Element> = mda_blocks()
            .iter()
            .enumerate()
            .map(|(i, block)| {
                let pts: Vec<_> = block.iter().map(|&(x, y)| point(x, y)).collect();
                Element::from_body(Body::from_points(&format!("element{i}"), &pts).unwrap())
            })
            .collect();
        let config = Configuration::new(elements);

        let counts = PanelCounts::Each(80);
        let paneled = config.panel_all(&counts).unwrap();
        let geometry = ConfigGeometry::from_paneled(&paneled).unwrap();

        assert_eq!(geometry.n_elements(), 3);
        // A closed element gains its closing node back; an open one does not.
        let restored = (0..3)
            .filter(|&element| paneled.element_is_closed(element))
            .count();
        assert!(restored > 0, "the fixture should have a sharp element");
        assert_eq!(geometry.total_nodes(), paneled.total_nodes() + restored);
        for element in 0..3 {
            let paneled_nodes = paneled.element_nodes(element);
            let closing = usize::from(paneled.element_is_closed(element));
            assert_eq!(
                geometry.element(element).n,
                paneled_nodes.len() + closing,
                "element {element}: node count"
            );
            // Nodes transfer unchanged from the paneling.
            for (node, paneled_node) in geometry.element_x(element).iter().zip(paneled_nodes) {
                assert_eq!(node.to_bits(), paneled_node.x.to_bits());
            }
            // And the restored node is the element's own first node.
            if closing == 1 {
                let geom = geometry.element(element);
                assert_eq!(
                    geometry.x()[geom.end() - 1].to_bits(),
                    geometry.x()[geom.start].to_bits()
                );
                assert_eq!(
                    geometry.y()[geom.end() - 1].to_bits(),
                    geometry.y()[geom.start].to_bits()
                );
            }
        }

        // The one-step constructor is the same thing.
        let direct = ConfigGeometry::from_configuration(&config, &counts).unwrap();
        assert_eq!(direct.total_nodes(), geometry.total_nodes());
        assert_bits_eq(direct.x(), geometry.x(), "from_configuration.x");
        assert_bits_eq(direct.y(), geometry.y(), "from_configuration.y");
        assert_bits_eq(
            direct.apanel(),
            geometry.apanel(),
            "from_configuration.apanel",
        );
    }

    #[test]
    fn an_empty_config_geometry_has_no_nodes() {
        let empty: [&[(f64, f64)]; 0] = [];
        let config = ConfigGeometry::from_element_points(&empty).unwrap();
        assert!(config.is_empty());
        assert_eq!(config.n_elements(), 0);
        assert_eq!(config.total_nodes(), 0);
        assert!(config.x().is_empty());
        assert!(config.try_element(0).is_none());
        assert!(config.sharp_te_bisector_control(0).is_none());
        assert!(config.te_coefficients_all().is_empty());
    }

    #[test]
    fn an_unusable_element_is_reported_rather_than_skipped() {
        let good = make_naca0012(80);
        let too_few: Vec<(f64, f64)> = good[..5].to_vec();
        let result = ConfigGeometry::from_element_points(&[good, too_few]);
        assert!(matches!(result, Err(InviscidError::InsufficientPoints(5))));
    }

    // =====================================================================
    // Leading and trailing edges of a real slat/main/flap section
    // =====================================================================

    /// The full 30P-30N section: 201, 221 and 242 coordinates for the slat, main
    /// element and flap.
    ///
    /// The trimmed copy in `testdata/` is the same three elements decimated to
    /// about a twentieth of the points, which loses the cove faces and the
    /// leading-edge resolution the landmark tests below turn on, so those use
    /// this one.
    fn full_mda_blocks() -> Vec<Vec<(f64, f64)>> {
        let path = repo_file("flexfoil-ui/public/airfoils/30p-30n.dat");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let mut blocks: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut current: Vec<(f64, f64)> = Vec::new();
        for line in text.lines() {
            let mut parts = line.trim().split_whitespace();
            let parsed = match (
                parts.next().and_then(|s| s.parse::<f64>().ok()),
                parts.next().and_then(|s| s.parse::<f64>().ok()),
            ) {
                (Some(x), Some(y)) => Some((x, y)),
                _ => None,
            };
            match parsed {
                Some(pair) => current.push(pair),
                None if !current.is_empty() => blocks.push(std::mem::take(&mut current)),
                None => {}
            }
        }
        if !current.is_empty() {
            blocks.push(current);
        }
        assert_eq!(
            blocks.iter().map(|b| b.len()).collect::<Vec<_>>(),
            vec![201, 221, 242],
            "expected the slat, main and flap blocks of the full 30P-30N"
        );
        blocks
    }

    /// The full 30P-30N as a `Configuration`, each element in the coordinates the
    /// file gives it.
    fn full_mda_configuration() -> rustfoil_core::Configuration {
        use rustfoil_core::{point, Body, Configuration, Element};
        let roles = ["slat", "main", "flap"];
        let elements: Vec<Element> = full_mda_blocks()
            .iter()
            .enumerate()
            .map(|(i, block)| {
                let pts: Vec<_> = block.iter().map(|&(x, y)| point(x, y)).collect();
                Element::from_body(Body::from_points(roles[i], &pts).unwrap())
            })
            .collect();
        Configuration::new(elements)
    }

    /// The node furthest from an element's trailing-edge midpoint, which brackets
    /// its leading edge to within one panel.
    fn furthest_node(x: &[f64], y: &[f64]) -> usize {
        let n = x.len();
        let xte = 0.5 * (x[0] + x[n - 1]);
        let yte = 0.5 * (y[0] + y[n - 1]);
        (0..n)
            .max_by(|&i, &j| {
                let di = (x[i] - xte).powi(2) + (y[i] - yte).powi(2);
                let dj = (x[j] - xte).powi(2) + (y[j] - yte).powi(2);
                di.partial_cmp(&dj).unwrap()
            })
            .unwrap()
    }

    /// A leading edge has to be at the front of the element, not on a cove face
    /// most of a chord downstream of it.
    ///
    /// The main element of this section has a straight vertical stretch where the
    /// flap tucks in — eight nodes sharing `x = 0.69993` — and the slat is
    /// deflected 30°, so on both of them the first place the contour turns back
    /// towards the trailing edge is not the leading edge. Checked on the raw
    /// coordinates and on the paneled contour at a spread of panel counts,
    /// because which element the local scan misses depends on the panel
    /// distribution.
    #[test]
    fn the_leading_edge_is_at_the_front_of_a_coved_element() {
        use rustfoil_core::paneling::PanelCounts;

        let blocks = full_mda_blocks();
        let raw = ConfigGeometry::from_element_points(&blocks).unwrap();
        let config = full_mda_configuration();

        let mut cases: Vec<(String, ConfigGeometry)> = vec![("raw coordinates".to_string(), raw)];
        for n in [60usize, 80, 100, 120, 140, 160, 200] {
            cases.push((
                format!("paneled Each({n})"),
                ConfigGeometry::from_configuration(&config, &PanelCounts::Each(n)).unwrap(),
            ));
        }

        for (what, geometry) in &cases {
            assert_eq!(geometry.n_elements(), 3, "{what}");
            for (element, role) in ["slat", "main", "flap"].iter().enumerate() {
                let geom = geometry.element(element);
                let x = geometry.element_x(element);
                let y = geometry.element_y(element);

                // The bracketing node, and the panels either side of it.
                let far = furthest_node(x, y);
                let s = geometry.element_s(element);
                let span = (s[(far + 1).min(geom.n - 1)] - s[far.saturating_sub(1)]).abs();

                let distance = ((geom.xle - x[far]).powi(2) + (geom.yle - y[far]).powi(2)).sqrt();
                assert!(
                    distance <= span,
                    "{what} {role}: leading edge ({}, {}) is {distance} from the furthest node \
                     ({}, {}), more than the {span} of arc length either side of it",
                    geom.xle,
                    geom.yle,
                    x[far],
                    y[far]
                );

                // And it is at the front of the element, not on its cove.
                let x_min = x.iter().copied().fold(f64::INFINITY, f64::min);
                let x_max = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                assert!(
                    geom.xle < x_min + 0.05 * (x_max - x_min),
                    "{what} {role}: leading edge at x = {} is not in the forward 5% of \
                     [{x_min}, {x_max}]",
                    geom.xle
                );
            }
        }
    }

    /// A sharp trailing edge is sharp on both construction routes.
    ///
    /// The slat and the main element of this section are closed contours and the
    /// flap is open, so one fixture covers both trailing-edge branches. The two
    /// routes differ in whether the closing point arrives with the coordinates:
    /// it does for the raw blocks of the file, and `from_paneled` restores it for
    /// the paneled contours, so both reach the convention the sharpness test
    /// needs.
    #[test]
    fn a_sharp_element_is_sharp_on_both_construction_routes() {
        use rustfoil_core::paneling::PanelCounts;

        let blocks = full_mda_blocks();
        let config = full_mda_configuration();
        let counts = PanelCounts::Each(120);
        let paneled = config.panel_all(&counts).unwrap();

        // The paneling's own verdict, which is the connectivity question.
        let closed: Vec<bool> = (0..3).map(|k| paneled.element_is_closed(k)).collect();
        assert_eq!(closed, vec![true, true, false], "slat, main, flap closure");

        for (what, geometry) in [
            (
                "raw coordinates",
                ConfigGeometry::from_element_points(&blocks).unwrap(),
            ),
            (
                "paneled",
                ConfigGeometry::from_configuration(&config, &counts).unwrap(),
            ),
        ] {
            for (element, role) in ["slat", "main", "flap"].iter().enumerate() {
                let geom = geometry.element(element);
                assert_eq!(
                    geom.sharp,
                    closed[element],
                    "{what} {role}: sharp = {}, but the contour is {}",
                    geom.sharp,
                    if closed[element] { "closed" } else { "open" }
                );

                if closed[element] {
                    // The two trailing-edge nodes are the same point, so there is
                    // no gap and the trailing-edge panel is the degenerate one.
                    assert_eq!(geom.dste, 0.0, "{what} {role}: dste");
                    assert_eq!(
                        geometry.apanel()[geom.end() - 1].to_bits(),
                        PI.to_bits(),
                        "{what} {role}: trailing-edge panel angle"
                    );
                    assert_eq!(
                        geometry.te_coefficients(element),
                        (1.0, 0.0),
                        "{what} {role}: te coefficients"
                    );
                    assert!(
                        geometry.sharp_te_bisector_control(element).is_some(),
                        "{what} {role}: bisector control point"
                    );
                } else {
                    assert!(
                        geom.dste > 1e-4 * geom.chord,
                        "{what} {role}: dste {} against chord {}",
                        geom.dste,
                        geom.chord
                    );
                    assert!(geometry.sharp_te_bisector_control(element).is_none());
                }
            }
        }
    }

    /// Landmark and chord values for the whole section, both routes, in one
    /// place: a slat and a main element whose leading edges are found by the
    /// global bracket, and a flap whose blunt trailing edge is measured across
    /// its own two corners.
    #[test]
    fn the_full_section_reports_a_landmark_per_element() {
        use rustfoil_core::paneling::PanelCounts;

        let config = full_mda_configuration();
        let raw = ConfigGeometry::from_element_points(&full_mda_blocks()).unwrap();
        let paneled = ConfigGeometry::from_configuration(&config, &PanelCounts::Each(120)).unwrap();

        // Element, xle, chord, sharp, tolerance on the two lengths. The
        // tolerances are the spread between the raw and the paneled contour, not
        // a claim about either one's precision.
        let expected = [
            ("slat", -0.0807, 0.1527, true),
            ("main", 0.0438, 0.8316, true),
            ("flap", 0.8735, 0.3003, false),
        ];
        for (what, geometry) in [("raw", &raw), ("paneled", &paneled)] {
            for (element, (role, xle, chord, sharp)) in expected.iter().enumerate() {
                let geom = geometry.element(element);
                assert!(
                    (geom.xle - xle).abs() < 5e-3,
                    "{what} {role}: xle {} against {xle}",
                    geom.xle
                );
                assert!(
                    (geom.chord - chord).abs() < 5e-3,
                    "{what} {role}: chord {} against {chord}",
                    geom.chord
                );
                assert_eq!(geom.sharp, *sharp, "{what} {role}: sharp");
            }
            // The main element is the largest, and it sets the D4 default.
            assert_eq!(
                geometry.main_element_index(),
                Some(1),
                "{what}: main element"
            );
            assert_eq!(
                geometry.default_ref_chord().to_bits(),
                geometry.element(1).chord.to_bits(),
                "{what}: reference chord"
            );
        }
    }

    /// Bracketing the leading-edge search leaves every single-element landmark
    /// exactly where it was.
    ///
    /// The values below were measured before the bracket was introduced, on the
    /// coordinate corpus in `testdata/` and on the analytic section. A
    /// conventional airfoil's first turn back towards the trailing edge *is* its
    /// leading edge, so the bracket agrees with the local scan, the iteration
    /// starts in the same place and lands in the same place, and these are raw
    /// bit patterns rather than tolerances because nothing about the arithmetic
    /// changed.
    #[test]
    fn improved_le_seeding_leaves_single_element_landmarks_bit_identical() {
        // (file, xle, yle, sle, chord)
        let pinned: [(&str, u64, u64, u64, u64); 5] = [
            (
                "naca0012.dat",
                0x0000000000000000,
                0x0000000000000000,
                0x3ff050471ff1e073,
                0x3ff0000000000000,
            ),
            (
                "naca2412.dat",
                0xbf1444200ec3adf8,
                0x3f59eb462913dccc,
                0x3ff06c2ef11e1d09,
                0x3ff00052605f7e89,
            ),
            (
                "naca0012_xfoil_paneled.dat",
                0xbe6bad1432c56e00,
                0x0000000000000000,
                0x3ff0505e655529ac,
                0x3ff000000dd68a19,
            ),
            (
                "naca0012_repaneled.dat",
                0x0000000000000000,
                0x0000000000000000,
                0x3ff01c98f7f3191b,
                0x3ff0000000000000,
            ),
            (
                "naca0012_buffer_real.dat",
                0x0000000000000000,
                0x0000000000000000,
                0x3ff05061b1e2a4fd,
                0x3ff0000000000000,
            ),
        ];

        for (file, xle, yle, sle, chord) in pinned {
            let blocks = read_dat_blocks(file);
            assert_eq!(blocks.len(), 1, "{file}: expected a single element");
            let geom = AirfoilGeometry::from_points(&blocks[0]).unwrap();
            for (name, actual, expected) in [
                ("xle", geom.xle.to_bits(), xle),
                ("yle", geom.yle.to_bits(), yle),
                ("sle", geom.sle.to_bits(), sle),
                ("chord", geom.chord.to_bits(), chord),
            ] {
                assert_eq!(
                    actual,
                    expected,
                    "{file}.{name}: {:#018x} vs {expected:#018x} ({})",
                    actual,
                    f64::from_bits(actual)
                );
            }
        }

        // And the analytic section, whose trailing edge is sharp.
        let geom = AirfoilGeometry::from_points(&make_naca0012(80)).unwrap();
        assert_eq!(geom.xle.to_bits(), 0x0000000000000000);
        assert_eq!(geom.yle.to_bits(), 0x0000000000000000);
        assert_eq!(geom.sle.to_bits(), 0x3ff0506853e7b5e7);
        assert_eq!(geom.chord.to_bits(), 0x3ff0000000000000);
    }

    /// The two chord definitions in the codebase, and how far apart they are.
    ///
    /// [`ElementGeometry::chord`] measures from the LEFIND leading edge and is
    /// the authority for anything aerodynamic.
    /// [`rustfoil_core::Element::chord`] measures from the contour's minimum-x
    /// node, which is what the geometric bookkeeping in rustfoil-core has
    /// available. This pins the gap so it cannot widen unnoticed.
    #[test]
    fn the_two_chord_definitions_disagree_by_a_pinned_amount() {
        use rustfoil_core::{point, Body, Element};

        let element_of = |points: &[(f64, f64)]| -> Element {
            let pts: Vec<_> = points.iter().map(|&(x, y)| point(x, y)).collect();
            Element::from_body(Body::from_points("element", &pts).unwrap())
        };

        // A conventional airfoil in its own coordinates: the minimum-x node is
        // within a rounded leading edge of the spline leading edge, so the two
        // agree to a few parts in 1e5.
        let single: [(&str, f64); 5] = [
            ("naca0012.dat", 0.0),
            ("naca2412.dat", 7.9e-5),
            ("naca0012_xfoil_paneled.dat", 2.6e-5),
            ("naca0012_repaneled.dat", 0.0),
            ("naca0012_buffer_real.dat", 0.0),
        ];
        for (file, bound) in single {
            let points = read_dat_blocks(file).remove(0);
            let inviscid = AirfoilGeometry::from_points(&points).unwrap().chord;
            let core = element_of(&points).chord();
            let relative = (core - inviscid).abs() / inviscid;
            assert!(
                relative <= bound.max(1e-15),
                "{file}: chords differ by {relative:e}, above the pinned {bound:e} \
                 (core {core}, inviscid {inviscid})"
            );
        }

        // An element whose deflection is baked into its coordinates is a
        // different matter: its minimum-x node is not near its leading edge, and
        // the disagreement is two orders larger.
        let deflected: [(&str, f64, f64); 3] = [
            ("slat", 0.02, 0.05),
            ("main", 0.0, 1e-4),
            ("flap", 5e-3, 0.02),
        ];
        for (block, (role, low, high)) in full_mda_blocks().iter().zip(deflected) {
            let inviscid = AirfoilGeometry::from_points(block).unwrap().chord;
            let core = element_of(block).chord();
            let relative = (core - inviscid).abs() / inviscid;
            assert!(
                relative >= low && relative <= high,
                "{role}: chords differ by {relative:e}, outside the pinned \
                 [{low:e}, {high:e}] (core {core}, inviscid {inviscid})"
            );
        }
    }
}
