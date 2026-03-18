//! XFOIL-faithful flap deflection (port of xgdes.f SUBROUTINE FLAP).
//!
//! Operates on a dense buffer airfoil, builds splines, uses SSS to find
//! exact break arc-lengths, rotates flap points, deletes folds, adds
//! arc-fill or corner points, runs SCHECK, and re-splines — matching
//! XFOIL 6.99 procedurally.

use crate::point::{point, Point};
use crate::xfoil_spline::Spline1D;

/// Deflect a flap on a buffer airfoil, matching XFOIL's GDES FLAP procedure.
///
/// * `pts` — buffer airfoil points (TE-upper → LE → TE-lower, Selig order)
/// * `hinge_x_frac` — hinge x as fraction of chord
/// * `hinge_y_frac` — hinge y as fraction of local thickness (0=lower, 1=upper)
/// * `deflection_deg` — flap deflection in degrees (positive = TE down)
///
/// Returns the modified buffer airfoil points ready for PANGEN repaneling.
pub fn xfoil_flap(
    pts: &[Point],
    hinge_x_frac: f64,
    hinge_y_frac: f64,
    deflection_deg: f64,
) -> Vec<Point> {
    let n = pts.len();
    if n < 6 || deflection_deg.abs() < 1e-10 {
        return pts.to_vec();
    }

    let le_idx = pts
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let x_min = pts[le_idx].x;
    let x_max = pts.iter().fold(f64::NEG_INFINITY, |a, p| a.max(p.x));
    let chord = (x_max - x_min).max(1e-10);
    let hinge_x = x_min + hinge_x_frac * chord;

    // ── Step 1: build spline of buffer (SCALC + SEGSPL) ──────────────
    let mut xb: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let mut yb: Vec<f64> = pts.iter().map(|p| p.y).collect();
    let mut sb = scalc(&xb, &yb);
    let nb = xb.len();

    let x_spl = Spline1D::new(&sb, &xb).unwrap();
    let y_spl = Spline1D::new(&sb, &yb).unwrap();

    // ── Step 2: GETXYF — find arc-length of hinge on upper/lower ─────
    let mut tops = sb[0] + (xb[0] - hinge_x);
    tops = sinvrt(tops, hinge_x, &x_spl);
    let mut bots = sb[nb - 1] - (xb[nb - 1] - hinge_x);
    bots = sinvrt(bots, hinge_x, &x_spl);

    let topy = y_spl.eval(tops);
    let boty = y_spl.eval(bots);
    let hinge_y = boty + hinge_y_frac.clamp(0.0, 1.0) * (topy - boty);

    // ── Step 3: determine which side folds / opens ───────────────────
    let rdef = deflection_deg.to_radians();
    let inside = is_inside(&xb, &yb, hinge_x, hinge_y);

    let (atop, abot) = if inside {
        (0.0_f64.max(-rdef), 0.0_f64.max(rdef))
    } else {
        let chx = x_spl.deriv(bots) - x_spl.deriv(tops);
        let chy = y_spl.deriv(bots) - y_spl.deriv(tops);
        let fvx = x_spl.eval(bots) + x_spl.eval(tops);
        let fvy = y_spl.eval(bots) + y_spl.eval(tops);
        let crsp = chx * (hinge_y - 0.5 * fvy) - chy * (hinge_x - 0.5 * fvx);
        if crsp > 0.0 {
            (0.0_f64.max(rdef), 0.0_f64.max(rdef))
        } else {
            (0.0_f64.max(-rdef), 0.0_f64.max(-rdef))
        }
    };

    // ── Step 4: SSS — find break arc-lengths ─────────────────────────
    let (st1, st2) = sss(tops, atop, hinge_x, hinge_y, &x_spl, &y_spl, 1);
    let (sb1, sb2) = sss(bots, abot, hinge_x, hinge_y, &x_spl, &y_spl, 2);

    // break-point coordinates
    let xt1 = x_spl.eval(st1);
    let yt1 = y_spl.eval(st1);
    let xt2 = x_spl.eval(st2);
    let yt2 = y_spl.eval(st2);
    let xb1 = x_spl.eval(sb1);
    let yb1 = y_spl.eval(sb1);
    let xb2 = x_spl.eval(sb2);
    let yb2 = y_spl.eval(sb2);

    // ── Step 5: find adjacent buffer-point indices ───────────────────
    let mut it1 = 1usize;
    let mut it2 = 1usize;
    let mut ib1 = 1usize;
    let mut ib2 = 1usize;
    for i in 0..nb - 1 {
        if sb[i] <= st1 && sb[i + 1] > st1 { it1 = i + 1; }
        if sb[i] < st2 && sb[i + 1] >= st2 { it2 = i; }
        if sb[i] <= sb1 && sb[i + 1] > sb1 { ib1 = i; }
        if sb[i] < sb2 && sb[i + 1] >= sb2 { ib2 = i + 1; }
    }

    let dsavg = (sb[nb - 1] - sb[0]) / (nb - 1) as f64;
    let sfrac = 1.0 / 3.0;

    // ── helper points for top break ──────────────────────────────────
    let mut lt1new = false;
    let mut lt2new = false;
    let mut xt1new = xt1;
    let mut yt1new = yt1;
    let mut xt2new = xt2;
    let mut yt2new = yt2;

    if atop != 0.0 {
        let st1q = st1 + sfrac * (sb[it1.min(nb - 1).saturating_add(1).min(nb - 1)] - st1);
        if sb[it1.min(nb - 1)] < st1q {
            xt1new = x_spl.eval(st1q);
            yt1new = y_spl.eval(st1q);
        } else {
            let st1p = st1 + sfrac * (sb[it1.min(nb - 1)] - st1);
            xt1new = x_spl.eval(st1p);
            yt1new = y_spl.eval(st1p);
            lt1new = true;
        }
        let it2q = it2.saturating_sub(1).max(0);
        let st2q = st2 + sfrac * (sb[it2q] - st2);
        if sb[it2.min(nb - 1)] > st2q {
            xt2new = x_spl.eval(st2q);
            yt2new = y_spl.eval(st2q);
        } else {
            let st2p = st2 + sfrac * (sb[it2.min(nb - 1)] - st2);
            xt2new = x_spl.eval(st2p);
            yt2new = y_spl.eval(st2p);
            lt2new = true;
        }
    }

    // ── helper points for bottom break ───────────────────────────────
    let mut lb1new = false;
    let mut lb2new = false;
    let mut xb1new = xb1;
    let mut yb1new = yb1;
    let mut xb2new = xb2;
    let mut yb2new = yb2;

    if abot != 0.0 {
        let sb1q = sb1 + sfrac * (sb[ib1.saturating_sub(1)] - sb1);
        if sb[ib1.min(nb - 1)] > sb1q {
            xb1new = x_spl.eval(sb1q);
            yb1new = y_spl.eval(sb1q);
        } else {
            let sb1p = sb1 + sfrac * (sb[ib1.min(nb - 1)] - sb1);
            xb1new = x_spl.eval(sb1p);
            yb1new = y_spl.eval(sb1p);
            lb1new = true;
        }
        let ib2q = (ib2 + 1).min(nb - 1);
        let sb2q = sb2 + sfrac * (sb[ib2q] - sb2);
        if sb[ib2.min(nb - 1)] < sb2q {
            xb2new = x_spl.eval(sb2q);
            yb2new = y_spl.eval(sb2q);
        } else {
            let sb2p = sb2 + sfrac * (sb[ib2.min(nb - 1)] - sb2);
            xb2new = x_spl.eval(sb2p);
            yb2new = y_spl.eval(sb2p);
            lb2new = true;
        }
    }

    let sind = rdef.sin();
    let cosd = rdef.cos();

    // ── Step 6: rotate flap points ───────────────────────────────────
    for i in 0..nb {
        if i >= it1 && i <= ib1 { continue; }
        let dx = xb[i] - hinge_x;
        let dy = yb[i] - hinge_y;
        xb[i] = hinge_x + dx * cosd + dy * sind;
        yb[i] = hinge_y - dx * sind + dy * cosd;
    }

    // ── Step 7: delete disappeared points on upper surface ───────────
    let idif_top = if it1 > it2 + 1 { it1 - it2 - 1 } else { 0 };
    if idif_top > 0 {
        xb.drain(it2 + 1..it2 + 1 + idif_top);
        yb.drain(it2 + 1..it2 + 1 + idif_top);
        it1 -= idif_top;
        ib1 -= idif_top;
        ib2 -= idif_top;
    }

    let idif_bot = if ib2 > ib1 + 1 { ib2 - ib1 - 1 } else { 0 };
    if idif_bot > 0 {
        xb.drain(ib1 + 1..ib1 + 1 + idif_bot);
        yb.drain(ib1 + 1..ib1 + 1 + idif_bot);
        ib2 -= idif_bot;
    }

    // ── Step 8: add new points at top break ──────────────────────────
    if atop == 0.0 {
        // Gap side: skip arc-fill (our corner-aware CubicSpline handles the
        // break; XFOIL's arc-fill creates x-reversals that CubicSpline can't
        // repanel without oscillation — XFOIL's own SEGSPL handles them fine)
    } else {
        // Fold side: place corner + helper points
        let npadd = 1 + lt1new as usize + lt2new as usize;
        for _ in 0..npadd {
            xb.insert(it1, 0.0);
            yb.insert(it1, 0.0);
        }
        it1 += npadd;
        ib1 += npadd;
        ib2 += npadd;

        if lt1new {
            xb[it1 - 1] = xt1new;
            yb[it1 - 1] = yt1new;
            xb[it1 - 2] = xt1;
            yb[it1 - 2] = yt1;
        } else {
            xb[it1] = xt1new;
            yb[it1] = yt1new;
            xb[it1 - 1] = xt1;
            yb[it1 - 1] = yt1;
        }

        let xbar = xt2new - hinge_x;
        let ybar = yt2new - hinge_y;
        let idx = if lt2new { it2 + 1 } else { it2 };
        if idx < xb.len() {
            xb[idx] = hinge_x + xbar * cosd + ybar * sind;
            yb[idx] = hinge_y - xbar * sind + ybar * cosd;
        }
    }

    // ── Step 9: add new points at bottom break ───────────────────────
    if abot == 0.0 {
        // Gap side: skip arc-fill (same reasoning as top)
    } else {
        let npadd = 1 + lb1new as usize + lb2new as usize;
        for _ in 0..npadd {
            xb.insert(ib2, 0.0);
            yb.insert(ib2, 0.0);
        }
        ib2 += npadd;

        if lb1new {
            xb[ib1 + 1] = xb1new;
            yb[ib1 + 1] = yb1new;
            xb[ib1 + 2] = xb1;
            yb[ib1 + 2] = yb1;
        } else {
            xb[ib1] = xb1new;
            yb[ib1] = yb1new;
            xb[ib1 + 1] = xb1;
            yb[ib1 + 1] = yb1;
        }

        let xbar = xb2new - hinge_x;
        let ybar = yb2new - hinge_y;
        let idx = if lb2new { ib2 - 1 } else { ib2 };
        if idx < xb.len() {
            xb[idx] = hinge_x + xbar * cosd + ybar * sind;
            yb[idx] = hinge_y - xbar * sind + ybar * cosd;
        }
    }

    // ── Step 10: SCHECK — remove splinter segments ───────────────────
    scheck(&mut xb, &mut yb, 0.2);

    // Convert back to points
    xb.iter().zip(yb.iter()).map(|(&x, &y)| point(x, y)).collect()
}

// ═══════════════════════════════════════════════════════════════════════
// Helper functions matching XFOIL's Fortran subroutines
// ═══════════════════════════════════════════════════════════════════════

fn scalc(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut s = vec![0.0; n];
    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
    s
}

/// XFOIL SINVRT — Newton inversion of spline to find s where x(s) = xi.
fn sinvrt(mut si: f64, xi: f64, spl: &Spline1D) -> f64 {
    let si_sav = si;
    for _ in 0..10 {
        let res = spl.eval(si) - xi;
        let resp = spl.deriv(si);
        if resp.abs() < 1e-30 { break; }
        let ds = -res / resp;
        si += ds;
        let stot = spl.s_range();
        if stot > 0.0 && (ds / stot).abs() < 1e-5 {
            return si;
        }
    }
    si_sav
}

/// XFOIL SSS — find break arc-lengths S1, S2 at flap surface break.
fn sss(
    ss: f64,
    del: f64,
    xbf: f64,
    ybf: f64,
    x_spl: &Spline1D,
    y_spl: &Spline1D,
    iside: usize,
) -> (f64, f64) {
    let eps = 1.0e-5;
    let stot = x_spl.s_range();
    let sind_half = (0.5 * del.abs()).sin();
    let ssgn = if iside == 1 { -1.0 } else { 1.0 };

    let rsq = (x_spl.eval(ss) - xbf).powi(2) + (y_spl.eval(ss) - ybf).powi(2);
    let mut s1 = ss - (sind_half * rsq.sqrt() + eps * stot) * ssgn;
    let mut s2 = ss + (sind_half * rsq.sqrt() + eps * stot) * ssgn;

    for _ in 0..10 {
        let x1 = x_spl.eval(s1);
        let x1p = x_spl.deriv(s1);
        let y1 = y_spl.eval(s1);
        let y1p = y_spl.deriv(s1);

        let x2 = x_spl.eval(s2);
        let x2p = x_spl.deriv(s2);
        let y2 = y_spl.eval(s2);
        let y2p = y_spl.deriv(s2);

        let r1 = ((x1 - xbf).powi(2) + (y1 - ybf).powi(2)).sqrt();
        let r2 = ((x2 - xbf).powi(2) + (y2 - ybf).powi(2)).sqrt();
        let rr = ((x1 - x2).powi(2) + (y1 - y2).powi(2)).sqrt();

        if r1 <= eps * stot || r2 <= eps * stot {
            return (ss, ss);
        }

        let r1_s1 = (x1p * (x1 - xbf) + y1p * (y1 - ybf)) / r1;
        let r2_s2 = (x2p * (x2 - xbf) + y2p * (y2 - ybf)) / r2;

        let (rs1, a11, a12, rs2, a21, a22);

        if sind_half > 0.01 {
            if rr == 0.0 { return (s1, s2); }

            let rr_s1 = (x1p * (x1 - x2) + y1p * (y1 - y2)) / rr;
            let rr_s2 = -(x2p * (x1 - x2) + y2p * (y1 - y2)) / rr;

            rs1 = ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) / rr - sind_half * r1;
            a11 = ((xbf - x1) * (-x1p) + (ybf - y1) * (-y1p)) / rr
                + ((-x1p) * (x2 - x1) + (-y1p) * (y2 - y1)) / rr
                - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rr_s1 / (rr * rr)
                - sind_half * r1_s1;
            a12 = ((xbf - x1) * x2p + (ybf - y1) * y2p) / rr
                - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rr_s2 / (rr * rr);
            rs2 = r1 - r2;
            a21 = r1_s1;
            a22 = -r2_s2;
        } else {
            rs1 = (r1 + r2) * sind_half + (s1 - s2) * ssgn;
            a11 = r1_s1 * sind_half + ssgn;
            a12 = r2_s2 * sind_half - ssgn;

            let (x1pp, y1pp) = d2val_1d_pair(s1, x_spl, y_spl);
            let (x2pp, y2pp) = d2val_1d_pair(s2, x_spl, y_spl);

            let xtot = x1 + x2 - 2.0 * xbf;
            let ytot = y1 + y2 - 2.0 * ybf;

            rs2 = xtot * (x1p + x2p) + ytot * (y1p + y2p);
            a21 = x1p * (x1p + x2p) + y1p * (y1p + y2p) + xtot * x1pp + ytot * y1pp;
            a22 = x2p * (x1p + x2p) + y2p * (y1p + y2p) + xtot * x2pp + ytot * y2pp;
        }

        let det = a11 * a22 - a12 * a21;
        if det.abs() < 1e-30 { break; }
        let mut ds1 = -(rs1 * a22 - a12 * rs2) / det;
        let mut ds2 = -(a11 * rs2 - rs1 * a21) / det;

        let lim = 0.01 * stot;
        ds1 = ds1.clamp(-lim, lim);
        ds2 = ds2.clamp(-lim, lim);

        s1 += ds1;
        s2 += ds2;
        if (ds1.abs() + ds2.abs()) < eps * stot { break; }
    }

    if del == 0.0 {
        let avg = 0.5 * (s1 + s2);
        (avg, avg)
    } else {
        (s1, s2)
    }
}

fn d2val_1d_pair(ss: f64, x_spl: &Spline1D, y_spl: &Spline1D) -> (f64, f64) {
    (x_spl.d2val(ss), y_spl.d2val(ss))
}

/// XFOIL INSIDE — test if (xf, yf) is inside the airfoil contour.
fn is_inside(x: &[f64], y: &[f64], xf: f64, yf: f64) -> bool {
    let n = x.len();
    let mut cross = 0i32;
    for i in 0..n {
        let j = (i + 1) % n;
        let (yi, yj) = (y[i] - yf, y[j] - yf);
        if (yi > 0.0) == (yj > 0.0) { continue; }
        let xi = x[i] + yi * (x[j] - x[i]) / (yi - yj) - xf;
        if xi > 0.0 { cross += 1; }
    }
    cross % 2 != 0
}

/// XFOIL SCHECK — remove splinter segments.
fn scheck(x: &mut Vec<f64>, y: &mut Vec<f64>, stol: f64) {
    loop {
        let mut changed = false;
        let mut i = 1;
        while i + 2 < x.len() {
            let dsm1 = ((x[i] - x[i - 1]).powi(2) + (y[i] - y[i - 1]).powi(2)).sqrt();
            let dsp1 = ((x[i + 1] - x[i]).powi(2) + (y[i + 1] - y[i]).powi(2)).sqrt();
            let dsp2 = ((x[i + 2] - x[i + 1]).powi(2) + (y[i + 2] - y[i + 1]).powi(2)).sqrt();

            if dsp1 == 0.0 { i += 1; continue; }

            if dsp1 < stol * dsm1 || dsp1 < stol * dsp2 {
                x[i] = 0.5 * (x[i] + x[i + 1]);
                y[i] = 0.5 * (y[i] + y[i + 1]);
                x.remove(i + 1);
                y.remove(i + 1);
                changed = true;
                continue; // re-check at same index
            }
            i += 1;
        }
        if !changed { break; }
    }
}
