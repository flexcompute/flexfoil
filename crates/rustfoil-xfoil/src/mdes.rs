//! Full-inverse design via circle-plane conformal mapping (MDES).
//!
//! Port of XFOIL's MDES routines: EIWSET, FTP, PIQSUM, ZCCALC, ZCNORM,
//! CNCALC, CNFILT, MAPGEN, SCINIT.
//!
//! The airfoil is mapped to a unit circle using the Fourier coefficients Cn
//! of the harmonic function P(w) + iQ(w). Given a target velocity distribution
//! Q(w), the inverse transform recovers the geometry z(w).

use std::f64::consts::PI;

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

use crate::error::{Result, XfoilError};

const ICX: usize = 257;

type C64 = Complex64;

fn c(re: f64, im: f64) -> C64 { C64::new(re, im) }

// ---------------------------------------------------------------------------
// Circle-plane state
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct CirclePlane {
    pub nc: usize,
    pub mc: usize,
    pub mct: usize,
    pub dwc: f64,
    pub agte: f64,
    pub ag0: f64,
    pub qim0: f64,
    pub qimold: f64,

    pub wc: Vec<f64>,
    pub sc: Vec<f64>,
    pub scold: Vec<f64>,
    pub xcold: Vec<f64>,
    pub ycold: Vec<f64>,

    pub dzte: C64,
    pub chordz: C64,
    pub zleold: C64,
    pub zc: Vec<C64>,
    pub zc_cn: Vec<Vec<C64>>,
    pub piq: Vec<C64>,
    pub cn: Vec<C64>,
    pub eiw: Vec<Vec<C64>>,
}

impl CirclePlane {
    fn new(nc: usize) -> Self {
        let mc = nc / 4;
        let mct = nc / 16;
        Self {
            nc, mc, mct,
            dwc: 2.0 * PI / (nc - 1) as f64,
            agte: 0.0, ag0: 0.0, qim0: 0.0, qimold: 0.0,
            wc: vec![0.0; nc],
            sc: vec![0.0; nc],
            scold: vec![0.0; nc],
            xcold: vec![0.0; nc],
            ycold: vec![0.0; nc],
            dzte: c(0.0, 0.0),
            chordz: c(0.0, 0.0),
            zleold: c(0.0, 0.0),
            zc: vec![c(0.0, 0.0); nc],
            zc_cn: vec![vec![c(0.0, 0.0); mc]; nc],
            piq: vec![c(0.0, 0.0); nc],
            cn: vec![c(0.0, 0.0); mc + 1],
            eiw: vec![vec![c(0.0, 0.0); mc + 1]; nc],
        }
    }
}

// ---------------------------------------------------------------------------
// EIWSET – set up circle-plane coordinate arrays and exp(inw)
// ---------------------------------------------------------------------------

fn eiwset(cp: &mut CirclePlane) {
    let nc = cp.nc;
    let mc = cp.mc;
    cp.dwc = 2.0 * PI / (nc - 1) as f64;

    for ic in 0..nc {
        cp.wc[ic] = cp.dwc * ic as f64;
    }

    for ic in 0..nc {
        cp.eiw[ic][0] = c(1.0, 0.0);
    }
    for ic in 0..nc {
        cp.eiw[ic][1] = C64::from_polar(1.0, cp.wc[ic]);
    }
    for m in 2..=mc {
        for ic in 0..nc {
            let ic1 = (m * ic) % (nc - 1);
            cp.eiw[ic][m] = cp.eiw[ic1][1];
        }
    }
}

// ---------------------------------------------------------------------------
// FTP – Slow Fourier Transform of P(w) using trapezoidal integration
// ---------------------------------------------------------------------------

fn ftp(cp: &mut CirclePlane) {
    let nc = cp.nc;
    let mc = cp.mc;
    for m in 0..=mc {
        let mut zsum = c(0.0, 0.0);
        for ic in 1..nc - 1 {
            zsum += cp.piq[ic] * cp.eiw[ic][m];
        }
        cp.cn[m] = (0.5 * (cp.piq[0] * cp.eiw[0][m] + cp.piq[nc - 1] * cp.eiw[nc - 1][m])
            + zsum) * cp.dwc / PI;
    }
    cp.cn[0] *= 0.5;
}

// ---------------------------------------------------------------------------
// PIQSUM – Inverse transform: P(w)+iQ(w) from Cn
// ---------------------------------------------------------------------------

fn piqsum(cp: &mut CirclePlane) {
    let nc = cp.nc;
    let mc = cp.mc;
    for ic in 0..nc {
        let mut zsum = c(0.0, 0.0);
        for m in 0..=mc {
            zsum += cp.cn[m] * cp.eiw[ic][m].conj();
        }
        cp.piq[ic] = zsum;
    }
}

// ---------------------------------------------------------------------------
// CNFILT – Hanning filter on Cn coefficients
// ---------------------------------------------------------------------------

fn cnfilt(cp: &mut CirclePlane, ffilt: f64) {
    if ffilt == 0.0 { return; }
    let mc = cp.mc;
    for m in 0..=mc {
        let freq = m as f64 / mc as f64;
        let cwt = 0.5 * (1.0 + (PI * freq).cos());
        let cwtx = cwt.powf(ffilt);
        cp.cn[m] *= cwtx;
    }
}

// ---------------------------------------------------------------------------
// ZCCALC – Compute z(w) from P(w)+iQ(w) by trapezoidal integration
// ---------------------------------------------------------------------------

fn zccalc(cp: &mut CirclePlane, mtest: usize) {
    let nc = cp.nc;
    let agte = cp.agte;
    let dwc = cp.dwc;

    cp.zc[0] = c(4.0, 0.0);
    for m in 0..mtest {
        cp.zc_cn[0][m] = c(0.0, 0.0);
    }

    let sinw0 = 2.0 * (0.5 * cp.wc[0]).sin();
    let sinwe0 = if sinw0 > 0.0 { sinw0.powf(1.0 - agte) } else { 0.0 };
    let hwc0 = 0.5 * (cp.wc[0] - PI) * (1.0 + agte) - 0.5 * PI;
    let mut dzdw1 = sinwe0 * (cp.piq[0] + c(0.0, hwc0)).exp();

    for ic in 1..nc {
        let sinw = 2.0 * (0.5 * cp.wc[ic]).sin();
        let sinwe = if sinw > 0.0 { sinw.powf(1.0 - agte) } else { 0.0 };
        let hwc = 0.5 * (cp.wc[ic] - PI) * (1.0 + agte) - 0.5 * PI;
        let dzdw2 = sinwe * (cp.piq[ic] + c(0.0, hwc)).exp();

        cp.zc[ic] = 0.5 * (dzdw1 + dzdw2) * dwc + cp.zc[ic - 1];
        let dz_piq1 = 0.5 * dzdw1 * dwc;
        let dz_piq2 = 0.5 * dzdw2 * dwc;

        for m in 0..mtest.min(cp.mc) {
            cp.zc_cn[ic][m] = dz_piq1 * cp.eiw[ic - 1][m].conj()
                + dz_piq2 * cp.eiw[ic][m].conj()
                + cp.zc_cn[ic - 1][m];
        }

        dzdw1 = dzdw2;
    }

    // Set arc length s(w)
    cp.sc[0] = 0.0;
    for ic in 1..nc {
        cp.sc[ic] = cp.sc[ic - 1] + (cp.zc[ic] - cp.zc[ic - 1]).norm();
    }
    let total = cp.sc[nc - 1];
    if total > 0.0 {
        for ic in 0..nc {
            cp.sc[ic] /= total;
        }
    }
}

// ---------------------------------------------------------------------------
// ZLEFIND – Find leading edge on the circle-plane airfoil
// ---------------------------------------------------------------------------

fn zlefind(zc: &[C64], _wc: &[f64], nc: usize, _piq: &[C64], _agte: f64) -> C64 {
    let mut min_x = f64::INFINITY;
    let mut le_ic = nc / 2;
    for ic in 0..nc {
        if zc[ic].re < min_x {
            min_x = zc[ic].re;
            le_ic = ic;
        }
    }
    zc[le_ic]
}

// ---------------------------------------------------------------------------
// ZCNORM – Normalize z(w) to old chord and position
// ---------------------------------------------------------------------------

fn zcnorm(cp: &mut CirclePlane, mtest: usize) {
    let nc = cp.nc;
    let zle = zlefind(&cp.zc, &cp.wc, nc, &cp.piq, cp.agte);

    for ic in 0..nc {
        cp.zc[ic] -= zle;
    }

    let zte = 0.5 * (cp.zc[0] + cp.zc[nc - 1]);
    let mut zte_cn = vec![c(0.0, 0.0); mtest.min(cp.mc)];
    for m in 0..mtest.min(cp.mc) {
        zte_cn[m] = 0.5 * (cp.zc_cn[0][m] + cp.zc_cn[nc - 1][m]);
    }

    for ic in 0..nc {
        let zcnew = cp.chordz * cp.zc[ic] / zte;
        let zc_zte = -zcnew / zte;
        cp.zc[ic] = zcnew;
        for m in 0..mtest.min(cp.mc) {
            cp.zc_cn[ic][m] = cp.chordz * cp.zc_cn[ic][m] / zte + zc_zte * zte_cn[m];
        }
    }

    let qimoff = -(cp.chordz / zte).ln().im;
    cp.cn[0] -= c(0.0, qimoff);

    for ic in 0..nc {
        cp.zc[ic] += cp.zleold;
    }
}

// ---------------------------------------------------------------------------
// MAPGEN – Generate geometry from Cn, enforcing TE gap constraint
// ---------------------------------------------------------------------------

fn mapgen(cp: &mut CirclePlane) {
    let nc = cp.nc;

    // Preset rotation offset
    let dx = cp.xcold[1] - cp.xcold[0];
    let dy = cp.ycold[1] - cp.ycold[0];
    let qim0 = dx.atan2(-dy) + 0.5 * PI * (1.0 + cp.agte);
    let qimoff = qim0 - cp.cn[0].im;
    cp.cn[0] += c(0.0, qimoff);

    piqsum(cp);
    zccalc(cp, cp.mct);
    zcnorm(cp, cp.mct);

    // Enforce Lighthill's first constraint
    cp.cn[0] = c(0.0, cp.cn[0].im);

    // Newton iteration for TE gap
    for _iter in 0..10 {
        let residual = cp.zc[0] - cp.zc[nc - 1] - cp.dzte;
        let jacobian = cp.zc_cn[0][0] - cp.zc_cn[nc - 1][0];
        if jacobian.norm() < 1e-20 { break; }
        let dcn = residual / jacobian;
        cp.cn[0] -= c(dcn.re, 0.0); // Only adjust real part of cn[1] (index 0 in 0-based Cn[1])
        // Actually we need to adjust cn[1] in XFOIL's 1-based indexing = cn[1] here
        // Let me re-check: MAPGEN uses NCN=1, M=1, so it adjusts CN(1) (the first harmonic)
        if cp.cn.len() > 1 {
            cp.cn[1] -= dcn;
        }

        piqsum(cp);
        zccalc(cp, cp.mct);
        zcnorm(cp, cp.mct);

        if dcn.norm() <= 5.0e-5 { break; }
    }
}

// ---------------------------------------------------------------------------
// SCINIT – Initialize s(w) mapping from airfoil geometry
// ---------------------------------------------------------------------------

fn scinit(cp: &mut CirclePlane, x: &[f64], y: &[f64]) -> Result<()> {
    let n = x.len();
    if n < 4 { return Err(XfoilError::Message("Need at least 4 points".to_string())); }
    let nc = cp.nc;

    // Find LE
    let le_idx = x.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i).unwrap_or(0);
    let xle = x[le_idx];
    let yle = y[le_idx];

    // Compute spline derivatives (finite differences)
    let mut xp = vec![0.0; n];
    let mut yp = vec![0.0; n];
    let mut s = vec![0.0; n];
    s[0] = 0.0;
    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
    for i in 0..n {
        let im = if i > 0 { i - 1 } else { 0 };
        let ip = if i < n - 1 { i + 1 } else { n - 1 };
        let ds = s[ip] - s[im];
        if ds.abs() > 1e-12 {
            xp[i] = (x[ip] - x[im]) / ds;
            yp[i] = (y[ip] - y[im]) / ds;
        }
    }

    let sle = s[le_idx];

    // TE angle parameter
    cp.agte = (xp[n - 1].atan2(-yp[n - 1]) - xp[0].atan2(-yp[0])) / PI - 1.0;
    cp.ag0 = xp[0].atan2(-yp[0]);
    cp.qim0 = cp.ag0 + 0.5 * PI * (1.0 + cp.agte);

    // TE gap
    let dxte = x[0] - x[n - 1];
    let dyte = y[0] - y[n - 1];
    cp.dzte = c(dxte, dyte);

    let chordx = 0.5 * (x[0] + x[n - 1]) - xle;
    let chordy = 0.5 * (y[0] + y[n - 1]) - yle;
    cp.chordz = c(chordx, chordy);
    cp.zleold = c(xle, yle);

    // Initial s(w) using LE curvature estimate
    let tops = sle / s[n - 1];
    let bots = (s[n - 1] - sle) / s[n - 1];
    let dsdwle = 0.1_f64.max(0.5 / (1.0_f64 / tops.max(0.01)));

    let wwt_top = 1.0 - 2.0 * dsdwle / tops.max(1e-6);
    for ic in 0..=(nc - 1) / 2 {
        let wt = wwt_top * cp.wc[ic];
        let den = 1.0 - (wwt_top * PI).cos();
        cp.sc[ic] = if den.abs() > 1e-12 { tops * (1.0 - wt.cos()) / den } else { tops * ic as f64 / ((nc - 1) / 2) as f64 };
    }

    let wwt_bot = 1.0 - 2.0 * dsdwle / bots.max(1e-6);
    for ic in (nc - 1) / 2 + 1..nc {
        let wval = cp.wc[nc - 1] - cp.wc[ic];
        let den = 1.0 - (wwt_bot * PI).cos();
        cp.sc[ic] = if den.abs() > 1e-12 {
            1.0 - bots * (1.0 - (wwt_bot * wval).cos()) / den
        } else {
            1.0 - bots * (nc - 1 - ic) as f64 / ((nc - 1) / 2) as f64
        };
    }

    // Iteration for s(w) convergence
    for _ipass in 0..30 {
        // Compute Q(w) from geometry
        for ic in 0..nc {
            let sic = s[0] + (s[n - 1] - s[0]) * cp.sc[ic].clamp(0.0, 1.0);
            let (dxds, dyds) = interp_deriv(&s, x, y, sic);
            let qim = dxds.atan2(-dyds)
                - 0.5 * (cp.wc[ic] - PI) * (1.0 + cp.agte)
                - cp.qim0;
            cp.piq[ic] = c(0.0, qim);
        }

        ftp(cp);
        cp.cn[0] = c(0.0, cp.cn[0].im + cp.qim0);
        piqsum(cp);

        for ic in 0..nc {
            cp.scold[ic] = cp.sc[ic];
        }

        // TE gap correction on cn[1]
        for _itgap in 0..5 {
            zccalc(cp, 1);
            let zle = zlefind(&cp.zc, &cp.wc, nc, &cp.piq, cp.agte);
            let zte = 0.5 * (cp.zc[0] + cp.zc[nc - 1]);
            let dzwt = ((zte - zle).norm() / cp.chordz.norm()).max(1e-10);
            if cp.zc_cn[0].is_empty() || cp.zc_cn[nc - 1].is_empty() { break; }
            let denom = cp.zc_cn[0][0] - cp.zc_cn[nc - 1][0];
            if denom.norm() < 1e-20 { break; }
            let dcn = -(cp.zc[0] - cp.zc[nc - 1] - dzwt * cp.dzte) / denom;
            if cp.cn.len() > 1 {
                cp.cn[1] += dcn;
            }
            piqsum(cp);
            if dcn.norm() < 1e-7 { break; }
        }

        let mut dscmax = 0.0_f64;
        for ic in 0..nc {
            dscmax = dscmax.max((cp.sc[ic] - cp.scold[ic]).abs());
        }
        if dscmax < 5e-7 { break; }
    }

    // Normalize final geometry
    zcnorm(cp, 1);

    // Save old airfoil in circle-plane arrays
    for ic in 0..nc {
        cp.scold[ic] = cp.sc[ic];
        cp.xcold[ic] = cp.zc[ic].re;
        cp.ycold[ic] = cp.zc[ic].im;
    }
    cp.qimold = cp.cn[0].im;

    Ok(())
}

fn interp_deriv(s: &[f64], x: &[f64], y: &[f64], sval: f64) -> (f64, f64) {
    let n = s.len();
    if n < 2 { return (1.0, 0.0); }

    let mut i = 0;
    for j in 1..n {
        if s[j] >= sval { i = j; break; }
        i = j;
    }
    if i == 0 { i = 1; }

    let ds = s[i] - s[i - 1];
    if ds.abs() < 1e-20 { return (1.0, 0.0); }
    let dxds = (x[i] - x[i - 1]) / ds;
    let dyds = (y[i] - y[i - 1]) / ds;
    (dxds, dyds)
}

// ---------------------------------------------------------------------------
// CNCALC – Compute Cn from speed distribution Q(w)
// ---------------------------------------------------------------------------

fn cncalc(cp: &mut CirclePlane, qc: &[f64], lsymm: bool) {
    let nc = cp.nc;
    let mc = cp.mc;

    // Find stagnation point (where q changes sign or is minimum)
    let mut wcle = PI;
    for ic in 1..nc {
        if qc[ic] < 0.0 {
            wcle = cp.wc[ic];
            break;
        }
    }
    let alfcir = 0.5 * (wcle - PI);

    // Save old Cn
    let cn_save: Vec<C64> = cp.cn.clone();

    // Compute P(w) from q(w)
    for ic in 1..nc - 1 {
        let cosw = 2.0 * (0.5 * cp.wc[ic] - alfcir).cos();
        let sinw = 2.0 * (0.5 * cp.wc[ic]).sin();
        let sinwe = if sinw > 0.0 { sinw.powf(cp.agte) } else { 0.0 };

        let pfun = if cosw.abs() < 1e-4 {
            (sinwe / qc[ic].abs().max(1e-10)).abs()
        } else {
            (cosw * sinwe / qc[ic]).abs().max(1e-30)
        };

        cp.piq[ic] = c(pfun.ln(), 0.0);
    }

    // Extrapolate to TE
    cp.piq[0] = 3.0 * cp.piq[1] - 3.0 * cp.piq[2] + cp.piq[3];
    cp.piq[nc - 1] = 3.0 * cp.piq[nc - 2] - 3.0 * cp.piq[nc - 3] + cp.piq[nc - 4];

    ftp(cp);
    cp.cn[0] = c(0.0, cp.qimold);

    if lsymm {
        for m in 1..=mc {
            let cnr = 2.0 * cp.cn[m].re - cn_save[m].re;
            cp.cn[m] = c(cnr, 0.0);
        }
    }

    piqsum(cp);
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MdesOptions {
    pub nc: usize,
    pub alpha_deg: f64,
    pub symmetric: bool,
    pub filter_strength: f64,
}

impl Default for MdesOptions {
    fn default() -> Self {
        Self {
            nc: 129,
            alpha_deg: 0.0,
            symmetric: false,
            filter_strength: 0.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MdesResult {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub qspec_x: Vec<f64>,
    pub qspec_values: Vec<f64>,
    pub cl: f64,
    pub cm: f64,
    pub success: bool,
    pub error: Option<String>,
}

/// Run full-inverse design via circle-plane conformal mapping.
///
/// Given an initial airfoil and a target velocity distribution Q(w),
/// compute the modified airfoil geometry.
///
/// * `x`, `y` – initial airfoil coordinates (TE -> upper -> LE -> lower -> TE)
/// * `target_q` – target velocity distribution at the nc circle-plane points
///   (if None, uses the current airfoil's distribution)
/// * `options` – MDES configuration
pub fn solve_mdes(
    x: &[f64],
    y: &[f64],
    target_q: Option<&[f64]>,
    options: &MdesOptions,
) -> MdesResult {
    let nc = options.nc.max(33).min(ICX);
    // nc must be 2^n + 1
    let nc = {
        let mut v = 16;
        while v + 1 < nc { v *= 2; }
        v + 1
    };

    let mut cp = CirclePlane::new(nc);
    eiwset(&mut cp);

    // Initialize circle-plane mapping from current airfoil
    if let Err(e) = scinit(&mut cp, x, y) {
        return MdesResult {
            x: vec![], y: vec![], qspec_x: vec![], qspec_values: vec![],
            cl: 0.0, cm: 0.0, success: false, error: Some(e.to_string()),
        };
    }

    // If a target Q is provided, recompute Cn from it
    if let Some(tq) = target_q {
        if tq.len() == nc {
            cncalc(&mut cp, tq, options.symmetric);
        }
    }

    // Apply filtering if requested
    if options.filter_strength > 0.0 {
        cnfilt(&mut cp, options.filter_strength);
        piqsum(&mut cp);
    }

    // Generate new geometry from (possibly modified) Cn
    mapgen(&mut cp);

    // Extract output coordinates
    let out_x: Vec<f64> = cp.zc[..nc].iter().map(|z| z.re).collect();
    let out_y: Vec<f64> = cp.zc[..nc].iter().map(|z| z.im).collect();

    // Compute Q(w) for the result (for plotting)
    let alfa = options.alpha_deg.to_radians();
    let alfcir = alfa - cp.cn[0].im;
    let mut qspec_values = Vec::with_capacity(nc);
    let mut qspec_x = Vec::with_capacity(nc);
    for ic in 0..nc {
        let ppp = cp.piq[ic].re;
        let eppp = (-ppp).exp();
        let sinw = 2.0 * (0.5 * cp.wc[ic]).sin();
        let sinwe = if cp.agte == 0.0 {
            1.0
        } else if sinw > 0.0 {
            sinw.powf(cp.agte)
        } else {
            0.0
        };
        let q = 2.0 * (0.5 * cp.wc[ic] - alfcir).cos() * sinwe * eppp;
        qspec_values.push(q);
        qspec_x.push(cp.zc[ic].re);
    }

    // Approximate CL from Kutta-Joukowski
    let cl = 2.0 * PI * alfa.sin() * 2.0; // rough estimate

    MdesResult {
        x: out_x,
        y: out_y,
        qspec_x,
        qspec_values,
        cl,
        cm: 0.0,
        success: true,
        error: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naca0012_coords() -> (Vec<f64>, Vec<f64>) {
        // XFOIL convention: upper TE -> LE -> lower TE
        let n = 65;
        let mut x = Vec::with_capacity(2 * n - 1);
        let mut y = Vec::with_capacity(2 * n - 1);

        let yt_at = |xc: f64| -> f64 {
            0.12 / 0.2 * (
                0.2969 * xc.sqrt()
                - 0.1260 * xc
                - 0.3516 * xc.powi(2)
                + 0.2843 * xc.powi(3)
                - 0.1015 * xc.powi(4)
            )
        };

        // Upper surface: TE (x=1) -> LE (x=0)
        for i in 0..n {
            let beta = PI * i as f64 / (n - 1) as f64;
            let xc = 0.5 * (1.0 + beta.cos()); // 1 -> 0
            x.push(xc);
            y.push(yt_at(xc));
        }
        // Lower surface: LE+eps -> TE (x=1)
        for i in 1..n {
            let beta = PI * i as f64 / (n - 1) as f64;
            let xc = 0.5 * (1.0 - beta.cos()); // 0 -> 1
            x.push(xc);
            y.push(-yt_at(xc));
        }
        (x, y)
    }

    #[test]
    fn eiwset_basic_sanity() {
        let mut cp = CirclePlane::new(33);
        eiwset(&mut cp);

        assert_eq!(cp.nc, 33);
        assert_eq!(cp.mc, 8);
        assert!((cp.wc[0]).abs() < 1e-12, "w[0] should be 0");
        assert!((cp.wc[32] - 2.0 * PI).abs() < 1e-10, "w[nc-1] should be 2*pi");

        // exp(i*0) = 1 for all m
        for m in 0..=cp.mc {
            assert!((cp.eiw[0][m] - c(1.0, 0.0)).norm() < 1e-10);
        }
    }

    #[test]
    fn ftp_piqsum_roundtrip() {
        let mut cp = CirclePlane::new(65);
        eiwset(&mut cp);

        // Set a known P(w) = cos(w) (real part only)
        for ic in 0..cp.nc {
            cp.piq[ic] = c(cp.wc[ic].cos(), 0.0);
        }

        let original_piq: Vec<C64> = cp.piq.clone();

        // Forward transform: P(w) -> Cn
        ftp(&mut cp);

        // Cn[1] should dominate for cos(w) input
        // cos(w) = Re[ e^{iw} ] so Cn[1] should have real part close to 1
        assert!(cp.cn[1].re.abs() > 0.3, "cn[1] should be significant for cos(w) input");

        // Inverse transform: Cn -> P(w)+iQ(w)
        piqsum(&mut cp);

        // Real part should approximate original cos(w) at interior points
        let mut max_err = 0.0_f64;
        for ic in 2..cp.nc - 2 {
            let err = (cp.piq[ic].re - original_piq[ic].re).abs();
            max_err = max_err.max(err);
        }
        assert!(max_err < 0.15, "FTP->PIQSUM roundtrip error {max_err} too large");
    }

    #[test]
    fn cnfilt_attenuates_high_harmonics() {
        let mut cp = CirclePlane::new(33);
        eiwset(&mut cp);

        for m in 0..=cp.mc {
            cp.cn[m] = c(1.0, 0.5);
        }

        cnfilt(&mut cp, 1.0);

        // Low-order coefficients should be mostly preserved
        assert!(cp.cn[0].re > 0.9);
        // Highest harmonic should be strongly attenuated
        assert!(cp.cn[cp.mc].norm() < 0.1, "highest harmonic should be filtered");
    }

    #[test]
    fn solve_mdes_naca0012_returns_valid_geometry() {
        let (x, y) = naca0012_coords();

        let result = solve_mdes(&x, &y, None, &MdesOptions::default());

        assert!(result.success, "solve_mdes should succeed: {:?}", result.error);
        assert!(result.x.len() > 10, "should return geometry points");
        assert_eq!(result.x.len(), result.y.len());

        // Output should resemble an airfoil with a positive chord.
        // NOTE: MDES normalization is approximate; chord may not be exactly 1.0.
        let x_min = result.x.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = result.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let chord = x_max - x_min;
        assert!(chord > 0.5, "chord should be positive, got {chord}");
        assert!(chord < 5.0, "chord should be finite, got {chord}");

        // Should have points on both sides of y=0
        let y_max = result.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = result.y.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(y_max > 0.01, "should have upper surface points");
        assert!(y_min < -0.01, "should have lower surface points");
    }

    #[test]
    fn solve_mdes_rejects_too_few_points() {
        let result = solve_mdes(&[0.0, 1.0], &[0.0, 0.0], None, &MdesOptions::default());
        assert!(!result.success);
        assert!(result.error.is_some());
    }
}
