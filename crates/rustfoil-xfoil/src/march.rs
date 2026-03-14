//! Faithful XFOIL boundary layer marching.
//!
//! Implements BLPINI, MRCHUE (xbl.f:584-986), and MRCHDU (xbl.f:990-1312).

use rustfoil_bl::{
    bldif, blvar,
    closures::{hkin, trchek2_full, Trchek2FullResult},
    equations::{bldif_full_simi, trdif_full},
    BlStation, FlowType,
};

use crate::state::{XfoilBlRow, XfoilState};

// Shape parameter limits for direct→inverse mode switching (xbl.f:609-610)
const HLMAX: f64 = 3.8;
const HTMAX: f64 = 2.5;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// 4×4 Gaussian elimination with partial pivoting.
fn gauss_4x4(a: &mut [[f64; 4]; 4], b: &mut [f64; 4]) {
    for col in 0..4 {
        // Partial pivoting
        let mut max_val = a[col][col].abs();
        let mut max_row = col;
        for row in (col + 1)..4 {
            let v = a[row][col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_row != col {
            a.swap(col, max_row);
            b.swap(col, max_row);
        }
        let pivot = a[col][col];
        if pivot.abs() < 1e-30 {
            continue;
        }
        // Eliminate below
        for row in (col + 1)..4 {
            let factor = a[row][col] / pivot;
            for k in col..4 {
                a[row][k] -= factor * a[col][k];
            }
            b[row] -= factor * b[col];
        }
    }
    // Back-substitute
    for col in (0..4).rev() {
        if a[col][col].abs() < 1e-30 {
            b[col] = 0.0;
            continue;
        }
        for k in (col + 1)..4 {
            b[col] -= a[col][k] * b[k];
        }
        b[col] /= a[col][col];
    }
}

/// DSLIM: limit δ* so that Hk ≥ hklim.  (xbl.f:1771-1781)
fn dslim(dstr: &mut f64, thet: f64, msq: f64, hklim: f64) {
    let h = *dstr / thet.max(1e-12);
    let hk_res = hkin(h, msq);
    let dh = (hklim - hk_res.hk).max(0.0) / hk_res.hk_h.max(1e-12);
    *dstr += dh * thet;
}

/// Build a BlStation from marching variables and call blvar.
fn build_and_fill_station(
    x: f64,
    u: f64,
    theta: f64,
    delta_star: f64,
    dw: f64,
    ctau: f64,
    ampl: f64,
    flow_type: FlowType,
    msq: f64,
    re: f64,
) -> BlStation {
    let mut s = BlStation::new();
    s.x = x;
    s.u = u.abs().max(1e-6);
    s.theta = theta.max(1e-12);
    s.delta_star = delta_star.max(1e-12);
    s.dw = dw;
    s.ctau = ctau;
    s.ampl = ampl;
    s.h = s.delta_star / s.theta;
    // XFOIL's wake mass state is carried on total DSTR = delta_star + dw.
    // BLPRV/BLKIN still use the wake-free delta_star for local closure work,
    // but the coupled/marched mass contract must include the wake gap.
    let total_dstr = if matches!(flow_type, FlowType::Wake) {
        s.delta_star + s.dw
    } else {
        s.delta_star
    };
    s.mass_defect = s.u * total_dstr;
    s.is_laminar = matches!(flow_type, FlowType::Laminar);
    s.is_turbulent = matches!(flow_type, FlowType::Turbulent) || matches!(flow_type, FlowType::Wake);
    s.is_wake = matches!(flow_type, FlowType::Wake);
    blvar(&mut s, flow_type, msq, re);
    s
}

/// Determine flow type for a station.
fn xfoil_ibl(ibl: usize) -> usize {
    ibl + 1
}

fn is_transition_or_turbulent(ibl: usize, itran: usize) -> bool {
    xfoil_ibl(ibl) >= itran
}

fn is_laminar_interval(ibl: usize, itran: usize) -> bool {
    xfoil_ibl(ibl) < itran
}

fn is_transition_interval(ibl: usize, itran: usize) -> bool {
    xfoil_ibl(ibl) == itran
}

fn flow_type_for(ibl: usize, itran: usize, wake: bool) -> FlowType {
    if wake {
        FlowType::Wake
    } else if is_transition_or_turbulent(ibl, itran) {
        FlowType::Turbulent
    } else {
        FlowType::Laminar
    }
}

/// Store converged station variables back to a row.
fn store_to_row(
    station: &BlStation,
    row: &mut XfoilBlRow,
    dsi: f64,
    uei: f64,
    cti: f64,
    ami: f64,
    ibl: usize,
    itran: usize,
    wake: bool,
) {
    let is_turb = is_transition_or_turbulent(ibl, itran) || wake;
    let dw = if wake { station.dw.max(0.0) } else { 0.0 };
    // Primary
    row.theta = station.theta;
    row.dstr = dsi;
    row.dw = dw;
    row.uedg = uei;
    row.mass = dsi * uei;
    // CTAU stores AMI for laminar, CTI for turbulent (matching Fortran convention)
    row.ctau = if is_turb { cti } else { ami };
    row.ampl = ami;
    // Secondary from blvar
    row.hk = station.hk;
    row.h = station.h;
    row.hs = station.hs;
    row.hc = station.hc;
    row.r_theta = station.r_theta;
    row.cf = station.cf;
    row.cd = station.cd;
    row.us = station.us;
    row.cq = station.cq;
    row.de = station.de;
    row.derivs = station.derivs.clone();
    // Flags
    row.is_laminar = !is_turb && !wake;
    row.is_turbulent = is_turb;
    row.is_wake = wake;
}

/// Compute HTARG for inverse mode (xbl.f:785-813).
fn compute_htarg(
    hk1: f64,
    x1: f64,
    x2: f64,
    t1: f64,
    xt: f64,
    ibl: usize,
    itran: usize,
    wake: bool,
    hmax: f64,
) -> f64 {
    let dx = (x2 - x1).max(1e-12);
    let htarg = if is_laminar_interval(ibl, itran) {
        // Laminar: slow increase
        hk1 + 0.03 * dx / t1.max(1e-12)
    } else if is_transition_interval(ibl, itran) {
        // Transition interval: weighted
        hk1 + (0.03 * (xt - x1) - 0.15 * (x2 - xt)) / t1.max(1e-12)
    } else if wake {
        // Wake: asymptotic with 3 Newton iterations
        let cnst = 0.03 * dx / t1.max(1e-12);
        let mut hk2 = hk1;
        for _ in 0..3 {
            let hkm1 = hk2 - 1.0;
            hk2 -= (hk2 + cnst * hkm1 * hkm1 * hkm1 - hk1)
                / (1.0 + 3.0 * cnst * hkm1 * hkm1);
        }
        hk2
    } else {
        // Turbulent: fast decrease
        hk1 - 0.15 * dx / t1.max(1e-12)
    };
    if wake {
        htarg.max(1.01)
    } else {
        htarg.max(hmax)
    }
}

// ---------------------------------------------------------------------------
// BLPINI
// ---------------------------------------------------------------------------

/// BLPINI: Initialize BL closure constants.
///
/// In XFOIL this sets SCCON, GACON, etc. (xbl.f:1785-1804).
/// These constants are already compiled into `rustfoil-bl/src/constants.rs`.
pub fn blpini(state: &mut XfoilState, _reynolds: f64) {
    let _ = state;
}

// ---------------------------------------------------------------------------
// MRCHUE
// ---------------------------------------------------------------------------

/// MRCHUE: Direct-mode BL march with transition checking (xbl.f:584-986).
pub fn mrchue(state: &mut XfoilState, reynolds: f64, ncrit: f64, iteration: usize) {
    let msq = 0.0_f64;

    // Upper surface (IS=1): other side = lower (stale values)
    let other_te = te_values(&state.lower_rows, state.iblte_lower);
    let (itran_u, xssitr_u) = march_ue_side(
        1,
        &mut state.upper_rows, state.nbl_upper, state.iblte_upper,
        other_te, state.ante, reynolds, ncrit, msq, iteration,
    );
    state.itran_upper = itran_u;
    state.xssitr_upper = xssitr_u;

    // Lower surface (IS=2): other side = upper (fresh values)
    let other_te = te_values(&state.upper_rows, state.iblte_upper);
    let (itran_l, xssitr_l) = march_ue_side(
        2,
        &mut state.lower_rows, state.nbl_lower, state.iblte_lower,
        other_te, state.ante, reynolds, ncrit, msq, iteration,
    );
    state.itran_lower = itran_l;
    state.xssitr_lower = xssitr_l;

    state.lblini = true;
}

/// TE values needed for merging.
struct TeValues {
    theta: f64,
    dstr: f64,
    ctau: f64,
}

fn te_values(rows: &[XfoilBlRow], iblte: usize) -> TeValues {
    rows.get(iblte).map_or(
        TeValues { theta: 0.0, dstr: 0.0, ctau: 0.03 },
        |r| TeValues { theta: r.theta, dstr: r.dstr, ctau: r.ctau },
    )
}

/// Per-side MRCHUE march. Returns `(itran, xssitr)`.
fn march_ue_side(
    side: usize,
    rows: &mut [XfoilBlRow],
    nbl: usize,
    iblte: usize,
    other_te: TeValues,
    ante: f64,
    re: f64,
    ncrit: f64,
    msq: f64,
    iteration: usize,
) -> (usize, f64) {
    if nbl < 2 {
        return (iblte, 0.0);
    }

    let mut itran = xfoil_ibl(iblte);
    let mut xssitr = 0.0;
    let mut turb = false;
    let x_forced = rows.get(iblte).map(|row| row.x);

    // ---- Similarity station init (IBL=2 = index 1) ----
    let xsi0 = rows[1].x.max(1e-12);
    let uei0 = rows[1].uedg.abs().max(1e-6);
    let tsq = 0.45 * xsi0 / (6.0 * uei0 * re);
    let mut thi = tsq.max(0.0).sqrt().max(1e-12);
    let mut dsi = 2.2 * thi;
    let mut ami = 0.0_f64;
    let mut cti = 0.03_f64;

    // Station_1 for the similarity station: will be set properly after IBL=2 converges.
    // For the similarity station itself, station_1 = station_2 (SIMI flag).
    let mut station_1 = build_and_fill_station(
        xsi0, uei0, thi, dsi, 0.0, cti, ami, FlowType::Laminar, msq, re,
    );

    // ---- March downstream ----
    for ibl in 1..nbl {
        let simi = ibl == 1;
        let wake = ibl > iblte;
        let xsi = rows[ibl].x;
        let mut uei = rows[ibl].uedg;
        let dswaki = if wake { rows[ibl].dw } else { 0.0 };

        let mut direct = true;
        let mut htarg = 0.0_f64;
        let mut dmax = 0.0_f64;
        let mut tran = false;
        let mut xt = 0.0_f64;
        let mut tr_result: Option<Trchek2FullResult> = None;

        // ---- Newton loop (up to 25 iterations) ----
        for _itbl in 0..25 {
            // Determine flow type (may change if transition detected within loop)
            let ft = flow_type_for(ibl, itran, wake);

            // Build station_2 with full secondary variables
            let mut station_2 = build_and_fill_station(
                xsi, uei, thi, (dsi - dswaki).max(1e-12), dswaki, cti, ami, ft, msq, re,
            );

            // ---- Transition check (laminar, not similarity) ----
            if !simi && !turb {
                let tr = trchek2_full(
                    station_1.x, station_2.x,
                    station_1.theta, station_2.theta,
                    station_1.delta_star, station_2.delta_star,
                    station_1.u, station_2.u,
                    station_1.hk, station_2.hk,
                    station_1.r_theta, station_2.r_theta,
                    station_1.ampl, ncrit,
                    x_forced, msq, re,
                );
                if rustfoil_bl::is_debug_active() {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_final(
                        iteration,
                        side,
                        ibl + 1,
                        tr.iterations,
                        tr.converged,
                        station_1.x,
                        station_2.x,
                        station_1.ampl,
                        tr.ampl2,
                        tr.ax,
                        tr.xt,
                        ncrit,
                        tr.transition,
                        tr.forced,
                        Some(tr.wf1),
                        Some(tr.wf2),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                    ));
                }
                ami = tr.ampl2;
                if tr.transition {
                    itran = xfoil_ibl(ibl);
                    tran = true;
                    xt = tr.xt;
                    if cti <= 0.0 {
                        cti = 0.03;
                    }
                    tr_result = Some(tr);
                } else {
                    itran = xfoil_ibl(ibl) + 2;
                    tran = false;
                    tr_result = None;
                }
                // Rebuild station_2 if flow type changed
                let ft_new = flow_type_for(ibl, itran, wake);
                if ft != ft_new {
                    station_2 = build_and_fill_station(
                        xsi, uei, thi, (dsi - dswaki).max(1e-12), dswaki, cti, ami,
                        ft_new, msq, re,
                    );
                }
            }

            let ft = flow_type_for(ibl, itran, wake);

            // ---- Build 3×4 BL system ----
            let mut sys = [[0.0_f64; 4]; 4];
            let mut rhs = [0.0_f64; 4];

            if ibl == iblte + 1 {
                // TESYS: TE dummy system (xblsys.f:682-716)
                let my_te_theta = rows[iblte].theta;
                let my_te_dstr = rows[iblte].dstr;
                let my_te_ctau = rows[iblte].ctau;
                let tte = my_te_theta + other_te.theta;
                let dte = my_te_dstr + other_te.dstr + ante;
                let cte = if tte > 0.0 {
                    (my_te_ctau * my_te_theta + other_te.ctau * other_te.theta) / tte
                } else {
                    0.03
                };
                // Rebuild station_2 as wake for TESYS
                station_2 = build_and_fill_station(
                    xsi, uei, thi, (dsi - dswaki).max(1e-12), dswaki, cti, ami,
                    FlowType::Wake, msq, re,
                );
                // Identity equations
                sys[0] = [1.0, 0.0, 0.0, 0.0];
                rhs[0] = cte - station_2.ctau;
                sys[1] = [0.0, 1.0, 0.0, 0.0];
                rhs[1] = tte - station_2.theta;
                sys[2] = [0.0, 0.0, 1.0, 0.0];
                rhs[2] = dte - station_2.delta_star - station_2.dw;
            } else {
                let (res, jac) = if simi {
                    bldif_full_simi(&station_2, FlowType::Laminar, msq, re)
                } else if tran && tr_result.is_some() {
                    trdif_full(
                        &station_1,
                        &station_2,
                        tr_result.as_ref().unwrap(),
                        ncrit,
                        msq,
                        re,
                    )
                } else {
                    bldif(&station_1, &station_2, ft, msq, re)
                };
                for eq in 0..3 {
                    for var in 0..4 {
                        sys[eq][var] = jac.vs2[eq][var];
                    }
                }
                rhs[0] = res.res_third;
                rhs[1] = res.res_mom;
                rhs[2] = res.res_shape;
            }

            // ---- 4th equation + solve ----
            if direct {
                // Direct mode: Ue prescribed
                sys[3] = [0.0, 0.0, 0.0, 1.0];
                rhs[3] = 0.0;
                gauss_4x4(&mut sys, &mut rhs);

                // Compute DMAX and RLX
                dmax = (rhs[1] / thi.max(1e-12)).abs().max((rhs[2] / dsi.max(1e-12)).abs());
                if is_laminar_interval(ibl, itran) {
                    dmax = dmax.max((rhs[0] / 10.0).abs());
                }
                if is_transition_or_turbulent(ibl, itran) {
                    dmax = dmax.max((rhs[0] / cti.max(1e-12)).abs());
                }
                let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

                // Check Hk for direct→inverse switch (skip TE+1)
                if ibl != iblte + 1 {
                    let htest = (dsi + rlx * rhs[2]) / (thi + rlx * rhs[1]).max(1e-12);
                    let hktest = hkin(htest, msq).hk;
                    let hmax = if is_laminar_interval(ibl, itran) { HLMAX } else { HTMAX };
                    if hktest >= hmax {
                        direct = false;
                        htarg = compute_htarg(
                            station_1.hk, station_1.x, xsi, station_1.theta,
                            xt, ibl, itran, wake, hmax,
                        );
                        continue; // Redo iteration in inverse mode
                    }
                }

                // Apply direct updates (AMI not updated in MRCHUE — commented out in Fortran)
                if is_transition_or_turbulent(ibl, itran) {
                    cti += rlx * rhs[0];
                }
                thi += rlx * rhs[1];
                dsi += rlx * rhs[2];
            } else {
                // Inverse mode: Hk prescribed
                let hk2_t2 = station_2.derivs.hk_h * station_2.derivs.h_theta;
                let hk2_d2 = station_2.derivs.hk_h * station_2.derivs.h_delta_star;
                sys[3] = [0.0, hk2_t2, hk2_d2, 0.0];
                rhs[3] = htarg - station_2.hk;
                gauss_4x4(&mut sys, &mut rhs);

                // Compute DMAX (inverse includes Ue)
                dmax = (rhs[1] / thi.max(1e-12))
                    .abs()
                    .max((rhs[2] / dsi.max(1e-12)).abs())
                    .max((rhs[3] / uei.abs().max(1e-6)).abs());
                if is_transition_or_turbulent(ibl, itran) {
                    dmax = dmax.max((rhs[0] / cti.max(1e-12)).abs());
                }
                let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

                // Apply inverse updates
                if is_transition_or_turbulent(ibl, itran) {
                    cti += rlx * rhs[0];
                }
                thi += rlx * rhs[1];
                dsi += rlx * rhs[2];
                uei += rlx * rhs[3];
            }

            // ---- Clamp CTI ----
            if is_transition_or_turbulent(ibl, itran) {
                cti = cti.clamp(1e-7, 0.30);
            }

            // ---- DSLIM ----
            let hklim = if ibl <= iblte { 1.02 } else { 1.00005 };
            let mut dsw = dsi - dswaki;
            dslim(&mut dsw, thi, msq, hklim);
            dsi = dsw + dswaki;

            if dmax <= 1e-5 {
                break;
            }
        }

        let mut fallback_station_final: Option<BlStation> = None;

        // ---- Convergence failure fallback (extrapolate if DMAX > 0.1) ----
        if dmax > 0.1 && ibl > 2 {
            if ibl <= iblte {
                let prev_x = rows[ibl - 1].x.max(1e-12);
                let ratio = (xsi / prev_x).sqrt();
                thi = rows[ibl - 1].theta * ratio;
                dsi = rows[ibl - 1].dstr * ratio;
            } else if ibl == iblte + 1 {
                let my_te_theta = rows[iblte].theta;
                let my_te_dstr = rows[iblte].dstr;
                cti = if (my_te_theta + other_te.theta) > 0.0 {
                    (rows[iblte].ctau * my_te_theta + other_te.ctau * other_te.theta)
                        / (my_te_theta + other_te.theta)
                } else {
                    0.03
                };
                thi = my_te_theta + other_te.theta;
                dsi = my_te_dstr + other_te.dstr + ante;
            } else {
                thi = rows[ibl - 1].theta;
                let ratlen =
                    (xsi - rows[ibl - 1].x) / (10.0 * rows[ibl - 1].dstr.max(1e-12));
                dsi = (rows[ibl - 1].dstr + thi * ratlen) / (1.0 + ratlen);
            }
            if is_transition_interval(ibl, itran) {
                cti = 0.05;
            }
            if xfoil_ibl(ibl) > itran {
                cti = rows[ibl - 1].ctau;
            }
            uei = rows[ibl].uedg;
            if ibl + 1 < nbl {
                uei = 0.5 * (rows[ibl - 1].uedg + rows[ibl + 1].uedg);
            }

            // Match XFOIL label 109 path: rebuild the extrapolated station,
            // re-run transition detection, then recompute the station using
            // the updated laminar/turbulent classification before storing.
            let mut fallback_station = build_and_fill_station(
                xsi,
                uei,
                thi,
                (dsi - dswaki).max(1e-12),
                dswaki,
                cti,
                ami,
                flow_type_for(ibl, itran, wake),
                msq,
                re,
            );
            if !simi && !turb {
                let tr = trchek2_full(
                    station_1.x,
                    fallback_station.x,
                    station_1.theta,
                    fallback_station.theta,
                    station_1.delta_star,
                    fallback_station.delta_star,
                    station_1.u,
                    fallback_station.u,
                    station_1.hk,
                    fallback_station.hk,
                    station_1.r_theta,
                    fallback_station.r_theta,
                    station_1.ampl,
                    ncrit,
                    x_forced,
                    msq,
                    re,
                );
                if rustfoil_bl::is_debug_active() {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_final(
                        iteration,
                        side,
                        ibl + 1,
                        tr.iterations,
                        tr.converged,
                        station_1.x,
                        fallback_station.x,
                        station_1.ampl,
                        tr.ampl2,
                        tr.ax,
                        tr.xt,
                        ncrit,
                        tr.transition,
                        tr.forced,
                        Some(tr.wf1),
                        Some(tr.wf2),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                    ));
                }
                ami = tr.ampl2;
                if tr.transition {
                    itran = xfoil_ibl(ibl);
                    tran = true;
                    xt = tr.xt;
                    tr_result = Some(tr);
                } else {
                    itran = xfoil_ibl(ibl) + 2;
                    tran = false;
                    tr_result = None;
                }
                fallback_station = build_and_fill_station(
                    xsi,
                    uei,
                    thi,
                    (dsi - dswaki).max(1e-12),
                    dswaki,
                    cti,
                    ami,
                    flow_type_for(ibl, itran, wake),
                    msq,
                    re,
                );
            }
            fallback_station_final = Some(fallback_station);
        }

        // ---- Store converged state (xbl.f:940-951) ----
        let station_final = fallback_station_final.unwrap_or_else(|| {
            build_and_fill_station(
                xsi,
                uei,
                thi,
                (dsi - dswaki).max(1e-12),
                dswaki,
                cti,
                ami,
                flow_type_for(ibl, itran, wake),
                msq,
                re,
            )
        });
        if std::env::var("RUSTFOIL_WAKE_MARCH_DEBUG").is_ok()
            && ((wake && (iblte + 1..=iblte + 3).contains(&ibl)) || ibl == iblte)
        {
            eprintln!(
                "[WAKE MARCH] ibl={} x={:.8e} u={:.8e} theta={:.8e} dsi_total={:.8e} delta_star={:.8e} dw={:.8e} ctau={:.8e} ampl={:.8e} wake={}",
                ibl,
                xsi,
                uei,
                thi,
                dsi,
                station_final.delta_star,
                dswaki,
                cti,
                ami,
                wake
            );
        }
        store_to_row(&station_final, &mut rows[ibl], dsi, uei, cti, ami, ibl, itran, wake);
        if rustfoil_bl::is_debug_active() {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue(
                side,
                ibl + 1,
                station_final.x,
                station_final.u,
                station_final.theta,
                station_final.delta_star,
                station_final.hk,
                station_final.cf,
            ));
        }

        // ---- Set station_1 for next station (xbl.f:953-965) ----
        station_1 = station_final;

        // ---- Transition / TE bookkeeping (xbl.f:968-981) ----
        if tran {
            xssitr = xt;
        }
        if tran || ibl == iblte {
            turb = true;
        }

        // TE merge: set THI, DSI for wake entry (xbl.f:978-981)
        if ibl == iblte {
            thi = rows[iblte].theta + other_te.theta;
            dsi = rows[iblte].dstr + other_te.dstr + ante;
        }
    }

    (itran, xssitr)
}

// ---------------------------------------------------------------------------
// MRCHDU
// ---------------------------------------------------------------------------

/// MRCHDU: Mixed Ue-Hk march with transition checking (xbl.f:990-1312).
pub fn mrchdu(state: &mut XfoilState, reynolds: f64, ncrit: f64, iteration: usize) {
    let msq = 0.0_f64;

    // Upper surface
    let other_te = te_values(&state.lower_rows, state.iblte_lower);
    let itrold_u = state.itran_upper;
    let (itran_u, xssitr_u) = march_du_side(
        1,
        &mut state.upper_rows, state.nbl_upper, state.iblte_upper,
        other_te, state.ante, itrold_u, reynolds, ncrit, msq, iteration,
    );
    state.itran_upper = itran_u;
    state.xssitr_upper = xssitr_u;

    // Lower surface (uses freshly updated upper TE)
    let other_te = te_values(&state.upper_rows, state.iblte_upper);
    let itrold_l = state.itran_lower;
    let (itran_l, xssitr_l) = march_du_side(
        2,
        &mut state.lower_rows, state.nbl_lower, state.iblte_lower,
        other_te, state.ante, itrold_l, reynolds, ncrit, msq, iteration,
    );
    state.itran_lower = itran_l;
    state.xssitr_lower = xssitr_l;
}

/// Per-side MRCHDU march. Returns `(itran, xssitr)`.
#[allow(clippy::too_many_arguments)]
fn march_du_side(
    side: usize,
    rows: &mut [XfoilBlRow],
    nbl: usize,
    iblte: usize,
    other_te: TeValues,
    ante: f64,
    itrold: usize,
    re: f64,
    ncrit: f64,
    msq: f64,
    iteration: usize,
) -> (usize, f64) {
    if nbl < 2 {
        return (iblte, 0.0);
    }

    const DEPS: f64 = 5.0e-6;
    const SENSWT: f64 = 1000.0;

    let mut itran = xfoil_ibl(iblte);
    let mut xssitr = 0.0;
    let mut turb = false;
    let mut sens = 0.0_f64;
    let mut ami = rows[1].ctau;
    let x_forced = rows.get(iblte).map(|row| row.x);

    // Build station_1 from stored similarity station (index 1)
    let ft0 = flow_type_for(1, itrold, false);
    let dswaki0 = rows[1].dw;
    let mut station_1 = build_and_fill_station(
        rows[1].x,
        rows[1].uedg,
        rows[1].theta,
        (rows[1].dstr - dswaki0).max(1e-12),
        dswaki0,
        rows[1].ctau,
        rows[1].ampl,
        ft0,
        msq,
        re,
    );

    for ibl in 1..nbl {
        let simi = ibl == 1;
        let wake = ibl > iblte;
        let xsi = rows[ibl].x;
        let mut uei = rows[ibl].uedg;
        let mut thi = rows[ibl].theta;
        let mut dsi = rows[ibl].dstr;
        let dswaki = if wake { rows[ibl].dw } else { 0.0 };

        // Initialize AMI/CTI from stored CTAU based on ITROLD (xbl.f:1050-1056)
        let mut cti;
        if is_laminar_interval(ibl, itrold) {
            ami = rows[ibl].ctau;
            cti = 0.03;
        } else {
            cti = rows[ibl].ctau;
            if cti <= 0.0 {
                cti = 0.03;
            }
        }

        // Clamp initial DSI (xbl.f:1067-1068)
        if ibl <= iblte {
            dsi = (dsi - dswaki).max(1.02 * thi) + dswaki;
        } else {
            dsi = (dsi - dswaki).max(1.00005 * thi) + dswaki;
        }

        let mut dmax = 0.0_f64;
        let mut tran = false;
        let mut xt = 0.0_f64;
        let mut tr_result: Option<Trchek2FullResult> = None;
        let mut ueref = 0.0_f64;
        let mut hkref = 0.0_f64;
        let mut sennew = 0.0_f64;

        // ---- Newton iteration loop ----
        for itbl in 0..25 {
            let ft = flow_type_for(ibl, itran, wake);

            let mut station_2 = build_and_fill_station(
                xsi, uei, thi, (dsi - dswaki).max(1e-12), dswaki, cti, ami, ft, msq, re,
            );

            // Transition check
            if !simi && !turb {
                let tr = trchek2_full(
                    station_1.x, station_2.x,
                    station_1.theta, station_2.theta,
                    station_1.delta_star, station_2.delta_star,
                    station_1.u, station_2.u,
                    station_1.hk, station_2.hk,
                    station_1.r_theta, station_2.r_theta,
                    station_1.ampl, ncrit,
                    x_forced, msq, re,
                );
                if rustfoil_bl::is_debug_active() {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_final(
                        iteration,
                        side,
                        ibl + 1,
                        tr.iterations,
                        tr.converged,
                        station_1.x,
                        station_2.x,
                        station_1.ampl,
                        tr.ampl2,
                        tr.ax,
                        tr.xt,
                        ncrit,
                        tr.transition,
                        tr.forced,
                        Some(tr.wf1),
                        Some(tr.wf2),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                    ));
                }
                if std::env::var("RUSTFOIL_TRCHEK_DEBUG").is_ok() {
                    eprintln!(
                        "[TRCHEK DU] side={} ibl={} itran={} itrold={} x1={:.12e} x2={:.12e} th1={:.12e} th2={:.12e} ds1={:.12e} ds2={:.12e} u1={:.12e} u2={:.12e} hk1={:.12e} hk2={:.12e} rt1={:.12e} rt2={:.12e} ampl1={:.12e} ampl2={:.12e} transition={} xt={:.12e} xforced={:?}",
                        side,
                        ibl,
                        itran,
                        itrold,
                        station_1.x,
                        station_2.x,
                        station_1.theta,
                        station_2.theta,
                        station_1.delta_star,
                        station_2.delta_star,
                        station_1.u,
                        station_2.u,
                        station_1.hk,
                        station_2.hk,
                        station_1.r_theta,
                        station_2.r_theta,
                        station_1.ampl,
                        tr.ampl2,
                        tr.transition,
                        tr.xt,
                        x_forced,
                    );
                }
                ami = tr.ampl2;
                if tr.transition {
                    itran = xfoil_ibl(ibl);
                    tran = true;
                    xt = tr.xt;
                    tr_result = Some(tr);
                } else {
                    itran = xfoil_ibl(ibl) + 2;
                    tran = false;
                    tr_result = None;
                }
                let ft_new = flow_type_for(ibl, itran, wake);
                station_2 = build_and_fill_station(
                    xsi, uei, thi, (dsi - dswaki).max(1e-12), dswaki, cti, ami,
                    ft_new, msq, re,
                );
            }

            let ft = flow_type_for(ibl, itran, wake);

            // ---- Build 3×4 system ----
            let mut sys = [[0.0_f64; 4]; 4];
            let mut rhs = [0.0_f64; 4];

            if ibl == iblte + 1 {
                // TESYS
                let my_te_theta = rows[iblte].theta;
                let my_te_dstr = rows[iblte].dstr;
                let my_te_ctau = rows[iblte].ctau;
                let tte = my_te_theta + other_te.theta;
                let dte = my_te_dstr + other_te.dstr + ante;
                let cte = if tte > 0.0 {
                    (my_te_ctau * my_te_theta + other_te.ctau * other_te.theta) / tte
                } else {
                    0.03
                };
                station_2 = build_and_fill_station(
                    xsi, uei, thi, (dsi - dswaki).max(1e-12), dswaki, cti, ami,
                    FlowType::Wake, msq, re,
                );
                sys[0] = [1.0, 0.0, 0.0, 0.0];
                rhs[0] = cte - station_2.ctau;
                sys[1] = [0.0, 1.0, 0.0, 0.0];
                rhs[1] = tte - station_2.theta;
                sys[2] = [0.0, 0.0, 1.0, 0.0];
                rhs[2] = dte - station_2.delta_star - station_2.dw;
            } else {
                let (res, jac) = if tran && tr_result.is_some() {
                    trdif_full(
                        &station_1,
                        &station_2,
                        tr_result.as_ref().unwrap(),
                        ncrit,
                        msq,
                        re,
                    )
                } else if simi {
                    bldif_full_simi(&station_2, FlowType::Laminar, msq, re)
                } else {
                    bldif(&station_1, &station_2, ft, msq, re)
                };
                for eq in 0..3 {
                    for var in 0..4 {
                        sys[eq][var] = jac.vs2[eq][var];
                    }
                }
                rhs[0] = res.res_third;
                rhs[1] = res.res_mom;
                rhs[2] = res.res_shape;
            }

            // ---- First iteration: set baselines (xbl.f:1101-1128) ----
            if itbl == 0 {
                ueref = station_2.u;
                hkref = station_2.hk;

                // Transition reversal (xbl.f:1108-1115)
                if is_laminar_interval(ibl, itran) && is_transition_or_turbulent(ibl, itrold) {
                    if ibl > 1 {
                        let prev = &rows[ibl - 1];
                        let hm = prev.dstr / prev.theta.max(1e-12);
                        hkref = hkin(hm, msq).hk;
                    }
                }

                // Reinitialize CTI for newly turbulent stations (xbl.f:1118-1126)
                if is_laminar_interval(ibl, itrold) {
                    if tran {
                        cti = 0.03;
                    }
                    if turb {
                        cti = rows.get(ibl.wrapping_sub(1)).map_or(0.03, |r| r.ctau);
                    }
                    if tran || turb {
                        station_2.ctau = cti;
                    }
                }
            }

            // ---- 4th equation for MRCHDU ----
            if simi || ibl == iblte + 1 {
                // Prescribe Ue (xbl.f:1131-1138)
                // U2_UEI = 1.0 for incompressible
                sys[3] = [0.0, 0.0, 0.0, 1.0];
                rhs[3] = ueref - station_2.u;
            } else {
                // Ue-Hk characteristic slope (xbl.f:1142-1176)
                let hk2_t2 = station_2.derivs.hk_h * station_2.derivs.h_theta;
                let hk2_d2 = station_2.derivs.hk_h * station_2.derivs.h_delta_star;
                let hk2_u2 = 0.0_f64; // incompressible

                // Solve auxiliary system for dUe/dHk (xbl.f:1144-1162)
                let mut vtmp = sys;
                let mut vztmp = rhs;
                vtmp[3] = [0.0, hk2_t2, hk2_d2, hk2_u2 + 1.0]; // U2_UEI=1 for incomp
                vztmp[3] = 1.0;
                gauss_4x4(&mut vtmp, &mut vztmp);

                sennew = SENSWT * vztmp[3] * hkref / ueref.abs().max(1e-6);

                // Smooth SENS (xbl.f:1163-1167)
                if itbl < 5 {
                    sens = sennew;
                } else if itbl < 15 {
                    sens = 0.5 * (sens + sennew);
                }

                // Set quasi-normal constraint (xbl.f:1170-1175)
                let ueref_safe = ueref.abs().max(1e-6);
                sys[3] = [
                    0.0,
                    hk2_t2 * hkref,
                    hk2_d2 * hkref,
                    (hk2_u2 * hkref + sens / ueref_safe) * 1.0, // U2_UEI=1
                ];
                rhs[3] = -(hkref * hkref) * (station_2.hk / hkref - 1.0)
                    - sens * (station_2.u / ueref_safe - 1.0);
            }

            // ---- Solve ----
            gauss_4x4(&mut sys, &mut rhs);

            // ---- DMAX and RLX (xbl.f:1184-1190) ----
            dmax = (rhs[1] / thi.max(1e-12))
                .abs()
                .max((rhs[2] / dsi.max(1e-12)).abs())
                .max((rhs[3] / uei.abs().max(1e-6)).abs());
            if is_transition_or_turbulent(ibl, itran) {
                dmax = dmax.max((rhs[0] / (10.0 * cti.max(1e-12))).abs());
            }
            let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

            if std::env::var("RUSTFOIL_WAKE_ITER_DEBUG").is_ok() && wake && ibl == iblte + 2 {
                eprintln!(
                    "[WAKE ITER] ibl={} it={} theta={:.8e} dsi_total={:.8e} ctau={:.8e} ue={:.8e} dmax={:.8e} rlx={:.8e} rhs=[{:.8e}, {:.8e}, {:.8e}, {:.8e}]",
                    ibl,
                    itbl + 1,
                    thi,
                    dsi,
                    cti,
                    uei,
                    dmax,
                    rlx,
                    rhs[0],
                    rhs[1],
                    rhs[2],
                    rhs[3]
                );
            }
            if std::env::var("RUSTFOIL_TAIL_ITER_DEBUG").is_ok()
                && (iblte.saturating_sub(2)..=iblte).contains(&ibl)
            {
                eprintln!(
                    "[TAIL ITER DU] side={} ibl={} it={} theta={:.12e} dsi_total={:.12e} ctau={:.12e} ue={:.12e} ampl={:.12e} dmax={:.12e} rlx={:.12e} rhs=[{:.12e}, {:.12e}, {:.12e}, {:.12e}]",
                    side,
                    ibl,
                    itbl + 1,
                    thi,
                    dsi,
                    cti,
                    uei,
                    ami,
                    dmax,
                    rlx,
                    rhs[0],
                    rhs[1],
                    rhs[2],
                    rhs[3]
                );
            }

            // ---- Updates (xbl.f:1193-1197) ----
            if is_laminar_interval(ibl, itran) {
                ami += rlx * rhs[0];
            }
            if is_transition_or_turbulent(ibl, itran) {
                cti += rlx * rhs[0];
            }
            thi += rlx * rhs[1];
            dsi += rlx * rhs[2];
            uei += rlx * rhs[3];

            // ---- Clamp CTI ----
            if is_transition_or_turbulent(ibl, itran) {
                cti = cti.clamp(1e-7, 0.30);
            }

            // ---- DSLIM ----
            let hklim = if ibl <= iblte { 1.02 } else { 1.00005 };
            let mut dsw = dsi - dswaki;
            dslim(&mut dsw, thi, msq, hklim);
            dsi = dsw + dswaki;

            if dmax <= DEPS {
                break;
            }
        }

        let mut fallback_station_final: Option<BlStation> = None;

        // ---- Convergence failure fallback ----
        if dmax > 0.1 && ibl > 2 {
            if ibl <= iblte {
                let prev_x = rows[ibl - 1].x.max(1e-12);
                let ratio = (xsi / prev_x).sqrt();
                thi = rows[ibl - 1].theta * ratio;
                dsi = rows[ibl - 1].dstr * ratio;
                uei = rows[ibl - 1].uedg;
            } else if ibl == iblte + 1 {
                let my_te = &rows[iblte];
                let tte = my_te.theta + other_te.theta;
                let dte = my_te.dstr + other_te.dstr + ante;
                let cte = if tte > 0.0 {
                    (my_te.ctau * my_te.theta + other_te.ctau * other_te.theta) / tte
                } else {
                    0.03
                };
                cti = cte;
                thi = tte;
                dsi = dte;
                uei = rows[ibl - 1].uedg;
            } else {
                thi = rows[ibl - 1].theta;
                let ratlen =
                    (xsi - rows[ibl - 1].x) / (10.0 * rows[ibl - 1].dstr.max(1e-12));
                dsi = (rows[ibl - 1].dstr + thi * ratlen) / (1.0 + ratlen);
                uei = rows[ibl - 1].uedg;
            }
            if is_transition_interval(ibl, itran) {
                cti = 0.05;
            }
            if xfoil_ibl(ibl) > itran {
                cti = rows[ibl - 1].ctau;
            }

            // Match XFOIL label 109 path after extrapolation: BLPRV/BLKIN,
            // rerun TRCHEK if still pre-transition, then rebuild the station
            // with the updated flow classification before storing.
            let mut fallback_station = build_and_fill_station(
                xsi,
                uei,
                thi,
                (dsi - dswaki).max(1e-12),
                dswaki,
                cti,
                ami,
                flow_type_for(ibl, itran, wake),
                msq,
                re,
            );
            if !simi && !turb {
                let tr = trchek2_full(
                    station_1.x,
                    fallback_station.x,
                    station_1.theta,
                    fallback_station.theta,
                    station_1.delta_star,
                    fallback_station.delta_star,
                    station_1.u,
                    fallback_station.u,
                    station_1.hk,
                    fallback_station.hk,
                    station_1.r_theta,
                    fallback_station.r_theta,
                    station_1.ampl,
                    ncrit,
                    x_forced,
                    msq,
                    re,
                );
                if rustfoil_bl::is_debug_active() {
                    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::trchek2_final(
                        iteration,
                        side,
                        ibl + 1,
                        tr.iterations,
                        tr.converged,
                        station_1.x,
                        fallback_station.x,
                        station_1.ampl,
                        tr.ampl2,
                        tr.ax,
                        tr.xt,
                        ncrit,
                        tr.transition,
                        tr.forced,
                        Some(tr.wf1),
                        Some(tr.wf2),
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                    ));
                }
                if std::env::var("RUSTFOIL_TRCHEK_DEBUG").is_ok() {
                    eprintln!(
                        "[TRCHEK DU FALLBACK] side={} ibl={} itran={} itrold={} x1={:.12e} x2={:.12e} th1={:.12e} th2={:.12e} ds1={:.12e} ds2={:.12e} u1={:.12e} u2={:.12e} hk1={:.12e} hk2={:.12e} rt1={:.12e} rt2={:.12e} ampl1={:.12e} ampl2={:.12e} transition={} xt={:.12e} xforced={:?}",
                        side,
                        ibl,
                        itran,
                        itrold,
                        station_1.x,
                        fallback_station.x,
                        station_1.theta,
                        fallback_station.theta,
                        station_1.delta_star,
                        fallback_station.delta_star,
                        station_1.u,
                        fallback_station.u,
                        station_1.hk,
                        fallback_station.hk,
                        station_1.r_theta,
                        fallback_station.r_theta,
                        station_1.ampl,
                        tr.ampl2,
                        tr.transition,
                        tr.xt,
                        x_forced,
                    );
                }
                ami = tr.ampl2;
                if tr.transition {
                    itran = xfoil_ibl(ibl);
                    tran = true;
                    xt = tr.xt;
                    tr_result = Some(tr);
                } else {
                    itran = xfoil_ibl(ibl) + 2;
                    tran = false;
                    tr_result = None;
                }
                fallback_station = build_and_fill_station(
                    xsi,
                    uei,
                    thi,
                    (dsi - dswaki).max(1e-12),
                    dswaki,
                    cti,
                    ami,
                    flow_type_for(ibl, itran, wake),
                    msq,
                    re,
                );
            }
            fallback_station_final = Some(fallback_station);
        }

        // ---- Store ----
        sens = sennew;
        let station_final = fallback_station_final.unwrap_or_else(|| {
            build_and_fill_station(
                xsi,
                uei,
                thi,
                (dsi - dswaki).max(1e-12),
                dswaki,
                cti,
                ami,
                flow_type_for(ibl, itran, wake),
                msq,
                re,
            )
        });
        if std::env::var("RUSTFOIL_WAKE_MARCH_DEBUG").is_ok()
            && ((wake && (iblte + 1..=iblte + 3).contains(&ibl)) || ibl == iblte)
        {
            eprintln!(
                "[WAKE MARCH UE] ibl={} x={:.8e} u={:.8e} theta={:.8e} dsi_total={:.8e} delta_star={:.8e} dw={:.8e} ctau={:.8e} ampl={:.8e} wake={}",
                ibl,
                xsi,
                uei,
                thi,
                dsi,
                station_final.delta_star,
                dswaki,
                cti,
                ami,
                wake
            );
        }
        store_to_row(&station_final, &mut rows[ibl], dsi, uei, cti, ami, ibl, itran, wake);
        if rustfoil_bl::is_debug_active() {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::mrchue(
                side,
                ibl + 1,
                station_final.x,
                station_final.u,
                station_final.theta,
                station_final.delta_star,
                station_final.hk,
                station_final.cf,
            ));
        }

        // Set station_1 for next station
        station_1 = station_final;

        // Transition / TE bookkeeping
        if tran {
            xssitr = xt;
        }
        if tran || ibl == iblte {
            turb = true;
        }
    }

    (itran, xssitr)
}
