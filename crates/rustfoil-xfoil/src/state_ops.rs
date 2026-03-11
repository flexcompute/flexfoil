use nalgebra::DMatrix;
use rustfoil_bl::state::BlStation;

use crate::state::{XfoilBlRow, XfoilState, XfoilSurface};

pub fn compute_arc_lengths(x: &[f64], y: &[f64]) -> Vec<f64> {
    let mut s = vec![0.0; x.len()];
    for i in 1..x.len() {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
    s
}

pub fn specal(state: &mut XfoilState, alpha_rad: f64) {
    state.alpha_rad = alpha_rad;
    let cosa = alpha_rad.cos();
    let sina = alpha_rad.sin();
    for i in 0..state.n_panel_nodes() {
        state.qinv[i] = cosa * state.qinvu_0[i] + sina * state.qinvu_90[i];
        state.qinv_a[i] = -sina * state.qinvu_0[i] + cosa * state.qinvu_90[i];
    }
}

pub fn stfind(state: &mut XfoilState) {
    let gamma = &state.gam;
    let n = gamma.len();
    let mut ist = n.saturating_div(2).saturating_sub(1);
    for i in 0..n.saturating_sub(1) {
        if gamma[i] >= 0.0 && gamma[i + 1] < 0.0 {
            ist = i;
            break;
        }
    }

    let dgam = gamma[ist + 1] - gamma[ist];
    let ds = state.panel_s[ist + 1] - state.panel_s[ist];
    let mut sst = if gamma[ist] < -gamma[ist + 1] {
        state.panel_s[ist] - ds * (gamma[ist] / dgam)
    } else {
        state.panel_s[ist + 1] - ds * (gamma[ist + 1] / dgam)
    };
    sst = sst.clamp(state.panel_s[ist] + 1.0e-7, state.panel_s[ist + 1] - 1.0e-7);

    state.ist = ist;
    state.sst = sst;
    state.sst_go = (sst - state.panel_s[ist + 1]) / dgam;
    state.sst_gp = (state.panel_s[ist] - sst) / dgam;
}

pub fn iblpan(state: &mut XfoilState) {
    state.ipan_upper.clear();
    state.ipan_lower.clear();
    state.ipan_upper.push(usize::MAX);
    state.ipan_lower.push(usize::MAX);

    for i in (0..=state.ist).rev() {
        state.ipan_upper.push(i);
    }
    state.iblte_upper = state.ipan_upper.len() - 1;
    state.nbl_upper = state.ipan_upper.len();

    for i in (state.ist + 1)..state.n_panel_nodes() {
        state.ipan_lower.push(i);
    }
    state.iblte_lower = state.ipan_lower.len() - 1;

    for iw in 0..state.wake_x.len() {
        state.ipan_lower.push(state.n_panel_nodes() + iw);
        state.ipan_upper.push(state.n_panel_nodes() + iw);
    }
    state.nbl_lower = state.ipan_lower.len();

    isys_from_ipan(state);
}

fn isys_from_ipan(state: &mut XfoilState) {
    state.isys_upper = vec![None; state.ipan_upper.len()];
    state.isys_lower = vec![None; state.ipan_lower.len()];

    let mut iv = 0usize;
    for ibl in 2..=state.nbl_upper {
        state.isys_upper[ibl - 1] = Some(iv);
        iv += 1;
    }
    for ibl in 2..=state.nbl_lower {
        state.isys_lower[ibl - 1] = Some(iv);
        iv += 1;
    }
}

pub fn xicalc(state: &mut XfoilState) {
    let xeps = 1.0e-7 * (state.panel_s.last().copied().unwrap_or(1.0) - state.panel_s[0]).abs();

    let mut upper_rows = Vec::with_capacity(state.ipan_upper.len());
    let mut lower_rows = Vec::with_capacity(state.ipan_lower.len());

    upper_rows.push(XfoilBlRow::default());
    lower_rows.push(XfoilBlRow::default());

    for &ipan in state.ipan_upper.iter().skip(1).take(state.iblte_upper) {
        let x = (state.sst - state.panel_s[ipan]).max(xeps);
        upper_rows.push(seed_row(x, state.panel_x[ipan], ipan, state.qinv[ipan].abs(), false));
    }

    for &ipan in state.ipan_lower.iter().skip(1).take(state.iblte_lower) {
        let x = (state.panel_s[ipan] - state.sst).max(xeps);
        lower_rows.push(seed_row(x, state.panel_x[ipan], ipan, (-state.qinv[ipan]).abs(), false));
    }

    let upper_te_x = upper_rows.last().map(|row| row.x).unwrap_or(0.0);
    let lower_te_x = lower_rows.last().map(|row| row.x).unwrap_or(0.0);

    for iw in 0..state.wake_x.len() {
        let arc = if iw == 0 {
            state.wake_s.get(iw).copied().unwrap_or(0.0)
        } else {
            state.wake_s[iw]
        };
        upper_rows.push(seed_row(
            upper_te_x + arc,
            state.wake_x[iw],
            state.n_panel_nodes() + iw,
            state.wake_qinv.get(iw).copied().unwrap_or(0.0).abs(),
            true,
        ));
        lower_rows.push(seed_row(
            lower_te_x + arc,
            state.wake_x[iw],
            state.n_panel_nodes() + iw,
            state.wake_qinv.get(iw).copied().unwrap_or(0.0).abs(),
            true,
        ));
    }

    set_stagnation_anchor(&mut upper_rows, state, XfoilSurface::Upper);
    set_stagnation_anchor(&mut lower_rows, state, XfoilSurface::Lower);
    state.upper_rows = upper_rows;
    state.lower_rows = lower_rows;
    update_wgap(state);
}

fn set_stagnation_anchor(rows: &mut [XfoilBlRow], state: &XfoilState, surface: XfoilSurface) {
    if rows.len() < 2 {
        return;
    }
    let next = &rows[1];
    let mut anchor = seed_row(0.0, next.x_coord, usize::MAX, 0.0, false);
    anchor.theta = 0.0;
    anchor.dstr = 0.0;
    anchor.mass = 0.0;
    anchor.ctau = 0.0;
    anchor.uedg = 0.0;
    anchor.uinv = 0.0;
    anchor.uinv_a = 0.0;
    anchor.panel_idx = usize::MAX;
    if surface == XfoilSurface::Lower {
        anchor.x_coord = state.panel_x[(state.ist + 1).min(state.panel_x.len() - 1)];
    }
    rows[0] = anchor;
}

fn seed_row(x: f64, x_coord: f64, panel_idx: usize, ue: f64, is_wake: bool) -> XfoilBlRow {
    let seed = BlStation::stagnation(x.max(1.0e-6), ue.max(0.01), 1.0e6);
    XfoilBlRow {
        x,
        x_coord,
        panel_idx,
        uedg: ue,
        uinv: ue,
        uinv_a: 0.0,
        theta: seed.theta,
        dstr: seed.delta_star,
        mass: ue * seed.delta_star,
        ctau: if is_wake { 0.03 } else { 0.0 },
        ampl: 0.0,
        hk: seed.hk,
        h: seed.h,
        hs: seed.hs,
        hc: seed.hc,
        r_theta: seed.r_theta,
        cf: if is_wake { 0.0 } else { 0.003 },
        cd: if is_wake { 0.0 } else { 0.001 },
        us: seed.us,
        cq: seed.cq,
        de: seed.de,
        dw: 0.0,
        is_laminar: !is_wake,
        is_turbulent: is_wake,
        is_wake,
        derivs: seed.derivs,
    }
}

pub fn uicalc(state: &mut XfoilState) {
    let n_panel_nodes = state.n_panel_nodes();
    let panel_qinv = &state.qinv;
    let wake_qinv = &state.wake_qinv;
    let panel_qinv_a = &state.qinv_a;
    let wake_qinv_a = &state.wake_qinv_a;
    let nbl_upper = state.nbl_upper;
    let nbl_lower = state.nbl_lower;

    if let Some(anchor) = state.upper_rows.first_mut() {
        anchor.uinv = 0.0;
        anchor.uinv_a = 0.0;
    }
    if let Some(anchor) = state.lower_rows.first_mut() {
        anchor.uinv = 0.0;
        anchor.uinv_a = 0.0;
    }

    for row in state.upper_rows.iter_mut().skip(1).take(nbl_upper.saturating_sub(1)) {
        row.uinv = signed_q_component_at(
            panel_qinv,
            wake_qinv,
            n_panel_nodes,
            row.panel_idx,
            XfoilSurface::Upper,
        );
        row.uinv_a = signed_q_component_at(
            panel_qinv_a,
            wake_qinv_a,
            n_panel_nodes,
            row.panel_idx,
            XfoilSurface::Upper,
        );
    }

    for row in state.lower_rows.iter_mut().skip(1).take(nbl_lower.saturating_sub(1)) {
        row.uinv = signed_q_component_at(
            panel_qinv,
            wake_qinv,
            n_panel_nodes,
            row.panel_idx,
            XfoilSurface::Lower,
        );
        row.uinv_a = signed_q_component_at(
            panel_qinv_a,
            wake_qinv_a,
            n_panel_nodes,
            row.panel_idx,
            XfoilSurface::Lower,
        );
    }
}

pub fn uedginit(state: &mut XfoilState) {
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        row.uedg = row.uinv;
    }

    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        row.uedg = row.uinv;
    }
}

fn signed_q_component_at(
    panel_values: &[f64],
    wake_values: &[f64],
    n_panel_nodes: usize,
    panel_idx: usize,
    surface: XfoilSurface,
) -> f64 {
    if panel_idx < panel_values.len() {
        return surface.vti() * panel_values[panel_idx];
    }
    let wake_idx = panel_idx.saturating_sub(n_panel_nodes);
    surface.vti() * wake_values.get(wake_idx).copied().unwrap_or(0.0)
}

pub fn qvfue(state: &mut XfoilState) {
    state.qvis = state
        .qinv
        .iter()
        .copied()
        .chain(state.wake_qinv.iter().copied())
        .collect();
    for row in state.upper_rows.iter().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        if row.panel_idx < state.qvis.len() {
            state.qvis[row.panel_idx] = row.uedg;
        }
    }
    for row in state.lower_rows.iter().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        if row.panel_idx < state.qvis.len() {
            state.qvis[row.panel_idx] = -row.uedg;
        }
    }
}

pub fn gamqv(state: &mut XfoilState) {
    for i in 0..state.n_panel_nodes() {
        state.gam[i] = state.qvis.get(i).copied().unwrap_or(0.0);
        state.gam_a[i] = state.qinv_a.get(i).copied().unwrap_or(0.0);
    }
}

pub fn stmove(state: &mut XfoilState) {
    let old_ist = state.ist;
    let old_upper = state.upper_rows.clone();
    let old_lower = state.lower_rows.clone();
    stfind(state);
    if state.ist != old_ist {
        iblpan(state);
        uicalc(state);
        xicalc(state);
        shift_stagnation_state(state, old_ist, &old_upper, &old_lower);
        apply_ue_floor(state);
    } else {
        xicalc(state);
        restore_rows_identity(state, &old_upper, &old_lower);
    }
    recompute_mass_from_state(state);
}

pub fn dsset(state: &mut XfoilState) {
    // XFOIL: DO IBL=2, NBL(IS). Upper has no wake, lower includes wake.
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        if row.uedg.abs() > 1.0e-12 {
            row.dstr = row.mass / row.uedg;
        }
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        if row.uedg.abs() > 1.0e-12 {
            row.dstr = row.mass / row.uedg;
        }
    }
}

pub fn ueset(state: &mut XfoilState) {
    // Build coupling vector respecting NBL bounds (XFOIL: JS=1 to NBL(JS)).
    // Upper side: NBL(1) = IBLTE(1), no wake.
    // Lower side: NBL(2) = IBLTE(2) + NW, includes wake.
    let mut coupling = Vec::new();
    for row in state.upper_rows.iter().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        coupling.push((XfoilSurface::Upper.vti(), row.panel_idx, row.mass));
    }
    for row in state.lower_rows.iter().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        coupling.push((XfoilSurface::Lower.vti(), row.panel_idx, row.mass));
    }

    let dij = state.dij.clone();

    // Outer loop also respects NBL bounds (XFOIL: IS=1 to NBL(IS)).
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        row.uedg = coupled_uedg(row.uinv, row.panel_idx, XfoilSurface::Upper.vti(), &coupling, &dij);
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        row.uedg = coupled_uedg(row.uinv, row.panel_idx, XfoilSurface::Lower.vti(), &coupling, &dij);
    }
}

fn coupled_uedg(
    uinv: f64,
    panel_idx: usize,
    vti: f64,
    coupling: &[(f64, usize, f64)],
    dij: &DMatrix<f64>,
) -> f64 {
    let mut ue = uinv;
    for (other_vti, other_panel, mass) in coupling {
        if panel_idx < dij.nrows() && *other_panel < dij.ncols() {
            ue -= vti * *other_vti * dij[(panel_idx, *other_panel)] * *mass;
        }
    }
    ue
}

fn update_wgap(state: &mut XfoilState) {
    state.wgap = vec![0.0; state.wake_x.len()];
    if state.sharp || state.wake_x.is_empty() || state.ante.abs() <= 1.0e-12 {
        return;
    }

    let telrat = 2.50;
    let xp0 = state.panel_xp.first().copied().unwrap_or(0.0);
    let yp0 = state.panel_yp.first().copied().unwrap_or(0.0);
    let xpn = state.panel_xp.last().copied().unwrap_or(0.0);
    let ypn = state.panel_yp.last().copied().unwrap_or(0.0);
    let denom = ((xp0 * xp0 + yp0 * yp0) * (xpn * xpn + ypn * ypn)).sqrt();
    if denom <= 1.0e-12 {
        return;
    }
    let crosp = (xp0 * ypn - yp0 * xpn) / denom;
    let mut dwdxte = crosp / (1.0 - crosp * crosp).max(1.0e-12).sqrt();
    dwdxte = dwdxte.clamp(-3.0 / telrat, 3.0 / telrat);
    let aa = 3.0 + telrat * dwdxte;
    let bb = -2.0 - telrat * dwdxte;

    let wake_start_x = state
        .lower_rows
        .get(state.iblte_lower)
        .map(|row| row.x)
        .unwrap_or(0.0);
    for iw in 0..state.wake_x.len() {
        let ibl = state.iblte_lower + 1 + iw;
        let xssi = state.lower_rows.get(ibl).map(|row| row.x).unwrap_or(wake_start_x);
        let zn = 1.0 - (xssi - wake_start_x) / (telrat * state.ante);
        state.wgap[iw] = if zn >= 0.0 {
            state.ante * (aa + bb * zn) * zn * zn
        } else {
            0.0
        };
    }
}

fn restore_rows_identity(state: &mut XfoilState, old_upper: &[XfoilBlRow], old_lower: &[XfoilBlRow]) {
    let upper_len = state.upper_rows.len().min(old_upper.len());
    let lower_len = state.lower_rows.len().min(old_lower.len());
    copy_row_window(&mut state.upper_rows, old_upper, upper_len);
    copy_row_window(&mut state.lower_rows, old_lower, lower_len);
}

fn shift_stagnation_state(state: &mut XfoilState, old_ist: usize, old_upper: &[XfoilBlRow], old_lower: &[XfoilBlRow]) {
    if state.ist > old_ist {
        let idif = state.ist - old_ist;
        for ibl in ((idif + 2)..=state.nbl_upper).rev() {
            copy_row_values(&mut state.upper_rows[ibl - 1], &old_upper[ibl - idif - 1]);
        }

        let seed_idx = idif + 1;
        let seed = state.upper_rows[seed_idx].clone();
        let dudx = if seed.x.abs() > 1.0e-12 { seed.uedg / seed.x } else { 0.0 };
        for ibl in (2..=idif + 1).rev() {
            copy_row_values(&mut state.upper_rows[ibl - 1], &seed);
            state.upper_rows[ibl - 1].uedg = dudx * state.upper_rows[ibl - 1].x;
        }

        for ibl in 2..=state.nbl_lower {
            copy_row_values(&mut state.lower_rows[ibl - 1], &old_lower[ibl + idif - 1]);
        }
    } else {
        let idif = old_ist - state.ist;
        for ibl in ((idif + 2)..=state.nbl_lower).rev() {
            copy_row_values(&mut state.lower_rows[ibl - 1], &old_lower[ibl - idif - 1]);
        }

        let seed_idx = idif + 1;
        let seed = state.lower_rows[seed_idx].clone();
        let dudx = if seed.x.abs() > 1.0e-12 { seed.uedg / seed.x } else { 0.0 };
        for ibl in (2..=idif + 1).rev() {
            copy_row_values(&mut state.lower_rows[ibl - 1], &seed);
            state.lower_rows[ibl - 1].uedg = dudx * state.lower_rows[ibl - 1].x;
        }

        for ibl in 2..=state.nbl_upper {
            copy_row_values(&mut state.upper_rows[ibl - 1], &old_upper[ibl + idif - 1]);
        }
    }
}

fn copy_row_window(dst: &mut [XfoilBlRow], src: &[XfoilBlRow], len: usize) {
    for idx in 1..len {
        copy_row_values(&mut dst[idx], &src[idx]);
    }
}

fn copy_row_values(dst: &mut XfoilBlRow, src: &XfoilBlRow) {
    dst.ctau = src.ctau;
    dst.theta = src.theta;
    dst.dstr = src.dstr;
    dst.uedg = src.uedg;
    dst.uinv = src.uinv;
    dst.uinv_a = src.uinv_a;
}

fn apply_ue_floor(state: &mut XfoilState) {
    const UEPS: f64 = 1.0e-7;
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        if row.uedg <= UEPS {
            row.uedg = UEPS;
            if row.panel_idx < state.qvis.len() {
                state.qvis[row.panel_idx] = XfoilSurface::Upper.vti() * UEPS;
            }
            if row.panel_idx < state.gam.len() {
                state.gam[row.panel_idx] = XfoilSurface::Upper.vti() * UEPS;
            }
        }
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        if row.uedg <= UEPS {
            row.uedg = UEPS;
            if row.panel_idx < state.qvis.len() {
                state.qvis[row.panel_idx] = XfoilSurface::Lower.vti() * UEPS;
            }
            if row.panel_idx < state.gam.len() {
                state.gam[row.panel_idx] = XfoilSurface::Lower.vti() * UEPS;
            }
        }
    }
}

fn recompute_mass_from_state(state: &mut XfoilState) {
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        row.mass = row.dstr * row.uedg;
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        row.mass = row.dstr * row.uedg;
    }
}
