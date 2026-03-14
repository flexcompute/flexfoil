use nalgebra::DMatrix;
use rustfoil_bl::state::BlStation;

use crate::canonical_state::{
    recompute_mass_from_state, rebuild_gam_from_qvis, rebuild_qvis_from_rows,
};
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
        state.gam[i] = cosa * state.gamu_0[i] + sina * state.gamu_90[i];
        state.gam_a[i] = -sina * state.gamu_0[i] + cosa * state.gamu_90[i];
    }
    state.wake_qinv = state
        .wake_qinvu_0
        .iter()
        .zip(state.wake_qinvu_90.iter())
        .map(|(&q0, &q90)| cosa * q0 + sina * q90)
        .collect();
    state.wake_qinv_a = state
        .wake_qinvu_0
        .iter()
        .zip(state.wake_qinvu_90.iter())
        .map(|(&q0, &q90)| -sina * q0 + cosa * q90)
        .collect();
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
    let s0 = state.panel_s[ist];
    let s1 = state.panel_s[ist + 1];
    let ds = s1 - s0;

    state.ist = ist;

    let interp_xy = |sst: f64| {
        let ds = s1 - s0;
        if ds.abs() <= 1.0e-20 {
            (state.panel_x[ist], state.panel_y[ist])
        } else {
            let t = (sst - s0) / ds;
            (
                state.panel_x[ist] + t * (state.panel_x[ist + 1] - state.panel_x[ist]),
                state.panel_y[ist] + t * (state.panel_y[ist + 1] - state.panel_y[ist]),
            )
        }
    };

    if dgam.abs() <= 1.0e-20 || !dgam.is_finite() {
        // Degenerate stagnation interpolation: keep the point finite and
        // disable XI_ULE sensitivity rather than propagating NaNs.
        state.sst = 0.5 * (s0 + s1);
        state.sst_go = 0.0;
        state.sst_gp = 0.0;
        if rustfoil_bl::is_debug_active() {
            let (xst, yst) = interp_xy(state.sst);
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::stfind(ist + 1, state.sst, xst, yst));
        }
        return;
    }

    let mut sst = if gamma[ist] < -gamma[ist + 1] {
        s0 - ds * (gamma[ist] / dgam)
    } else {
        s1 - ds * (gamma[ist + 1] / dgam)
    };
    sst = sst.clamp(s0 + 1.0e-7, s1 - 1.0e-7);

    state.sst = sst;
    state.sst_go = (sst - s1) / dgam;
    state.sst_gp = (s0 - sst) / dgam;

    if rustfoil_bl::is_debug_active() {
        let (xst, yst) = interp_xy(sst);
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::stfind(ist + 1, sst, xst, yst));
    }
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

    // Keep upper-surface wake pointers only as plotting ghosts. The canonical
    // BL/Newton owner is the lower side, matching XFOIL's NBL/ISYS contract.
    for iw in 0..state.wake_x.len() {
        state.ipan_lower.push(state.n_panel_nodes() + iw);
        state.ipan_upper.push(state.n_panel_nodes() + iw);
    }
    state.nbl_lower = state.ipan_lower.len();

    isys_from_ipan(state);

    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::iblpan(
            state.ist + 1,
            state.nbl_upper,
            state.nbl_lower,
            state.iblte_upper + 1,
            state.iblte_lower + 1,
        ));
    }
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

    upper_rows.push(state.upper_rows.first().cloned().unwrap_or_default());
    lower_rows.push(state.lower_rows.first().cloned().unwrap_or_default());

    for (ibl, &ipan) in state.ipan_upper.iter().enumerate().skip(1).take(state.iblte_upper) {
        let x = (state.sst - state.panel_s[ipan]).max(xeps);
        upper_rows.push(refresh_or_seed_row(
            state.upper_rows.get(ibl),
            x,
            state.panel_x[ipan],
            state.panel_y[ipan],
            ipan,
            signed_q_component_at(
                &state.qinv,
                &state.wake_qinv,
                state.n_panel_nodes(),
                ipan,
                XfoilSurface::Upper,
            ),
            signed_q_component_at(
                &state.qinv_a,
                &state.wake_qinv_a,
                state.n_panel_nodes(),
                ipan,
                XfoilSurface::Upper,
            ),
            false,
        ));
    }

    for (ibl, &ipan) in state.ipan_lower.iter().enumerate().skip(1).take(state.iblte_lower) {
        let x = (state.panel_s[ipan] - state.sst).max(xeps);
        lower_rows.push(refresh_or_seed_row(
            state.lower_rows.get(ibl),
            x,
            state.panel_x[ipan],
            state.panel_y[ipan],
            ipan,
            signed_q_component_at(
                &state.qinv,
                &state.wake_qinv,
                state.n_panel_nodes(),
                ipan,
                XfoilSurface::Lower,
            ),
            signed_q_component_at(
                &state.qinv_a,
                &state.wake_qinv_a,
                state.n_panel_nodes(),
                ipan,
                XfoilSurface::Lower,
            ),
            false,
        ));
    }

    let upper_te_x = upper_rows.last().map(|row| row.x).unwrap_or(0.0);
    let lower_te_x = lower_rows.last().map(|row| row.x).unwrap_or(0.0);
    let mut upper_wake_arc = upper_te_x;
    let mut lower_wake_arc = lower_te_x;

    for iw in 0..state.wake_x.len() {
        // XFOIL's XICALC does not reuse the target wake spacing array here.
        // It pins the first wake station to the TE xi and then accumulates the
        // actual geometric distance between wake nodes.
        if iw > 0 {
            let dx = state.wake_x[iw] - state.wake_x[iw - 1];
            let dy = state.wake_y[iw] - state.wake_y[iw - 1];
            let dxssi = (dx * dx + dy * dy).sqrt();
            upper_wake_arc += dxssi;
            lower_wake_arc += dxssi;
        }
        let wake_panel = state.n_panel_nodes() + iw;
        let upper_ibl = state.iblte_upper + 1 + iw;
        let lower_ibl = state.iblte_lower + 1 + iw;
        upper_rows.push(refresh_or_seed_row(
            state.upper_rows.get(upper_ibl),
            upper_wake_arc,
            state.wake_x[iw],
            state.wake_y[iw],
            wake_panel,
            signed_q_component_at(
                &state.qinv,
                &state.wake_qinv,
                state.n_panel_nodes(),
                wake_panel,
                XfoilSurface::Upper,
            ),
            signed_q_component_at(
                &state.qinv_a,
                &state.wake_qinv_a,
                state.n_panel_nodes(),
                wake_panel,
                XfoilSurface::Upper,
            ),
            true,
        ));
        lower_rows.push(refresh_or_seed_row(
            state.lower_rows.get(lower_ibl),
            lower_wake_arc,
            state.wake_x[iw],
            state.wake_y[iw],
            wake_panel,
            signed_q_component_at(
                &state.qinv,
                &state.wake_qinv,
                state.n_panel_nodes(),
                wake_panel,
                XfoilSurface::Lower,
            ),
            signed_q_component_at(
                &state.qinv_a,
                &state.wake_qinv_a,
                state.n_panel_nodes(),
                wake_panel,
                XfoilSurface::Lower,
            ),
            true,
        ));
    }

    set_stagnation_anchor(&mut upper_rows, state, XfoilSurface::Upper);
    set_stagnation_anchor(&mut lower_rows, state, XfoilSurface::Lower);
    state.upper_rows = upper_rows;
    state.lower_rows = lower_rows;
    update_wgap(state);
}

fn refresh_or_seed_row(
    existing: Option<&XfoilBlRow>,
    x: f64,
    x_coord: f64,
    y_coord: f64,
    panel_idx: usize,
    uinv: f64,
    uinv_a: f64,
    is_wake: bool,
) -> XfoilBlRow {
    let mut row = existing
        .cloned()
        .unwrap_or_else(|| seed_row(x, x_coord, y_coord, panel_idx, uinv, is_wake));
    row.x = x;
    row.x_coord = x_coord;
    row.y_coord = y_coord;
    row.panel_idx = panel_idx;
    row.uinv = uinv;
    row.uinv_a = uinv_a;
    row.is_wake = is_wake;
    if is_wake {
        row.is_laminar = false;
        row.is_turbulent = true;
    }
    row
}

fn set_stagnation_anchor(rows: &mut [XfoilBlRow], state: &XfoilState, surface: XfoilSurface) {
    if rows.len() < 2 {
        return;
    }
    let next = &rows[1];
    let mut anchor = seed_row(0.0, next.x_coord, next.y_coord, usize::MAX, 0.0, false);
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
        anchor.y_coord = state.panel_y[(state.ist + 1).min(state.panel_y.len() - 1)];
    }
    rows[0] = anchor;
}

fn seed_row(x: f64, x_coord: f64, y_coord: f64, panel_idx: usize, ue: f64, is_wake: bool) -> XfoilBlRow {
    let seed = BlStation::stagnation(x.max(1.0e-6), ue.max(0.01), 1.0e6);
    XfoilBlRow {
        x,
        x_coord,
        y_coord,
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
    let nbl_upper = state.nbl_upper.min(state.upper_rows.len()).min(state.ipan_upper.len());
    let nbl_lower = state.nbl_lower.min(state.lower_rows.len()).min(state.ipan_lower.len());

    if let Some(anchor) = state.upper_rows.first_mut() {
        anchor.uinv = 0.0;
        anchor.uinv_a = 0.0;
    }
    if let Some(anchor) = state.lower_rows.first_mut() {
        anchor.uinv = 0.0;
        anchor.uinv_a = 0.0;
    }

    for ibl in 1..nbl_upper {
        let panel_idx = state.ipan_upper[ibl];
        let row = &mut state.upper_rows[ibl];
        row.uinv = signed_q_component_at(
            panel_qinv,
            wake_qinv,
            n_panel_nodes,
            panel_idx,
            XfoilSurface::Upper,
        );
        row.uinv_a = signed_q_component_at(
            panel_qinv_a,
            wake_qinv_a,
            n_panel_nodes,
            panel_idx,
            XfoilSurface::Upper,
        );
    }

    for ibl in 1..nbl_lower {
        let panel_idx = state.ipan_lower[ibl];
        let row = &mut state.lower_rows[ibl];
        row.uinv = signed_q_component_at(
            panel_qinv,
            wake_qinv,
            n_panel_nodes,
            panel_idx,
            XfoilSurface::Lower,
        );
        row.uinv_a = signed_q_component_at(
            panel_qinv_a,
            wake_qinv_a,
            n_panel_nodes,
            panel_idx,
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
    rebuild_qvis_from_rows(state);
}

pub fn gamqv(state: &mut XfoilState) {
    rebuild_gam_from_qvis(state);
}

pub fn stmove(state: &mut XfoilState) {
    let old_ist = state.ist;
    maybe_debug_stmove_rows(state, "before_stfind", old_ist, old_ist);
    stfind(state);
    maybe_debug_stmove_rows(state, "after_stfind", old_ist, state.ist);
    if state.ist != old_ist {
        let old_upper = state.upper_rows.clone();
        let old_lower = state.lower_rows.clone();
        adjust_transition_indices(state, old_ist);
        iblpan(state);
        resize_row_storage_for_topology(state);
        uicalc(state);
        xicalc_geometry_only(state);
        maybe_debug_stmove_rows(state, "after_xicalc", old_ist, state.ist);
        shift_stagnation_state(state, old_ist, &old_upper, &old_lower);
        maybe_debug_stmove_rows(state, "after_shift", old_ist, state.ist);
        apply_ue_floor(state);
        refresh_shifted_row_views(state);
        recompute_mass_from_state(state);
        maybe_debug_stmove_rows(state, "after_mass", old_ist, state.ist);
    } else {
        xicalc_geometry_only(state);
        maybe_debug_stmove_rows(state, "after_xicalc", old_ist, state.ist);
    }
}

fn resize_row_storage_for_topology(state: &mut XfoilState) {
    let upper_len = state.ipan_upper.len();
    let lower_len = state.ipan_lower.len();
    if state.upper_rows.len() < upper_len {
        state.upper_rows.resize(upper_len, XfoilBlRow::default());
    } else {
        state.upper_rows.truncate(upper_len);
    }
    if state.lower_rows.len() < lower_len {
        state.lower_rows.resize(lower_len, XfoilBlRow::default());
    } else {
        state.lower_rows.truncate(lower_len);
    }
}

fn xicalc_geometry_only(state: &mut XfoilState) {
    let xeps = 1.0e-7 * (state.panel_s.last().copied().unwrap_or(1.0) - state.panel_s[0]).abs();

    if let Some(anchor) = state.upper_rows.first_mut() {
        anchor.x = 0.0;
        anchor.panel_idx = usize::MAX;
    }
    if let Some(anchor) = state.lower_rows.first_mut() {
        anchor.x = 0.0;
        anchor.panel_idx = usize::MAX;
        let idx = (state.ist + 1).min(state.panel_x.len().saturating_sub(1));
        anchor.x_coord = state.panel_x[idx];
        anchor.y_coord = state.panel_y[idx];
    }

    for (ibl, &ipan) in state.ipan_upper.iter().enumerate().skip(1).take(state.iblte_upper) {
        if let Some(row) = state.upper_rows.get_mut(ibl) {
            row.x = (state.sst - state.panel_s[ipan]).max(xeps);
            row.x_coord = state.panel_x[ipan];
            row.y_coord = state.panel_y[ipan];
            row.panel_idx = ipan;
        }
    }

    for (ibl, &ipan) in state.ipan_lower.iter().enumerate().skip(1).take(state.iblte_lower) {
        if let Some(row) = state.lower_rows.get_mut(ibl) {
            row.x = (state.panel_s[ipan] - state.sst).max(xeps);
            row.x_coord = state.panel_x[ipan];
            row.y_coord = state.panel_y[ipan];
            row.panel_idx = ipan;
        }
    }

    let upper_te_x = state
        .upper_rows
        .get(state.iblte_upper)
        .map(|row| row.x)
        .unwrap_or(0.0);
    let lower_te_x = state
        .lower_rows
        .get(state.iblte_lower)
        .map(|row| row.x)
        .unwrap_or(0.0);
    let mut upper_wake_arc = upper_te_x;
    let mut lower_wake_arc = lower_te_x;

    for iw in 0..state.wake_x.len() {
        if iw > 0 {
            let dx = state.wake_x[iw] - state.wake_x[iw - 1];
            let dy = state.wake_y[iw] - state.wake_y[iw - 1];
            let dxssi = (dx * dx + dy * dy).sqrt();
            upper_wake_arc += dxssi;
            lower_wake_arc += dxssi;
        }

        let wake_panel = state.n_panel_nodes() + iw;
        let upper_ibl = state.iblte_upper + 1 + iw;
        let lower_ibl = state.iblte_lower + 1 + iw;

        if let Some(row) = state.upper_rows.get_mut(upper_ibl) {
            row.x = upper_wake_arc;
            row.x_coord = state.wake_x[iw];
            row.y_coord = state.wake_y[iw];
            row.panel_idx = wake_panel;
        }
        if let Some(row) = state.lower_rows.get_mut(lower_ibl) {
            row.x = lower_wake_arc;
            row.x_coord = state.wake_x[iw];
            row.y_coord = state.wake_y[iw];
            row.panel_idx = wake_panel;
        }
    }

    update_wgap(state);
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
    // XFOIL's UESET consumes MASS as the coupling source term; keep the
    // cached row mass coherent with the current row-space BL state first.
    recompute_mass_from_state(state);

    let upper_before: Vec<f64> = state
        .upper_rows
        .iter()
        .skip(1)
        .take(state.nbl_upper.saturating_sub(1))
        .map(|row| row.uedg)
        .collect();
    let lower_before: Vec<f64> = state
        .lower_rows
        .iter()
        .skip(1)
        .take(state.nbl_lower.saturating_sub(1))
        .map(|row| row.uedg)
        .collect();
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
    let n_panel_nodes = state.n_panel_nodes();

    // Outer loop also respects NBL bounds (XFOIL: IS=1 to NBL(IS)).
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        row.uedg = coupled_uedg(
            row.uinv,
            row.panel_idx,
            XfoilSurface::Upper.vti(),
            n_panel_nodes,
            &coupling,
            &dij,
        );
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        row.uedg = coupled_uedg(
            row.uinv,
            row.panel_idx,
            XfoilSurface::Lower.vti(),
            n_panel_nodes,
            &coupling,
            &dij,
        );
    }

    if rustfoil_bl::is_debug_active() {
        let upper_surface = rustfoil_bl::UesetSurfaceData {
            ue_inviscid: state
                .upper_rows
                .iter()
                .skip(1)
                .take(state.nbl_upper.saturating_sub(1))
                .map(|row| row.uinv)
                .collect(),
            mass_defect: state
                .upper_rows
                .iter()
                .skip(1)
                .take(state.nbl_upper.saturating_sub(1))
                .map(|row| row.mass)
                .collect(),
            ue_before: upper_before,
            ue_after: state
                .upper_rows
                .iter()
                .skip(1)
                .take(state.nbl_upper.saturating_sub(1))
                .map(|row| row.uedg)
                .collect(),
        };
        let lower_surface = rustfoil_bl::UesetSurfaceData {
            ue_inviscid: state
                .lower_rows
                .iter()
                .skip(1)
                .take(state.nbl_lower.saturating_sub(1))
                .map(|row| row.uinv)
                .collect(),
            mass_defect: state
                .lower_rows
                .iter()
                .skip(1)
                .take(state.nbl_lower.saturating_sub(1))
                .map(|row| row.mass)
                .collect(),
            ue_before: lower_before,
            ue_after: state
                .lower_rows
                .iter()
                .skip(1)
                .take(state.nbl_lower.saturating_sub(1))
                .map(|row| row.uedg)
                .collect(),
        };
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::ueset(
            state.iterations + 1,
            state.nbl_upper,
            state.nbl_lower,
            upper_surface,
            lower_surface,
        ));
    }
}

fn coupled_uedg(
    uinv: f64,
    panel_idx: usize,
    vti: f64,
    n_panel_nodes: usize,
    coupling: &[(f64, usize, f64)],
    dij: &DMatrix<f64>,
) -> f64 {
    let mut ue = uinv;
    let debug_panels = matches!(
        panel_idx,
        0 | 69 | 70 | 98 | 102 | 103 | 104 | 146 | 147 | 148
    );
    let debug_first = std::env::var("RUSTFOIL_UESET_STATE_DEBUG").is_ok()
        && panel_idx != usize::MAX
        && (vti > 0.0 || (vti < 0.0 && debug_panels));
    let mut debug_upper = 0.0;
    let mut debug_lower_airfoil = 0.0;
    let mut debug_lower_wake = 0.0;
    for (other_vti, other_panel, mass) in coupling {
        if panel_idx < dij.nrows() && *other_panel < dij.ncols() {
            let contrib = -vti * *other_vti * dij[(panel_idx, *other_panel)] * *mass;
            ue += contrib;
            if debug_first {
                if *other_vti > 0.0 {
                    debug_upper += contrib;
                } else if *other_panel < n_panel_nodes {
                    debug_lower_airfoil += contrib;
                } else {
                    debug_lower_wake += contrib;
                }
            }
        }
    }
    if debug_first && debug_panels {
        eprintln!(
            "[STATE UESET] panel={} vti={:.0} uinv={:.8e} upper={:.8e} lower_airfoil={:.8e} lower_wake={:.8e} uedg={:.8e}",
            panel_idx, vti, uinv, debug_upper, debug_lower_airfoil, debug_lower_wake, ue
        );
        if std::env::var("RUSTFOIL_UESET_DETAIL_DEBUG").is_ok() {
            for (other_vti, other_panel, mass) in coupling {
                if *other_vti > 0.0
                    && panel_idx < dij.nrows()
                    && *other_panel < dij.ncols()
                    && debug_panels
                {
                    let dij_val = dij[(panel_idx, *other_panel)];
                    let contrib = -vti * *other_vti * dij_val * *mass;
                    if contrib.abs() > 1.0e-5 {
                        eprintln!(
                            "[STATE UESET DETAIL UPPER] panel={} other_panel={} dij={:.12e} mass={:.12e} contrib={:.12e}",
                            panel_idx, other_panel, dij_val, mass, contrib
                        );
                    }
                }
                if *other_vti < 0.0 && *other_panel < n_panel_nodes
                    && panel_idx < dij.nrows() && *other_panel < dij.ncols()
                {
                    let dij_val = dij[(panel_idx, *other_panel)];
                    let contrib = -vti * *other_vti * dij_val * *mass;
                    if contrib.abs() > 1.0e-5 && debug_panels {
                        eprintln!(
                            "[STATE UESET DETAIL LOWER] panel={} other_panel={} dij={:.8e} mass={:.8e} contrib={:.8e}",
                            panel_idx, other_panel, dij_val, mass, contrib
                        );
                    }
                }
                if *other_vti < 0.0 && *other_panel >= n_panel_nodes
                    && panel_idx < dij.nrows() && *other_panel < dij.ncols()
                    && debug_panels
                {
                    let dij_val = dij[(panel_idx, *other_panel)];
                    let contrib = -vti * *other_vti * dij_val * *mass;
                    if contrib.abs() > 1.0e-5 {
                        eprintln!(
                            "[STATE UESET DETAIL WAKE] panel={} other_panel={} dij={:.8e} mass={:.8e} contrib={:.8e}",
                            panel_idx, other_panel, dij_val, mass, contrib
                        );
                    }
                }
            }
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
        let wgap = if zn >= 0.0 {
            state.ante * (aa + bb * zn) * zn * zn
        } else {
            0.0
        };
        state.wgap[iw] = wgap;

        if let Some(row) = state.lower_rows.get_mut(ibl) {
            let wake_free_dstr = if row.is_wake {
                (row.dstr - row.dw).max(1.0e-12)
            } else {
                row.dstr.max(1.0e-12)
            };
            row.dw = wgap;
            row.dstr = wake_free_dstr + row.dw;
        }
        if let Some(row) = state.upper_rows.get_mut(state.iblte_upper + 1 + iw) {
            let wake_free_dstr = if row.is_wake {
                (row.dstr - row.dw).max(1.0e-12)
            } else {
                row.dstr.max(1.0e-12)
            };
            row.dw = wgap;
            row.dstr = wake_free_dstr + row.dw;
        }
    }
}

fn shift_stagnation_state(state: &mut XfoilState, old_ist: usize, old_upper: &[XfoilBlRow], old_lower: &[XfoilBlRow]) {
    if state.ist > old_ist {
        let idif = state.ist - old_ist;
        for ibl in ((idif + 2)..=state.nbl_upper).rev() {
            copy_row_values(&mut state.upper_rows[ibl - 1], &old_upper[ibl - idif - 1]);
        }

        let seed_idx = idif + 2;
        let seed = state.upper_rows[seed_idx - 1].clone();
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

        let seed_idx = idif + 2;
        let seed = state.lower_rows[seed_idx - 1].clone();
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

fn copy_row_values(dst: &mut XfoilBlRow, src: &XfoilBlRow) {
    dst.ctau = src.ctau;
    dst.theta = src.theta;
    dst.dstr = src.dstr;
    dst.uedg = src.uedg;
}

fn adjust_transition_indices(state: &mut XfoilState, old_ist: usize) {
    if state.ist > old_ist {
        let idif = state.ist - old_ist;
        state.itran_upper = state.itran_upper.saturating_add(idif);
        state.itran_lower = state.itran_lower.saturating_sub(idif);
    } else {
        let idif = old_ist - state.ist;
        state.itran_upper = state.itran_upper.saturating_sub(idif);
        state.itran_lower = state.itran_lower.saturating_add(idif);
    }
}

fn refresh_shifted_row_views(state: &mut XfoilState) {
    refresh_surface_row_views(&mut state.upper_rows, state.nbl_upper, state.iblte_upper, state.itran_upper);
    refresh_surface_row_views(&mut state.lower_rows, state.nbl_lower, state.iblte_lower, state.itran_lower);
}

fn refresh_surface_row_views(rows: &mut [XfoilBlRow], nbl: usize, iblte: usize, itran: usize) {
    let len = nbl.min(rows.len());
    for (ibl, row) in rows.iter_mut().enumerate().take(len).skip(1) {
        let wake = ibl > iblte;
        let turbulent = wake || ibl + 1 >= itran;
        if !wake && row.dw.abs() > 0.0 {
            // Wake rows store total DSTR = delta_star + dw. If a stagnation move
            // reclassifies that row back onto the airfoil, remove the carried wake gap.
            row.dstr = (row.dstr - row.dw).max(1.0e-12);
            row.dw = 0.0;
        }
        row.is_wake = wake;
        row.is_turbulent = turbulent;
        row.is_laminar = !turbulent;
        if row.is_laminar {
            row.ampl = row.ctau.max(0.0);
        }
    }
}

fn maybe_debug_stmove_rows(state: &XfoilState, phase: &str, old_ist: usize, new_ist: usize) {
    if std::env::var("RUSTFOIL_STMOVE_DEBUG").is_err() {
        return;
    }
    eprintln!(
        "[STMOVE DEBUG] phase={} old_ist={} new_ist={} nbl_upper={} iblte_upper={} nbl_lower={} iblte_lower={}",
        phase,
        old_ist + 1,
        new_ist + 1,
        state.nbl_upper,
        state.iblte_upper,
        state.nbl_lower,
        state.iblte_lower
    );
    debug_stmove_surface_rows("upper", &state.upper_rows, state.nbl_upper, phase);
    debug_stmove_surface_rows("lower", &state.lower_rows, state.nbl_lower, phase);
}

fn debug_stmove_surface_rows(label: &str, rows: &[XfoilBlRow], active_len: usize, phase: &str) {
    let mut ibls: Vec<usize> = (1..6.min(rows.len())).collect();
    let tail_limit = active_len.min(rows.len());
    let tail_start = tail_limit.saturating_sub(6).max(1);
    for ibl in tail_start..tail_limit {
        if !ibls.contains(&ibl) {
            ibls.push(ibl);
        }
    }
    for ibl in ibls {
        let row = &rows[ibl];
        eprintln!(
            "[STMOVE DEBUG] phase={} {} ibl={} panel={} x={:.8e} uinv={:.8e} uedg={:.8e} theta={:.8e} dstr={:.8e} mass={:.8e} dw={:.8e} wake={} lam={} turb={}",
            phase,
            label,
            ibl + 1,
            row.panel_idx,
            row.x,
            row.uinv,
            row.uedg,
            row.theta,
            row.dstr,
            row.mass,
            row.dw,
            row.is_wake,
            row.is_laminar,
            row.is_turbulent
        );
    }
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

