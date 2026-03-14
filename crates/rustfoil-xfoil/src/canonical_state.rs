use rustfoil_bl::{blvar, BlStation, FlowType};

use crate::state::{XfoilBlRow, XfoilState, XfoilSurface};

pub(crate) fn flow_type_for_station_index(ibl: usize, iblte: usize, itran: usize) -> FlowType {
    if ibl > iblte {
        FlowType::Wake
    } else if ibl + 1 >= itran {
        FlowType::Turbulent
    } else {
        FlowType::Laminar
    }
}

pub(crate) fn rows_to_stations(
    rows: &[XfoilBlRow],
    iblte: usize,
    itran: usize,
    msq: f64,
    reynolds: f64,
) -> Vec<BlStation> {
    rows.iter()
        .enumerate()
        .map(|(ibl, row)| {
            let flow_type = flow_type_for_station_index(ibl, iblte, itran);
            let is_wake = matches!(flow_type, FlowType::Wake);
            let is_turbulent = matches!(flow_type, FlowType::Turbulent | FlowType::Wake);
            let is_laminar = matches!(flow_type, FlowType::Laminar);
            let mut station = BlStation::new();
            let total_dstr = row.dstr.max(1.0e-12);
            station.x = row.x;
            station.x_coord = row.x_coord;
            station.panel_idx = row.panel_idx;
            station.u = row.uedg.abs().max(1.0e-9);
            station.theta = row.theta.max(1.0e-12);
            station.dw = row.dw.max(0.0);
            station.delta_star = if is_wake {
                (total_dstr - station.dw).max(1.0e-12)
            } else {
                total_dstr
            };
            station.ctau = if is_laminar {
                0.03
            } else {
                row.ctau.max(1.0e-7)
            };
            station.ampl = row.ampl.max(0.0);
            station.mass_defect = row.mass.max(1.0e-12);
            station.is_laminar = is_laminar;
            station.is_turbulent = is_turbulent;
            station.is_wake = is_wake;
            blvar(&mut station, flow_type, msq, reynolds);
            station
        })
        .collect()
}

pub(crate) fn sync_rows_from_stations(rows: &mut [XfoilBlRow], stations: &[BlStation]) {
    for (row, station) in rows.iter_mut().zip(stations.iter()) {
        row.uedg = station.u;
        row.theta = station.theta;
        row.dw = if station.is_wake {
            station.dw.max(0.0)
        } else {
            0.0
        };
        row.dstr = if station.is_wake {
            station.delta_star + station.dw
        } else {
            station.delta_star
        };
        row.ampl = station.ampl;
        row.ctau = if station.is_laminar {
            station.ampl
        } else {
            station.ctau
        };
        row.h = station.h;
        row.hk = station.hk;
        row.hs = station.hs;
        row.hc = station.hc;
        row.r_theta = station.r_theta;
        row.cf = station.cf;
        row.cd = station.cd;
        row.us = station.us;
        row.cq = station.cq;
        row.de = station.de;
        row.mass = station.mass_defect;
        row.is_laminar = station.is_laminar;
        row.is_turbulent = station.is_turbulent;
        row.is_wake = station.is_wake;
        row.derivs = station.derivs.clone();
    }
}

pub(crate) fn mirror_upper_wake_from_lower(state: &mut XfoilState) {
    let upper_wake_start = state.nbl_upper;
    let lower_wake_start = state.iblte_lower + 1;
    for iw in 0..state.wake_x.len() {
        let upper_idx = upper_wake_start + iw;
        let lower_idx = lower_wake_start + iw;
        if upper_idx >= state.upper_rows.len() || lower_idx >= state.lower_rows.len() {
            break;
        }
        let source = state.lower_rows[lower_idx].clone();
        let panel_idx = state.upper_rows[upper_idx].panel_idx;
        let x_coord = state.upper_rows[upper_idx].x_coord;
        state.upper_rows[upper_idx] = source;
        state.upper_rows[upper_idx].panel_idx = panel_idx;
        state.upper_rows[upper_idx].x_coord = x_coord;
        state.upper_rows[upper_idx].is_wake = true;
        state.upper_rows[upper_idx].is_turbulent = true;
        state.upper_rows[upper_idx].is_laminar = false;
    }
}

pub(crate) fn rebuild_qvis_from_rows(state: &mut XfoilState) {
    let n_total_nodes = state.n_total_nodes();
    state.qvis.clear();
    state.qvis.resize(n_total_nodes, 0.0);

    for ibl in 1..state.nbl_upper.min(state.upper_rows.len()).min(state.ipan_upper.len()) {
        let panel_idx = state.ipan_upper[ibl];
        if panel_idx < state.qvis.len() {
            state.qvis[panel_idx] = XfoilSurface::Upper.vti() * state.upper_rows[ibl].uedg;
        }
    }
    for ibl in 1..state.nbl_lower.min(state.lower_rows.len()).min(state.ipan_lower.len()) {
        let panel_idx = state.ipan_lower[ibl];
        if panel_idx < state.qvis.len() {
            state.qvis[panel_idx] = XfoilSurface::Lower.vti() * state.lower_rows[ibl].uedg;
        }
    }
}

pub(crate) fn rebuild_gam_from_qvis(state: &mut XfoilState) {
    state.gam.resize(state.n_panel_nodes(), 0.0);
    state.gam_a.resize(state.n_panel_nodes(), 0.0);
    for i in 0..state.n_panel_nodes() {
        state.gam[i] = state.qvis.get(i).copied().unwrap_or(0.0);
        state.gam_a[i] = state.qinv_a.get(i).copied().unwrap_or(0.0);
    }
}

pub(crate) fn recompute_mass_from_state(state: &mut XfoilState) {
    for row in state.upper_rows.iter_mut().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        row.mass = row.dstr * row.uedg;
    }
    for row in state.lower_rows.iter_mut().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        row.mass = row.dstr * row.uedg;
    }
}

pub(crate) fn sync_state_views_from_rows(state: &mut XfoilState) {
    mirror_upper_wake_from_lower(state);
    rebuild_qvis_from_rows(state);
    rebuild_gam_from_qvis(state);
}
