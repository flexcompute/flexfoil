use rustfoil_bl::{blvar, BlStation, FlowType};

use crate::state::{XfoilBlRow, XfoilState, XfoilSurface};

pub(crate) fn flow_type_for_row(row: &XfoilBlRow) -> FlowType {
    if row.is_wake {
        FlowType::Wake
    } else if row.is_turbulent {
        FlowType::Turbulent
    } else {
        FlowType::Laminar
    }
}

pub(crate) fn rows_to_stations(rows: &[XfoilBlRow], msq: f64, reynolds: f64) -> Vec<BlStation> {
    rows.iter()
        .map(|row| {
            let mut station = BlStation::new();
            let total_dstr = row.dstr.max(1.0e-12);
            station.x = row.x;
            station.x_coord = row.x_coord;
            station.panel_idx = row.panel_idx;
            station.u = row.uedg.abs().max(1.0e-9);
            station.theta = row.theta.max(1.0e-12);
            station.dw = row.dw.max(0.0);
            station.delta_star = if row.is_wake {
                (total_dstr - station.dw).max(1.0e-12)
            } else {
                total_dstr
            };
            station.ctau = if row.is_laminar {
                0.03
            } else {
                row.ctau.max(1.0e-7)
            };
            station.ampl = row.ampl.max(0.0);
            station.mass_defect = row.mass.max(1.0e-12);
            station.is_laminar = row.is_laminar;
            station.is_turbulent = row.is_turbulent;
            station.is_wake = row.is_wake;
            blvar(&mut station, flow_type_for_row(row), msq, reynolds);
            station
        })
        .collect()
}

pub(crate) fn sync_rows_from_stations(rows: &mut [XfoilBlRow], stations: &[BlStation]) {
    for (row, station) in rows.iter_mut().zip(stations.iter()) {
        row.uedg = station.u;
        row.theta = station.theta;
        row.dstr = station.delta_star;
        row.mass = station.mass_defect;
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
    let n_panel_nodes = state.n_panel_nodes();
    let n_total_nodes = state.n_total_nodes();
    state.qvis.clear();
    state.qvis.resize(n_total_nodes, 0.0);

    for row in state.upper_rows.iter().skip(1).take(state.nbl_upper.saturating_sub(1)) {
        if !row.is_wake && row.panel_idx < state.qvis.len() {
            state.qvis[row.panel_idx] = XfoilSurface::Upper.vti() * row.uedg;
        }
    }
    for row in state.lower_rows.iter().skip(1).take(state.nbl_lower.saturating_sub(1)) {
        if !row.is_wake && row.panel_idx < state.qvis.len() {
            state.qvis[row.panel_idx] = XfoilSurface::Lower.vti() * row.uedg;
        } else if row.is_wake {
            let wake_idx = row.panel_idx.saturating_sub(n_panel_nodes);
            let qvis_idx = n_panel_nodes + wake_idx;
            if qvis_idx < state.qvis.len() {
                state.qvis[qvis_idx] = row.uedg.abs();
            }
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
