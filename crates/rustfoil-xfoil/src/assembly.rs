use rustfoil_bl::FlowType;
use rustfoil_coupling::global_newton::GlobalNewtonSystem;

use crate::{
    canonical_state::{flow_type_for_station_index, recompute_mass_from_state, rows_to_stations},
    config::OperatingMode,
    march::{mrchdu, mrchue},
    state::{XfoilBlRow, XfoilState},
    state_ops::ueset,
};

#[derive(Debug, Clone)]
pub struct AssemblyState {
    pub system: GlobalNewtonSystem,
}

pub fn setbl(
    state: &mut XfoilState,
    reynolds: f64,
    ncrit: f64,
    mach: f64,
    iteration: usize,
) -> AssemblyState {
    if !state.lblini {
        mrchue(state, reynolds, ncrit, iteration);
        if std::env::var("RUSTFOIL_SETBL_HANDOFF_DEBUG").is_ok() {
            let start = state.nbl_upper.saturating_sub(3);
            for ibl in start..state.nbl_upper {
                if let Some(row) = state.upper_rows.get(ibl) {
                    eprintln!(
                        "[SETBL AFTER MRCHUE] upper[{ibl}] panel={} x={:.12e} ue={:.12e} theta={:.12e} dstr={:.12e} mass={:.12e} ampl={:.12e} ctau={:.12e}",
                        row.panel_idx,
                        row.x,
                        row.uedg,
                        row.theta,
                        row.dstr,
                        row.mass,
                        row.ampl,
                        row.ctau,
                    );
                }
            }
        }
    }
    mrchdu(state, reynolds, ncrit, iteration);

    // XFOIL re-enters SETBL with MASS coherent to the current BL state.
    // Keep the row cache aligned before building canonical stations or the
    // mass-coupled UESET reference state used for forced-change assembly.
    recompute_mass_from_state(state);

    state.allocate_newton_state();
    state.usav_upper = state.upper_rows.iter().map(|row| row.uedg).collect();
    state.usav_lower = state.lower_rows.iter().map(|row| row.uedg).collect();

    let upper_rows = &state.upper_rows[..state.nbl_upper.min(state.upper_rows.len())];
    let lower_rows = &state.lower_rows[..state.nbl_lower.min(state.lower_rows.len())];

    let upper_stations = rows_to_stations(
        upper_rows,
        state.iblte_upper,
        state.itran_upper,
        mach * mach,
        reynolds,
    );
    let lower_stations = rows_to_stations(
        lower_rows,
        state.iblte_lower,
        state.itran_lower,
        mach * mach,
        reynolds,
    );
    let upper_ue_current: Vec<f64> = upper_rows.iter().map(|row| row.uedg).collect();
    let lower_ue_current: Vec<f64> = lower_rows.iter().map(|row| row.uedg).collect();
    let upper_flow = rows_to_flow_types(upper_rows, state.iblte_upper, state.itran_upper);
    let lower_flow = rows_to_flow_types(lower_rows, state.iblte_lower, state.itran_lower);
    let upper_ue_inv: Vec<f64> = upper_rows.iter().map(|row| row.uinv).collect();
    let lower_ue_inv: Vec<f64> = lower_rows.iter().map(|row| row.uinv).collect();
    let mut mass_state = state.clone();
    ueset(&mut mass_state);
    let upper_ue_mass: Vec<f64> = mass_state.upper_rows[..state.nbl_upper.min(mass_state.upper_rows.len())]
        .iter()
        .map(|row| row.uedg)
        .collect();
    let lower_ue_mass: Vec<f64> = mass_state.lower_rows[..state.nbl_lower.min(mass_state.lower_rows.len())]
        .iter()
        .map(|row| row.uedg)
        .collect();
    if std::env::var("RUSTFOIL_SETBL_DEBUG").is_ok() {
        for (name, stations) in [("upper", &upper_stations), ("lower", &lower_stations)] {
            for (ibl, station) in stations.iter().enumerate().take(6) {
                eprintln!(
                    "[SETBL DEBUG] iter={} {name}[{ibl}] panel={} x={:.6e} u={:.6e} theta={:.6e} dstar={:.6e} mass={:.6e} ampl={:.6e} ctau={:.6e} wake={} turb={}",
                    iteration,
                    station.panel_idx,
                    station.x,
                    station.u,
                    station.theta,
                    station.delta_star,
                    station.mass_defect,
                    station.ampl,
                    station.ctau,
                    station.is_wake,
                    station.is_turbulent,
                );
            }
        }
        for ibl in 0..6.min(upper_stations.len()) {
            let station = &upper_stations[ibl];
            eprintln!(
                "[SETBL DEBUG EXACT] iter={} upper[{ibl}] panel={} x={:.12e} u={:.12e} theta={:.12e} dstar={:.12e} mass={:.12e} ampl={:.12e} ctau={:.12e}",
                iteration,
                station.panel_idx,
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.mass_defect,
                station.ampl,
                station.ctau,
            );
        }
        let upper_airfoil_start = upper_stations.len().saturating_sub(6);
        for ibl in upper_airfoil_start..upper_stations.len() {
            let station = &upper_stations[ibl];
            let row = upper_rows.get(ibl);
            eprintln!(
                "[SETBL DEBUG EXACT TAIL] iter={} upper[{ibl}] panel={} x={:.12e} u={:.12e} theta={:.12e} dstar={:.12e} mass={:.12e} dw={:.12e} ampl={:.12e} ctau={:.12e}",
                iteration,
                station.panel_idx,
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.mass_defect,
                row.map(|r| r.dw).unwrap_or(0.0),
                station.ampl,
                station.ctau,
            );
        }
        for ibl in 28..32.min(upper_stations.len()) {
            let station = &upper_stations[ibl];
            eprintln!(
                "[SETBL DEBUG] iter={} upper_transition[{ibl}] panel={} x={:.6e} u={:.6e} theta={:.6e} dstar={:.6e} mass={:.6e} ampl={:.6e} ctau={:.6e} wake={} turb={}",
                iteration,
                station.panel_idx,
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.mass_defect,
                station.ampl,
                station.ctau,
                station.is_wake,
                station.is_turbulent,
            );
        }
        for ibl in 60..65.min(lower_stations.len()) {
            let station = &lower_stations[ibl];
            eprintln!(
                "[SETBL DEBUG] iter={} lower_tail[{ibl}] panel={} x={:.6e} u={:.6e} theta={:.6e} dstar={:.6e} mass={:.6e} dw={:.6e} ampl={:.6e} ctau={:.6e} wake={} turb={}",
                iteration,
                station.panel_idx,
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.mass_defect,
                station.dw,
                station.ampl,
                station.ctau,
                station.is_wake,
                station.is_turbulent,
            );
        }
        for ibl in 54..59.min(lower_stations.len()) {
            let station = &lower_stations[ibl];
            eprintln!(
                "[SETBL DEBUG] iter={} lower_mid[{ibl}] panel={} x={:.6e} u={:.6e} theta={:.6e} dstar={:.6e} mass={:.6e} dw={:.6e} ampl={:.6e} ctau={:.6e} wake={} turb={}",
                iteration,
                station.panel_idx,
                station.x,
                station.u,
                station.theta,
                station.delta_star,
                station.mass_defect,
                station.dw,
                station.ampl,
                station.ctau,
                station.is_wake,
                station.is_turbulent,
            );
        }
        for (name, values) in [("upper_mass_ue", &upper_ue_mass), ("lower_mass_ue", &lower_ue_mass)] {
            for (ibl, value) in values.iter().enumerate().take(6) {
                eprintln!("[SETBL DEBUG] iter={} {name}[{ibl}] ue={value:.6e}", iteration);
            }
        }
    }
    let (upper_ue_operating, lower_ue_operating) = match state.operating_mode {
        OperatingMode::PrescribedAlpha => (
            vec![0.0; upper_rows.len()],
            vec![0.0; lower_rows.len()],
        ),
        OperatingMode::PrescribedCl { .. } => (
            upper_rows.iter().map(|row| row.uinv_a).collect(),
            lower_rows.iter().map(|row| row.uinv_a).collect(),
        ),
    };

    let mut system = GlobalNewtonSystem::new(
        state.nbl_upper,
        state.nbl_lower,
        state.iblte_upper,
        state.iblte_lower,
    );
    system.ante = state.ante;
    system.set_stagnation_derivs(state.sst_go, state.sst_gp);
    system.build_global_system(
        &upper_stations,
        &lower_stations,
        &upper_ue_current,
        &lower_ue_current,
        &upper_flow,
        &lower_flow,
        &state.dij,
        &upper_ue_inv,
        &lower_ue_inv,
        &upper_ue_mass,
        &lower_ue_mass,
        &upper_ue_operating,
        &lower_ue_operating,
        ncrit,
        mach * mach,
        reynolds,
        iteration,
    );
    system.emit_setbl_debug(iteration);

    copy_system_to_state(state, &system);
    AssemblyState { system }
}

fn rows_to_flow_types(rows: &[XfoilBlRow], iblte: usize, itran: usize) -> Vec<FlowType> {
    rows.iter()
        .enumerate()
        .skip(1)
        .map(|(ibl, _)| flow_type_for_station_index(ibl, iblte, itran))
        .collect()
}

fn copy_system_to_state(state: &mut XfoilState, system: &GlobalNewtonSystem) {
    state.nsys = system.nsys;
    state.vz = system.vz;
    for iv in 0..=system.nsys {
        for eq in 0..3 {
            state.va[iv][eq][0] = system.va[iv][eq][0];
            state.va[iv][eq][1] = system.va[iv][eq][1];
            state.vb[iv][eq][0] = system.vb[iv][eq][0];
            state.vb[iv][eq][1] = system.vb[iv][eq][1];
            state.vdel[iv][eq][0] = system.vdel[iv][eq];
            state.vdel[iv][eq][1] = system.vdel_operating[iv][eq];
        }
        for jv in 0..=system.nsys {
            state.vm[iv][jv] = system.vm[iv][jv];
        }
    }
    state.due1 = system.dule1;
    state.due2 = system.dule2;
}
