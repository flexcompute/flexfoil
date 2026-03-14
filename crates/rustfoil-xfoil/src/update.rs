use rustfoil_bl::{blvar, BlStation, FlowType};
use rustfoil_coupling::global_newton::{apply_global_updates, preview_global_update_ue, GlobalNewtonSystem};

use crate::{
    assembly::AssemblyState,
    canonical_state::{rows_to_stations, sync_rows_from_stations, sync_state_views_from_rows},
    config::OperatingMode,
    forces::compute_panel_forces_from_gamma,
    solve::SolveState,
    state::XfoilState,
};

pub fn update(
    state: &mut XfoilState,
    assembly: &AssemblyState,
    solve: &SolveState,
    mach: f64,
    reynolds: f64,
    iteration: usize,
) {
    let msq = mach * mach;
    let upper_len = state.nbl_upper.min(state.upper_rows.len());
    let lower_len = state.nbl_lower.min(state.lower_rows.len());
    let mut upper_stations = rows_to_stations(
        &state.upper_rows[..upper_len],
        state.iblte_upper,
        state.itran_upper,
        msq,
        reynolds,
    );
    let mut lower_stations = rows_to_stations(
        &state.lower_rows[..lower_len],
        state.iblte_lower,
        state.itran_lower,
        msq,
        reynolds,
    );
    let upper_before = upper_stations.clone();
    let lower_before = lower_stations.clone();
    let upper_ue_inv: Vec<f64> = state.upper_rows[..upper_len].iter().map(|row| row.uinv).collect();
    let lower_ue_inv: Vec<f64> = state.lower_rows[..lower_len].iter().map(|row| row.uinv).collect();
    let (upper_ue_operating, lower_ue_operating) = match state.operating_mode {
        OperatingMode::PrescribedAlpha => (
            vec![0.0; upper_len],
            vec![0.0; lower_len],
        ),
        OperatingMode::PrescribedCl { .. } => (
            state.upper_rows[..upper_len].iter().map(|row| row.uinv_a).collect(),
            state.lower_rows[..lower_len].iter().map(|row| row.uinv_a).collect(),
        ),
    };

    let preview = preview_global_update_ue(
        &upper_stations,
        &lower_stations,
        &solve.state_deltas,
        Some(&solve.operating_deltas),
        &assembly.system,
        Some(&state.dij),
        Some(&upper_ue_inv),
        Some(&lower_ue_inv),
        Some(&upper_ue_operating),
        Some(&lower_ue_operating),
    );
    let (cl_new, cl_a, cl_ac, q_new, q_ac) = compute_trial_cl_terms(
        state,
        &upper_stations,
        &lower_stations,
        &preview.upper_u_new,
        &preview.lower_u_new,
        &preview.upper_u_ac,
        &preview.lower_u_ac,
    );
    let dac = compute_operating_correction(
        state.operating_mode,
        state.cl,
        cl_new,
        cl_a,
        cl_ac,
    );
    let initial_relaxation = relaxation_limit_from_operating_correction(state.operating_mode, dac);
    let update_result = apply_global_updates(
        &mut upper_stations,
        &mut lower_stations,
        &solve.state_deltas,
        Some(&solve.operating_deltas),
        &assembly.system,
        iteration,
        initial_relaxation,
        Some(&state.dij),
        Some(&upper_ue_inv),
        Some(&lower_ue_inv),
        Some(&upper_ue_operating),
        Some(&lower_ue_operating),
        dac,
        msq,
        reynolds,
    );

    refresh_secondary_vars(&mut upper_stations, msq, reynolds);
    refresh_secondary_vars(&mut lower_stations, msq, reynolds);
    emit_update_debug_events(
        &assembly.system,
        &upper_before,
        &lower_before,
        &upper_stations,
        &lower_stations,
        &preview.upper_u_new,
        &preview.lower_u_new,
        &preview.upper_u_ac,
        &preview.lower_u_ac,
        solve,
        dac,
        update_result.relaxation_used,
        iteration,
    );
    sync_rows_from_stations(&mut state.upper_rows[..upper_len], &upper_stations);
    sync_rows_from_stations(&mut state.lower_rows[..lower_len], &lower_stations);
    sync_state_views_from_rows(state);
    if std::env::var("RUSTFOIL_UPDATE_DEBUG").is_ok() {
        eprintln!(
            "[UPDATE DEBUG] cl_new={:.6e} cl_a={:.6e} cl_ac={:.6e} dac={:.6e} rlx={:.6e}",
            cl_new, cl_a, cl_ac, dac, update_result.relaxation_used
        );
        for (name, rows) in [("upper", &state.upper_rows), ("lower", &state.lower_rows)] {
            for (ibl, row) in rows.iter().enumerate().take(6) {
                eprintln!(
                    "[UPDATE DEBUG] {name}[{ibl}] u={:.6e} theta={:.6e} dstar={:.6e} mass={:.6e} ctau={:.6e} ampl={:.6e}",
                    row.uedg,
                    row.theta,
                    row.dstr,
                    row.mass,
                    row.ctau,
                    row.ampl,
                );
            }
        }
    }

    state.dac = dac;
    if let OperatingMode::PrescribedCl { .. } = state.operating_mode {
        state.alpha_rad += dac;
    }
    state.rlx = update_result.relaxation_used;
    state.rmsbl = update_result.rms_change;
    state.rmxbl = update_result.max_change;
    state.residual = update_result.rms_change;
    state.cl_new = cl_new;
    state.cl_a = cl_a;
    state.cl_ac = cl_ac;
    state.cl_ms = 0.0;
    state.q_new = q_new;
    state.q_ac = q_ac;
    state.u_new = flatten_surface_values(
        &assembly.system,
        &update_result.upper_u_new,
        &update_result.lower_u_new,
    );
    state.u_ac = flatten_surface_values(
        &assembly.system,
        &update_result.upper_u_ac,
        &update_result.lower_u_ac,
    );
}

fn refresh_secondary_vars(stations: &mut [BlStation], msq: f64, reynolds: f64) {
    for station in stations.iter_mut().skip(1) {
        blvar(
            station,
            if station.is_wake {
                FlowType::Wake
            } else if station.is_turbulent {
                FlowType::Turbulent
            } else {
                FlowType::Laminar
            },
            msq,
            reynolds,
        );
    }
}

fn flatten_surface_values(
    system: &GlobalNewtonSystem,
    upper: &[f64],
    lower: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; system.nsys + 1];
    for ibl in 1..upper.len() {
        let iv = system.to_global(0, ibl);
        if iv < out.len() {
            out[iv] = upper[ibl];
        }
    }
    for ibl in 1..lower.len() {
        let iv = system.to_global(1, ibl);
        if iv < out.len() {
            out[iv] = lower[ibl];
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn emit_update_debug_events(
    system: &GlobalNewtonSystem,
    upper_before: &[BlStation],
    lower_before: &[BlStation],
    upper_after: &[BlStation],
    lower_after: &[BlStation],
    upper_u_new: &[f64],
    lower_u_new: &[f64],
    upper_u_ac: &[f64],
    lower_u_ac: &[f64],
    solve: &SolveState,
    dac: f64,
    relaxation: f64,
    iteration: usize,
) {
    if !rustfoil_bl::is_debug_active() {
        return;
    }

    emit_update_debug_surface(
        system,
        1,
        upper_before,
        upper_after,
        upper_u_new,
        upper_u_ac,
        solve,
        dac,
        relaxation,
        iteration,
    );
    emit_update_debug_surface(
        system,
        2,
        lower_before,
        lower_after,
        lower_u_new,
        lower_u_ac,
        solve,
        dac,
        relaxation,
        iteration,
    );
}

#[allow(clippy::too_many_arguments)]
fn emit_update_debug_surface(
    system: &GlobalNewtonSystem,
    side: usize,
    before: &[BlStation],
    after: &[BlStation],
    u_new: &[f64],
    u_ac: &[f64],
    solve: &SolveState,
    dac: f64,
    relaxation: f64,
    iteration: usize,
) {
    let surface = side - 1;
    for ibl in 1..before.len() {
        let iv = system.to_global(surface, ibl);
        if iv >= solve.state_deltas.len() || iv >= solve.operating_deltas.len() {
            break;
        }

        let before_station = &before[ibl];
        let after_station = &after[ibl];
        let raw_dctau = solve.state_deltas[iv][0] - dac * solve.operating_deltas[iv][0];
        let raw_dtheta = solve.state_deltas[iv][1] - dac * solve.operating_deltas[iv][1];
        let raw_dmass = solve.state_deltas[iv][2] - dac * solve.operating_deltas[iv][2];
        let raw_due = u_new
            .get(ibl)
            .zip(u_ac.get(ibl))
            .map(|(&new_u, &ac_u)| new_u - dac * ac_u - before_station.u)
            .unwrap_or(0.0);
        let ctau_before = if before_station.is_laminar {
            before_station.ampl
        } else {
            before_station.ctau
        };
        let ctau_after = if after_station.is_laminar {
            after_station.ampl
        } else {
            after_station.ctau
        };

        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::update(
            iteration,
            side,
            ibl + 1,
            relaxation * raw_dctau,
            relaxation * raw_dtheta,
            relaxation * raw_dmass,
            relaxation * raw_due,
            relaxation,
        ));
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::update_detailed(
            iteration,
            side,
            ibl + 1,
            ctau_before,
            before_station.theta,
            before_station.delta_star,
            before_station.u,
            before_station.mass_defect,
            before_station.h,
            before_station.hk,
            raw_dctau,
            raw_dtheta,
            raw_dmass,
            raw_due,
            relaxation,
            ctau_after,
            after_station.theta,
            after_station.delta_star,
            after_station.u,
            after_station.mass_defect,
            after_station.h,
            after_station.hk,
        ));
    }
}

fn build_panel_values_from_surface_values(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    upper_values: &[f64],
    lower_values: &[f64],
    n_panels: usize,
) -> Vec<f64> {
    let mut panel_values = vec![0.0; n_panels];
    for (station, &value) in upper_stations.iter().zip(upper_values.iter()).skip(1) {
        if station.panel_idx < n_panels && !station.is_wake {
            panel_values[station.panel_idx] = value;
        }
    }
    for (station, &value) in lower_stations.iter().zip(lower_values.iter()).skip(1) {
        if station.panel_idx < n_panels && !station.is_wake {
            panel_values[station.panel_idx] = -value;
        }
    }
    if let Some(last) = panel_values.last_mut() {
        *last = 0.0;
    }
    panel_values
}

fn compute_trial_cl_terms(
    state: &XfoilState,
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    upper_u_new: &[f64],
    lower_u_new: &[f64],
    upper_u_ac: &[f64],
    lower_u_ac: &[f64],
) -> (f64, f64, f64, Vec<f64>, Vec<f64>) {
    let qnew = build_panel_values_from_surface_values(
        upper_stations,
        lower_stations,
        upper_u_new,
        lower_u_new,
        state.panel_x.len(),
    );
    let qac = build_panel_values_from_surface_values(
        upper_stations,
        lower_stations,
        upper_u_ac,
        lower_u_ac,
        state.panel_x.len(),
    );
    let (cl_new, _) = compute_panel_forces_from_gamma(&state.panel_x, &state.panel_y, &qnew, state.alpha_rad);
    let eps = 1.0e-6;
    let (cl_plus, _) = compute_panel_forces_from_gamma(
        &state.panel_x,
        &state.panel_y,
        &qnew,
        state.alpha_rad + eps,
    );
    let (cl_minus, _) = compute_panel_forces_from_gamma(
        &state.panel_x,
        &state.panel_y,
        &qnew,
        state.alpha_rad - eps,
    );
    let cl_a = (cl_plus - cl_minus) / (2.0 * eps);
    let cl_ac = if qac.iter().any(|value| value.abs() > 0.0) {
        let qnew_operating: Vec<f64> = qnew
            .iter()
            .zip(qac.iter())
            .map(|(base, operating)| base + eps * operating)
            .collect();
        let (cl_operating, _) = compute_panel_forces_from_gamma(
            &state.panel_x,
            &state.panel_y,
            &qnew_operating,
            state.alpha_rad,
        );
        (cl_operating - cl_new) / eps
    } else {
        0.0
    };
    (cl_new, cl_a, cl_ac, qnew, qac)
}

fn compute_operating_correction(
    mode: OperatingMode,
    current_cl: f64,
    cl_new: f64,
    cl_a: f64,
    cl_ac: f64,
) -> f64 {
    match mode {
        OperatingMode::PrescribedAlpha => {
            let denominator = 1.0 - cl_ac;
            if denominator.abs() < 1.0e-12 {
                0.0
            } else {
                (cl_new - current_cl) / denominator
            }
        }
        OperatingMode::PrescribedCl { target_cl } => {
            let denominator = -(cl_ac + cl_a);
            if denominator.abs() < 1.0e-12 {
                0.0
            } else {
                const DALMAX: f64 = 0.5_f64.to_radians();
                const DALMIN: f64 = -0.5_f64.to_radians();
                ((cl_new - target_cl) / denominator).clamp(DALMIN, DALMAX)
            }
        }
    }
}

fn relaxation_limit_from_operating_correction(mode: OperatingMode, dac: f64) -> f64 {
    if dac.abs() <= 1.0e-12 {
        return 1.0;
    }

    let (max_delta, min_delta) = match mode {
        OperatingMode::PrescribedAlpha => (0.5, -0.5),
        OperatingMode::PrescribedCl { .. } => (0.5_f64.to_radians(), -0.5_f64.to_radians()),
    };

    let mut rlx = 1.0;
    if rlx * dac > max_delta {
        rlx = max_delta / dac;
    }
    if rlx * dac < min_delta {
        rlx = min_delta / dac;
    }
    rlx.clamp(0.0, 1.0)
}
