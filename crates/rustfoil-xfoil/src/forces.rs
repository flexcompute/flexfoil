use rustfoil_bl::BlStation;

use crate::canonical_state::rows_to_stations;
use crate::state::XfoilState;
use rustfoil_solver::viscous::config::ViscousSolverConfig;
use rustfoil_solver::viscous::forces::compute_forces_from_canonical_state;
use rustfoil_solver::viscous::state::{CanonicalBlRow, XfoilLikeViscousState};

pub fn compute_panel_forces_from_gamma(
    panel_x: &[f64],
    panel_y: &[f64],
    gamma: &[f64],
    alpha: f64,
) -> (f64, f64) {
    let cp: Vec<f64> = gamma.iter().map(|&g| 1.0 - g * g).collect();
    let cosa = alpha.cos();
    let sina = alpha.sin();
    let x_ref = 0.25;

    let mut cl = 0.0;
    let mut cm = 0.0;
    for i in 0..panel_x.len().saturating_sub(1) {
        let ip = (i + 1) % panel_x.len();
        let dx = panel_x[ip] - panel_x[i];
        let dy = panel_y[ip] - panel_y[i];
        let dx_wind = dx * cosa + dy * sina;
        let dy_wind = dy * cosa - dx * sina;
        let cp_avg = 0.5 * (cp[i] + cp[ip]);
        let dg = cp[ip] - cp[i];
        cl += cp_avg * dx_wind;

        let x_mid = 0.5 * (panel_x[i] + panel_x[ip]);
        let y_mid = 0.5 * (panel_y[i] + panel_y[ip]);
        let ax = (x_mid - x_ref) * cosa + y_mid * sina;
        let ay = y_mid * cosa - (x_mid - x_ref) * sina;
        cm -= cp_avg * (ax * dx_wind + ay * dy_wind);
        cm -= dg * dx_wind * dx_wind / 12.0;
        cm -= dg * dy_wind * dy_wind / 12.0;
    }
    (cl, cm)
}

pub fn update_force_state(state: &mut XfoilState, mach: f64, reynolds: f64) {
    let (cl, cm) =
        compute_panel_forces_from_gamma(&state.panel_x, &state.panel_y, &state.gam, state.alpha_rad);
    let msq = mach * mach;
    let upper_len = state.nbl_upper.min(state.upper_rows.len());
    let lower_len = state.nbl_lower.min(state.lower_rows.len());
    let upper_stations = rows_to_stations(&state.upper_rows[..upper_len], msq, reynolds);
    let lower_stations = rows_to_stations(&state.lower_rows[..lower_len], msq, reynolds);
    let canonical_state = canonical_force_state(state, &upper_stations, &lower_stations);
    let drag_forces = compute_forces_from_canonical_state(
        &canonical_state,
        &state.panel_x,
        &state.panel_y,
        state.alpha_rad,
        &ViscousSolverConfig::default(),
    );

    state.cl = cl;
    state.cm = cm;
    state.cdf = drag_forces.cd_friction;
    state.cdp = drag_forces.cd_pressure;
    state.cd = drag_forces.cd;
}

fn canonical_force_state(
    state: &XfoilState,
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
) -> XfoilLikeViscousState {
    let mut canonical = XfoilLikeViscousState::new(state.n_panel_nodes());
    canonical.gam = state.gam.clone();
    canonical.gam_a = state.gam_a.clone();
    canonical.upper_rows = upper_stations.iter().map(canonical_row_from_station).collect();
    canonical.lower_rows = lower_stations.iter().map(canonical_row_from_station).collect();
    canonical
}

fn canonical_row_from_station(station: &BlStation) -> CanonicalBlRow {
    CanonicalBlRow {
        x: station.x,
        x_coord: station.x_coord,
        panel_idx: station.panel_idx,
        uedg: station.u,
        theta: station.theta,
        dstr: station.delta_star,
        ctau_or_ampl: if station.is_laminar && !station.is_turbulent {
            station.ampl
        } else {
            station.ctau
        },
        ctau: station.ctau,
        ampl: station.ampl,
        mass: station.mass_defect,
        h: station.h,
        hk: station.hk,
        hs: station.hs,
        hc: station.hc,
        r_theta: station.r_theta,
        cf: station.cf,
        cd: station.cd,
        us: station.us,
        cq: station.cq,
        de: station.de,
        dw: station.dw,
        is_laminar: station.is_laminar,
        is_turbulent: station.is_turbulent,
        is_wake: station.is_wake,
        derivs: station.derivs.clone(),
    }
}
