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

pub fn update_force_state(state: &mut XfoilState) {
    let (cl, cm) =
        compute_panel_forces_from_gamma(&state.panel_x, &state.panel_y, &state.gam, state.alpha_rad);
    let canonical_state = canonical_force_state(state);
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

fn canonical_force_state(state: &XfoilState) -> XfoilLikeViscousState {
    let mut canonical = XfoilLikeViscousState::new(state.n_panel_nodes());
    canonical.gam = state.gam.clone();
    canonical.gam_a = state.gam_a.clone();
    canonical.upper_rows = state.upper_rows[..state.nbl_upper.min(state.upper_rows.len())]
        .iter()
        .map(canonical_row_from_xfoil)
        .collect();
    canonical.lower_rows = state.lower_rows[..state.nbl_lower.min(state.lower_rows.len())]
        .iter()
        .map(canonical_row_from_xfoil)
        .collect();
    canonical
}

fn canonical_row_from_xfoil(row: &crate::state::XfoilBlRow) -> CanonicalBlRow {
    CanonicalBlRow {
        x: row.x,
        x_coord: row.x_coord,
        panel_idx: row.panel_idx,
        uedg: row.uedg,
        theta: row.theta,
        dstr: row.dstr,
        ctau_or_ampl: if row.is_laminar && !row.is_turbulent {
            row.ampl
        } else {
            row.ctau
        },
        ctau: row.ctau,
        ampl: row.ampl,
        mass: row.mass,
        h: row.h,
        hk: row.hk,
        hs: row.hs,
        hc: row.hc,
        r_theta: row.r_theta,
        cf: row.cf,
        cd: row.cd,
        us: row.us,
        cq: row.cq,
        de: row.de,
        dw: row.dw,
        is_laminar: row.is_laminar,
        is_turbulent: row.is_turbulent,
        is_wake: row.is_wake,
        derivs: row.derivs.clone(),
    }
}
