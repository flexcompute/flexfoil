use crate::state::XfoilState;

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
    let cdf_upper = surface_friction_drag(&state.upper_rows);
    let cdf_lower = surface_friction_drag(&state.lower_rows);
    let cdf = cdf_upper + cdf_lower;
    let cdp = wake_pressure_drag(state);

    state.cl = cl;
    state.cm = cm;
    state.cdf = cdf;
    state.cdp = cdp;
    state.cd = cdf + cdp;
}

fn surface_friction_drag(rows: &[crate::state::XfoilBlRow]) -> f64 {
    rows.windows(2)
        .map(|pair| {
            let dx = (pair[1].x_coord - pair[0].x_coord).abs().max(1.0e-6);
            0.5 * (pair[0].cf + pair[1].cf) * dx
        })
        .sum::<f64>()
        .max(1.0e-6)
}

fn wake_pressure_drag(state: &XfoilState) -> f64 {
    match (state.upper_rows.last(), state.lower_rows.last()) {
        (Some(upper), Some(lower)) => (upper.dstr + lower.dstr).abs() * 1.0e-3,
        _ => 0.0,
    }
}
