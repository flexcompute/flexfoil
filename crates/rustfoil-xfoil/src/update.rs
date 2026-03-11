use crate::{
    config::OperatingMode,
    solve::SolveState,
    state::XfoilState,
};

pub fn update(state: &mut XfoilState, solve: &SolveState) {
    state.residual = solve.rms;
    state.rlx = 1.0;
    state.u_new.clear();
    state.u_ac.clear();

    let mut delta_iter = solve.deltas.iter();
    for row in state.upper_rows.iter_mut().skip(1) {
        let delta = delta_iter.next().copied().unwrap_or(0.0);
        row.mass = (row.mass + delta).max(0.0);
        if row.uedg.abs() > 1.0e-12 {
            row.dstr = row.mass / row.uedg;
        }
        state.u_new.push(row.uedg);
        state.u_ac.push(0.0);
    }
    for row in state.lower_rows.iter_mut().skip(1) {
        let delta = delta_iter.next().copied().unwrap_or(0.0);
        row.mass = (row.mass + delta).max(0.0);
        if row.uedg.abs() > 1.0e-12 {
            row.dstr = row.mass / row.uedg;
        }
        state.u_new.push(row.uedg);
        state.u_ac.push(0.0);
    }

    state.dac = match state.operating_mode {
        OperatingMode::PrescribedAlpha => 0.0,
        OperatingMode::PrescribedCl { target_cl } => 0.1 * (target_cl - state.cl),
    };
    if let OperatingMode::PrescribedCl { .. } = state.operating_mode {
        state.alpha_rad += state.dac;
    }
}
