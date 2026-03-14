use rustfoil_coupling::global_newton::{emit_blsolv_solution_debug, solve_global_system};

use crate::{assembly::AssemblyState, state::XfoilState};

#[derive(Debug, Clone)]
pub struct SolveState {
    pub state_deltas: Vec<[f64; 3]>,
    pub operating_deltas: Vec<[f64; 3]>,
    pub rms: f64,
    pub max: f64,
}

pub fn blsolv(state: &mut XfoilState, assembly: &mut AssemblyState, iteration: usize) -> SolveState {
    let rms = assembly.system.rms_residual();
    let max = assembly.system.max_residual();
    let result = solve_global_system(&mut assembly.system);
    emit_blsolv_solution_debug(iteration, &result.state_deltas);
    if std::env::var("RUSTFOIL_SOLVE_DEBUG").is_ok() {
        for iv in 1..=5.min(result.state_deltas.len().saturating_sub(1)) {
            eprintln!(
                "[SOLVE DEBUG] iv={iv} state={:?} operating={:?}",
                result.state_deltas[iv],
                result.operating_deltas[iv]
            );
        }
    }

    for iv in 0..=assembly.system.nsys.min(state.nsys) {
        for eq in 0..3 {
            state.vdel[iv][eq][0] = result.state_deltas[iv][eq];
            state.vdel[iv][eq][1] = result.operating_deltas[iv][eq];
        }
    }

    SolveState {
        state_deltas: result.state_deltas,
        operating_deltas: result.operating_deltas,
        rms,
        max,
    }
}
