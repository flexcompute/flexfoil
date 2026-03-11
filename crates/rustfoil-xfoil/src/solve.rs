use crate::assembly::AssemblyState;

#[derive(Debug, Clone)]
pub struct SolveState {
    pub deltas: Vec<f64>,
    pub rms: f64,
    pub max: f64,
}

pub fn blsolv(assembly: &AssemblyState) -> SolveState {
    let deltas = assembly.residuals.iter().map(|&r| -0.25 * r).collect::<Vec<_>>();
    let rms = if assembly.residuals.is_empty() {
        0.0
    } else {
        (assembly
            .residuals
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            / assembly.residuals.len() as f64)
            .sqrt()
    };
    let max = assembly.residuals.iter().copied().fold(0.0, f64::max);
    SolveState { deltas, rms, max }
}
