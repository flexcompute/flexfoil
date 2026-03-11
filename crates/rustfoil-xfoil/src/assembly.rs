use crate::state::XfoilState;

#[derive(Debug, Clone)]
pub struct AssemblyState {
    pub residuals: Vec<f64>,
}

pub fn setbl(state: &XfoilState) -> AssemblyState {
    let residuals = state
        .upper_rows
        .iter()
        .chain(state.lower_rows.iter())
        .skip(1)
        .map(|row| (row.mass - row.uedg * row.dstr).abs())
        .collect();
    AssemblyState { residuals }
}
