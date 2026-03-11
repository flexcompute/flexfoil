#[derive(Debug, Clone)]
pub struct XfoilViscousResult {
    pub alpha_deg: f64,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub x_tr_upper: f64,
    pub x_tr_lower: f64,
    pub converged: bool,
    pub iterations: usize,
    pub residual: f64,
    pub cd_friction: f64,
    pub cd_pressure: f64,
    pub x_separation: Option<f64>,
}
