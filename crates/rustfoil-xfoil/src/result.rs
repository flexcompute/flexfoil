use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum QdesTargetKind {
    EdgeVelocity,
    PressureCoefficient,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QdesTarget {
    pub x: Vec<f64>,
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QdesSpec {
    pub operating_point: crate::oper::AlphaSpec,
    pub target_kind: QdesTargetKind,
    pub upper: Option<QdesTarget>,
    pub lower: Option<QdesTarget>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceDistribution {
    pub x: Vec<f64>,
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QdesIterationSnapshot {
    pub iteration: usize,
    pub rms_error: f64,
    pub max_error: f64,
    pub geometry_delta_norm: f64,
    pub oper_result: XfoilViscousResult,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QdesResult {
    pub input_coords: Vec<(f64, f64)>,
    pub output_coords: Vec<(f64, f64)>,
    pub paneled_coords: Vec<(f64, f64)>,
    pub oper_result: XfoilViscousResult,
    pub converged: bool,
    pub iterations: usize,
    pub rms_error: f64,
    pub max_error: f64,
    pub target_kind: QdesTargetKind,
    pub target_upper: Option<QdesTarget>,
    pub target_lower: Option<QdesTarget>,
    pub achieved_upper: Option<SurfaceDistribution>,
    pub achieved_lower: Option<SurfaceDistribution>,
    pub history: Vec<QdesIterationSnapshot>,
}
