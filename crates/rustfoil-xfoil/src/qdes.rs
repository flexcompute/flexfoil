use std::f64::consts::PI;

use nalgebra::{DMatrix, DVector};

use rustfoil_core::{Body, CubicSpline, PanelingParams, Point};

use crate::{
    config::XfoilOptions,
    error::{Result, XfoilError},
    oper::{build_state_from_coords, coords_from_body, solve_operating_point_from_state},
    result::{
        QdesIterationSnapshot, QdesResult, QdesSpec, QdesTarget, QdesTargetKind,
        SurfaceDistribution, XfoilViscousResult,
    },
    state::XfoilBlRow,
};

#[derive(Debug, Clone)]
pub struct QdesOptions {
    pub xfoil_options: XfoilOptions,
    pub outer_iterations: usize,
    pub target_tolerance: f64,
    pub basis_count_per_side: usize,
    pub finite_difference_step: f64,
    pub update_damping: f64,
    pub max_basis_delta: f64,
    pub panel_count: usize,
    pub paneling: PanelingParams,
}

impl Default for QdesOptions {
    fn default() -> Self {
        Self {
            xfoil_options: XfoilOptions::default(),
            outer_iterations: 6,
            target_tolerance: 5.0e-3,
            basis_count_per_side: 4,
            finite_difference_step: 2.5e-3,
            update_damping: 0.6,
            max_basis_delta: 1.5e-2,
            panel_count: 160,
            paneling: PanelingParams::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct InverseDesignSession {
    name: String,
    design_points: Vec<Point>,
    spec: QdesSpec,
    options: QdesOptions,
    history: Vec<QdesIterationSnapshot>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceId {
    Upper,
    Lower,
}

#[derive(Debug, Clone)]
struct BasisFunction {
    surface: SurfaceId,
    center: f64,
    width: f64,
}

#[derive(Debug, Clone)]
struct ErrorStats {
    max_error: f64,
    rms_error: f64,
    error_vector: Vec<f64>,
    achieved_upper: Option<SurfaceDistribution>,
    achieved_lower: Option<SurfaceDistribution>,
}

impl QdesSpec {
    fn validate(&self) -> Result<()> {
        if self.upper.is_none() && self.lower.is_none() {
            return Err(XfoilError::Message(
                "QDES requires at least one upper or lower surface target".to_string(),
            ));
        }
        if let Some(target) = &self.upper {
            validate_target("upper", target)?;
        }
        if let Some(target) = &self.lower {
            validate_target("lower", target)?;
        }
        Ok(())
    }
}

impl InverseDesignSession {
    pub fn new(name: impl Into<String>, coords: &[(f64, f64)], spec: QdesSpec, options: QdesOptions) -> Result<Self> {
        spec.validate()?;
        options
            .xfoil_options
            .validate()
            .map_err(|msg| XfoilError::Message(msg.to_string()))?;
        if coords.len() < 4 {
            return Err(XfoilError::Message(
                "QDES requires at least four coordinate points".to_string(),
            ));
        }
        let design_points = coords
            .iter()
            .map(|&(x, y)| Point::new(x, y))
            .collect::<Vec<_>>();
        Ok(Self {
            name: name.into(),
            design_points,
            spec,
            options,
            history: Vec::new(),
        })
    }

    pub fn run(mut self) -> Result<QdesResult> {
        let input_coords = self
            .design_points
            .iter()
            .map(|point| (point.x, point.y))
            .collect::<Vec<_>>();
        let basis = build_basis_functions(&self.spec, self.options.basis_count_per_side);

        let mut last_paneled = Vec::new();
        let mut last_result: Option<XfoilViscousResult> = None;
        let mut last_stats: Option<ErrorStats> = None;

        for iteration in 0..=self.options.outer_iterations {
            let paneled_points = repanel_points(
                &self.design_points,
                self.options.panel_count,
                &self.options.paneling,
            )?;
            let paneled_coords = paneled_points
                .iter()
                .map(|point| (point.x, point.y))
                .collect::<Vec<_>>();

            let (mut state, factorized) = build_state_from_coords(
                &self.name,
                &paneled_coords,
                self.spec.operating_point,
                &self.options.xfoil_options,
            )?;
            let oper_result =
                solve_operating_point_from_state(&mut state, &factorized, &self.options.xfoil_options)?;
            let stats = error_stats_from_state(&self.spec, &state);
            let geometry_delta_norm = if let Some(previous) = self.history.last() {
                previous.geometry_delta_norm
            } else {
                0.0
            };
            self.history.push(QdesIterationSnapshot {
                iteration,
                rms_error: stats.rms_error,
                max_error: stats.max_error,
                geometry_delta_norm,
                oper_result: oper_result.clone(),
            });

            last_paneled = paneled_coords;
            last_result = Some(oper_result.clone());
            last_stats = Some(stats.clone());

            if stats.max_error <= self.options.target_tolerance || iteration == self.options.outer_iterations {
                return Ok(QdesResult {
                    input_coords,
                    output_coords: self
                        .design_points
                        .iter()
                        .map(|point| (point.x, point.y))
                        .collect(),
                    paneled_coords: last_paneled,
                    oper_result,
                    converged: stats.max_error <= self.options.target_tolerance,
                    iterations: iteration,
                    rms_error: stats.rms_error,
                    max_error: stats.max_error,
                    target_kind: self.spec.target_kind,
                    target_upper: self.spec.upper.clone(),
                    target_lower: self.spec.lower.clone(),
                    achieved_upper: stats.achieved_upper,
                    achieved_lower: stats.achieved_lower,
                    history: self.history,
                });
            }

            if basis.is_empty() {
                break;
            }

            let coefficients = self.solve_basis_step(&basis, &stats.error_vector)?;
            let geometry_delta_norm = apply_basis_update(
                &mut self.design_points,
                &basis,
                &coefficients,
                self.options.max_basis_delta,
            );
            if let Some(snapshot) = self.history.last_mut() {
                snapshot.geometry_delta_norm = geometry_delta_norm;
            }
            if geometry_delta_norm <= 1.0e-10 {
                break;
            }
        }

        let final_result = last_result.ok_or_else(|| {
            XfoilError::Message("QDES did not produce an operating-point result".to_string())
        })?;
        let final_stats = last_stats.ok_or_else(|| {
            XfoilError::Message("QDES did not produce target-error diagnostics".to_string())
        })?;
        Ok(QdesResult {
            input_coords,
            output_coords: self
                .design_points
                .iter()
                .map(|point| (point.x, point.y))
                .collect(),
            paneled_coords: last_paneled,
            oper_result: final_result,
            converged: final_stats.max_error <= self.options.target_tolerance,
            iterations: self.history.len().saturating_sub(1),
            rms_error: final_stats.rms_error,
            max_error: final_stats.max_error,
            target_kind: self.spec.target_kind,
            target_upper: self.spec.upper.clone(),
            target_lower: self.spec.lower.clone(),
            achieved_upper: final_stats.achieved_upper,
            achieved_lower: final_stats.achieved_lower,
            history: self.history,
        })
    }

    fn solve_basis_step(&self, basis: &[BasisFunction], current_error: &[f64]) -> Result<Vec<f64>> {
        let m = current_error.len();
        let n = basis.len();
        if m == 0 || n == 0 {
            return Ok(Vec::new());
        }

        let mut sensitivity = DMatrix::<f64>::zeros(m, n);
        for (col, basis_function) in basis.iter().enumerate() {
            let mut perturbed = self.design_points.clone();
            let step = self.options.finite_difference_step.max(1.0e-5);
            apply_basis_update(
                &mut perturbed,
                std::slice::from_ref(basis_function),
                &[step],
                step,
            );
            let paneled_points =
                repanel_points(&perturbed, self.options.panel_count, &self.options.paneling)?;
            let paneled_coords = paneled_points
                .iter()
                .map(|point| (point.x, point.y))
                .collect::<Vec<_>>();
            let stats = match build_state_from_coords(
                &self.name,
                &paneled_coords,
                self.spec.operating_point,
                &self.options.xfoil_options,
            )
            .and_then(|(mut state, factorized)| {
                solve_operating_point_from_state(&mut state, &factorized, &self.options.xfoil_options)?;
                Ok(error_stats_from_state(&self.spec, &state))
            }) {
                Ok(stats) => stats,
                Err(_) => continue,
            };
            for row in 0..m.min(stats.error_vector.len()) {
                sensitivity[(row, col)] = (stats.error_vector[row] - current_error[row]) / step;
            }
        }

        let s_t = sensitivity.transpose();
        let lhs =
            &s_t * &sensitivity + DMatrix::<f64>::identity(n, n) * (1.0e-6 / self.options.update_damping.max(1.0e-3));
        let rhs = -(&s_t * DVector::from_column_slice(current_error));
        let Some(solution) = lhs.lu().solve(&rhs) else {
            return Ok(vec![0.0; n]);
        };
        Ok(solution
            .iter()
            .map(|value| (value * self.options.update_damping).clamp(
                -self.options.max_basis_delta,
                self.options.max_basis_delta,
            ))
            .collect())
    }
}

pub fn solve_body_qdes(body: &Body, spec: QdesSpec, options: QdesOptions) -> Result<QdesResult> {
    let coords = coords_from_body(body);
    solve_coords_qdes(&body.name, &coords, spec, options)
}

pub fn solve_coords_qdes(
    name: &str,
    coords: &[(f64, f64)],
    spec: QdesSpec,
    options: QdesOptions,
) -> Result<QdesResult> {
    InverseDesignSession::new(name, coords, spec, options)?.run()
}

fn validate_target(surface: &str, target: &QdesTarget) -> Result<()> {
    if target.x.len() != target.values.len() || target.x.is_empty() {
        return Err(XfoilError::Message(format!(
            "{surface} target x/value arrays must have the same non-zero length"
        )));
    }
    for pair in target.x.windows(2) {
        if pair[1] < pair[0] {
            return Err(XfoilError::Message(format!(
                "{surface} target x locations must be monotonic"
            )));
        }
    }
    Ok(())
}

fn repanel_points(points: &[Point], panel_count: usize, paneling: &PanelingParams) -> Result<Vec<Point>> {
    let spline = CubicSpline::from_points(points)
        .map_err(|err| XfoilError::Geometry(err.to_string()))?;
    Ok(spline.resample_xfoil(panel_count.max(8), paneling))
}

fn build_basis_functions(spec: &QdesSpec, basis_count_per_side: usize) -> Vec<BasisFunction> {
    if basis_count_per_side == 0 {
        return Vec::new();
    }
    let width = (0.75 / basis_count_per_side as f64).max(0.12);
    let mut basis = Vec::new();
    for surface in [SurfaceId::Lower, SurfaceId::Upper] {
        let target_present = match surface {
            SurfaceId::Upper => spec.upper.is_some(),
            SurfaceId::Lower => spec.lower.is_some(),
        };
        if !target_present {
            continue;
        }
        for idx in 0..basis_count_per_side {
            let center = (idx + 1) as f64 / (basis_count_per_side + 1) as f64;
            basis.push(BasisFunction {
                surface,
                center,
                width,
            });
        }
    }
    basis
}

fn apply_basis_update(
    points: &mut [Point],
    basis: &[BasisFunction],
    coefficients: &[f64],
    coefficient_limit: f64,
) -> f64 {
    if basis.is_empty() || coefficients.is_empty() || points.len() < 4 {
        return 0.0;
    }
    let le_idx = leading_edge_index(points);
    let x_min = points[le_idx].x;
    let x_max = points
        .iter()
        .fold(f64::NEG_INFINITY, |acc, point| acc.max(point.x));
    let chord = (x_max - x_min).abs().max(1.0e-6);
    let mut delta_norm_sq = 0.0;

    // Points are in XFOIL order: upper TE → LE → lower TE.
    // Indices 1..le_idx are the UPPER surface; le_idx+1..n-1 are LOWER.
    for idx in 1..le_idx {
        let eta = ((points[idx].x - x_min) / chord).clamp(0.0, 1.0);
        let delta = basis_delta(eta, SurfaceId::Upper, basis, coefficients, coefficient_limit);
        points[idx].y += delta;
        delta_norm_sq += delta * delta;
    }
    for idx in (le_idx + 1)..points.len().saturating_sub(1) {
        let eta = ((points[idx].x - x_min) / chord).clamp(0.0, 1.0);
        let delta = basis_delta(eta, SurfaceId::Lower, basis, coefficients, coefficient_limit);
        points[idx].y += delta;
        delta_norm_sq += delta * delta;
    }
    delta_norm_sq.sqrt()
}

fn basis_delta(
    eta: f64,
    surface: SurfaceId,
    basis: &[BasisFunction],
    coefficients: &[f64],
    coefficient_limit: f64,
) -> f64 {
    let envelope = (PI * eta).sin().powi(2);
    if envelope <= 1.0e-12 {
        return 0.0;
    }
    basis
        .iter()
        .zip(coefficients.iter().copied().chain(std::iter::repeat(0.0)))
        .filter(|(basis_function, _)| basis_function.surface == surface)
        .map(|(basis_function, coefficient)| {
            let local = ((eta - basis_function.center) / basis_function.width).powi(2);
            coefficient.clamp(-coefficient_limit, coefficient_limit) * envelope * (-0.5 * local).exp()
        })
        .sum()
}

fn leading_edge_index(points: &[Point]) -> usize {
    points
        .iter()
        .enumerate()
        .min_by(|(_, lhs), (_, rhs)| lhs.x.partial_cmp(&rhs.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(idx, _)| idx)
        .unwrap_or(0)
}

fn error_stats_from_state(spec: &QdesSpec, state: &crate::state::XfoilState) -> ErrorStats {
    let upper_actual = surface_distribution(&state.upper_rows, state.iblte_upper, spec.target_kind);
    let lower_actual = surface_distribution(&state.lower_rows, state.iblte_lower, spec.target_kind);

    let mut error_vector = Vec::new();
    let achieved_upper = sampled_distribution(spec.upper.as_ref(), &upper_actual, &mut error_vector);
    let achieved_lower = sampled_distribution(spec.lower.as_ref(), &lower_actual, &mut error_vector);

    let rms_error = if error_vector.is_empty() {
        0.0
    } else {
        (error_vector.iter().map(|value| value * value).sum::<f64>() / error_vector.len() as f64).sqrt()
    };
    let max_error = error_vector
        .iter()
        .fold(0.0_f64, |acc, value| acc.max(value.abs()));

    ErrorStats {
        max_error,
        rms_error,
        error_vector,
        achieved_upper,
        achieved_lower,
    }
}

fn surface_distribution(
    rows: &[XfoilBlRow],
    iblte: usize,
    target_kind: QdesTargetKind,
) -> SurfaceDistribution {
    if rows.is_empty() {
        return SurfaceDistribution {
            x: Vec::new(),
            values: Vec::new(),
        };
    }
    let end = iblte.min(rows.len().saturating_sub(1));
    let mut samples = rows[..=end]
        .iter()
        .map(|row| {
            let value = match target_kind {
                QdesTargetKind::EdgeVelocity => row.uedg,
                QdesTargetKind::PressureCoefficient => 1.0 - row.uedg * row.uedg,
            };
            (row.x_coord, value)
        })
        .collect::<Vec<_>>();
    samples.sort_by(|lhs, rhs| lhs.0.partial_cmp(&rhs.0).unwrap_or(std::cmp::Ordering::Equal));
    SurfaceDistribution {
        x: samples.iter().map(|(x, _)| *x).collect(),
        values: samples.iter().map(|(_, value)| *value).collect(),
    }
}

fn sampled_distribution(
    target: Option<&QdesTarget>,
    actual: &SurfaceDistribution,
    error_vector: &mut Vec<f64>,
) -> Option<SurfaceDistribution> {
    let target = target?;
    let mut sampled = Vec::with_capacity(target.x.len());
    for (x, target_value) in target.x.iter().copied().zip(target.values.iter().copied()) {
        let achieved = interpolate_clamped(&actual.x, &actual.values, x);
        sampled.push(achieved);
        error_vector.push(achieved - target_value);
    }
    Some(SurfaceDistribution {
        x: target.x.clone(),
        values: sampled,
    })
}

fn interpolate_clamped(x: &[f64], y: &[f64], target: f64) -> f64 {
    if x.is_empty() || y.is_empty() {
        return 0.0;
    }
    if target <= x[0] {
        return y[0];
    }
    for idx in 1..x.len().min(y.len()) {
        if target <= x[idx] {
            let dx = x[idx] - x[idx - 1];
            if dx.abs() <= 1.0e-12 {
                return y[idx];
            }
            let weight = (target - x[idx - 1]) / dx;
            return y[idx - 1] + weight * (y[idx] - y[idx - 1]);
        }
    }
    *y.last().unwrap_or(&0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oper::AlphaSpec;
    use crate::result::{QdesSpec, QdesTarget, QdesTargetKind};

    fn naca0012_coords() -> Vec<(f64, f64)> {
        let n = 65;
        let mut pts = Vec::new();
        // Upper TE -> LE
        for i in (0..n).rev() {
            let beta = std::f64::consts::PI * i as f64 / (n - 1) as f64;
            let x = 0.5 * (1.0 - beta.cos());
            let yt = 0.12 / 0.2 * (
                0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2)
                + 0.2843 * x.powi(3) - 0.1015 * x.powi(4)
            );
            pts.push((x, yt));
        }
        // Lower LE+1 -> TE
        for i in 1..n {
            let beta = std::f64::consts::PI * i as f64 / (n - 1) as f64;
            let x = 0.5 * (1.0 - beta.cos());
            let yt = 0.12 / 0.2 * (
                0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2)
                + 0.2843 * x.powi(3) - 0.1015 * x.powi(4)
            );
            pts.push((x, -yt));
        }
        pts
    }

    #[test]
    fn qdes_validate_requires_target() {
        let spec = QdesSpec {
            operating_point: AlphaSpec::AlphaDeg(0.0),
            target_kind: QdesTargetKind::PressureCoefficient,
            upper: None,
            lower: None,
        };
        assert!(spec.validate().is_err());
    }

    #[test]
    fn qdes_validate_rejects_mismatched_lengths() {
        let spec = QdesSpec {
            operating_point: AlphaSpec::AlphaDeg(0.0),
            target_kind: QdesTargetKind::PressureCoefficient,
            upper: Some(QdesTarget {
                x: vec![0.1, 0.5],
                values: vec![0.5], // wrong length
            }),
            lower: None,
        };
        assert!(spec.validate().is_err());
    }

    #[test]
    fn qdes_validate_accepts_valid_spec() {
        let spec = QdesSpec {
            operating_point: AlphaSpec::AlphaDeg(0.0),
            target_kind: QdesTargetKind::PressureCoefficient,
            upper: Some(QdesTarget {
                x: vec![0.1, 0.3, 0.5, 0.7],
                values: vec![-0.5, -0.3, -0.1, 0.1],
            }),
            lower: None,
        };
        assert!(spec.validate().is_ok());
    }

    #[test]
    fn interpolate_clamped_at_endpoints() {
        let x = vec![0.0, 0.5, 1.0];
        let y = vec![1.0, 2.0, 3.0];

        assert!((interpolate_clamped(&x, &y, -0.5) - 1.0).abs() < 1e-10);
        assert!((interpolate_clamped(&x, &y, 1.5) - 3.0).abs() < 1e-10);
        assert!((interpolate_clamped(&x, &y, 0.25) - 1.5).abs() < 1e-10);
    }

    #[test]
    fn leading_edge_index_finds_minimum_x() {
        let pts = vec![
            Point::new(1.0, 0.05),
            Point::new(0.5, 0.06),
            Point::new(0.0, 0.0),
            Point::new(0.5, -0.06),
            Point::new(1.0, -0.05),
        ];
        assert_eq!(leading_edge_index(&pts), 2);
    }

    #[test]
    fn basis_functions_cover_targeted_surfaces() {
        let spec = QdesSpec {
            operating_point: AlphaSpec::AlphaDeg(0.0),
            target_kind: QdesTargetKind::PressureCoefficient,
            upper: Some(QdesTarget {
                x: vec![0.1, 0.5, 0.9],
                values: vec![-0.5, -0.2, 0.1],
            }),
            lower: None,
        };
        let basis = build_basis_functions(&spec, 4);
        assert_eq!(basis.len(), 4); // 4 per targeted surface
        for b in &basis {
            assert!(matches!(b.surface, SurfaceId::Upper));
        }
    }

    #[test]
    fn solve_coords_qdes_runs_on_naca0012() {
        let coords = naca0012_coords();
        let spec = QdesSpec {
            operating_point: AlphaSpec::AlphaDeg(0.0),
            target_kind: QdesTargetKind::PressureCoefficient,
            upper: Some(QdesTarget {
                x: vec![0.1, 0.3, 0.5, 0.7, 0.9],
                values: vec![-0.4, -0.2, -0.1, 0.0, 0.1],
            }),
            lower: None,
        };
        let options = QdesOptions {
            outer_iterations: 2,
            ..Default::default()
        };

        let result = solve_coords_qdes("NACA0012", &coords, spec, options);
        assert!(result.is_ok(), "QDES should not error: {:?}", result.err());

        let qr = result.unwrap();
        assert!(!qr.output_coords.is_empty());
        assert!(qr.iterations <= 3);
        assert!(qr.rms_error.is_finite());
    }

    #[test]
    fn session_rejects_too_few_coords() {
        let spec = QdesSpec {
            operating_point: AlphaSpec::AlphaDeg(0.0),
            target_kind: QdesTargetKind::PressureCoefficient,
            upper: Some(QdesTarget { x: vec![0.5], values: vec![0.0] }),
            lower: None,
        };
        let result = InverseDesignSession::new("test", &[(0.0, 0.0), (1.0, 0.0)], spec, QdesOptions::default());
        assert!(result.is_err());
    }
}
