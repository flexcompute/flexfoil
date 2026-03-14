mod support;

use rustfoil_xfoil::{
    solve_coords_qdes, AlphaSpec, QdesOptions, QdesSpec, QdesTarget, QdesTargetKind, SurfaceDistribution,
    XfoilOptions,
};

use support::{qdes_fortran, xfoil_binary_qdes_distributions, FortranSurfaceDistribution};

fn pairs(xs: &[f64], ys: &[f64]) -> Vec<(f64, f64)> {
    xs.iter().copied().zip(ys.iter().copied()).collect()
}

fn upper_surface(coords: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let le_idx = coords
        .iter()
        .enumerate()
        .min_by(|(_, lhs), (_, rhs)| lhs.0.partial_cmp(&rhs.0).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);
    let mut upper = coords[..=le_idx].to_vec();
    upper.sort_by(|lhs, rhs| lhs.0.partial_cmp(&rhs.0).unwrap());
    upper
}

fn max_abs_pair_delta(lhs: &[(f64, f64)], rhs: &[(f64, f64)]) -> f64 {
    lhs.iter()
        .zip(rhs.iter())
        .map(|(lhs_pt, rhs_pt)| (lhs_pt.0 - rhs_pt.0).abs().max((lhs_pt.1 - rhs_pt.1).abs()))
        .fold(0.0_f64, f64::max)
}

fn max_abs_y_delta(lhs: &[(f64, f64)], rhs: &[(f64, f64)]) -> f64 {
    lhs.iter()
        .zip(rhs.iter())
        .map(|(lhs_pt, rhs_pt)| (lhs_pt.1 - rhs_pt.1).abs())
        .fold(0.0_f64, f64::max)
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

fn max_sample_delta(
    actual: &SurfaceDistribution,
    reference: &FortranSurfaceDistribution,
    min_x: f64,
    max_x: f64,
) -> f64 {
    actual
        .x
        .iter()
        .copied()
        .zip(actual.values.iter().copied())
        .filter(|(x, _)| *x >= min_x && *x <= max_x)
        .map(|(x, value)| (value - interpolate_clamped(&reference.x, &reference.values, x)).abs())
        .fold(0.0_f64, f64::max)
}

fn sample_pairs(
    actual: &SurfaceDistribution,
    reference: &FortranSurfaceDistribution,
    min_x: f64,
    max_x: f64,
    count: usize,
) -> Vec<(f64, f64, f64)> {
    actual
        .x
        .iter()
        .copied()
        .zip(actual.values.iter().copied())
        .filter(|(x, _)| *x >= min_x && *x <= max_x)
        .take(count)
        .map(|(x, value)| (x, value, interpolate_clamped(&reference.x, &reference.values, x)))
        .collect()
}

#[test]
fn qdes_noop_matches_xfoil_driver_geometry() {
    let reference = qdes_fortran(0.0, 2);
    let base_coords = pairs(&reference.base_x, &reference.base_y);
    let spec = QdesSpec {
        operating_point: AlphaSpec::AlphaDeg(reference.alpha_deg),
        target_kind: QdesTargetKind::EdgeVelocity,
        upper: Some(QdesTarget {
            x: reference.target_upper_x.clone(),
            values: reference.target_upper_q.clone(),
        }),
        lower: Some(QdesTarget {
            x: reference.target_lower_x.clone(),
            values: reference.target_lower_q.clone(),
        }),
    };
    let result = solve_coords_qdes(
        "qdes-noop",
        &base_coords,
        spec,
        QdesOptions {
            xfoil_options: XfoilOptions {
                reynolds: 1.0e6,
                mach: 0.0,
                ncrit: 9.0,
                max_iterations: 4,
                ..Default::default()
            },
            outer_iterations: reference.iterations,
            basis_count_per_side: 4,
            panel_count: 80,
            target_tolerance: 1.0e-6,
            ..Default::default()
        },
    )
    .expect("run rust qdes noop parity case");

    let xfoil_coords = pairs(&reference.output_x, &reference.output_y);
    let xfoil_baseline_delta = max_abs_pair_delta(&base_coords, &xfoil_coords);
    assert!(
        xfoil_baseline_delta < 1.0e-6,
        "XFOIL noop reference drifted from the input geometry: {xfoil_baseline_delta}"
    );
    let max_delta = max_abs_pair_delta(&result.output_coords, &xfoil_coords);
    assert!(
        max_delta < 2.0e-2,
        "Rust/XFOIL noop geometry diverged too much: {max_delta}"
    );

    let (xfoil_upper, xfoil_lower) =
        xfoil_binary_qdes_distributions("noop", &xfoil_coords, reference.alpha_deg);
    let min_compare_x = 5.0e-2;
    let max_compare_x = 9.5e-1;
    let upper_delta = max_sample_delta(
        result.achieved_upper.as_ref().expect("Rust noop achieved upper distribution"),
        &xfoil_upper,
        min_compare_x,
        max_compare_x,
    );
    let lower_delta = max_sample_delta(
        result.achieved_lower.as_ref().expect("Rust noop achieved lower distribution"),
        &xfoil_lower,
        min_compare_x,
        max_compare_x,
    );
    assert!(
        upper_delta < 1.5e-1,
        "Rust/XFOIL noop upper distribution mismatch too large: {upper_delta}; samples={:?}",
        sample_pairs(
            result.achieved_upper.as_ref().expect("Rust noop achieved upper distribution"),
            &xfoil_upper,
            min_compare_x,
            max_compare_x,
            6,
        )
    );
    assert!(
        lower_delta < 1.5e-1,
        "Rust/XFOIL noop lower distribution mismatch too large: {lower_delta}"
    );
}

#[test]
fn qdes_upper_surface_case_tracks_xfoil_driver_shape() {
    let reference = qdes_fortran(0.02, 2);
    let base_coords = pairs(&reference.base_x, &reference.base_y);
    let spec = QdesSpec {
        operating_point: AlphaSpec::AlphaDeg(reference.alpha_deg),
        target_kind: QdesTargetKind::EdgeVelocity,
        upper: Some(QdesTarget {
            x: reference.target_upper_x.clone(),
            values: reference.target_upper_q.clone(),
        }),
        lower: Some(QdesTarget {
            x: reference.target_lower_x.clone(),
            values: reference.target_lower_q.clone(),
        }),
    };
    let result = solve_coords_qdes(
        "qdes-parity",
        &base_coords,
        spec,
        QdesOptions {
            xfoil_options: XfoilOptions {
                reynolds: 1.0e6,
                mach: 0.0,
                ncrit: 9.0,
                max_iterations: 4,
                ..Default::default()
            },
            outer_iterations: 4,
            basis_count_per_side: 6,
            panel_count: 80,
            update_damping: 0.5,
            max_basis_delta: 1.0e-2,
            target_tolerance: 2.0e-3,
            ..Default::default()
        },
    )
    .expect("run rust qdes parity case");

    let rust_upper = upper_surface(&result.output_coords);
    let xfoil_upper = upper_surface(&pairs(&reference.output_x, &reference.output_y));
    assert_eq!(rust_upper.len(), xfoil_upper.len(), "upper surface sample count mismatch");
    let max_y_delta = max_abs_y_delta(&rust_upper, &xfoil_upper);
    assert!(
        max_y_delta < 2.5e-2,
        "Rust/XFOIL upper-surface shape mismatch too large: {max_y_delta}"
    );

    let (xfoil_upper_dist, xfoil_lower_dist) =
        xfoil_binary_qdes_distributions("perturbed", &pairs(&reference.output_x, &reference.output_y), reference.alpha_deg);
    let min_compare_x = 5.0e-2;
    let max_compare_x = 9.5e-1;
    let upper_q_delta = max_sample_delta(
        result.achieved_upper.as_ref().expect("Rust achieved upper distribution"),
        &xfoil_upper_dist,
        min_compare_x,
        max_compare_x,
    );
    let lower_q_delta = max_sample_delta(
        result.achieved_lower.as_ref().expect("Rust achieved lower distribution"),
        &xfoil_lower_dist,
        min_compare_x,
        max_compare_x,
    );
    assert!(
        upper_q_delta < 1.5e-1,
        "Rust/XFOIL upper speed-distribution mismatch too large: {upper_q_delta}; samples={:?}",
        sample_pairs(
            result.achieved_upper.as_ref().expect("Rust achieved upper distribution"),
            &xfoil_upper_dist,
            min_compare_x,
            max_compare_x,
            6,
        )
    );
    assert!(
        lower_q_delta < 3.0e-1,
        "Rust/XFOIL lower speed-distribution mismatch too large: {lower_q_delta}"
    );
}
