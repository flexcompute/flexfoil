mod support;

use rustfoil_xfoil::{
    build_state_from_coords, solve_coords_qdes, solve_operating_point_from_state, AlphaSpec,
    QdesOptions, QdesSpec, QdesTarget, QdesTargetKind, SurfaceDistribution, XfoilOptions,
};

use rustfoil_core::{CubicSpline, PanelingParams, Point};

use support::naca0012_coords;

fn cp_distribution_from_coords(coords: &[(f64, f64)], alpha_deg: f64) -> (SurfaceDistribution, SurfaceDistribution) {
    let options = XfoilOptions {
        reynolds: 1.0e6,
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 4,
        ..Default::default()
    };
    let (mut state, factorized) =
        build_state_from_coords("test-airfoil", coords, AlphaSpec::AlphaDeg(alpha_deg), &options)
            .expect("build faithful state");
    solve_operating_point_from_state(&mut state, &factorized, &options).expect("solve faithful state");
    (
        cp_distribution(&state.upper_rows, state.iblte_upper),
        cp_distribution(&state.lower_rows, state.iblte_lower),
    )
}

fn cp_distribution(
    rows: &[rustfoil_xfoil::XfoilBlRow],
    iblte: usize,
) -> SurfaceDistribution {
    let end = iblte.min(rows.len().saturating_sub(1));
    let mut samples = rows[..=end]
        .iter()
        .map(|row| (row.x_coord, 1.0 - row.uedg * row.uedg))
        .collect::<Vec<_>>();
    samples.sort_by(|lhs, rhs| lhs.0.partial_cmp(&rhs.0).unwrap());
    SurfaceDistribution {
        x: samples.iter().map(|(x, _)| *x).collect(),
        values: samples.iter().map(|(_, value)| *value).collect(),
    }
}

fn build_cp_spec(coords: &[(f64, f64)], alpha_deg: f64) -> QdesSpec {
    let (upper, lower) = cp_distribution_from_coords(coords, alpha_deg);
    QdesSpec {
        operating_point: AlphaSpec::AlphaDeg(alpha_deg),
        target_kind: QdesTargetKind::PressureCoefficient,
        upper: Some(QdesTarget {
            x: upper.x,
            values: upper.values,
        }),
        lower: Some(QdesTarget {
            x: lower.x,
            values: lower.values,
        }),
    }
}

fn camber_shift(coords: &[(f64, f64)], magnitude: f64) -> Vec<(f64, f64)> {
    let mut shifted = coords.to_vec();
    let interior_len = shifted.len().saturating_sub(2);
    for point in shifted.iter_mut().skip(1).take(interior_len) {
        let x = point.0.clamp(0.0, 1.0);
        point.1 += magnitude * x * (1.0 - x);
    }
    shifted
}

fn max_target_error(result: &rustfoil_xfoil::QdesResult) -> f64 {
    result.max_error
}

fn repanel_coords(coords: &[(f64, f64)], panels: usize) -> Vec<(f64, f64)> {
    let points = coords
        .iter()
        .map(|&(x, y)| Point::new(x, y))
        .collect::<Vec<_>>();
    CubicSpline::from_points(&points)
        .expect("build spline for test repanel")
        .resample_xfoil(panels, &PanelingParams::default())
        .into_iter()
        .map(|point| (point.x, point.y))
        .collect()
}

#[test]
fn qdes_matches_self_target_without_shape_change() {
    let coords = naca0012_coords(80);
    let seed_spec = build_cp_spec(&repanel_coords(&coords, 80), 4.0);
    let seed_options = QdesOptions {
        xfoil_options: XfoilOptions {
            reynolds: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 4,
            ..Default::default()
        },
        outer_iterations: 0,
        basis_count_per_side: 2,
        panel_count: 80,
        ..Default::default()
    };
    let seed_result =
        solve_coords_qdes("naca0012-seed", &coords, seed_spec, seed_options).expect("seed qdes result");
    let spec = QdesSpec {
        operating_point: AlphaSpec::AlphaDeg(4.0),
        target_kind: QdesTargetKind::PressureCoefficient,
        upper: seed_result.achieved_upper.clone().map(|surface| QdesTarget {
            x: surface.x,
            values: surface.values,
        }),
        lower: seed_result.achieved_lower.clone().map(|surface| QdesTarget {
            x: surface.x,
            values: surface.values,
        }),
    };
    let options = QdesOptions {
        xfoil_options: XfoilOptions {
            reynolds: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 4,
            ..Default::default()
        },
        outer_iterations: 2,
        basis_count_per_side: 2,
        panel_count: 80,
        target_tolerance: 1.0e-8,
        ..Default::default()
    };

    let result = solve_coords_qdes("naca0012", &coords, spec, options).expect("run qdes self target");

    assert!(result.max_error < 1.0e-8, "self target max error = {}", result.max_error);
    assert!(result.rms_error < 1.0e-8, "self target rms error = {}", result.rms_error);
    assert!(result.converged);
}

#[test]
fn qdes_reduces_target_error_for_nearby_geometry() {
    let base_coords = naca0012_coords(80);
    let target_coords = camber_shift(&base_coords, 0.01);
    let spec = build_cp_spec(&target_coords, 4.0);

    let baseline = solve_coords_qdes(
        "baseline",
        &base_coords,
        spec.clone(),
        QdesOptions {
            xfoil_options: XfoilOptions {
                reynolds: 1.0e6,
                mach: 0.0,
                ncrit: 9.0,
                max_iterations: 4,
                ..Default::default()
            },
            outer_iterations: 0,
            basis_count_per_side: 3,
            panel_count: 80,
            ..Default::default()
        },
    )
    .expect("baseline qdes evaluation");

    let improved = solve_coords_qdes(
        "improved",
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
            outer_iterations: 2,
            basis_count_per_side: 3,
            panel_count: 80,
            update_damping: 0.45,
            max_basis_delta: 8.0e-3,
            ..Default::default()
        },
    )
    .expect("improved qdes evaluation");

    assert!(
        max_target_error(&improved) < max_target_error(&baseline),
        "expected improved max error < baseline max error, got {} vs {}",
        improved.max_error,
        baseline.max_error
    );
    assert_ne!(improved.output_coords, improved.input_coords);
}
