mod support;

use rustfoil_testkit::timing::{assert_ratio_within, benchmark_closure, BenchmarkConfig};
use rustfoil_xfoil::wake_panel::{qdcalc, qwcalc, xywake};

use support::{
    assert_close_scalar, assert_close_slice, build_qdcalc_state, qdcalc_fortran, rust_qdcalc_output,
    MICROBENCH_MAX_RATIO,
};

#[test]
fn qdcalc_matches_fortran() {
    let reference = qdcalc_fortran();
    let actual = rust_qdcalc_output();

    assert_eq!(actual.matrix_size, reference.matrix_size);
    assert_close_slice("qdcalc.dij_flat", &actual.dij_flat, &reference.dij_flat, 5.0e-6);

    let n = actual.matrix_size;
    let nw = actual.wake_x.len();
    let te_row = n - nw - 1;
    let first_wake_row = te_row + 1;
    for j in 0..n {
        assert_close_scalar(
            "qdcalc.first_wake_row_equals_te",
            actual.dij_flat[n * first_wake_row + j],
            actual.dij_flat[n * te_row + j],
            5.0e-6,
        );
    }
}

#[test]
fn xywake_matches_fortran_when_port_is_more_literal() {
    let reference = qdcalc_fortran();
    let (mut state, factorized) = build_qdcalc_state();
    xywake(&mut state, &factorized, 1.0);
    assert_close_slice("xywake.wake_x", &state.wake_x, &reference.wake_x, 1.0e-4);
    assert_close_slice("xywake.wake_y", &state.wake_y, &reference.wake_y, 1.0e-4);
    assert_close_slice("xywake.wake_s", &state.wake_s, &reference.wake_s, 1.0e-10);
    assert_close_slice("xywake.wake_nx", &state.wake_nx, &reference.wake_nx, 1.0e-4);
    assert_close_slice("xywake.wake_ny", &state.wake_ny, &reference.wake_ny, 1.0e-4);
    assert_close_slice("xywake.wake_apanel", &state.wake_apanel, &reference.wake_apanel, 1.0e-4);
}

#[test]
fn qwcalc_matches_fortran_when_port_is_more_literal() {
    let reference = qdcalc_fortran();
    let (mut state, factorized) = build_qdcalc_state();
    xywake(&mut state, &factorized, 1.0);
    qwcalc(&mut state, &factorized);
    assert_close_slice(
        "qwcalc.wake_qinvu_0.downstream",
        &state.wake_qinvu_0[1..],
        &reference.wake_qinvu_0[1..],
        5.0e-6,
    );
    assert_close_slice(
        "qwcalc.wake_qinvu_90.downstream",
        &state.wake_qinvu_90[1..],
        &reference.wake_qinvu_90[1..],
        5.0e-6,
    );
    assert_close_slice(
        "qwcalc.wake_qinv.downstream",
        &state.wake_qinv[1..],
        &reference.wake_qinv[1..],
        5.0e-6,
    );
    assert_close_slice(
        "qwcalc.wake_qinv_a.downstream",
        &state.wake_qinv_a[1..],
        &reference.wake_qinv_a[1..],
        5.0e-6,
    );
}

#[test]
fn specal_recombines_wake_basis_without_rebuilding_wake() {
    let (mut state, _) = build_qdcalc_state();
    let wake_x = state.wake_x.clone();
    let wake_y = state.wake_y.clone();
    let wake_qinvu_0 = state.wake_qinvu_0.clone();
    let wake_qinvu_90 = state.wake_qinvu_90.clone();
    let alpha = 7.0_f64.to_radians();
    rustfoil_xfoil::state_ops::specal(&mut state, alpha);

    assert_eq!(state.wake_x, wake_x);
    assert_eq!(state.wake_y, wake_y);
    assert_eq!(state.wake_qinvu_0, wake_qinvu_0);
    assert_eq!(state.wake_qinvu_90, wake_qinvu_90);

    let cosa = alpha.cos();
    let sina = alpha.sin();
    let expected_qinv: Vec<f64> = wake_qinvu_0
        .iter()
        .zip(wake_qinvu_90.iter())
        .map(|(&q0, &q90)| cosa * q0 + sina * q90)
        .collect();
    let expected_qinv_a: Vec<f64> = wake_qinvu_0
        .iter()
        .zip(wake_qinvu_90.iter())
        .map(|(&q0, &q90)| -sina * q0 + cosa * q90)
        .collect();

    assert_close_slice("specal.wake_qinv", &state.wake_qinv, &expected_qinv, 1.0e-12);
    assert_close_slice("specal.wake_qinv_a", &state.wake_qinv_a, &expected_qinv_a, 1.0e-12);
}

#[test]
#[ignore = "diagnostic for wake fixture geometry"]
fn dump_qdcalc_geometry_comparison() {
    let reference = qdcalc_fortran();
    let actual = rust_qdcalc_output();
    println!(
        "xle rust={:.12e} fortran={:.12e} diff={:.12e}",
        actual.xle,
        reference.xle,
        actual.xle - reference.xle
    );
    println!(
        "yle rust={:.12e} fortran={:.12e} diff={:.12e}",
        actual.yle,
        reference.yle,
        actual.yle - reference.yle
    );
    println!(
        "xte rust={:.12e} fortran={:.12e} diff={:.12e}",
        actual.xte,
        reference.xte,
        actual.xte - reference.xte
    );
    println!(
        "yte rust={:.12e} fortran={:.12e} diff={:.12e}",
        actual.yte,
        reference.yte,
        actual.yte - reference.yte
    );
    println!(
        "chord rust={:.12e} fortran={:.12e} diff={:.12e}",
        actual.chord,
        reference.chord,
        actual.chord - reference.chord
    );
}

#[test]
#[ignore = "release-mode wake microbenchmark gate"]
fn perf_qdcalc_vs_fortran() {
    assert!(
        !cfg!(debug_assertions),
        "run perf gates with `cargo test --release -- --ignored`"
    );
    let reference = &qdcalc_fortran().perf;
    let (mut state, factorized) = build_qdcalc_state();
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            state.ladij = false;
            state.lwdij = false;
            qdcalc(&mut state, &factorized).expect("qdcalc");
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "qdcalc");
}
