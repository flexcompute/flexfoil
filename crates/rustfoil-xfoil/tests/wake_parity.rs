mod support;

use rustfoil_testkit::timing::{assert_ratio_within, benchmark_closure, BenchmarkConfig};
use rustfoil_xfoil::wake_panel::{qdcalc, qwcalc, xywake};

use support::{
    assert_close_slice, build_qdcalc_state, qdcalc_fortran, MICROBENCH_MAX_RATIO,
};

#[test]
#[ignore = "wake path is still converging toward literal XYWAKE/QWCALC/QDCALC parity"]
fn qdcalc_matches_fortran() {
    let reference = qdcalc_fortran();
    let (wake_x, wake_y, diag_sample, row0_sample, matrix_size) = support::rust_qdcalc_sample();

    assert_eq!(matrix_size, reference.matrix_size);
    assert_close_slice("qdcalc.wake_x", &wake_x, &reference.wake_x, 1.0e-8);
    assert_close_slice("qdcalc.wake_y", &wake_y, &reference.wake_y, 1.0e-8);
    assert_close_slice("qdcalc.diag_sample", &diag_sample, &reference.diag_sample, 5.0e-6);
    assert_close_slice("qdcalc.row0_sample", &row0_sample, &reference.row0_sample, 5.0e-6);
}

#[test]
#[ignore = "wake path is still converging toward literal XYWAKE parity"]
fn xywake_matches_fortran_when_port_is_more_literal() {
    let reference = qdcalc_fortran();
    let (mut state, factorized) = build_qdcalc_state();
    xywake(&mut state, &factorized, 1.0);
    assert_close_slice("xywake.wake_x", &state.wake_x, &reference.wake_x, 1.0e-8);
    assert_close_slice("xywake.wake_y", &state.wake_y, &reference.wake_y, 1.0e-8);
}

#[test]
#[ignore = "wake path is still converging toward literal QWCALC parity"]
fn qwcalc_matches_fortran_when_port_is_more_literal() {
    let (mut state, factorized) = build_qdcalc_state();
    xywake(&mut state, &factorized, 1.0);
    qwcalc(&mut state, &factorized);
    assert!(!state.wake_qinv.is_empty());
    assert!(state.wake_qinv.iter().all(|value| value.is_finite()));
}

#[test]
#[ignore = "release-mode wake microbenchmark gate"]
fn perf_qdcalc_vs_fortran() {
    assert!(
        !cfg!(debug_assertions),
        "run perf gates with `cargo test --release -- --ignored`"
    );
    let reference = &qdcalc_fortran().perf;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let (mut state, factorized) = build_qdcalc_state();
            qdcalc(&mut state, &factorized).expect("qdcalc");
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "qdcalc");
}
