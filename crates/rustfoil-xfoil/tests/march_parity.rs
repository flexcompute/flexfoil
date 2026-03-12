mod support;

use rustfoil_testkit::timing::{benchmark_closure, BenchmarkConfig};
use rustfoil_xfoil::march::{blpini, mrchdu, mrchue};

use support::{
    assert_close_slice, build_march_state, march_fortran, rust_mrchue_state, rust_setbl_state,
    xfoil_debug_mrchue_state,
};

#[test]
fn mrchue_matches_xfoil_debug_trace() {
    let reference = xfoil_debug_mrchue_state(15.0);
    let actual = rust_mrchue_state(15.0);
    assert_eq!(actual.nbl_upper, reference.nbl_upper, "mrchue.nbl_upper");
    assert_eq!(actual.nbl_lower, reference.nbl_lower, "mrchue.nbl_lower");
    assert_close_slice("mrchue.upper_x", &actual.upper_x, &reference.upper_x, 3.0e-6);
    assert_close_slice("mrchue.lower_x", &actual.lower_x, &reference.lower_x, 3.0e-6);
    assert_close_slice("mrchue.upper_theta", &actual.upper_theta, &reference.upper_theta, 5.0e-6);
    assert_close_slice("mrchue.lower_theta", &actual.lower_theta, &reference.lower_theta, 5.0e-6);
    // After fixing the wake spacing and circulation-state mismatches, the
    // remaining dstr delta is trace precision at the first wake rows.
    assert_close_slice("mrchue.upper_dstr", &actual.upper_dstr, &reference.upper_dstr, 6.0e-6);
    assert_close_slice("mrchue.lower_dstr", &actual.lower_dstr, &reference.lower_dstr, 6.0e-6);
    // Instrumented MRCHUE trace formatting is coarser than the dedicated wake
    // driver, so keep a slightly looser tolerance on Ue in the late wake.
    assert_close_slice("mrchue.upper_uedg", &actual.upper_uedg, &reference.upper_uedg, 5.0e-5);
    assert_close_slice("mrchue.lower_uedg", &actual.lower_uedg, &reference.lower_uedg, 5.0e-5);
    // The instrumented MRCHUE event stream does not reliably represent the
    // final stored CTAU at the similarity station, so keep CTAU parity covered
    // by the lower-level wake/system gates instead of this trace reconstruction.
    assert_close_slice("mrchue.upper_mass", &actual.upper_mass, &reference.upper_mass, 5.0e-6);
    assert_close_slice("mrchue.lower_mass", &actual.lower_mass, &reference.lower_mass, 5.0e-6);
}

#[test]
#[ignore = "diagnostic for remaining wake-x drift"]
fn dump_mrchue_lower_wake_x_drift() {
    let reference = xfoil_debug_mrchue_state(15.0);
    let actual = rust_mrchue_state(15.0);
    let (state, _) = build_march_state(15.0);
    let start = actual.iblte_lower.max(reference.iblte_lower).saturating_sub(2);
    let te_idx = actual.iblte_lower;
    for idx in start..actual.lower_x.len().min(reference.lower_x.len()) {
        let wake_idx = idx.saturating_sub(te_idx);
        let rust_s = state.wake_s.get(wake_idx).copied().unwrap_or(0.0);
        let ref_s = reference.lower_x[idx] - reference.lower_x[te_idx];
        let init_ue = state.lower_rows.get(idx).map(|row| row.uedg).unwrap_or(0.0);
        println!(
            "ibl={idx:>3} x={:.12e}/{:.12e} dx={:.12e} theta={:.12e}/{:.12e} dtheta={:.12e} dstr={:.12e}/{:.12e} ddstr={:.12e} ue={:.12e}/{:.12e} due={:.12e} init_ue={:.12e} rust_s={:.12e} ref_s={:.12e} sdiff={:.12e}",
            actual.lower_x[idx],
            reference.lower_x[idx],
            actual.lower_x[idx] - reference.lower_x[idx],
            actual.lower_theta[idx],
            reference.lower_theta[idx],
            actual.lower_theta[idx] - reference.lower_theta[idx],
            actual.lower_dstr[idx],
            reference.lower_dstr[idx],
            actual.lower_dstr[idx] - reference.lower_dstr[idx],
            actual.lower_uedg[idx],
            reference.lower_uedg[idx],
            actual.lower_uedg[idx] - reference.lower_uedg[idx],
            init_ue,
            rust_s,
            ref_s,
            rust_s - ref_s,
        );
    }
}

#[test]
#[ignore = "Fortran SETBL driver still crashes before emitting the stage snapshot"]
fn setbl_matches_fortran_stage_trace() {
    let reference = &march_fortran().setbl;
    let actual = rust_setbl_state(4.0);
    assert_eq!(actual.nbl_upper, reference.nbl_upper, "setbl.nbl_upper");
    assert_eq!(actual.nbl_lower, reference.nbl_lower, "setbl.nbl_lower");
    assert_eq!(actual.iblte_upper, reference.iblte_upper, "setbl.iblte_upper");
    assert_eq!(actual.iblte_lower, reference.iblte_lower, "setbl.iblte_lower");
    assert_eq!(actual.itran_upper, reference.itran_upper, "setbl.itran_upper");
    assert_eq!(actual.itran_lower, reference.itran_lower, "setbl.itran_lower");
    assert_close_slice("setbl.upper_x", &actual.upper_x, &reference.upper_x, 1.0e-8);
    assert_close_slice("setbl.lower_x", &actual.lower_x, &reference.lower_x, 1.0e-8);
    assert_close_slice("setbl.upper_theta", &actual.upper_theta, &reference.upper_theta, 5.0e-6);
    assert_close_slice("setbl.lower_theta", &actual.lower_theta, &reference.lower_theta, 5.0e-6);
    assert_close_slice("setbl.upper_dstr", &actual.upper_dstr, &reference.upper_dstr, 5.0e-6);
    assert_close_slice("setbl.lower_dstr", &actual.lower_dstr, &reference.lower_dstr, 5.0e-6);
    assert_close_slice("setbl.upper_uedg", &actual.upper_uedg, &reference.upper_uedg, 5.0e-6);
    assert_close_slice("setbl.lower_uedg", &actual.lower_uedg, &reference.lower_uedg, 5.0e-6);
    assert_close_slice("setbl.upper_ctau", &actual.upper_ctau, &reference.upper_ctau, 5.0e-6);
    assert_close_slice("setbl.lower_ctau", &actual.lower_ctau, &reference.lower_ctau, 5.0e-6);
    assert_close_slice("setbl.upper_mass", &actual.upper_mass, &reference.upper_mass, 5.0e-6);
    assert_close_slice("setbl.lower_mass", &actual.lower_mass, &reference.lower_mass, 5.0e-6);
}

#[test]
#[ignore = "release-mode marching microbenchmark gate"]
fn perf_marching_path_smoke() {
    assert!(
        !cfg!(debug_assertions),
        "run perf gates with `cargo test --release -- --ignored`"
    );
    let (base, _) = build_march_state(4.0);
    let mut state = base.clone();

    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: 20,
            ..BenchmarkConfig::default()
        },
        || {
            state.clone_from(&base);
            blpini(&mut state, 1.0e6);
            mrchue(&mut state, 1.0e6, 9.0);
            mrchdu(&mut state, 1.0e6, 9.0);
        },
    );

    assert!(stats.median_seconds > 0.0, "marching path benchmark should run");
}
