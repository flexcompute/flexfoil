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
    assert_close_slice("mrchue.upper_x", &actual.upper_x, &reference.upper_x, 1.0e-6);
    assert_close_slice("mrchue.lower_x", &actual.lower_x, &reference.lower_x, 1.0e-6);
    assert_close_slice("mrchue.upper_theta", &actual.upper_theta, &reference.upper_theta, 5.0e-6);
    assert_close_slice("mrchue.lower_theta", &actual.lower_theta, &reference.lower_theta, 5.0e-6);
    assert_close_slice("mrchue.upper_dstr", &actual.upper_dstr, &reference.upper_dstr, 5.0e-6);
    assert_close_slice("mrchue.lower_dstr", &actual.lower_dstr, &reference.lower_dstr, 5.0e-6);
    assert_close_slice("mrchue.upper_uedg", &actual.upper_uedg, &reference.upper_uedg, 5.0e-6);
    assert_close_slice("mrchue.lower_uedg", &actual.lower_uedg, &reference.lower_uedg, 5.0e-6);
    assert_close_slice("mrchue.upper_ctau", &actual.upper_ctau, &reference.upper_ctau, 5.0e-6);
    assert_close_slice("mrchue.lower_ctau", &actual.lower_ctau, &reference.lower_ctau, 5.0e-6);
    assert_close_slice("mrchue.upper_mass", &actual.upper_mass, &reference.upper_mass, 5.0e-6);
    assert_close_slice("mrchue.lower_mass", &actual.lower_mass, &reference.lower_mass, 5.0e-6);
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
