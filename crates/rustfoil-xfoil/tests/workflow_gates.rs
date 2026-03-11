mod support;

use rustfoil_testkit::timing::{assert_ratio_within, benchmark_closure, BenchmarkConfig};

use support::{assert_close_scalar, run_workflow_case, workflow_fortran, WORKFLOW_MAX_RATIO};

#[test]
#[ignore = "faithful workflow parity gate; enable as the VISCAL port tightens"]
fn viscal_style_oper_point_matches_fortran() {
    let alpha_deg = 4.0;
    let rust = run_workflow_case(alpha_deg);
    let reference = workflow_fortran(alpha_deg);
    let case = reference
        .cases
        .iter()
        .find(|case| (case.alpha_deg - alpha_deg).abs() < 1.0e-9)
        .expect("Fortran case for requested alpha");

    assert_close_scalar("cl", rust.cl, case.cl, 5.0e-2);
    assert_close_scalar("cd", rust.cd, case.cd, 5.0e-3);
    assert_close_scalar("cm", rust.cm, case.cm, 5.0e-2);
    assert_eq!(rust.converged, case.converged);
}

#[test]
#[ignore = "release-mode end-to-end workflow gate"]
fn perf_viscal_style_oper_point_vs_fortran() {
    assert!(
        !cfg!(debug_assertions),
        "run perf gates with `cargo test --release -- --ignored`"
    );
    let alpha_deg = 4.0;
    let rust_stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: 10,
            warmup_runs: 1,
            sample_runs: 5,
        },
        || {
            let _ = run_workflow_case(alpha_deg);
        },
    );
    let fortran_stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: 10,
            warmup_runs: 1,
            sample_runs: 5,
        },
        || {
            let _ = workflow_fortran(alpha_deg);
        },
    );
    assert_ratio_within(
        &rust_stats,
        fortran_stats.median_seconds,
        WORKFLOW_MAX_RATIO,
        "viscal_oper_point",
    );
}
