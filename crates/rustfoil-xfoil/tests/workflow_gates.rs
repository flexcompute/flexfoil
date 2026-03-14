mod support;

use rustfoil_testkit::timing::{assert_ratio_within, benchmark_closure, BenchmarkConfig};
use rustfoil_testkit::fortran_runner::run_xfoil_instrumented;

use std::path::PathBuf;

use support::{
    assert_close_scalar, run_workflow_case, xfoil_binary_oper_point, xfoil_binary_reference,
    WORKFLOW_MAX_RATIO,
};

fn workflow_case_against_xfoil_binary(alpha_deg: f64) -> (rustfoil_xfoil::result::XfoilViscousResult, support::XfoilBinaryReference) {
    let rust = run_workflow_case(alpha_deg);
    let reference = xfoil_binary_oper_point(alpha_deg);
    (rust, reference)
}

#[test]
fn viscal_style_oper_point_matches_xfoil_binary() {
    let alpha_deg = 15.0;
    let rust = run_workflow_case(alpha_deg);
    let case = xfoil_binary_reference();

    assert_close_scalar("cl", rust.cl, case.cl, 7.5e-2);
    assert_close_scalar("cd", rust.cd, case.cd, 7.5e-3);
    assert_close_scalar("cm", rust.cm, case.cm, 7.5e-2);
}

#[test]
fn viscal_workflow_matrix_matches_xfoil_binary() {
    let mut failures = Vec::new();
    for alpha_deg in [-4.0, 0.0, 4.0, 8.0, 12.0, 15.0] {
        let (rust, reference) = workflow_case_against_xfoil_binary(alpha_deg);
        println!(
            "alpha={alpha_deg:>5.1} rust(cl={:>8.4}, cd={:>8.5}, cm={:>8.4}, xtr_u={:>6.3}, xtr_l={:>6.3}, conv={}, iter={}) xfoil(cl={:>8.4}, cd={:>8.5}, cm={:>8.4}, xtr_u={:>6.3}, xtr_l={:>6.3}, conv={}, iter={})",
            rust.cl,
            rust.cd,
            rust.cm,
            rust.x_tr_upper,
            rust.x_tr_lower,
            rust.converged,
            rust.iterations,
            reference.cl,
            reference.cd,
            reference.cm,
            reference.x_tr_upper,
            reference.x_tr_lower,
            reference.converged,
            reference.iterations,
        );

        if (rust.cl - reference.cl).abs() > 9.0e-2 {
            failures.push(format!("cl@{alpha_deg:.1}"));
        }
        if (rust.cd - reference.cd).abs() > 1.0e-2 {
            failures.push(format!("cd@{alpha_deg:.1}"));
        }
        if (rust.cm - reference.cm).abs() > 9.0e-2 {
            failures.push(format!("cm@{alpha_deg:.1}"));
        }
        if (rust.x_tr_upper - reference.x_tr_upper).abs() > 1.5e-1 {
            failures.push(format!("x_tr_upper@{alpha_deg:.1}"));
        }
        if (rust.x_tr_lower - reference.x_tr_lower).abs() > 1.5e-1 {
            failures.push(format!("x_tr_lower@{alpha_deg:.1}"));
        }
    }
    assert!(failures.is_empty(), "workflow matrix mismatches: {}", failures.join(", "));
}

#[test]
#[ignore = "release-mode end-to-end workflow gate"]
fn perf_viscal_style_oper_point_vs_xfoil_binary() {
    assert!(
        !cfg!(debug_assertions),
        "run perf gates with `cargo test --release -- --ignored`"
    );
    let alpha_deg = 4.0;
    let rust_stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: 3,
            warmup_runs: 1,
            sample_runs: 5,
        },
        || {
            let _ = run_workflow_case(alpha_deg);
        },
    );
    let xfoil_workdir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("target")
        .join("xfoil-binary-perf");
    std::fs::create_dir_all(&xfoil_workdir).expect("create xfoil perf dir");
    let commands =
        format!("NACA 0012\nOPER\nVISC 1000000\nMACH 0\nITER 1\nALFA {alpha_deg:.1}\n\nQUIT\n");
    let xfoil_stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: 3,
            warmup_runs: 1,
            sample_runs: 5,
        },
        || {
            let _ = run_xfoil_instrumented(&commands, &xfoil_workdir).expect("run xfoil binary");
        },
    );
    assert_ratio_within(
        &rust_stats,
        xfoil_stats.median_seconds,
        WORKFLOW_MAX_RATIO,
        "viscal_oper_point",
    );
}
