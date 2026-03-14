mod support;

use rustfoil_testkit::timing::{assert_ratio_within, benchmark_closure, BenchmarkConfig};
use rustfoil_xfoil::state_ops::{dsset, gamqv, iblpan, qvfue, stfind, uicalc, ueset, xicalc};

use support::{build_state_topology_fixture, prepare_state_topology_state, state_topology_fortran, MICROBENCH_MAX_RATIO};

fn release_only() {
    assert!(
        !cfg!(debug_assertions),
        "run perf gates with `cargo test --release -- --ignored`"
    );
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_stfind_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.stfind;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = build_state_topology_fixture();
            stfind(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "stfind");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_iblpan_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.iblpan;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = build_state_topology_fixture();
            stfind(&mut state);
            iblpan(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "iblpan");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_xicalc_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.xicalc;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = build_state_topology_fixture();
            stfind(&mut state);
            iblpan(&mut state);
            xicalc(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "xicalc");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_uicalc_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.uicalc;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = build_state_topology_fixture();
            stfind(&mut state);
            iblpan(&mut state);
            xicalc(&mut state);
            uicalc(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "uicalc");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_qvfue_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.qvfue;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = prepare_state_topology_state();
            qvfue(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "qvfue");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_gamqv_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.gamqv;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = prepare_state_topology_state();
            qvfue(&mut state);
            gamqv(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "gamqv");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_ueset_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.ueset;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = prepare_state_topology_state();
            ueset(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "ueset");
}

#[test]
#[ignore = "release-mode microbenchmark gate"]
fn perf_dsset_vs_fortran() {
    release_only();
    let reference = &state_topology_fortran().perf.dsset;
    let stats = benchmark_closure(
        BenchmarkConfig {
            inner_loops: reference.inner_loops,
            ..BenchmarkConfig::default()
        },
        || {
            let mut state = prepare_state_topology_state();
            ueset(&mut state);
            dsset(&mut state);
        },
    );
    assert_ratio_within(&stats, reference.median_seconds, MICROBENCH_MAX_RATIO, "dsset");
}
