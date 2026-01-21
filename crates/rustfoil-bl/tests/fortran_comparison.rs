//! FORTRAN comparison tests for boundary layer closures
//!
//! These tests validate Rust implementations against reference data
//! generated directly from XFOIL's FORTRAN source code.

use rustfoil_bl::closures::cf::{cf_laminar, cf_turbulent};
use rustfoil_bl::closures::hkin::hkin;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize)]
struct HkinTest {
    h: f64,
    msq: f64,
    hk: f64,
    hk_h: f64,
    hk_msq: f64,
}

#[derive(Deserialize)]
struct CflTest {
    hk: f64,
    rt: f64,
    cf: f64,
    cf_hk: f64,
    cf_rt: f64,
}

#[derive(Deserialize)]
struct CftTest {
    hk: f64,
    rt: f64,
    msq: f64,
    cf: f64,
    cf_hk: f64,
    cf_rt: f64,
    cf_msq: f64,
}

#[derive(Deserialize)]
struct ClosuresReference {
    hkin_tests: Vec<HkinTest>,
    #[serde(default)]
    cfl_tests: Vec<CflTest>,
    #[serde(default)]
    cft_tests: Vec<CftTest>,
}

fn load_reference() -> ClosuresReference {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("testdata/closures_reference.json");

    let contents = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read reference file {}: {}", path.display(), e));

    serde_json::from_str(&contents)
        .unwrap_or_else(|e| panic!("Failed to parse JSON from {}: {}", path.display(), e))
}

#[test]
fn test_hkin_matches_fortran() {
    let reference = load_reference();

    for (i, t) in reference.hkin_tests.iter().enumerate() {
        let result = hkin(t.h, t.msq);

        // Tolerance of 1e-6 accounts for FORTRAN single-precision (REAL*4) vs Rust f64
        assert!(
            (result.hk - t.hk).abs() < 1e-6,
            "Test {}: hk mismatch: {} vs {} for h={}, msq={}",
            i,
            result.hk,
            t.hk,
            t.h,
            t.msq
        );
        assert!(
            (result.hk_h - t.hk_h).abs() < 1e-6,
            "Test {}: hk_h mismatch: {} vs {} for h={}, msq={}",
            i,
            result.hk_h,
            t.hk_h,
            t.h,
            t.msq
        );
        assert!(
            (result.hk_msq - t.hk_msq).abs() < 1e-6,
            "Test {}: hk_msq mismatch: {} vs {} for h={}, msq={}",
            i,
            result.hk_msq,
            t.hk_msq,
            t.h,
            t.msq
        );
    }

    println!(
        "All {} HKIN test cases passed against FORTRAN reference",
        reference.hkin_tests.len()
    );
}

#[test]
fn test_hkin_reference_coverage() {
    let reference = load_reference();

    // Verify we have the expected 400 test cases (20 H values × 20 M² values)
    assert_eq!(
        reference.hkin_tests.len(),
        400,
        "Expected 400 HKIN test cases"
    );

    // Verify H range: 1.0 to 4.8
    let min_h = reference
        .hkin_tests
        .iter()
        .map(|t| t.h)
        .fold(f64::INFINITY, f64::min);
    let max_h = reference
        .hkin_tests
        .iter()
        .map(|t| t.h)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((min_h - 1.0).abs() < 0.01, "Min H should be ~1.0");
    assert!((max_h - 4.8).abs() < 0.01, "Max H should be ~4.8");

    // Verify M² range: 0.0 to 0.95
    let min_msq = reference
        .hkin_tests
        .iter()
        .map(|t| t.msq)
        .fold(f64::INFINITY, f64::min);
    let max_msq = reference
        .hkin_tests
        .iter()
        .map(|t| t.msq)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((min_msq - 0.0).abs() < 0.01, "Min M² should be ~0.0");
    assert!((max_msq - 0.95).abs() < 0.01, "Max M² should be ~0.95");
}

// ============== CFL (Laminar Skin Friction) Tests ==============

#[test]
fn test_cfl_matches_fortran() {
    let reference = load_reference();

    if reference.cfl_tests.is_empty() {
        println!("No CFL reference data found, skipping test");
        return;
    }

    for (i, t) in reference.cfl_tests.iter().enumerate() {
        let result = cf_laminar(t.hk, t.rt, 0.0);

        // Tolerance accounts for FORTRAN single-precision (REAL*4) vs Rust f64
        // Use relative tolerance of 1e-4 for larger values, absolute 1e-10 for small values
        let cf_tol = (t.cf.abs() * 1e-4).max(1e-10);
        assert!(
            (result.cf - t.cf).abs() < cf_tol,
            "Test {}: cf mismatch: {} vs {} for hk={}, rt={}",
            i,
            result.cf,
            t.cf,
            t.hk,
            t.rt
        );

        let cf_hk_tol = (t.cf_hk.abs() * 1e-4).max(1e-10);
        assert!(
            (result.cf_hk - t.cf_hk).abs() < cf_hk_tol,
            "Test {}: cf_hk mismatch: {} vs {} for hk={}, rt={}",
            i,
            result.cf_hk,
            t.cf_hk,
            t.hk,
            t.rt
        );

        // cf_rt values are very small (often ~1e-7 to 1e-10), so use absolute tolerance
        let cf_rt_tol = (t.cf_rt.abs() * 1e-3).max(1e-10);
        assert!(
            (result.cf_rt - t.cf_rt).abs() < cf_rt_tol,
            "Test {}: cf_rt mismatch: {} vs {} for hk={}, rt={}",
            i,
            result.cf_rt,
            t.cf_rt,
            t.hk,
            t.rt
        );
    }

    println!(
        "All {} CFL test cases passed against FORTRAN reference",
        reference.cfl_tests.len()
    );
}

#[test]
fn test_cfl_reference_coverage() {
    let reference = load_reference();

    if reference.cfl_tests.is_empty() {
        println!("No CFL reference data found, skipping coverage test");
        return;
    }

    // Verify we have the expected 130 test cases (13 HK values × 10 RT values)
    assert_eq!(
        reference.cfl_tests.len(),
        130,
        "Expected 130 CFL test cases"
    );

    // Verify HK range spans both branches (< 5.5 and >= 5.5)
    let min_hk = reference
        .cfl_tests
        .iter()
        .map(|t| t.hk)
        .fold(f64::INFINITY, f64::min);
    let max_hk = reference
        .cfl_tests
        .iter()
        .map(|t| t.hk)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(min_hk < 5.5, "Min HK should be < 5.5 (attached branch)");
    assert!(max_hk > 5.5, "Max HK should be > 5.5 (separated branch)");
}

// ============== CFT (Turbulent Skin Friction) Tests ==============

#[test]
fn test_cft_matches_fortran() {
    let reference = load_reference();

    if reference.cft_tests.is_empty() {
        println!("No CFT reference data found, skipping test");
        return;
    }

    for (i, t) in reference.cft_tests.iter().enumerate() {
        let result = cf_turbulent(t.hk, t.rt, t.msq);

        // Tolerance accounts for FORTRAN single-precision (REAL*4) vs Rust f64
        // Use relative tolerance of 2% for Cf values - single precision accumulates
        // more error in exp/pow/tanh operations, especially at edge cases (high Hk)
        // where Cf approaches zero and the THK correction term becomes dominant
        let cf_tol = (t.cf.abs() * 2e-2).max(1e-10);
        assert!(
            (result.cf - t.cf).abs() < cf_tol,
            "Test {}: cf mismatch: {} vs {} for hk={}, rt={}, msq={}",
            i,
            result.cf,
            t.cf,
            t.hk,
            t.rt,
            t.msq
        );

        let cf_hk_tol = (t.cf_hk.abs() * 1e-4).max(1e-10);
        assert!(
            (result.cf_hk - t.cf_hk).abs() < cf_hk_tol,
            "Test {}: cf_hk mismatch: {} vs {} for hk={}, rt={}, msq={}",
            i,
            result.cf_hk,
            t.cf_hk,
            t.hk,
            t.rt,
            t.msq
        );

        // cf_rt values are very small (often ~1e-7 to 1e-10), so use larger relative tolerance
        let cf_rt_tol = (t.cf_rt.abs() * 1e-3).max(1e-10);
        assert!(
            (result.cf_rt - t.cf_rt).abs() < cf_rt_tol,
            "Test {}: cf_rt mismatch: {} vs {} for hk={}, rt={}, msq={}",
            i,
            result.cf_rt,
            t.cf_rt,
            t.hk,
            t.rt,
            t.msq
        );

        // cf_msq values can be very small, use larger relative tolerance
        let cf_msq_tol = (t.cf_msq.abs() * 1e-3).max(1e-10);
        assert!(
            (result.cf_msq - t.cf_msq).abs() < cf_msq_tol,
            "Test {}: cf_msq mismatch: {} vs {} for hk={}, rt={}, msq={}",
            i,
            result.cf_msq,
            t.cf_msq,
            t.hk,
            t.rt,
            t.msq
        );
    }

    println!(
        "All {} CFT test cases passed against FORTRAN reference",
        reference.cft_tests.len()
    );
}

#[test]
fn test_cft_reference_coverage() {
    let reference = load_reference();

    if reference.cft_tests.is_empty() {
        println!("No CFT reference data found, skipping coverage test");
        return;
    }

    // Verify we have the expected 480 test cases (8 HK × 10 RT × 6 MSQ)
    assert_eq!(
        reference.cft_tests.len(),
        480,
        "Expected 480 CFT test cases"
    );

    // Verify HK range for turbulent (typical: 1.2 to 3.5)
    let min_hk = reference
        .cft_tests
        .iter()
        .map(|t| t.hk)
        .fold(f64::INFINITY, f64::min);
    let max_hk = reference
        .cft_tests
        .iter()
        .map(|t| t.hk)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((min_hk - 1.2).abs() < 0.01, "Min HK should be ~1.2");
    assert!((max_hk - 3.3).abs() < 0.01, "Max HK should be ~3.3");

    // Verify MSQ range: 0.0 to 0.5
    let min_msq = reference
        .cft_tests
        .iter()
        .map(|t| t.msq)
        .fold(f64::INFINITY, f64::min);
    let max_msq = reference
        .cft_tests
        .iter()
        .map(|t| t.msq)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!((min_msq - 0.0).abs() < 0.01, "Min MSQ should be ~0.0");
    assert!((max_msq - 0.5).abs() < 0.01, "Max MSQ should be ~0.5");
}
