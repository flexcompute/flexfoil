//! Inviscid solver debug tests for comparing RustFoil vs XFOIL.
//!
//! These tests load exact XFOIL panel coordinates and compare intermediate
//! values (geometry, influence coefficients, GAMU solutions) to identify
//! discrepancies in the inviscid solver.

use rustfoil_inviscid::{AirfoilGeometry, FlowConditions, InviscidSolver};
use rustfoil_inviscid::influence::{psilin, psilin_single_panel, psilin_debug};
use rustfoil_inviscid::system::build_system_matrix;
use serde::Deserialize;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

// ============================================================================
// Reference Data Structures
// ============================================================================

#[derive(Debug, Deserialize)]
struct XfoilReference {
    test: String,
    source: String,
    n_panels: usize,
    airfoil: String,
    data: Vec<serde_json::Value>,
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Load airfoil coordinates from XFOIL .dat file format.
fn load_xfoil_dat(filename: &str) -> Vec<(f64, f64)> {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.pop(); // crates/
    path.pop(); // project root
    path.push(filename);
    
    let content = fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Failed to read file: {:?}", path));
    
    let mut points = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                points.push((x, y));
            }
        }
    }
    
    points
}

/// Load reference JSON data.
fn load_reference_json() -> XfoilReference {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.pop();
    path.pop();
    path.push("testdata");
    path.push("inviscid_ref.json");
    
    let content = fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Failed to read reference file: {:?}", path));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse reference file: {}", e))
}

/// Extract specific subroutine data from reference.
fn get_subroutine_data<'a>(ref_data: &'a XfoilReference, name: &str) -> Option<&'a serde_json::Value> {
    ref_data.data.iter().find(|v| {
        v.get("subroutine").and_then(|s| s.as_str()) == Some(name)
    })
}

/// Print comparison table header.
fn print_comparison_header(title: &str) {
    println!("\n{}", "=".repeat(80));
    println!("{}", title);
    println!("{}", "=".repeat(80));
}

// ============================================================================
// Phase 1: Panel Geometry Validation
// ============================================================================

#[test]
fn test_phase1_exact_xfoil_geometry() {
    print_comparison_header("Phase 1: Panel Geometry Validation (Exact XFOIL Coordinates)");
    
    // Load exact XFOIL paneled coordinates
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    println!("Loaded {} points from XFOIL paneled file", points.len());
    
    // Create geometry
    let geom = AirfoilGeometry::from_points(&points)
        .expect("Failed to create geometry from XFOIL coordinates");
    
    println!("\n--- Basic Geometry Info ---");
    println!("  N nodes:    {}", geom.n);
    println!("  Chord:      {:.6}", geom.chord);
    println!("  XLE, YLE:   ({:.6}, {:.6})", geom.xle, geom.yle);
    println!("  DSTE:       {:.6e}", geom.dste);
    println!("  Sharp TE:   {}", geom.sharp);
    
    // Load reference data
    let ref_data = load_reference_json();
    
    // Compare TECALC values
    if let Some(tecalc) = get_subroutine_data(&ref_data, "TECALC") {
        let ref_ante = tecalc["ANTE"].as_f64().unwrap();
        let ref_aste = tecalc["ASTE"].as_f64().unwrap();
        let ref_dste = tecalc["DSTE"].as_f64().unwrap();
        let ref_sharp = tecalc["SHARP"].as_bool().unwrap();
        
        println!("\n--- TE Geometry Comparison ---");
        println!("  {:12} {:>16} {:>16} {:>12}", "Parameter", "RustFoil", "XFOIL", "Ratio");
        println!("  {:12} {:16.10e} {:16.10e} {:12.6}", "ANTE", geom.ante, ref_ante, 
            if ref_ante.abs() > 1e-15 { geom.ante / ref_ante } else { f64::NAN });
        println!("  {:12} {:16.10e} {:16.10e} {:12.6}", "ASTE", geom.aste, ref_aste,
            if ref_aste.abs() > 1e-15 { geom.aste / ref_aste } else { f64::NAN });
        println!("  {:12} {:16.10e} {:16.10e} {:12.6}", "DSTE", geom.dste, ref_dste, geom.dste / ref_dste);
        println!("  {:12} {:>16} {:>16}", "SHARP", geom.sharp, ref_sharp);
        
        // Assertions with tolerance
        assert!((geom.dste - ref_dste).abs() / ref_dste < 0.01, 
            "DSTE mismatch: {} vs {}", geom.dste, ref_dste);
    }
    
    // Compare normals (NCALC)
    if let Some(ncalc) = get_subroutine_data(&ref_data, "NCALC") {
        let ref_nx: Vec<f64> = ncalc["NX"].as_array().unwrap()
            .iter().map(|v| v.as_f64().unwrap()).collect();
        let ref_ny: Vec<f64> = ncalc["NY"].as_array().unwrap()
            .iter().map(|v| v.as_f64().unwrap()).collect();
        
        println!("\n--- Normal Vectors (NX, NY) First 10 nodes ---");
        println!("  {:>4} {:>14} {:>14} {:>14} {:>14}", "i", "NX_rust", "NX_xfoil", "NY_rust", "NY_xfoil");
        
        let mut nx_max_err: f64 = 0.0;
        let mut ny_max_err: f64 = 0.0;
        
        for i in 0..10.min(geom.n) {
            let nx_err = (geom.nx[i] - ref_nx[i]).abs();
            let ny_err = (geom.ny[i] - ref_ny[i]).abs();
            nx_max_err = nx_max_err.max(nx_err);
            ny_max_err = ny_max_err.max(ny_err);
            
            println!("  {:4} {:14.10} {:14.10} {:14.10} {:14.10}", 
                i, geom.nx[i], ref_nx[i], geom.ny[i], ref_ny[i]);
        }
        
        println!("\n  Max NX error: {:.6e}", nx_max_err);
        println!("  Max NY error: {:.6e}", ny_max_err);
        
        // Check all normals
        for i in 0..geom.n.min(ref_nx.len()) {
            let nx_err = (geom.nx[i] - ref_nx[i]).abs();
            let ny_err = (geom.ny[i] - ref_ny[i]).abs();
            nx_max_err = nx_max_err.max(nx_err);
            ny_max_err = ny_max_err.max(ny_err);
        }
        
        assert!(nx_max_err < 0.001, "NX max error too large: {:.6e}", nx_max_err);
        assert!(ny_max_err < 0.001, "NY max error too large: {:.6e}", ny_max_err);
    }
    
    // Compare panel angles (APCALC)
    if let Some(apcalc) = get_subroutine_data(&ref_data, "APCALC") {
        let ref_apanel: Vec<f64> = apcalc["APANEL"].as_array().unwrap()
            .iter().map(|v| v.as_f64().unwrap()).collect();
        
        println!("\n--- Panel Angles (APANEL) First 10 and Last 5 ---");
        println!("  {:>4} {:>14} {:>14} {:>14}", "i", "APANEL_rust", "APANEL_xfoil", "Diff");
        
        for i in 0..10.min(geom.n) {
            let diff = geom.apanel[i] - ref_apanel[i];
            println!("  {:4} {:14.10} {:14.10} {:14.6e}", i, geom.apanel[i], ref_apanel[i], diff);
        }
        println!("  ...");
        for i in (geom.n - 5)..geom.n {
            let diff = geom.apanel[i] - ref_apanel[i];
            println!("  {:4} {:14.10} {:14.10} {:14.6e}", i, geom.apanel[i], ref_apanel[i], diff);
        }
        
        let mut apanel_max_err: f64 = 0.0;
        for i in 0..geom.n.min(ref_apanel.len()) {
            let err = (geom.apanel[i] - ref_apanel[i]).abs();
            apanel_max_err = apanel_max_err.max(err);
        }
        println!("\n  Max APANEL error: {:.6e}", apanel_max_err);
        
        assert!(apanel_max_err < 0.01, "APANEL max error too large: {:.6e}", apanel_max_err);
    }
    
    println!("\n[PASS] Phase 1: Geometry validation passed");
}

// ============================================================================
// Phase 2: Influence Coefficient Debug
// ============================================================================

#[test]
fn test_phase2_influence_coefficients() {
    print_comparison_header("Phase 2: Influence Coefficient Analysis");
    
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let geom = AirfoilGeometry::from_points(&points).unwrap();
    
    println!("\n--- PSILIN Analysis for Field Point i=1 ---");
    println!("Analyzing influence from first 5 panels on node 1");
    println!();
    
    // Compute PSILIN for node 1 (second node)
    let i = 1;
    let xi = geom.x[i];
    let yi = geom.y[i];
    
    println!("Field point: i={}, (xi, yi) = ({:.6}, {:.6})", i, xi, yi);
    println!();
    
    // Analyze individual panel contributions
    println!("{:>4} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}", 
        "JO", "X1", "X2", "YY", "G1", "G2", "PSIS", "PSID");
    
    for jo in 0..5.min(geom.n - 1) {
        let contrib = psilin_single_panel(&geom, i, jo);
        println!("{:4} {:12.6} {:12.6} {:12.6} {:12.6} {:12.6} {:12.6} {:12.6}",
            jo, contrib.x1, contrib.x2, contrib.yy, contrib.g1, contrib.g2, contrib.psis, contrib.psid);
    }
    
    // Full PSILIN result
    let result = psilin(&geom, i, xi, yi);
    
    println!("\n--- DZDG coefficients (first 10) ---");
    println!("{:>4} {:>16}", "j", "DZDG[j]");
    for j in 0..10.min(geom.n) {
        println!("{:4} {:16.10e}", j, result.dzdg[j]);
    }
    
    // Self-influence at diagonal
    println!("\n--- Diagonal (self-influence) coefficients ---");
    for i in 0..5 {
        let result = psilin(&geom, i, geom.x[i], geom.y[i]);
        println!("  A[{},{}] = {:16.10e}", i, i, result.dzdg[i]);
    }
}

#[test]
fn test_phase2_detailed_debug() {
    print_comparison_header("Phase 2b: Detailed PSILIN Debug Output");
    
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let geom = AirfoilGeometry::from_points(&points).unwrap();
    
    // Test for field point i=0 (upper trailing edge)
    let i = 0;
    let debug = psilin_debug(&geom, i, geom.x[i], geom.y[i]);
    
    println!("\n--- Full Debug for i={} (upper TE) ---", i);
    println!("Field point: ({:.10}, {:.10})", debug.xi, debug.yi);
    
    // Print first few panels
    println!("\n{:>4} {:>4} {:>12} {:>12} {:>12} {:>12} {:>12} {:>14} {:>14}",
        "JO", "JP", "X1", "X2", "YY", "T1", "T2", "DZDG[jo]", "DZDG[jp]");
    
    for panel in debug.panels.iter().take(10) {
        println!("{:4} {:4} {:12.6} {:12.6} {:12.6} {:12.6} {:12.6} {:14.10e} {:14.10e}",
            panel.jo, panel.jp, 
            panel.contrib.x1, panel.contrib.x2, panel.contrib.yy,
            panel.contrib.t1, panel.contrib.t2,
            panel.dzdg_jo, panel.dzdg_jp);
    }
    
    // Print panels around diagonal (i-1, i, i+1)
    println!("\n--- Panels affecting diagonal (self-influence) ---");
    for panel in &debug.panels {
        if panel.jo == i || panel.jp == i || panel.jo == i.saturating_sub(1) || panel.jp == (i + 1) % geom.n {
            println!("Panel jo={}, jp={}: PSIS={:.10}, PSID={:.10}, dzdg_jo={:.10e}, dzdg_jp={:.10e}",
                panel.jo, panel.jp, panel.contrib.psis, panel.contrib.psid, panel.dzdg_jo, panel.dzdg_jp);
        }
    }
    
    // Final DZDG at diagonal
    println!("\n--- Final DZDG[{}] (diagonal) = {:.10e}", i, debug.dzdg[i]);
    
    // Test for field point i=1
    println!("\n\n--- Full Debug for i=1 ---");
    let debug1 = psilin_debug(&geom, 1, geom.x[1], geom.y[1]);
    println!("Field point: ({:.10}, {:.10})", debug1.xi, debug1.yi);
    println!("DZDG[1] (diagonal) = {:.10e}", debug1.dzdg[1]);
    
    // Compare with panels jo=0 and jo=1
    for panel in &debug1.panels {
        if panel.jo <= 2 {
            println!("Panel jo={}, jp={}: X1={:.6}, X2={:.6}, YY={:.6}, G1={:.6}, G2={:.6}",
                panel.jo, panel.jp, 
                panel.contrib.x1, panel.contrib.x2, panel.contrib.yy,
                panel.contrib.g1, panel.contrib.g2);
            println!("  PSIS={:.10}, PSID={:.10}", panel.contrib.psis, panel.contrib.psid);
            println!("  dzdg[{}] += {:.10e}, dzdg[{}] += {:.10e}", 
                panel.jo, panel.dzdg_jo, panel.jp, panel.dzdg_jp);
        }
    }
}

// ============================================================================
// Phase 3: GAMU Base Solution Comparison
// ============================================================================

#[test]
fn test_phase3_gamu_comparison() {
    print_comparison_header("Phase 3: GAMU Base Solution Comparison");
    
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let ref_data = load_reference_json();
    
    // Get XFOIL GAMU values
    let ggcalc = get_subroutine_data(&ref_data, "GGCALC")
        .expect("No GGCALC data in reference");
    
    let ref_gamu_0: Vec<f64> = ggcalc["GAMU_0"].as_array().unwrap()
        .iter().map(|v| v.as_f64().unwrap()).collect();
    let ref_gamu_90: Vec<f64> = ggcalc["GAMU_90"].as_array().unwrap()
        .iter().map(|v| v.as_f64().unwrap()).collect();
    
    println!("XFOIL reference: {} GAMU_0 values, {} GAMU_90 values", 
        ref_gamu_0.len(), ref_gamu_90.len());
    
    // Solve with RustFoil
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&points).expect("Factorization failed");
    
    println!("\n--- GAMU_0 Comparison (alpha=0 solution) ---");
    println!("{:>4} {:>16} {:>16} {:>12}", "i", "RustFoil", "XFOIL", "Ratio");
    
    let n_compare = ref_gamu_0.len().min(factorized.gamu_0.len());
    let mut gamu0_max_err: f64 = 0.0;
    let mut gamu0_sum_sq: f64 = 0.0;
    
    for i in 0..n_compare {
        let ratio = if ref_gamu_0[i].abs() > 1e-15 {
            factorized.gamu_0[i] / ref_gamu_0[i]
        } else {
            f64::NAN
        };
        let err = (factorized.gamu_0[i] - ref_gamu_0[i]).abs();
        gamu0_max_err = gamu0_max_err.max(err);
        gamu0_sum_sq += (factorized.gamu_0[i] - ref_gamu_0[i]).powi(2);
        
        if i < 20 {
            println!("{:4} {:16.10} {:16.10} {:12.6}", i, factorized.gamu_0[i], ref_gamu_0[i], ratio);
        }
    }
    let gamu0_rms = (gamu0_sum_sq / n_compare as f64).sqrt();
    println!("\n  GAMU_0 Max error: {:.6e}", gamu0_max_err);
    println!("  GAMU_0 RMS error: {:.6e}", gamu0_rms);
    
    println!("\n--- GAMU_90 Comparison (alpha=90 solution) - THE KEY DIAGNOSTIC ---");
    println!("{:>4} {:>16} {:>16} {:>12}", "i", "RustFoil", "XFOIL", "Ratio");
    
    let n_compare_90 = ref_gamu_90.len().min(factorized.gamu_90.len());
    let mut gamu90_max_err: f64 = 0.0;
    let mut gamu90_sum_sq: f64 = 0.0;
    let mut gamu90_ratios = Vec::new();
    
    for i in 0..n_compare_90 {
        let ratio = if ref_gamu_90[i].abs() > 1e-10 {
            factorized.gamu_90[i] / ref_gamu_90[i]
        } else {
            f64::NAN
        };
        if ratio.is_finite() {
            gamu90_ratios.push(ratio);
        }
        
        let err = (factorized.gamu_90[i] - ref_gamu_90[i]).abs();
        gamu90_max_err = gamu90_max_err.max(err);
        gamu90_sum_sq += (factorized.gamu_90[i] - ref_gamu_90[i]).powi(2);
        
        if i < 20 {
            println!("{:4} {:16.10} {:16.10} {:12.6}", i, factorized.gamu_90[i], ref_gamu_90[i], ratio);
        }
    }
    let gamu90_rms = (gamu90_sum_sq / n_compare_90 as f64).sqrt();
    
    let avg_ratio = if !gamu90_ratios.is_empty() {
        gamu90_ratios.iter().sum::<f64>() / gamu90_ratios.len() as f64
    } else {
        f64::NAN
    };
    
    println!("\n  GAMU_90 Max error: {:.6e}", gamu90_max_err);
    println!("  GAMU_90 RMS error: {:.6e}", gamu90_rms);
    println!("  GAMU_90 Average ratio: {:.4} (should be ~1.0)", avg_ratio);
    
    // Key diagnostic: If ratio is ~0.3, we have the bug described in Obsidian
    if avg_ratio < 0.5 {
        println!("\n  ⚠️  WARNING: GAMU_90 ratio is {:.4}, suggesting the bug persists!", avg_ratio);
        println!("      Expected: ratio ~1.0");
        println!("      Actual:   ratio ~{:.2}", avg_ratio);
    }
    
    // Check Kutta condition
    let n = factorized.gamu_0.len();
    let kutta_0 = factorized.gamu_0[0] + factorized.gamu_0[n - 1];
    let kutta_90 = factorized.gamu_90[0] + factorized.gamu_90[n - 1];
    
    println!("\n--- Kutta Condition Check ---");
    println!("  γ₀[0] + γ₀[N-1] = {:.6e} (should be ~0)", kutta_0);
    println!("  γ₉₀[0] + γ₉₀[N-1] = {:.6e} (should be ~0)", kutta_90);
    
    assert!(kutta_0.abs() < 1e-8, "Kutta violated for alpha=0: {}", kutta_0);
    assert!(kutta_90.abs() < 1e-8, "Kutta violated for alpha=90: {}", kutta_90);
}

// ============================================================================
// Phase 4: Multi-Alpha Validation
// ============================================================================

#[test]
fn test_phase4_multi_alpha_cl() {
    print_comparison_header("Phase 4: CL vs Alpha Comparison");
    
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let ref_data = load_reference_json();
    
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&points).expect("Factorization failed");
    
    // Get reference CLCALC at alpha=0 (for future comparison)
    let _ref_cl_0: Option<f64> = ref_data.data.iter()
        .filter(|ev| ev.get("subroutine").and_then(|s| s.as_str()) == Some("CLCALC"))
        .filter_map(|ev| ev.get("CL").and_then(|v| v.as_f64()))
        .find(|&cl| cl.abs() < 0.01);
    
    println!("\n{:>8} {:>12} {:>12} {:>12}", "α (deg)", "CL_rust", "CL_xfoil*", "Diff");
    
    let alphas = [0.0, 2.0, 4.0, 6.0, 8.0];
    for &alpha_deg in &alphas {
        let flow = FlowConditions::with_alpha_deg(alpha_deg);
        let solution = factorized.solve_alpha(&flow);
        
        // Approximate expected CL (thin airfoil theory: CL = 2π·α)
        let cl_theory = 2.0 * PI * alpha_deg.to_radians();
        
        println!("{:8.1} {:12.6} {:12.6} {:12.6}", 
            alpha_deg, solution.cl, cl_theory, solution.cl - cl_theory);
    }
    
    // Lift curve slope
    let flow_0 = FlowConditions::with_alpha_deg(0.0);
    let flow_8 = FlowConditions::with_alpha_deg(8.0);
    let cl_0 = factorized.solve_alpha(&flow_0).cl;
    let cl_8 = factorized.solve_alpha(&flow_8).cl;
    let cl_alpha = (cl_8 - cl_0) / 8.0_f64.to_radians();
    
    println!("\n  Lift curve slope: {:.4}/rad (expected ~{:.4})", cl_alpha, 2.0 * PI);
    println!("  Error: {:.1}%", (cl_alpha - 2.0 * PI).abs() / (2.0 * PI) * 100.0);
    
    // Check that CL at alpha=0 is near zero for symmetric airfoil
    assert!(cl_0.abs() < 0.1, "CL at alpha=0 should be near 0 for symmetric airfoil, got {}", cl_0);
}

// ============================================================================
// Phase 5: Matrix System Analysis
// ============================================================================

#[test]
fn test_phase5_matrix_analysis() {
    print_comparison_header("Phase 5: Influence Matrix Analysis");
    
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let geom = AirfoilGeometry::from_points(&points).unwrap();
    let ref_data = load_reference_json();
    
    // Build system matrix
    let (a_matrix, rhs_0, rhs_90) = build_system_matrix(&geom);
    
    println!("Matrix size: {}x{}", a_matrix.nrows(), a_matrix.ncols());
    
    // Compare diagonal
    if let Some(ggcalc) = get_subroutine_data(&ref_data, "GGCALC") {
        let ref_diag: Vec<f64> = ggcalc["AIJ_diagonal"].as_array().unwrap()
            .iter().map(|v| v.as_f64().unwrap()).collect();
        
        println!("\n--- Matrix Diagonal Comparison ---");
        println!("{:>4} {:>16} {:>16} {:>12}", "i", "A[i,i]_rust", "A[i,i]_xfoil", "Ratio");
        
        for i in 0..ref_diag.len().min(20) {
            let rust_val = a_matrix[(i, i)];
            let ratio = if ref_diag[i].abs() > 1e-15 { rust_val / ref_diag[i] } else { f64::NAN };
            println!("{:4} {:16.10e} {:16.10e} {:12.6}", i, rust_val, ref_diag[i], ratio);
        }
    }
    
    // Check RHS vectors
    println!("\n--- RHS Vectors (first 10) ---");
    println!("{:>4} {:>16} {:>16} {:>16}", "i", "y[i]", "RHS_0[i]", "RHS_90[i]");
    for i in 0..10 {
        println!("{:4} {:16.10} {:16.10} {:16.10}", 
            i, geom.y[i], rhs_0[i], rhs_90[i]);
    }
    
    // Verify RHS is -y for alpha=0 and +x for alpha=90
    println!("\n  Checking RHS construction:");
    println!("    RHS_0[0] = {:.10} should equal -y[0] = {:.10}", rhs_0[0], -geom.y[0]);
    println!("    RHS_90[0] = {:.10} should equal x[0] = {:.10}", rhs_90[0], geom.x[0]);
}

// ============================================================================
// Full Panel Comparison
// ============================================================================

/// Load XFOIL dump file (from DUMP command in OPER menu)
fn load_xfoil_dump(filename: &str) -> Option<Vec<XfoilDumpRow>> {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.pop();
    path.pop();
    path.push(filename);
    
    let content = fs::read_to_string(&path).ok()?;
    let mut rows = Vec::new();
    
    for line in content.lines().skip(1) {  // Skip header
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 5 {
            rows.push(XfoilDumpRow {
                s: parts[0].parse().unwrap_or(0.0),
                x: parts[1].parse().unwrap_or(0.0),
                y: parts[2].parse().unwrap_or(0.0),
                ue_vinf: parts[3].parse().unwrap_or(0.0),  // This is gamma (vortex strength)
            });
        }
    }
    Some(rows)
}

#[derive(Debug)]
struct XfoilDumpRow {
    s: f64,
    x: f64,
    y: f64,
    ue_vinf: f64,  // γ = Ue/Vinf (edge velocity / freestream)
}

#[test]
fn test_full_panel_comparison() {
    print_comparison_header("Full Panel Comparison: Gamma, Velocity, Cp");
    
    // Load XFOIL reference data
    let xfoil_data = load_xfoil_dump("naca0012_alpha0_dump.txt");
    if xfoil_data.is_none() {
        println!("  [SKIP] XFOIL dump file not found");
        return;
    }
    let xfoil_data = xfoil_data.unwrap();
    
    // Compute RustFoil solution
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&points).expect("Factorization failed");
    let flow = FlowConditions::with_alpha_deg(0.0);
    let solution = factorized.solve_alpha(&flow);
    
    println!("\nNACA 0012, α=0°, {} panels", points.len());
    println!("\n{:>4} {:>10} {:>10} {:>12} {:>12} {:>10} {:>12} {:>12} {:>10}",
        "i", "x", "y", "γ_rust", "γ_xfoil", "γ_ratio", "Cp_rust", "Cp_xfoil", "Cp_diff");
    println!("{}", "-".repeat(110));
    
    let n = solution.gamma.len().min(xfoil_data.len());
    let mut gamma_max_err: f64 = 0.0;
    let mut gamma_rms: f64 = 0.0;
    let mut cp_max_err: f64 = 0.0;
    let mut cp_rms: f64 = 0.0;
    
    for i in 0..n {
        let gamma_rust = solution.gamma[i];
        let gamma_xfoil = xfoil_data[i].ue_vinf;
        let gamma_ratio = if gamma_xfoil.abs() > 1e-10 { gamma_rust / gamma_xfoil } else { f64::NAN };
        let gamma_err = (gamma_rust - gamma_xfoil).abs();
        gamma_max_err = gamma_max_err.max(gamma_err);
        gamma_rms += (gamma_rust - gamma_xfoil).powi(2);
        
        // Cp = 1 - (Ue/Vinf)^2 = 1 - γ²
        let cp_rust = solution.cp[i];
        let cp_xfoil = 1.0 - gamma_xfoil * gamma_xfoil;
        let cp_diff = cp_rust - cp_xfoil;
        cp_max_err = cp_max_err.max(cp_diff.abs());
        cp_rms += cp_diff.powi(2);
        
        // Print every 10th row, plus first 5 and last 5
        if i < 5 || i >= n - 5 || i % 10 == 0 {
            println!("{:4} {:10.6} {:10.6} {:12.6} {:12.6} {:10.6} {:12.6} {:12.6} {:10.6e}",
                i, xfoil_data[i].x, xfoil_data[i].y, 
                gamma_rust, gamma_xfoil, gamma_ratio,
                cp_rust, cp_xfoil, cp_diff);
        }
        if i == 5 { println!("  ..."); }
        if i == n - 6 { println!("  ..."); }
    }
    
    gamma_rms = (gamma_rms / n as f64).sqrt();
    cp_rms = (cp_rms / n as f64).sqrt();
    
    println!("{}", "-".repeat(110));
    println!("\n--- Statistics ---");
    println!("  Gamma max error: {:.6e}", gamma_max_err);
    println!("  Gamma RMS error: {:.6e}", gamma_rms);
    println!("  Cp max error:    {:.6e}", cp_max_err);
    println!("  Cp RMS error:    {:.6e}", cp_rms);
    
    // Also show α=4° comparison if available
    let xfoil_data_4 = load_xfoil_dump("naca0012_alpha4_dump.txt");
    if let Some(xfoil_4) = xfoil_data_4 {
        let flow_4 = FlowConditions::with_alpha_deg(4.0);
        let sol_4 = factorized.solve_alpha(&flow_4);
        
        println!("\n\n--- NACA 0012, α=4° ---");
        println!("{:>4} {:>10} {:>10} {:>12} {:>12} {:>10}",
            "i", "x", "y", "γ_rust", "γ_xfoil", "γ_ratio");
        println!("{}", "-".repeat(70));
        
        let n4 = sol_4.gamma.len().min(xfoil_4.len());
        let mut g4_max_err: f64 = 0.0;
        
        for i in 0..n4 {
            let gr = sol_4.gamma[i];
            let gx = xfoil_4[i].ue_vinf;
            let ratio = if gx.abs() > 1e-10 { gr / gx } else { f64::NAN };
            let err = (gr - gx).abs();
            g4_max_err = g4_max_err.max(err);
            
            if i < 5 || i >= n4 - 5 || i % 20 == 0 {
                println!("{:4} {:10.6} {:10.6} {:12.6} {:12.6} {:10.6}",
                    i, xfoil_4[i].x, xfoil_4[i].y, gr, gx, ratio);
            }
        }
        println!("\n  α=4° Gamma max error: {:.6e}", g4_max_err);
    }
    
    // Assertions
    assert!(gamma_max_err < 0.001, "Gamma max error too large: {}", gamma_max_err);
    assert!(cp_max_err < 0.002, "Cp max error too large: {}", cp_max_err);
}

#[test]
fn test_full_panel_all_values() {
    print_comparison_header("ALL Panel Values: Gamma, Velocity, Cp at α=0°");
    
    // Load XFOIL reference data
    let xfoil_data = load_xfoil_dump("naca0012_alpha0_dump.txt");
    if xfoil_data.is_none() {
        println!("  [SKIP] XFOIL dump file not found");
        return;
    }
    let xfoil_data = xfoil_data.unwrap();
    
    // Compute RustFoil solution
    let points = load_xfoil_dat("naca0012_xfoil_paneled.dat");
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&points).expect("Factorization failed");
    let flow = FlowConditions::with_alpha_deg(0.0);
    let solution = factorized.solve_alpha(&flow);
    
    println!("\nNACA 0012, α=0°, 160 panels - ALL VALUES");
    println!("\n{:>4} {:>10} {:>10} {:>12} {:>12} {:>12} {:>12}",
        "i", "x", "y", "γ_rust", "γ_xfoil", "Cp_rust", "Cp_xfoil");
    println!("{}", "=".repeat(85));
    
    let n = solution.gamma.len().min(xfoil_data.len());
    
    for i in 0..n {
        let gamma_rust = solution.gamma[i];
        let gamma_xfoil = xfoil_data[i].ue_vinf;
        let cp_rust = solution.cp[i];
        let cp_xfoil = 1.0 - gamma_xfoil * gamma_xfoil;
        
        println!("{:4} {:10.6} {:10.6} {:12.6} {:12.6} {:12.6} {:12.6}",
            i, xfoil_data[i].x, xfoil_data[i].y, 
            gamma_rust, gamma_xfoil, cp_rust, cp_xfoil);
    }
    println!("{}", "=".repeat(85));
}

// ============================================================================
// Phase 6: Multi-Foil Validation
// ============================================================================

#[test]
fn test_phase6_multi_foil() {
    print_comparison_header("Phase 6: Multi-Foil Validation");
    
    // Test NACA 2412 (cambered)
    let naca2412 = load_xfoil_dat("naca2412_xfoil_paneled.dat");
    if !naca2412.is_empty() {
        let solver = InviscidSolver::new();
        let factorized = solver.factorize(&naca2412).expect("2412 factorization failed");
        
        println!("\n--- NACA 2412 ---");
        let alphas = [-2.0, 0.0, 2.0, 4.0, 6.0];
        println!("{:>8} {:>12} {:>12}", "α (deg)", "CL", "CM");
        
        for &alpha in &alphas {
            let flow = FlowConditions::with_alpha_deg(alpha);
            let sol = factorized.solve_alpha(&flow);
            println!("{:8.1} {:12.6} {:12.6}", alpha, sol.cl, sol.cm);
        }
        
        // Check CL at alpha=0 (should be non-zero for cambered airfoil)
        let flow_0 = FlowConditions::with_alpha_deg(0.0);
        let cl_0 = factorized.solve_alpha(&flow_0).cl;
        println!("\n  CL at α=0: {:.4} (should be ~0.26 for 2% camber)", cl_0);
        
        // Check zero-lift alpha (should be ~-2° for 2% camber)
        // Linear interpolation between α=-2 and α=0
        let flow_m2 = FlowConditions::with_alpha_deg(-2.0);
        let cl_m2 = factorized.solve_alpha(&flow_m2).cl;
        let alpha_0l = -2.0 + 2.0 * (-cl_m2) / (cl_0 - cl_m2);
        println!("  Zero-lift α: {:.2}° (should be ~-2°)", alpha_0l);
    } else {
        println!("  Skipping NACA 2412 (buffer file not found)");
    }
    
    // Test NACA 4412 (high camber)
    let naca4412 = load_xfoil_dat("naca4412_xfoil_paneled.dat");
    if !naca4412.is_empty() {
        let solver = InviscidSolver::new();
        let factorized = solver.factorize(&naca4412).expect("4412 factorization failed");
        
        println!("\n--- NACA 4412 ---");
        let alphas = [-4.0, 0.0, 4.0, 8.0];
        println!("{:>8} {:>12} {:>12}", "α (deg)", "CL", "CM");
        
        for &alpha in &alphas {
            let flow = FlowConditions::with_alpha_deg(alpha);
            let sol = factorized.solve_alpha(&flow);
            println!("{:8.1} {:12.6} {:12.6}", alpha, sol.cl, sol.cm);
        }
        
        // Check CL at alpha=0 (should be ~0.5 for 4% camber)
        let flow_0 = FlowConditions::with_alpha_deg(0.0);
        let cl_0 = factorized.solve_alpha(&flow_0).cl;
        println!("\n  CL at α=0: {:.4} (should be ~0.52 for 4% camber)", cl_0);
    } else {
        println!("  Skipping NACA 4412 (buffer file not found)");
    }
    
    println!("\n[PASS] Phase 6: Multi-foil validation");
}

// ============================================================================
// Run all phases
// ============================================================================

#[test]
fn test_run_all_debug_phases() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                    INVISCID SOLVER DEBUG SUITE                               ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    
    // Run tests manually to get output
    test_phase1_exact_xfoil_geometry();
    test_phase2_influence_coefficients();
    test_phase3_gamu_comparison();
    test_phase4_multi_alpha_cl();
    test_phase5_matrix_analysis();
    
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                         DEBUG SUITE COMPLETE                                 ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
}
