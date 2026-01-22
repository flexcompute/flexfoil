//! Matrix comparison test - compares XFOIL's VS2 matrix with RustFoil's build_4x4_system
//!
//! This test loads the XFOIL debug trace and compares:
//! 1. The 4x4 Jacobian matrix (VS2) before GAUSS
//! 2. The RHS vector (VSREZ) before GAUSS
//! 3. The solution vector (VSREZ) after GAUSS

use rustfoil_bl::equations::{bldif, blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::solve::{build_4x4_system, solve_4x4};
use serde::Deserialize;
use std::fs;

/// Load the XFOIL debug trace
fn load_debug_trace() -> Option<serde_json::Value> {
    let paths = [
        "testdata/newton_debug_trace.json",
        "../testdata/newton_debug_trace.json",
        "../../testdata/newton_debug_trace.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(data) = serde_json::from_str(&content) {
                return Some(data);
            }
        }
    }
    None
}

/// Load the station reference data
fn load_station_reference() -> Option<serde_json::Value> {
    let paths = [
        "testdata/mrchue_iterations.json",
        "../testdata/mrchue_iterations.json",
        "../../testdata/mrchue_iterations.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(data) = serde_json::from_str(&content) {
                return Some(data);
            }
        }
    }
    None
}

/// Extract VS2_BEFORE events from debug trace
fn extract_vs2_before_events(trace: &serde_json::Value) -> Vec<serde_json::Value> {
    let events = trace["events"].as_array().unwrap();
    events
        .iter()
        .filter(|e| e["subroutine"].as_str() == Some("VS2_BEFORE"))
        .cloned()
        .collect()
}

/// Extract VSREZ_AFTER events from debug trace
fn extract_vsrez_after_events(trace: &serde_json::Value) -> Vec<serde_json::Value> {
    let events = trace["events"].as_array().unwrap();
    events
        .iter()
        .filter(|e| e["subroutine"].as_str() == Some("VSREZ_AFTER"))
        .cloned()
        .collect()
}

/// Parse a 4x4 matrix from JSON
fn parse_matrix_4x4(json: &serde_json::Value) -> [[f64; 4]; 4] {
    let mut m = [[0.0; 4]; 4];
    if let Some(rows) = json.as_array() {
        for (i, row) in rows.iter().enumerate() {
            if let Some(cols) = row.as_array() {
                for (j, val) in cols.iter().enumerate() {
                    if i < 4 && j < 4 {
                        m[i][j] = val.as_f64().unwrap_or(0.0);
                    }
                }
            }
        }
    }
    m
}

/// Parse a 4-element vector from JSON
fn parse_vector_4(json: &serde_json::Value) -> [f64; 4] {
    let mut v = [0.0; 4];
    if let Some(arr) = json.as_array() {
        for (i, val) in arr.iter().enumerate() {
            if i < 4 {
                v[i] = val.as_f64().unwrap_or(0.0);
            }
        }
    }
    v
}

#[test]
fn test_matrix_comparison_ibl2_iter1() {
    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let station_ref = match load_station_reference() {
        Some(r) => r,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let vs2_before_events = extract_vs2_before_events(&trace);
    let vsrez_after_events = extract_vsrez_after_events(&trace);

    println!("\n=== Matrix Comparison: IBL=2, Iteration 1 ===\n");

    // Find the first event for side=1, ibl=2, newton_iter=1
    let xfoil_before = vs2_before_events
        .iter()
        .find(|e| {
            e["side"].as_i64() == Some(1)
                && e["ibl"].as_i64() == Some(2)
                && e["newton_iter"].as_i64() == Some(1)
        });

    let xfoil_after = vsrez_after_events
        .iter()
        .find(|e| {
            e["side"].as_i64() == Some(1)
                && e["ibl"].as_i64() == Some(2)
                && e["newton_iter"].as_i64() == Some(1)
        });

    let xfoil_before = match xfoil_before {
        Some(e) => e,
        None => {
            eprintln!("VS2_BEFORE event not found for IBL=2 iter=1");
            return;
        }
    };

    let xfoil_after = match xfoil_after {
        Some(e) => e,
        None => {
            eprintln!("VSREZ_AFTER event not found for IBL=2 iter=1");
            return;
        }
    };

    // Parse XFOIL matrices
    let xfoil_vs2 = parse_matrix_4x4(&xfoil_before["VS2_4x4"]);
    let xfoil_rhs = parse_vector_4(&xfoil_before["VSREZ_rhs"]);
    let xfoil_solution = parse_vector_4(&xfoil_after["VSREZ_solution"]);

    println!("XFOIL VS2 Matrix (before GAUSS):");
    for i in 0..4 {
        println!(
            "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
            xfoil_vs2[i][0], xfoil_vs2[i][1], xfoil_vs2[i][2], xfoil_vs2[i][3]
        );
    }
    println!("\nXFOIL RHS (before GAUSS):");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        xfoil_rhs[0], xfoil_rhs[1], xfoil_rhs[2], xfoil_rhs[3]
    );
    println!("\nXFOIL Solution (after GAUSS):");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        xfoil_solution[0], xfoil_solution[1], xfoil_solution[2], xfoil_solution[3]
    );

    // Now build RustFoil's system for comparison
    // For IBL=2, we need to use the stagnation initialization
    let re = station_ref["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;
    let stations = &station_ref["sides"]["1"]["stations"];
    let s0 = &stations[0]; // IBL=2

    // Build stagnation point (IBL=1) as "prev"
    let mut prev = BlStation::new();
    // Use Thwaites formula for stagnation initialization
    let x_init = s0["x"].as_f64().unwrap();
    let ue_init = s0["Ue"].as_f64().unwrap();
    
    // Thwaites: theta ~ sqrt(0.45 * x / (6 * Ue * Re))
    let theta_init = (0.45 * x_init / (6.0 * ue_init * re)).sqrt();
    prev.x = x_init * 0.5; // Approximate stagnation point
    prev.u = ue_init * 0.5;
    prev.theta = theta_init;
    prev.delta_star = 2.2 * theta_init;
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    // Build "curr" station (IBL=2) with initial guess
    let mut curr = BlStation::new();
    curr.x = s0["x"].as_f64().unwrap();
    curr.u = s0["Ue"].as_f64().unwrap();
    curr.theta = s0["initial"]["theta"].as_f64().unwrap();
    curr.delta_star = s0["initial"]["delta_star"].as_f64().unwrap();
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== RustFoil Computation ===\n");
    println!("prev: x={:.6e}, Ue={:.6e}, θ={:.6e}, δ*={:.6e}", prev.x, prev.u, prev.theta, prev.delta_star);
    println!("curr: x={:.6e}, Ue={:.6e}, θ={:.6e}, δ*={:.6e}", curr.x, curr.u, curr.theta, curr.delta_star);

    // Compute residuals and Jacobian
    let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);

    // Build 4x4 system
    let (rust_a, rust_b) = build_4x4_system(
        &jac.vs2,
        &[res.res_third, res.res_mom, res.res_shape],
        true, // direct mode
        curr.hk / curr.theta,
        -curr.hk / curr.delta_star,
        0.0,
        4.0, // hlmax
        curr.hk,
    );

    println!("\nRustFoil VS2 Matrix:");
    for i in 0..4 {
        println!(
            "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
            rust_a[i][0], rust_a[i][1], rust_a[i][2], rust_a[i][3]
        );
    }
    println!("\nRustFoil RHS:");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        rust_b[0], rust_b[1], rust_b[2], rust_b[3]
    );

    // Solve
    let rust_solution = solve_4x4(&rust_a, &rust_b);
    println!("\nRustFoil Solution:");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        rust_solution[0], rust_solution[1], rust_solution[2], rust_solution[3]
    );

    // Compare matrices element by element
    println!("\n=== Element-by-Element Comparison ===\n");
    println!("Matrix A (Jacobian):");
    for i in 0..4 {
        for j in 0..4 {
            let xfoil_val = xfoil_vs2[i][j];
            let rust_val = rust_a[i][j];
            let diff = rust_val - xfoil_val;
            let rel_err = if xfoil_val.abs() > 1e-30 {
                (diff / xfoil_val).abs() * 100.0
            } else {
                if rust_val.abs() > 1e-30 { 100.0 } else { 0.0 }
            };
            
            let sign_match = (xfoil_val >= 0.0) == (rust_val >= 0.0);
            let sign_str = if sign_match { "  " } else { "!!" };
            
            if rel_err > 1.0 || !sign_match {
                println!(
                    "  A[{}][{}]: XFOIL={:+12.5e}, Rust={:+12.5e}, err={:6.1}% {}",
                    i, j, xfoil_val, rust_val, rel_err, sign_str
                );
            }
        }
    }

    println!("\nRHS (b vector):");
    for i in 0..4 {
        let xfoil_val = xfoil_rhs[i];
        let rust_val = rust_b[i];
        let diff = rust_val - xfoil_val;
        let rel_err = if xfoil_val.abs() > 1e-30 {
            (diff / xfoil_val).abs() * 100.0
        } else {
            if rust_val.abs() > 1e-30 { 100.0 } else { 0.0 }
        };
        
        let sign_match = (xfoil_val >= 0.0) == (rust_val >= 0.0);
        let sign_str = if sign_match { "  " } else { "!!" };
        
        println!(
            "  b[{}]: XFOIL={:+12.5e}, Rust={:+12.5e}, err={:6.1}% {}",
            i, xfoil_val, rust_val, rel_err, sign_str
        );
    }

    println!("\nSolution (dx vector):");
    for i in 0..4 {
        let xfoil_val = xfoil_solution[i];
        let rust_val = rust_solution[i];
        let diff = rust_val - xfoil_val;
        let rel_err = if xfoil_val.abs() > 1e-30 {
            (diff / xfoil_val).abs() * 100.0
        } else {
            if rust_val.abs() > 1e-30 { 100.0 } else { 0.0 }
        };
        
        let sign_match = (xfoil_val >= 0.0) == (rust_val >= 0.0);
        let sign_str = if sign_match { "  " } else { "!! SIGN MISMATCH" };
        
        println!(
            "  dx[{}]: XFOIL={:+12.5e}, Rust={:+12.5e}, err={:6.1}% {}",
            i, xfoil_val, rust_val, rel_err, sign_str
        );
    }

    // Summary
    println!("\n=== Summary ===");
    println!("dx[1] = dθ: XFOIL={:+.6e}, Rust={:+.6e}", xfoil_solution[1], rust_solution[1]);
    println!("dx[2] = dδ*: XFOIL={:+.6e}, Rust={:+.6e}", xfoil_solution[2], rust_solution[2]);
    
    let theta_sign_ok = (xfoil_solution[1] >= 0.0) == (rust_solution[1] >= 0.0);
    let dstar_sign_ok = (xfoil_solution[2] >= 0.0) == (rust_solution[2] >= 0.0);
    
    if theta_sign_ok && dstar_sign_ok {
        println!("\nSigns match for θ and δ* updates.");
    } else {
        println!("\n*** SIGN MISMATCH DETECTED ***");
        if !theta_sign_ok {
            println!("  - θ update has wrong sign!");
        }
        if !dstar_sign_ok {
            println!("  - δ* update has wrong sign!");
        }
    }
}

#[test]
fn test_matrix_comparison_ibl3_iter1() {
    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let station_ref = match load_station_reference() {
        Some(r) => r,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let vs2_before_events = extract_vs2_before_events(&trace);
    let vsrez_after_events = extract_vsrez_after_events(&trace);

    println!("\n=== Matrix Comparison: IBL=3, Iteration 1 ===\n");

    // Find event for side=1, ibl=3, newton_iter=1
    let xfoil_before = vs2_before_events
        .iter()
        .find(|e| {
            e["side"].as_i64() == Some(1)
                && e["ibl"].as_i64() == Some(3)
                && e["newton_iter"].as_i64() == Some(1)
        });

    let xfoil_after = vsrez_after_events
        .iter()
        .find(|e| {
            e["side"].as_i64() == Some(1)
                && e["ibl"].as_i64() == Some(3)
                && e["newton_iter"].as_i64() == Some(1)
        });

    let xfoil_before = match xfoil_before {
        Some(e) => e,
        None => {
            eprintln!("VS2_BEFORE event not found for IBL=3 iter=1");
            return;
        }
    };

    let xfoil_after = match xfoil_after {
        Some(e) => e,
        None => {
            eprintln!("VSREZ_AFTER event not found for IBL=3 iter=1");
            return;
        }
    };

    // Parse XFOIL matrices
    let xfoil_vs2 = parse_matrix_4x4(&xfoil_before["VS2_4x4"]);
    let xfoil_rhs = parse_vector_4(&xfoil_before["VSREZ_rhs"]);
    let xfoil_solution = parse_vector_4(&xfoil_after["VSREZ_solution"]);

    println!("XFOIL VS2 Matrix (before GAUSS):");
    for i in 0..4 {
        println!(
            "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
            xfoil_vs2[i][0], xfoil_vs2[i][1], xfoil_vs2[i][2], xfoil_vs2[i][3]
        );
    }
    println!("\nXFOIL RHS (before GAUSS):");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        xfoil_rhs[0], xfoil_rhs[1], xfoil_rhs[2], xfoil_rhs[3]
    );
    println!("\nXFOIL Solution (after GAUSS): dx = [d_ampl, d_theta, d_dstar, d_ue]");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        xfoil_solution[0], xfoil_solution[1], xfoil_solution[2], xfoil_solution[3]
    );

    // Build RustFoil's system
    let re = station_ref["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;
    let stations = &station_ref["sides"]["1"]["stations"];
    let s0 = &stations[0]; // IBL=2 (prev)
    let s1 = &stations[1]; // IBL=3 (curr)

    // Build "prev" station (IBL=2 converged)
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    // Build "curr" station (IBL=3 initial)
    let mut curr = BlStation::new();
    curr.x = s1["x"].as_f64().unwrap();
    curr.u = s1["Ue"].as_f64().unwrap();
    // For first iteration, current station starts with prev station's values
    curr.theta = s1["initial"]["theta"].as_f64().unwrap();
    curr.delta_star = s1["initial"]["delta_star"].as_f64().unwrap();
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== RustFoil Computation ===\n");
    println!("prev (IBL=2): x={:.6e}, Ue={:.6e}, θ={:.6e}, δ*={:.6e}", prev.x, prev.u, prev.theta, prev.delta_star);
    println!("curr (IBL=3): x={:.6e}, Ue={:.6e}, θ={:.6e}, δ*={:.6e}", curr.x, curr.u, curr.theta, curr.delta_star);
    println!("curr Hk={:.4}, Rθ={:.4}", curr.hk, curr.r_theta);

    // Compute residuals and Jacobian
    let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);

    println!("\nResiduals from bldif:");
    println!("  res_third={:+.6e}, res_mom={:+.6e}, res_shape={:+.6e}", 
        res.res_third, res.res_mom, res.res_shape);

    // Build 4x4 system
    let (rust_a, rust_b) = build_4x4_system(
        &jac.vs2,
        &[res.res_third, res.res_mom, res.res_shape],
        true, // direct mode
        curr.hk / curr.theta,
        -curr.hk / curr.delta_star,
        0.0,
        4.0,
        curr.hk,
    );

    println!("\nRustFoil VS2 Matrix:");
    for i in 0..4 {
        println!(
            "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
            rust_a[i][0], rust_a[i][1], rust_a[i][2], rust_a[i][3]
        );
    }
    println!("\nRustFoil RHS:");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        rust_b[0], rust_b[1], rust_b[2], rust_b[3]
    );

    // Solve
    let rust_solution = solve_4x4(&rust_a, &rust_b);
    println!("\nRustFoil Solution:");
    println!(
        "  [{:+12.5e}, {:+12.5e}, {:+12.5e}, {:+12.5e}]",
        rust_solution[0], rust_solution[1], rust_solution[2], rust_solution[3]
    );

    // Compare
    println!("\n=== Element-by-Element Comparison ===\n");
    
    println!("Matrix A (rows where significant difference exists):");
    for i in 0..4 {
        for j in 0..4 {
            let xfoil_val = xfoil_vs2[i][j];
            let rust_val = rust_a[i][j];
            let rel_err = if xfoil_val.abs() > 1e-30 {
                ((rust_val - xfoil_val) / xfoil_val).abs() * 100.0
            } else {
                if rust_val.abs() > 1e-30 { 100.0 } else { 0.0 }
            };
            
            let sign_match = (xfoil_val >= 0.0) == (rust_val >= 0.0);
            
            if rel_err > 5.0 || !sign_match {
                println!(
                    "  A[{}][{}]: XFOIL={:+12.5e}, Rust={:+12.5e}, err={:6.1}%{}",
                    i, j, xfoil_val, rust_val, rel_err,
                    if sign_match { "" } else { " !! SIGN" }
                );
            }
        }
    }

    println!("\nRHS comparison:");
    for i in 0..4 {
        let xfoil_val = xfoil_rhs[i];
        let rust_val = rust_b[i];
        let rel_err = if xfoil_val.abs() > 1e-30 {
            ((rust_val - xfoil_val) / xfoil_val).abs() * 100.0
        } else {
            if rust_val.abs() > 1e-30 { 100.0 } else { 0.0 }
        };
        let sign_match = (xfoil_val >= 0.0) == (rust_val >= 0.0);
        
        println!(
            "  b[{}]: XFOIL={:+12.5e}, Rust={:+12.5e}, err={:6.1}%{}",
            i, xfoil_val, rust_val, rel_err,
            if sign_match { "" } else { " !! SIGN" }
        );
    }

    println!("\nSolution comparison:");
    for i in 0..4 {
        let xfoil_val = xfoil_solution[i];
        let rust_val = rust_solution[i];
        let rel_err = if xfoil_val.abs() > 1e-30 {
            ((rust_val - xfoil_val) / xfoil_val).abs() * 100.0
        } else {
            if rust_val.abs() > 1e-30 { 100.0 } else { 0.0 }
        };
        let sign_match = (xfoil_val >= 0.0) == (rust_val >= 0.0);
        
        let var_name = match i {
            0 => "d_ampl",
            1 => "d_theta",
            2 => "d_dstar",
            3 => "d_ue",
            _ => "?",
        };
        
        println!(
            "  dx[{}] ({}): XFOIL={:+12.5e}, Rust={:+12.5e}, err={:6.1}%{}",
            i, var_name, xfoil_val, rust_val, rel_err,
            if sign_match { "" } else { " !! SIGN MISMATCH !!" }
        );
    }

    // Final verdict
    println!("\n=== VERDICT ===");
    let theta_sign_ok = (xfoil_solution[1] >= 0.0) == (rust_solution[1] >= 0.0);
    let dstar_sign_ok = (xfoil_solution[2] >= 0.0) == (rust_solution[2] >= 0.0);
    
    if theta_sign_ok && dstar_sign_ok {
        println!("Signs match for θ and δ* updates.");
    } else {
        println!("*** SIGN MISMATCH DETECTED ***");
        if !theta_sign_ok {
            println!("  θ update: XFOIL={:+.2e}, Rust={:+.2e} - WRONG SIGN", 
                xfoil_solution[1], rust_solution[1]);
        }
        if !dstar_sign_ok {
            println!("  δ* update: XFOIL={:+.2e}, Rust={:+.2e} - WRONG SIGN", 
                xfoil_solution[2], rust_solution[2]);
        }
    }
}
