//! Trace a single Newton step (station 1) to debug

use rustfoil_bl::equations::{blvar, bldif, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::solve::{build_4x4_system, solve_4x4};
use serde::Deserialize;
use std::fs;

fn load_reference() -> Option<serde_json::Value> {
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

/// Extract VSREZ_AFTER event for specific station/iteration
fn get_xfoil_solution(trace: &serde_json::Value, side: i64, ibl: i64, iter: i64) -> Option<[f64; 4]> {
    let events = trace["events"].as_array()?;
    for e in events {
        if e["subroutine"].as_str() == Some("VSREZ_AFTER")
            && e["side"].as_i64() == Some(side)
            && e["ibl"].as_i64() == Some(ibl)
            && e["newton_iter"].as_i64() == Some(iter)
        {
            let sol = e["VSREZ_solution"].as_array()?;
            let mut v = [0.0; 4];
            for (i, val) in sol.iter().enumerate() {
                if i < 4 {
                    v[i] = val.as_f64().unwrap_or(0.0);
                }
            }
            return Some(v);
        }
    }
    None
}

#[test]
fn test_single_newton_step() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;

    let stations = &ref_data["sides"]["1"]["stations"];
    let s0 = &stations[0];
    let s1 = &stations[1];

    // Build prev station (IBL=2 converged = XFOIL station 0 final)
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    // Build curr station (IBL=3 initial = use prev values as initial guess)
    let x_new = s1["x"].as_f64().unwrap();
    let ue_new = s1["Ue"].as_f64().unwrap();
    
    let mut curr = BlStation::new();
    curr.x = x_new;
    curr.u = ue_new;
    curr.theta = prev.theta;  // Initial guess
    curr.delta_star = prev.delta_star;
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== Single Newton Step (Station 0 → Station 1) ===");
    println!("\nPrev (IBL=2 final): x={:.6}, θ={:.6e}, δ*={:.6e}", prev.x, prev.theta, prev.delta_star);
    println!("Curr (IBL=3 init):  x={:.6}, θ={:.6e}, δ*={:.6e}", curr.x, curr.theta, curr.delta_star);
    println!("\nXFOIL target: θ={:.6e}", s1["final"]["theta"].as_f64().unwrap());

    // Compute residuals and Jacobian
    let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);
    
    println!("\n=== Residuals ===");
    println!("res_third = {:.6e} (XFOIL ~0 for laminar)", res.res_third);
    println!("res_mom = {:.6e} (XFOIL IBL=3 iter 1: -0.130)", res.res_mom);
    println!("res_shape = {:.6e}", res.res_shape);

    println!("\n=== Jacobian VS2 ===");
    for (i, row) in jac.vs2.iter().enumerate() {
        println!("  Row {}: [{:.4e}, {:.4e}, {:.4e}, {:.4e}, {:.4e}]",
            i, row[0], row[1], row[2], row[3], row[4]);
    }

    // Build 4x4 system
    let hlmax = 3.8;
    let (a, b) = build_4x4_system(
        &jac.vs2,
        &[res.res_third, res.res_mom, res.res_shape],
        true, // direct mode
        curr.hk / curr.theta,  // hk2_t
        -curr.hk / curr.delta_star, // hk2_d
        0.0, // hk2_u
        hlmax,
        curr.hk,
    );

    println!("\n=== 4x4 System ===");
    for (i, row) in a.iter().enumerate() {
        println!("  A[{}] = [{:.4e}, {:.4e}, {:.4e}, {:.4e}]",
            i, row[0], row[1], row[2], row[3]);
    }
    println!("  b = [{:.4e}, {:.4e}, {:.4e}, {:.4e}]", b[0], b[1], b[2], b[3]);

    // Solve
    let dx = solve_4x4(&a, &b);
    println!("\n=== Solution ===");
    println!("dx = [{:.4e}, {:.4e}, {:.4e}, {:.4e}]", dx[0], dx[1], dx[2], dx[3]);
    println!("(dx[1] = dθ, dx[2] = dδ*)");

    // Compute updates
    let dmax = (dx[1] / curr.theta.max(1e-12)).abs()
        .max((dx[2] / curr.delta_star.max(1e-12)).abs())
        .max((dx[0] / 10.0).abs());
    let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };
    
    println!("\n=== Update ===");
    println!("dmax = {:.4e}", dmax);
    println!("rlx = {:.4}", rlx);
    println!("θ_new = θ + rlx*dθ = {:.6e} + {:.4}*{:.6e} = {:.6e}",
        curr.theta, rlx, dx[1], curr.theta + rlx * dx[1]);
    println!("δ*_new = δ* + rlx*dδ* = {:.6e} + {:.4}*{:.6e} = {:.6e}",
        curr.delta_star, rlx, dx[2], curr.delta_star + rlx * dx[2]);

    // XFOIL's first iteration result
    let theta_after_iter1 = s1["initial"]["theta"].as_f64().unwrap();  // This is after IBL=2 converged
    // For IBL=3 iter 1: input θ = 4.121e-5, output θ = 3.622e-5
    println!("\nXFOIL IBL=3 iter 1: θ changes from 4.121e-5 to 3.622e-5 (dθ = -5.0e-6)");
    println!("Our Newton: θ changes from {:.6e} to {:.6e} (dθ = {:.6e})",
        curr.theta, curr.theta + rlx * dx[1], rlx * dx[1]);
}

#[test]
fn test_compare_with_xfoil_trace() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let trace = match load_debug_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let re = ref_data["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;

    let stations = &ref_data["sides"]["1"]["stations"];
    let s0 = &stations[0];
    let s1 = &stations[1];

    // Build prev station (IBL=2 converged)
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    // Build curr station (IBL=3 initial)
    let mut curr = BlStation::new();
    curr.x = s1["x"].as_f64().unwrap();
    curr.u = s1["Ue"].as_f64().unwrap();
    curr.theta = s1["initial"]["theta"].as_f64().unwrap();
    curr.delta_star = s1["initial"]["delta_star"].as_f64().unwrap();
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== Compare RustFoil vs XFOIL Newton Solution (IBL=3, iter=1) ===\n");

    // Get XFOIL's solution from debug trace
    let xfoil_dx = match get_xfoil_solution(&trace, 1, 3, 1) {
        Some(dx) => dx,
        None => {
            eprintln!("Could not find XFOIL solution for side=1, ibl=3, iter=1");
            return;
        }
    };

    // Compute RustFoil's solution
    let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);
    let (a, b) = build_4x4_system(
        &jac.vs2,
        &[res.res_third, res.res_mom, res.res_shape],
        true,
        curr.hk / curr.theta,
        -curr.hk / curr.delta_star,
        0.0,
        4.0,
        curr.hk,
    );
    let rust_dx = solve_4x4(&a, &b);

    println!("XFOIL dx: [{:+.6e}, {:+.6e}, {:+.6e}, {:+.6e}]", 
        xfoil_dx[0], xfoil_dx[1], xfoil_dx[2], xfoil_dx[3]);
    println!("Rust  dx: [{:+.6e}, {:+.6e}, {:+.6e}, {:+.6e}]", 
        rust_dx[0], rust_dx[1], rust_dx[2], rust_dx[3]);

    println!("\n--- Solution Comparison ---");
    let names = ["d_ctau/ampl", "d_theta", "d_dstar", "d_ue"];
    let mut sign_errors = 0;
    for i in 0..4 {
        let xfoil_val = xfoil_dx[i];
        let rust_val = rust_dx[i];
        let sign_ok = (xfoil_val >= 0.0) == (rust_val >= 0.0);
        if !sign_ok && xfoil_val.abs() > 1e-15 {
            sign_errors += 1;
        }
        let rel_err = if xfoil_val.abs() > 1e-30 {
            ((rust_val - xfoil_val) / xfoil_val).abs() * 100.0
        } else {
            0.0
        };
        println!("  {:12}: XFOIL={:+.4e}, Rust={:+.4e}, err={:5.1}%{}", 
            names[i], xfoil_val, rust_val, rel_err,
            if sign_ok { "" } else { " !! SIGN" });
    }

    if sign_errors > 0 {
        println!("\n*** SIGN MISMATCH: {} variable(s) have wrong sign! ***", sign_errors);
        println!("This indicates a bug in the Jacobian computation.");
    }
}
