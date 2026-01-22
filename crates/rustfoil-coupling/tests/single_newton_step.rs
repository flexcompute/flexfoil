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
