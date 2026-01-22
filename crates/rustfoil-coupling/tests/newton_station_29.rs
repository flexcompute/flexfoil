//! Newton solver test at divergence point (IBL 31, x=0.1276)
//!
//! This test isolates the Newton solve at the exact station where theta diverges
//! from XFOIL in the march_fixed_ue function.
//!
//! Key finding: At IBL 31, Hk jumps to 5.39 which exceeds hlmax=3.8.
//! The Hk clamping logic should handle this, but may be causing the divergence.

use rustfoil_bl::equations::{bldif, blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::solve::{build_4x4_system, solve_4x4};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    mach: f64,
    ncrit: f64,
}

#[derive(Debug, Deserialize)]
struct StationValues {
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
    #[serde(default, rename = "Hk")]
    hk: f64,
    #[serde(default, rename = "Rtheta")]
    r_theta: f64,
}

#[derive(Debug, Deserialize)]
struct StationRef {
    ibl: usize,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
    n_iterations: usize,
    converged: bool,
    initial: StationValues,
    #[serde(rename = "final")]
    final_values: StationValues,
}

#[derive(Debug, Deserialize)]
struct SideData {
    stations: Vec<StationRef>,
}

#[derive(Debug, Deserialize)]
struct MrchueReference {
    metadata: Metadata,
    sides: std::collections::HashMap<String, SideData>,
}

fn load_reference() -> Option<MrchueReference> {
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

/// Test Newton solve at IBL 31 (index 29) where theta diverges
/// 
/// XFOIL data at this station:
/// - x = 0.1276
/// - Hk_final = 5.39 (exceeds hlmax=3.8)
/// - theta_final = 2.218e-4
/// - ampl_final = 7.88
#[test]
fn test_newton_at_divergence_point() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data.metadata.reynolds;
    let msq = ref_data.metadata.mach * ref_data.metadata.mach;
    let stations = &ref_data.sides.get("1").unwrap().stations;
    
    // IBL 30 is at index 28, IBL 31 is at index 29
    let prev_ref = &stations[28];  // IBL 30
    let curr_ref = &stations[29];  // IBL 31 - divergence point
    
    println!("\n=== Newton Solve at Divergence Point (IBL {}) ===", curr_ref.ibl);
    println!("x = {:.4}, Ue = {:.4}", curr_ref.x, curr_ref.ue);
    println!("\nXFOIL reference:");
    println!("  prev (IBL {}): theta={:.6e}, delta_star={:.6e}, Hk={:.4}", 
             prev_ref.ibl, prev_ref.final_values.theta, prev_ref.final_values.delta_star, prev_ref.final_values.hk);
    println!("  curr (IBL {}): theta={:.6e}, delta_star={:.6e}, Hk={:.4}",
             curr_ref.ibl, curr_ref.final_values.theta, curr_ref.final_values.delta_star, curr_ref.final_values.hk);
    
    // Build previous station with XFOIL's converged values
    let mut prev = BlStation::new();
    prev.x = prev_ref.x;
    prev.u = prev_ref.ue;
    prev.theta = prev_ref.final_values.theta;
    prev.delta_star = prev_ref.final_values.delta_star;
    prev.ampl = prev_ref.final_values.ampl;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);
    
    println!("\nPrev station after blvar: Hk={:.4}, Rt={:.2}", prev.hk, prev.r_theta);
    
    // Build current station with XFOIL's initial guess
    let mut curr = BlStation::new();
    curr.x = curr_ref.x;
    curr.u = curr_ref.ue;
    curr.theta = curr_ref.initial.theta;
    curr.delta_star = curr_ref.initial.delta_star;
    curr.ampl = curr_ref.initial.ampl;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);
    
    println!("Curr initial: theta={:.6e}, delta_star={:.6e}, Hk={:.4}", 
             curr.theta, curr.delta_star, curr.hk);
    
    // XFOIL's Hk limits
    let hlmax = 3.8;
    let htmax = 2.5;
    let tolerance = 0.1;
    let max_iter = 25;
    
    println!("\n=== Newton Iterations ===");
    println!("{:>4} {:>12} {:>12} {:>8} {:>8} {:>8} {:>10}",
             "Iter", "theta", "delta_star", "Hk", "dmax", "rlx", "converged");
    
    for iter in 0..max_iter {
        // Compute residuals and Jacobian
        let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);
        
        if iter == 0 {
            println!("\nDebug iteration 0:");
            println!("  Residuals: res_third={:.4e}, res_mom={:.4e}, res_shape={:.4e}",
                     res.res_third, res.res_mom, res.res_shape);
            println!("  VS2[0]: [{:.4e}, {:.4e}, {:.4e}, {:.4e}]",
                     jac.vs2[0][0], jac.vs2[0][1], jac.vs2[0][2], jac.vs2[0][3]);
            println!("  VS2[1]: [{:.4e}, {:.4e}, {:.4e}, {:.4e}]",
                     jac.vs2[1][0], jac.vs2[1][1], jac.vs2[1][2], jac.vs2[1][3]);
            println!("  VS2[2]: [{:.4e}, {:.4e}, {:.4e}, {:.4e}]",
                     jac.vs2[2][0], jac.vs2[2][1], jac.vs2[2][2], jac.vs2[2][3]);
        }
        
        // Check if Hk exceeds limit (approaching separation)
        let hklim = hlmax;
        let is_direct = curr.hk < hklim;
        
        // Build 4x4 system
        let hk_t = curr.hk / curr.theta.max(1e-12);
        let hk_d = -curr.hk / curr.delta_star.max(1e-12);
        
        let (a, b) = build_4x4_system(
            &jac.vs2,
            &[res.res_third, res.res_mom, res.res_shape],
            is_direct,
            hk_t,
            hk_d,
            0.0,
            hklim,
            curr.hk,
        );
        
        if iter == 0 {
            println!("  4x4 A diagonal: [{:.4e}, {:.4e}, {:.4e}, {:.4e}]",
                     a[0][0], a[1][1], a[2][2], a[3][3]);
            println!("  4x4 b: [{:.4e}, {:.4e}, {:.4e}, {:.4e}]", b[0], b[1], b[2], b[3]);
        }
        
        // Solve
        // Variable order: [dS/dAmpl, dθ, dδ*, dUe]
        let deltas = solve_4x4(&a, &b);
        
        if iter == 0 {
            println!("  Deltas: [dAmpl={:.4e}, dθ={:.4e}, dδ*={:.4e}, dUe={:.4e}]", 
                     deltas[0], deltas[1], deltas[2], deltas[3]);
            println!("  d_theta = {:.4e}, curr.theta = {:.4e}", deltas[1], curr.theta);
            println!("  Ratio d_theta/theta = {:.2}", deltas[1]/curr.theta);
        }
        
        // Find max delta for relaxation
        // Variable order from build_4x4_system: [dS/dAmpl, dθ, dδ*, dUe]
        let d_ampl = deltas[0];
        let d_theta = deltas[1];
        let d_delta = deltas[2];
        let d_ue = deltas[3];
        let _ = d_ampl; // Unused in laminar direct mode
        let _ = d_ue;   // Zero in direct mode
        
        let dmax = d_theta.abs().max(d_delta.abs());
        
        // Relaxation
        let rlx = if dmax > 3.0 { 3.0 / dmax } else { 1.0 };
        
        // Apply updates
        let theta_new = curr.theta + rlx * d_theta;
        let delta_new = curr.delta_star + rlx * d_delta;
        
        // Hk clamping (matches march.rs logic)
        let dsw = delta_new;
        let dslim = hklim * theta_new;
        let delta_clamped = if dsw > dslim { dslim } else { dsw };
        
        // Store old values for comparison
        let theta_old = curr.theta;
        let delta_old = curr.delta_star;
        let hk_old = curr.hk;
        
        curr.theta = theta_new;
        curr.delta_star = delta_clamped;
        
        // Recompute secondary variables
        blvar(&mut curr, FlowType::Laminar, msq, re);
        
        let converged = dmax <= tolerance;
        
        println!("{:4} {:12.6e} {:12.6e} {:8.4} {:8.4} {:8.4} {:>10}",
                 iter, curr.theta, curr.delta_star, curr.hk, dmax, rlx,
                 if converged { "YES" } else { "" });
        
        if delta_clamped != dsw {
            println!("      ^ Hk clamped: dsw={:.6e} -> dslim={:.6e}", dsw, dslim);
        }
        
        if converged {
            break;
        }
    }
    
    println!("\n=== Final Comparison ===");
    println!("                {:>12} {:>12} {:>8}", "theta", "delta_star", "Hk");
    println!("XFOIL:          {:12.6e} {:12.6e} {:8.4}", 
             curr_ref.final_values.theta, curr_ref.final_values.delta_star, curr_ref.final_values.hk);
    println!("RustFoil:       {:12.6e} {:12.6e} {:8.4}", 
             curr.theta, curr.delta_star, curr.hk);
    
    let theta_error = (curr.theta - curr_ref.final_values.theta).abs() / curr_ref.final_values.theta * 100.0;
    let delta_error = (curr.delta_star - curr_ref.final_values.delta_star).abs() / curr_ref.final_values.delta_star * 100.0;
    
    println!("\nErrors: theta={:.1}%, delta_star={:.1}%", theta_error, delta_error);
    
    // XFOIL's Hk is 5.39, which exceeds hlmax=3.8
    // This means XFOIL is NOT clamping Hk at this station
    // But our code might be clamping it, causing the divergence
    if curr_ref.final_values.hk > hlmax {
        println!("\nNOTE: XFOIL's Hk ({:.2}) exceeds hlmax ({:.2})", 
                 curr_ref.final_values.hk, hlmax);
        println!("This suggests XFOIL is NOT using Hk clamping at this station.");
    }
}

/// Run march_fixed_ue and trace what happens at each station
#[test]
fn test_march_with_trace() {
    use rustfoil_coupling::march::{march_fixed_ue, MarchConfig};
    
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data.metadata.reynolds;
    let msq = ref_data.metadata.mach * ref_data.metadata.mach;
    let stations = &ref_data.sides.get("1").unwrap().stations;
    
    let x: Vec<f64> = stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = stations.iter().map(|s| s.ue).collect();
    
    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 1e-5,
        debug_trace: true,  // Enable tracing
        ..Default::default()
    };
    
    let result = march_fixed_ue(&x, &ue, re, msq, &config);
    
    println!("\n=== March Trace Around Divergence Point ===");
    println!("{:>4} {:>8} {:>12} {:>12} {:>8} {:>8} {:>8} {:>8}",
             "Idx", "x", "theta_xfoil", "theta_rust", "Hk_xfoil", "Hk_rust", "Ue_xfoil", "Ue_rust");
    println!("{}", "-".repeat(88));
    
    for (i, (xfoil_st, rust_st)) in stations.iter().zip(result.stations.iter()).enumerate() {
        if i >= 25 && i <= 35 {  // Around divergence point
            let theta_err = (rust_st.theta - xfoil_st.final_values.theta).abs() 
                          / xfoil_st.final_values.theta * 100.0;
            
            let marker = if theta_err > 10.0 { " <-- DIVERGE" } else { "" };
            
            println!("{:4} {:8.4} {:12.6e} {:12.6e} {:8.4} {:8.4} {:8.4} {:8.4}{}",
                     i, xfoil_st.x, 
                     xfoil_st.final_values.theta, rust_st.theta,
                     xfoil_st.final_values.hk, rust_st.hk,
                     xfoil_st.ue, rust_st.u,
                     marker);
        }
    }
    
    println!("\nTransition:");
    println!("  XFOIL: IBL 32 at x ≈ 0.1408");
    if let Some(x_tr) = result.x_transition {
        println!("  RustFoil: x = {:.4}", x_tr);
    } else {
        println!("  RustFoil: No transition detected");
    }
}
