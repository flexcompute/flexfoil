//! Theta evolution comparison test
//!
//! Compares θ and δ* station-by-station against XFOIL to identify divergence.

use rustfoil_coupling::march::{march_fixed_ue, MarchConfig};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct MrchueReference {
    metadata: Metadata,
    sides: std::collections::HashMap<String, SideData>,
}

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    mach: f64,
    ncrit: f64,
}

#[derive(Debug, Deserialize)]
struct SideData {
    stations: Vec<StationRef>,
}

#[derive(Debug, Deserialize)]
struct StationRef {
    ibl: usize,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
    initial: InitialValues,
    #[serde(rename = "final")]
    final_values: FinalValues,
}

#[derive(Debug, Deserialize)]
struct InitialValues {
    theta: f64,
    delta_star: f64,
}

#[derive(Debug, Deserialize)]
struct FinalValues {
    theta: f64,
    delta_star: f64,
    #[serde(default)]
    ampl: Option<f64>,
    #[serde(rename = "Hk", default)]
    hk: Option<f64>,
    #[serde(rename = "Rtheta", default)]
    rtheta: Option<f64>,
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

#[test]
fn test_theta_evolution_station_by_station() {
    let ref_data = match load_reference() {
        Some(data) => data,
        None => {
            eprintln!("Skipping test: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data.metadata.reynolds;
    let msq = ref_data.metadata.mach * ref_data.metadata.mach;

    let side_data = match ref_data.sides.get("1") {
        Some(data) => data,
        None => {
            eprintln!("Skipping test: side 1 data not found");
            return;
        }
    };

    // Extract x and Ue arrays
    let x: Vec<f64> = side_data.stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = side_data.stations.iter().map(|s| s.ue).collect();

    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 0.1, // Match XFOIL's convergence criterion
        ..Default::default()
    };

    let result = march_fixed_ue(&x, &ue, re, msq, &config);

    // Print transition info
    println!("\n=== Transition Info ===");
    println!("RustFoil transition: x={:?}, station={:?}", 
             result.x_transition, result.transition_index);
    println!("Ncrit = {}", ref_data.metadata.ncrit);
    
    // Find XFOIL transition from ampl values
    for (i, station) in side_data.stations.iter().enumerate() {
        let ampl = station.final_values.ampl.unwrap_or(0.0);
        if ampl >= ref_data.metadata.ncrit {
            println!("XFOIL transition: x={:.6}, station={} (ampl={:.4})", 
                     station.x, i, ampl);
            break;
        }
    }
    
    // Analyze H growth pattern between XFOIL and us
    println!("\n=== H Growth Pattern Analysis ===");
    println!("{:>4} {:>10} {:>8} {:>8} {:>10} {:>8} {:>8} {:>10}",
             "Idx", "x", "H_xfoil", "H_rust", "H_err%", "dH_xf%", "dH_rs%", "dH_err");
    
    for i in 15..25.min(result.stations.len()) {
        let xfoil = &side_data.stations[i];
        let ours = &result.stations[i];
        
        let h_xfoil = xfoil.final_values.delta_star / xfoil.final_values.theta;
        let h_rust = ours.h;
        let h_err = (h_rust - h_xfoil) / h_xfoil * 100.0;
        
        // H growth from previous station
        let (dh_xf, dh_rs) = if i > 0 {
            let xf_prev = &side_data.stations[i-1];
            let h_xf_prev = xf_prev.final_values.delta_star / xf_prev.final_values.theta;
            let rs_prev = &result.stations[i-1];
            ((h_xfoil - h_xf_prev) / h_xf_prev * 100.0, (h_rust - rs_prev.h) / rs_prev.h * 100.0)
        } else {
            (0.0, 0.0)
        };
        
        let marker = if h_err.abs() > 1.0 { " <--" } else { "" };
        println!("{:4} {:10.4} {:8.4} {:8.4} {:+10.1} {:+8.2} {:+8.2} {:+10.2}{}",
                 i, xfoil.x, h_xfoil, h_rust, h_err, dh_xf, dh_rs, dh_rs - dh_xf, marker);
    }
    
    println!("\nKey insight: Our H growth per station is LOWER than XFOIL's starting at station 18.");
    println!("This causes H to diverge over time, making Rcrit higher and blocking amplification.");
    
    // Check what happens at station 18 specifically
    println!("\n=== Station 18 Analysis ===");
    let i = 18;
    let xfoil = &side_data.stations[i];
    let xfoil_prev = &side_data.stations[i-1];
    let ours = &result.stations[i];
    let ours_prev = &result.stations[i-1];
    
    // XFOIL inputs to Newton solve
    let x_new = xfoil.x;
    let ue_new = xfoil.ue;
    
    println!("Inputs (should match):");
    println!("  x[18] = {:.6}", x_new);
    println!("  Ue[18] = {:.6}", ue_new);
    println!("  prev θ: XFOIL={:.6e}, Rust={:.6e}", xfoil_prev.final_values.theta, ours_prev.theta);
    println!("  prev δ*: XFOIL={:.6e}, Rust={:.6e}", xfoil_prev.final_values.delta_star, ours_prev.delta_star);
    
    println!("\nOutputs:");
    println!("  θ: XFOIL={:.6e}, Rust={:.6e}", xfoil.final_values.theta, ours.theta);
    println!("  δ*: XFOIL={:.6e}, Rust={:.6e}", xfoil.final_values.delta_star, ours.delta_star);
    println!("  H: XFOIL={:.4}, Rust={:.4}", 
             xfoil.final_values.delta_star / xfoil.final_values.theta, ours.h);
    
    // What dH does each get?
    let dh_xfoil = xfoil.final_values.delta_star / xfoil.final_values.theta 
                   - xfoil_prev.final_values.delta_star / xfoil_prev.final_values.theta;
    let dh_rust = ours.h - ours_prev.h;
    println!("  ΔH: XFOIL={:+.6}, Rust={:+.6}", dh_xfoil, dh_rust);
    println!("  θ growth: XFOIL={:+.1}%, Rust={:+.1}%",
             (xfoil.final_values.theta / xfoil_prev.final_values.theta - 1.0) * 100.0,
             (ours.theta / ours_prev.theta - 1.0) * 100.0);
    println!("  δ* growth: XFOIL={:+.1}%, Rust={:+.1}%",
             (xfoil.final_values.delta_star / xfoil_prev.final_values.delta_star - 1.0) * 100.0,
             (ours.delta_star / ours_prev.delta_star - 1.0) * 100.0);

    // Print detailed comparison at key stations
    println!("\n=== Detailed BL State Comparison (stations 15-20) ===");
    println!("NOTE: Reference Hk/Rtheta values are INCORRECT (from different source)");
    println!("Our values should match Re*Ue*theta and H=dstar/theta\n");
    
    for i in 15..21.min(result.stations.len()) {
        let xfoil = &side_data.stations[i];
        let ours = &result.stations[i];
        let xfoil_ampl = xfoil.final_values.ampl.unwrap_or(0.0);
        
        // Compute correct H and Rtheta from XFOIL's theta/dstar
        let h_correct = xfoil.final_values.delta_star / xfoil.final_values.theta;
        let rtheta_correct = 1e6 * xfoil.ue * xfoil.final_values.theta;
        
        println!("Station {} (x={:.6}):", i, xfoil.x);
        println!("  theta:  XFOIL={:.6e}  Rust={:.6e}  err={:+.1}%", 
                 xfoil.final_values.theta, ours.theta,
                 (ours.theta - xfoil.final_values.theta) / xfoil.final_values.theta * 100.0);
        println!("  dstar:  XFOIL={:.6e}  Rust={:.6e}  err={:+.1}%",
                 xfoil.final_values.delta_star, ours.delta_star,
                 (ours.delta_star - xfoil.final_values.delta_star) / xfoil.final_values.delta_star * 100.0);
        println!("  H:      correct={:.4}       Rust={:.4}       err={:+.1}%",
                 h_correct, ours.h, (ours.h - h_correct) / h_correct * 100.0);
        println!("  Hk:     correct={:.4}       Rust={:.4}       (should be ~H at M=0)",
                 h_correct, ours.hk);
        println!("  Rtheta: correct={:.1}       Rust={:.1}       err={:+.1}%",
                 rtheta_correct, ours.r_theta, (ours.r_theta - rtheta_correct) / rtheta_correct * 100.0);
        println!("  ampl:   XFOIL={:.6e}  Rust={:.6e}",
                 xfoil_ampl, ours.ampl);
        println!();
    }

    println!("\n=== θ Evolution Comparison (Upper Surface) ===");
    println!("{:>4} {:>10} {:>12} {:>12} {:>8} {:>12} {:>12} {:>8}",
        "Idx", "x", "θ_xfoil", "θ_rust", "θ_err%", "δ*_xfoil", "δ*_rust", "δ*_err%");

    let mut first_divergence_idx: Option<usize> = None;
    let mut cumulative_theta_error = 0.0;

    for (i, (xfoil_station, our_station)) in side_data.stations.iter()
        .zip(result.stations.iter())
        .enumerate()
    {
        let xfoil_theta = xfoil_station.final_values.theta;
        let rust_theta = our_station.theta;
        let xfoil_dstar = xfoil_station.final_values.delta_star;
        let rust_dstar = our_station.delta_star;

        let theta_err = if xfoil_theta > 1e-12 {
            (rust_theta - xfoil_theta) / xfoil_theta * 100.0
        } else {
            0.0
        };

        let dstar_err = if xfoil_dstar > 1e-12 {
            (rust_dstar - xfoil_dstar) / xfoil_dstar * 100.0
        } else {
            0.0
        };

        cumulative_theta_error += theta_err.abs();

        // Mark first significant divergence (>10%)
        if first_divergence_idx.is_none() && theta_err.abs() > 10.0 {
            first_divergence_idx = Some(i);
        }

        // Print first 25 stations and any with large error
        if i < 25 || theta_err.abs() > 20.0 {
            let marker = if theta_err.abs() > 10.0 { " ***" } else { "" };
            println!(
                "{:4} {:10.4} {:12.6e} {:12.6e} {:+8.1} {:12.6e} {:12.6e} {:+8.1}{}",
                i, xfoil_station.x, 
                xfoil_theta, rust_theta, theta_err,
                xfoil_dstar, rust_dstar, dstar_err,
                marker
            );
        }
    }

    let n = side_data.stations.len().min(result.stations.len());
    let avg_error = cumulative_theta_error / n as f64;

    println!("\n--- Summary ---");
    println!("Stations compared: {}", n);
    println!("Average |θ error|: {:.1}%", avg_error);
    if let Some(idx) = first_divergence_idx {
        println!("First divergence (>10%): station {} at x={:.4}", 
            idx, side_data.stations[idx].x);
    }
}

#[test]
fn test_initial_theta_comparison() {
    let ref_data = match load_reference() {
        Some(data) => data,
        None => {
            eprintln!("Skipping test: mrchue_iterations.json not found");
            return;
        }
    };

    let side_data = match ref_data.sides.get("1") {
        Some(data) => data,
        None => return,
    };

    // Check initial values at first few stations
    println!("\n=== XFOIL Initial Values (First 5 Stations) ===");
    println!("{:>4} {:>10} {:>12} {:>12} {:>8}",
        "IBL", "x", "θ_initial", "δ*_initial", "H_init");

    for station in side_data.stations.iter().take(5) {
        let h_init = station.initial.delta_star / station.initial.theta;
        println!("{:4} {:10.4} {:12.6e} {:12.6e} {:8.4}",
            station.ibl, station.x,
            station.initial.theta, station.initial.delta_star, h_init);
    }

    // The first station's initial θ is crucial - it's the stagnation point value
    let first = &side_data.stations[0];
    println!("\n=== Stagnation Point Analysis ===");
    println!("First station IBL={}, x={:.6}", first.ibl, first.x);
    println!("  Initial θ = {:.6e}", first.initial.theta);
    println!("  Initial δ* = {:.6e}", first.initial.delta_star);
    println!("  Initial H = {:.4}", first.initial.delta_star / first.initial.theta);
    println!("  Final θ = {:.6e}", first.final_values.theta);
    println!("  Final δ* = {:.6e}", first.final_values.delta_star);
}
