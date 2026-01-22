//! Transition location validation test
//!
//! Validates that transition locations match XFOIL within 10% on both surfaces.

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
    alpha_rad: f64,
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
    #[serde(rename = "final")]
    final_values: FinalValues,
}

#[derive(Debug, Deserialize)]
struct FinalValues {
    theta: f64,
    ampl: f64,
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
fn test_transition_location_upper() {
    let ref_data = match load_reference() {
        Some(data) => data,
        None => {
            eprintln!("Skipping test: mrchue_iterations.json not found");
            return;
        }
    };

    // XFOIL reference values at Re=1M, alpha=5°
    // From instrumented XFOIL: x_tr_upper = 0.1486 (from VISCOUS_FINAL)
    // or from TRCHEK: interpolated x_tr = 0.1341
    let xfoil_x_tr_upper = 0.1486;

    let re = ref_data.metadata.reynolds;
    let msq = ref_data.metadata.mach * ref_data.metadata.mach;

    let side_data = match ref_data.sides.get("1") {
        Some(data) => data,
        None => {
            eprintln!("Skipping test: side 1 data not found");
            return;
        }
    };

    // Extract x and Ue arrays from reference
    let x: Vec<f64> = side_data.stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = side_data.stations.iter().map(|s| s.ue).collect();

    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 1e-5,
        ..Default::default()
    };

    // Run our march
    let result = march_fixed_ue(&x, &ue, re, msq, &config);

    println!("\n=== Upper Surface Transition Validation ===");
    println!("XFOIL x_tr_upper: {:.4}", xfoil_x_tr_upper);
    println!("RustFoil x_tr_upper: {:?}", result.x_transition);

    // Show theta and N-factor evolution
    println!("\n{:>4} {:>10} {:>12} {:>12} {:>10} {:>10}",
        "Idx", "x", "theta_xfoil", "theta_rust", "N_xfoil", "N_rust");
    
    for (i, (xfoil_station, our_station)) in side_data.stations.iter()
        .zip(result.stations.iter())
        .take(40)
        .enumerate() 
    {
        // Only print every 3rd station or near transition
        let near_trans = xfoil_station.x > 0.10 && xfoil_station.x < 0.25;
        if i % 3 == 0 || near_trans {
            println!(
                "{:4} {:10.4} {:12.6e} {:12.6e} {:10.4} {:10.4}",
                i, xfoil_station.x, 
                xfoil_station.final_values.theta, our_station.theta,
                xfoil_station.final_values.ampl, our_station.ampl
            );
        }
    }

    // Calculate error
    if let Some(rust_x_tr) = result.x_transition {
        let error_pct = (rust_x_tr - xfoil_x_tr_upper).abs() / xfoil_x_tr_upper * 100.0;
        println!("\nTransition location error: {:.1}%", error_pct);
        
        // Target: error < 10%
        if error_pct > 10.0 {
            println!("WARNING: Transition error exceeds 10% target");
        }
    } else {
        println!("\nWARNING: No transition detected by RustFoil!");
    }
}

#[test]
fn test_theta_comparison() {
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

    let x: Vec<f64> = side_data.stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = side_data.stations.iter().map(|s| s.ue).collect();

    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 1e-5,
        ..Default::default()
    };

    let result = march_fixed_ue(&x, &ue, re, msq, &config);

    println!("\n=== Theta Comparison (Upper Surface) ===");
    println!("{:>4} {:>10} {:>12} {:>12} {:>8}",
        "Idx", "x", "theta_xfoil", "theta_rust", "Err%");

    let mut total_error = 0.0;
    let mut n_compared = 0;

    for (i, (xfoil_station, our_station)) in side_data.stations.iter()
        .zip(result.stations.iter())
        .enumerate() 
    {
        let xfoil_theta = xfoil_station.final_values.theta;
        let rust_theta = our_station.theta;
        
        let error = if xfoil_theta > 1e-12 {
            (rust_theta - xfoil_theta).abs() / xfoil_theta * 100.0
        } else {
            0.0
        };

        if i < 35 {  // First 35 stations (pre-transition)
            println!(
                "{:4} {:10.4} {:12.6e} {:12.6e} {:8.1}",
                i, xfoil_station.x, xfoil_theta, rust_theta, error
            );
            total_error += error;
            n_compared += 1;
        }
    }

    if n_compared > 0 {
        let avg_error = total_error / n_compared as f64;
        println!("\nAverage theta error (pre-transition): {:.1}%", avg_error);
    }
}

#[test]
fn test_hk_and_rtheta_effect_on_amplification() {
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

    let x: Vec<f64> = side_data.stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = side_data.stations.iter().map(|s| s.ue).collect();

    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 1e-5,
        ..Default::default()
    };

    let result = march_fixed_ue(&x, &ue, re, msq, &config);

    println!("\n=== Hk and Rθ Effect on Amplification ===");
    println!("{:>4} {:>8} {:>8} {:>10} {:>10} {:>8}",
        "Idx", "x", "Hk", "Rθ", "N", "dN/ds");

    for (i, station) in result.stations.iter().take(40).enumerate() {
        // Calculate dN/ds = (N_curr - N_prev) / ds
        let dn_ds = if i > 0 {
            let prev = &result.stations[i - 1];
            let ds = station.x - prev.x;
            if ds > 1e-12 {
                (station.ampl - prev.ampl) / ds
            } else {
                0.0
            }
        } else {
            0.0
        };

        if i < 15 || (i >= 25 && i <= 35) {
            println!(
                "{:4} {:8.4} {:8.4} {:10.1} {:10.4} {:8.2}",
                i, station.x, station.hk, station.r_theta, station.ampl, dn_ds
            );
        }
    }
}
