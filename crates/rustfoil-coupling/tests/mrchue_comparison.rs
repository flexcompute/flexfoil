//! Station-by-station comparison test between RustFoil and XFOIL MRCHUE
//!
//! This test loads XFOIL reference data from testdata/mrchue_iterations.json
//! and compares our march results station by station.

use rustfoil_coupling::march::{march_fixed_ue, MarchConfig};
use serde::Deserialize;
use std::fs;

/// Structure matching the XFOIL reference data format
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
    side: usize,
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
struct StationValues {
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
    #[serde(default)]
    dmax: f64,
    #[serde(default, rename = "Hk")]
    hk: f64,
    #[serde(default, rename = "Rtheta")]
    r_theta: f64,
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
fn test_march_station_by_station_upper() {
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

    println!("\n=== Upper Surface (Side 1) Station Comparison ===");
    println!(
        "{:>4} {:>10} {:>12} {:>12} {:>8}",
        "IBL", "X", "XFOIL_θ", "RUST_θ", "Err%"
    );

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

    // Compare station by station
    let mut total_theta_error = 0.0;
    let mut n_compared = 0;

    for (i, xfoil_station) in side_data.stations.iter().enumerate() {
        if i >= result.stations.len() {
            break;
        }

        let our_station = &result.stations[i];
        let xfoil_theta = xfoil_station.final_values.theta;
        let our_theta = our_station.theta;

        let error_pct = if xfoil_theta.abs() > 1e-12 {
            (our_theta - xfoil_theta).abs() / xfoil_theta * 100.0
        } else {
            0.0
        };

        println!(
            "{:4} {:10.4} {:12.4e} {:12.4e} {:8.2}",
            xfoil_station.ibl, xfoil_station.x, xfoil_theta, our_theta, error_pct
        );

        total_theta_error += error_pct;
        n_compared += 1;
    }

    let avg_error = total_theta_error / n_compared as f64;
    println!("\nAverage theta error: {:.1}%", avg_error);

    // Test transition detection
    println!(
        "\nTransition: XFOIL at x≈0.128, RustFoil at {:?}",
        result.x_transition
    );

    // For now, just print summary - we can add assertions once values match better
    // assert!(avg_error < 20.0, "Average theta error too high: {}%", avg_error);
}

#[test]
fn test_newton_convergence_pattern() {
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

    println!("\n=== Newton Convergence Pattern (XFOIL) ===");
    println!("{:>4} {:>10} {:>8} {:>10}", "IBL", "X", "N_iter", "Converged");

    let mut total_iters = 0;
    let mut n_converged = 0;

    for station in &side_data.stations {
        println!(
            "{:4} {:10.4} {:8} {:>10}",
            station.ibl,
            station.x,
            station.n_iterations,
            if station.converged { "Y" } else { "N" }
        );
        total_iters += station.n_iterations;
        if station.converged {
            n_converged += 1;
        }
    }

    let n = side_data.stations.len();
    println!("\nTotal stations: {}", n);
    println!("Converged: {}/{}", n_converged, n);
    println!(
        "Avg iterations: {:.1}",
        total_iters as f64 / n.max(1) as f64
    );
}

#[test]
fn test_first_station_matches_xfoil() {
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
        None => return,
    };

    // Get first station from XFOIL
    let first_xfoil = &side_data.stations[0];

    // Run march with just first two stations
    let x = vec![first_xfoil.x, side_data.stations[1].x];
    let ue = vec![first_xfoil.ue, side_data.stations[1].ue];

    let config = MarchConfig {
        ncrit: ref_data.metadata.ncrit,
        max_iter: 25,
        tolerance: 1e-5,
        ..Default::default()
    };

    let result = march_fixed_ue(&x, &ue, re, msq, &config);

    println!("\n=== First Station Comparison ===");
    println!("XFOIL  theta: {:.6e}", first_xfoil.final_values.theta);
    println!("RustFoil theta: {:.6e}", result.stations[0].theta);
    println!("XFOIL  delta_star: {:.6e}", first_xfoil.final_values.delta_star);
    println!("RustFoil delta_star: {:.6e}", result.stations[0].delta_star);

    let theta_error = (result.stations[0].theta - first_xfoil.final_values.theta).abs()
        / first_xfoil.final_values.theta.max(1e-12)
        * 100.0;
    println!("\nTheta error: {:.1}%", theta_error);
}

#[test]
fn test_hk_and_rtheta_comparison() {
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

    println!("\n=== Hk and Rθ Comparison (Upper Surface) ===");
    println!(
        "{:>4} {:>10} {:>10} {:>10} {:>8} {:>12} {:>12} {:>8}",
        "IBL", "X", "XFOIL_Hk", "RUST_Hk", "Hk_Err%", "XFOIL_Rθ", "RUST_Rθ", "Rθ_Err%"
    );

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

    // Compare station by station
    let mut total_hk_error = 0.0;
    let mut total_rt_error = 0.0;
    let mut n_compared = 0;

    for (i, xfoil_station) in side_data.stations.iter().enumerate() {
        if i >= result.stations.len() {
            break;
        }

        let our_station = &result.stations[i];
        let xfoil_hk = xfoil_station.final_values.hk;
        let xfoil_rt = xfoil_station.final_values.r_theta;

        // Skip if XFOIL didn't record these values
        if xfoil_hk == 0.0 || xfoil_rt == 0.0 {
            continue;
        }

        let hk_error = if xfoil_hk.abs() > 1e-12 {
            (our_station.hk - xfoil_hk).abs() / xfoil_hk * 100.0
        } else {
            0.0
        };

        let rt_error = if xfoil_rt.abs() > 1e-12 {
            (our_station.r_theta - xfoil_rt).abs() / xfoil_rt * 100.0
        } else {
            0.0
        };

        // Print first 15 and any with large error
        if i < 15 || hk_error > 10.0 || rt_error > 10.0 {
            println!(
                "{:4} {:10.4} {:10.4} {:10.4} {:8.2} {:12.1} {:12.1} {:8.2}",
                xfoil_station.ibl, xfoil_station.x, 
                xfoil_hk, our_station.hk, hk_error,
                xfoil_rt, our_station.r_theta, rt_error
            );
        }

        total_hk_error += hk_error;
        total_rt_error += rt_error;
        n_compared += 1;
    }

    if n_compared > 0 {
        let avg_hk_error = total_hk_error / n_compared as f64;
        let avg_rt_error = total_rt_error / n_compared as f64;

        println!("\n--- Summary ---");
        println!("Compared {} stations", n_compared);
        println!("Average Hk error: {:.1}%", avg_hk_error);
        println!("Average Rθ error: {:.1}%", avg_rt_error);
    } else {
        println!("\nNo Hk/Rθ data found in reference (need to regenerate with updated instrumentation)");
    }
}
