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
        tolerance: 1e-5,
        ..Default::default()
    };

    let result = march_fixed_ue(&x, &ue, re, msq, &config);

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
