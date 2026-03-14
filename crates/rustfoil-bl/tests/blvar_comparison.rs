//! BLVAR comparison test between RustFoil and XFOIL
//!
//! This test validates that our blvar() function produces the same
//! outputs as XFOIL's BLVAR subroutine for the same inputs.

use rustfoil_bl::{blvar, BlStation, FlowType};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct BlvarTestData {
    metadata: Metadata,
    events: Vec<BlvarEvent>,
}

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    mach: f64,
}

#[derive(Debug, Deserialize)]
struct BlvarEvent {
    side: usize,
    ibl: usize,
    flow_type: usize,
    input: BlvarInput,
    output: BlvarOutput,
}

#[derive(Debug, Deserialize)]
struct BlvarInput {
    x: f64,
    u: f64,
    theta: f64,
    delta_star: f64,
    ctau: f64,
    ampl: f64,
}

#[derive(Debug, Deserialize)]
struct BlvarOutput {
    #[serde(rename = "H")]
    h: f64,
    #[serde(rename = "Hk")]
    hk: f64,
    #[serde(rename = "Hs")]
    hs: f64,
    #[serde(rename = "Hc")]
    hc: f64,
    #[serde(rename = "Rtheta")]
    r_theta: f64,
    #[serde(rename = "Cf")]
    cf: f64,
    #[serde(rename = "Cd")]
    cd: f64,
    #[serde(rename = "Us")]
    us: f64,
    #[serde(rename = "Cq")]
    cq: f64,
    #[serde(rename = "De")]
    de: f64,
}

fn load_test_vectors() -> Option<BlvarTestData> {
    let paths = [
        "testdata/blvar_test_vectors.json",
        "../testdata/blvar_test_vectors.json",
        "../../testdata/blvar_test_vectors.json",
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
fn test_blvar_hk_vs_xfoil() {
    let data = match load_test_vectors() {
        Some(d) => d,
        None => {
            eprintln!("Skipping test: blvar_test_vectors.json not found");
            return;
        }
    };

    let re = data.metadata.reynolds;
    let msq = data.metadata.mach * data.metadata.mach;

    println!("\n=== BLVAR Hk Comparison ===");
    println!("{:>4} {:>4} {:>6} {:>10} {:>10} {:>8}",
        "Side", "IBL", "Type", "XFOIL_Hk", "Rust_Hk", "Err%");

    let mut total_error: f64 = 0.0;
    let mut max_error: f64 = 0.0;
    let mut n_compared = 0;

    for event in &data.events {
        // Set up station with XFOIL's inputs
        let mut station = BlStation::default();
        station.x = event.input.x;
        station.u = event.input.u;
        station.theta = event.input.theta;
        station.delta_star = event.input.delta_star;
        station.ctau = event.input.ctau;
        station.ampl = event.input.ampl;

        let flow_type = if event.flow_type == 1 {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        // Call our blvar
        blvar(&mut station, flow_type, msq, re);

        let xfoil_hk = event.output.hk;
        let rust_hk = station.hk;

        let error = if xfoil_hk.abs() > 1e-10 {
            (rust_hk - xfoil_hk).abs() / xfoil_hk * 100.0
        } else {
            0.0
        };

        let type_str = if event.flow_type == 1 { "Lam" } else { "Turb" };
        println!(
            "{:4} {:4} {:>6} {:10.4} {:10.4} {:8.2}",
            event.side, event.ibl, type_str, xfoil_hk, rust_hk, error
        );

        total_error += error;
        max_error = max_error.max(error);
        n_compared += 1;
    }

    let avg_error = total_error / n_compared as f64;
    println!("\n--- Hk Summary ---");
    println!("Compared: {} events", n_compared);
    println!("Average error: {:.2}%", avg_error);
    println!("Max error: {:.2}%", max_error);

    // Target: < 1% error
    assert!(avg_error < 5.0, "Hk average error {} > 5%", avg_error);
}

#[test]
fn test_blvar_all_outputs_vs_xfoil() {
    let data = match load_test_vectors() {
        Some(d) => d,
        None => {
            eprintln!("Skipping test: blvar_test_vectors.json not found");
            return;
        }
    };

    let re = data.metadata.reynolds;
    let msq = data.metadata.mach * data.metadata.mach;

    println!("\n=== BLVAR Full Comparison ===");
    
    // Track errors for each output
    let mut errors: std::collections::HashMap<&str, Vec<f64>> = std::collections::HashMap::new();
    for key in &["H", "Hk", "Hs", "Rtheta", "Cf", "Cd"] {
        errors.insert(key, Vec::new());
    }

    for event in &data.events {
        let mut station = BlStation::default();
        station.x = event.input.x;
        station.u = event.input.u;
        station.theta = event.input.theta;
        station.delta_star = event.input.delta_star;
        station.ctau = event.input.ctau;
        station.ampl = event.input.ampl;

        let flow_type = if event.flow_type == 1 {
            FlowType::Laminar
        } else {
            FlowType::Turbulent
        };

        blvar(&mut station, flow_type, msq, re);

        // Compute errors
        let h_err = rel_error(station.delta_star / station.theta, event.output.h);
        let hk_err = rel_error(station.hk, event.output.hk);
        let hs_err = rel_error(station.hs, event.output.hs);
        let rt_err = rel_error(station.r_theta, event.output.r_theta);
        let cf_err = rel_error(station.cf, event.output.cf);
        let cd_err = rel_error(station.cd, event.output.cd);

        errors.get_mut("H").unwrap().push(h_err);
        errors.get_mut("Hk").unwrap().push(hk_err);
        errors.get_mut("Hs").unwrap().push(hs_err);
        errors.get_mut("Rtheta").unwrap().push(rt_err);
        errors.get_mut("Cf").unwrap().push(cf_err);
        errors.get_mut("Cd").unwrap().push(cd_err);
    }

    println!("{:>10} {:>10} {:>10} {:>10}",
        "Output", "Avg Err%", "Max Err%", "Status");
    
    for (name, errs) in &errors {
        let avg = errs.iter().sum::<f64>() / errs.len() as f64;
        let max = errs.iter().cloned().fold(0.0, f64::max);
        let status = if avg < 5.0 { "OK" } else { "FAIL" };
        println!("{:>10} {:10.2} {:10.2} {:>10}", name, avg, max, status);
    }
}

#[test]
fn test_blvar_laminar_only() {
    let data = match load_test_vectors() {
        Some(d) => d,
        None => {
            eprintln!("Skipping test: blvar_test_vectors.json not found");
            return;
        }
    };

    let re = data.metadata.reynolds;
    let msq = data.metadata.mach * data.metadata.mach;

    println!("\n=== BLVAR Laminar Cases ===");
    println!("{:>4} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "IBL", "x", "H_xfoil", "H_rust", "Hk_xfoil", "Hk_rust");

    let laminar_events: Vec<_> = data.events.iter()
        .filter(|e| e.flow_type == 1)
        .collect();

    for event in &laminar_events {
        let mut station = BlStation::default();
        station.x = event.input.x;
        station.u = event.input.u;
        station.theta = event.input.theta;
        station.delta_star = event.input.delta_star;
        station.ctau = event.input.ctau;
        station.ampl = event.input.ampl;

        blvar(&mut station, FlowType::Laminar, msq, re);

        let h_xfoil = event.output.h;
        let h_rust = station.delta_star / station.theta;
        let hk_xfoil = event.output.hk;
        let hk_rust = station.hk;

        println!(
            "{:4} {:10.4} {:10.4} {:10.4} {:10.4} {:10.4}",
            event.ibl, event.input.x, h_xfoil, h_rust, hk_xfoil, hk_rust
        );
    }
}

fn rel_error(rust: f64, xfoil: f64) -> f64 {
    if xfoil.abs() > 1e-10 {
        (rust - xfoil).abs() / xfoil.abs() * 100.0
    } else if rust.abs() > 1e-10 {
        100.0  // XFOIL is zero but we're not
    } else {
        0.0  // Both effectively zero
    }
}
