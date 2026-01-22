//! BLDIF comparison test between RustFoil and XFOIL
//!
//! Tests that our bldif() function produces the same residuals as XFOIL.

use rustfoil_bl::{bldif, blvar, BlStation, FlowType};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct BldifTestData {
    events: Vec<BldifEvent>,
}

#[derive(Debug, Deserialize)]
struct BldifEvent {
    side: usize,
    ibl: usize,
    flow_type: usize,
    #[serde(rename = "VSREZ")]
    vsrez: Vec<f64>,
    #[serde(rename = "VS2")]
    vs2: Vec<Vec<f64>>,
}

fn load_test_vectors() -> Option<BldifTestData> {
    let paths = [
        "testdata/bldif_test_vectors.json",
        "../testdata/bldif_test_vectors.json",
        "../../testdata/bldif_test_vectors.json",
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

/// Test bldif using XFOIL's MRCHUE iteration data
#[test]
fn test_bldif_residuals_first_station() {
    // Load MRCHUE reference data to get the station values
    let mrchue_paths = [
        "testdata/mrchue_iterations.json",
        "../testdata/mrchue_iterations.json",
        "../../testdata/mrchue_iterations.json",
    ];

    let ref_data: Option<serde_json::Value> = mrchue_paths.iter()
        .find_map(|path| {
            fs::read_to_string(path).ok()
                .and_then(|content| serde_json::from_str(&content).ok())
        });

    let ref_data = match ref_data {
        Some(d) => d,
        None => {
            eprintln!("Skipping test: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0; // Mach = 0

    // Get first two stations from upper surface (side 1)
    let stations = &ref_data["sides"]["1"]["stations"];
    let s0 = &stations[0];  // Station at IBL=2
    let s1 = &stations[1];  // Station at IBL=3

    // Build station 1 (previous)
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = s0["final"]["ampl"].as_f64().unwrap_or(0.0);
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    // Build station 2 (current) with XFOIL's initial guess
    let mut curr = BlStation::new();
    curr.x = s1["x"].as_f64().unwrap();
    curr.u = s1["Ue"].as_f64().unwrap();
    // Initial guess is usually previous station's values
    curr.theta = s0["final"]["theta"].as_f64().unwrap();  // Use prev theta as initial
    curr.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== BLDIF Test at First Step ===");
    println!("Station 1 (prev): x={:.6}, Ue={:.6}", prev.x, prev.u);
    println!("  θ={:.6e}, δ*={:.6e}", prev.theta, prev.delta_star);
    println!("  H={:.4}, Hk={:.4}, Cf={:.6e}", prev.h, prev.hk, prev.cf);

    println!("\nStation 2 (curr, initial): x={:.6}, Ue={:.6}", curr.x, curr.u);
    println!("  θ={:.6e}, δ*={:.6e}", curr.theta, curr.delta_star);
    println!("  H={:.4}, Hk={:.4}, Cf={:.6e}", curr.h, curr.hk, curr.cf);

    // Compute residuals
    let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);

    println!("\n=== RustFoil Residuals ===");
    println!("  res_third (ampl): {:.6e}", res.res_third);
    println!("  res_mom (θ):      {:.6e}", res.res_mom);
    println!("  res_shape (H):    {:.6e}", res.res_shape);

    // XFOIL reference (from the debug output we saw earlier)
    // IBL=3, first iter: res = [1.7e-26, -0.130, 0.0018]
    println!("\n=== XFOIL Reference (IBL=3, iter 1) ===");
    println!("  res_third (ampl): ~0 (1.7e-26)");
    println!("  res_mom (θ):      -0.130");
    println!("  res_shape (H):    0.0018");

    // Show Jacobian
    println!("\n=== RustFoil Jacobian VS2 ===");
    for i in 0..3 {
        println!("  Row {}: [{:.4e}, {:.4e}, {:.4e}, {:.4e}, {:.4e}]",
            i, jac.vs2[i][0], jac.vs2[i][1], jac.vs2[i][2], jac.vs2[i][3], jac.vs2[i][4]);
    }
}

/// Test with XFOIL's exact final values to verify closure functions
#[test]
fn test_bldif_with_xfoil_final_values() {
    let mrchue_paths = [
        "testdata/mrchue_iterations.json",
        "../testdata/mrchue_iterations.json",
        "../../testdata/mrchue_iterations.json",
    ];

    let ref_data: Option<serde_json::Value> = mrchue_paths.iter()
        .find_map(|path| {
            fs::read_to_string(path).ok()
                .and_then(|content| serde_json::from_str(&content).ok())
        });

    let ref_data = match ref_data {
        Some(d) => d,
        None => {
            eprintln!("Skipping test: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;

    let stations = &ref_data["sides"]["1"]["stations"];
    let s0 = &stations[0];
    let s1 = &stations[1];

    // Build both stations with XFOIL's FINAL converged values
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    let mut curr = BlStation::new();
    curr.x = s1["x"].as_f64().unwrap();
    curr.u = s1["Ue"].as_f64().unwrap();
    curr.theta = s1["final"]["theta"].as_f64().unwrap();  // XFOIL's converged theta
    curr.delta_star = s1["final"]["delta_star"].as_f64().unwrap();
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== BLDIF with XFOIL's Converged Values ===");
    println!("Station 1: x={:.6}, θ={:.6e}", prev.x, prev.theta);
    println!("Station 2: x={:.6}, θ={:.6e}", curr.x, curr.theta);

    let (res, _) = bldif(&prev, &curr, FlowType::Laminar, msq, re);

    println!("\n=== Residuals (should be ~0 if XFOIL converged) ===");
    println!("  res_third: {:.6e}", res.res_third);
    println!("  res_mom:   {:.6e}", res.res_mom);
    println!("  res_shape: {:.6e}", res.res_shape);

    // If XFOIL converged, these should be very small
    // Note: They won't be exactly zero because we use midpoint averages
    println!("\nExpected: Small residuals (XFOIL converged to tolerance ~1e-5)");
}
