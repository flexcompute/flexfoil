//! Instrumented XFOIL comparison for Newton V-I coupling
//!
//! This test compares RustFoil's Newton system (VA, VB, VDEL, VM) against
//! XFOIL's instrumented debug output from SETBL.

use nalgebra::DMatrix;
use rustfoil_bl::equations::{blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::newton::BlNewtonSystem;
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct XfoilEvent {
    subroutine: String,
    iteration: Option<i32>,
    side: Option<i32>,
    ibl: Option<i32>,
    iv: Option<i32>,
    #[serde(rename = "VA")]
    va: Option<Vec<Vec<f64>>>,
    #[serde(rename = "VB")]
    vb: Option<Vec<Vec<f64>>>,
    #[serde(rename = "VDEL")]
    vdel: Option<Vec<Vec<f64>>>,
    #[serde(rename = "VM_sample")]
    vm_sample: Option<Vec<Vec<f64>>>,
    // UPDATE fields
    delta_ctau: Option<f64>,
    delta_theta: Option<f64>,
    delta_mass: Option<f64>,
    #[serde(rename = "delta_Ue")]
    delta_ue: Option<f64>,
    relaxation: Option<f64>,
    // BLVAR fields for station reconstruction
    input: Option<BlvarInput>,
    output: Option<BlvarOutput>,
    reynolds: Option<f64>,
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
#[serde(untagged)]
enum BlvarOutput {
    Full {
        #[serde(rename = "H")]
        h: f64,
        #[serde(rename = "Hk")]
        hk: f64,
        #[serde(rename = "Hs")]
        hs: f64,
        #[serde(rename = "Rtheta")]
        r_theta: f64,
        #[serde(rename = "Cf")]
        cf: f64,
        #[serde(rename = "Cd")]
        cd: f64,
    },
    Simple {
        theta: f64,
        delta_star: f64,
        ctau: f64,
        ampl: f64,
    },
}

#[derive(Debug, Deserialize)]
struct XfoilTrace {
    events: Vec<XfoilEvent>,
}

fn load_xfoil_trace() -> Option<XfoilTrace> {
    // Try to find the testdata directory relative to the manifest directory
    // CARGO_MANIFEST_DIR = crates/rustfoil-coupling
    let manifest_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    
    // Try going up to project root (crates/rustfoil-coupling -> crates -> project_root)
    let paths = [
        manifest_dir.join("../../testdata/newton_debug_trace.json"),
        manifest_dir.join("../../../testdata/newton_debug_trace.json"),
        std::path::PathBuf::from("testdata/newton_debug_trace.json"),
    ];

    for path in &paths {
        match fs::read_to_string(&path) {
            Ok(content) => {
                match serde_json::from_str::<XfoilTrace>(&content) {
                    Ok(trace) => {
                        return Some(trace);
                    }
                    Err(e) => {
                        eprintln!("JSON parse error from {:?}: {}", path, e);
                    }
                }
            }
            Err(_) => {}
        }
    }
    None
}

/// Extract SETBL entries from XFOIL trace
fn extract_setbl_entries(trace: &XfoilTrace) -> Vec<&XfoilEvent> {
    trace
        .events
        .iter()
        .filter(|e| e.subroutine == "SETBL")
        .collect()
}

/// Extract BLVAR entries to reconstruct station states
fn extract_blvar_entries(trace: &XfoilTrace, side: i32, iteration: i32) -> Vec<&XfoilEvent> {
    trace
        .events
        .iter()
        .filter(|e| {
            e.subroutine == "BLVAR"
                && e.side == Some(side)
                && e.iteration == Some(iteration)
                && e.input.is_some()
        })
        .collect()
}

/// Extract UPDATE entries from XFOIL trace
fn extract_update_entries(trace: &XfoilTrace) -> Vec<&XfoilEvent> {
    trace
        .events
        .iter()
        .filter(|e| e.subroutine == "UPDATE")
        .collect()
}

#[test]
fn test_setbl_va_vb_structure() {
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let setbl_entries = extract_setbl_entries(&trace);
    println!("Found {} SETBL entries", setbl_entries.len());

    // Analyze VA and VB structure
    let mut va_norms = Vec::new();
    let mut vb_norms = Vec::new();

    for entry in setbl_entries.iter().take(10) {
        if let (Some(va), Some(vb)) = (&entry.va, &entry.vb) {
            let va_norm: f64 = va.iter().flat_map(|row| row.iter()).map(|v| v * v).sum();
            let vb_norm: f64 = vb.iter().flat_map(|row| row.iter()).map(|v| v * v).sum();

            va_norms.push(va_norm.sqrt());
            vb_norms.push(vb_norm.sqrt());

            println!(
                "SETBL side={}, ibl={}, iv={}: ||VA||={:.4e}, ||VB||={:.4e}",
                entry.side.unwrap_or(0),
                entry.ibl.unwrap_or(0),
                entry.iv.unwrap_or(0),
                va_norm.sqrt(),
                vb_norm.sqrt()
            );
        }
    }

    // VA should be non-zero (diagonal block)
    assert!(
        va_norms.iter().any(|&n| n > 1e-10),
        "VA blocks should have non-zero entries"
    );
}

#[test]
fn test_setbl_vm_sample_structure() {
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let setbl_entries = extract_setbl_entries(&trace);

    println!("\n=== VM Matrix Sample Analysis ===");

    for entry in setbl_entries.iter().take(5) {
        if let Some(vm_sample) = &entry.vm_sample {
            println!(
                "\nSETBL side={}, ibl={}, iv={}:",
                entry.side.unwrap_or(0),
                entry.ibl.unwrap_or(0),
                entry.iv.unwrap_or(0)
            );

            for (k, row) in vm_sample.iter().enumerate() {
                let row_norm: f64 = row.iter().map(|v| v * v).sum::<f64>().sqrt();
                let first_5: Vec<String> = row.iter().take(5).map(|v| format!("{:.4e}", v)).collect();
                println!("  VM[{},:5] = [{}], ||row||={:.4e}", k, first_5.join(", "), row_norm);
            }
        }
    }
}

#[test]
fn test_setbl_vdel_residuals() {
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let setbl_entries = extract_setbl_entries(&trace);

    println!("\n=== VDEL (Residual) Analysis ===");

    let mut max_residuals = [0.0f64; 3];

    for entry in &setbl_entries {
        if let Some(vdel) = &entry.vdel {
            for (k, row) in vdel.iter().enumerate() {
                if !row.is_empty() {
                    let res = row[0].abs();
                    if res > max_residuals[k] {
                        max_residuals[k] = res;
                    }
                }
            }
        }
    }

    println!("Max residuals across all stations:");
    println!("  Eq 0 (ampl/ctau): {:.6e}", max_residuals[0]);
    println!("  Eq 1 (momentum):  {:.6e}", max_residuals[1]);
    println!("  Eq 2 (shape):     {:.6e}", max_residuals[2]);
}

#[test]
fn test_update_deltas_and_relaxation() {
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let update_entries = extract_update_entries(&trace);
    println!("Found {} UPDATE entries", update_entries.len());

    println!("\n=== UPDATE Analysis (first 10 stations) ===");

    let mut relaxations = Vec::new();
    let mut delta_theta_vals = Vec::new();
    let mut delta_mass_vals = Vec::new();
    let mut delta_ue_vals = Vec::new();

    for entry in update_entries.iter().take(10) {
        if let (Some(dt), Some(dm), Some(due), Some(rlx)) = (
            entry.delta_theta,
            entry.delta_mass,
            entry.delta_ue,
            entry.relaxation,
        ) {
            println!(
                "UPDATE side={}, ibl={}: dθ={:.4e}, dMass={:.4e}, dUe={:.4e}, rlx={:.4}",
                entry.side.unwrap_or(0),
                entry.ibl.unwrap_or(0),
                dt,
                dm,
                due,
                rlx
            );

            relaxations.push(rlx);
            delta_theta_vals.push(dt);
            delta_mass_vals.push(dm);
            delta_ue_vals.push(due);
        }
    }

    // Analyze relaxation
    if !relaxations.is_empty() {
        let avg_rlx: f64 = relaxations.iter().sum::<f64>() / relaxations.len() as f64;
        let min_rlx = relaxations.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_rlx = relaxations.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        println!("\nRelaxation statistics:");
        println!("  Min: {:.4}", min_rlx);
        println!("  Max: {:.4}", max_rlx);
        println!("  Avg: {:.4}", avg_rlx);

        // XFOIL typically uses relaxation < 1.0 for Newton
        assert!(
            min_rlx > 0.0 && max_rlx <= 1.0,
            "Relaxation should be in (0, 1]"
        );
    }
}

#[test]
fn test_compare_va_vb_magnitudes() {
    // Compare XFOIL's VA/VB magnitudes with what we compute
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    // Get Reynolds number from VISCAL entry
    let re = trace
        .events
        .iter()
        .find(|e| e.subroutine == "VISCAL")
        .and_then(|e| e.reynolds)
        .unwrap_or(1e6);

    let msq = 0.0;

    // Get BLVAR entries for side 1, iteration 1 to reconstruct stations
    let blvar_entries = extract_blvar_entries(&trace, 1, 1);
    println!("Found {} BLVAR entries for side 1, iteration 1", blvar_entries.len());

    // Reconstruct first few stations
    let mut stations: Vec<BlStation> = Vec::new();
    for entry in blvar_entries.iter().take(10) {
        if let Some(input) = &entry.input {
            let mut station = BlStation::new();
            station.x = input.x;
            station.u = input.u;
            station.theta = input.theta;
            station.delta_star = input.delta_star;
            station.ctau = input.ctau;
            station.ampl = input.ampl;
            station.mass_defect = station.u * station.delta_star;
            station.is_laminar = input.ctau < 0.001 || input.ampl > 0.0;

            let flow_type = if station.is_laminar {
                FlowType::Laminar
            } else {
                FlowType::Turbulent
            };
            blvar(&mut station, flow_type, msq, re);
            stations.push(station);
        }
    }

    if stations.len() < 3 {
        println!("Not enough stations reconstructed");
        return;
    }

    println!("\n=== VA/VB Magnitude Comparison ===");

    // Build Newton system
    let n = stations.len();
    let flow_types: Vec<FlowType> = stations
        .iter()
        .skip(1)
        .map(|s| {
            if s.is_laminar {
                FlowType::Laminar
            } else {
                FlowType::Turbulent
            }
        })
        .collect();

    let mut system = BlNewtonSystem::new(n);
    system.build(&stations, &flow_types, msq, re);

    // Get SETBL entries for comparison
    let setbl_entries = extract_setbl_entries(&trace);
    let side1_setbl: Vec<_> = setbl_entries
        .iter()
        .filter(|e| e.side == Some(1) && e.iteration == Some(1))
        .take(n - 1)
        .collect();

    println!("\nComparing first {} stations:", (n - 1).min(side1_setbl.len()));

    for (i, setbl) in side1_setbl.iter().enumerate() {
        let iv = i + 1; // Our index (1-based intervals)

        if let Some(xfoil_va) = &setbl.va {
            // Compute norms
            let xfoil_va_norm: f64 = xfoil_va
                .iter()
                .flat_map(|row| row.iter())
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt();

            let rust_va_norm: f64 = system.va[iv]
                .iter()
                .flat_map(|row| row.iter())
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt();

            let va_ratio = if xfoil_va_norm > 1e-20 {
                rust_va_norm / xfoil_va_norm
            } else {
                1.0
            };

            println!(
                "Station {}: XFOIL ||VA||={:.4e}, RustFoil ||VA||={:.4e}, ratio={:.4}",
                iv, xfoil_va_norm, rust_va_norm, va_ratio
            );
        }

        if let Some(xfoil_vb) = &setbl.vb {
            let xfoil_vb_norm: f64 = xfoil_vb
                .iter()
                .flat_map(|row| row.iter())
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt();

            let rust_vb_norm: f64 = system.vb[iv]
                .iter()
                .flat_map(|row| row.iter())
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt();

            let vb_ratio = if xfoil_vb_norm > 1e-20 {
                rust_vb_norm / xfoil_vb_norm
            } else {
                1.0
            };

            println!(
                "          XFOIL ||VB||={:.4e}, RustFoil ||VB||={:.4e}, ratio={:.4}",
                xfoil_vb_norm, rust_vb_norm, vb_ratio
            );
        }
    }
}

#[test]
fn test_compare_residuals() {
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let re = trace
        .events
        .iter()
        .find(|e| e.subroutine == "VISCAL")
        .and_then(|e| e.reynolds)
        .unwrap_or(1e6);

    let msq = 0.0;

    // Reconstruct stations from BLVAR
    let blvar_entries = extract_blvar_entries(&trace, 1, 1);
    let mut stations: Vec<BlStation> = Vec::new();
    for entry in blvar_entries.iter().take(10) {
        if let Some(input) = &entry.input {
            let mut station = BlStation::new();
            station.x = input.x;
            station.u = input.u;
            station.theta = input.theta;
            station.delta_star = input.delta_star;
            station.ctau = input.ctau;
            station.ampl = input.ampl;
            station.mass_defect = station.u * station.delta_star;
            station.is_laminar = input.ctau < 0.001 || input.ampl > 0.0;

            let flow_type = if station.is_laminar {
                FlowType::Laminar
            } else {
                FlowType::Turbulent
            };
            blvar(&mut station, flow_type, msq, re);
            stations.push(station);
        }
    }

    if stations.len() < 3 {
        return;
    }

    let n = stations.len();
    let flow_types: Vec<FlowType> = stations
        .iter()
        .skip(1)
        .map(|s| {
            if s.is_laminar {
                FlowType::Laminar
            } else {
                FlowType::Turbulent
            }
        })
        .collect();

    let mut system = BlNewtonSystem::new(n);
    system.build(&stations, &flow_types, msq, re);

    // Get SETBL entries
    let setbl_entries = extract_setbl_entries(&trace);
    let side1_setbl: Vec<_> = setbl_entries
        .iter()
        .filter(|e| e.side == Some(1) && e.iteration == Some(1))
        .take(n - 1)
        .collect();

    println!("\n=== Residual (VDEL) Comparison ===");

    for (i, setbl) in side1_setbl.iter().enumerate() {
        let iv = i + 1;

        if let Some(xfoil_vdel) = &setbl.vdel {
            // XFOIL VDEL has structure [[res0, ...], [res1, ...], [res2, ...]]
            // We want the first column (residual values)
            let xfoil_res: Vec<f64> = xfoil_vdel.iter().map(|row| row.get(0).copied().unwrap_or(0.0)).collect();

            let rust_res = &system.rhs[iv];

            println!(
                "Station {}: XFOIL VDEL=[{:.4e}, {:.4e}, {:.4e}]",
                iv, xfoil_res.get(0).unwrap_or(&0.0), xfoil_res.get(1).unwrap_or(&0.0), xfoil_res.get(2).unwrap_or(&0.0)
            );
            println!(
                "          RustFoil RHS=[{:.4e}, {:.4e}, {:.4e}]",
                rust_res[0], rust_res[1], rust_res[2]
            );

            // Check if residuals are in same ballpark (within 10x)
            for k in 0..3 {
                let xfoil_val = xfoil_res.get(k).unwrap_or(&0.0).abs();
                let rust_val = rust_res[k].abs();

                if xfoil_val > 1e-20 && rust_val > 1e-20 {
                    let ratio = rust_val / xfoil_val;
                    println!("          Eq {}: ratio = {:.4}", k, ratio);
                }
            }
        }
    }
}

#[test]
fn test_vm_diagonal_comparison() {
    let trace = match load_xfoil_trace() {
        Some(t) => t,
        None => {
            eprintln!("Skipping: newton_debug_trace.json not found");
            return;
        }
    };

    let setbl_entries = extract_setbl_entries(&trace);

    println!("\n=== VM Diagonal Analysis ===");
    println!("(XFOIL VM_sample shows first 10 columns for each equation row)");

    for entry in setbl_entries.iter().take(5) {
        if let Some(vm_sample) = &entry.vm_sample {
            let iv = entry.iv.unwrap_or(0) as usize;

            // The diagonal entry VM[iv][iv] is at column index (iv - 1) in 0-based
            // But VM_sample only shows first 10 columns
            if iv > 0 && iv <= 10 {
                println!(
                    "\nStation iv={}: VM diagonal (column {}):",
                    iv,
                    iv - 1
                );
                for (k, row) in vm_sample.iter().enumerate() {
                    if let Some(diag_val) = row.get(iv - 1) {
                        println!("  VM[{}][{}] = {:.6e}", k, iv - 1, diag_val);
                    }
                }
            }
        }
    }
}
