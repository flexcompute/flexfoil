//! Transition trace comparison test
//!
//! Compares RustFoil's transition detection step-by-step against XFOIL's
//! intermediate values from the instrumented debug output.

use rustfoil_bl::closures::{amplification_rate, axset};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    #[allow(dead_code)]
    mach: f64,
    ncrit: f64,
}

#[derive(Debug, Deserialize)]
struct Station {
    ibl: i32,
    #[allow(dead_code)]
    side: i32,
    x1: Option<f64>,
    x2: Option<f64>,
    ampl1: Option<f64>,
    ampl2: Option<f64>,
    #[serde(rename = "Hk1")]
    hk1: Option<f64>,
    #[serde(rename = "Hk2")]
    hk2: Option<f64>,
    #[serde(rename = "Rt1")]
    rt1: Option<f64>,
    #[serde(rename = "Rt2")]
    rt2: Option<f64>,
    ax_final: Option<f64>,
    transition: Option<bool>,
}

#[derive(Debug, Deserialize)]
struct Side {
    stations: Vec<Station>,
}

#[derive(Debug, Deserialize)]
struct TransitionData {
    metadata: Metadata,
    sides: std::collections::HashMap<String, Side>,
}

fn load_transition_trace() -> Option<TransitionData> {
    let paths = [
        "testdata/transition_trace.json",
        "../testdata/transition_trace.json",
        "../../testdata/transition_trace.json",
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

// Also load mrchue_iterations.json which has actual theta values
#[derive(Debug, Deserialize)]
struct MrchueMetadata {
    reynolds: f64,
    ncrit: f64,
}

#[derive(Debug, Deserialize)]
struct FinalValues {
    theta: f64,
    delta_star: f64,
    ampl: f64,
    #[serde(rename = "Hk", default)]
    hk: f64,
    #[serde(rename = "Rtheta", default)]
    rtheta: f64,
}

#[derive(Debug, Deserialize)]
struct MrchueStation {
    ibl: i32,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
    #[serde(rename = "final")]
    final_values: FinalValues,
}

#[derive(Debug, Deserialize)]
struct MrchueSide {
    stations: Vec<MrchueStation>,
}

#[derive(Debug, Deserialize)]
struct MrchueData {
    metadata: MrchueMetadata,
    sides: std::collections::HashMap<String, MrchueSide>,
}

fn load_mrchue_data() -> Option<MrchueData> {
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

/// Test DAMPL (amplification_rate) against XFOIL values
#[test]
fn test_dampl_vs_xfoil() {
    let data = match load_transition_trace() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load transition trace data");
            return;
        }
    };
    
    let stations = &data.sides.get("1").unwrap().stations;
    
    println!("\n=== DAMPL Comparison (Upper Surface) ===");
    println!("{:>4} {:>8} {:>10} {:>12} {:>12} {:>10}",
             "IBL", "Hk2", "Rt2", "AX_XFOIL", "AX_RUST", "Error%");
    println!("{}", "-".repeat(66));
    
    let mut total_dampl_tests = 0;
    let mut dampl_matches = 0;
    
    for st in stations.iter().take(30) {
        let hk2 = st.hk2.unwrap_or(0.0);
        let rt2 = st.rt2.unwrap_or(0.0);
        
        if hk2 == 0.0 || rt2 == 0.0 {
            continue;
        }
        
        // We need theta to call amplification_rate
        // Estimate from Rt = Re * Ue * theta, assuming Ue ≈ 1.7 for upper surface
        // theta ≈ Rt / (Re * Ue) = Rt / (1e6 * 1.7)
        let re = data.metadata.reynolds;
        let ue_estimate = 1.7;
        let theta = rt2 / (re * ue_estimate);
        
        let rust_result = amplification_rate(hk2, theta, rt2);
        
        // Estimate XFOIL's raw DAMPL result
        // When in subcritical region, DAMPL returns 0
        // When supercritical, it returns the amplification rate
        // We can check if ax_final is primarily from DAX or from DAMPL
        
        let ax_xfoil = st.ax_final.unwrap_or(0.0);
        let ax_rust = rust_result.ax;
        
        // Check if both are in subcritical region (AX ≈ 0 or very small)
        let both_subcritical = ax_rust < 1e-6 && ax_xfoil < 1e-6;
        let both_supercritical = ax_rust > 0.01 && ax_xfoil > 0.01;
        
        total_dampl_tests += 1;
        
        if both_subcritical || both_supercritical {
            dampl_matches += 1;
        }
        
        let error = if ax_xfoil.abs() > 1e-10 {
            (ax_rust - ax_xfoil) / ax_xfoil * 100.0
        } else if ax_rust.abs() < 1e-10 {
            0.0  // Both effectively zero
        } else {
            f64::NAN
        };
        
        let status = if both_subcritical { "sub" } 
                    else if both_supercritical { "sup" } 
                    else { "???" };
        
        println!("{:4} {:8.4} {:10.2} {:12.4e} {:12.4e} {:10.1} {}",
                 st.ibl, hk2, rt2, ax_xfoil, ax_rust, error, status);
    }
    
    println!();
    println!("DAMPL region agreement: {}/{} stations", dampl_matches, total_dampl_tests);
}

/// Test AXSET against XFOIL values
#[test]
fn test_axset_vs_xfoil() {
    let data = match load_transition_trace() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load transition trace data");
            return;
        }
    };
    
    let ncrit = data.metadata.ncrit;
    let re = data.metadata.reynolds;
    let stations = &data.sides.get("1").unwrap().stations;
    
    println!("\n=== AXSET Comparison (Upper Surface) ===");
    println!("{:>4} {:>10} {:>10} {:>12} {:>12} {:>10}",
             "IBL", "ampl1", "ax_XFOIL", "ax_RUST", "dN_XFOIL", "dN_RUST");
    println!("{}", "-".repeat(70));
    
    let ue_estimate = 1.7;
    
    for st in stations.iter().take(25) {
        let hk1 = st.hk1.unwrap_or(0.0);
        let hk2 = st.hk2.unwrap_or(0.0);
        let rt1 = st.rt1.unwrap_or(0.0);
        let rt2 = st.rt2.unwrap_or(0.0);
        let ampl1 = st.ampl1.unwrap_or(0.0);
        let ampl2 = st.ampl2.unwrap_or(0.0);
        let x1 = st.x1.unwrap_or(0.0);
        let x2 = st.x2.unwrap_or(0.0);
        
        if hk1 == 0.0 || hk2 == 0.0 {
            continue;
        }
        
        // Estimate theta from Rt
        let th1 = rt1 / (re * ue_estimate);
        let th2 = rt2 / (re * ue_estimate);
        
        // Call RustFoil's axset
        let rust_result = axset(hk1, th1, rt1, ampl1, hk2, th2, rt2, ampl1, ncrit);
        
        let ax_xfoil = st.ax_final.unwrap_or(0.0);
        let ax_rust = rust_result.ax;
        
        let dx = x2 - x1;
        let dn_xfoil = ampl2 - ampl1;
        let dn_rust = ax_rust * dx;
        
        println!("{:4} {:10.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e}",
                 st.ibl, ampl1, ax_xfoil, ax_rust, dn_xfoil, dn_rust);
    }
}

/// Test N-factor evolution - compare accumulated N vs XFOIL
/// Uses actual theta values from mrchue_iterations.json
#[test]
fn test_nfactor_evolution() {
    let mrchue = match load_mrchue_data() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load mrchue_iterations.json");
            return;
        }
    };
    
    let ncrit = mrchue.metadata.ncrit;
    let stations = &mrchue.sides.get("1").unwrap().stations;
    
    println!("\n=== N-Factor Evolution Comparison (Upper Surface) ===");
    println!("Using actual theta values from mrchue_iterations.json");
    println!();
    println!("{:>4} {:>8} {:>12} {:>12} {:>12} {:>10}",
             "IBL", "x", "N_XFOIL", "N_RUST", "dN_XFOIL", "dN_RUST");
    println!("{}", "-".repeat(70));
    
    // Accumulate N-factor as we march
    let mut n_rust = 0.0;
    let mut prev_x = 0.0;
    let mut prev_hk = 0.0;
    let mut prev_th = 0.0;
    let mut prev_rt = 0.0;
    let mut prev_ampl = 0.0;
    
    let mut transition_xfoil = None;
    let mut transition_rust = None;
    
    for (i, st) in stations.iter().enumerate() {
        let hk2 = st.final_values.hk;
        let rt2 = st.final_values.rtheta;
        let th2 = st.final_values.theta;
        let x2 = st.x;
        let ampl2_xfoil = st.final_values.ampl;
        
        if hk2 == 0.0 || th2 == 0.0 {
            continue;
        }
        
        // Compute dN using axset
        let dx = x2 - prev_x;
        
        if i > 0 && dx > 0.0 {
            let hk1 = prev_hk;
            let th1 = prev_th;
            let rt1 = prev_rt;
            
            // Use previous N for ampl1 (RustFoil's accumulated value)
            let axset_result = axset(hk1, th1, rt1, n_rust, hk2, th2, rt2, n_rust, ncrit);
            let dn_rust = axset_result.ax * dx;
            n_rust += dn_rust;
            
            let dn_xfoil = ampl2_xfoil - prev_ampl;
            
            // Check for transition
            if n_rust >= ncrit && transition_rust.is_none() {
                transition_rust = Some((st.ibl, x2));
            }
            if ampl2_xfoil >= ncrit && transition_xfoil.is_none() {
                transition_xfoil = Some((st.ibl, x2));
            }
            
            let n_error = if ampl2_xfoil.abs() > 1e-10 {
                (n_rust - ampl2_xfoil) / ampl2_xfoil * 100.0
            } else if n_rust.abs() < 1e-10 {
                0.0
            } else {
                f64::NAN
            };
            
            let marker = if n_rust >= ncrit && transition_rust == Some((st.ibl, x2)) { " <-- RUST transition" }
                        else if ampl2_xfoil >= ncrit && transition_xfoil == Some((st.ibl, x2)) { " <-- XFOIL transition" }
                        else { "" };
            
            if i < 35 || n_rust >= ncrit * 0.8 || ampl2_xfoil >= ncrit * 0.8 {
                println!("{:4} {:8.4} {:12.4e} {:12.4e} {:12.4e} {:10.4e} {:6.1}%{}",
                         st.ibl, x2, ampl2_xfoil, n_rust, dn_xfoil, dn_rust, n_error, marker);
            }
        }
        
        prev_x = x2;
        prev_hk = hk2;
        prev_th = th2;
        prev_rt = rt2;
        prev_ampl = ampl2_xfoil;
    }
    
    println!();
    if let Some((ibl, x)) = transition_xfoil {
        println!("XFOIL transition: IBL {} at x = {:.4}", ibl, x);
    }
    if let Some((ibl, x)) = transition_rust {
        println!("RustFoil transition: IBL {} at x = {:.4}", ibl, x);
    }
    
    // Check that transition locations are close
    if let (Some((ibl_xfoil, _)), Some((ibl_rust, _))) = (transition_xfoil, transition_rust) {
        let diff = (ibl_xfoil - ibl_rust).abs();
        println!("\nTransition station difference: {} stations", diff);
        assert!(diff <= 2, "Transition location should be within 2 stations of XFOIL");
    }
}

/// Test critical Reynolds number (Rcrit) calculation
#[test]
fn test_rcrit_calculation() {
    let data = match load_transition_trace() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load transition trace data");
            return;
        }
    };
    
    let stations = &data.sides.get("1").unwrap().stations;
    
    println!("\n=== Rcrit Analysis (Upper Surface) ===");
    println!("Shows when Rtheta exceeds Rcrit (amplification onset)");
    println!();
    println!("{:>4} {:>8} {:>10} {:>12} {:>12} {:>8}",
             "IBL", "Hk", "Rtheta", "Rcrit", "Rth/Rcrit", "Status");
    println!("{}", "-".repeat(66));
    
    let dgr = 0.08_f64;  // XFOIL's DGR constant
    
    for st in stations.iter().take(30) {
        let hk = st.hk2.unwrap_or(0.0);
        let rt = st.rt2.unwrap_or(0.0);
        
        if hk <= 1.0 || rt == 0.0 {
            continue;
        }
        
        // Compute Rcrit using DAMPL correlation
        let hmi = 1.0 / (hk - 1.0);
        let aa = 2.492 * hmi.powf(0.43);
        let bb = (14.0 * hmi - 9.24).tanh();
        let log_rcrit = aa + 0.7 * (bb + 1.0);
        let rcrit = 10.0_f64.powf(log_rcrit);
        
        // Compute the ramp region boundaries
        let rcrit_lower = 10.0_f64.powf(log_rcrit - dgr);  // Below this: AX = 0
        let rcrit_upper = 10.0_f64.powf(log_rcrit + dgr);  // Above this: AX = full
        
        let ratio = rt / rcrit;
        
        let status = if rt < rcrit_lower {
            "subcrit"
        } else if rt < rcrit_upper {
            "ramp"
        } else {
            "supercrit"
        };
        
        let marker = if status == "supercrit" { "*" } else { "" };
        
        println!("{:4} {:8.4} {:10.2} {:12.2} {:12.4} {:>8} {}",
                 st.ibl, hk, rt, rcrit, ratio, status, marker);
    }
    
    println!();
    println!("Note: * marks stations where amplification rate is fully active");
}
