//! Detailed amplification debugging test
//!
//! Compares RustFoil's amplification_rate and axset calculations against
//! XFOIL reference data station-by-station to find where differences arise.

use rustfoil_bl::closures::{amplification_rate, axset};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct Metadata {
    reynolds: f64,
    mach: f64,
    ncrit: f64,
}

#[derive(Debug, Deserialize)]
struct FinalValues {
    theta: f64,
    delta_star: f64,
    #[serde(default)]
    ampl: f64,
    #[serde(rename = "Hk", default)]
    hk: f64,
    #[serde(rename = "Rtheta", default)]
    rtheta: f64,
}

#[derive(Debug, Deserialize)]
struct Station {
    ibl: i32,
    x: f64,
    #[serde(rename = "Ue")]
    ue: f64,
    #[serde(rename = "final")]
    final_values: FinalValues,
}

#[derive(Debug, Deserialize)]
struct Side {
    stations: Vec<Station>,
}

#[derive(Debug, Deserialize)]
struct ReferenceData {
    metadata: Metadata,
    sides: std::collections::HashMap<String, Side>,
}

fn load_reference() -> Option<ReferenceData> {
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
fn test_amplification_rate_comparison() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load reference data");
            return;
        }
    };
    
    let ncrit = ref_data.metadata.ncrit;
    let stations = &ref_data.sides.get("1").unwrap().stations;
    
    println!("\n=== Amplification Rate (DAMPL) Comparison ===");
    println!("{:>4} {:>10} {:>10} {:>12} {:>12} {:>12} {:>10}",
             "IBL", "Hk", "theta", "Rtheta", "AX_XFOIL", "AX_RUST", "Error%");
    println!("{}", "-".repeat(76));
    
    // For each station, compute amplification rate and compare
    // We can estimate XFOIL's AX from the N-factor differences
    for i in 1..stations.len().min(35) {
        let prev = &stations[i - 1];
        let curr = &stations[i];
        
        let hk = curr.final_values.hk;
        let theta = curr.final_values.theta;
        let rtheta = curr.final_values.rtheta;
        
        // Compute RustFoil's amplification rate
        let rust_result = amplification_rate(hk, theta, rtheta);
        
        // Estimate XFOIL's AX from dN/dx = (N2 - N1) / (x2 - x1)
        let dn = curr.final_values.ampl - prev.final_values.ampl;
        let dx = curr.x - prev.x;
        let ax_xfoil_approx = if dx > 0.0 { dn / dx } else { 0.0 };
        
        let error_pct = if ax_xfoil_approx.abs() > 1e-10 {
            ((rust_result.ax - ax_xfoil_approx) / ax_xfoil_approx * 100.0)
        } else if rust_result.ax.abs() > 1e-10 {
            999.9  // XFOIL has 0, we don't
        } else {
            0.0
        };
        
        println!("{:4} {:10.4} {:10.6e} {:12.2} {:12.4e} {:12.4e} {:>10.1}",
                 curr.ibl, hk, theta, rtheta, ax_xfoil_approx, rust_result.ax, error_pct);
    }
}

#[test]
fn test_axset_comparison() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load reference data");
            return;
        }
    };
    
    let ncrit = ref_data.metadata.ncrit;
    let stations = &ref_data.sides.get("1").unwrap().stations;
    
    println!("\n=== AXSET (RMS Averaged) Comparison ===");
    println!("{:>4} {:>12} {:>12} {:>12} {:>12} {:>10}",
             "IBL", "dN(XFOIL)", "AXSET*dx", "dN_diff", "AMPL_err%", "AX1/AX2");
    println!("{}", "-".repeat(76));
    
    let mut our_ampl = 0.0;
    
    for i in 1..stations.len().min(35) {
        let prev = &stations[i - 1];
        let curr = &stations[i];
        
        let hk1 = prev.final_values.hk;
        let t1 = prev.final_values.theta;
        let rt1 = prev.final_values.rtheta;
        let ampl1_xfoil = prev.final_values.ampl;
        
        let hk2 = curr.final_values.hk;
        let t2 = curr.final_values.theta;
        let rt2 = curr.final_values.rtheta;
        let ampl2_xfoil = curr.final_values.ampl;
        
        let dx = curr.x - prev.x;
        
        // Compute our AXSET
        let axset_result = axset(hk1, t1, rt1, our_ampl, hk2, t2, rt2, our_ampl, ncrit);
        
        // Update our amplitude
        let dn_ours = axset_result.ax * dx;
        our_ampl += dn_ours;
        
        let dn_xfoil = ampl2_xfoil - ampl1_xfoil;
        let dn_diff = dn_ours - dn_xfoil;
        
        let ampl_err_pct = if ampl2_xfoil.abs() > 1e-10 {
            (our_ampl - ampl2_xfoil) / ampl2_xfoil * 100.0
        } else if our_ampl.abs() > 1e-10 {
            999.9
        } else {
            0.0
        };
        
        println!("{:4} {:12.6e} {:12.6e} {:12.6e} {:12.1} {:>10}",
                 curr.ibl, dn_xfoil, dn_ours, dn_diff, ampl_err_pct, 
                 format!("{:.2e}/{:.2e}", axset_result.ax1, axset_result.ax2));
    }
    
    println!("\nFinal: XFOIL ampl = {:.6e}, Our ampl = {:.6e}", 
             stations[34].final_values.ampl, our_ampl);
}

#[test]
fn test_rcrit_calculation() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            println!("Warning: Could not load reference data");
            return;
        }
    };
    
    let stations = &ref_data.sides.get("1").unwrap().stations;
    
    println!("\n=== Rcrit vs Rtheta Analysis ===");
    println!("{:>4} {:>8} {:>12} {:>12} {:>12} {:>12}",
             "IBL", "Hk", "Rtheta", "Rcrit", "Rth/Rcrit", "AX");
    println!("{}", "-".repeat(72));
    
    for i in 0..stations.len().min(30) {
        let st = &stations[i];
        let hk = st.final_values.hk;
        let rtheta = st.final_values.rtheta;
        let theta = st.final_values.theta;
        
        // Compute Rcrit from Hk (matching XFOIL's DAMPL correlation)
        let hmi = 1.0 / (hk - 1.0);
        let aa = 2.492 * hmi.powf(0.43);
        let bb = (14.0 * hmi - 9.24).tanh();
        let log_rcrit = aa + 0.7 * (bb + 1.0);
        let rcrit = 10.0_f64.powf(log_rcrit);
        
        let ratio = rtheta / rcrit;
        
        // Compute AX
        let ax_result = amplification_rate(hk, theta, rtheta);
        
        let marker = if ratio > 1.0 { "*" } else { "" };
        
        println!("{:4} {:8.4} {:12.2} {:12.2} {:12.4} {:12.4e} {}",
                 st.ibl, hk, rtheta, rcrit, ratio, ax_result.ax, marker);
    }
}
