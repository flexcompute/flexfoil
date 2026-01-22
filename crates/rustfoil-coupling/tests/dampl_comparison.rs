//! Compare our DAMPL against XFOIL test vectors

use rustfoil_bl::closures::amplification_rate;
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct DamplVector {
    #[serde(rename = "Hk")]
    hk: f64,
    theta: f64,
    #[serde(rename = "Rtheta")]
    rtheta: f64,
    #[serde(rename = "Ax")]
    ax: f64,
    #[serde(rename = "Ax_Hk")]
    ax_hk: f64,
    #[serde(rename = "Ax_theta")]
    ax_theta: f64,
    #[serde(rename = "Ax_Rt")]
    ax_rt: f64,
}

fn load_vectors() -> Vec<DamplVector> {
    let paths = [
        "testdata/dampl_test_vectors.json",
        "../testdata/dampl_test_vectors.json",
        "../../testdata/dampl_test_vectors.json",
    ];
    
    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(vectors) = serde_json::from_str(&content) {
                return vectors;
            }
        }
    }
    vec![]
}

#[test]
fn test_dampl_against_xfoil_vectors() {
    let vectors = load_vectors();
    if vectors.is_empty() {
        println!("Warning: Could not load DAMPL test vectors");
        return;
    }
    
    println!("\n=== DAMPL Comparison with XFOIL Test Vectors ===");
    println!("{:>8} {:>10} {:>10} {:>12} {:>12} {:>10}",
             "Hk", "theta", "Rtheta", "XFOIL_Ax", "Rust_Ax", "Error%");
    println!("{}", "-".repeat(70));
    
    let mut total_error = 0.0;
    let mut max_error = 0.0;
    let mut n_compared = 0;
    
    for v in &vectors {
        let result = amplification_rate(v.hk, v.theta, v.rtheta);
        
        let error_pct = if v.ax.abs() > 1e-10 {
            ((result.ax - v.ax) / v.ax).abs() * 100.0
        } else if result.ax.abs() > 1e-10 {
            999.9
        } else {
            0.0
        };
        
        // Only print interesting cases (non-zero or large error)
        if v.ax.abs() > 1e-10 || result.ax.abs() > 1e-10 {
            println!("{:8.4} {:10.4e} {:10.2} {:12.4e} {:12.4e} {:10.1}",
                     v.hk, v.theta, v.rtheta, v.ax, result.ax, error_pct);
            
            total_error += error_pct;
            if error_pct > max_error {
                max_error = error_pct;
            }
            n_compared += 1;
        }
    }
    
    println!("\nSummary:");
    println!("  Vectors with non-zero Ax: {}", n_compared);
    println!("  Max error: {:.1}%", max_error);
    if n_compared > 0 {
        println!("  Avg error: {:.1}%", total_error / n_compared as f64);
    }
    
    // Check some specific cases in detail
    println!("\n=== Detailed trace for Hk=2.6660, Rt=139.45 ===");
    let hk = 2.6660_f64;
    let theta = 8.8441e-5_f64;
    let rt = 139.45_f64;
    
    // Compute Rcrit
    let hmi = 1.0 / (hk - 1.0);
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let grcrit = aa + 0.7 * (bb + 1.0);
    let rcrit = 10.0_f64.powf(grcrit);
    
    println!("  HMI = {:.6}", hmi);
    println!("  AA = {:.6}", aa);
    println!("  BB = {:.6}", bb);
    println!("  GRCRIT = {:.6}", grcrit);
    println!("  Rcrit = {:.2}", rcrit);
    println!("  log10(Rtheta) = {:.6}", rt.log10());
    println!("  GR - GRCRIT = {:.6}", rt.log10() - grcrit);
    
    let dgr = 0.08;
    let rnorm = (rt.log10() - (grcrit - dgr)) / (2.0 * dgr);
    println!("  RNORM = {:.6}", rnorm);
    
    let rfac = if rnorm >= 1.0 {
        1.0
    } else {
        3.0 * rnorm.powi(2) - 2.0 * rnorm.powi(3)
    };
    println!("  RFAC = {:.6}", rfac);
    
    let result = amplification_rate(hk, theta, rt);
    println!("  Our Ax = {:.6e}", result.ax);
    println!("  XFOIL Ax = {:.6e}", 2.1432e-2);
}
