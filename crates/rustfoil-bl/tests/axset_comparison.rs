//! AXSET function comparison test
//!
//! Compares our axset function against XFOIL's AXSET output
//! to verify the RMS averaging implementation matches.

use rustfoil_bl::closures::axset;
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct AxsetEvent {
    #[serde(rename = "Hk1")]
    hk1: f64,
    theta1: f64,
    #[serde(rename = "Rtheta1")]
    rtheta1: f64,
    ampl1: f64,
    #[serde(rename = "Hk2")]
    hk2: f64,
    theta2: f64,
    #[serde(rename = "Rtheta2")]
    rtheta2: f64,
    ampl2: f64,
    #[serde(rename = "Ncrit")]
    ncrit: f64,
    #[serde(rename = "Ax1")]
    ax1: f64,
    #[serde(rename = "Ax2")]
    ax2: f64,
    #[serde(rename = "Axa_rms")]
    axa_rms: f64,
    #[serde(rename = "Dax")]
    dax: f64,
    #[serde(rename = "Ax_final")]
    ax_final: f64,
}

fn load_axset_events() -> Option<Vec<AxsetEvent>> {
    let paths = [
        "testdata/axset_test_vectors.json",
        "../testdata/axset_test_vectors.json",
        "../../testdata/axset_test_vectors.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(events) = serde_json::from_str::<Vec<AxsetEvent>>(&content) {
                return Some(events);
            }
        }
    }
    None
}

#[test]
fn test_axset_vs_xfoil() {
    let events = match load_axset_events() {
        Some(e) => e,
        None => {
            eprintln!("Skipping test: AXSET test data not found");
            return;
        }
    };

    println!("\n=== AXSET Comparison: RustFoil vs XFOIL ===\n");
    println!("Total test vectors: {}", events.len());

    let mut total_ax1_error = 0.0;
    let mut total_ax2_error = 0.0;
    let mut total_rms_error = 0.0;
    let mut total_dax_error = 0.0;
    let mut total_final_error = 0.0;
    let mut n_compared = 0;

    println!("\n{:>4} {:>8} {:>8} {:>12} {:>12} {:>12} {:>8}",
        "Idx", "Hk1", "Hk2", "XFOIL_Ax", "RUST_Ax", "Diff", "Err%");

    for (i, e) in events.iter().enumerate() {
        let result = axset(
            e.hk1, e.theta1, e.rtheta1, e.ampl1,
            e.hk2, e.theta2, e.rtheta2, e.ampl2,
            e.ncrit,
        );

        // Compare individual values
        let ax1_err = if e.ax1.abs() > 1e-20 {
            (result.ax1 - e.ax1).abs() / e.ax1.abs() * 100.0
        } else if result.ax1.abs() < 1e-20 {
            0.0
        } else {
            100.0
        };

        let ax2_err = if e.ax2.abs() > 1e-20 {
            (result.ax2 - e.ax2).abs() / e.ax2.abs() * 100.0
        } else if result.ax2.abs() < 1e-20 {
            0.0
        } else {
            100.0
        };

        let rms_err = if e.axa_rms.abs() > 1e-20 {
            (result.axa_rms - e.axa_rms).abs() / e.axa_rms.abs() * 100.0
        } else if result.axa_rms.abs() < 1e-20 {
            0.0
        } else {
            100.0
        };

        let dax_err = if e.dax.abs() > 1e-20 {
            (result.dax - e.dax).abs() / e.dax.abs() * 100.0
        } else if result.dax.abs() < 1e-20 {
            0.0
        } else {
            100.0
        };

        let final_err = if e.ax_final.abs() > 1e-20 {
            (result.ax - e.ax_final).abs() / e.ax_final.abs() * 100.0
        } else if result.ax.abs() < 1e-20 {
            0.0
        } else {
            100.0
        };

        total_ax1_error += ax1_err;
        total_ax2_error += ax2_err;
        total_rms_error += rms_err;
        total_dax_error += dax_err;
        total_final_error += final_err;
        n_compared += 1;

        // Print first 20 or any with large error
        if i < 20 || final_err > 5.0 {
            println!(
                "{:4} {:8.4} {:8.4} {:12.4e} {:12.4e} {:12.4e} {:8.2}",
                i, e.hk1, e.hk2, e.ax_final, result.ax, result.ax - e.ax_final, final_err
            );
        }
    }

    if n_compared > 0 {
        let avg_ax1_error = total_ax1_error / n_compared as f64;
        let avg_ax2_error = total_ax2_error / n_compared as f64;
        let avg_rms_error = total_rms_error / n_compared as f64;
        let avg_dax_error = total_dax_error / n_compared as f64;
        let avg_final_error = total_final_error / n_compared as f64;

        println!("\n--- Summary ---");
        println!("  Compared {} test vectors", n_compared);
        println!("  Average AX1 error: {:.2}%", avg_ax1_error);
        println!("  Average AX2 error: {:.2}%", avg_ax2_error);
        println!("  Average RMS error: {:.2}%", avg_rms_error);
        println!("  Average DAX error: {:.2}%", avg_dax_error);
        println!("  Average Final Ax error: {:.2}%", avg_final_error);

        // Assert reasonable accuracy
        assert!(
            avg_final_error < 5.0,
            "Average Ax_final error too high: {:.2}%",
            avg_final_error
        );
    }
}

#[test]
fn test_axset_rms_formula() {
    // Verify the RMS formula matches XFOIL exactly
    // RMS = sqrt(0.5*(ax1² + ax2²))
    
    let events = match load_axset_events() {
        Some(e) => e,
        None => {
            eprintln!("Skipping test: AXSET test data not found");
            return;
        }
    };

    println!("\n=== RMS Formula Verification ===\n");
    
    for (i, e) in events.iter().take(10).enumerate() {
        // Compute expected RMS from XFOIL's AX1 and AX2
        let axsq = 0.5 * (e.ax1 * e.ax1 + e.ax2 * e.ax2);
        let expected_rms = if axsq > 0.0 { axsq.sqrt() } else { 0.0 };
        
        let diff = (expected_rms - e.axa_rms).abs();
        
        println!(
            "[{}] XFOIL: AX1={:.4e}, AX2={:.4e}, RMS={:.4e}, computed={:.4e}, diff={:.4e}",
            i, e.ax1, e.ax2, e.axa_rms, expected_rms, diff
        );
        
        assert!(
            diff < 1e-15,
            "RMS formula mismatch at {}: XFOIL={:.4e}, computed={:.4e}",
            i, e.axa_rms, expected_rms
        );
    }
}
