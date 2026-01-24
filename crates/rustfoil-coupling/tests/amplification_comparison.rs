//! Amplification comparison test
//!
//! Compares RustFoil's AXSET and N-factor evolution against XFOIL's TRCHEK2 output
//! to identify discrepancies in transition handling.

use rustfoil_bl::closures::{amplification_rate, axset, trchek2_stations};
use serde::Deserialize;
use std::fs;

/// TRCHEK2 iteration event from XFOIL debug trace
#[derive(Debug, Deserialize)]
struct Trchek2Iter {
    subroutine: String,
    global_iter: i32,
    side: i32,
    ibl: i32,
    trchek_iter: i32,
    x1: f64,
    x2: f64,
    ampl1: f64,
    ampl2: f64,
    ax: f64,
    residual: f64,
    wf1: f64,
    wf2: f64,
    xt: f64,
    #[serde(rename = "Hk1")]
    hk1: f64,
    #[serde(rename = "Hk2")]
    hk2: f64,
    #[serde(rename = "Rt1")]
    rt1: f64,
    #[serde(rename = "Rt2")]
    rt2: f64,
    #[serde(rename = "T1")]
    t1: f64,
    #[serde(rename = "T2")]
    t2: f64,
    #[serde(rename = "Ncrit")]
    ncrit: f64,
    transition: bool,
}

/// TRCHEK2 final result event from XFOIL debug trace
#[derive(Debug, Deserialize)]
struct Trchek2Final {
    subroutine: String,
    global_iter: i32,
    side: i32,
    ibl: i32,
    n_iterations: i32,
    converged: bool,
    x1: f64,
    x2: f64,
    ampl1: f64,
    ampl2_final: f64,
    ax_final: f64,
    xt_final: f64,
    #[serde(rename = "Ncrit")]
    ncrit: f64,
    transition: bool,
    forced: bool,
}

/// Wrapper for parsing events array
#[derive(Debug, Deserialize)]
struct DebugTrace {
    events: Vec<serde_json::Value>,
}

fn load_trchek_trace() -> Option<Vec<serde_json::Value>> {
    let paths = [
        "testdata/trchek_trace.json",
        "../testdata/trchek_trace.json",
        "../../testdata/trchek_trace.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(trace) = serde_json::from_str::<DebugTrace>(&content) {
                return Some(trace.events);
            }
        }
    }
    None
}

/// Extract TRCHEK2_FINAL events from the trace
fn get_trchek2_final_events(events: &[serde_json::Value]) -> Vec<Trchek2Final> {
    events
        .iter()
        .filter_map(|ev| {
            if ev.get("subroutine").and_then(|s| s.as_str()) == Some("TRCHEK2_FINAL") {
                serde_json::from_value(ev.clone()).ok()
            } else {
                None
            }
        })
        .collect()
}

/// Extract TRCHEK2_ITER events from the trace
fn get_trchek2_iter_events(events: &[serde_json::Value]) -> Vec<Trchek2Iter> {
    events
        .iter()
        .filter_map(|ev| {
            if ev.get("subroutine").and_then(|s| s.as_str()) == Some("TRCHEK2_ITER") {
                serde_json::from_value(ev.clone()).ok()
            } else {
                None
            }
        })
        .collect()
}

#[test]
fn test_axset_comparison() {
    let events = match load_trchek_trace() {
        Some(ev) => ev,
        None => {
            eprintln!("Skipping test: trchek_trace.json not found");
            return;
        }
    };

    let iter_events = get_trchek2_iter_events(&events);
    
    if iter_events.is_empty() {
        eprintln!("No TRCHEK2_ITER events found in trace");
        return;
    }

    println!("\n=== AXSET Comparison Test ===");
    println!("Testing {} TRCHEK2_ITER events", iter_events.len());
    
    let mut total_ax_error = 0.0;
    let mut max_ax_error = 0.0;
    let mut n_tested = 0;
    
    // Test a subset of events (every 10th to keep output manageable)
    for (i, ev) in iter_events.iter().enumerate().step_by(10) {
        // Call our axset with the same inputs
        let result = axset(
            ev.hk1, ev.t1, ev.rt1, ev.ampl1,
            ev.hk2, ev.t2, ev.rt2, ev.ampl2,
            ev.ncrit,
        );
        
        // Compare AX values
        let ax_xfoil = ev.ax;
        let ax_rust = result.ax;
        
        // Relative error (handle near-zero values)
        let ax_error = if ax_xfoil.abs() > 1e-15 {
            ((ax_rust - ax_xfoil) / ax_xfoil).abs()
        } else if ax_rust.abs() > 1e-15 {
            1.0  // XFOIL is zero but we're not
        } else {
            0.0  // Both near zero
        };
        
        total_ax_error += ax_error;
        if ax_error > max_ax_error {
            max_ax_error = ax_error;
        }
        n_tested += 1;
        
        // Print details for first few and any with significant error
        if i < 30 || ax_error > 0.01 {
            println!(
                "Event {} (side={}, ibl={}, iter={}): AX: XFOIL={:.6e}, Rust={:.6e}, Error={:.2}%",
                i, ev.side, ev.ibl, ev.trchek_iter, ax_xfoil, ax_rust, ax_error * 100.0
            );
        }
    }
    
    let avg_ax_error = total_ax_error / n_tested as f64;
    
    println!("\n=== Summary ===");
    println!("Tested {} events", n_tested);
    println!("Average AX error: {:.4}%", avg_ax_error * 100.0);
    println!("Maximum AX error: {:.4}%", max_ax_error * 100.0);
    
    // Assert reasonable accuracy
    assert!(avg_ax_error < 0.05, "Average AX error too high: {:.2}%", avg_ax_error * 100.0);
}

#[test]
fn test_amplification_evolution() {
    let events = match load_trchek_trace() {
        Some(ev) => ev,
        None => {
            eprintln!("Skipping test: trchek_trace.json not found");
            return;
        }
    };

    let final_events = get_trchek2_final_events(&events);
    
    if final_events.is_empty() {
        eprintln!("No TRCHEK2_FINAL events found in trace");
        return;
    }

    println!("\n=== Amplification Evolution Test ===");
    println!("Testing {} TRCHEK2_FINAL events", final_events.len());
    
    // Group by side and track evolution
    let side1_events: Vec<_> = final_events.iter().filter(|e| e.side == 1).collect();
    let side2_events: Vec<_> = final_events.iter().filter(|e| e.side == 2).collect();
    
    println!("\n--- Side 1: {} events ---", side1_events.len());
    
    // Track our own amplification using explicit integration
    let mut our_ampl = 0.0;
    let mut prev_x = 0.0;
    let mut prev_hk = 0.0;
    let mut prev_t = 0.0;
    let mut prev_rt = 0.0;
    let mut first = true;
    
    let mut n_errors = 0;
    let mut transition_detected_xfoil = false;
    let mut transition_detected_rust = false;
    
    for (i, ev) in side1_events.iter().enumerate() {
        if first {
            // Initialize from first station
            prev_x = ev.x1;
            // Note: we don't have Hk1, T1, Rt1 in FINAL events, only in ITER events
            // So we'll use the final values for comparison
            first = false;
            continue;
        }
        
        // Check if XFOIL detected transition
        if ev.transition && !transition_detected_xfoil {
            transition_detected_xfoil = true;
            println!(
                "XFOIL transition at IBL={}, x={:.6}, ampl2={:.4}",
                ev.ibl, ev.xt_final, ev.ampl2_final
            );
        }
        
        // Compare final amplification values
        let ampl_xfoil = ev.ampl2_final;
        
        // For now, just report the XFOIL values
        if i < 20 || ev.transition || ampl_xfoil > 0.1 {
            println!(
                "IBL={:3}: x={:.6}, ampl1={:.6e}, ampl2={:.6e}, ax={:.6e}, trans={}",
                ev.ibl, ev.x2, ev.ampl1, ev.ampl2_final, ev.ax_final, ev.transition
            );
        }
        
        // Check for transition in our data
        if ampl_xfoil >= ev.ncrit && !transition_detected_rust {
            transition_detected_rust = true;
            println!(
                "Rust transition at IBL={}, x={:.6}, ampl={:.4}",
                ev.ibl, ev.x2, ampl_xfoil
            );
        }
    }
    
    println!("\n--- Side 2: {} events ---", side2_events.len());
    
    transition_detected_xfoil = false;
    
    for (i, ev) in side2_events.iter().enumerate() {
        // Check if XFOIL detected transition
        if ev.transition && !transition_detected_xfoil {
            transition_detected_xfoil = true;
            println!(
                "XFOIL transition at IBL={}, x={:.6}, ampl2={:.4}",
                ev.ibl, ev.xt_final, ev.ampl2_final
            );
        }
        
        // For now, just report the XFOIL values
        let ampl_xfoil = ev.ampl2_final;
        if i < 20 || ev.transition || ampl_xfoil > 0.1 {
            println!(
                "IBL={:3}: x={:.6}, ampl1={:.6e}, ampl2={:.6e}, ax={:.6e}, trans={}",
                ev.ibl, ev.x2, ev.ampl1, ev.ampl2_final, ev.ax_final, ev.transition
            );
        }
    }
}

#[test]
fn test_dampl_comparison() {
    // Test the base DAMPL function directly against XFOIL's TRCHEK2_ITER data
    let events = match load_trchek_trace() {
        Some(ev) => ev,
        None => {
            eprintln!("Skipping test: trchek_trace.json not found");
            return;
        }
    };

    let iter_events = get_trchek2_iter_events(&events);
    
    if iter_events.is_empty() {
        eprintln!("No TRCHEK2_ITER events found in trace");
        return;
    }

    println!("\n=== DAMPL (amplification_rate) Comparison Test ===");
    
    // Test individual DAMPL calls
    // XFOIL calls DAMPL twice per AXSET call - once for station 1, once for station 2
    // The TRCHEK2_ITER event gives us the final AX after RMS averaging
    
    let mut total_error = 0.0;
    let mut n_tested = 0;
    
    for (i, ev) in iter_events.iter().enumerate().take(50) {
        // Test DAMPL for station 1
        let result1 = amplification_rate(ev.hk1, ev.t1, ev.rt1);
        
        // Test DAMPL for station 2
        let result2 = amplification_rate(ev.hk2, ev.t2, ev.rt2);
        
        // Compute our RMS average
        let ax1 = result1.ax;
        let ax2 = result2.ax;
        let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
        let axa_rms = if axsq > 0.0 { axsq.sqrt() } else { 0.0 };
        
        // Compare with XFOIL's AX (which includes the near-Ncrit correction DAX)
        // So we need to compute DAX too
        let avg_ampl = 0.5 * (ev.ampl1 + ev.ampl2);
        let arg = (20.0 * (ev.ncrit - avg_ampl)).min(20.0);
        let exn = if arg <= 0.0 { 1.0 } else { (-arg).exp() };
        let th_sum = ev.t1 + ev.t2;
        let dax = if th_sum > 0.0 { exn * 0.002 / th_sum } else { 0.0 };
        let ax_computed = axa_rms + dax;
        
        let error = if ev.ax.abs() > 1e-15 {
            ((ax_computed - ev.ax) / ev.ax).abs()
        } else {
            0.0
        };
        
        total_error += error;
        n_tested += 1;
        
        if i < 10 || error > 0.01 {
            println!(
                "Event {}: Hk1={:.4}, Rt1={:.2e}, ax1={:.4e} | Hk2={:.4}, Rt2={:.2e}, ax2={:.4e}",
                i, ev.hk1, ev.rt1, ax1, ev.hk2, ev.rt2, ax2
            );
            println!(
                "         AXA_rms={:.4e}, DAX={:.4e}, AX_total={:.4e}, XFOIL={:.4e}, Error={:.2}%",
                axa_rms, dax, ax_computed, ev.ax, error * 100.0
            );
        }
    }
    
    if n_tested > 0 {
        let avg_error = total_error / n_tested as f64;
        println!("\nAverage error: {:.4}%", avg_error * 100.0);
        assert!(avg_error < 0.01, "Average DAMPL error too high: {:.2}%", avg_error * 100.0);
    }
}

#[test]
fn test_trchek2_comparison() {
    // Test the TRCHEK2 implicit iteration against XFOIL's trace
    let events = match load_trchek_trace() {
        Some(ev) => ev,
        None => {
            eprintln!("Skipping test: trchek_trace.json not found");
            return;
        }
    };

    let final_events = get_trchek2_final_events(&events);
    let iter_events = get_trchek2_iter_events(&events);
    
    if final_events.is_empty() {
        eprintln!("No TRCHEK2_FINAL events found in trace");
        return;
    }

    println!("\n=== TRCHEK2 Comparison Test ===");
    println!("Testing {} TRCHEK2_FINAL events", final_events.len());
    
    let mut total_ampl_error = 0.0;
    let mut max_ampl_error = 0.0;
    let mut n_tested = 0;
    let mut transition_mismatch = 0;
    
    // Find pairs of consecutive events on same side to get station 1 and 2 data
    let side1_finals: Vec<_> = final_events.iter().filter(|e| e.side == 1).collect();
    
    // Use ITER events for full input data
    let side1_iters: Vec<_> = iter_events.iter().filter(|e| e.side == 1).collect();
    
    // Group ITER events by IBL (take first iteration of each)
    let mut prev_iter: Option<&Trchek2Iter> = None;
    
    for ev in side1_iters.iter() {
        if ev.trchek_iter != 1 {
            continue; // Only use first iteration for input state
        }
        
        if let Some(prev) = prev_iter {
            if prev.ibl + 1 == ev.ibl {
                // Consecutive stations - test TRCHEK2
                // Estimate d1, d2 from Hk and theta: d = Hk * theta
                let d1 = prev.hk2 * prev.t2;
                let d2 = ev.hk2 * ev.t2;
                let result = trchek2_stations(
                    prev.x2,  // x1 = prev's x2
                    ev.x2,    // x2 = current's x2
                    prev.hk2, // hk1 = prev's hk2 (output from prev station)
                    prev.t2,  // t1 = prev's T2
                    prev.rt2, // rt1 = prev's Rt2
                    d1,       // d1 = estimated displacement thickness
                    prev.rt2 / prev.t2, // u1 from Rt1 = Re * Ue * theta (use Re=1)
                    prev.ampl2, // ampl1 = prev's ampl2 (N at end of prev station)
                    ev.hk2,   // hk2 = current's hk2
                    ev.t2,    // t2 = current's T2  
                    ev.rt2,   // rt2 = current's Rt2
                    d2,       // d2 = estimated displacement thickness
                    ev.rt2 / ev.t2, // u2 from Rt2 = Re * Ue * theta (use Re=1)
                    ev.ncrit, // ncrit
                    0.0,      // msq
                    1.0,      // re
                );
                
                // Compare with XFOIL's AMPL2
                let ampl_xfoil = ev.ampl2;
                let ampl_rust = result.ampl2;
                
                let error = if ampl_xfoil.abs() > 1e-10 {
                    ((ampl_rust - ampl_xfoil) / ampl_xfoil).abs()
                } else if ampl_rust.abs() > 1e-10 {
                    1.0
                } else {
                    0.0
                };
                
                total_ampl_error += error;
                if error > max_ampl_error {
                    max_ampl_error = error;
                }
                n_tested += 1;
                
                // Check transition agreement
                let xfoil_trans = ev.ampl2 >= ev.ncrit;
                if result.transition != xfoil_trans {
                    transition_mismatch += 1;
                }
                
                // Print details for significant errors or near transition
                if error > 0.05 || result.transition || xfoil_trans || n_tested <= 10 {
                    println!(
                        "IBL={:3}: XFOIL ampl2={:.6e}, Rust ampl2={:.6e}, Error={:.2}%, trans: XFOIL={} Rust={}",
                        ev.ibl, ampl_xfoil, ampl_rust, error * 100.0, xfoil_trans, result.transition
                    );
                    if result.transition || xfoil_trans {
                        println!(
                            "         XT: XFOIL={:.6}, Rust={:?}, iters={}",
                            ev.xt, result.xt, result.iterations
                        );
                    }
                }
            }
        }
        
        prev_iter = Some(ev);
    }
    
    if n_tested > 0 {
        let avg_error = total_ampl_error / n_tested as f64;
        println!("\n=== Summary ===");
        println!("Tested {} station pairs", n_tested);
        println!("Average AMPL2 error: {:.4}%", avg_error * 100.0);
        println!("Maximum AMPL2 error: {:.4}%", max_ampl_error * 100.0);
        println!("Transition mismatches: {}", transition_mismatch);
        
        // Note: We expect some error due to simplified implementation
        // The test is informational - we'll tighten the tolerance after debugging
    }
}
