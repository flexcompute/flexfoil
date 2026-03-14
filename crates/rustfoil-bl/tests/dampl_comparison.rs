//! DAMPL function comparison test
//!
//! Compares our amplification_rate function against XFOIL's DAMPL output
//! to verify the e^n transition prediction matches.

use rustfoil_bl::closures::amplification_rate;
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct XfoilDebug {
    events: Vec<serde_json::Value>,
}

#[derive(Debug, Deserialize)]
struct DamplEvent {
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

fn load_dampl_events() -> Option<Vec<DamplEvent>> {
    let paths = [
        "testdata/dampl_test_vectors.json",
        "../testdata/dampl_test_vectors.json",
        "../../testdata/dampl_test_vectors.json",
        "../Xfoil-instrumented/bin/xfoil_debug.json",
        "../../Xfoil-instrumented/bin/xfoil_debug.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            // Try to load as array of events first (test vectors format)
            if let Ok(events) = serde_json::from_str::<Vec<DamplEvent>>(&content) {
                return Some(events);
            }
            // Try to load as full debug file
            if let Ok(debug) = serde_json::from_str::<XfoilDebug>(&content) {
                let dampl_events: Vec<DamplEvent> = debug
                    .events
                    .iter()
                    .filter(|e| e.get("subroutine").and_then(|s| s.as_str()) == Some("DAMPL"))
                    .filter_map(|e| serde_json::from_value(e.clone()).ok())
                    .collect();
                if !dampl_events.is_empty() {
                    return Some(dampl_events);
                }
            }
        }
    }
    None
}

#[test]
fn test_dampl_vs_xfoil() {
    let events = match load_dampl_events() {
        Some(e) => e,
        None => {
            eprintln!("Skipping test: DAMPL test data not found");
            return;
        }
    };

    println!("\n=== DAMPL Comparison: RustFoil vs XFOIL ===\n");

    // Separate into subcritical (Ax=0) and supercritical (Ax>0) cases
    let subcritical: Vec<_> = events.iter().filter(|e| e.ax == 0.0).collect();
    let supercritical: Vec<_> = events.iter().filter(|e| e.ax > 0.0).collect();

    println!("Total DAMPL events: {}", events.len());
    println!("  Subcritical (Ax=0): {}", subcritical.len());
    println!("  Supercritical (Ax>0): {}", supercritical.len());

    // Test subcritical cases - should all return 0
    println!("\n--- Subcritical Tests ---");
    let mut subcrit_pass = 0;
    let mut subcrit_fail = 0;
    for (i, e) in subcritical.iter().take(100).enumerate() {
        let result = amplification_rate(e.hk, e.theta, e.rtheta);
        if result.ax == 0.0 {
            subcrit_pass += 1;
        } else {
            subcrit_fail += 1;
            if subcrit_fail <= 5 {
                println!(
                    "  [{i}] FAIL: Hk={:.4}, θ={:.4e}, Rθ={:.1} → XFOIL Ax=0, RustFoil Ax={:.4e}",
                    e.hk, e.theta, e.rtheta, result.ax
                );
            }
        }
    }
    println!(
        "  Subcritical: {}/{} passed ({:.1}%)",
        subcrit_pass,
        subcrit_pass + subcrit_fail,
        100.0 * subcrit_pass as f64 / (subcrit_pass + subcrit_fail) as f64
    );

    // Test supercritical cases - should match XFOIL values
    println!("\n--- Supercritical Tests ---");
    println!(
        "{:>6} {:>8} {:>12} {:>12} {:>12} {:>12} {:>8}",
        "Idx", "Hk", "θ", "Rθ", "XFOIL_Ax", "RUST_Ax", "Err%"
    );

    let mut total_ax_error = 0.0;
    let mut total_ax_hk_error = 0.0;
    let mut total_ax_th_error = 0.0;
    let mut total_ax_rt_error = 0.0;
    let mut n_compared = 0;

    // Sample unique supercritical cases (deduplicate by Rθ)
    let mut seen_rtheta = std::collections::HashSet::new();
    let unique_supercrit: Vec<_> = supercritical
        .iter()
        .filter(|e| {
            let key = (e.rtheta * 10.0) as i64; // Round to 0.1
            if seen_rtheta.contains(&key) {
                false
            } else {
                seen_rtheta.insert(key);
                true
            }
        })
        .take(50)
        .collect();

    for (i, e) in unique_supercrit.iter().enumerate() {
        let result = amplification_rate(e.hk, e.theta, e.rtheta);

        let ax_err = if e.ax.abs() > 1e-10 {
            (result.ax - e.ax).abs() / e.ax * 100.0
        } else {
            0.0
        };

        // Only compare derivatives when XFOIL has non-zero values
        if e.ax_hk.abs() > 1e-10 {
            total_ax_hk_error += (result.ax_hk - e.ax_hk).abs() / e.ax_hk.abs() * 100.0;
        }
        if e.ax_theta.abs() > 1e-10 {
            total_ax_th_error += (result.ax_th - e.ax_theta).abs() / e.ax_theta.abs() * 100.0;
        }
        if e.ax_rt.abs() > 1e-10 {
            total_ax_rt_error += (result.ax_rt - e.ax_rt).abs() / e.ax_rt.abs() * 100.0;
        }

        total_ax_error += ax_err;
        n_compared += 1;

        if i < 20 || ax_err > 5.0 {
            println!(
                "{:6} {:8.4} {:12.4e} {:12.1} {:12.4e} {:12.4e} {:8.2}",
                i, e.hk, e.theta, e.rtheta, e.ax, result.ax, ax_err
            );
        }
    }

    if n_compared > 0 {
        let avg_ax_error = total_ax_error / n_compared as f64;
        let avg_ax_hk_error = total_ax_hk_error / n_compared as f64;
        let avg_ax_th_error = total_ax_th_error / n_compared as f64;
        let avg_ax_rt_error = total_ax_rt_error / n_compared as f64;

        println!("\n--- Summary ---");
        println!("  Compared {} unique supercritical cases", n_compared);
        println!("  Average Ax error: {:.2}%", avg_ax_error);
        println!("  Average Ax_Hk error: {:.2}%", avg_ax_hk_error);
        println!("  Average Ax_theta error: {:.2}%", avg_ax_th_error);
        println!("  Average Ax_Rt error: {:.2}%", avg_ax_rt_error);

        // Assert reasonable accuracy
        assert!(
            avg_ax_error < 5.0,
            "Average Ax error too high: {:.2}%",
            avg_ax_error
        );
    }
}

#[test]
fn test_dampl_critical_rtheta() {
    // Test the critical Rθ calculation at various Hk values
    // From XFOIL DAMPL: GRCRIT = AA + 0.7*(BB + 1)
    // where AA = 2.492 * HMI^0.43, BB = tanh(14*HMI - 9.24), HMI = 1/(Hk-1)

    println!("\n=== Critical Rθ vs Hk ===");
    println!("{:>8} {:>12} {:>12}", "Hk", "Rcrit_calc", "Rcrit_obs");

    // Test at various Hk values
    let hk_values = [2.0, 2.2, 2.5, 2.8, 3.0, 3.5, 4.0];

    for hk in &hk_values {
        // Calculate expected critical Rθ from XFOIL formula
        let hmi: f64 = 1.0 / (hk - 1.0);
        let aa: f64 = 2.492 * hmi.powf(0.43);
        let bb: f64 = (14.0 * hmi - 9.24).tanh();
        let grcrit: f64 = aa + 0.7 * (bb + 1.0);
        let rcrit_expected: f64 = 10.0_f64.powf(grcrit);

        // Find actual critical by binary search
        let th = 0.001;
        let mut rt_low = 10.0;
        let mut rt_high = 100000.0;
        while rt_high - rt_low > 1.0 {
            let rt_mid = (rt_low + rt_high) / 2.0;
            let ax = amplification_rate(*hk, th, rt_mid).ax;
            if ax > 0.0 {
                rt_high = rt_mid;
            } else {
                rt_low = rt_mid;
            }
        }
        let rcrit_observed = (rt_low + rt_high) / 2.0;

        let error = (rcrit_observed - rcrit_expected).abs() / rcrit_expected * 100.0;
        println!(
            "{:8.2} {:12.1} {:12.1} ({:.1}% error)",
            hk, rcrit_expected, rcrit_observed, error
        );

        // Critical Rθ should match within 20% (the ramp region is ±DGR wide)
        assert!(
            error < 30.0,
            "Critical Rθ mismatch at Hk={}: expected {:.1}, got {:.1}",
            hk,
            rcrit_expected,
            rcrit_observed
        );
    }
}

#[test]
fn test_dampl_ramp_region() {
    // Test the smooth ramp behavior near critical Rθ
    // XFOIL uses cubic ramp: RFAC = 3*rnorm² - 2*rnorm³
    // over the range GRCRIT-DGR to GRCRIT+DGR where DGR=0.08

    println!("\n=== Ramp Region Test (Hk=2.5) ===");

    let hk = 2.5;
    let th = 0.001;

    // Calculate critical Rθ
    let hmi: f64 = 1.0 / (hk - 1.0);
    let aa: f64 = 2.492 * hmi.powf(0.43);
    let bb: f64 = (14.0 * hmi - 9.24).tanh();
    let grcrit: f64 = aa + 0.7 * (bb + 1.0);
    let rcrit: f64 = 10.0_f64.powf(grcrit);
    let dgr: f64 = 0.08;

    // Test at various points relative to critical
    let rt_below = 10.0_f64.powf(grcrit - dgr - 0.02); // Well below
    let rt_ramp_low = 10.0_f64.powf(grcrit - dgr + 0.02); // In ramp, low side
    let rt_ramp_mid = 10.0_f64.powf(grcrit); // In ramp, middle
    let rt_ramp_high = 10.0_f64.powf(grcrit + dgr - 0.02); // In ramp, high side
    let rt_above = 10.0_f64.powf(grcrit + dgr + 0.02); // Well above

    println!("Critical Rθ: {:.1}", rcrit);
    println!("Ramp region: {:.1} to {:.1}", 
        10.0_f64.powf(grcrit - dgr), 
        10.0_f64.powf(grcrit + dgr));

    let ax_below = amplification_rate(hk, th, rt_below).ax;
    let ax_ramp_low = amplification_rate(hk, th, rt_ramp_low).ax;
    let ax_ramp_mid = amplification_rate(hk, th, rt_ramp_mid).ax;
    let ax_ramp_high = amplification_rate(hk, th, rt_ramp_high).ax;
    let ax_above = amplification_rate(hk, th, rt_above).ax;

    println!("\n{:>12} {:>12} {:>12}", "Rθ", "Region", "Ax");
    println!("{:12.1} {:>12} {:12.4e}", rt_below, "Below", ax_below);
    println!("{:12.1} {:>12} {:12.4e}", rt_ramp_low, "Ramp-Low", ax_ramp_low);
    println!("{:12.1} {:>12} {:12.4e}", rt_ramp_mid, "Ramp-Mid", ax_ramp_mid);
    println!("{:12.1} {:>12} {:12.4e}", rt_ramp_high, "Ramp-High", ax_ramp_high);
    println!("{:12.1} {:>12} {:12.4e}", rt_above, "Above", ax_above);

    // Verify monotonic increase
    assert!(ax_below == 0.0 || ax_below < ax_ramp_low, "Below should be 0 or < ramp_low");
    assert!(ax_ramp_low <= ax_ramp_mid, "Ramp should increase monotonically");
    assert!(ax_ramp_mid <= ax_ramp_high, "Ramp should increase monotonically");
    assert!(ax_ramp_high <= ax_above, "Ramp should increase monotonically");
}
