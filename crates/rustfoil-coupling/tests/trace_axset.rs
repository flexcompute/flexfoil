//! Detailed AXSET trace for stations 17-18-19
//!
//! Traces through the AXSET calculation to find where our
//! N-factor integration differs from XFOIL.

use rustfoil_bl::closures::{amplification_rate, axset};

fn compute_rcrit(hk: f64) -> f64 {
    let hmi = 1.0 / (hk - 1.0);
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let grcrit = aa + 0.7 * (bb + 1.0);
    10.0_f64.powf(grcrit)
}

#[test]
fn test_trace_axset_17_to_19() {
    // XFOIL data from mrchue_iterations.json
    let stations = [
        (17, 0.034705, 2.6225_f64, 6.125481e-5_f64, 118.60_f64, 1.771892e-9_f64),
        (18, 0.037666, 2.6815, 6.800512e-5, 130.49, 1.866317e-9),
        (19, 0.040885, 2.7424, 7.533508e-5, 142.93, 9.204327e-4),
    ];
    let ncrit = 9.0;
    
    println!("\n=== Rcrit Analysis ===");
    for (ibl, x, hk, theta, rtheta, ampl) in &stations {
        let rcrit = compute_rcrit(*hk);
        let ratio = rtheta / rcrit;
        let in_ramp = ratio > 0.6 && ratio < 1.0;  // Approximate ramp region
        println!("IBL={}: Rtheta={:.2}, Rcrit={:.2}, ratio={:.4} {}",
                 ibl, rtheta, rcrit, ratio, if in_ramp { "[RAMP]" } else if ratio >= 1.0 { "[SUPER]" } else { "" });
    }
    
    println!("\n=== DAMPL at each station ===");
    for (ibl, _x, hk, theta, rtheta, _ampl) in &stations {
        let result = amplification_rate(*hk, *theta, *rtheta);
        println!("IBL={}: AX={:.6e}", ibl, result.ax);
    }
    
    println!("\n=== AXSET Trace: Station 17 -> 18 ===");
    let (ibl1, x1, hk1, t1, rt1, ampl1_xfoil) = stations[0];
    let (ibl2, x2, hk2, t2, rt2, ampl2_xfoil) = stations[1];
    
    let dx = x2 - x1;
    println!("dx = {:.6}", dx);
    
    // Test with ampl1 = XFOIL's value (nearly 0)
    let axset_result = axset(hk1, t1, rt1, ampl1_xfoil, hk2, t2, rt2, ampl1_xfoil, ncrit);
    println!("\nUsing ampl1 = XFOIL's value ({:.6e}):", ampl1_xfoil);
    println!("  AX1 = {:.6e}", axset_result.ax1);
    println!("  AX2 = {:.6e}", axset_result.ax2);
    println!("  AXA_RMS = {:.6e}", axset_result.axa_rms);
    println!("  DAX = {:.6e}", axset_result.dax);
    println!("  AX (final) = {:.6e}", axset_result.ax);
    
    let dn_ours = axset_result.ax * dx;
    let dn_xfoil = ampl2_xfoil - ampl1_xfoil;
    println!("\n  Our dN = AX * dx = {:.6e}", dn_ours);
    println!("  XFOIL dN = {:.6e}", dn_xfoil);
    println!("  Ratio (ours/xfoil) = {:.1}", dn_ours / dn_xfoil);
    
    println!("\n=== AXSET Trace: Station 18 -> 19 ===");
    let (ibl1, x1, hk1, t1, rt1, ampl1_xfoil) = stations[1];
    let (ibl2, x2, hk2, t2, rt2, ampl2_xfoil) = stations[2];
    
    let dx = x2 - x1;
    println!("dx = {:.6}", dx);
    
    // Test with ampl1 = XFOIL's value at station 18
    let axset_result = axset(hk1, t1, rt1, ampl1_xfoil, hk2, t2, rt2, ampl1_xfoil, ncrit);
    println!("\nUsing ampl1 = XFOIL's value ({:.6e}):", ampl1_xfoil);
    println!("  AX1 = {:.6e}", axset_result.ax1);
    println!("  AX2 = {:.6e}", axset_result.ax2);
    println!("  AXA_RMS = {:.6e}", axset_result.axa_rms);
    println!("  DAX = {:.6e}", axset_result.dax);
    println!("  AX (final) = {:.6e}", axset_result.ax);
    
    let dn_ours = axset_result.ax * dx;
    let dn_xfoil = ampl2_xfoil - ampl1_xfoil;
    println!("\n  Our dN = AX * dx = {:.6e}", dn_ours);
    println!("  XFOIL dN = {:.6e}", dn_xfoil);
    println!("  Ratio (ours/xfoil) = {:.1}", dn_ours / dn_xfoil);
    
    // Now check what happens if we START with ampl=0 and accumulate
    println!("\n=== Accumulated N-factor starting from 0 ===");
    let mut our_ampl = 0.0;
    for i in 1..stations.len() {
        let (_, x1, hk1, t1, rt1, _) = stations[i-1];
        let (ibl2, x2, hk2, t2, rt2, ampl2_xfoil) = stations[i];
        let dx = x2 - x1;
        
        let axset_result = axset(hk1, t1, rt1, our_ampl, hk2, t2, rt2, our_ampl, ncrit);
        let dn = axset_result.ax * dx;
        our_ampl += dn;
        
        println!("IBL={}: XFOIL_ampl={:.6e}, our_ampl={:.6e}, ratio={:.1}",
                 ibl2, ampl2_xfoil, our_ampl, our_ampl / ampl2_xfoil.max(1e-20));
    }
}
