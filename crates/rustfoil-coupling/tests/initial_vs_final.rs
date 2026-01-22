//! Compare using INITIAL vs FINAL Newton values for TRCHEK2

use rustfoil_bl::closures::{amplification_rate, axset};

fn compute_rcrit(hk: f64) -> f64 {
    let hmi = 1.0 / (hk - 1.0);
    let aa = 2.492 * hmi.powf(0.43);
    let bb = (14.0 * hmi - 9.24).tanh();
    let grcrit = aa + 0.7 * (bb + 1.0);
    10.0_f64.powf(grcrit)
}

#[test]
fn test_initial_vs_final_values() {
    // Data from mrchue_iterations.json
    // Station 17 FINAL (which becomes station 18 INITIAL)
    let st17_final = (6.1254812e-5_f64, 1.5714112e-4_f64, 118.60_f64);  // theta, delta_star, Rt
    let st17_h = st17_final.1 / st17_final.0;  // H = delta_star / theta
    let st17_hk = st17_h;  // At M=0, Hk = H
    
    // Station 18 FINAL (after Newton converges)
    let st18_final = (6.8005125e-5_f64, 1.7834412e-4_f64, 130.49_f64);  // theta, delta_star, Rt  
    let st18_h_final = st18_final.1 / st18_final.0;
    let st18_hk_final = st18_h_final;
    
    // Station 18 INITIAL = Station 17 FINAL (before Newton runs)
    let st18_initial = st17_final;
    let st18_h_initial = st18_initial.1 / st18_initial.0;
    let st18_hk_initial = st18_h_initial;
    
    println!("\n=== Station 17 FINAL ===");
    println!("theta={:.6e}, delta_star={:.6e}", st17_final.0, st17_final.1);
    println!("H = {:.4}, Hk = {:.4}", st17_h, st17_hk);
    let rcrit17 = compute_rcrit(st17_hk);
    println!("Rcrit = {:.1}, Rt = {:.1}, Rt/Rcrit = {:.4}", 
             rcrit17, st17_final.2, st17_final.2 / rcrit17);
    let ax17 = amplification_rate(st17_hk, st17_final.0, st17_final.2);
    println!("AX = {:.6e}", ax17.ax);
    
    println!("\n=== Station 18 INITIAL (before Newton) ===");
    println!("theta={:.6e}, delta_star={:.6e}", st18_initial.0, st18_initial.1);
    println!("H = {:.4}, Hk = {:.4}", st18_h_initial, st18_hk_initial);
    // For station 18 initial, Rt would be similar (Re * Ue * theta) 
    // but Ue might be slightly different. Approximate as same Rt for now.
    let rcrit18_init = compute_rcrit(st18_hk_initial);
    let rt18_init = st17_final.2;  // Approximately same
    println!("Rcrit = {:.1}, Rt ≈ {:.1}, Rt/Rcrit = {:.4}", 
             rcrit18_init, rt18_init, rt18_init / rcrit18_init);
    let ax18_init = amplification_rate(st18_hk_initial, st18_initial.0, rt18_init);
    println!("AX = {:.6e}", ax18_init.ax);
    
    println!("\n=== Station 18 FINAL (after Newton) ===");
    println!("theta={:.6e}, delta_star={:.6e}", st18_final.0, st18_final.1);
    println!("H = {:.4}, Hk = {:.4}", st18_h_final, st18_hk_final);
    let rcrit18_final = compute_rcrit(st18_hk_final);
    println!("Rcrit = {:.1}, Rt = {:.1}, Rt/Rcrit = {:.4}", 
             rcrit18_final, st18_final.2, st18_final.2 / rcrit18_final);
    let ax18_final = amplification_rate(st18_hk_final, st18_final.0, st18_final.2);
    println!("AX = {:.6e}", ax18_final.ax);
    
    println!("\n=== AXSET Comparison ===");
    let ncrit = 9.0;
    let ampl1 = 1.77e-9;
    let dx = 0.002961;
    
    println!("\nUsing INITIAL values (st17_final + st18_initial):");
    let axset_init = axset(st17_hk, st17_final.0, st17_final.2, ampl1,
                           st18_hk_initial, st18_initial.0, rt18_init, ampl1, ncrit);
    println!("  AX1={:.6e}, AX2={:.6e}, RMS={:.6e}", axset_init.ax1, axset_init.ax2, axset_init.axa_rms);
    println!("  DAX={:.6e}, Final AX={:.6e}", axset_init.dax, axset_init.ax);
    let dn_init = axset_init.ax * dx;
    println!("  dN = {:.6e}", dn_init);
    
    println!("\nUsing FINAL values (st17_final + st18_final):");
    let axset_final = axset(st17_hk, st17_final.0, st17_final.2, ampl1,
                            st18_hk_final, st18_final.0, st18_final.2, ampl1, ncrit);
    println!("  AX1={:.6e}, AX2={:.6e}, RMS={:.6e}", axset_final.ax1, axset_final.ax2, axset_final.axa_rms);
    println!("  DAX={:.6e}, Final AX={:.6e}", axset_final.dax, axset_final.ax);
    let dn_final = axset_final.ax * dx;
    println!("  dN = {:.6e}", dn_final);
    
    println!("\n=== Summary ===");
    println!("XFOIL actual dN: 9.4e-11");
    println!("Our dN with INITIAL values: {:.4e}", dn_init);
    println!("Our dN with FINAL values: {:.4e}", dn_final);
    println!("\nRatio (final/xfoil): {:.0}", dn_final / 9.4e-11);
    println!("Ratio (init/xfoil): {:.1}", dn_init / 9.4e-11);
}
