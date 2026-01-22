//! Check midpoint Cf (CFM) computation

use rustfoil_bl::equations::{blvar, bldif, FlowType};
use rustfoil_bl::state::BlStation;

#[test]
fn test_cfm_computation() {
    let re = 1e6;
    let msq = 0.0;

    // Station 1 (prev, IBL=2 final)
    let mut s1 = BlStation::new();
    s1.x = 0.002310;
    s1.u = 0.115336;
    s1.theta = 4.121333e-5;
    s1.delta_star = 9.188543e-5;
    s1.ampl = 0.0;
    s1.ctau = 0.03;
    s1.is_laminar = true;
    blvar(&mut s1, FlowType::Laminar, msq, re);

    // Station 2 (curr, IBL=3 initial)
    let mut s2 = BlStation::new();
    s2.x = 0.004704;
    s2.u = 0.258499;
    s2.theta = 4.121333e-5;
    s2.delta_star = 9.188543e-5;
    s2.ampl = 0.0;
    s2.ctau = 0.03;
    s2.is_laminar = true;
    blvar(&mut s2, FlowType::Laminar, msq, re);

    // Call bldif to get residuals
    let (res, jac) = bldif(&s1, &s2, FlowType::Laminar, msq, re);

    println!("\n=== Station values ===");
    println!("Station 1: Hk={:.4}, Rθ={:.2}, Cf={:.6e}", s1.hk, s1.r_theta, s1.cf);
    println!("Station 2: Hk={:.4}, Rθ={:.2}, Cf={:.6e}", s2.hk, s2.r_theta, s2.cf);

    println!("\n=== bldif Residuals ===");
    println!("res_mom = {:.6}", res.res_mom);
    println!("XFOIL res_mom = -0.130");

    // Compute our own CFX for comparison
    // Our blmid averages Hk and Rθ, then computes Cf at those midpoints
    let hka = 0.5 * (s1.hk + s2.hk);
    let rta = 0.5 * (s1.r_theta + s2.r_theta);
    
    // Laminar Cf closure at midpoint
    // This is what blmid should compute
    let cfm = laminar_cf_approx(hka, rta);
    
    println!("\n=== Midpoint analysis ===");
    println!("HKA (avg Hk) = {:.4}", hka);
    println!("RTA (avg Rθ) = {:.2}", rta);
    println!("Simple avg Cf = {:.6e}", 0.5*(s1.cf + s2.cf));
    println!("CFM at midpoint Hk,Rθ ≈ {:.6e}", cfm);

    // Compute CFX
    let xa = 0.5 * (s1.x + s2.x);
    let ta = 0.5 * (s1.theta + s2.theta);
    let cfx_simple = 0.5 * (0.5*(s1.cf + s2.cf)) * xa / ta 
        + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);
    let cfx_midpt = 0.5 * cfm * xa / ta 
        + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);

    println!("\n=== CFX comparison ===");
    println!("CFX (simple avg) = {:.4}", cfx_simple);
    println!("CFX (midpoint Cf) = {:.4}", cfx_midpt);

    // Compute residual
    let xlog = (s2.x / s1.x).ln();
    let ulog = (s2.u / s1.u).ln();
    let tlog = (s2.theta / s1.theta).ln();
    let ha = 0.5 * (s1.h + s2.h);
    let btmp = ha + 2.0;

    let res_simple = -(tlog + btmp * ulog - xlog * 0.5 * cfx_simple);
    let res_midpt = -(tlog + btmp * ulog - xlog * 0.5 * cfx_midpt);

    println!("\n=== Residual comparison ===");
    println!("res_mom (simple avg) = {:.4}", res_simple);
    println!("res_mom (midpoint Cf) = {:.4}", res_midpt);
    println!("res_mom (bldif) = {:.4}", res.res_mom);
    println!("res_mom (XFOIL) = -0.130");
}

// Approximate laminar Cf formula
fn laminar_cf_approx(hk: f64, rt: f64) -> f64 {
    // Falkner-Skan based formula from XFOIL
    let fc = (hk - 1.0).powi(2) / ((hk + 1.0) * (hk - 0.9).max(0.1));
    0.3 * fc / rt.max(0.1).sqrt()
}
