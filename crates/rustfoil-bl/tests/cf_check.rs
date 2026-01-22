//! Check Cf values computed by blvar

use rustfoil_bl::equations::{blvar, FlowType};
use rustfoil_bl::state::BlStation;

#[test]
fn test_cf_values_at_ibl3_setup() {
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
    s2.theta = 4.121333e-5;  // Initial guess = prev theta
    s2.delta_star = 9.188543e-5;  // Initial guess = prev delta_star
    s2.ampl = 0.0;
    s2.ctau = 0.03;
    s2.is_laminar = true;
    blvar(&mut s2, FlowType::Laminar, msq, re);

    println!("\n=== Station 1 (IBL=2 final, prev) ===");
    println!("  x={:.6}, Ue={:.6}", s1.x, s1.u);
    println!("  θ={:.6e}, δ*={:.6e}", s1.theta, s1.delta_star);
    println!("  H={:.4}, Hk={:.4}", s1.h, s1.hk);
    println!("  Rθ={:.2}", s1.r_theta);
    println!("  Cf={:.6e}", s1.cf);

    println!("\n=== Station 2 (IBL=3 initial, curr) ===");
    println!("  x={:.6}, Ue={:.6}", s2.x, s2.u);
    println!("  θ={:.6e}, δ*={:.6e}", s2.theta, s2.delta_star);
    println!("  H={:.4}, Hk={:.4}", s2.h, s2.hk);
    println!("  Rθ={:.2}", s2.r_theta);
    println!("  Cf={:.6e}", s2.cf);

    // Compute CFX term
    let xa = 0.5 * (s1.x + s2.x);
    let ta = 0.5 * (s1.theta + s2.theta);
    let cfm = 0.5 * (s1.cf + s2.cf);
    let cfx = 0.5 * cfm * xa / ta 
        + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);

    println!("\n=== CFX computation ===");
    println!("  CFM = 0.5*(Cf1+Cf2) = {:.6e}", cfm);
    println!("  xa={:.6}, ta={:.6e}", xa, ta);
    println!("  CFX = {:.4}", cfx);

    // Compute residual
    let xlog = (s2.x / s1.x).ln();
    let ulog = (s2.u / s1.u).ln();
    let tlog = (s2.theta / s1.theta).ln();
    let ha = 0.5 * (s1.h + s2.h);
    let btmp = ha + 2.0;

    let res_mom = -(tlog + btmp * ulog - xlog * 0.5 * cfx);

    println!("\n=== Residual computation ===");
    println!("  xlog={:.4}, ulog={:.4}, tlog={:.4e}", xlog, ulog, tlog);
    println!("  HA={:.4}, BTMP={:.4}", ha, btmp);
    println!("  res_mom = -({:.4e} + {:.4}*{:.4} - {:.4}*0.5*{:.4})", 
        tlog, btmp, ulog, xlog, cfx);
    println!("  res_mom = {:.4}", res_mom);
    println!("\n  XFOIL IBL=3 res_mom = -0.130");
}
