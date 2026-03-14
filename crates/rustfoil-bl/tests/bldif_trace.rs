//! Trace through bldif computation step by step

use rustfoil_bl::equations::{blvar, FlowType};
use rustfoil_bl::closures::cf::cf_laminar;
use rustfoil_bl::state::BlStation;

#[test]
fn test_bldif_trace() {
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

    println!("\n=== BLDIF Trace ===");
    println!("\nStation 1: x={:.6}, u={:.6}, θ={:.6e}, δ*={:.6e}",
        s1.x, s1.u, s1.theta, s1.delta_star);
    println!("  H={:.4}, Hk={:.4}, Hs={:.4}, Cf={:.6e}", s1.h, s1.hk, s1.hs, s1.cf);

    println!("\nStation 2: x={:.6}, u={:.6}, θ={:.6e}, δ*={:.6e}",
        s2.x, s2.u, s2.theta, s2.delta_star);
    println!("  H={:.4}, Hk={:.4}, Hs={:.4}, Cf={:.6e}", s2.h, s2.hk, s2.hs, s2.cf);

    // Compute log ratios
    let xlog = (s2.x / s1.x).ln();
    let ulog = (s2.u / s1.u).ln();
    let tlog = (s2.theta / s1.theta).ln();
    let hlog = (s2.hs / s1.hs).ln();

    println!("\n=== Log ratios ===");
    println!("XLOG = ln({:.6}/{:.6}) = {:.6}", s2.x, s1.x, xlog);
    println!("ULOG = ln({:.6}/{:.6}) = {:.6}", s2.u, s1.u, ulog);
    println!("TLOG = ln({:.6e}/{:.6e}) = {:.6e}", s2.theta, s1.theta, tlog);
    println!("HLOG = ln({:.4}/{:.4}) = {:.6}", s2.hs, s1.hs, hlog);

    // Compute CFM at midpoint
    let hka = 0.5 * (s1.hk + s2.hk);
    let rta = 0.5 * (s1.r_theta + s2.r_theta);
    let cfm_result = cf_laminar(hka, rta, msq);
    let cfm = cfm_result.cf;

    println!("\n=== Midpoint Cf (BLMID) ===");
    println!("HKA = 0.5*({:.4} + {:.4}) = {:.4}", s1.hk, s2.hk, hka);
    println!("RTA = 0.5*({:.2} + {:.2}) = {:.2}", s1.r_theta, s2.r_theta, rta);
    println!("CFM = cf_laminar({:.4}, {:.2}) = {:.6e}", hka, rta, cfm);

    // Compute CFX
    let xa = 0.5 * (s1.x + s2.x);
    let ta = 0.5 * (s1.theta + s2.theta);
    let cfx = 0.5 * cfm * xa / ta 
        + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);

    println!("\n=== CFX computation ===");
    println!("XA = {:.6}", xa);
    println!("TA = {:.6e}", ta);
    println!("Term1 = 0.5*{:.6e}*{:.6}/{:.6e} = {:.4}", cfm, xa, ta, 0.5*cfm*xa/ta);
    println!("Term2a = 0.25*{:.6e}*{:.6}/{:.6e} = {:.4}", s1.cf, s1.x, s1.theta, 0.25*s1.cf*s1.x/s1.theta);
    println!("Term2b = 0.25*{:.6e}*{:.6}/{:.6e} = {:.4}", s2.cf, s2.x, s2.theta, 0.25*s2.cf*s2.x/s2.theta);
    println!("CFX = {:.4}", cfx);

    // Compute residual
    let ha = 0.5 * (s1.h + s2.h);
    let btmp = ha + 2.0;

    println!("\n=== Momentum residual ===");
    println!("HA = {:.4}", ha);
    println!("BTMP = HA + 2 = {:.4}", btmp);
    println!("REZT = TLOG + BTMP*ULOG - XLOG*0.5*CFX");
    println!("     = {:.4e} + {:.4}*{:.4} - {:.4}*0.5*{:.4}", tlog, btmp, ulog, xlog, cfx);
    let rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
    println!("     = {:.4e} + {:.4} - {:.4}", tlog, btmp * ulog, xlog * 0.5 * cfx);
    println!("     = {:.4}", rezt);
    println!("res_mom = -REZT = {:.4}", -rezt);
    println!("\nXFOIL IBL=3 res_mom = -0.130");
}
