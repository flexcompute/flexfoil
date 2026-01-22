//! Detailed DAMPL trace for station 18
//!
//! Traces through every step of the amplification_rate calculation
//! to find where it differs from XFOIL.

use rustfoil_bl::closures::amplification_rate;

#[test]
fn test_trace_station18_dampl() {
    // Station 18 values from XFOIL
    let hk: f64 = 2.6815;
    let th: f64 = 6.800512e-5;
    let rt: f64 = 130.49;
    
    println!("\n=== DAMPL Trace for Station 18 ===");
    println!("Input: Hk={}, th={:.6e}, Rt={}", hk, th, rt);
    
    // Step 1: HMI
    let hmi: f64 = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;
    println!("\nStep 1: HMI = 1/(Hk-1) = {:.6}", hmi);
    println!("        HMI_HK = {:.6}", hmi_hk);
    
    // Step 2: AA (log10 critical Rtheta correlation)
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;
    println!("\nStep 2: AA = 2.492 * HMI^0.43 = {:.6}", aa);
    println!("        AA_HK = {:.6}", aa_hk);
    
    // Step 3: BB (tanh term)
    let bb_arg: f64 = 14.0 * hmi - 9.24;
    let bb = bb_arg.tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
    println!("\nStep 3: BB_arg = 14*HMI - 9.24 = {:.6}", bb_arg);
    println!("        BB = tanh(BB_arg) = {:.6}", bb);
    println!("        BB_HK = {:.6}", bb_hk);
    
    // Step 4: GRCRIT (log10 of critical Rtheta)
    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;
    println!("\nStep 4: GRCRIT = AA + 0.7*(BB+1) = {:.6}", grcrit);
    println!("        GRC_HK = {:.6}", grc_hk);
    println!("        Rcrit = 10^GRCRIT = {:.2}", 10.0_f64.powf(grcrit));
    
    // Step 5: GR (log10 of actual Rtheta)
    let gr = rt.log10();
    let gr_rt = 1.0 / (std::f64::consts::LN_10 * rt);
    println!("\nStep 5: GR = log10(Rt) = {:.6}", gr);
    println!("        GR_RT = {:.6}", gr_rt);
    
    // Step 6: Check subcritical
    const DGR: f64 = 0.08;
    let gr_minus_grcrit = gr - grcrit;
    println!("\nStep 6: GR - GRCRIT = {:.6}", gr_minus_grcrit);
    println!("        DGR = {}", DGR);
    println!("        GRCRIT - DGR = {:.6}", grcrit - DGR);
    
    if gr < grcrit - DGR {
        println!("        ==> GR < GRCRIT-DGR: SUBCRITICAL, AX should be 0");
    } else {
        println!("        ==> GR >= GRCRIT-DGR: IN RAMP or SUPERCRITICAL");
    }
    
    // Step 7: RNORM (normalized position in ramp)
    let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);
    println!("\nStep 7: RNORM = (GR - (GRCRIT-DGR)) / (2*DGR) = {:.6}", rnorm);
    println!("        RN_HK = {:.6}", rn_hk);
    println!("        RN_RT = {:.6}", rn_rt);
    
    // Step 8: RFAC (cubic ramp function)
    let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
        println!("\n        RNORM >= 1.0: FULLY SUPERCRITICAL");
        (1.0, 0.0, 0.0)
    } else {
        let rfac = 3.0 * rnorm.powi(2) - 2.0 * rnorm.powi(3);
        let rfac_rn = 6.0 * rnorm - 6.0 * rnorm.powi(2);
        println!("\n        RNORM < 1.0: IN RAMP REGION");
        println!("        RFAC = 3*RNORM^2 - 2*RNORM^3 = {:.6}", rfac);
        println!("        RFAC_RN = {:.6}", rfac_rn);
        (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
    };
    println!("Step 8: RFAC = {:.6}", rfac);
    println!("        RFAC_HK = {:.6}", rfac_hk);
    println!("        RFAC_RT = {:.6}", rfac_rt);
    
    // Step 9: Envelope slope correlation
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;
    let ex = (-arg * arg).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);
    println!("\nStep 9: ARG = 3.87*HMI - 2.52 = {:.6}", arg);
    println!("        EX = exp(-ARG^2) = {:.6}", ex);
    
    // Step 10: DADR
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;
    println!("\nStep 10: DADR = 0.028*(Hk-1) - 0.0345*EX = {:.6}", dadr);
    println!("         DADR_HK = {:.6}", dadr_hk);
    
    // Step 11: AF (new m(H) correlation)
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi.powi(2) + 3.0 * hmi.powi(3);
    let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi.powi(2);
    let af_hk = af_hmi * hmi_hk;
    println!("\nStep 11: AF = -0.05 + 2.7*HMI - 5.5*HMI^2 + 3.0*HMI^3 = {:.6}", af);
    println!("         AF_HK = {:.6}", af_hk);
    
    // Step 12: Final AX
    let ax_base = af * dadr / th;
    let ax = ax_base * rfac;
    println!("\nStep 12: AX_BASE = AF * DADR / TH = {:.6}", ax_base);
    println!("         AX = AX_BASE * RFAC = {:.6}", ax);
    
    // Call our function to verify
    let result = amplification_rate(hk, th, rt);
    println!("\n=== Verification ===");
    println!("Our amplification_rate() returns: AX = {:.6e}", result.ax);
    
    // What XFOIL gives
    println!("\n=== XFOIL Reference ===");
    println!("XFOIL N-factor at station 17: 1.771892e-9");
    println!("XFOIL N-factor at station 18: 1.866317e-9");
    println!("XFOIL dN = 9.4e-11 over dx ≈ 0.003");
    println!("XFOIL effective AX ≈ 3e-8");
    println!("\nOur AX = {:.6e} is {} times larger!", 
             result.ax, result.ax / 3e-8);
}
