//! Direct comparison test between RustFoil and XFOIL at transition
//!
//! This test uses exact XFOIL values from the debug output to validate
//! RustFoil's TRDIF implementation.

use rustfoil_bl::closures::{trchek2_full, trchek2_stations};
use rustfoil_bl::equations::{blvar, trdif, trdif_full, FlowType};
use rustfoil_bl::state::BlStation;

/// XFOIL transition station exact values
/// From NACA 0012, alpha=4°, Re=1e6 debug output
#[test]
fn test_xfoil_transition_comparison() {
    // ===================================================================
    // XFOIL Station 37 (laminar, just before transition)
    // From xfoil_debug.json MRCHUE at ibl=37
    // ===================================================================
    let mut s37 = BlStation::new();
    s37.x = 0.22295522;
    s37.u = 1.3352973;
    s37.theta = 0.32148343e-3;
    s37.delta_star = 0.10817976e-2;
    s37.ctau = 0.03;  // Laminar default
    s37.ampl = 8.2344856;
    s37.is_laminar = true;
    s37.is_turbulent = false;
    
    let re = 1e6;
    let msq = 0.0;
    let ncrit = 9.0;
    
    blvar(&mut s37, FlowType::Laminar, msq, re);
    
    println!("\n=== XFOIL Station 37 (laminar, pre-transition) ===");
    println!("x = {:.6}", s37.x);
    println!("theta = {:.6e} (XFOIL: 3.2148343e-4)", s37.theta);
    println!("delta* = {:.6e} (XFOIL: 1.0817976e-3)", s37.delta_star);
    println!("Hk = {:.4} (XFOIL: 3.365)", s37.hk);
    println!("Rθ = {:.1} (XFOIL: 429.3)", s37.r_theta);
    println!("N = {:.4} (XFOIL: 8.234)", s37.ampl);
    
    // ===================================================================
    // XFOIL Station 38 (transition station)
    // This is where transition occurs
    // ===================================================================
    let x38 = 0.23851869;
    
    // IMPORTANT: Use the BLDIF U2 value from XFOIL debug (first Newton iteration),
    // NOT the final converged MRCHUE Ue value!
    // XFOIL BLDIF debug shows U2 = 1.3244895 for turbulent part
    let u38_bldif = 1.3244895;  // U2 during first Newton iteration in TRDIF
    let _u38_final = 1.2998873;  // Final converged Ue (for reference only)
    
    // What XFOIL produces after Newton converges:
    let xfoil_theta_38 = 0.37622254e-3;
    let xfoil_dstar_38 = 0.16022328e-2;
    let xfoil_ctau_38 = 0.067784877;
    let xfoil_hk_38 = 4.2587370;
    
    println!("\n=== XFOIL Station 38 EXPECTED (transition) ===");
    println!("x = {:.6}", x38);
    println!("theta = {:.6e}", xfoil_theta_38);
    println!("delta* = {:.6e}", xfoil_dstar_38);
    println!("Hk = {:.4}", xfoil_hk_38);
    println!("ctau = {:.6}", xfoil_ctau_38);
    
    // ===================================================================
    // XFOIL TRCHEK2 result at transition
    // ===================================================================
    // XT = 0.23752098 (transition point within interval)
    // WF1 = 0.064106, WF2 = 0.935894
    let xfoil_xt = 0.23752098;
    let xfoil_wf1 = 0.064106111;
    let xfoil_wf2 = 0.93589389;
    
    println!("\n=== XFOIL TRCHEK2 result ===");
    println!("XT = {:.6}", xfoil_xt);
    println!("WF1 = {:.6}, WF2 = {:.6}", xfoil_wf1, xfoil_wf2);
    
    // ===================================================================
    // Run RustFoil TRCHEK2
    // ===================================================================
    // Use exact XFOIL first-iteration values (before Newton convergence)
    let u2 = u38_bldif;  // BLDIF U2 value, not final converged Ue
    let t2_guess = s37.theta;  // Same as station 37 (initial guess)
    let d2_guess = s37.delta_star;  // Same as station 37 (initial guess)
    let rt2_guess = re * u2 * t2_guess;
    
    let trchek_simple = trchek2_stations(
        s37.x, x38,
        s37.hk, s37.theta, s37.r_theta, s37.delta_star, s37.u, s37.ampl,
        s37.hk, t2_guess, rt2_guess, d2_guess, u2,  // Use s37 values as initial for s38
        ncrit, msq, re,
    );
    
    println!("\n=== RustFoil trchek2_stations ===");
    println!("transition = {}", trchek_simple.transition);
    println!("XT = {:?}", trchek_simple.xt);
    println!("ampl2 = {:.4}", trchek_simple.ampl2);
    
    // ===================================================================
    // Run RustFoil trchek2_full
    // ===================================================================
    let trchek_full = trchek2_full(
        s37.x, x38,
        s37.theta, t2_guess,
        s37.delta_star, d2_guess,
        s37.u, u2,
        s37.hk, s37.hk,  // Hk guess
        s37.r_theta, re * u2 * t2_guess,
        s37.ampl,
        ncrit,
        None,
        msq,
        re,
    );
    
    println!("\n=== RustFoil trchek2_full ===");
    println!("transition = {}", trchek_full.transition);
    println!("XT = {:.6} (XFOIL: {:.6})", trchek_full.xt, xfoil_xt);
    println!("WF1 = {:.6} (XFOIL: {:.6})", trchek_full.wf1, xfoil_wf1);
    println!("WF2 = {:.6} (XFOIL: {:.6})", trchek_full.wf2, xfoil_wf2);
    
    println!("\n=== Transition point derivatives (used in chain rule) ===");
    println!("TT (theta at XT):");
    println!("  tt_t1 = {:.6} (expect ~WF1 = {:.6})", trchek_full.tt_t1, xfoil_wf1);
    println!("  tt_t2 = {:.6} (expect ~WF2 = {:.6})", trchek_full.tt_t2, xfoil_wf2);
    println!("  tt_a1 = {:.6e}", trchek_full.tt_a1);
    println!("DT (delta* at XT):");
    println!("  dt_d1 = {:.6}", trchek_full.dt_d1);
    println!("  dt_d2 = {:.6}", trchek_full.dt_d2);
    println!("XT derivatives (from implicit function theorem):");
    println!("  xt_t1 = {:.6e}", trchek_full.xt_t1);
    println!("  xt_t2 = {:.6e}", trchek_full.xt_t2);
    println!("  xt_a1 = {:.6e}", trchek_full.xt_a1);
    
    // The XT_T2 contribution to VS2[0][1] is:
    // jac_turb.vs1[0][4] * xt_t2
    // XFOIL VS1[0][4] = -0.326 (derivative w.r.t. X at transition point)
    // So contribution = -0.326 * 101 ≈ -33
    // This is consistent with VS2[0][1] going from 0.779 to something negative
    let vs1_0_4_approx = -0.326;  // XFOIL value
    let xt_t2_contribution = vs1_0_4_approx * trchek_full.xt_t2;
    println!("\nChain rule contribution to VS2[0][1]:");
    println!("  VS1[0][4] ≈ {:.3} (derivative w.r.t. X at XT)", vs1_0_4_approx);
    println!("  VS1[0][4] * xt_t2 = {:.3} * {:.2} = {:.2}", 
             vs1_0_4_approx, trchek_full.xt_t2, xt_t2_contribution);
    
    if trchek_full.transition {
        let xt_error = (trchek_full.xt - xfoil_xt).abs() / xfoil_xt * 100.0;
        println!("XT error = {:.2}%", xt_error);
        
        // ===================================================================
        // Run RustFoil TRDIF
        // ===================================================================
        // Create station 2 with exact XFOIL first-iteration values
        // CRITICAL: Use BLDIF U2 (1.3244895), NOT final Ue (1.2998873)
        let mut s38 = BlStation::new();
        s38.x = x38;
        s38.u = u38_bldif;  // BLDIF U2 from XFOIL debug
        s38.theta = s37.theta;  // Same as T1 (T2 = T1 in first iteration)
        s38.delta_star = s37.delta_star;  // Same as D1 (D2 = D1 in first iteration)
        s38.ctau = 0.03;  // S2 from XFOIL debug: initial turbulent ctau
        s38.ampl = trchek_full.ampl2;
        s38.is_laminar = false;
        s38.is_turbulent = true;
        
        blvar(&mut s38, FlowType::Turbulent, msq, re);
        
        println!("\n=== RustFoil TRDIF (single call, no iteration) ===");
        
        // Simple trdif
        let (res_simple, jac_simple) = trdif(&s37, &s38, trchek_full.xt, msq, re);
        println!("Simple TRDIF residuals:");
        println!("  third = {:.6e}", res_simple.res_third);
        println!("  mom   = {:.6e}", res_simple.res_mom);
        println!("  shape = {:.6e}", res_simple.res_shape);
        
        // Full trdif
        let (res_full, jac_full) = trdif_full(&s37, &s38, &trchek_full, msq, re);
        println!("\nFull TRDIF residuals:");
        println!("  third = {:.6e}", res_full.res_third);
        println!("  mom   = {:.6e}", res_full.res_mom);
        println!("  shape = {:.6e}", res_full.res_shape);
        
        // Compare Jacobians
        println!("\nJacobian VS1[0] (shear-lag row):");
        println!("  Simple: {:?}", jac_simple.vs1[0]);
        println!("  Full:   {:?}", jac_full.vs1[0]);
        
        println!("\nJacobian VS2 (derivatives w.r.t. station 2 - used in Newton):");
        println!("  Row 0 (shear-lag):");
        println!("    Simple: {:?}", jac_simple.vs2[0]);
        println!("    Full:   {:?}", jac_full.vs2[0]);
        println!("  Row 1 (momentum):");
        println!("    Simple: {:?}", jac_simple.vs2[1]);
        println!("    Full:   {:?}", jac_full.vs2[1]);
        println!("  Row 2 (shape):");
        println!("    Simple: {:?}", jac_simple.vs2[2]);
        println!("    Full:   {:?}", jac_full.vs2[2]);
        
        // XFOIL turbulent BLDIF VS2 (before TRDIF transformation):
        // [-0.1588, 0.7792, 0.3742, -0.00351, 0.3255]  (row 0)
        // [0.0, 3109.3, 0.6893, 4.051, -0.3362]        (row 1)
        // [-0.1165, 143.2, -57.39, -1.799, -5.777]     (row 2)
        //
        // XFOIL combined TRDIF residual: [-1.63e-3, 4.88e-2, 3.62e-3]
        // (row 0 from turb, rows 1-2 are sums)
        
        println!("\n=== XFOIL BLDIF (turbulent part, raw): ===");
        println!("VS1 Row 0: [0.114, 0.779, 0.374, 0.00352, -0.326]");
        println!("VS2 Row 0: [-0.159, 0.779, 0.374, -0.00351, 0.326]");
        println!("VS1 Row 1: [0.0, -3112, 0.686, -4.05, 0.336]");
        println!("VS2 Row 1: [0.0, 3109.3, 0.689, 4.05, -0.336]");
        println!("VS1 Row 2: [-0.154, -110, 53.7, 1.80, 5.78]");
        println!("VS2 Row 2: [-0.117, 143.2, -57.4, -1.80, -5.78]");
        
        // The chain rule transformation for bt2[0][1] is:
        // bt2[0][1] = vs2[0][1] + vs1[0][0]*st_t2 + vs1[0][1]*tt_t2 + vs1[0][2]*dt_t2 + vs1[0][3]*ut_t2 + vs1[0][4]*xt_t2
        // Expected: ~0.779 + 0.114*st_t2 + 0.779*tt_t2 + 0.374*dt_t2 + ...
        // If tt_t2 ≈ 0.936, this should add ~0.73, giving ~1.5, NOT -63.9
        println!("\n=== Chain rule analysis ===");
        println!("The large difference in VS2[0][1] suggests tt_t2 or other derivatives");
        println!("are much larger than simple interpolation weights (WF1, WF2).");
        println!("This happens because XT itself depends on station 2 variables.");
    }
    
    // ===================================================================
    // Key insight: The theta error likely comes from:
    // 1. Different initial conditions for Newton iteration
    // 2. Different handling of the TRDIF residuals/Jacobian
    // 3. The fact that XFOIL uses interpolated values at XT, not station 2 values
    // ===================================================================
    println!("\n=== Analysis ===");
    println!("XFOIL theta growth: {:.4e} -> {:.4e} ({:.1}% increase)", 
             s37.theta, xfoil_theta_38, 
             (xfoil_theta_38 - s37.theta) / s37.theta * 100.0);
    println!("XFOIL Hk change: {:.3} -> {:.3}", s37.hk, xfoil_hk_38);
    println!("\nKey observation:");
    println!("XFOIL uses T1=T2={:.6e} for TRDIF (interpolated at XT)", s37.theta);
    println!("This is because WF1 ≈ 0.064, so TT ≈ T1 at the transition point.");
}

/// Test to understand what XFOIL's BLDIF calls look like during TRDIF
#[test]
fn test_xfoil_trdif_bldif_calls() {
    println!("\n=== XFOIL TRDIF BLDIF Calls Analysis ===\n");
    
    // From xfoil_debug.json, BLDIF calls during TRDIF at IBL 38:
    
    // LAMINAR BLDIF (call_id 3139):
    // X1 = 0.22295522, X2 = 0.23752098 (XT)
    // T1 = T2 = 0.32148343e-3
    // D1 = D2 = 0.10817976e-2
    // S1 = 0.03, S2 = 0.0 (laminar at transition point)
    // A1 = 8.234, A2 = 9.0 (Ncrit)
    // VSREZ = [2.14e-5, 4.57e-2, -9.04e-4, 0]
    
    println!("LAMINAR part (X1 -> XT):");
    println!("  X1 = 0.22295522, X2 = 0.23752098 (XT)");
    println!("  T1 = T2 = 3.2148e-4 (SAME - interpolated)");
    println!("  D1 = D2 = 1.0818e-3 (SAME - interpolated)");
    println!("  VSREZ = [2.14e-5, 4.57e-2, -9.04e-4, 0]");
    
    // TURBULENT BLDIF (call_id 3145):
    // X1 = 0.23752098 (XT), X2 = 0.23851869
    // T1 = T2 = 0.32148343e-3 (still same!)
    // D1 = D2 = 0.10817976e-2
    // S1 = 0.03968823 (computed ST at transition)
    // S2 = 0.03
    // VSREZ = [-1.63e-3, 3.14e-3, 4.52e-3, 0]
    
    println!("\nTURBULENT part (XT -> X2):");
    println!("  X1 = 0.23752098 (XT), X2 = 0.23851869");
    println!("  T1 = T2 = 3.2148e-4 (SAME - not Newton-iterated yet)");
    println!("  D1 = D2 = 1.0818e-3");
    println!("  S1 = 0.03969 (computed ctau at transition)");
    println!("  VSREZ = [-1.63e-3, 3.14e-3, 4.52e-3, 0]");
    
    println!("\nCOMBINED residuals (row 0 from turb, rows 1-2 summed):");
    println!("  Row 0 (shear-lag): -1.63e-3");
    println!("  Row 1 (momentum):  4.57e-2 + 3.14e-3 = 4.88e-2");
    println!("  Row 2 (shape):     -9.04e-4 + 4.52e-3 = 3.62e-3");
    
    println!("\n=== KEY INSIGHT ===");
    println!("XFOIL's TRDIF uses the SAME theta/delta_star for both");
    println!("laminar and turbulent BLDIF calls (the interpolated values at XT).");
    println!("The Newton iteration THEN converges these to the final values.");
    println!("This is fundamentally different from re-solving with different theta.");
}
