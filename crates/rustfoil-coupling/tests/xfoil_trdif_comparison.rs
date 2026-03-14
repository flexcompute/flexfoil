//! Direct comparison of RustFoil TRDIF with XFOIL instrumented output
//!
//! Uses exact values from xfoil_debug.json (Jan 25 2026 run) at transition station IBL 38

use rustfoil_bl::closures::trchek2_full;
use rustfoil_bl::equations::{blvar, trdif_full, FlowType};
use rustfoil_bl::state::BlStation;

/// Compare RustFoil TRDIF output with XFOIL's transformed Jacobian
///
/// XFOIL conditions from xfoil_debug.json at IBL 38:
/// - x1 = 0.22560078, x2 = 0.24114425
/// - T1 = T2 = 0.32419833e-3
/// - U1 = 1.3334211, U2 = 1.3227175  
/// - ampl1 = 8.3759789
/// - Ncrit = 9.0
/// - xt = 0.23734327
/// - WF1 = 0.24453843, WF2 = 0.75546157
#[test]
fn test_trdif_jacobian_vs_xfoil() {
    let re = 1e6;
    let msq = 0.0;
    let ncrit = 9.0;
    
    // Station 37 (previous laminar station) - exact XFOIL values
    let x1 = 0.22560078;
    let t1 = 0.32419833e-3;
    let d1 = t1 * 3.3825711;  // H = Hk ≈ 3.38 at low Mach
    let u1 = 1.3334211;
    let ampl1 = 8.3759789;
    
    // Station 38 (transition station) - exact XFOIL values
    let x2 = 0.24114425;
    let t2 = 0.32419833e-3;  // Same as T1 (initial Newton iteration)
    let d2 = t2 * 3.3825711;  // Same H as station 1
    let u2 = 1.3227175;  // XFOIL TRCHEK2 U2 value (NOT final Ue)
    
    // Expected XFOIL values
    let xfoil_xt = 0.23734327;
    let xfoil_wf1 = 0.24453843;
    let xfoil_wf2 = 0.75546157;
    
    // Create station 1
    let mut s1 = BlStation::new();
    s1.x = x1;
    s1.u = u1;
    s1.theta = t1;
    s1.delta_star = d1;
    s1.ctau = 0.03;
    s1.ampl = ampl1;
    s1.is_laminar = true;
    blvar(&mut s1, FlowType::Laminar, msq, re);
    
    println!("\n=== Station 1 (IBL 37) ===");
    println!("x = {:.6}, u = {:.6}", s1.x, s1.u);
    println!("theta = {:.6e}, delta* = {:.6e}", s1.theta, s1.delta_star);
    println!("Hk = {:.4} (XFOIL: 3.3826)", s1.hk);
    println!("Rθ = {:.1} (XFOIL: 432.3)", s1.r_theta);
    println!("ampl = {:.4} (XFOIL: 8.376)", s1.ampl);
    
    // Create station 2 (initial guess, same as station 1 for theta/delta*)
    let mut s2 = BlStation::new();
    s2.x = x2;
    s2.u = u2;
    s2.theta = t2;
    s2.delta_star = d2;
    s2.ctau = 0.03;  // Initial turbulent ctau
    s2.ampl = 9.2;   // Near Ncrit
    s2.is_laminar = false;
    s2.is_turbulent = true;
    blvar(&mut s2, FlowType::Turbulent, msq, re);
    
    println!("\n=== Station 2 (IBL 38) - Initial ===");
    println!("x = {:.6}, u = {:.6}", s2.x, s2.u);
    println!("theta = {:.6e}, delta* = {:.6e}", s2.theta, s2.delta_star);
    println!("Hk = {:.4}", s2.hk);
    println!("Rθ = {:.1}", s2.r_theta);
    
    // Run trchek2_full
    let tr = trchek2_full(
        s1.x, s2.x,
        s1.theta, s2.theta,
        s1.delta_star, s2.delta_star,
        s1.u, s2.u,
        s1.hk, s2.hk,
        s1.r_theta, s2.r_theta,
        s1.ampl, ncrit, None, msq, re,
    );
    
    println!("\n=== TRCHEK2 Result ===");
    println!("transition = {}", tr.transition);
    println!("XT = {:.6} (XFOIL: {:.6})", tr.xt, xfoil_xt);
    println!("WF1 = {:.6} (XFOIL: {:.6})", tr.wf1, xfoil_wf1);
    println!("WF2 = {:.6} (XFOIL: {:.6})", tr.wf2, xfoil_wf2);
    
    if tr.transition {
        let xt_err = (tr.xt - xfoil_xt).abs() / xfoil_xt * 100.0;
        let wf1_err = (tr.wf1 - xfoil_wf1).abs() / xfoil_wf1.max(0.01) * 100.0;
        println!("XT error = {:.2}%, WF1 error = {:.2}%", xt_err, wf1_err);
        
        // Run TRDIF
        let (res, jac) = trdif_full(&s1, &s2, &tr, msq, re);
        
        // XFOIL IBL 37 reference values (call 3148 in xfoil_debug.json)
        // These are for the SAME station conditions as our test!
        let xfoil_vsrez = [-0.2624109e-2, 0.4825283e-1, 0.1663855e-1];
        let xfoil_vs2_0 = [-0.1703540, -0.4596360e2, 0.1340944e2, -0.8422569e-3, 0.3276345];
        let xfoil_vs2_1 = [0.0, 0.3065808e4, 0.1001512e2, 0.4071182e1, -0.3133783];
        let xfoil_vs2_2 = [-0.4438060, 0.7283655e3, -0.1713198e3, -0.1798596e1, -0.5803800e1];
        
        println!("\n=== TRDIF Comparison (vs XFOIL IBL 37, call 3148) ===");
        println!("VSREZ (residuals):");
        println!("  RustFoil: [{:.6e}, {:.6e}, {:.6e}]", res.res_third, res.res_mom, res.res_shape);
        println!("  XFOIL:    [{:.6e}, {:.6e}, {:.6e}]", xfoil_vsrez[0], xfoil_vsrez[1], xfoil_vsrez[2]);
        
        // Check residual match
        let res_err_0 = (res.res_third - xfoil_vsrez[0]).abs();
        let res_err_1 = (res.res_mom - xfoil_vsrez[1]).abs();
        let res_err_2 = (res.res_shape - xfoil_vsrez[2]).abs();
        println!("  Errors:   [{:.2e}, {:.2e}, {:.2e}]", res_err_0, res_err_1, res_err_2);
        
        println!("\nVS2[0] (shear-lag row):");
        println!("  RustFoil: [{:.5}, {:.4}, {:.4}, {:.6e}, {:.5}]", 
                 jac.vs2[0][0], jac.vs2[0][1], jac.vs2[0][2], jac.vs2[0][3], jac.vs2[0][4]);
        println!("  XFOIL:    [{:.5}, {:.4}, {:.4}, {:.6e}, {:.5}]", 
                 xfoil_vs2_0[0], xfoil_vs2_0[1], xfoil_vs2_0[2], xfoil_vs2_0[3], xfoil_vs2_0[4]);
        
        println!("\nVS2[1] (momentum row):");
        println!("  RustFoil: [{:.3}, {:.2}, {:.4}, {:.4}, {:.5}]", 
                 jac.vs2[1][0], jac.vs2[1][1], jac.vs2[1][2], jac.vs2[1][3], jac.vs2[1][4]);
        println!("  XFOIL:    [{:.3}, {:.2}, {:.4}, {:.4}, {:.5}]", 
                 xfoil_vs2_1[0], xfoil_vs2_1[1], xfoil_vs2_1[2], xfoil_vs2_1[3], xfoil_vs2_1[4]);
        
        println!("\nVS2[2] (shape row):");
        println!("  RustFoil: [{:.5}, {:.3}, {:.3}, {:.5}, {:.5}]", 
                 jac.vs2[2][0], jac.vs2[2][1], jac.vs2[2][2], jac.vs2[2][3], jac.vs2[2][4]);
        println!("  XFOIL:    [{:.5}, {:.3}, {:.3}, {:.5}, {:.5}]", 
                 xfoil_vs2_2[0], xfoil_vs2_2[1], xfoil_vs2_2[2], xfoil_vs2_2[3], xfoil_vs2_2[4]);
        
        // Verify VS2 matches
        let vs2_max_err = [
            (jac.vs2[0][1] - xfoil_vs2_0[1]).abs() / xfoil_vs2_0[1].abs(),
            (jac.vs2[1][1] - xfoil_vs2_1[1]).abs() / xfoil_vs2_1[1].abs(),
            (jac.vs2[2][1] - xfoil_vs2_2[1]).abs() / xfoil_vs2_2[1].abs(),
        ];
        println!("\nVS2[:,1] (theta column) relative errors: [{:.4}%, {:.4}%, {:.4}%]",
                 vs2_max_err[0] * 100.0, vs2_max_err[1] * 100.0, vs2_max_err[2] * 100.0);
        
        // Assert residuals match (tolerance 1e-7 for floating point differences)
        assert!((res.res_third - xfoil_vsrez[0]).abs() < 1e-7, "Shear-lag residual mismatch");
        assert!((res.res_mom - xfoil_vsrez[1]).abs() < 1e-7, "Momentum residual mismatch");
        assert!((res.res_shape - xfoil_vsrez[2]).abs() < 1e-7, "Shape residual mismatch");
        
        // Assert VS2 theta column matches (key for Newton convergence)
        assert!(vs2_max_err[0] < 0.001, "VS2[0][1] mismatch");
        assert!(vs2_max_err[1] < 0.001, "VS2[1][1] mismatch");
        assert!(vs2_max_err[2] < 0.001, "VS2[2][1] mismatch");
        
        println!("\n✓ TRDIF Jacobian matches XFOIL within tolerance!");
    }
}
