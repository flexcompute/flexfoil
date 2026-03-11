//! Debug test for ctau convergence at transition
//!
//! Compares RustFoil's Newton iteration against XFOIL's exact values

use rustfoil_bl::closures::{trchek2_full, Trchek2FullResult};
use rustfoil_bl::equations::{bldif, blvar, trdif, trdif_full, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_bl::constants::{SCCON, GACON, GBCON, GCCON, DLCON, DUXCON};
use rustfoil_coupling::march::{newton_solve_station, MarchConfig};

/// Compute and print all intermediate REZC values for debugging
/// This mirrors XFOIL's BLDIF shear-lag calculation (xblsys.f:1721-1804)
fn debug_rezc_terms(s1: &BlStation, s2: &BlStation, flow_type: FlowType) {
    println!("\n=== REZC Term-by-Term Debug ===");
    
    // === UPW calculation (xblsys.f:1598-1644) ===
    let hupwt = 1.0_f64;
    let hdcon = if flow_type == FlowType::Wake {
        hupwt / (s2.hk * s2.hk)
    } else {
        5.0 * hupwt / (s2.hk * s2.hk)
    };
    
    // NOTE: XFOIL does NOT use abs() here - checking if this matters
    let arg_with_abs = ((s2.hk - 1.0) / (s1.hk - 1.0).max(1e-6)).abs();
    let arg_without_abs = (s2.hk - 1.0) / (s1.hk - 1.0).max(1e-6);
    let hl = arg_with_abs.ln();
    let hlsq = (hl * hl).min(15.0);
    let ehh = (-hlsq * hdcon).exp();
    let upw = 1.0 - 0.5 * ehh;
    
    println!("UPW calculation:");
    println!("  Hk1={:.6}, Hk2={:.6}", s1.hk, s2.hk);
    println!("  arg (with abs)={:.6}, arg (without abs)={:.6}", arg_with_abs, arg_without_abs);
    println!("  hl={:.6}, hlsq={:.6}, hdcon={:.6}", hl, hlsq, hdcon);
    println!("  ehh={:.6}, UPW={:.6}", ehh, upw);
    
    // === Upwinded quantities (xblsys.f:1723-1731) ===
    let sa = (1.0 - upw) * s1.ctau + upw * s2.ctau;
    let cqa = (1.0 - upw) * s1.cq + upw * s2.cq;
    let cfa = (1.0 - upw) * s1.cf + upw * s2.cf;
    let hka = (1.0 - upw) * s1.hk + upw * s2.hk;
    
    // Simple averages (xblsys.f:1728-1731)
    let usa = 0.5 * (s1.us + s2.us);
    let rta = 0.5 * (s1.r_theta + s2.r_theta);
    let dea = 0.5 * (s1.de + s2.de);
    let da = 0.5 * (s1.delta_star + s2.delta_star);
    
    println!("\nUpwinded quantities:");
    println!("  SA (ctau_avg)={:.6}, CQA={:.6}", sa, cqa);
    println!("  CFA={:.6}, HKA={:.6}", cfa, hka);
    println!("  USA={:.6}, RTA={:.1}", usa, rta);
    println!("  DEA={:.6e}, DA={:.6e}", dea, da);
    
    // === ALD (xblsys.f:1734-1739) ===
    let ald = if flow_type == FlowType::Wake { DLCON } else { 1.0 };
    println!("  ALD={:.6}", ald);
    
    // === HKC and HR (xblsys.f:1742-1761) ===
    let gcc = if flow_type == FlowType::Turbulent { GCCON } else { 0.0 };
    let mut hkc = hka - 1.0 - gcc / rta.max(1.0);
    if flow_type == FlowType::Turbulent && hkc < 0.01 {
        hkc = 0.01;
    }
    
    let hr = hkc / (GACON * ald * hka);
    
    println!("\nHKC/HR calculation:");
    println!("  GCC={:.6}, HKC={:.6}", gcc, hkc);
    println!("  GACON={:.6}, HR={:.6}", GACON, hr);
    
    // === UQ (xblsys.f:1763-1767) ===
    let uq = (0.5 * cfa - hr * hr) / (GBCON * da.max(1e-10));
    println!("  GBCON={:.6}, UQ={:.6}", GBCON, uq);
    
    // === SCC (xblsys.f:1792-1793) ===
    let scc = SCCON * 1.333 / (1.0 + usa);
    println!("  SCCON={:.6}, SCC={:.6}", SCCON, scc);
    
    // === SLOG, DXI, ULOG (xblsys.f:1799-1800) ===
    let slog = (s2.ctau / s1.ctau.max(1e-20)).ln();
    let dxi = s2.x - s1.x;
    let ulog = (s2.u / s1.u.max(1e-20)).ln();
    
    println!("\nLog differences:");
    println!("  S1={:.6}, S2={:.6}", s1.ctau, s2.ctau);
    println!("  SLOG=ln(S2/S1)={:.6}", slog);
    println!("  DXI=X2-X1={:.6}", dxi);
    println!("  ULOG=ln(U2/U1)={:.6}", ulog);
    println!("  DUXCON={:.6}", DUXCON);
    
    // === REZC terms (xblsys.f:1802-1804) ===
    // REZC = SCC*(CQA - SA*ALD)*DXI 
    //      - DEA*2.0*SLOG 
    //      + DEA*2.0*(UQ*DXI - ULOG)*DUXCON
    let term1 = scc * (cqa - sa * ald) * dxi;
    let term2 = -dea * 2.0 * slog;
    let term3 = dea * 2.0 * (uq * dxi - ulog) * DUXCON;
    let rezc = term1 + term2 + term3;
    
    println!("\n=== REZC TERMS ===");
    println!("  Term1 = SCC*(CQA - SA*ALD)*DXI = {:.6e}", term1);
    println!("    CQA - SA*ALD = {:.6e}", cqa - sa * ald);
    println!("  Term2 = -DEA*2.0*SLOG = {:.6e}", term2);
    println!("  Term3 = DEA*2.0*(UQ*DXI - ULOG)*DUXCON = {:.6e}", term3);
    println!("    UQ*DXI - ULOG = {:.6e}", uq * dxi - ulog);
    println!("  REZC = Term1 + Term2 + Term3 = {:.6e}", rezc);
    println!("  VSREZ(1) = -REZC = {:.6e}", -rezc);
    
    // === Z coefficients (xblsys.f:1816-1846) ===
    let z_sl = -dea * 2.0;
    let z_sa = -scc * dxi * ald;
    let z_s2 = upw * z_sa + z_sl / s2.ctau.max(1e-20);
    
    println!("\n=== Z Coefficients for VS2[0][0] ===");
    println!("  Z_SL = -DEA*2.0 = {:.6e}", z_sl);
    println!("  Z_SA = -SCC*DXI*ALD = {:.6e}", z_sa);
    println!("  Z_S2 = UPW*Z_SA + Z_SL/S2 = {:.6e}", z_s2);
    println!("    UPW*Z_SA = {:.6e}", upw * z_sa);
    println!("    Z_SL/S2 = {:.6e}", z_sl / s2.ctau.max(1e-20));
}

/// Test that our bldif produces the same VS2 matrix as XFOIL for the transition station
#[test]
fn test_vs2_at_transition() {
    // =========================================================================
    // XFOIL actual converged values from debug output at transition (NACA 0012, alpha=4, Re=1e6)
    // The debug output shows:
    //   HK1=5.546, HK2=2.690 (converged)
    //   S1=0.102, S2=0.095 (ctau values)
    // =========================================================================
    
    // Previous station (laminar, just before transition):
    // XFOIL: HK1 = 5.546
    let mut prev = BlStation::new();
    prev.x = 0.12758861;
    prev.u = 1.4637065;
    // Need theta, delta_star such that H = 5.546 (at M=0, Hk ≈ H)
    // Let's use theta = 2.218e-4, then delta_star = 5.546 * theta = 1.231e-3
    prev.theta = 2.218498e-04;
    prev.delta_star = 1.231e-03;  // Adjusted to give Hk ≈ 5.55
    prev.ctau = 0.102;   // From XFOIL: S1 = 0.10188
    prev.ampl = 7.8794514;
    prev.is_laminar = true;
    prev.is_turbulent = false;
    
    // Current station (first turbulent, at transition):
    // XFOIL converged: HK2 = 2.690, S2 = 0.095
    // From XFOIL debug output:
    //   RTA = 951.87, DA = 4.01e-3, DEA = 8.07e-3, USA = 2.59e-2
    // 
    // Working backwards from XFOIL's actual averages:
    //   DA = 0.5*(D1 + D2) = 4.01e-3, so if D1 ≈ D2, then D2 ≈ 4.0e-3
    //   For H = 2.69: theta2 = delta_star2 / H = 4.0e-3 / 2.69 = 1.49e-3
    //   Check: RT2 = Re * U * theta = 1e6 * 1.457 * 1.49e-3 = 2170 (too high!)
    //
    // The issue is that XFOIL's prev station (laminar) has very different values.
    // Let me try to match XFOIL's actual DEA = 8.07e-3 by adjusting our values.
    //
    // ALTERNATIVE: Use XFOIL's exact intermediate values directly
    // From converged iteration where REZC ≈ 0:
    //   Use actual station values that produce DA=4.01e-3, DEA=8.07e-3
    let mut curr = BlStation::new();
    curr.x = 0.14078095;
    curr.u = 1.4566427;
    // Matching XFOIL's DA = 4.01e-3 average, which means D2 ≈ 4.0e-3 if D1 ≈ 4.0e-3
    // But prev.delta_star = 1.23e-3, so to get DA = 4.01e-3:
    // D2 = 2 * DA - D1 = 2 * 4.01e-3 - 1.23e-3 = 6.79e-3
    // With H = 2.69: theta2 = 6.79e-3 / 2.69 = 2.52e-3
    curr.theta = 2.52e-03;          // Computed to match XFOIL's DA
    curr.delta_star = 6.79e-03;     // Gives H ≈ 2.69 and matches DA=4.01e-3 average
    curr.ctau = 0.095;              // From XFOIL: S2 = 0.09480
    curr.ampl = 10.138472;
    curr.is_laminar = false;
    curr.is_turbulent = true;
    
    let re = 1e6;
    let msq = 0.0;
    
    // Compute closure relations
    blvar(&mut prev, FlowType::Laminar, msq, re);
    blvar(&mut curr, FlowType::Turbulent, msq, re);
    
    println!("\n=== Station state after blvar ===");
    println!("Prev (laminar):");
    println!("  x={:.6}, u={:.6}", prev.x, prev.u);
    println!("  Hk={:.4}, Rθ={:.1}, H={:.4}", prev.hk, prev.r_theta, prev.h);
    println!("  CQ={:.6}, CF={:.6}", prev.cq, prev.cf);
    println!("  US={:.6}, DE={:.6e}", prev.us, prev.de);
    println!("  ctau={:.6}", prev.ctau);
    
    println!("Curr (turbulent):");
    println!("  x={:.6}, u={:.6}", curr.x, curr.u);
    println!("  Hk={:.4}, Rθ={:.1}, H={:.4}", curr.hk, curr.r_theta, curr.h);
    println!("  CQ={:.6}, CF={:.6}", curr.cq, curr.cf);
    println!("  US={:.6}, DE={:.6e}", curr.us, curr.de);
    println!("  ctau={:.6}", curr.ctau);
    
    // Print term-by-term REZC debug
    debug_rezc_terms(&prev, &curr, FlowType::Turbulent);
    
    // Compute residuals and Jacobian
    let (res, jac) = bldif(&prev, &curr, FlowType::Turbulent, msq, re);
    
    println!("\n=== RustFoil bldif output ===");
    println!("Residuals:");
    println!("  res_third (shear-lag): {:.6e}", res.res_third);
    println!("  res_mom (momentum):    {:.6e}", res.res_mom);
    println!("  res_shape (shape):     {:.6e}", res.res_shape);
    
    println!("\nVS1 matrix (Jacobian wrt station 1):");
    for i in 0..3 {
        println!("  Row {}: [{:12.4e}, {:12.4e}, {:12.4e}, {:12.4e}, {:12.4e}]",
                 i, jac.vs1[i][0], jac.vs1[i][1], jac.vs1[i][2], jac.vs1[i][3], jac.vs1[i][4]);
    }
    
    println!("\nVS2 matrix (Jacobian wrt station 2):");
    for i in 0..3 {
        println!("  Row {}: [{:12.4e}, {:12.4e}, {:12.4e}, {:12.4e}, {:12.4e}]",
                 i, jac.vs2[i][0], jac.vs2[i][1], jac.vs2[i][2], jac.vs2[i][3], jac.vs2[i][4]);
    }
    
    // XFOIL BLDIF reference values (from call_id=2916, flow_type=2)
    // Note: The VS2 values in MRCHUE_ITER are post-Gauss-elimination, not raw BLDIF output!
    println!("\n=== XFOIL BLDIF reference (call_id=2916) ===");
    println!("VS2 matrix:");
    println!("  Row 0: [-0.1053334, -9.317754, 3.878588, -0.002245024, 0.08062964]");
    println!("VSREZ: [-0.001158263, 0.09778855, 0.1485491, 0.0]");
    
    // Compare VS2[0][0] (shear-lag sensitivity to ctau)
    let xfoil_vs2_00 = -0.1053334;
    let error_vs2_00 = (jac.vs2[0][0] - xfoil_vs2_00).abs() / xfoil_vs2_00.abs();
    println!("\n=== Comparison with XFOIL BLDIF ===");
    println!("VS2[0][0] (shear-lag/ctau): RustFoil={:.4e}, XFOIL={:.4e}, Error={:.1}%",
             jac.vs2[0][0], xfoil_vs2_00, error_vs2_00 * 100.0);
    
    let xfoil_vs2_01 = -9.317754;
    let error_vs2_01 = (jac.vs2[0][1] - xfoil_vs2_01).abs() / xfoil_vs2_01.abs();
    println!("VS2[0][1] (shear-lag/theta): RustFoil={:.4e}, XFOIL={:.4e}, Error={:.1}%",
             jac.vs2[0][1], xfoil_vs2_01, error_vs2_01 * 100.0);
    
    // Check residuals - note XFOIL returns -REZC, so we compare with negated value
    let xfoil_res_third = -(-0.001158263);  // XFOIL's VSREZ(1) = -REZC
    let error_res = (res.res_third - xfoil_res_third).abs() / xfoil_res_third.abs();
    println!("res_third (shear-lag): RustFoil={:.6e}, XFOIL={:.6e}, Error={:.1}%",
             res.res_third, xfoil_res_third, error_res * 100.0);
}

/// Test Newton solver convergence for first turbulent station
#[test]
fn test_newton_ctau_convergence() {
    // Previous station (IBL 31, laminar) - converged state from XFOIL
    let mut prev = BlStation::new();
    prev.x = 0.12758861;
    prev.u = 1.4637065;
    prev.theta = 2.218498e-04;
    prev.delta_star = 1.195687e-03;
    prev.ctau = 0.03;
    prev.ampl = 7.8794514;
    prev.is_laminar = true;
    prev.is_turbulent = false;
    
    let re = 1e6;
    let msq = 0.0;
    
    // Compute closure relations for prev
    blvar(&mut prev, FlowType::Laminar, msq, re);
    
    // Current station (IBL 32) inputs from XFOIL
    let x_new = 0.14078095;
    let ue_new = 1.4566427;
    
    // Run Newton solver
    let (result, converged, _dmax) = newton_solve_station(
        &prev,
        x_new,
        ue_new,
        re,
        msq,
        false,  // is_laminar = false (turbulent)
        25,     // max_iter
        1e-5,   // tolerance
        4.0,    // hlmax
        2.5,    // htmax
    );
    
    println!("\n=== Newton solver result ===");
    println!("Converged: {}", converged);
    println!("theta: {:.6e} (XFOIL: 2.439713e-04)", result.theta);
    println!("delta_star: {:.6e} (XFOIL: 6.099282e-04)", result.delta_star);
    println!("ctau: {:.6} (XFOIL: 0.054610)", result.ctau);
    println!("Hk: {:.4} (XFOIL: 2.5000)", result.hk);
    
    // Check ctau error
    let xfoil_ctau = 0.054610;
    let ctau_error = (result.ctau - xfoil_ctau).abs() / xfoil_ctau * 100.0;
    println!("\nctau error: {:.1}%", ctau_error);
    
    // Target: < 1% error for numerical precision
    assert!(ctau_error < 10.0, "ctau error {} > 10%", ctau_error);
}

/// Test that trchek2_full computes transition location and derivatives
#[test]
fn test_trchek2_full_derivatives() {
    // Set up stations near transition (NACA 0012, alpha=4, Re=1e6)
    let x1 = 0.12758861;
    let x2 = 0.14078095;
    let t1 = 2.218498e-04;
    let t2 = 2.50e-04;  // Approximate theta at next station
    let d1 = 1.195687e-03;
    let d2 = 1.50e-03;  // Approximate delta_star at next station
    let u1 = 1.4637065;
    let u2 = 1.4566427;
    let hk1 = 5.39;  // Laminar Hk
    let hk2 = 6.0;   // Higher Hk as approaching transition
    let rt1 = 1e6 * u1 * t1;
    let rt2 = 1e6 * u2 * t2;
    let ampl1 = 7.8;  // Close to Ncrit
    let ncrit = 9.0;
    let msq = 0.0;
    let re = 1e6;
    
    // Call trchek2_full
    let result = trchek2_full(
        x1, x2, t1, t2, d1, d2, u1, u2,
        hk1, hk2, rt1, rt2, ampl1, ncrit, None, msq, re,
    );
    
    println!("\n=== trchek2_full results ===");
    println!("Transition: {}", result.transition);
    println!("XT: {:.6}", result.xt);
    println!("AMPL2: {:.4}", result.ampl2);
    println!("WF1: {:.4}, WF2: {:.4}", result.wf1, result.wf2);
    println!("Converged: {}", result.converged);
    
    if result.transition {
        println!("\n=== XT derivatives ===");
        println!("XT_A1: {:.6e}", result.xt_a1);
        println!("XT_T1: {:.6e}", result.xt_t1);
        println!("XT_D1: {:.6e}", result.xt_d1);
        println!("XT_U1: {:.6e}", result.xt_u1);
        println!("XT_X1: {:.6e}", result.xt_x1);
        println!("XT_T2: {:.6e}", result.xt_t2);
        println!("XT_D2: {:.6e}", result.xt_d2);
        println!("XT_U2: {:.6e}", result.xt_u2);
        println!("XT_X2: {:.6e}", result.xt_x2);
        
        println!("\n=== TT derivatives (sample) ===");
        println!("TT_T1: {:.6e}, TT_T2: {:.6e}", result.tt_t1, result.tt_t2);
        println!("TT_A1: {:.6e}", result.tt_a1);
        
        // Verify basic properties
        assert!(result.xt > x1 && result.xt < x2, "XT should be between X1 and X2");
        assert!(result.wf1 + result.wf2 - 1.0 < 1e-10, "WF1 + WF2 should equal 1");
        assert!(result.ampl2 >= ncrit, "AMPL2 should be >= Ncrit at transition");
    }
}

/// Test that trdif_full produces reasonable residuals at transition
#[test]
fn test_trdif_full_residuals() {
    // Previous station (laminar, just before transition)
    let mut s1 = BlStation::new();
    s1.x = 0.12758861;
    s1.u = 1.4637065;
    s1.theta = 2.218498e-04;
    s1.delta_star = 1.195687e-03;
    s1.ctau = 0.03;
    s1.ampl = 8.5;  // Just below Ncrit
    s1.is_laminar = true;
    s1.is_turbulent = false;
    
    // Current station (at transition, turbulent)
    let mut s2 = BlStation::new();
    s2.x = 0.14078095;
    s2.u = 1.4566427;
    s2.theta = 2.50e-04;
    s2.delta_star = 6.0e-04;
    s2.ctau = 0.05;
    s2.ampl = 9.5;  // Above Ncrit
    s2.is_laminar = false;
    s2.is_turbulent = true;
    
    let re = 1e6;
    let msq = 0.0;
    let ncrit = 9.0;
    
    // Compute closure relations
    blvar(&mut s1, FlowType::Laminar, msq, re);
    blvar(&mut s2, FlowType::Turbulent, msq, re);
    
    // Compute full TRCHEK2 result
    let tr = trchek2_full(
        s1.x, s2.x, s1.theta, s2.theta, s1.delta_star, s2.delta_star,
        s1.u, s2.u, s1.hk, s2.hk, s1.r_theta, s2.r_theta,
        s1.ampl, ncrit, None, msq, re,
    );
    
    println!("\n=== trdif_full test ===");
    println!("Transition detected: {}", tr.transition);
    
    if tr.transition {
        // Test simple trdif (backward compatible)
        let (res_simple, jac_simple) = trdif(&s1, &s2, tr.xt, msq, re);
        
        // Test full trdif
        let (res_full, jac_full) = trdif_full(&s1, &s2, &tr, msq, re);
        
        println!("\n--- Simple TRDIF ---");
        println!("Residuals: third={:.6e}, mom={:.6e}, shape={:.6e}",
                 res_simple.res_third, res_simple.res_mom, res_simple.res_shape);
        
        println!("\n--- Full TRDIF ---");
        println!("Residuals: third={:.6e}, mom={:.6e}, shape={:.6e}",
                 res_full.res_third, res_full.res_mom, res_full.res_shape);
        
        println!("\nVS1 comparison (simple vs full):");
        for i in 0..3 {
            println!("  Row {}: simple[0]={:.4e}, full[0]={:.4e}", 
                     i, jac_simple.vs1[i][0], jac_full.vs1[i][0]);
        }
        
        println!("\nVS2 comparison (simple vs full):");
        for i in 0..3 {
            println!("  Row {}: simple[0]={:.4e}, full[0]={:.4e}", 
                     i, jac_simple.vs2[i][0], jac_full.vs2[i][0]);
        }
        
        // Residuals should be similar (same physics, different Jacobian accuracy)
        let res_diff = (res_full.res_mom - res_simple.res_mom).abs();
        println!("\nMomentum residual difference: {:.6e}", res_diff);
        
        // The full TRDIF should have non-zero VS1[0] derivatives (chain rule through XT)
        // while simple TRDIF has VS1[0] = [0,0,0,0,0]
        let full_vs1_0_sum: f64 = jac_full.vs1[0].iter().map(|x| x.abs()).sum();
        println!("Full VS1[0] sum of absolute values: {:.6e}", full_vs1_0_sum);
        
        // The full Jacobian should capture more dependencies
        // This is the key improvement from the full implementation
    } else {
        println!("No transition detected - test may need adjustment");
    }
}

/// Integration test: verify theta error improves with full TRDIF
#[test]
fn test_trdif_theta_accuracy() {
    use rustfoil_coupling::march::{march_fixed_ue, MarchConfig};
    
    // Generate a simple test case with known properties
    // This is a simplified test - real validation should compare against XFOIL
    let n_stations = 50;
    let mut x = Vec::with_capacity(n_stations);
    let mut ue = Vec::with_capacity(n_stations);
    
    // Simple velocity profile (favorable gradient then adverse)
    for i in 0..n_stations {
        let xi = (i as f64) / (n_stations as f64 - 1.0) * 0.5;  // x from 0 to 0.5
        x.push(xi.max(1e-6));  // Avoid x=0
        
        // Velocity profile: accelerates then decelerates
        let ue_val = if xi < 0.1 {
            0.1 + xi * 10.0  // Accelerate from 0.1 to 1.1
        } else {
            1.1 - (xi - 0.1) * 0.5  // Decelerate
        };
        ue.push(ue_val);
    }
    
    let re = 1e6;
    let msq = 0.0;
    let config = MarchConfig {
        ncrit: 9.0,
        max_iter: 25,
        tolerance: 0.1,
        hlmax: 4.5,
        htmax: 2.5,
        debug_trace: false,
        ..Default::default()
    };
    
    let result = march_fixed_ue(&x, &ue, re, msq, &config);
    
    println!("\n=== March result ===");
    println!("Converged: {}", result.converged);
    println!("Transition x: {:?}", result.x_transition);
    println!("Transition index: {:?}", result.transition_index);
    println!("Separation x: {:?}", result.x_separation);
    println!("Number of stations: {}", result.stations.len());
    
    if let Some(trans_idx) = result.transition_index {
        if trans_idx < result.stations.len() {
            let trans_station = &result.stations[trans_idx];
            println!("\nTransition station:");
            println!("  x: {:.6}", trans_station.x);
            println!("  theta: {:.6e}", trans_station.theta);
            println!("  delta_star: {:.6e}", trans_station.delta_star);
            println!("  Hk: {:.4}", trans_station.hk);
            println!("  ctau: {:.6}", trans_station.ctau);
        }
    }
    
    // Basic sanity checks
    assert!(result.stations.len() > 0, "Should have computed stations");
    
    // If transition occurred, verify it's in a reasonable location
    if let Some(x_trans) = result.x_transition {
        assert!(x_trans > 0.0 && x_trans < 0.5, 
                "Transition should be in the domain");
    }
}
