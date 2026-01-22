//! Check IBL=3 step residual

use rustfoil_bl::equations::{blvar, bldif, FlowType};
use rustfoil_bl::state::BlStation;

#[test]
fn test_ibl3_step_residual() {
    let re = 1e6;
    let msq = 0.0;

    // Station 1 (prev, IBL=2 converged = our station 1)
    let mut s1 = BlStation::new();
    s1.x = 0.004704;
    s1.u = 0.258499;
    s1.theta = 3.667001e-5;  // IBL=2 final
    s1.delta_star = 8.084212e-5;  // IBL=2 final
    s1.ampl = 0.0;
    s1.ctau = 0.03;
    s1.is_laminar = true;
    blvar(&mut s1, FlowType::Laminar, msq, re);

    // Station 2 (curr, IBL=3 initial = our station 2)
    let mut s2 = BlStation::new();
    s2.x = 0.006961;
    s2.u = 0.416645;
    s2.theta = 3.667001e-5;  // Initial guess = prev final
    s2.delta_star = 8.084212e-5;
    s2.ampl = 0.0;
    s2.ctau = 0.03;
    s2.is_laminar = true;
    blvar(&mut s2, FlowType::Laminar, msq, re);

    // Call bldif
    let (res, _jac) = bldif(&s1, &s2, FlowType::Laminar, msq, re);

    println!("\n=== IBL=3 step (Station 1 → Station 2) ===");
    println!("Prev (IBL=2 final): x={:.6}, θ={:.6e}", s1.x, s1.theta);
    println!("Curr (IBL=3 init):  x={:.6}, θ={:.6e}", s2.x, s2.theta);
    
    println!("\n=== Residuals ===");
    println!("res_mom (RustFoil) = {:.4}", res.res_mom);
    println!("res_mom (XFOIL IBL=3) = -0.130");
    
    let error = (res.res_mom - (-0.130)).abs();
    let pct_error = error / 0.130 * 100.0;
    println!("\nError: {:.4} ({:.1}%)", error, pct_error);
}
