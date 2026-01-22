// Run with: cargo test --package rustfoil-coupling --lib -- debug_station_31 --nocapture

#[cfg(test)]
mod tests {
    use rustfoil_bl::equations::{blvar, bldif, FlowType};
    use rustfoil_bl::state::BlStation;

    #[test]
    fn debug_station_31() {
        // Station 30 (previous station) - XFOIL final values
        let mut s1 = BlStation::default();
        s1.x = 0.115205;
        s1.u = 1.503356;
        s1.theta = 2.073123e-04;
        s1.delta_star = 7.458165e-04;
        s1.ctau = 0.03;
        s1.ampl = 6.1372;
        
        // Station 31 - XFOIL initial (= station 30 final)
        let mut s2 = BlStation::default();
        s2.x = 0.127589;
        s2.u = 1.482947;
        s2.theta = 2.073123e-04;  // Initial guess = previous final
        s2.delta_star = 7.458165e-04;
        s2.ctau = 0.03;
        s2.ampl = 6.1372;
        
        let msq = 0.0;
        let re = 1e6;
        
        // Compute secondary variables
        blvar(&mut s1, FlowType::Laminar, msq, re);
        blvar(&mut s2, FlowType::Laminar, msq, re);
        
        println!("=== Station 30 (s1) - blvar output ===");
        println!("  H = {:.4}", s1.h);
        println!("  Hk = {:.4}", s1.hk);
        println!("  Hs = {:.4}", s1.hs);
        println!("  Cf = {:.6}", s1.cf);
        println!("  Cd = {:.6}", s1.cd);
        println!("  Rtheta = {:.2}", s1.r_theta);
        println!();
        
        println!("=== Station 31 (s2) initial - blvar output ===");
        println!("  H = {:.4}", s2.h);
        println!("  Hk = {:.4}", s2.hk);
        println!("  Hs = {:.4}", s2.hs);
        println!("  Cf = {:.6}", s2.cf);
        println!("  Cd = {:.6}", s2.cd);
        println!("  Rtheta = {:.2}", s2.r_theta);
        println!();
        
        // Compute residuals and Jacobian
        let (res, jac) = bldif(&s1, &s2, FlowType::Laminar, msq, re);
        
        println!("=== bldif output ===");
        println!("Residuals:");
        println!("  res_third (ampl) = {:.6e}", res.res_third);
        println!("  res_mom = {:.6e}", res.res_mom);
        println!("  res_shape = {:.6e}", res.res_shape);
        println!();
        
        println!("Jacobian VS2:");
        println!("  Row 0: [{:12.4}, {:12.4}, {:12.4}, {:12.4}, {:12.4}]",
            jac.vs2[0][0], jac.vs2[0][1], jac.vs2[0][2], jac.vs2[0][3], jac.vs2[0][4]);
        println!("  Row 1: [{:12.4}, {:12.4}, {:12.4}, {:12.4}, {:12.4}]",
            jac.vs2[1][0], jac.vs2[1][1], jac.vs2[1][2], jac.vs2[1][3], jac.vs2[1][4]);
        println!("  Row 2: [{:12.4}, {:12.4}, {:12.4}, {:12.4}, {:12.4}]",
            jac.vs2[2][0], jac.vs2[2][1], jac.vs2[2][2], jac.vs2[2][3], jac.vs2[2][4]);
        println!();
        
        println!("Key comparisons with XFOIL:");
        println!("                        XFOIL      RustFoil");
        println!("  VS2[1][1] (∂mom/∂θ):  4792.37    {:.2}", jac.vs2[1][1]);
        println!("  VS2[2][1] (∂shape/∂θ): 371.69    {:.2}", jac.vs2[2][1]);
        println!("  VS2[2][2] (∂shape/∂δ*): -52.99   {:.2}", jac.vs2[2][2]);
        println!();
        
        // Check the derivatives stored in s2
        println!("=== Derivative chain for δ* ===");
        println!("  h_delta_star = {:.4e}", s2.derivs.h_delta_star);
        println!("  hk_h = {:.4}", s2.derivs.hk_h);
        println!("  hs_hk = {:.4}", s2.derivs.hs_hk);
        println!("  cf_hk = {:.6}", s2.derivs.cf_hk);
        println!("  cd_hk = {:.6}", s2.derivs.cd_hk);
    }
}
