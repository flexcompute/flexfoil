use rustfoil_inviscid::{InviscidSolver, FlowConditions};
use std::f64::consts::PI;

#[test]
fn test_cl_slope() {
    let content = std::fs::read_to_string("../../Xfoil-instrumented/bin/naca0012_xfoil.dat")
        .expect("Failed to read file");
    
    let coords: Vec<(f64, f64)> = content
        .lines()
        .skip(1)
        .filter_map(|line| {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                Some((parts[0].parse::<f64>().ok()?, parts[1].parse::<f64>().ok()?))
            } else { None }
        })
        .collect();
    
    println!("\nLoaded {} points", coords.len());
    
    let solver = InviscidSolver::new();
    let factorized = solver.factorize(&coords).expect("Factorization failed");
    
    println!("\n{:>6} {:>12} {:>12}", "Alpha", "CL", "CLα (per rad)");
    println!("{}", "-".repeat(34));
    
    let alphas = [0.0, 1.0, 2.0, 4.0];
    let mut prev_cl = 0.0;
    let mut prev_alpha = 0.0;
    
    for alpha_deg in alphas {
        let flow = FlowConditions::with_alpha_deg(alpha_deg);
        let solution = factorized.solve_alpha(&flow);
        
        if alpha_deg > 0.0 {
            let dalpha = (alpha_deg - prev_alpha).to_radians();
            let dcl = solution.cl - prev_cl;
            let cl_alpha = dcl / dalpha;
            println!("{:>6.1} {:>12.4} {:>12.4}", alpha_deg, solution.cl, cl_alpha);
        } else {
            println!("{:>6.1} {:>12.4} {:>12}", alpha_deg, solution.cl, "-");
        }
        
        prev_cl = solution.cl;
        prev_alpha = alpha_deg;
    }
    
    println!("\nThin airfoil theory: CLα = 2π = {:.4}", 2.0 * PI);
    
    // Compute CLα from 0 to 4 degrees
    let cl0 = factorized.solve_alpha(&FlowConditions::with_alpha_deg(0.0)).cl;
    let cl4 = factorized.solve_alpha(&FlowConditions::with_alpha_deg(4.0)).cl;
    let cl_alpha = (cl4 - cl0) / 4.0_f64.to_radians();
    
    println!("\nOverall CLα (0-4°): {:.4}", cl_alpha);
    println!("Ratio to 2π: {:.3}", cl_alpha / (2.0 * PI));
    
    // Check Kutta condition at α=2°
    let sol2 = factorized.solve_alpha(&FlowConditions::with_alpha_deg(2.0));
    let n = sol2.gamma.len();
    println!("\nKutta condition at α=2°: γ[0] + γ[{}] = {:.6e}", n-1, sol2.gamma[0] + sol2.gamma[n-1]);
    
    // Check gamma at TE
    println!("γ[0] = {:.6}", sol2.gamma[0]);
    println!("γ[{}] = {:.6}", n-1, sol2.gamma[n-1]);
}
