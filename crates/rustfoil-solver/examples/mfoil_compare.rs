//! Comparison between mfoil-style solver and other coupling methods.
//!
//! Usage: cargo run --release --example mfoil_compare

use rustfoil_core::{naca::naca4, Body, CubicSpline, PanelingParams};
use rustfoil_solver::inviscid::FlowConditions;
use rustfoil_solver::viscous::{ViscousSolver, ViscousConfig};
use rustfoil_solver::mfoil_bl::{MfoilSolver, MfoilConfig};

fn main() {
    // Test case: NACA 0012 at Re=1e6, alpha=5°
    let designation = 12u32;
    let reynolds = 1_000_000.0;
    let alpha = 5.0;
    
    // Generate airfoil
    let buffer = naca4(designation, Some(123));
    let spline = CubicSpline::from_points(&buffer).expect("Failed to create spline");
    let params = PanelingParams::default();
    let paneled = spline.resample_xfoil(160, &params);
    let airfoil = Body::from_points("NACA0012", &paneled).expect("Failed to create body");
    
    println!("=========================================");
    println!("Mfoil-style BL Solver Comparison");
    println!("NACA 00{} at Re = {:.0e}, α = {:.1}°", designation, reynolds, alpha);
    println!("=========================================\n");
    
    // Test 1: Transpiration coupling (current default)
    {
        println!("--- Transpiration Coupling ---");
        let config = ViscousConfig::with_transpiration(reynolds);
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        println!("  Cl = {:.6}", result.cl);
        println!("  Cd = {:.6}", result.cd);
        println!("  Cm = {:.6}", result.cm);
        println!("  xtr_upper = {:.4}", result.x_tr_upper);
        println!("  xtr_lower = {:.4}", result.x_tr_lower);
        println!("  converged = {}, iterations = {}", result.converged, result.iterations);
        println!();
    }
    
    // Test 2: Full Newton coupling
    {
        println!("--- Full Newton Coupling ---");
        let config = ViscousConfig::newton_xfoil(reynolds, 9.0);
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        println!("  Cl = {:.6}", result.cl);
        println!("  Cd = {:.6}", result.cd);
        println!("  Cm = {:.6}", result.cm);
        println!("  xtr_upper = {:.4}", result.x_tr_upper);
        println!("  xtr_lower = {:.4}", result.x_tr_lower);
        println!("  converged = {}, iterations = {}", result.converged, result.iterations);
        println!();
    }
    
    // Test 3: Mfoil Newton coupling via ViscousSolver
    {
        println!("--- Mfoil Newton (via ViscousSolver) ---");
        let config = ViscousConfig::mfoil_style(reynolds, 9.0);
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        println!("  Cl = {:.6}", result.cl);
        println!("  Cd = {:.6}", result.cd);
        println!("  Cm = {:.6}", result.cm);
        println!("  xtr_upper = {:.4}", result.x_tr_upper);
        println!("  xtr_lower = {:.4}", result.x_tr_lower);
        println!("  converged = {}, iterations = {}", result.converged, result.iterations);
        println!();
    }
    
    // Test 4: Direct MfoilSolver usage
    {
        println!("--- MfoilSolver Direct ---");
        let coords: Vec<[f64; 2]> = airfoil.panels().iter()
            .map(|p| [p.p1.x, p.p1.y])
            .collect();
        
        let mfoil_config = MfoilConfig {
            reynolds,
            ncrit: 9.0,
            max_iter: 100,
            rtol: 1e-10,
            verbose: 1,
            ..Default::default()
        };
        let mut mfoil_solver = MfoilSolver::new(mfoil_config);
        mfoil_solver.set_airfoil(&coords);
        let result = mfoil_solver.solve(alpha, true);
        
        println!("  Cl = {:.6}", result.cl);
        println!("  Cd = {:.6}", result.cd);
        println!("  Cm = {:.6}", result.cm);
        println!("  xtr_upper = {:.4}", result.xtr_upper);
        println!("  xtr_lower = {:.4}", result.xtr_lower);
        println!("  converged = {}, iterations = {}", result.converged, result.iterations);
        println!();
    }
    
    println!("=========================================");
    println!("Note: For accurate comparison with mfoil.py,");
    println!("run scripts/xfoil_comparison.py");
    println!("=========================================");
}
