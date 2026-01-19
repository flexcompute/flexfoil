//! Polar generation utility for XFOIL comparison.
//!
//! Usage: cargo run --release --example polar_gen -- <naca> <reynolds> <alpha_start> <alpha_end> <alpha_step>
//!
//! Example: cargo run --release --example polar_gen -- 0012 3000000 -4 14 1

use rustfoil_core::{naca::naca4, Body, CubicSpline, PanelingParams};
use rustfoil_solver::inviscid::FlowConditions;
use rustfoil_solver::viscous::{ViscousSolver, ViscousConfig};
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        eprintln!("Usage: polar_gen <naca> <reynolds> <alpha_start> <alpha_end> <alpha_step>");
        eprintln!("Example: polar_gen 0012 3000000 -4 14 1");
        std::process::exit(1);
    }
    
    let designation: u32 = args[1].parse().expect("Invalid NACA designation");
    let reynolds: f64 = args[2].parse().expect("Invalid Reynolds number");
    let alpha_start: f64 = args[3].parse().expect("Invalid alpha_start");
    let alpha_end: f64 = args[4].parse().expect("Invalid alpha_end");
    let alpha_step: f64 = args[5].parse().expect("Invalid alpha_step");
    
    // Generate airfoil with XFOIL-style paneling
    let buffer = naca4(designation, Some(123));
    let spline = CubicSpline::from_points(&buffer).expect("Failed to create spline");
    let params = PanelingParams::default();
    let paneled = spline.resample_xfoil(160, &params);
    let airfoil = Body::from_points("test", &paneled).expect("Failed to create body");
    
    // Setup viscous solver
    let config = ViscousConfig {
        reynolds,
        n_crit: 9.0,
        ..Default::default()
    };
    let solver = ViscousSolver::new(config);
    
    // Run polar
    let mut alpha = alpha_start;
    let mut alphas = Vec::new();
    let mut cls = Vec::new();
    let mut cds = Vec::new();
    let mut cms = Vec::new();
    let mut xtr_tops = Vec::new();
    let mut xtr_bots = Vec::new();
    
    while alpha <= alpha_end + 0.001 {
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        alphas.push(alpha);
        cls.push(result.cl);
        cds.push(result.cd);
        cms.push(result.cm);
        xtr_tops.push(result.x_tr_upper);
        xtr_bots.push(result.x_tr_lower);
        
        alpha += alpha_step;
    }
    
    // Output JSON
    println!("{{");
    println!("  \"alpha\": {:?},", alphas);
    println!("  \"cl\": {:?},", cls);
    println!("  \"cd\": {:?},", cds);
    println!("  \"cm\": {:?},", cms);
    println!("  \"xtr_top\": {:?},", xtr_tops);
    println!("  \"xtr_bot\": {:?}", xtr_bots);
    println!("}}");
}
