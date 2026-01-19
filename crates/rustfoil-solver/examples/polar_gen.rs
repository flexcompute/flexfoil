//! Polar generation utility for XFOIL comparison.
//!
//! Usage: cargo run --release --example polar_gen -- <naca> <reynolds> <alpha_start> <alpha_end> <alpha_step>
//!
//! Example: cargo run --release --example polar_gen -- 0012 3000000 -4 14 1

use rustfoil_core::{naca::naca4, Body, CubicSpline, PanelingParams};
use rustfoil_solver::inviscid::FlowConditions;
use rustfoil_solver::viscous::{ViscousSolver, ViscousConfig};
use serde_json::json;
use std::env;

fn sanitize(values: &[f64]) -> Vec<Option<f64>> {
    values
        .iter()
        .map(|v| if v.is_finite() { Some(*v) } else { None })
        .collect()
}

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
    
    // Setup viscous solver (Full Newton, XFOIL-style)
    let mut config = ViscousConfig::newton_xfoil(reynolds, 9.0);
    config.max_iterations = 30;
    config.tolerance = 1e-5;
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
    let output = json!({
        "alpha": sanitize(&alphas),
        "cl": sanitize(&cls),
        "cd": sanitize(&cds),
        "cm": sanitize(&cms),
        "xtr_top": sanitize(&xtr_tops),
        "xtr_bot": sanitize(&xtr_bots),
    });
    println!("{}", serde_json::to_string_pretty(&output).expect("Failed to serialize JSON"));
}
