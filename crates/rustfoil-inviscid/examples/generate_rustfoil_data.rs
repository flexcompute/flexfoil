//! Generate RustFoil data for comparison with XFOIL
//! Run with: cargo run --release --example generate_rustfoil_data

use rustfoil_inviscid::{FlowConditions, InviscidSolver};
use serde::Serialize;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

#[derive(Serialize)]
struct FoilData {
    foil: String,
    alphas: HashMap<i32, AlphaData>,
}

#[derive(Serialize)]
struct AlphaData {
    s: Vec<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
    gamma: Vec<f64>,
    cp: Vec<f64>,
    cl: f64,
    cm: f64,
}

fn load_dat_file(filename: &str) -> Option<Vec<(f64, f64)>> {
    let content = fs::read_to_string(filename).ok()?;
    let mut points = Vec::new();
    
    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                points.push((x, y));
            }
        }
    }
    
    if points.is_empty() { None } else { Some(points) }
}

fn main() {
    let foils = vec!["naca0012", "naca2412", "naca4412"];
    let alphas: Vec<i32> = vec![-4, -2, 0, 2, 4, 6, 8, 10];
    
    let mut all_data: HashMap<String, FoilData> = HashMap::new();
    
    for foil_name in &foils {
        println!("Processing {}...", foil_name);
        
        let dat_file = format!("{}_xfoil_paneled.dat", foil_name);
        let points = match load_dat_file(&dat_file) {
            Some(p) => p,
            None => {
                println!("  Could not load {}", dat_file);
                continue;
            }
        };
        
        let solver = InviscidSolver::new();
        let factorized = match solver.factorize(&points) {
            Ok(f) => f,
            Err(e) => {
                println!("  Factorization failed: {:?}", e);
                continue;
            }
        };
        
        let geom = factorized.geometry();
        
        let mut foil_data = FoilData {
            foil: foil_name.to_string(),
            alphas: HashMap::new(),
        };
        
        for &alpha in &alphas {
            let flow = FlowConditions::with_alpha_deg(alpha as f64);
            let solution = factorized.solve_alpha(&flow);
            
            let alpha_data = AlphaData {
                s: geom.s.clone(),
                x: geom.x.clone(),
                y: geom.y.clone(),
                gamma: solution.gamma.clone(),
                cp: solution.cp.clone(),
                cl: solution.cl,
                cm: solution.cm,
            };
            
            foil_data.alphas.insert(alpha, alpha_data);
            println!("  {} α={:+}°: CL={:.4}", foil_name, alpha, solution.cl);
        }
        
        all_data.insert(foil_name.to_string(), foil_data);
    }
    
    // Save to JSON
    let output_path = Path::new("scripts/comparison_data/rustfoil_data.json");
    let json = serde_json::to_string_pretty(&all_data).expect("Failed to serialize");
    fs::write(output_path, json).expect("Failed to write file");
    
    println!("\nSaved to {}", output_path.display());
}
