//! Debug curvature computation for paneling comparison.
//!
//! Usage: cargo run -p rustfoil-core --example debug_curvature

use rustfoil_core::spline::{CubicSpline, PanelingParams};
use rustfoil_core::xfoil_spline::XfoilSpline;
use rustfoil_core::point;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

fn load_dat_file(path: &Path) -> Vec<(f64, f64)> {
    let file = File::open(path).expect("Failed to open file");
    let reader = BufReader::new(file);
    
    let mut coords = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                coords.push((x, y));
            }
        }
    }
    coords
}

fn main() {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().parent().unwrap();
    
    for foil in &["naca0012", "naca2412", "naca4412"] {
        println!("\n============================================================");
        println!("=== {} ===", foil.to_uppercase());
        println!("============================================================");
        
        let buffer_path = workspace.join(format!("{}_buffer_xfoil.dat", foil));
        if !buffer_path.exists() {
            println!("  Buffer not found, skipping");
            continue;
        }
        
        let coords = load_dat_file(&buffer_path);
        let points: Vec<_> = coords.iter().map(|(x, y)| point::point(*x, *y)).collect();
        
        let n = coords.len();
        println!("Buffer: {} points", n);
        
        // Create XfoilSpline
        let xfoil_spline = XfoilSpline::from_points(&points).expect("Failed to create spline");
        
        // Compute arc lengths
        let s_max = xfoil_spline.total_arc_length();
        let sbref = s_max / 2.0;
        
        // Find LE
        let s_le = xfoil_spline.lefind();
        let le_pt = xfoil_spline.seval(s_le);
        let cv_le = xfoil_spline.curvature(s_le).abs() * sbref;
        
        println!("s_max: {:.8}", s_max);
        println!("s_le:  {:.8} ({:.4}%)", s_le, 100.0 * s_le / s_max);
        println!("LE point: ({:.8}, {:.8})", le_pt.x, le_pt.y);
        println!("cv_le: {:.4}", cv_le);
        
        // Sample curvature at buffer points
        println!("\nCurvature at key buffer points (near LE):");
        let mut s_buffer = vec![0.0f64];
        for i in 1..n {
            let ds = ((coords[i].0 - coords[i-1].0).powi(2) + 
                      (coords[i].1 - coords[i-1].1).powi(2)).sqrt();
            s_buffer.push(s_buffer[i-1] + ds);
        }
        
        // Find buffer LE index
        let le_idx = s_buffer.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((**a - s_le).abs()).partial_cmp(&((**b - s_le).abs())).unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(n / 2);
        
        println!("Buffer LE index: {}", le_idx);
        
        for i in le_idx.saturating_sub(5)..=(le_idx + 5).min(n - 1) {
            let cv = xfoil_spline.curvature(s_buffer[i]).abs() * sbref;
            let marker = if i == le_idx { " <-- LE" } else { "" };
            println!("  [{:3}] s={:.6}, x={:.6}, y={:+.6}, cv={:.4}{}",
                     i, s_buffer[i], coords[i].0, coords[i].1, cv, marker);
        }
        
        // Check symmetry of curvature around LE
        println!("\nCurvature symmetry check (offset from LE):");
        for offset in 1..=5 {
            let s_minus = s_le - 0.01 * offset as f64;
            let s_plus = s_le + 0.01 * offset as f64;
            let cv_minus = xfoil_spline.curvature(s_minus).abs() * sbref;
            let cv_plus = xfoil_spline.curvature(s_plus).abs() * sbref;
            println!("  offset={}: cv(-{:.2})={:.4}, cv(+{:.2})={:.4}, diff={:.4}",
                     offset, 0.01 * offset as f64, cv_minus, 
                     0.01 * offset as f64, cv_plus, cv_plus - cv_minus);
        }
    }
}
