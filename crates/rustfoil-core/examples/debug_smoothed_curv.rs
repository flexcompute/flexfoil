//! Debug smoothed curvature for paneling comparison.
//!
//! This dumps the smoothed/normalized curvature values used in the Newton iteration.

use rustfoil_core::spline::{CubicSpline, PanelingParams};
use rustfoil_core::xfoil_spline::{XfoilSpline, Spline1D};
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
        
        let buffer_path = workspace.join(format!("{}_buffer_xfoil.dat", foil));
        if !buffer_path.exists() {
            println!("  Buffer not found, skipping");
            continue;
        }
        
        let coords = load_dat_file(&buffer_path);
        let points: Vec<_> = coords.iter().map(|(x, y)| point::point(*x, *y)).collect();
        let n = coords.len();
        
        // Create XfoilSpline and compute curvature
        let xfoil_spline = XfoilSpline::from_points(&points).expect("Failed to create spline");
        let s_max = xfoil_spline.total_arc_length();
        let sbref = s_max / 2.0;
        let s_le = xfoil_spline.lefind();
        
        // Compute arc lengths
        let mut s_buffer = vec![0.0f64];
        for i in 1..n {
            let ds = ((coords[i].0 - coords[i-1].0).powi(2) + 
                      (coords[i].1 - coords[i-1].1).powi(2)).sqrt();
            s_buffer.push(s_buffer[i-1] + ds);
        }
        
        // Compute raw curvature
        let mut curv_buffer: Vec<f64> = s_buffer.iter()
            .map(|&s| xfoil_spline.curvature(s).abs() * sbref)
            .collect();
        
        // Find buffer LE index
        let le_idx = s_buffer.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((**a - s_le).abs()).partial_cmp(&((**b - s_le).abs())).unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(n / 2);
        
        // Compute LE curvature
        let cv_le = xfoil_spline.curvature(s_le).abs() * sbref;
        
        // Average curvature near LE
        let nk = 3;
        let mut cv_sum = 0.0;
        for k in -nk..=nk {
            let frac = k as f64 / nk as f64;
            let sbk = s_le + frac * sbref / cv_le.max(20.0);
            let cvk = xfoil_spline.curvature(sbk).abs() * sbref;
            cv_sum += cvk;
        }
        let cv_avg = (cv_sum / (2 * nk + 1) as f64).max(10.0);
        
        // Compute smoothing length
        let smool = (1.0 / cv_avg.max(20.0)).max(0.25 / (160.0 / 2.0));
        
        println!("s_le: {:.6}, cv_le: {:.2}, cv_avg: {:.2}", s_le, cv_le, cv_avg);
        println!("smool: {:.6}", smool);
        
        // Print raw curvature near LE
        println!("\nRaw curvature near LE (buffer indices):");
        for i in le_idx.saturating_sub(3)..=(le_idx + 3).min(n - 1) {
            let marker = if i == le_idx { " <-- LE" } else { "" };
            println!("  [{:3}] s={:.6}, raw_cv={:.4}{}",
                     i, s_buffer[i], curv_buffer[i], marker);
        }
        
        // Note: We can't easily dump the smoothed curvature without modifying the library,
        // but we can compute what the raw curvature looks like at key points
        println!("\nCurvature asymmetry around s_le:");
        for offset in [0.005, 0.01, 0.02, 0.03] {
            let cv_minus = xfoil_spline.curvature(s_le - offset).abs() * sbref;
            let cv_plus = xfoil_spline.curvature(s_le + offset).abs() * sbref;
            let ratio = cv_plus / cv_minus;
            println!("  offset={:.3}: cv(-)={:.2}, cv(+)={:.2}, ratio={:.3}",
                     offset, cv_minus, cv_plus, ratio);
        }
    }
}
