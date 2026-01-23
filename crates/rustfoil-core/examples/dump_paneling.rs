//! Dump XFOIL-style paneling coordinates for comparison with XFOIL.
//!
//! Usage: cargo run -p rustfoil-core --example dump_paneling

use rustfoil_core::spline::{CubicSpline, PanelingParams};
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

fn save_dat_file(path: &Path, coords: &[(f64, f64)], header: &str) {
    let mut file = File::create(path).expect("Failed to create file");
    writeln!(file, "{}", header).expect("Failed to write header");
    for (x, y) in coords {
        writeln!(file, "  {:.8}  {:.8}", x, y).expect("Failed to write coords");
    }
    println!("Saved {} points to {:?}", coords.len(), path);
}

fn main() {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().parent().unwrap();
    
    let foils = ["naca0012", "naca2412", "naca4412"];
    
    for foil in &foils {
        println!("\n=== {} ===", foil.to_uppercase());
        
        // Load buffer coordinates (input to paneling)
        let buffer_path = workspace.join(format!("{}_buffer.dat", foil));
        if !buffer_path.exists() {
            // Try the real buffer file
            let alt_path = workspace.join(format!("{}_buffer_real.dat", foil));
            if !alt_path.exists() {
                println!("  Buffer file not found, skipping");
                continue;
            }
        }
        
        // Prefer the XFOIL-generated buffer file (from PCOP command)
        let buffer_path = if workspace.join(format!("{}_buffer_xfoil.dat", foil)).exists() {
            workspace.join(format!("{}_buffer_xfoil.dat", foil))
        } else if workspace.join(format!("{}_buffer_real.dat", foil)).exists() {
            workspace.join(format!("{}_buffer_real.dat", foil))
        } else {
            workspace.join(format!("{}_buffer.dat", foil))
        };
        
        let buffer_coords = load_dat_file(&buffer_path);
        println!("  Loaded {} buffer points from {:?}", buffer_coords.len(), buffer_path);
        
        // Load XFOIL paneled output for comparison (prefer freshly generated)
        let xfoil_path = if workspace.join(format!("{}_paneled_xfoil.dat", foil)).exists() {
            workspace.join(format!("{}_paneled_xfoil.dat", foil))
        } else {
            workspace.join(format!("{}_xfoil_paneled.dat", foil))
        };
        let xfoil_coords = if xfoil_path.exists() {
            let coords = load_dat_file(&xfoil_path);
            println!("  Loaded {} XFOIL paneled points", coords.len());
            Some(coords)
        } else {
            println!("  XFOIL paneled file not found");
            None
        };
        
        // Convert to points for CubicSpline
        let points: Vec<_> = buffer_coords.iter()
            .map(|(x, y)| rustfoil_core::point(*x, *y))
            .collect();
        
        // Create spline and resample with XFOIL-style paneling
        let spline = CubicSpline::from_points(&points).expect("Failed to create spline");
        
        let n_panels = xfoil_coords.as_ref().map(|c| c.len()).unwrap_or(160);
        let params = PanelingParams::default();
        let repaneled = spline.resample_xfoil(n_panels, &params);
        
        println!("  Generated {} RustFoil paneled points", repaneled.len());
        
        // Save RustFoil output
        let rustfoil_path = workspace.join(format!("{}_rustfoil_paneled.dat", foil));
        let rustfoil_coords: Vec<_> = repaneled.iter().map(|p| (p.x, p.y)).collect();
        save_dat_file(&rustfoil_path, &rustfoil_coords, &format!("{} (RustFoil paneled)", foil));
        
        // Compare if XFOIL data available
        if let Some(xfoil) = &xfoil_coords {
            if xfoil.len() == rustfoil_coords.len() {
                let mut max_err = 0.0f64;
                let mut sum_sq = 0.0f64;
                let mut max_idx = 0;
                
                for (i, ((xf_x, xf_y), (rf_x, rf_y))) in xfoil.iter().zip(rustfoil_coords.iter()).enumerate() {
                    let err = ((xf_x - rf_x).powi(2) + (xf_y - rf_y).powi(2)).sqrt();
                    sum_sq += err * err;
                    if err > max_err {
                        max_err = err;
                        max_idx = i;
                    }
                }
                
                let rms = (sum_sq / xfoil.len() as f64).sqrt();
                println!("\n  Comparison with XFOIL:");
                println!("    RMS error:  {:.2e}", rms);
                println!("    Max error:  {:.2e} at point {}", max_err, max_idx);
                
                // Show first few and last few points
                println!("\n  First 5 points:");
                println!("    {:>4}  {:>12} {:>12}  {:>12} {:>12}  {:>10}", 
                         "idx", "xfoil_x", "xfoil_y", "rustfoil_x", "rustfoil_y", "error");
                for i in 0..5.min(xfoil.len()) {
                    let err = ((xfoil[i].0 - rustfoil_coords[i].0).powi(2) + 
                               (xfoil[i].1 - rustfoil_coords[i].1).powi(2)).sqrt();
                    println!("    {:>4}  {:>12.8} {:>12.8}  {:>12.8} {:>12.8}  {:>10.2e}",
                             i, xfoil[i].0, xfoil[i].1, 
                             rustfoil_coords[i].0, rustfoil_coords[i].1, err);
                }
                
                // LE region
                let le_idx = rustfoil_coords.iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.0.partial_cmp(&b.0).unwrap())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                
                println!("\n  LE region (idx {}-{}):", le_idx.saturating_sub(2), le_idx + 2);
                for i in le_idx.saturating_sub(2)..=(le_idx + 2).min(xfoil.len() - 1) {
                    let err = ((xfoil[i].0 - rustfoil_coords[i].0).powi(2) + 
                               (xfoil[i].1 - rustfoil_coords[i].1).powi(2)).sqrt();
                    println!("    {:>4}  {:>12.8} {:>12.8}  {:>12.8} {:>12.8}  {:>10.2e}",
                             i, xfoil[i].0, xfoil[i].1, 
                             rustfoil_coords[i].0, rustfoil_coords[i].1, err);
                }
            }
        }
    }
    
    println!("\n\nOutput files saved. Use these for interactive comparison.");
}
