//! Debug tests for lefind algorithm comparison.

#[cfg(test)]
mod lefind_debug {
    use crate::xfoil_spline::XfoilSpline;
    use crate::point;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    fn load_buffer(path: &Path) -> Vec<(f64, f64)> {
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

    #[test]
    fn debug_lefind_initial_guess() {
        let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().parent().unwrap();
        
        for foil in &["naca0012", "naca2412", "naca4412"] {
            let buffer_path = workspace.join(format!("{}_buffer_xfoil.dat", foil));
            if !buffer_path.exists() {
                println!("Skipping {}: buffer not found", foil);
                continue;
            }
            
            let coords = load_buffer(&buffer_path);
            let n = coords.len();
            
            println!("\n=== {} ({} points) ===", foil.to_uppercase(), n);
            
            // Build spline
            let points: Vec<_> = coords.iter().map(|(x, y)| point::point(*x, *y)).collect();
            let spline = XfoilSpline::from_points(&points).expect("Failed to create spline");
            
            // TE point
            let x_te = 0.5 * (coords[0].0 + coords[n-1].0);
            let y_te = 0.5 * (coords[0].1 + coords[n-1].1);
            println!("TE: ({:.8}, {:.8})", x_te, y_te);
            
            // Find initial guess - trace through the loop
            let mut i_le_initial = n / 2;
            let mut found = false;
            
            println!("\nInitial guess search (looking for DOTP < 0):");
            for i in 2..n-2 {
                let dx_te = coords[i].0 - x_te;
                let dy_te = coords[i].1 - y_te;
                let dx = coords[i + 1].0 - coords[i].0;
                let dy = coords[i + 1].1 - coords[i].1;
                let dotp = dx_te * dx + dy_te * dy;
                
                // Print around the transition
                if !found && dotp < 0.0 {
                    // Print a few before and after
                    for j in i.saturating_sub(3)..=(i+3).min(n-3) {
                        let dx_te_j = coords[j].0 - x_te;
                        let dy_te_j = coords[j].1 - y_te;
                        let dx_j = coords[j + 1].0 - coords[j].0;
                        let dy_j = coords[j + 1].1 - coords[j].1;
                        let dotp_j = dx_te_j * dx_j + dy_te_j * dy_j;
                        let marker = if j == i { " <-- BREAK" } else { "" };
                        println!("  i={:3}: x={:.6}, y={:+.6}, dotp={:+.6e}{}",
                                 j, coords[j].0, coords[j].1, dotp_j, marker);
                    }
                    i_le_initial = i;
                    found = true;
                    break;
                }
            }
            
            if !found {
                println!("  WARNING: Loop completed without finding DOTP < 0!");
                println!("  Default i_le = {}", i_le_initial);
            }
            
            // Now run Newton iteration
            let s_le = spline.lefind();
            let le_pt = spline.seval(s_le);
            
            // Find closest index in original coords
            let mut closest_idx = 0;
            let mut min_dist = f64::MAX;
            for (i, (x, y)) in coords.iter().enumerate() {
                let d = (x - le_pt.x).powi(2) + (y - le_pt.y).powi(2);
                if d < min_dist {
                    min_dist = d;
                    closest_idx = i;
                }
            }
            
            println!("\nResults:");
            println!("  Initial guess index: {}", i_le_initial);
            println!("  Final s_le: {:.8}", s_le);
            println!("  Final LE point: ({:.8}, {:.8})", le_pt.x, le_pt.y);
            println!("  Closest coord index: {}", closest_idx);
            println!("  Min x index: {}", coords.iter().enumerate()
                     .min_by(|(_, a), (_, b)| a.0.partial_cmp(&b.0).unwrap())
                     .map(|(i, _)| i).unwrap_or(0));
        }
    }
}
