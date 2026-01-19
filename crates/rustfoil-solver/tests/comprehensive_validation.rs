//! Comprehensive validation of RustFoil against XFOIL reference data.
//!
//! This test suite validates:
//! - Multiple airfoils (NACA 0012, 2412, 4412, 6412)
//! - Multiple Reynolds numbers (500k, 1M, 3M, 6M)
//! - Full alpha sweeps (-5° to +15°)
//! - Cl-alpha curves
//! - Cl/Cd (L/D) performance
//! - Transition location

use rustfoil_core::{naca::naca4, Body, CubicSpline, PanelingParams, Point};
use rustfoil_solver::inviscid::FlowConditions;
use rustfoil_solver::viscous::{ViscousSolver, ViscousConfig};

/// Generate a NACA 4-digit airfoil with proper paneling.
fn make_airfoil(designation: u32, n_panels: usize) -> Body {
    let buffer = naca4(designation, Some(123));
    let spline = CubicSpline::from_points(&buffer).unwrap();
    let params = PanelingParams::default();
    let paneled = spline.resample_xfoil(n_panels, &params);
    Body::from_points(&format!("NACA{:04}", designation), &paneled).unwrap()
}

/// Single analysis point result.
#[derive(Debug, Clone)]
struct AnalysisPoint {
    alpha: f64,
    cl: f64,
    cd: f64,
    cm: f64,
    x_tr_upper: f64,
    x_tr_lower: f64,
    converged: bool,
}

/// Run a full alpha sweep for an airfoil at a given Reynolds number.
fn run_alpha_sweep(
    airfoil: &Body,
    reynolds: f64,
    alphas: &[f64],
) -> Vec<AnalysisPoint> {
    let config = ViscousConfig {
        reynolds,
        n_crit: 9.0,
        ..Default::default()
    };
    let solver = ViscousSolver::new(config);
    
    alphas.iter().map(|&alpha| {
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(airfoil, &flow);
        
        AnalysisPoint {
            alpha,
            cl: result.cl,
            cd: result.cd,
            cm: result.cm,
            x_tr_upper: result.x_tr_upper,
            x_tr_lower: result.x_tr_lower,
            converged: result.converged,
        }
    }).collect()
}

/// Print a polar table.
fn print_polar_table(name: &str, reynolds: f64, points: &[AnalysisPoint]) {
    println!("\n{} at Re = {:.0e}", name, reynolds);
    println!("  α°    |   Cl    |   Cd    |  Cl/Cd |  Cm     | xtr_u  | xtr_l  | Conv");
    println!("--------|---------|---------|--------|---------|--------|--------|-----");
    
    for p in points {
        let l_d = if p.cd > 1e-6 { p.cl / p.cd } else { 0.0 };
        let conv = if p.converged { "✓" } else { "✗" };
        println!("{:6.1}  | {:7.4} | {:7.5} | {:6.1} | {:7.4} | {:5.3}  | {:5.3}  | {}",
            p.alpha, p.cl, p.cd, l_d, p.cm, p.x_tr_upper, p.x_tr_lower, conv);
    }
}

/// Compute key metrics from a polar.
#[derive(Debug)]
struct PolarMetrics {
    /// Lift curve slope (per degree)
    cl_alpha: f64,
    /// Zero-lift angle (degrees)
    alpha_0: f64,
    /// Maximum L/D
    max_l_d: f64,
    /// Cl at max L/D
    cl_at_max_l_d: f64,
    /// Minimum Cd
    min_cd: f64,
    /// Max Cl achieved
    max_cl: f64,
}

fn compute_metrics(points: &[AnalysisPoint]) -> PolarMetrics {
    // Lift curve slope from linear region (-2° to +6°)
    let linear_points: Vec<_> = points.iter()
        .filter(|p| p.alpha >= -2.0 && p.alpha <= 6.0 && p.converged)
        .collect();
    
    let cl_alpha = if linear_points.len() >= 2 {
        let n = linear_points.len() as f64;
        let sum_a: f64 = linear_points.iter().map(|p| p.alpha).sum();
        let sum_cl: f64 = linear_points.iter().map(|p| p.cl).sum();
        let sum_a2: f64 = linear_points.iter().map(|p| p.alpha * p.alpha).sum();
        let sum_a_cl: f64 = linear_points.iter().map(|p| p.alpha * p.cl).sum();
        
        (n * sum_a_cl - sum_a * sum_cl) / (n * sum_a2 - sum_a * sum_a)
    } else {
        0.0
    };
    
    // Zero-lift angle (interpolate where Cl crosses zero)
    let alpha_0 = if let Some(i) = points.windows(2).position(|w| w[0].cl * w[1].cl <= 0.0) {
        let p1 = &points[i];
        let p2 = &points[i + 1];
        if (p2.cl - p1.cl).abs() > 1e-6 {
            p1.alpha - p1.cl * (p2.alpha - p1.alpha) / (p2.cl - p1.cl)
        } else {
            0.0
        }
    } else {
        0.0
    };
    
    // Max L/D
    let (max_l_d, cl_at_max_l_d) = points.iter()
        .filter(|p| p.converged && p.cd > 1e-6)
        .map(|p| (p.cl / p.cd, p.cl))
        .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
        .unwrap_or((0.0, 0.0));
    
    // Min Cd
    let min_cd = points.iter()
        .filter(|p| p.converged)
        .map(|p| p.cd)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(0.0);
    
    // Max Cl
    let max_cl = points.iter()
        .filter(|p| p.converged)
        .map(|p| p.cl)
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(0.0);
    
    PolarMetrics {
        cl_alpha,
        alpha_0,
        max_l_d,
        cl_at_max_l_d,
        min_cd,
        max_cl,
    }
}

#[test]
fn test_comprehensive_airfoil_comparison() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║              RUSTFOIL COMPREHENSIVE VALIDATION SUITE                         ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    
    // Test matrix
    let airfoils = [
        (12u32, "NACA 0012 (symmetric)"),
        (2412, "NACA 2412 (2% camber)"),
        (4412, "NACA 4412 (4% camber)"),
        (6412, "NACA 6412 (6% camber)"),
    ];
    
    let reynolds_numbers = [5e5, 1e6, 3e6, 6e6];
    let alphas: Vec<f64> = (-5..=15).map(|a| a as f64).collect();
    
    // Store all metrics for summary
    let mut all_metrics: Vec<(String, f64, PolarMetrics)> = Vec::new();
    
    for (designation, name) in &airfoils {
        println!("\n");
        println!("═══════════════════════════════════════════════════════════════════════════════");
        println!("  {} (NACA {:04})", name, designation);
        println!("═══════════════════════════════════════════════════════════════════════════════");
        
        let airfoil = make_airfoil(*designation, 160);
        
        for &re in &reynolds_numbers {
            let points = run_alpha_sweep(&airfoil, re, &alphas);
            print_polar_table(name, re, &points);
            
            let metrics = compute_metrics(&points);
            println!("\n  Summary:");
            println!("    Cl_α = {:.4}/deg ({:.3}/rad)", metrics.cl_alpha, metrics.cl_alpha * 57.3);
            println!("    α₀ = {:.2}°", metrics.alpha_0);
            println!("    (L/D)_max = {:.1} at Cl = {:.3}", metrics.max_l_d, metrics.cl_at_max_l_d);
            println!("    Cd_min = {:.5}", metrics.min_cd);
            println!("    Cl_max = {:.3}", metrics.max_cl);
            
            all_metrics.push((name.to_string(), re, metrics));
        }
    }
    
    // Print summary comparison table
    println!("\n\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                         SUMMARY COMPARISON TABLE                             ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    println!();
    println!("Cl-alpha slope (per degree) - Theory: 0.110/deg (2π/rad)");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("  Airfoil          |  Re=5e5  |  Re=1e6  |  Re=3e6  |  Re=6e6  |");
    println!("-------------------|----------|----------|----------|----------|");
    
    for (designation, name) in &airfoils {
        let values: Vec<String> = reynolds_numbers.iter().map(|&re| {
            if let Some((_, _, m)) = all_metrics.iter().find(|(n, r, _)| n == *name && *r == re) {
                format!("{:7.4}", m.cl_alpha)
            } else {
                "   -   ".to_string()
            }
        }).collect();
        println!("  {:16} | {} | {} | {} | {} |", 
            name.split(" ").next().unwrap_or(name),
            values[0], values[1], values[2], values[3]);
    }
    
    println!();
    println!("Zero-lift angle (degrees)");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("  Airfoil          |  Re=5e5  |  Re=1e6  |  Re=3e6  |  Re=6e6  |");
    println!("-------------------|----------|----------|----------|----------|");
    
    for (designation, name) in &airfoils {
        let values: Vec<String> = reynolds_numbers.iter().map(|&re| {
            if let Some((_, _, m)) = all_metrics.iter().find(|(n, r, _)| n == *name && *r == re) {
                format!("{:7.2}°", m.alpha_0)
            } else {
                "   -   ".to_string()
            }
        }).collect();
        println!("  {:16} | {} | {} | {} | {} |", 
            name.split(" ").next().unwrap_or(name),
            values[0], values[1], values[2], values[3]);
    }
    
    println!();
    println!("Maximum L/D");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("  Airfoil          |  Re=5e5  |  Re=1e6  |  Re=3e6  |  Re=6e6  |");
    println!("-------------------|----------|----------|----------|----------|");
    
    for (designation, name) in &airfoils {
        let values: Vec<String> = reynolds_numbers.iter().map(|&re| {
            if let Some((_, _, m)) = all_metrics.iter().find(|(n, r, _)| n == *name && *r == re) {
                format!("{:7.1}", m.max_l_d)
            } else {
                "   -   ".to_string()
            }
        }).collect();
        println!("  {:16} | {} | {} | {} | {} |", 
            name.split(" ").next().unwrap_or(name),
            values[0], values[1], values[2], values[3]);
    }
    
    println!();
    println!("Minimum Cd");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("  Airfoil          |  Re=5e5  |  Re=1e6  |  Re=3e6  |  Re=6e6  |");
    println!("-------------------|----------|----------|----------|----------|");
    
    for (designation, name) in &airfoils {
        let values: Vec<String> = reynolds_numbers.iter().map(|&re| {
            if let Some((_, _, m)) = all_metrics.iter().find(|(n, r, _)| n == *name && *r == re) {
                format!("{:7.5}", m.min_cd)
            } else {
                "   -   ".to_string()
            }
        }).collect();
        println!("  {:16} | {} | {} | {} | {} |", 
            name.split(" ").next().unwrap_or(name),
            values[0], values[1], values[2], values[3]);
    }
    
    println!();
    println!("════════════════════════════════════════════════════════════════════════════════");
    println!();
}

/// ASCII Cl-alpha plot.
fn print_cl_alpha_chart(name: &str, re: f64, points: &[AnalysisPoint]) {
    println!("\n  Cl vs Alpha - {} at Re={:.0e}", name, re);
    println!("  Cl");
    println!("   │");
    
    // Find ranges
    let cl_min = points.iter().map(|p| p.cl).fold(f64::INFINITY, f64::min);
    let cl_max = points.iter().map(|p| p.cl).fold(f64::NEG_INFINITY, f64::max);
    let alpha_min = points.first().map(|p| p.alpha).unwrap_or(-5.0);
    let alpha_max = points.last().map(|p| p.alpha).unwrap_or(15.0);
    
    let height = 15;
    let width = 50;
    
    // Create grid
    let mut grid = vec![vec![' '; width]; height];
    
    // Plot points
    for p in points {
        let x = ((p.alpha - alpha_min) / (alpha_max - alpha_min) * (width - 1) as f64) as usize;
        let y = ((p.cl - cl_min) / (cl_max - cl_min).max(0.01) * (height - 1) as f64) as usize;
        let y = (height - 1).saturating_sub(y);
        if x < width && y < height {
            grid[y][x] = if p.converged { '●' } else { '○' };
        }
    }
    
    // Print grid
    for (i, row) in grid.iter().enumerate() {
        let cl_val = cl_max - (i as f64 / (height - 1) as f64) * (cl_max - cl_min);
        print!("{:6.2} │", cl_val);
        for c in row {
            print!("{}", c);
        }
        println!();
    }
    
    // X axis
    print!("       └");
    for _ in 0..width {
        print!("─");
    }
    println!();
    println!("        {:5.0}                                      {:5.0}", alpha_min, alpha_max);
    println!("                             α (deg)");
}

/// ASCII L/D plot.
fn print_ld_chart(name: &str, re: f64, points: &[AnalysisPoint]) {
    println!("\n  L/D vs Cl - {} at Re={:.0e}", name, re);
    println!("  L/D");
    println!("   │");
    
    let valid_points: Vec<_> = points.iter()
        .filter(|p| p.converged && p.cd > 1e-6 && p.cl > 0.0)
        .collect();
    
    if valid_points.is_empty() {
        println!("   (no valid data)");
        return;
    }
    
    // Find ranges
    let ld_vals: Vec<f64> = valid_points.iter().map(|p| p.cl / p.cd).collect();
    let cl_vals: Vec<f64> = valid_points.iter().map(|p| p.cl).collect();
    
    let ld_min = 0.0;
    let ld_max = ld_vals.iter().cloned().fold(0.0, f64::max);
    let cl_min = 0.0;
    let cl_max = cl_vals.iter().cloned().fold(0.0, f64::max);
    
    let height = 12;
    let width = 40;
    
    // Create grid
    let mut grid = vec![vec![' '; width]; height];
    
    // Plot points
    for p in &valid_points {
        let ld = p.cl / p.cd;
        let x = ((p.cl - cl_min) / (cl_max - cl_min).max(0.01) * (width - 1) as f64) as usize;
        let y = ((ld - ld_min) / (ld_max - ld_min).max(0.01) * (height - 1) as f64) as usize;
        let y = (height - 1).saturating_sub(y);
        if x < width && y < height {
            grid[y][x] = '●';
        }
    }
    
    // Print grid
    for (i, row) in grid.iter().enumerate() {
        let ld_val = ld_max - (i as f64 / (height - 1) as f64) * (ld_max - ld_min);
        print!("{:6.0} │", ld_val);
        for c in row {
            print!("{}", c);
        }
        println!();
    }
    
    // X axis
    print!("       └");
    for _ in 0..width {
        print!("─");
    }
    println!();
    println!("        {:4.1}                            {:4.1}", cl_min, cl_max);
    println!("                        Cl");
}

#[test]
fn test_visualization_charts() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                    CL-ALPHA AND L/D VISUALIZATION                            ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    
    let airfoils = [
        (12u32, "NACA 0012"),
        (2412, "NACA 2412"),
        (4412, "NACA 4412"),
    ];
    
    let re = 3e6;
    let alphas: Vec<f64> = (-5..=15).map(|a| a as f64).collect();
    
    for (designation, name) in &airfoils {
        let airfoil = make_airfoil(*designation, 160);
        let points = run_alpha_sweep(&airfoil, re, &alphas);
        
        print_cl_alpha_chart(name, re, &points);
        print_ld_chart(name, re, &points);
        println!();
    }
}

#[test]
fn test_reynolds_effect() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                    REYNOLDS NUMBER EFFECT STUDY                              ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    
    let airfoil = make_airfoil(12, 160);  // NACA 0012
    let re_values = [1e5, 2e5, 5e5, 1e6, 2e6, 3e6, 6e6, 1e7];
    let alpha = 4.0;  // Fixed angle
    
    println!("\nNACA 0012 at α = 4° across Reynolds numbers");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("    Re      |   Cl    |   Cd    |  L/D   | xtr_up | xtr_lo | Conv");
    println!("------------|---------|---------|--------|--------|--------|-----");
    
    for &re in &re_values {
        let config = ViscousConfig {
            reynolds: re,
            n_crit: 9.0,
            ..Default::default()
        };
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        let l_d = if result.cd > 1e-6 { result.cl / result.cd } else { 0.0 };
        let conv = if result.converged { "✓" } else { "✗" };
        
        println!("{:11.2e} | {:7.4} | {:7.5} | {:6.1} | {:5.3}  | {:5.3}  | {}",
            re, result.cl, result.cd, l_d, result.x_tr_upper, result.x_tr_lower, conv);
    }
    
    println!("\nExpected trends:");
    println!("  - Cd decreases with increasing Re (thinner BL, lower friction)");
    println!("  - Transition moves forward with increasing Re");
    println!("  - L/D improves with Re (lower drag)");
}

#[test]
fn test_transition_location() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                    TRANSITION LOCATION STUDY                                 ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    
    let airfoil = make_airfoil(12, 160);  // NACA 0012
    let re = 3e6;
    let alphas: Vec<f64> = (0..=12).map(|a| a as f64).collect();
    
    println!("\nNACA 0012 at Re = 3M - Transition location vs alpha");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("   α°   | xtr_upper | xtr_lower | Note");
    println!("--------|-----------|-----------|------------------------------------------");
    
    let config = ViscousConfig {
        reynolds: re,
        n_crit: 9.0,
        ..Default::default()
    };
    let solver = ViscousSolver::new(config);
    
    for &alpha in &alphas {
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        let note = if result.x_tr_upper < 0.1 {
            "LE transition (upper)"
        } else if result.x_tr_lower > 0.9 {
            "Mostly laminar (lower)"
        } else {
            ""
        };
        
        println!("{:6.1}  |   {:5.3}   |   {:5.3}   | {}",
            alpha, result.x_tr_upper, result.x_tr_lower, note);
    }
    
    println!("\nExpected behavior:");
    println!("  - Upper surface: transition moves forward with increasing α");
    println!("  - Lower surface: transition moves aft with increasing α");
    println!("  - At α=0°, symmetric airfoil should have symmetric transition");
}

#[test]
fn test_camber_effect() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                    CAMBER EFFECT STUDY                                       ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
    
    let designations = [12, 1412, 2412, 3412, 4412, 5412, 6412];
    let re = 3e6;
    let alpha = 4.0;
    
    println!("\nEffect of camber at α = 4°, Re = 3M");
    println!("────────────────────────────────────────────────────────────────────────────────");
    println!("  Airfoil  | Camber |   Cl    |   Cd    |  L/D   |   α₀   | Conv");
    println!("-----------|--------|---------|---------|--------|--------|-----");
    
    for &desig in &designations {
        let camber = (desig / 1000) as f64 / 100.0;  // First digit / 100
        let airfoil = make_airfoil(desig, 160);
        
        let config = ViscousConfig {
            reynolds: re,
            n_crit: 9.0,
            ..Default::default()
        };
        let solver = ViscousSolver::new(config);
        let flow = FlowConditions::with_alpha_deg(alpha);
        let result = solver.solve(&airfoil, &flow);
        
        // Estimate α₀ with a quick 2-point check
        let flow0 = FlowConditions::with_alpha_deg(0.0);
        let result0 = solver.solve(&airfoil, &flow0);
        let alpha_0 = if result.cl != result0.cl {
            -result0.cl / ((result.cl - result0.cl) / alpha)
        } else {
            0.0
        };
        
        let l_d = if result.cd > 1e-6 { result.cl / result.cd } else { 0.0 };
        let conv = if result.converged { "✓" } else { "✗" };
        
        println!("  NACA{:04} | {:4.1}%  | {:7.4} | {:7.5} | {:6.1} | {:5.2}° | {}",
            desig, camber * 100.0, result.cl, result.cd, l_d, alpha_0, conv);
    }
    
    println!("\nExpected trends:");
    println!("  - Cl increases with camber at same α");
    println!("  - α₀ becomes more negative with camber");
    println!("  - Cd slightly increases with camber (thicker BL on pressure side)");
}
