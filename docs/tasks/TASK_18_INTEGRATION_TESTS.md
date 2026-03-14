# Task 18: Integration Tests and XFOIL Comparison

## Objective
Create end-to-end tests comparing RustFoil with XFOIL.

## Prerequisites
- All previous tasks completed (01-17)

## Deliverables

### Create rustfoil-solver/tests/xfoil_comparison.rs
```rust
//! End-to-end comparison with XFOIL binary
//!
//! These tests run both RustFoil and XFOIL on the same cases
//! and verify results match within acceptable tolerance.

use rustfoil_core::Body;
use rustfoil_solver::viscous::{solve_viscous, ViscousSolverConfig};
use std::process::Command;
use std::io::Write;
use tempfile::NamedTempFile;

/// Helper to run XFOIL and extract results
struct XfoilRunner {
    binary_path: String,
}

impl XfoilRunner {
    fn new(path: &str) -> Self {
        Self { binary_path: path.to_string() }
    }
    
    fn analyze(&self, airfoil: &str, alpha: f64, re: f64) -> XfoilResult {
        // Create XFOIL command file
        let commands = format!(
            "load {}\npane\noper\nvisc {}\nalfa {}\n\nquit\n",
            airfoil, re, alpha
        );
        
        let mut cmd_file = NamedTempFile::new().unwrap();
        write!(cmd_file, "{}", commands).unwrap();
        
        let output = Command::new(&self.binary_path)
            .stdin(std::fs::File::open(cmd_file.path()).unwrap())
            .output()
            .expect("Failed to run XFOIL");
        
        // Parse output for CL, CD, CM
        let stdout = String::from_utf8_lossy(&output.stdout);
        parse_xfoil_output(&stdout)
    }
}

fn parse_xfoil_output(output: &str) -> XfoilResult {
    // Parse XFOIL's output format
    // Looking for lines like:
    //   CL =   0.5432   CM = -0.0123   CD =  0.00654
    
    let mut result = XfoilResult::default();
    
    for line in output.lines() {
        if line.contains("CL =") {
            // Parse CL, CM, CD from this line
            // TODO: Implement proper parsing
        }
    }
    
    result
}

#[derive(Default)]
struct XfoilResult {
    cl: f64,
    cd: f64,
    cm: f64,
    x_tr_upper: f64,
    x_tr_lower: f64,
}

#[test]
#[ignore] // Run manually: cargo test xfoil -- --ignored
fn test_naca0012_vs_xfoil() {
    let xfoil = XfoilRunner::new("../Xfoil/bin/xfoil");
    
    // Test cases: NACA 0012 at Re=3e6
    let test_alphas = [-4.0, 0.0, 4.0, 8.0, 12.0];
    let re = 3e6;
    
    let body = Body::from_naca4(12, 160);
    let config = ViscousSolverConfig::with_reynolds(re);
    
    for alpha in test_alphas {
        println!("\n=== Testing alpha = {}° ===", alpha);
        
        // Run RustFoil
        let rust_result = solve_viscous(&body, alpha, &config)
            .expect("RustFoil failed");
        
        // Run XFOIL
        let xfoil_result = xfoil.analyze("naca0012.dat", alpha, re);
        
        // Compare
        let cl_diff = (rust_result.cl - xfoil_result.cl).abs();
        let cd_diff = (rust_result.cd - xfoil_result.cd).abs();
        
        println!("CL: Rust={:.4}, XFOIL={:.4}, diff={:.4}",
            rust_result.cl, xfoil_result.cl, cl_diff);
        println!("CD: Rust={:.5}, XFOIL={:.5}, diff={:.5}",
            rust_result.cd, xfoil_result.cd, cd_diff);
        
        // Tolerances
        assert!(cl_diff < 0.01, "CL difference too large: {}", cl_diff);
        assert!(cd_diff < 0.001, "CD difference too large: {}", cd_diff);
    }
}

#[test]
fn test_blasius_flat_plate() {
    // Analytical test: flat plate boundary layer (Blasius solution)
    // H = 2.59 for laminar BL
    
    // This doesn't need XFOIL, just validates BL equations
    let x = vec![0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
    let ue = vec![1.0; x.len()]; // Constant edge velocity
    let re = 1e6;
    
    let config = rustfoil_coupling::march::MarchConfig::default();
    let result = rustfoil_coupling::march::march_fixed_ue(
        &x, &ue, re, 0.0, &config
    );
    
    // Check shape factor approaches Blasius value
    for station in &result.stations {
        if station.is_laminar {
            // Blasius: H ≈ 2.59
            assert!(station.h > 2.3 && station.h < 2.9,
                "H={} outside expected range", station.h);
        }
    }
}

#[test]
fn test_transition_location() {
    // Test that transition occurs at reasonable location
    let body = Body::from_naca4(12, 160);
    let config = ViscousSolverConfig {
        reynolds: 3e6,
        ncrit: 9.0,
        ..Default::default()
    };
    
    let result = solve_viscous(&body, 0.0, &config).unwrap();
    
    // At Re=3e6, α=0°, transition should be roughly mid-chord
    assert!(result.x_tr_upper > 0.1 && result.x_tr_upper < 0.8,
        "Transition at unexpected location: {}", result.x_tr_upper);
}

#[test]
fn test_separation_detection() {
    // At high angle of attack, should detect separation
    let body = Body::from_naca4(12, 160);
    let config = ViscousSolverConfig::with_reynolds(1e6);
    
    // At 16°, NACA 0012 should be stalled
    let result = solve_viscous(&body, 16.0, &config);
    
    // Solution may not converge due to separation
    // That's expected behavior at post-stall conditions
    match result {
        Ok(r) => {
            // If it converged, CL should be reduced from inviscid
            println!("Converged at stall: CL={}", r.cl);
        }
        Err(_) => {
            println!("Non-convergence at stall (expected)");
        }
    }
}
```

### Create testdata for reference comparisons
Generate reference data by running XFOIL:
```bash
cd /Users/harry/flexfoil-boundary-layer/testdata
./generate_xfoil_reference.sh
```

## Verification
```bash
# Run unit tests
cargo test -p rustfoil-bl
cargo test -p rustfoil-coupling

# Run integration tests (excluding XFOIL comparison)
cargo test -p rustfoil-solver

# Run XFOIL comparison (requires XFOIL binary)
cargo test -p rustfoil-solver xfoil -- --ignored
```

## Final Checklist

- [ ] All closure functions match FORTRAN to 1e-10
- [ ] BL march produces Blasius solution for flat plate
- [ ] Transition prediction gives reasonable x_tr
- [ ] CL matches XFOIL within 0.01
- [ ] CD matches XFOIL within 0.001
- [ ] CLI commands work correctly
- [ ] Parallel polar generation works

## Project Complete!

After this task, the viscous solver implementation is complete. 
Run comprehensive validation and document any limitations.

---

## Documentation Requirements

Also ensure that you update Docusaurus with progress.

Explain what tests were for, what they show, and how they passed/failed/worked and consequences.
