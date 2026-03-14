//! Newton iteration debugging test

use rustfoil_bl::equations::{bldif, blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::solve::{build_4x4_system, solve_4x4};
use serde::Deserialize;
use std::fs;

fn load_reference() -> Option<serde_json::Value> {
    let paths = [
        "testdata/mrchue_iterations.json",
        "../testdata/mrchue_iterations.json",
        "../../testdata/mrchue_iterations.json",
    ];

    for path in &paths {
        if let Ok(content) = fs::read_to_string(path) {
            if let Ok(data) = serde_json::from_str(&content) {
                return Some(data);
            }
        }
    }
    None
}

#[test]
fn test_newton_iteration_detail() {
    let ref_data = match load_reference() {
        Some(d) => d,
        None => {
            eprintln!("Skipping: mrchue_iterations.json not found");
            return;
        }
    };

    let re = ref_data["metadata"]["reynolds"].as_f64().unwrap();
    let msq = 0.0;

    let stations = &ref_data["sides"]["1"]["stations"];
    let s0 = &stations[0]; // IBL=2 (first station)
    let s1 = &stations[1]; // IBL=3 (second station)

    // Build "prev" station with XFOIL's converged values
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    // Build "curr" station with initial guess (prev station's values)
    let mut curr = BlStation::new();
    curr.x = s1["x"].as_f64().unwrap();
    curr.u = s1["Ue"].as_f64().unwrap();
    curr.theta = prev.theta; // Initial guess
    curr.delta_star = prev.delta_star;
    curr.ampl = 0.0;
    curr.ctau = 0.03;
    curr.is_laminar = true;
    blvar(&mut curr, FlowType::Laminar, msq, re);

    println!("\n=== Newton Iteration at Station 2 (IBL=3) ===");
    println!("Previous station: x={:.6}, θ={:.6e}", prev.x, prev.theta);
    println!("Current station:  x={:.6}, Ue={:.6}", curr.x, curr.u);
    println!("\nXFOIL target: θ={:.6e}", s1["final"]["theta"].as_f64().unwrap());

    // Manual Newton iteration
    for iter in 0..5 {
        println!("\n--- Iteration {} ---", iter);
        println!("Current θ={:.6e}, δ*={:.6e}, H={:.4}, Hk={:.4}", 
            curr.theta, curr.delta_star, curr.h, curr.hk);

        // Compute residuals and Jacobian
        let (res, jac) = bldif(&prev, &curr, FlowType::Laminar, msq, re);
        
        println!("Residuals: res_ampl={:.4e}, res_mom={:.4e}, res_shape={:.4e}",
            res.res_third, res.res_mom, res.res_shape);

        // Build 4x4 system
        let hlmax = 4.0;
        let (a, b) = build_4x4_system(
            &jac.vs2,
            &[res.res_third, res.res_mom, res.res_shape],
            true, // direct mode
            curr.hk / curr.theta,  // hk2_t
            -curr.hk / curr.delta_star, // hk2_d
            0.0, // hk2_u
            hlmax, // htarg
            curr.hk, // hk_current
        );

        // Print system matrix
        println!("\n4x4 System:");
        println!("A[0] = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", a[0][0], a[0][1], a[0][2], a[0][3]);
        println!("A[1] = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", a[1][0], a[1][1], a[1][2], a[1][3]);
        println!("A[2] = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", a[2][0], a[2][1], a[2][2], a[2][3]);
        println!("A[3] = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", a[3][0], a[3][1], a[3][2], a[3][3]);
        println!("b = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", b[0], b[1], b[2], b[3]);

        // Solve
        let dx = solve_4x4(&a, &b);
        println!("\nSolution dx = [{:+.4e}, {:+.4e}, {:+.4e}, {:+.4e}]", dx[0], dx[1], dx[2], dx[3]);

        // Compute dmax
        let dmax = (dx[1] / curr.theta.max(1e-12)).abs()
            .max((dx[2] / curr.delta_star.max(1e-12)).abs())
            .max((dx[0] / 10.0).abs());
        println!("dmax = {:.4e}", dmax);

        // Relaxation
        let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };
        println!("rlx = {:.4}", rlx);

        // Apply updates
        let theta_new = curr.theta + rlx * dx[1];
        let dstar_new = curr.delta_star + rlx * dx[2];
        
        println!("\nUpdate: θ += {:.4e}, δ* += {:.4e}", rlx * dx[1], rlx * dx[2]);
        println!("New θ={:.6e}, δ*={:.6e}", theta_new.max(1e-12), dstar_new.max(1e-12));

        // Check convergence
        if dmax <= 1e-5 {
            println!("\nConverged!");
            break;
        }

        // Update station for next iteration
        curr.theta = theta_new.max(1e-12);
        curr.delta_star = dstar_new.max(1e-12);
        blvar(&mut curr, FlowType::Laminar, msq, re);
    }

    println!("\n=== Final comparison ===");
    println!("RustFoil θ = {:.6e}", curr.theta);
    println!("XFOIL θ    = {:.6e}", s1["final"]["theta"].as_f64().unwrap());
}
