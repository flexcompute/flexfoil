//! Test march starting from XFOIL's converged first station values

use rustfoil_bl::equations::{blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::march::newton_solve_station;
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
fn test_march_from_xfoil_initial() {
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
    
    println!("\n=== March from XFOIL Initial Values ===\n");
    println!(" Idx          x      θ_xfoil       θ_rust   θ_err%");

    // Start from XFOIL's converged station 0 (IBL=2)
    let s0 = &stations[0];
    let mut prev = BlStation::new();
    prev.x = s0["x"].as_f64().unwrap();
    prev.u = s0["Ue"].as_f64().unwrap();
    prev.theta = s0["final"]["theta"].as_f64().unwrap();  // Use XFOIL's converged value!
    prev.delta_star = s0["final"]["delta_star"].as_f64().unwrap();
    prev.ampl = 0.0;
    prev.ctau = 0.03;
    prev.is_laminar = true;
    blvar(&mut prev, FlowType::Laminar, msq, re);

    println!("   0    {:7.4}  {:.6e}  {:.6e}    {:+5.1}",
        prev.x, prev.theta, prev.theta, 0.0);

    // March through subsequent stations
    for i in 1..10.min(stations.as_array().unwrap().len()) {
        let s = &stations[i];
        let x_new = s["x"].as_f64().unwrap();
        let ue_new = s["Ue"].as_f64().unwrap();
        let theta_xfoil = s["final"]["theta"].as_f64().unwrap();

        // Solve using Newton
        let (station, converged, _dmax) = newton_solve_station(
            &prev,
            x_new,
            ue_new,
            re,
            msq,
            true, // is_laminar
            25, // max_iter
            1e-5, // tolerance
            3.8, // hlmax
            2.5, // htmax
        );

        let err = (station.theta - theta_xfoil) / theta_xfoil * 100.0;
        let flag = if err.abs() > 10.0 { " ***" } else { "" };
        println!("  {:2}    {:7.4}  {:.6e}  {:.6e}  {:+6.1}{}",
            i, x_new, theta_xfoil, station.theta, err, flag);

        // Update prev for next iteration
        prev = station;
    }
}
