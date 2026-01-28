//! Check station to DIJ mapping

use rustfoil_solver::viscous::{
    setup_from_coords, ViscousSolverConfig,
    extract_surface_xfoil, compute_arc_lengths,
};
use std::fs;

#[test]
fn test_station_to_panel_mapping() {
    let content = fs::read_to_string("../../Xfoil-instrumented/bin/naca0012_xfoil.dat")
        .expect("Failed to read coordinate file");
    let coords: Vec<(f64, f64)> = content
        .lines()
        .skip(1)
        .filter_map(|line| {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                Some((parts[0].parse().ok()?, parts[1].parse().ok()?))
            } else { None }
        })
        .collect();
    
    let result = setup_from_coords(&coords, 2.0).expect("Setup failed");
    let dij = &result.setup.dij;
    
    let full_arc = compute_arc_lengths(&result.node_x, &result.node_y);
    let ist = result.ist;
    let sst = result.sst;
    let ue_stag = result.setup.ue_inviscid[ist].abs().max(0.01);
    
    let (upper_arc, upper_x, _, upper_ue) = extract_surface_xfoil(
        ist, sst, ue_stag, &full_arc, &result.node_x, &result.node_y, 
        &result.setup.ue_inviscid, true
    );
    
    let config = ViscousSolverConfig::with_reynolds(1e6);
    // Use initialize_surface_stations_with_panel_idx to properly set panel indices
    use rustfoil_solver::viscous::initialize_surface_stations_with_panel_idx;
    let upper_stations = initialize_surface_stations_with_panel_idx(
        &upper_arc, &upper_ue, &upper_x, ist, true, config.reynolds);
    
    println!("\n=== Station to Panel Mapping (Upper Surface) ===");
    println!("Stagnation index (ist): {}", ist);
    println!("DIJ matrix size: {}x{}", dij.nrows(), dij.ncols());
    println!("Number of upper stations: {}", upper_stations.len());
    
    println!("\nFirst 10 stations:");
    for (i, s) in upper_stations.iter().take(10).enumerate() {
        println!("  Station {}: panel_idx={}, x/c={:.4}, Ue={:.4}", 
            i, s.panel_idx, s.x_coord, s.u);
    }
    
    println!("\nLast 5 stations:");
    for (i, s) in upper_stations.iter().rev().take(5).rev().enumerate() {
        let idx = upper_stations.len() - 5 + i;
        println!("  Station {}: panel_idx={}, x/c={:.4}, Ue={:.4}", 
            idx, s.panel_idx, s.x_coord, s.u);
    }
    
    // Check if panel_idx values are valid DIJ indices
    let valid_indices: Vec<bool> = upper_stations.iter()
        .map(|s| s.panel_idx < dij.nrows())
        .collect();
    let all_valid = valid_indices.iter().all(|&v| v);
    println!("\nAll panel_idx values valid for DIJ? {}", all_valid);
    
    if !all_valid {
        for (i, &valid) in valid_indices.iter().enumerate() {
            if !valid {
                println!("  Invalid: station {} has panel_idx={}, but DIJ has {} rows",
                    i, upper_stations[i].panel_idx, dij.nrows());
            }
        }
    }
}
