//! Check if viscous solver updates Ue values from inviscid

use rustfoil_solver::viscous::{
    setup_from_coords, ViscousSolverConfig, 
    extract_surface_xfoil, compute_arc_lengths,
    solve_viscous_two_surfaces,
};
use std::fs;

#[test]
fn test_ue_change_from_viscous() {
    // Load coordinates
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
    
    let alpha_deg = 2.0;
    println!("\n=== Testing Ue change at α = {:.1}° ===", alpha_deg);
    
    let result = setup_from_coords(&coords, alpha_deg).expect("Setup failed");
    
    let full_arc = compute_arc_lengths(&result.node_x, &result.node_y);
    let ist = result.ist;
    let sst = result.sst;
    let ue_stag = result.setup.ue_inviscid[ist].abs().max(0.01);
    
    let (upper_arc, upper_x, _, upper_ue_inviscid) = extract_surface_xfoil(
        ist, sst, ue_stag, &full_arc, &result.node_x, &result.node_y, 
        &result.setup.ue_inviscid, true
    );
    let (lower_arc, lower_x, _, lower_ue_inviscid) = extract_surface_xfoil(
        ist, sst, ue_stag, &full_arc, &result.node_x, &result.node_y,
        &result.setup.ue_inviscid, false
    );
    
    let config = ViscousSolverConfig::with_reynolds(1e6).with_max_iterations(20);
    
    // CRITICAL: Use initialize_surface_stations_with_panel_idx to properly map
    // BL stations to panel indices for VI coupling
    use rustfoil_solver::viscous::initialize_surface_stations_with_panel_idx;
    let mut upper_stations = initialize_surface_stations_with_panel_idx(
        &upper_arc, &upper_ue_inviscid, &upper_x, ist, true, config.reynolds);
    let mut lower_stations = initialize_surface_stations_with_panel_idx(
        &lower_arc, &lower_ue_inviscid, &lower_x, ist, false, config.reynolds);
    
    // Save inviscid Ue for comparison
    let upper_ue_before: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
    let lower_ue_before: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();
    
    // Run viscous solver
    let visc_result = solve_viscous_two_surfaces(
        &mut upper_stations,
        &mut lower_stations,
        &upper_ue_inviscid,
        &lower_ue_inviscid,
        &result.setup.dij,
        &config,
    ).expect("Viscous solve failed");
    
    // Get viscous Ue after solving
    let upper_ue_after: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
    let lower_ue_after: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();
    
    // Compute max relative change
    let mut max_upper_change = 0.0f64;
    for i in 0..upper_ue_before.len() {
        if upper_ue_before[i].abs() > 0.01 {
            let rel_change = ((upper_ue_after[i] - upper_ue_before[i]) / upper_ue_before[i]).abs();
            max_upper_change = max_upper_change.max(rel_change);
        }
    }
    
    let mut max_lower_change = 0.0f64;
    for i in 0..lower_ue_before.len() {
        if lower_ue_before[i].abs() > 0.01 {
            let rel_change = ((lower_ue_after[i] - lower_ue_before[i]) / lower_ue_before[i]).abs();
            max_lower_change = max_lower_change.max(rel_change);
        }
    }
    
    println!("\nInviscid CL from panel method: {:.4}", result.inviscid.cl);
    println!("Viscous CL from circulation: {:.4}", visc_result.cl);
    println!("CL difference: {:.4} ({:.1}%)", 
        visc_result.cl - result.inviscid.cl,
        (visc_result.cl - result.inviscid.cl) / result.inviscid.cl.abs() * 100.0);
    
    println!("\nMax upper surface Ue relative change: {:.4}%", max_upper_change * 100.0);
    println!("Max lower surface Ue relative change: {:.4}%", max_lower_change * 100.0);
    
    // Sample some Ue values
    println!("\nUpper surface Ue (inviscid -> viscous):");
    let sample_points = [0, upper_ue_before.len()/4, upper_ue_before.len()/2, upper_ue_before.len()-1];
    for &i in &sample_points {
        if i < upper_ue_before.len() {
            let x = upper_arc.get(i).unwrap_or(&0.0);
            println!("  i={}: x/c={:.3}, Ue: {:.4} -> {:.4} (change: {:.2}%)",
                i, x, upper_ue_before[i], upper_ue_after[i],
                (upper_ue_after[i] - upper_ue_before[i]) / upper_ue_before[i].abs().max(0.001) * 100.0);
        }
    }
}
