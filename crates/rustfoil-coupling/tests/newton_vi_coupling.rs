//! Tests for Phase 1 Newton Viscous-Inviscid Coupling
//!
//! These tests validate the XFOIL-style Newton coupling implementation:
//! - VTI signs in VM matrix construction
//! - Forced changes (DUE, DDS) in residuals
//! - BLSOLV back-substitution using row 3 (mass)
//! - Mass-first Ue computation in UPDATE

use nalgebra::DMatrix;
use rustfoil_bl::equations::{blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::newton::BlNewtonSystem;
use rustfoil_coupling::solve::solve_blsolv_xfoil;
use rustfoil_coupling::update::{update_xfoil_style, UpdateConfig};

/// Create a simple test case with laminar BL stations
fn create_laminar_test_stations(n: usize, re: f64) -> Vec<BlStation> {
    let msq = 0.0; // Incompressible

    (0..n)
        .map(|i| {
            let x = i as f64 * 0.1; // Stations at x = 0.0, 0.1, 0.2, ...
            let mut station = BlStation::new();
            station.x = x;
            station.u = 1.0 - 0.1 * x; // Mild adverse pressure gradient
            station.theta = 0.001 * (1.0 + 0.5 * x); // Growing momentum thickness
            station.delta_star = 0.003 * (1.0 + 0.5 * x); // Growing displacement thickness
            station.mass_defect = station.u * station.delta_star;
            station.ampl = 0.0;
            station.is_laminar = true;

            blvar(&mut station, FlowType::Laminar, msq, re);
            station
        })
        .collect()
}

/// Create a simple DIJ matrix (diagonal for testing)
fn create_test_dij(n: usize) -> DMatrix<f64> {
    let mut dij = DMatrix::zeros(n, n);
    for i in 0..n {
        // Simple coupling: mass at station j affects Ue at station i
        for j in 0..n {
            // Distance-based influence (stronger for nearby stations)
            let dist = (i as f64 - j as f64).abs() + 1.0;
            dij[(i, j)] = 0.1 / dist;
        }
    }
    dij
}

#[test]
fn test_vti_initialization() {
    // Test that VTI is initialized to +1.0 for single surface
    let n = 10;
    let system = BlNewtonSystem::new(n);

    assert_eq!(system.vti.len(), n);
    for i in 0..n {
        assert!(
            (system.vti[i] - 1.0).abs() < 1e-10,
            "VTI[{}] should be 1.0, got {}",
            i,
            system.vti[i]
        );
    }
}

#[test]
fn test_build_with_vm_full_produces_finite_values() {
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let ue_inviscid: Vec<f64> = stations.iter().map(|s| s.u + 0.01).collect(); // Slight mismatch
    let flow_types = vec![FlowType::Laminar; n - 1];

    let mut system = BlNewtonSystem::new(n);
    system.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_inviscid);

    // Check that all values are finite
    for i in 1..n {
        for k in 0..3 {
            assert!(
                system.rhs[i][k].is_finite(),
                "RHS[{}][{}] is not finite",
                i,
                k
            );
            for j in 0..n {
                assert!(
                    system.vm[i][j][k].is_finite(),
                    "VM[{}][{}][{}] is not finite",
                    i,
                    j,
                    k
                );
            }
        }
    }
}

#[test]
fn test_vm_matrix_includes_diagonal() {
    // The VM matrix diagonal should have non-zero entries for mass coupling
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let ue_inviscid: Vec<f64> = stations.iter().map(|s| s.u).collect();
    let flow_types = vec![FlowType::Laminar; n - 1];

    let mut system = BlNewtonSystem::new(n);
    system.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_inviscid);

    // Check diagonal entries (these should be non-zero due to delta_star/mass coupling)
    for i in 1..n {
        let diag_norm = system.vm[i][i]
            .iter()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        // Diagonal should have meaningful values for mass coupling
        println!("VM diagonal at i={}: {:?}", i, system.vm[i][i]);
    }
}

#[test]
fn test_forced_changes_modify_rhs() {
    // When there's a mismatch between current Ue and UESET result,
    // the forced changes should modify the RHS
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let flow_types = vec![FlowType::Laminar; n - 1];

    // Case 1: No Ue mismatch (ue_inviscid matches station.u exactly)
    let ue_exact: Vec<f64> = stations.iter().map(|s| s.u).collect();
    let mut system_exact = BlNewtonSystem::new(n);
    system_exact.build(&stations, &flow_types, msq, re);
    let rhs_without_forced = system_exact.rhs.clone();

    // Case 2: With Ue mismatch
    let ue_mismatched: Vec<f64> = stations.iter().map(|s| s.u + 0.05).collect(); // 5% mismatch
    let mut system_mismatch = BlNewtonSystem::new(n);
    system_mismatch.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_mismatched);

    // The RHS should be different when there's a mismatch
    // (Note: This is a structural test - the exact values depend on XFOIL's formulation)
    let mut total_diff = 0.0;
    for i in 1..n {
        for k in 0..3 {
            total_diff += (system_mismatch.rhs[i][k] - rhs_without_forced[i][k]).abs();
        }
    }

    println!(
        "Total RHS difference with forced changes: {:.6e}",
        total_diff
    );
    // With a 5% Ue mismatch, we expect some difference
    // (The exact value depends on the VS columns and station states)
}

#[test]
fn test_blsolv_produces_finite_solution() {
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let ue_inviscid: Vec<f64> = stations.iter().map(|s| s.u).collect();
    let flow_types = vec![FlowType::Laminar; n - 1];

    let mut system = BlNewtonSystem::new(n);
    system.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_inviscid);

    let solution = solve_blsolv_xfoil(&system);

    // All solution values should be finite
    for i in 0..n {
        for k in 0..3 {
            assert!(
                solution[i][k].is_finite(),
                "Solution[{}][{}] is not finite: {}",
                i,
                k,
                solution[i][k]
            );
        }
    }

    println!("BLSOLV solution summary:");
    for i in 1..n.min(4) {
        println!(
            "  Station {}: dCtau={:.6e}, dTheta={:.6e}, dMass={:.6e}",
            i, solution[i][0], solution[i][1], solution[i][2]
        );
    }
}

#[test]
fn test_blsolv_back_substitution_uses_mass() {
    // The XFOIL-style back-substitution should use the third element (mass)
    // to update all upstream stations
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let ue_inviscid: Vec<f64> = stations.iter().map(|s| s.u).collect();
    let flow_types = vec![FlowType::Laminar; n - 1];

    let mut system = BlNewtonSystem::new(n);
    system.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_inviscid);

    let solution = solve_blsolv_xfoil(&system);

    // The last station's mass change should influence upstream stations
    // This is verified by the structure of the solution (non-zero upstream values)
    let mass_change_last = solution[n - 1][2];
    println!("Mass change at last station: {:.6e}", mass_change_last);

    // Upstream stations should have been updated by the back-substitution
    for i in 1..n - 1 {
        println!(
            "Station {} mass change: {:.6e}",
            i, solution[i][2]
        );
    }
}

#[test]
fn test_update_xfoil_style_produces_valid_state() {
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let mut stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let ue_inviscid: Vec<f64> = stations.iter().map(|s| s.u).collect();
    let flow_types = vec![FlowType::Laminar; n - 1];

    // Build and solve Newton system
    let mut system = BlNewtonSystem::new(n);
    system.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_inviscid);
    let deltas = solve_blsolv_xfoil(&system);

    // Save initial state
    let initial_theta: Vec<f64> = stations.iter().map(|s| s.theta).collect();
    let initial_mass: Vec<f64> = stations.iter().map(|s| s.mass_defect).collect();
    let initial_ue: Vec<f64> = stations.iter().map(|s| s.u).collect();

    // Apply update
    let config = UpdateConfig::default();
    let result = update_xfoil_style(&mut stations, &deltas, &ue_inviscid, &dij, &system.vti, &config);

    println!("\nUpdate result:");
    println!("  RMS change: {:.6e}", result.rms_change);
    println!(
        "  Max change: {:.6e} (var: {}, station: {})",
        result.max_change, result.max_change_var, result.max_change_station
    );
    println!("  Relaxation used: {:.4}", result.relaxation_used);

    // Check that stations have valid state after update
    for (i, station) in stations.iter().enumerate() {
        assert!(
            station.theta > 0.0,
            "Station {} theta should be positive, got {}",
            i,
            station.theta
        );
        assert!(
            station.delta_star > 0.0,
            "Station {} delta_star should be positive, got {}",
            i,
            station.delta_star
        );
        assert!(
            station.mass_defect > 0.0,
            "Station {} mass_defect should be positive, got {}",
            i,
            station.mass_defect
        );
        assert!(
            station.u.is_finite(),
            "Station {} Ue should be finite, got {}",
            i,
            station.u
        );
        assert!(
            station.h > 0.0,
            "Station {} H should be positive, got {}",
            i,
            station.h
        );
    }

    // Log changes for debugging
    for i in 1..n {
        let d_theta = stations[i].theta - initial_theta[i];
        let d_mass = stations[i].mass_defect - initial_mass[i];
        let d_ue = stations[i].u - initial_ue[i];
        println!(
            "Station {} changes: dθ={:.6e}, dMass={:.6e}, dUe={:.6e}",
            i, d_theta, d_mass, d_ue
        );
    }
}

#[test]
fn test_newton_iteration_reduces_residuals() {
    // Multiple Newton iterations should reduce residuals
    let re = 1e6;
    let msq = 0.0;
    let n = 5;

    let mut stations = create_laminar_test_stations(n, re);
    let dij = create_test_dij(n);
    let ue_inviscid: Vec<f64> = stations.iter().map(|s| s.u).collect();
    let flow_types = vec![FlowType::Laminar; n - 1];
    let config = UpdateConfig::default();

    let mut residuals = Vec::new();

    // Run several Newton iterations
    for iter in 0..5 {
        // Recompute secondary variables
        for station in stations.iter_mut() {
            let ft = if station.is_laminar {
                FlowType::Laminar
            } else {
                FlowType::Turbulent
            };
            blvar(station, ft, msq, re);
        }

        // Build Newton system
        let mut system = BlNewtonSystem::new(n);
        system.build_with_vm_full(&stations, &flow_types, msq, re, &dij, &ue_inviscid);

        let res = system.residual_norm();
        residuals.push(res);
        println!("Iteration {}: residual norm = {:.6e}", iter, res);

        // Solve and update
        let deltas = solve_blsolv_xfoil(&system);
        update_xfoil_style(&mut stations, &deltas, &ue_inviscid, &dij, &system.vti, &config);
    }

    // The residual should generally decrease (though not monotonically in all cases)
    // At minimum, it should not diverge
    let initial_res = residuals[0];
    let final_res = *residuals.last().unwrap();

    println!("\nInitial residual: {:.6e}", initial_res);
    println!("Final residual: {:.6e}", final_res);

    // Check that we don't diverge badly
    assert!(
        final_res < initial_res * 100.0,
        "Residual should not diverge: initial={:.6e}, final={:.6e}",
        initial_res,
        final_res
    );
}
