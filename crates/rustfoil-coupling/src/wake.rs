//! Wake marching and trailing edge combination.
//!
//! This module implements wake handling for viscous-inviscid coupling:
//! - Combining upper and lower surface BL at the trailing edge
//! - Marching the wake downstream using BLDIF(FlowType::Wake)
//!
//! # XFOIL Reference
//! - TESYS: xbl.f lines 236-276 (TE system setup)
//! - Wake marching: xbl.f lines 607-927 (same as MRCHUE but with wake flag)

use rustfoil_bl::equations::{bldif, blvar, FlowType};
use rustfoil_bl::state::BlStation;

/// Combine upper and lower trailing edge stations into wake initial condition.
///
/// At the trailing edge, the boundary layers from both surfaces merge into
/// a single wake. The momentum and displacement thicknesses add, while
/// the shear stress coefficient is a momentum-weighted average.
///
/// # Arguments
/// * `upper_te` - Upper surface trailing edge station
/// * `lower_te` - Lower surface trailing edge station
///
/// # Returns
/// Initial wake station with combined properties
///
/// # XFOIL Reference
/// xbl.f TESYS (lines 236-276):
/// ```fortran
/// TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
/// DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
/// CTE = (CTAU(IBLTE(1),1)*THET(IBLTE(1),1) + CTAU(IBLTE(2),2)*THET(IBLTE(2),2)) / TTE
/// ```
pub fn combine_te_for_wake(upper_te: &BlStation, lower_te: &BlStation) -> BlStation {
    let mut wake = BlStation::default();

    // Momentum thickness adds (conservation of momentum)
    wake.theta = upper_te.theta + lower_te.theta;

    // Displacement thickness adds (both surfaces contribute)
    // Note: ANTE (base thickness for blunt TE) is not included here
    wake.delta_star = upper_te.delta_star + lower_te.delta_star;

    // Shape factor from combined values
    wake.h = if wake.theta > 1e-12 {
        wake.delta_star / wake.theta
    } else {
        1.5 // Default wake H
    };

    // Shear stress coefficient is momentum-weighted average
    wake.ctau = if wake.theta > 1e-12 {
        (upper_te.ctau * upper_te.theta + lower_te.ctau * lower_te.theta) / wake.theta
    } else {
        0.5 * (upper_te.ctau + lower_te.ctau)
    };

    // Edge velocity is average (should be similar at TE)
    wake.u = 0.5 * (upper_te.u.abs() + lower_te.u.abs());

    // Mass defect
    wake.mass_defect = wake.u * wake.delta_star;

    // Wake flags
    wake.is_wake = true;
    wake.is_turbulent = true;
    wake.is_laminar = false;

    // Skin friction is zero in wake (no wall)
    wake.cf = 0.0;

    // Position at trailing edge (x/c = 1.0)
    wake.x = 1.0;
    wake.x_coord = 1.0;

    wake
}

/// Generate wake station positions downstream of trailing edge.
///
/// Creates arc lengths for wake stations extending from TE to downstream.
/// XFOIL uses wake panels from potential flow; we use a simplified spacing.
///
/// # Arguments
/// * `n_wake` - Number of wake stations
/// * `wake_length` - Total wake length in chord lengths
///
/// # Returns
/// Vector of arc lengths for wake stations
pub fn generate_wake_positions(n_wake: usize, wake_length: f64) -> Vec<f64> {
    if n_wake == 0 {
        return Vec::new();
    }

    // Geometric spacing - finer near TE, coarser downstream
    let mut positions = Vec::with_capacity(n_wake);
    let ratio = 1.2_f64; // Growth ratio

    let mut x = 1.0; // Start at TE (x/c = 1)
    let mut dx = wake_length / (n_wake as f64 * 2.0); // Initial spacing

    for _ in 0..n_wake {
        x += dx;
        positions.push(x);
        dx *= ratio;
    }

    // Normalize to fit within wake_length
    let x_max = positions.last().copied().unwrap_or(1.0);
    let scale = (1.0 + wake_length) / x_max;
    for p in &mut positions {
        *p = 1.0 + (*p - 1.0) * scale;
    }

    positions
}

/// Estimate wake edge velocity at downstream position.
///
/// Wake velocity recovers toward freestream as x → ∞.
/// Uses simple potential flow approximation.
///
/// # Arguments
/// * `x` - Position (x/c, where TE is at 1.0)
/// * `ue_te` - Edge velocity at trailing edge
///
/// # Returns
/// Estimated edge velocity at position x
pub fn wake_edge_velocity(x: f64, ue_te: f64) -> f64 {
    // Velocity recovers toward 1.0 (freestream)
    // Simple exponential recovery
    let x_wake = (x - 1.0).max(0.0);
    let recovery_rate = 2.0; // Empirical
    ue_te + (1.0 - ue_te) * (1.0 - (-recovery_rate * x_wake).exp())
}

/// Compute wake inverse mode Hk target.
///
/// In the wake, the shape factor Hk approaches 1.0 asymptotically
/// (fully developed wake). XFOIL uses a backward Euler iteration.
///
/// # Arguments
/// * `hk_prev` - Shape factor at previous station
/// * `dx` - Step size
/// * `theta_prev` - Momentum thickness at previous station
///
/// # Returns
/// Target Hk for inverse mode
///
/// # XFOIL Reference
/// xbl.f lines 735-746
pub fn wake_hk_target(hk_prev: f64, dx: f64, theta_prev: f64) -> f64 {
    let theta_safe = theta_prev.max(1e-12);
    let const_term = 0.03 * dx / theta_safe;

    // Backward Euler for Hk → 1.0
    // Solve: hk2 = hk1 - const * (hk2 - 1)^3
    // Using Newton iteration (one step approximation)
    let hk2 = hk_prev
        - (hk_prev + const_term * (hk_prev - 1.0).powi(3) - hk_prev)
            / (1.0 + 3.0 * const_term * (hk_prev - 1.0).powi(2));

    // Wake floor - Hk cannot go below 1.00005
    hk2.max(1.00005)
}

/// March wake from trailing edge downstream.
///
/// Uses simplified wake marching that solves the BL equations with
/// FlowType::Wake (no wall friction, wake dissipation model).
///
/// # Arguments
/// * `initial` - Initial wake station (from combine_te_for_wake)
/// * `wake_x` - Arc length positions for wake stations
/// * `wake_ue` - Edge velocities at wake positions
/// * `re` - Reynolds number
/// * `msq` - Mach number squared
///
/// # Returns
/// Vector of wake BlStations
///
/// # XFOIL Reference
/// xbl.f MRCHUE with WAKE=.TRUE. (lines 607-927)
pub fn march_wake(
    initial: &BlStation,
    wake_x: &[f64],
    wake_ue: &[f64],
    re: f64,
    msq: f64,
) -> Vec<BlStation> {
    if wake_x.is_empty() || wake_x.len() != wake_ue.len() {
        return vec![initial.clone()];
    }

    let mut stations = Vec::with_capacity(wake_x.len() + 1);
    stations.push(initial.clone());

    let mut prev = initial.clone();
    
    for i in 0..wake_x.len() {
        let x_new = wake_x[i];
        let ue_new = wake_ue[i].abs().max(0.01);
        
        // Solve wake station using simplified direct mode
        let mut station = solve_wake_station(&prev, x_new, ue_new, re, msq);
        
        // Compute secondary variables
        blvar(&mut station, FlowType::Wake, msq, re);
        
        stations.push(station.clone());
        prev = station;
    }

    stations
}

/// Solve for a single wake station using simplified Newton iteration.
///
/// Wake mode uses:
/// - Cf = 0 (no wall)
/// - Wake dissipation model
/// - Hk → 1.0 asymptotic behavior
fn solve_wake_station(
    prev: &BlStation,
    x_new: f64,
    ue_new: f64,
    re: f64,
    msq: f64,
) -> BlStation {
    let mut station = BlStation::default();
    station.x = x_new;
    station.x_coord = x_new;
    station.u = ue_new;
    station.is_wake = true;
    station.is_turbulent = true;
    station.is_laminar = false;
    station.cf = 0.0; // No wall friction in wake
    
    // Initialize from previous
    station.theta = prev.theta;
    station.delta_star = prev.delta_star;
    station.ctau = prev.ctau;
    
    // Simple step: use momentum equation with Cf = 0
    // dθ/dx = -(H + 2 - M²) * θ * (dUe/dx) / Ue
    let dx = (x_new - prev.x).max(1e-10);
    let due_dx = (ue_new - prev.u) / dx;
    let h = prev.h.clamp(1.0, 4.0);
    
    // Wake momentum growth (simplified)
    let dtheta_dx = -(h + 2.0 - msq) * prev.theta * due_dx / prev.u.max(0.01);
    station.theta = (prev.theta + dtheta_dx * dx).max(1e-8);
    
    // Wake Hk approaches 1.0 asymptotically
    let hk_target = wake_hk_target(prev.hk, dx, prev.theta);
    
    // Compute delta_star from target Hk
    // Hk ≈ H for low Mach, so δ* = Hk * θ
    station.delta_star = hk_target * station.theta;
    
    // Update shape factors
    station.h = station.delta_star / station.theta.max(1e-12);
    station.hk = hk_target;
    
    // Ctau decays in wake
    let ctau_decay = 0.9_f64; // Exponential decay factor
    station.ctau = (prev.ctau * ctau_decay).max(0.001);
    
    // Mass defect
    station.mass_defect = station.u * station.delta_star;
    
    // Rtheta
    station.r_theta = re * station.u * station.theta;
    
    station
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_combine_te_for_wake() {
        let mut upper = BlStation::default();
        upper.theta = 0.002;
        upper.delta_star = 0.005;
        upper.ctau = 0.03;
        upper.u = 0.95;

        let mut lower = BlStation::default();
        lower.theta = 0.0018;
        lower.delta_star = 0.004;
        lower.ctau = 0.025;
        lower.u = 0.93;

        let wake = combine_te_for_wake(&upper, &lower);

        // Theta adds
        assert!((wake.theta - 0.0038).abs() < 1e-10);
        // Delta star adds
        assert!((wake.delta_star - 0.009).abs() < 1e-10);
        // Wake flags
        assert!(wake.is_wake);
        assert!(wake.is_turbulent);
        assert!(!wake.is_laminar);
        // Cf is zero
        assert!((wake.cf - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_generate_wake_positions() {
        let positions = generate_wake_positions(5, 1.0);
        assert_eq!(positions.len(), 5);
        // All positions should be > 1.0 (after TE)
        assert!(positions.iter().all(|&x| x > 1.0));
        // Should end at approximately 2.0 (1.0 + wake_length)
        assert!(positions.last().unwrap() <= &2.0);
    }

    #[test]
    fn test_wake_hk_target() {
        // Starting Hk in wake
        let hk_prev = 1.3;
        let dx = 0.1;
        let theta = 0.003;

        let hk_target = wake_hk_target(hk_prev, dx, theta);

        // Should be less than hk_prev (approaching 1.0)
        assert!(hk_target < hk_prev);
        // Should be greater than 1.00005 (floor)
        assert!(hk_target >= 1.00005);
    }
}
