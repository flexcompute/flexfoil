//! Wake marching and trailing edge combination.
//!
//! This module implements wake handling for viscous-inviscid coupling:
//! - Combining upper and lower surface BL at the trailing edge
//! - Marching the wake downstream using BLDIF(FlowType::Wake)
//!
//! # XFOIL Reference
//! - TESYS: xbl.f lines 236-276 (TE system setup)
//! - Wake marching: xbl.f lines 607-927 (same as MRCHUE but with wake flag)

use crate::solve::{invert_3x3, multiply_3x3_vec};
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
/// Wake velocity recovery follows XFOIL-style physics:
/// - Near-wake: Ue stays close to TE value (wake mixing region)
/// - Far-wake: Ue gradually recovers toward freestream
///
/// This is critical for Squire-Young CD computation - wrong Ue evolution
/// causes θ to decrease (via momentum conservation) instead of staying
/// approximately constant, leading to ~50% CD overestimation.
///
/// # Arguments
/// * `x` - Position (x/c, where TE is at 1.0)
/// * `ue_te` - Edge velocity at trailing edge
///
/// # Returns
/// Estimated edge velocity at position x
///
/// # XFOIL Reference
/// In XFOIL, wake Ue comes from the inviscid wake panel solution and recovers
/// algebraically toward the freestream. We use the same qualitative 1/x-style
/// recovery law here until the wake-point velocity extraction is fully wired in.
pub fn wake_edge_velocity(x: f64, ue_te: f64) -> f64 {
    let x_wake = (x - 1.0).max(0.0);

    if x_wake <= 1.0e-12 {
        return ue_te.max(0.01);
    }

    let deficit_te = (1.0 - ue_te).max(0.0);
    let recovered = 1.0 - deficit_te / (1.0 + 2.0 * x_wake + 0.5 * x_wake * x_wake);
    recovered.clamp(ue_te.min(1.0).max(0.01), 1.0)
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

/// Solve for a single wake station using the full wake BLDIF equations.
///
/// This mirrors XFOIL's wake MRCHUE path much more closely than the older
/// empirical wake evolution. The wake still receives its prescribed Ue profile
/// from the caller, but the state march itself is now a Newton solve on
/// `BLDIF(FlowType::Wake)` rather than an ad hoc decay model.
fn solve_wake_station(
    prev: &BlStation,
    x_new: f64,
    ue_new: f64,
    re: f64,
    msq: f64,
) -> BlStation {
    let mut station = BlStation::default();
    station.x = x_new;
    station.x_coord = 1.0;
    station.u = ue_new;
    station.is_wake = true;
    station.is_turbulent = true;
    station.is_laminar = false;

    let ue_prev = prev.u.max(0.01);
    let ue_curr = ue_new.max(0.01);
    let ue_ratio = ue_prev / ue_curr;

    station.theta = (prev.theta * ue_ratio).max(1e-8);
    station.delta_star = (prev.delta_star * ue_ratio).max(station.theta * 1.00005);
    station.ctau = prev.ctau.max(1e-4);

    let max_iter = 20;
    let tolerance = 1.0e-10;

    for _ in 0..max_iter {
        blvar(&mut station, FlowType::Wake, msq, re);
        let (residuals, jacobian) = bldif(prev, &station, FlowType::Wake, msq, re);

        let a = [
            [jacobian.vs2[0][0], jacobian.vs2[0][1], jacobian.vs2[0][2]],
            [jacobian.vs2[1][0], jacobian.vs2[1][1], jacobian.vs2[1][2]],
            [jacobian.vs2[2][0], jacobian.vs2[2][1], jacobian.vs2[2][2]],
        ];
        let rhs = [residuals.res_third, residuals.res_mom, residuals.res_shape];
        let delta = multiply_3x3_vec(&invert_3x3(&a), &rhs);

        let mut rlx: f64 = 1.0;
        if delta[0] < 0.0 {
            rlx = rlx.min(0.5 * station.ctau.max(1.0e-4) / (-delta[0]).max(1.0e-12));
        }
        if delta[1] < 0.0 {
            rlx = rlx.min(0.5 * station.theta.max(1.0e-8) / (-delta[1]).max(1.0e-12));
        }
        if delta[2] < 0.0 {
            rlx = rlx.min(
                0.5 * station.delta_star.max(1.0e-8) / (-delta[2]).max(1.0e-12),
            );
        }
        rlx = rlx.clamp(0.05, 1.0);

        station.ctau = (station.ctau + rlx * delta[0]).clamp(1.0e-6, 0.25);
        station.theta = (station.theta + rlx * delta[1]).max(1.0e-8);
        station.delta_star = (station.delta_star + rlx * delta[2]).max(station.theta * 1.00005);

        let max_rel = [
            delta[0].abs() / station.ctau.max(1.0e-4),
            delta[1].abs() / station.theta.max(1.0e-8),
            delta[2].abs() / station.delta_star.max(1.0e-8),
        ]
        .into_iter()
        .fold(0.0, f64::max);

        if max_rel < tolerance {
            break;
        }
    }

    blvar(&mut station, FlowType::Wake, msq, re);
    station.mass_defect = station.u * station.delta_star;
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
