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
/// In XFOIL, wake Ue comes from inviscid panel solution which naturally
/// shows slow recovery. Here we approximate this behavior.
pub fn wake_edge_velocity(x: f64, ue_te: f64) -> f64 {
    let x_wake = (x - 1.0).max(0.0);
    
    // XFOIL wake physics:
    // 1. Near wake (x < ~1.5c): Ue stays close to TE value
    //    - Wake mixing creates velocity deficit region
    //    - Recovery is very slow due to wake spreading
    // 2. Far wake (x > ~1.5c): Ue slowly approaches freestream
    //    - Asymptotic approach to 1.0 as x → ∞
    //
    // Use delayed exponential recovery to match XFOIL behavior:
    // - Very slow recovery in near wake
    // - Faster (but still slow) recovery in far wake
    
    // XFOIL wake velocity recovery:
    // From inviscid potential flow, wake velocity recovers as:
    // Ue ≈ 1 - c₁/x + c₂/x² + ...
    // For practical purposes, use exponential recovery that reaches
    // ~0.99 by x = 2-3c downstream of TE.
    //
    // At TE (x=1), Ue = ue_te ≈ 0.8
    // At far wake (x=3), Ue should approach 0.99
    
    let near_wake_length = 0.2; // Short near-wake region
    let recovery_rate = 1.5;    // Faster recovery to reach ~0.99 by x=3
    
    if x_wake < near_wake_length {
        // Near wake: Ue stays close to TE value
        ue_te
    } else {
        // Far wake: exponential recovery toward freestream
        let x_far = x_wake - near_wake_length;
        let recovery_factor = 1.0 - (-recovery_rate * x_far).exp();
        ue_te + (1.0 - ue_te) * recovery_factor
    }
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

/// Solve for a single wake station using momentum conservation.
///
/// Wake physics (XFOIL-style):
/// - Cf = 0 (no wall)
/// - Momentum approximately conserved: θ * Ue ≈ constant
/// - Hk → 1.0 asymptotic behavior (thin wake limit)
/// - Ctau decays due to turbulent dissipation
///
/// # Key insight for Squire-Young
/// The Squire-Young formula CD = 2θ * Ue^((5+H)/2) is sensitive to far-wake
/// values. With momentum conservation (θ * Ue ≈ const):
/// - If Ue stays near TE value (~0.8), θ stays near TE value (~0.007)
/// - Far downstream with Ue → 1, θ → θ_te * Ue_te (slightly smaller)
/// - Hk approaches 1.0, so exponent (5+H)/2 → 3.0
///
/// This is critical: at TE, Hk≈2.5 gives exponent≈3.75, while far-wake
/// Hk≈1.04 gives exponent≈3.02. The ~50% CD overestimate comes from using
/// TE values instead of far-wake values.
fn solve_wake_station(
    prev: &BlStation,
    x_new: f64,
    ue_new: f64,
    re: f64,
    _msq: f64,
) -> BlStation {
    let mut station = BlStation::default();
    station.x = x_new;
    station.x_coord = x_new;
    station.u = ue_new;
    station.is_wake = true;
    station.is_turbulent = true;
    station.is_laminar = false;
    station.cf = 0.0; // No wall friction in wake
    
    let dx = (x_new - prev.x).max(1e-10);
    let ue_prev = prev.u.max(0.01);
    let ue_curr = ue_new.max(0.01);
    
    // === Wake momentum evolution with dissipation ===
    // XFOIL data shows θ*Ue² is NOT strictly conserved in the wake.
    // Analysis of XFOIL wake evolution shows:
    //   - Combined TE: θ=0.0066, Ue=0.77, θ*Ue²=0.0039
    //   - Far wake:    θ=0.0026, Ue=0.99, θ*Ue²=0.0025
    // This is a 35% reduction in θ*Ue², indicating dissipation.
    //
    // The wake momentum equation with Cf=0 is:
    //   dθ/dx = -θ * (H + 2 - M²) * (dUe/dx) / Ue + Cd
    // where Cd is the wake dissipation coefficient.
    //
    // XFOIL's wake dissipation leads to faster θ decay than pure momentum
    // conservation. Empirically, a decay rate of ~1.0 per chord length
    // matches XFOIL's wake evolution.
    
    // Compute nominal momentum-conserved theta
    let ue_ratio = ue_prev / ue_curr;
    let theta_conserved = prev.theta * ue_ratio * ue_ratio;
    
    // Apply incremental wake dissipation (per step, not total distance)
    // XFOIL wake at α=0° shows additional θ decay beyond momentum conservation.
    // Target: Total additional dissipation factor ≈ 0.69 over ~2.5c wake length
    // exp(-rate * 2.5) = 0.69  =>  rate = -ln(0.69)/2.5 ≈ 0.15
    let dissipation_rate = 0.15; // Per chord length (calibrated for 2.5c wake)
    let dissipation_factor = (-dissipation_rate * dx).exp();
    
    station.theta = theta_conserved * dissipation_factor;
    station.theta = station.theta.max(1e-8);
    
    // === Shape factor evolution ===
    // Wake Hk approaches 1.0 asymptotically (thin wake limit)
    // Use backward Euler: Hk2 = Hk1 - const * (Hk1 - 1)³
    let hk_target = wake_hk_target(prev.hk, dx, prev.theta);
    station.hk = hk_target;
    
    // Compute delta_star from Hk
    // For low Mach: H ≈ Hk, so δ* = H * θ ≈ Hk * θ
    station.delta_star = station.hk * station.theta;
    station.h = station.delta_star / station.theta.max(1e-12);
    
    // === Ctau decay ===
    // Turbulent shear stress coefficient decays in wake due to dissipation
    // Use exponential decay with length scale
    let decay_length = 1.0; // Characteristic length for ctau decay
    let ctau_decay = (-dx / decay_length).exp();
    station.ctau = (prev.ctau * ctau_decay).max(0.001);
    
    // === Derived quantities ===
    station.mass_defect = station.u * station.delta_star;
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
