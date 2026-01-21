//! Solution update with limiting for viscous-inviscid coupling
//!
//! This module implements XFOIL's UPDATE and UESET routines which:
//! - Apply Newton correction deltas to BL variables with under-relaxation
//! - Compute edge velocities from inviscid solution plus mass defect influence
//! - Limit changes to prevent instability and non-physical states
//!
//! XFOIL Reference: xbl.f UPDATE (line 1253), xpanel.f UESET (line 1758)

use nalgebra::DMatrix;
use rustfoil_bl::closures::hkin::hkin;
use rustfoil_bl::state::BlStation;

/// Configuration for update limiting
///
/// These parameters control the under-relaxation and limiting applied
/// during Newton iterations to ensure stability and convergence.
#[derive(Debug, Clone)]
pub struct UpdateConfig {
    /// Maximum relative change in θ (momentum thickness)
    pub max_theta_change: f64,
    /// Maximum relative change in δ* (displacement thickness)
    pub max_delta_star_change: f64,
    /// Maximum relative change in Cτ (shear stress coefficient)
    pub max_ctau_change: f64,
    /// Maximum absolute change in Ue (normalized by 0.25)
    pub max_ue_change: f64,
    /// Maximum relative increase (DHI in XFOIL)
    pub max_relative_increase: f64,
    /// Maximum relative decrease (DLO in XFOIL)
    pub max_relative_decrease: f64,
    /// Maximum Cτ value to prevent blow-up
    pub ctau_max: f64,
    /// Relaxation factor (0-1), 1.0 = no relaxation
    pub relaxation: f64,
    /// Minimum kinematic shape factor for attached flow
    pub hk_min_attached: f64,
    /// Minimum kinematic shape factor in wake
    pub hk_min_wake: f64,
}

impl Default for UpdateConfig {
    fn default() -> Self {
        Self {
            max_theta_change: 0.3,
            max_delta_star_change: 0.3,
            max_ctau_change: 0.3,
            max_ue_change: 0.25,
            // XFOIL values: DHI = 1.5, DLO = -0.5
            max_relative_increase: 1.5,
            max_relative_decrease: -0.5,
            ctau_max: 0.25,
            relaxation: 1.0,
            // XFOIL: HKLIM = 1.02 for attached, 1.00005 for wake
            hk_min_attached: 1.02,
            hk_min_wake: 1.00005,
        }
    }
}

/// Result of update operation including convergence metrics
#[derive(Debug, Clone)]
pub struct UpdateResult {
    /// RMS of normalized changes
    pub rms_change: f64,
    /// Maximum normalized change
    pub max_change: f64,
    /// Variable with maximum change ('T' = theta, 'D' = delta*, 'C' = ctau, 'U' = Ue)
    pub max_change_var: char,
    /// Station index with maximum change
    pub max_change_station: usize,
    /// Under-relaxation factor actually used
    pub relaxation_used: f64,
}

/// Apply Newton updates to BL stations with limiting
///
/// This implements XFOIL's UPDATE subroutine which:
/// 1. Computes the change in each BL variable from Newton deltas
/// 2. Normalizes changes to determine if under-relaxation is needed
/// 3. Applies changes with relaxation factor
/// 4. Limits δ* to prevent Hk from dropping below physical limits
///
/// # Arguments
/// * `stations` - BL stations to update
/// * `deltas` - Newton solution `[Δcτ/Δn, Δθ, Δmass]` at each station
/// * `ue_new` - New edge velocities after update
/// * `config` - Update limiting configuration
///
/// # Returns
/// Update metrics for convergence checking
///
/// # Reference
/// XFOIL xbl.f UPDATE (line 1253)
pub fn update_stations(
    stations: &mut [BlStation],
    deltas: &[[f64; 3]],
    ue_new: &[f64],
    config: &UpdateConfig,
) -> UpdateResult {
    let n = stations.len();
    if n == 0 || deltas.len() != n || ue_new.len() != n {
        return UpdateResult {
            rms_change: 0.0,
            max_change: 0.0,
            max_change_var: ' ',
            max_change_station: 0,
            relaxation_used: 1.0,
        };
    }

    let dhi = config.max_relative_increase;
    let dlo = config.max_relative_decrease;

    // First pass: compute all changes and determine required relaxation
    let mut rlx = config.relaxation;
    let mut rms_sum = 0.0;
    let mut max_change = 0.0_f64;
    let mut max_change_var = ' ';
    let mut max_change_station = 0;

    // Pre-compute changes to determine relaxation factor
    let mut changes: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(n);

    for i in 0..n {
        let station = &stations[i];
        let delta = &deltas[i];

        // Extract Newton deltas
        // delta[0] = Δcτ (turbulent) or Δn (laminar)
        // delta[1] = Δθ
        // delta[2] = Δmass = Δ(Ue * δ*)
        let dctau = delta[0];
        let dthet = delta[1];
        let dmass = delta[2];

        // Compute Ue change
        let duedg = ue_new[i] - station.u;

        // Compute δ* change from mass defect change
        // mass = Ue * δ*, so Δδ* = (Δmass - δ* * ΔUe) / Ue
        let ddstr = if station.u.abs() > 1e-12 {
            (dmass - station.delta_star * duedg) / station.u
        } else {
            0.0
        };

        // Normalize changes
        let dn1 = if station.is_laminar {
            dctau / 10.0 // Amplification factor change scaled
        } else if station.ctau.abs() > 1e-12 {
            dctau / station.ctau
        } else {
            dctau / 0.01
        };

        let dn2 = if station.theta.abs() > 1e-12 {
            dthet / station.theta
        } else {
            0.0
        };

        let dn3 = if station.delta_star.abs() > 1e-12 {
            ddstr / station.delta_star
        } else {
            0.0
        };

        let dn4 = duedg.abs() / config.max_ue_change;

        // Accumulate for RMS
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

        // Track maximum change
        if dn1.abs() > max_change.abs() {
            max_change = dn1;
            max_change_var = if station.is_laminar { 'n' } else { 'C' };
            max_change_station = i;
        }
        if dn2.abs() > max_change.abs() {
            max_change = dn2;
            max_change_var = 'T';
            max_change_station = i;
        }
        if dn3.abs() > max_change.abs() {
            max_change = dn3;
            max_change_var = 'D';
            max_change_station = i;
        }
        if dn4.abs() > max_change.abs() {
            max_change = duedg;
            max_change_var = 'U';
            max_change_station = i;
        }

        // Check if under-relaxation needed for each variable
        let rdn1 = rlx * dn1;
        if rdn1 > dhi {
            rlx = dhi / dn1;
        }
        if rdn1 < dlo {
            rlx = dlo / dn1;
        }

        let rdn2 = rlx * dn2;
        if rdn2 > dhi {
            rlx = dhi / dn2;
        }
        if rdn2 < dlo {
            rlx = dlo / dn2;
        }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi {
            rlx = dhi / dn3;
        }
        if rdn3 < dlo {
            rlx = dlo / dn3;
        }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi {
            rlx = dhi / dn4;
        }
        if rdn4 < dlo {
            rlx = dlo / dn4;
        }

        changes.push((dctau, dthet, ddstr, duedg));
    }

    // Compute RMS change
    let rms_change = (rms_sum / (4.0 * n as f64)).sqrt();

    // Second pass: apply updates with computed relaxation factor
    for i in 0..n {
        let (dctau, dthet, ddstr, duedg) = changes[i];

        // Apply relaxed updates
        if stations[i].is_laminar {
            stations[i].ampl += rlx * dctau;
            stations[i].ampl = stations[i].ampl.max(0.0);
        } else {
            stations[i].ctau += rlx * dctau;
            // Clamp Cτ to prevent blow-up
            stations[i].ctau = stations[i].ctau.clamp(0.0, config.ctau_max);
        }

        stations[i].theta += rlx * dthet;
        stations[i].delta_star += rlx * ddstr;
        stations[i].u += rlx * duedg;

        // Limit δ* to prevent Hk from dropping below physical minimum
        let hk_limit = if stations[i].is_wake {
            config.hk_min_wake
        } else {
            config.hk_min_attached
        };

        limit_delta_star_for_hk(&mut stations[i], hk_limit, 0.0);

        // Update mass defect (nonlinear update)
        stations[i].mass_defect = stations[i].delta_star * stations[i].u;

        // Update shape factor
        if stations[i].theta.abs() > 1e-12 {
            stations[i].h = stations[i].delta_star / stations[i].theta;
        }
    }

    // Ensure no negative Ue "islands" (separate pass to avoid borrow conflicts)
    for i in 1..n {
        if stations[i - 1].u > 0.0 && stations[i].u <= 0.0 {
            stations[i].u = stations[i - 1].u;
            stations[i].mass_defect = stations[i].delta_star * stations[i].u;
        }
    }

    UpdateResult {
        rms_change,
        max_change,
        max_change_var,
        max_change_station,
        relaxation_used: rlx,
    }
}

/// Limit δ* to prevent Hk from dropping below a minimum value
///
/// This implements XFOIL's DSLIM subroutine which prevents the boundary
/// layer from reaching an unphysical separated state by limiting the
/// kinematic shape factor.
///
/// # Arguments
/// * `station` - BL station to potentially modify
/// * `hk_min` - Minimum allowed kinematic shape factor
/// * `msq` - Local Mach number squared
///
/// # Reference
/// XFOIL xbl.f DSLIM (line 1564)
pub fn limit_delta_star_for_hk(station: &mut BlStation, hk_min: f64, msq: f64) {
    if station.theta.abs() < 1e-12 {
        return;
    }

    let h = station.delta_star / station.theta;
    let result = hkin(h, msq);

    if result.hk < hk_min && result.hk_h.abs() > 1e-12 {
        // Increase δ* to bring Hk up to minimum
        let dh = (hk_min - result.hk) / result.hk_h;
        station.delta_star += dh * station.theta;
    }
}

/// Limit a change to prevent instability
///
/// Clamps the change to be within ±max_relative of the current value.
///
/// # Arguments
/// * `delta` - Proposed change
/// * `current` - Current value
/// * `max_relative` - Maximum relative change (e.g., 0.3 for 30%)
///
/// # Returns
/// Limited change value
#[inline]
pub fn limit_change(delta: f64, current: f64, max_relative: f64) -> f64 {
    let max_abs = max_relative * current.abs();
    delta.clamp(-max_abs, max_abs)
}

/// Set edge velocities from inviscid solution plus mass defect influence
///
/// This implements XFOIL's UESET subroutine which computes the edge
/// velocity at each BL station as:
///
///   Ue = Ue_inviscid + Σ DIJ * (Ue * δ*)
///
/// The DIJ matrix represents the influence of mass defect (source strength)
/// at station j on the edge velocity at station i.
///
/// # Arguments
/// * `stations` - BL stations (read mass_defect, write u)
/// * `ue_inviscid` - Edge velocities from inviscid solver
/// * `dij` - Mass defect influence matrix
/// * `vti` - Sign array (+1 or -1) for upper/lower surface direction
///
/// # Reference
/// XFOIL xpanel.f UESET (line 1758)
pub fn set_edge_velocities(
    stations: &mut [BlStation],
    ue_inviscid: &[f64],
    dij: &DMatrix<f64>,
    vti: Option<&[f64]>,
) {
    let n = stations.len();
    if n == 0 || ue_inviscid.len() != n {
        return;
    }

    // Default VTI to 1.0 if not provided (all same direction)
    let default_vti = vec![1.0; n];
    let vti = vti.unwrap_or(&default_vti);

    // XFOIL formula:
    // UE_M = -VTI(i) * VTI(j) * DIJ(i,j)
    // DUI = sum_j UE_M * MASS(j)
    // UEDG(i) = UINV(i) + DUI

    for i in 0..n {
        let mut dui = 0.0;

        for j in 0..n {
            // Mass defect at station j
            let mass_j = stations[j].mass_defect;

            // Influence coefficient with sign correction
            let ue_m = -vti[i] * vti[j] * dij[(i, j)];

            dui += ue_m * mass_j;
        }

        stations[i].u = ue_inviscid[i] + dui;
    }

    // Update mass defect after Ue change
    for station in stations.iter_mut() {
        station.mass_defect = station.u * station.delta_star;
    }
}

/// Compute new edge velocities without modifying stations
///
/// This is useful for the UPDATE routine which needs to compute new Ue
/// values before deciding on relaxation.
///
/// # Arguments
/// * `stations` - BL stations (read mass_defect)
/// * `ue_inviscid` - Edge velocities from inviscid solver
/// * `dij` - Mass defect influence matrix
/// * `mass_deltas` - Changes to mass defect from Newton solution (optional)
/// * `vti` - Sign array for upper/lower surface direction
///
/// # Returns
/// Vector of new edge velocities
pub fn compute_new_edge_velocities(
    stations: &[BlStation],
    ue_inviscid: &[f64],
    dij: &DMatrix<f64>,
    mass_deltas: Option<&[f64]>,
    vti: Option<&[f64]>,
) -> Vec<f64> {
    let n = stations.len();
    if n == 0 || ue_inviscid.len() != n {
        return Vec::new();
    }

    let default_vti = vec![1.0; n];
    let vti = vti.unwrap_or(&default_vti);

    let mut ue_new = vec![0.0; n];

    for i in 0..n {
        let mut dui = 0.0;

        for j in 0..n {
            // Total mass at station j (current + delta if provided)
            let mass_j = if let Some(deltas) = mass_deltas {
                stations[j].mass_defect + deltas[j]
            } else {
                stations[j].mass_defect
            };

            let ue_m = -vti[i] * vti[j] * dij[(i, j)];
            dui += ue_m * mass_j;
        }

        ue_new[i] = ue_inviscid[i] + dui;
    }

    ue_new
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_limit_change_large_positive() {
        // Large positive change should be limited
        let limited = limit_change(0.5, 1.0, 0.3);
        assert!((limited - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_limit_change_large_negative() {
        // Large negative change should be limited
        let limited = limit_change(-0.5, 1.0, 0.3);
        assert!((limited - (-0.3)).abs() < 1e-10);
    }

    #[test]
    fn test_limit_change_small() {
        // Small change should pass through
        let limited = limit_change(0.1, 1.0, 0.3);
        assert!((limited - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_limit_change_zero_current() {
        // Zero current value should give zero limit
        let limited = limit_change(0.5, 0.0, 0.3);
        assert_eq!(limited, 0.0);
    }

    #[test]
    fn test_limit_change_negative_current() {
        // Negative current value should still limit based on magnitude
        let limited = limit_change(0.5, -1.0, 0.3);
        assert!((limited - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_update_config_default() {
        let config = UpdateConfig::default();

        // Check XFOIL default values
        assert!((config.max_relative_increase - 1.5).abs() < 1e-10);
        assert!((config.max_relative_decrease - (-0.5)).abs() < 1e-10);
        assert!((config.ctau_max - 0.25).abs() < 1e-10);
        assert!((config.hk_min_attached - 1.02).abs() < 1e-10);
        assert!((config.hk_min_wake - 1.00005).abs() < 1e-10);
    }

    #[test]
    fn test_update_stations_empty() {
        let mut stations: Vec<BlStation> = vec![];
        let deltas: Vec<[f64; 3]> = vec![];
        let ue_new: Vec<f64> = vec![];
        let config = UpdateConfig::default();

        let result = update_stations(&mut stations, &deltas, &ue_new, &config);

        assert_eq!(result.rms_change, 0.0);
        assert_eq!(result.relaxation_used, 1.0);
    }

    #[test]
    fn test_update_stations_small_change() {
        // Create a single station
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].u = 1.0;
        stations[0].ctau = 0.1;
        stations[0].is_laminar = false;

        // Small delta that shouldn't trigger limiting
        let deltas = vec![[0.001, 0.0001, 0.0001]];
        let ue_new = vec![1.001];
        let config = UpdateConfig::default();

        let result = update_stations(&mut stations, &deltas, &ue_new, &config);

        // Should use full relaxation (no limiting needed)
        assert!(result.relaxation_used > 0.9);
        assert!(result.rms_change > 0.0);
    }

    #[test]
    fn test_update_stations_large_change_triggers_limiting() {
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].u = 1.0;
        stations[0].ctau = 0.1;
        stations[0].is_laminar = false;

        // Large delta that should trigger limiting
        let deltas = vec![[0.5, 0.05, 0.1]]; // 500% change in ctau
        let ue_new = vec![1.5]; // 50% change in Ue
        let config = UpdateConfig::default();

        let result = update_stations(&mut stations, &deltas, &ue_new, &config);

        // Should have reduced relaxation
        assert!(
            result.relaxation_used < 1.0,
            "Large change should trigger under-relaxation"
        );
    }

    #[test]
    fn test_update_stations_ctau_clamping() {
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].u = 1.0;
        stations[0].ctau = 0.2;
        stations[0].is_laminar = false;

        // Delta that would push ctau above maximum
        let deltas = vec![[0.2, 0.0, 0.0]];
        let ue_new = vec![1.0];
        let config = UpdateConfig::default();

        update_stations(&mut stations, &deltas, &ue_new, &config);

        // ctau should be clamped to maximum
        assert!(
            stations[0].ctau <= config.ctau_max,
            "ctau should be clamped to max"
        );
    }

    #[test]
    fn test_update_stations_laminar() {
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].u = 1.0;
        stations[0].ampl = 5.0;
        stations[0].is_laminar = true;

        let deltas = vec![[1.0, 0.0001, 0.0001]]; // Amplification change
        let ue_new = vec![1.0];
        let config = UpdateConfig::default();

        update_stations(&mut stations, &deltas, &ue_new, &config);

        // Amplification should have increased
        assert!(stations[0].ampl > 5.0, "Amplification should increase");
        assert!(
            stations[0].ampl >= 0.0,
            "Amplification should remain non-negative"
        );
    }

    #[test]
    fn test_limit_delta_star_for_hk() {
        let mut station = BlStation::new();
        station.theta = 0.01;
        station.delta_star = 0.0101; // H = 1.01, which gives Hk ≈ 1.01 at M=0
        station.u = 1.0;

        // Try to limit Hk to 1.02
        limit_delta_star_for_hk(&mut station, 1.02, 0.0);

        // δ* should have increased to bring Hk up to 1.02
        let h_new = station.delta_star / station.theta;
        assert!(
            h_new >= 1.02 - 1e-10,
            "H should be at least 1.02, got {}",
            h_new
        );
    }

    #[test]
    fn test_set_edge_velocities_basic() {
        let mut stations = vec![BlStation::new(); 3];
        for (i, station) in stations.iter_mut().enumerate() {
            station.u = 1.0;
            station.delta_star = 0.01;
            station.mass_defect = station.u * station.delta_star;
            station.x = i as f64 * 0.1;
        }

        let ue_inviscid = vec![1.0, 1.05, 1.02];
        let dij = DMatrix::zeros(3, 3); // No influence for simplicity

        set_edge_velocities(&mut stations, &ue_inviscid, &dij, None);

        // With zero DIJ, Ue should equal inviscid values
        for (i, station) in stations.iter().enumerate() {
            assert!(
                (station.u - ue_inviscid[i]).abs() < 1e-10,
                "With zero DIJ, Ue should equal Ue_inviscid"
            );
        }
    }

    #[test]
    fn test_set_edge_velocities_with_influence() {
        let mut stations = vec![BlStation::new(); 2];
        stations[0].mass_defect = 0.01;
        stations[1].mass_defect = 0.02;

        let ue_inviscid = vec![1.0, 1.0];

        // Simple influence matrix
        let mut dij = DMatrix::zeros(2, 2);
        dij[(0, 1)] = -0.1; // Station 1 influences station 0
        dij[(1, 0)] = 0.1; // Station 0 influences station 1

        set_edge_velocities(&mut stations, &ue_inviscid, &dij, None);

        // Check that Ue has been modified by the influence
        // With default vti=1.0, UE_M = -1*1*DIJ = -DIJ
        // So station 0: Ue = 1.0 + (-1)*(-0.1)*0.02 = 1.0 + 0.002 = 1.002
        // Station 1: Ue = 1.0 + (-1)*(0.1)*0.01 = 1.0 - 0.001 = 0.999
        assert!(
            (stations[0].u - 1.002).abs() < 1e-10,
            "Station 0 Ue incorrect: {}",
            stations[0].u
        );
        assert!(
            (stations[1].u - 0.999).abs() < 1e-10,
            "Station 1 Ue incorrect: {}",
            stations[1].u
        );
    }

    #[test]
    fn test_compute_new_edge_velocities() {
        let stations = vec![BlStation::new(); 2];
        let ue_inviscid = vec![1.0, 1.1];
        let dij = DMatrix::zeros(2, 2);

        let ue_new = compute_new_edge_velocities(&stations, &ue_inviscid, &dij, None, None);

        assert_eq!(ue_new.len(), 2);
        assert!((ue_new[0] - 1.0).abs() < 1e-10);
        assert!((ue_new[1] - 1.1).abs() < 1e-10);
    }

    #[test]
    fn test_compute_new_edge_velocities_with_deltas() {
        let mut stations = vec![BlStation::new(); 2];
        stations[0].mass_defect = 0.01;
        stations[1].mass_defect = 0.01;

        let ue_inviscid = vec![1.0, 1.0];
        let dij = DMatrix::zeros(2, 2);
        let mass_deltas = vec![0.001, 0.002];

        let ue_new =
            compute_new_edge_velocities(&stations, &ue_inviscid, &dij, Some(&mass_deltas), None);

        // With zero DIJ, result should still be inviscid values
        assert_eq!(ue_new.len(), 2);
        assert!((ue_new[0] - 1.0).abs() < 1e-10);
        assert!((ue_new[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_update_result_fields() {
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].u = 1.0;
        stations[0].ctau = 0.1;
        stations[0].is_laminar = false;

        let deltas = vec![[0.01, 0.001, 0.001]];
        let ue_new = vec![1.01];
        let config = UpdateConfig::default();

        let result = update_stations(&mut stations, &deltas, &ue_new, &config);

        // Verify all result fields are populated
        assert!(result.rms_change >= 0.0);
        assert!(result.relaxation_used > 0.0 && result.relaxation_used <= 1.0);
        assert!(result.max_change_station < stations.len());
        // max_change_var should be one of the valid characters
        assert!(
            result.max_change_var == 'T'
                || result.max_change_var == 'D'
                || result.max_change_var == 'C'
                || result.max_change_var == 'U'
                || result.max_change_var == 'n'
        );
    }

    #[test]
    fn test_mass_defect_updated_after_update() {
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].u = 1.0;
        stations[0].mass_defect = 0.02; // u * delta_star
        stations[0].is_laminar = false;
        stations[0].ctau = 0.1;

        let deltas = vec![[0.0, 0.001, 0.005]];
        let ue_new = vec![1.05];
        let config = UpdateConfig::default();

        update_stations(&mut stations, &deltas, &ue_new, &config);

        // Mass defect should be updated to u * delta_star
        let expected_mass = stations[0].u * stations[0].delta_star;
        assert!(
            (stations[0].mass_defect - expected_mass).abs() < 1e-12,
            "Mass defect should equal u * delta_star after update"
        );
    }

    #[test]
    fn test_shape_factor_updated() {
        let mut stations = vec![BlStation::new()];
        stations[0].theta = 0.01;
        stations[0].delta_star = 0.02;
        stations[0].h = 2.0;
        stations[0].u = 1.0;
        stations[0].is_laminar = false;
        stations[0].ctau = 0.1;

        let deltas = vec![[0.0, 0.001, 0.003]];
        let ue_new = vec![1.0];
        let config = UpdateConfig::default();

        update_stations(&mut stations, &deltas, &ue_new, &config);

        // H should be updated to delta_star / theta
        let expected_h = stations[0].delta_star / stations[0].theta;
        assert!(
            (stations[0].h - expected_h).abs() < 1e-10,
            "Shape factor H should be updated"
        );
    }

    #[test]
    fn test_negative_ue_prevention() {
        let mut stations = vec![BlStation::new(); 3];
        for station in stations.iter_mut() {
            station.theta = 0.01;
            station.delta_star = 0.02;
            station.u = 1.0;
            station.is_laminar = false;
            station.ctau = 0.1;
        }

        // Try to make middle station negative while neighbors are positive
        let deltas = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let ue_new = vec![1.0, -0.5, 1.0]; // Negative in middle
        let config = UpdateConfig::default();

        update_stations(&mut stations, &deltas, &ue_new, &config);

        // Middle station should not go negative if previous is positive
        assert!(
            stations[1].u > 0.0,
            "Ue should not go negative between positive neighbors"
        );
    }
}
