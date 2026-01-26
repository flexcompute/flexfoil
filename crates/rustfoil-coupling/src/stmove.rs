//! Stagnation point relocation during viscous-inviscid iteration.
//!
//! This module implements XFOIL's STMOVE algorithm, which relocates
//! the stagnation point when the circulation distribution changes
//! during Newton iteration.
//!
//! # XFOIL Reference
//! - STFIND: xpanel.f lines 1357-1392
//! - STMOVE: xpanel.f lines 1628-1750

use rustfoil_bl::state::BlStation;

/// Find stagnation point from circulation (gamma) sign change.
///
/// XFOIL's STFIND algorithm: finds where gamma changes from positive
/// to negative, which indicates the stagnation point.
///
/// # Arguments
/// * `gamma` - Circulation distribution around airfoil
/// * `s` - Arc length array
///
/// # Returns
/// `(ist, sst)` - Panel index and exact arc length of stagnation point
///
/// # XFOIL Reference
/// xpanel.f STFIND (lines 1357-1392)
pub fn find_stagnation_by_gamma(gamma: &[f64], s: &[f64]) -> Option<(usize, f64)> {
    if gamma.len() < 2 || gamma.len() != s.len() {
        return None;
    }

    // Find where gamma changes from >= 0 to < 0
    for i in 0..gamma.len() - 1 {
        if gamma[i] >= 0.0 && gamma[i + 1] < 0.0 {
            // Found sign change - interpolate exact location
            let dgam = gamma[i + 1] - gamma[i];
            let ds = s[i + 1] - s[i];

            // Interpolate to minimize roundoff (XFOIL lines 1378-1382)
            let sst = if gamma[i] < -gamma[i + 1] {
                s[i] - ds * (gamma[i] / dgam)
            } else {
                s[i + 1] - ds * (gamma[i + 1] / dgam)
            };

            // Tweak if falls exactly on a node (very unlikely)
            let sst = if sst <= s[i] {
                s[i] + 1.0e-7
            } else if sst >= s[i + 1] {
                s[i + 1] - 1.0e-7
            } else {
                sst
            };

            return Some((i, sst));
        }
    }

    // No sign change found - stagnation at midpoint (fallback)
    Some((gamma.len() / 2, s[gamma.len() / 2]))
}

/// Shift BL stations downstream (increase indices).
///
/// Used when stagnation moves to add more stations at the beginning.
/// Stations are shifted so that station[i] gets the value of station[i-idif].
///
/// # Arguments
/// * `stations` - BL stations to shift
/// * `idif` - Number of positions to shift
fn shift_stations_downstream(stations: &mut [BlStation], idif: usize) {
    if idif == 0 || idif >= stations.len() {
        return;
    }

    // Shift from end to beginning to avoid overwriting
    for i in (idif..stations.len()).rev() {
        let src_idx = i - idif;
        // Copy primary BL variables
        stations[i].theta = stations[src_idx].theta;
        stations[i].delta_star = stations[src_idx].delta_star;
        stations[i].ctau = stations[src_idx].ctau;
        stations[i].ampl = stations[src_idx].ampl;
        stations[i].u = stations[src_idx].u;
        // Copy secondary variables
        stations[i].h = stations[src_idx].h;
        stations[i].hk = stations[src_idx].hk;
        stations[i].hs = stations[src_idx].hs;
        stations[i].r_theta = stations[src_idx].r_theta;
        stations[i].cf = stations[src_idx].cf;
        stations[i].mass_defect = stations[src_idx].mass_defect;
        // Copy flags
        stations[i].is_laminar = stations[src_idx].is_laminar;
        stations[i].is_turbulent = stations[src_idx].is_turbulent;
    }
}

/// Shift BL stations upstream (decrease indices).
///
/// Used when stagnation moves to remove stations from the beginning.
/// Stations are shifted so that station[i] gets the value of station[i+idif].
///
/// # Arguments
/// * `stations` - BL stations to shift
/// * `idif` - Number of positions to shift
fn shift_stations_upstream(stations: &mut [BlStation], idif: usize) {
    if idif == 0 || idif >= stations.len() {
        return;
    }

    // Shift from beginning to end to avoid overwriting
    for i in 0..stations.len() - idif {
        let src_idx = i + idif;
        // Copy primary BL variables
        stations[i].theta = stations[src_idx].theta;
        stations[i].delta_star = stations[src_idx].delta_star;
        stations[i].ctau = stations[src_idx].ctau;
        stations[i].ampl = stations[src_idx].ampl;
        stations[i].u = stations[src_idx].u;
        // Copy secondary variables
        stations[i].h = stations[src_idx].h;
        stations[i].hk = stations[src_idx].hk;
        stations[i].hs = stations[src_idx].hs;
        stations[i].r_theta = stations[src_idx].r_theta;
        stations[i].cf = stations[src_idx].cf;
        stations[i].mass_defect = stations[src_idx].mass_defect;
        // Copy flags
        stations[i].is_laminar = stations[src_idx].is_laminar;
        stations[i].is_turbulent = stations[src_idx].is_turbulent;
    }
}

/// Interpolate BL variables in gap between old and new stagnation.
///
/// When stagnation moves, there are stations between the old and new
/// positions that need interpolated values. XFOIL uses a linear
/// velocity gradient: DUDX = Ue/x.
///
/// # Arguments
/// * `stations` - BL stations (already shifted)
/// * `idif` - Number of gap stations to interpolate
/// * `reference_idx` - Index of first valid station after gap
///
/// # XFOIL Reference
/// xpanel.f lines 1675-1681, 1706-1718
fn interpolate_gap_stations(stations: &mut [BlStation], idif: usize, reference_idx: usize) {
    if idif == 0 || reference_idx >= stations.len() {
        return;
    }

    // Get reference values from first valid station
    let ref_station = &stations[reference_idx];
    let ref_ue = ref_station.u;
    let ref_x = ref_station.x.max(1e-12);
    let ref_theta = ref_station.theta;
    let ref_dstar = ref_station.delta_star;
    let ref_ctau = ref_station.ctau;
    let ref_h = ref_station.h;

    // DUDX = Ue / x (velocity gradient at stagnation)
    let dudx = ref_ue / ref_x;

    // Interpolate gap stations
    for i in 1..=idif {
        let gap_idx = reference_idx - idif + i - 1;
        if gap_idx >= stations.len() {
            continue;
        }

        let station = &mut stations[gap_idx];
        let x = station.x.max(1e-12);

        // Linear velocity profile: Ue = DUDX * x
        station.u = dudx * x;
        station.u = station.u.max(1e-6); // Ensure positive

        // Keep BL variables from reference (XFOIL copies from first valid station)
        station.theta = ref_theta;
        station.delta_star = ref_dstar;
        station.ctau = ref_ctau;
        station.h = ref_h;

        // Update mass defect
        station.mass_defect = station.u * station.delta_star;
    }
}

/// Relocate stagnation point based on updated circulation.
///
/// This is the main STMOVE implementation. It:
/// 1. Finds the new stagnation point from gamma sign change
/// 2. If stagnation moved, shifts BL arrays appropriately
/// 3. Interpolates BL variables in the gap
///
/// # Arguments
/// * `upper_stations` - Upper surface BL stations
/// * `lower_stations` - Lower surface BL stations
/// * `gamma` - Updated circulation distribution
/// * `arc_lengths` - Full arc length array
/// * `old_ist` - Previous stagnation panel index
///
/// # Returns
/// `Some(new_ist)` if stagnation moved, `None` if unchanged
///
/// # XFOIL Reference
/// xpanel.f STMOVE (lines 1628-1750)
pub fn stmove(
    upper_stations: &mut [BlStation],
    lower_stations: &mut [BlStation],
    gamma: &[f64],
    arc_lengths: &[f64],
    old_ist: usize,
) -> Option<usize> {
    // 1. Find new stagnation from gamma sign change
    let (new_ist, _new_sst) = find_stagnation_by_gamma(gamma, arc_lengths)?;

    if new_ist == old_ist {
        return None; // No change needed
    }

    let idif = (new_ist as i32 - old_ist as i32).abs() as usize;

    // Limit maximum shift to avoid wild oscillations
    let idif = idif.min(3);

    if idif == 0 {
        return None;
    }

    // 2. Shift BL arrays based on direction
    if new_ist > old_ist {
        // Stagnation moved downstream - more points on upper surface
        // Upper: shift downstream (add stations at beginning)
        // Lower: shift upstream (remove stations from beginning)
        shift_stations_downstream(upper_stations, idif);
        shift_stations_upstream(lower_stations, idif);

        // 3. Interpolate gap stations
        if idif + 1 < upper_stations.len() {
            interpolate_gap_stations(upper_stations, idif, idif + 1);
        }
    } else {
        // Stagnation moved upstream - more points on lower surface
        // Upper: shift upstream (remove stations from beginning)
        // Lower: shift downstream (add stations at beginning)
        shift_stations_upstream(upper_stations, idif);
        shift_stations_downstream(lower_stations, idif);

        // 3. Interpolate gap stations
        if idif + 1 < lower_stations.len() {
            interpolate_gap_stations(lower_stations, idif, idif + 1);
        }
    }

    Some(new_ist)
}

/// Update transition indices after stagnation move.
///
/// When stagnation moves, the BL station indices shift, so transition
/// locations need adjustment.
///
/// # Arguments
/// * `x_tr` - Current transition x location
/// * `old_ist` - Old stagnation index
/// * `new_ist` - New stagnation index
/// * `is_upper` - True for upper surface
///
/// # Returns
/// Adjusted transition location (or same if no change needed)
pub fn adjust_transition_for_stmove(
    x_tr: f64,
    old_ist: usize,
    new_ist: usize,
    _is_upper: bool,
) -> f64 {
    // Transition x-location doesn't change, only the station index would
    // For now, keep the physical x location the same
    // The marching will find the correct station index
    if old_ist == new_ist {
        return x_tr;
    }

    // Physical location stays the same - the BL march will handle it
    x_tr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_stagnation_by_gamma() {
        // Gamma positive on one side, negative on other
        let gamma = vec![0.5, 0.2, 0.0, -0.1, -0.3];
        let s = vec![0.0, 0.25, 0.5, 0.75, 1.0];

        let result = find_stagnation_by_gamma(&gamma, &s);
        assert!(result.is_some());

        let (ist, sst) = result.unwrap();
        assert_eq!(ist, 2); // Sign change between index 2 and 3
        assert!(sst > 0.5 && sst < 0.75);
    }

    #[test]
    fn test_find_stagnation_no_sign_change() {
        // All positive - fallback to midpoint
        let gamma = vec![0.5, 0.3, 0.2, 0.1, 0.05];
        let s = vec![0.0, 0.25, 0.5, 0.75, 1.0];

        let result = find_stagnation_by_gamma(&gamma, &s);
        assert!(result.is_some());
        let (ist, _) = result.unwrap();
        assert_eq!(ist, 2); // Midpoint
    }

    #[test]
    fn test_shift_downstream() {
        let mut stations: Vec<BlStation> = (0..5)
            .map(|i| {
                let mut s = BlStation::default();
                s.theta = i as f64 * 0.001;
                s.u = 0.5 + i as f64 * 0.1;
                s
            })
            .collect();

        shift_stations_downstream(&mut stations, 2);

        // Stations 2,3,4 should have values from 0,1,2
        assert!((stations[2].theta - 0.0).abs() < 1e-10);
        assert!((stations[3].theta - 0.001).abs() < 1e-10);
        assert!((stations[4].theta - 0.002).abs() < 1e-10);
    }

    #[test]
    fn test_shift_upstream() {
        let mut stations: Vec<BlStation> = (0..5)
            .map(|i| {
                let mut s = BlStation::default();
                s.theta = i as f64 * 0.001;
                s.u = 0.5 + i as f64 * 0.1;
                s
            })
            .collect();

        shift_stations_upstream(&mut stations, 2);

        // Stations 0,1,2 should have values from 2,3,4
        assert!((stations[0].theta - 0.002).abs() < 1e-10);
        assert!((stations[1].theta - 0.003).abs() < 1e-10);
        assert!((stations[2].theta - 0.004).abs() < 1e-10);
    }
}
