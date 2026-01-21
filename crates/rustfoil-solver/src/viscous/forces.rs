//! Aerodynamic force computation from viscous solution.
//!
//! This module computes lift, drag, and moment coefficients from the
//! boundary layer solution. The key components are:
//!
//! - **Friction drag (Cd_f)**: Integrated skin friction coefficient along surface
//! - **Pressure drag (Cd_p)**: Momentum deficit at wake trailing edge (Squire-Young)
//! - **Total drag**: Cd = Cd_f + Cd_p
//!
//! # XFOIL Reference
//! - CDCALC: xfoil.f line 1172

use rustfoil_bl::state::BlStation;

use super::config::ViscousSolverConfig;

/// Aerodynamic force coefficients computed from viscous solution.
#[derive(Debug, Clone, Copy, Default)]
pub struct AeroForces {
    /// Lift coefficient
    pub cl: f64,
    /// Total drag coefficient (friction + pressure)
    pub cd: f64,
    /// Moment coefficient (about quarter-chord)
    pub cm: f64,
    /// Pressure drag coefficient (from wake momentum deficit)
    pub cd_pressure: f64,
    /// Friction drag coefficient (integrated skin friction)
    pub cd_friction: f64,
}

impl AeroForces {
    /// Compute L/D ratio.
    pub fn lift_to_drag(&self) -> f64 {
        if self.cd.abs() > 1e-12 {
            self.cl / self.cd
        } else {
            f64::INFINITY
        }
    }

    /// Check if forces are physically reasonable.
    pub fn is_valid(&self) -> bool {
        self.cd >= 0.0
            && self.cd_friction >= 0.0
            && self.cd_pressure.is_finite()
            && self.cl.is_finite()
            && self.cm.is_finite()
    }
}

/// Compute aerodynamic forces from BL solution.
///
/// # Friction Drag
/// Friction drag is computed by integrating the skin friction coefficient
/// along the airfoil surface:
///
/// ```text
/// Cd_f = ∫ Cf · cos(θ) ds ≈ Σ Cf_avg · Δs
/// ```
///
/// For small surface angles, cos(θ) ≈ 1.
///
/// # Pressure Drag (Squire-Young Formula)
/// Pressure drag is computed from the momentum deficit at the wake trailing edge
/// using the Squire-Young formula:
///
/// ```text
/// Cd_p = 2 · θ_wake · (Ue_wake)^((5+H)/2)
/// ```
///
/// This accounts for the momentum lost due to the boundary layer and wake.
///
/// # Arguments
/// * `stations` - BL stations from solve_viscous()
/// * `config` - Solver configuration (for Reynolds number, etc.)
///
/// # Returns
/// `AeroForces` struct with all coefficients.
///
/// # Note
/// CL and CM are typically taken from the inviscid solution, as viscous
/// corrections to lift are small for attached flow. After merge with flexfoil,
/// this function will accept the inviscid solution to extract CL/CM.
///
/// # XFOIL Reference
/// XFOIL xfoil.f CDCALC (line 1172)
pub fn compute_forces(stations: &[BlStation], _config: &ViscousSolverConfig) -> AeroForces {
    if stations.len() < 2 {
        return AeroForces::default();
    }

    // === Friction Drag ===
    // Integrate Cf along the surface
    let cd_friction: f64 = stations
        .windows(2)
        .map(|pair| {
            let s1 = &pair[0];
            let s2 = &pair[1];

            // Arc length step
            let ds = (s2.x - s1.x).abs();

            // Average skin friction
            let cf_avg = 0.5 * (s1.cf.abs() + s2.cf.abs());

            // For friction drag, we need Cf to be positive
            // Negative Cf indicates separation; use 0 for drag contribution
            let cf_eff = if s1.cf > 0.0 && s2.cf > 0.0 {
                cf_avg
            } else {
                0.0
            };

            cf_eff * ds
        })
        .sum();

    // === Pressure Drag (Squire-Young) ===
    // Find the wake trailing edge (last station or last wake station)
    let wake_station = stations
        .iter()
        .rev()
        .find(|s| s.is_wake)
        .unwrap_or_else(|| stations.last().unwrap());

    // Squire-Young formula for pressure drag
    // Cd_p = 2 * theta * Ue^((5+H)/2)
    let h_wake = wake_station.h.clamp(1.0, 4.0);
    let ue_wake = wake_station.u.abs().max(0.01);
    let theta_wake = wake_station.theta;

    let exponent = (5.0 + h_wake) / 2.0;
    let cd_pressure = 2.0 * theta_wake * ue_wake.powf(exponent);

    // Total drag
    let cd = cd_friction + cd_pressure;

    // CL and CM would come from inviscid solution
    // For now, estimate CL from circulation (placeholder)
    // TODO: After merge, get CL/CM from InviscidSolution
    let cl = 0.0; // Placeholder - set by caller or from inviscid
    let cm = 0.0; // Placeholder

    AeroForces {
        cl,
        cd,
        cm,
        cd_pressure,
        cd_friction,
    }
}

/// Compute friction drag only (useful for debugging).
///
/// # Arguments
/// * `stations` - BL stations
///
/// # Returns
/// Integrated friction drag coefficient.
pub fn compute_friction_drag(stations: &[BlStation]) -> f64 {
    if stations.len() < 2 {
        return 0.0;
    }

    stations
        .windows(2)
        .map(|pair| {
            let ds = (pair[1].x - pair[0].x).abs();
            let cf_avg = 0.5 * (pair[0].cf.max(0.0) + pair[1].cf.max(0.0));
            cf_avg * ds
        })
        .sum()
}

/// Compute pressure drag using Squire-Young formula.
///
/// # Arguments
/// * `theta` - Momentum thickness at wake trailing edge
/// * `ue` - Edge velocity at wake trailing edge
/// * `h` - Shape factor at wake trailing edge
///
/// # Returns
/// Pressure drag coefficient.
pub fn squire_young_drag(theta: f64, ue: f64, h: f64) -> f64 {
    let h_clamped = h.clamp(1.0, 4.0);
    let exponent = (5.0 + h_clamped) / 2.0;
    2.0 * theta * ue.abs().powf(exponent)
}

/// Alternative drag computation using momentum deficit integration.
///
/// This integrates the momentum deficit along the wake, which should
/// give the same result as Squire-Young for a properly resolved wake.
///
/// # Arguments
/// * `wake_stations` - Stations in the wake region
///
/// # Returns
/// Profile drag coefficient from momentum integration.
pub fn compute_wake_drag(wake_stations: &[BlStation]) -> f64 {
    if wake_stations.is_empty() {
        return 0.0;
    }

    // Take momentum deficit at wake end
    let last = wake_stations.last().unwrap();
    2.0 * last.theta
}

/// Estimate transition drag penalty.
///
/// The drag increase due to transition from laminar to turbulent flow
/// can be estimated from the change in skin friction.
///
/// # Arguments
/// * `stations` - BL stations
/// * `x_transition` - Transition location (x/c)
///
/// # Returns
/// Estimated drag penalty due to transition.
pub fn transition_drag_penalty(stations: &[BlStation], x_transition: f64) -> f64 {
    // Find stations around transition
    let laminar: Vec<_> = stations
        .iter()
        .filter(|s| s.x < x_transition && s.is_laminar)
        .collect();
    let turbulent: Vec<_> = stations
        .iter()
        .filter(|s| s.x > x_transition && s.is_turbulent)
        .collect();

    if laminar.is_empty() || turbulent.is_empty() {
        return 0.0;
    }

    // Average Cf in each region
    let cf_lam_avg: f64 = laminar.iter().map(|s| s.cf).sum::<f64>() / laminar.len() as f64;
    let cf_turb_avg: f64 = turbulent.iter().map(|s| s.cf).sum::<f64>() / turbulent.len() as f64;

    // Drag penalty is the extra Cf integrated over the turbulent region
    let x_turb_start = turbulent.first().unwrap().x;
    let x_turb_end = turbulent.last().unwrap().x;
    let delta_cf = cf_turb_avg - cf_lam_avg;

    delta_cf.max(0.0) * (x_turb_end - x_turb_start)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_stations(n: usize) -> Vec<BlStation> {
        (0..n)
            .map(|i| {
                let x = i as f64 / (n - 1) as f64;
                let mut s = BlStation::new();
                s.x = x;
                s.u = 1.0;
                s.theta = 0.001 * (1.0 + x); // Growing theta
                s.delta_star = 2.5 * s.theta;
                s.h = 2.5;
                s.cf = 0.003 * (1.0 - 0.5 * x); // Decreasing Cf
                s.is_laminar = x < 0.3;
                s.is_turbulent = x >= 0.3;
                s.is_wake = x > 0.95;
                s
            })
            .collect()
    }

    #[test]
    fn test_compute_forces() {
        let stations = create_test_stations(50);
        let config = ViscousSolverConfig::default();

        let forces = compute_forces(&stations, &config);

        // Drag should be positive
        assert!(forces.cd > 0.0);
        assert!(forces.cd_friction > 0.0);
        assert!(forces.cd_pressure >= 0.0);

        // Friction + pressure = total
        assert!((forces.cd - forces.cd_friction - forces.cd_pressure).abs() < 1e-10);
    }

    #[test]
    fn test_compute_friction_drag() {
        let stations = create_test_stations(50);

        let cd_f = compute_friction_drag(&stations);

        assert!(cd_f > 0.0);
        assert!(cd_f < 1.0); // Reasonable magnitude
    }

    #[test]
    fn test_squire_young_drag() {
        // Test Squire-Young formula
        let cd = squire_young_drag(0.001, 1.0, 2.5);

        // Cd = 2 * 0.001 * 1.0^(7.5/2) = 0.002
        assert!((cd - 0.002).abs() < 1e-10);
    }

    #[test]
    fn test_squire_young_ue_dependence() {
        // Higher Ue should give higher drag (more momentum)
        let cd_low = squire_young_drag(0.001, 0.5, 2.0);
        let cd_high = squire_young_drag(0.001, 1.0, 2.0);

        assert!(cd_high > cd_low);
    }

    #[test]
    fn test_aero_forces_lift_to_drag() {
        let forces = AeroForces {
            cl: 1.0,
            cd: 0.01,
            cm: -0.1,
            cd_pressure: 0.002,
            cd_friction: 0.008,
        };

        assert!((forces.lift_to_drag() - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_aero_forces_lift_to_drag_zero_drag() {
        let forces = AeroForces {
            cl: 1.0,
            cd: 0.0,
            ..Default::default()
        };

        assert!(forces.lift_to_drag().is_infinite());
    }

    #[test]
    fn test_aero_forces_is_valid() {
        let valid = AeroForces {
            cl: 1.0,
            cd: 0.01,
            cm: -0.1,
            cd_pressure: 0.002,
            cd_friction: 0.008,
        };
        assert!(valid.is_valid());

        let invalid = AeroForces {
            cd: -0.01, // Negative drag is unphysical
            ..Default::default()
        };
        assert!(!invalid.is_valid());
    }

    #[test]
    fn test_transition_drag_penalty() {
        let stations = create_test_stations(50);
        let x_trans = 0.3;

        let penalty = transition_drag_penalty(&stations, x_trans);

        // Penalty should be non-negative
        assert!(penalty >= 0.0);
    }

    #[test]
    fn test_compute_wake_drag() {
        let mut stations = create_test_stations(50);
        // Mark last few as wake
        for s in stations.iter_mut().skip(45) {
            s.is_wake = true;
        }

        let wake: Vec<_> = stations.iter().filter(|s| s.is_wake).cloned().collect();
        let cd = compute_wake_drag(&wake);

        // Should be 2*theta at last station
        let last = wake.last().unwrap();
        assert!((cd - 2.0 * last.theta).abs() < 1e-10);
    }

    #[test]
    fn test_empty_stations() {
        let empty: Vec<BlStation> = vec![];
        let config = ViscousSolverConfig::default();

        let forces = compute_forces(&empty, &config);

        assert_eq!(forces.cd, 0.0);
        assert_eq!(forces.cd_friction, 0.0);
    }
}
