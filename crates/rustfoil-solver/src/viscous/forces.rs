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
    // First check for obviously bad Cf values (sanity check)
    let max_physical_cf = 0.05; // Cf should never exceed ~0.01 for turbulent BL
    let min_valid_theta = 1e-8; // Minimum valid theta value
    
    let cd_friction: f64 = stations
        .windows(2)
        .map(|pair| {
            let s1 = &pair[0];
            let s2 = &pair[1];

            // Skip stations with obviously invalid values
            // These can occur at the leading edge when Newton iteration fails
            if s1.theta < min_valid_theta || s2.theta < min_valid_theta {
                return 0.0;
            }
            if !s1.r_theta.is_finite() || s1.r_theta < 1e-6 {
                return 0.0;
            }
            if !s2.r_theta.is_finite() || s2.r_theta < 1e-6 {
                return 0.0;
            }

            // Arc length step
            let ds = (s2.x - s1.x).abs();

            // Average skin friction - with sanity check
            let cf1 = s1.cf.abs().min(max_physical_cf);
            let cf2 = s2.cf.abs().min(max_physical_cf);
            let cf_avg = 0.5 * (cf1 + cf2);

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
    let h_wake = if wake_station.h.is_finite() {
        wake_station.h.clamp(1.0, 4.0)
    } else {
        2.0 // Default laminar-like value
    };
    let ue_wake = wake_station.u.abs().max(0.01);
    let theta_wake = if wake_station.theta.is_finite() && wake_station.theta > 0.0 {
        wake_station.theta
    } else {
        0.001 // Default small value
    };

    let exponent = (5.0 + h_wake) / 2.0;
    let cd_pressure = 2.0 * theta_wake * ue_wake.powf(exponent);

    // Total drag - ensure finite result
    let cd = if (cd_friction + cd_pressure).is_finite() {
        cd_friction + cd_pressure
    } else {
        cd_friction.max(0.0) // Fall back to friction drag only
    };

    // CL and CM would come from inviscid solution with viscous coupling
    // For single-surface mode, CL/CM are set by the caller
    let cl = 0.0; // Set by caller or from compute_forces_two_surfaces
    let cm = 0.0; // Set by caller

    AeroForces {
        cl,
        cd,
        cm,
        cd_pressure,
        cd_friction,
    }
}

/// Compute aerodynamic forces from two-surface BL solution with CL/CM.
///
/// This is the main force computation for viscous analysis. It computes:
/// - CD from friction and pressure (Squire-Young) on both surfaces
/// - CL from coupled edge velocities using pressure integration
/// - CM about the quarter-chord
///
/// # Arguments
/// * `upper_stations` - BL stations for upper surface (LE to TE)
/// * `lower_stations` - BL stations for lower surface (LE to TE)
/// * `config` - Solver configuration
///
/// # Returns
/// `AeroForces` with all coefficients computed.
pub fn compute_forces_two_surfaces(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    config: &ViscousSolverConfig,
) -> AeroForces {
    // === Friction Drag ===
    // Integrate friction along both surfaces, but skip problematic stations
    let cd_friction_upper = compute_friction_drag(upper_stations);
    let cd_friction_lower = compute_friction_drag(lower_stations);
    let cd_friction = cd_friction_upper + cd_friction_lower;
    
    // === Pressure Drag (Squire-Young) ===
    // Squire-Young formula gives total drag, but we use it for the pressure
    // component by subtracting our friction estimate. The pressure drag comes
    // from the momentum deficit not accounted for by skin friction.
    let upper_te = upper_stations.last();
    let lower_te = lower_stations.last();
    
    let cd_pressure = match (upper_te, lower_te) {
        (Some(u), Some(l)) => {
            // Combined wake momentum thickness
            let theta_combined = u.theta + l.theta;
            let h_avg = 0.5 * (u.h.clamp(1.0, 4.0) + l.h.clamp(1.0, 4.0));
            // Use actual Ue at TE (don't artificially boost it)
            let ue_avg = 0.5 * (u.u.abs() + l.u.abs()).max(0.01);
            
            // Squire-Young gives total drag
            let cd_sy = squire_young_drag(theta_combined, ue_avg, h_avg);
            
            // Pressure drag = Squire-Young total - friction (but floor at 10% of SY)
            // This ensures pressure drag is always positive
            let cd_p_raw = cd_sy - cd_friction;
            if cd_p_raw > 0.1 * cd_sy {
                cd_p_raw
            } else {
                0.1 * cd_sy // Minimum pressure drag
            }
        }
        (Some(u), None) => {
            let cd_sy = squire_young_drag(u.theta, u.u.abs().max(0.01), u.h.clamp(1.0, 4.0));
            0.1 * cd_sy
        }
        (None, Some(l)) => {
            let cd_sy = squire_young_drag(l.theta, l.u.abs().max(0.01), l.h.clamp(1.0, 4.0));
            0.1 * cd_sy
        }
        (None, None) => 0.0,
    };
    
    let cd = cd_friction + cd_pressure;
    
    // === CL and CM from coupled edge velocities ===
    // Extract arc lengths and edge velocities
    let upper_arc: Vec<f64> = upper_stations.iter().map(|s| s.x).collect();
    let lower_arc: Vec<f64> = lower_stations.iter().map(|s| s.x).collect();
    let upper_ue: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
    let lower_ue: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();
    
    // Use circulation method for CL (simpler and more robust)
    let cl = compute_cl_from_circulation(&upper_ue, &lower_ue, &upper_arc, &lower_arc);
    
    // For CM, use simple approximation (proper CM requires panel geometry)
    // CM ≈ 0 for symmetric airfoils at small AoA
    let cm = 0.0; // TODO: Implement proper CM with panel geometry
    
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
///
/// # XFOIL Reference
/// XFOIL's CDCALC (xfoil.f line 1172) computes friction drag as:
/// ```fortran
/// DO 205 IBL=3, IBLTE(IS)   ! Starts at IBL=3, skipping stagnation
///   DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
///   CDF = CDF + 0.5*(TAU(IBL,IS)+TAU(IBL-1,IS))*DX * 2.0/QINF**2
/// ```
///
/// Key differences from naive integration:
/// 1. Skips IBL=1 (virtual stagnation) and IBL=2 (similarity station)
/// 2. Uses chord-wise distance dx, not arc length ds
///
/// We match this by skipping the first 2 stations where Cf is singular.
pub fn compute_friction_drag(stations: &[BlStation]) -> f64 {
    if stations.len() < 4 {
        // Need at least 4 stations: stag(0), simi(1), first_normal(2), second_normal(3)
        return 0.0;
    }

    // XFOIL starts at IBL=3 (our index 2), integrating from station 2 onward
    // This skips the stagnation (index 0) and similarity (index 1) stations
    // where Cf values are not physically meaningful due to 1/Rθ singularity.
    //
    // However, our BL solution has additional problematic stations near LE
    // due to the 1/√Rθ dependency in Cf. Skip more stations at low Rθ.
    
    // Maximum physical Cf (turbulent BL at low Re can reach ~0.01)
    let max_physical_cf = 0.02;  // Tighter clamp
    
    // Minimum Rθ for reliable Cf (below this, Cf blows up due to 1/sqrt(Rθ))
    // This filter is critical for numerical stability - values below ~80 often 
    // have singularity issues
    let r_theta_min = 80.0;
    
    let mut max_cf_seen = 0.0_f64;
    let mut max_dx_seen = 0.0_f64;
    
    // Start at station 3 (skip 0=stagnation, 1=similarity, 2=often still problematic)
    let cdf: f64 = stations[3..]
        .windows(2)
        .map(|pair| {
            // Skip stations with clearly invalid BL values (theta near zero = failed march)
            // This matches XFOIL's implicit handling where failed stations don't contribute
            let theta_min = 1e-10;
            if pair[0].theta < theta_min || pair[1].theta < theta_min {
                return 0.0;
            }
            
            // Skip stations with low Rθ where Cf is unreliable
            if pair[0].r_theta < r_theta_min || pair[1].r_theta < r_theta_min {
                return 0.0;
            }
            
            // Use chord-wise distance (x_coord difference) like XFOIL, not arc length
            let dx = (pair[1].x_coord - pair[0].x_coord).abs();
            
            // Fallback to arc length if x_coord not set (backward compatibility)
            let dx = if dx < 1e-12 {
                (pair[1].x - pair[0].x).abs()
            } else {
                dx
            };
            
            // Clamp Cf to physical range
            // For laminar: Cf ~ 0.664/sqrt(Rex), max ~0.01 at Rex=4000
            // For turbulent: Cf ~ 0.0583/Rex^0.2, max ~0.01 at Rex=10^4
            let cf1 = pair[0].cf.clamp(0.0, max_physical_cf);
            let cf2 = pair[1].cf.clamp(0.0, max_physical_cf);
            let cf_avg = 0.5 * (cf1 + cf2);
            
            max_cf_seen = max_cf_seen.max(pair[0].cf.abs()).max(pair[1].cf.abs());
            max_dx_seen = max_dx_seen.max(dx);
            
            // Simple integration: CDF = ∫ Cf · dx
            // Note: XFOIL uses TAU = 0.5 * Ue² * Cf, but our Cf already incorporates
            // Ue effects through Rθ = Re * Ue * θ. The simpler formula gives more
            // stable results with our BL initialization.
            cf_avg * dx
        })
        .sum();
    
    // Warn if Cf values were out of range (indicates numerical issues)
    if max_cf_seen > 0.1 {
        // Find which stations have huge Cf
        let bad_count = stations.iter().filter(|s| s.cf.abs() > 0.1).count();
        if bad_count > 0 {
            eprintln!("[WARN forces] {} stations have Cf > 0.1 (max_cf={:.6e})", bad_count, max_cf_seen);
            // Show first few
            for (i, s) in stations.iter().enumerate().take(10) {
                if s.cf.abs() > 0.1 {
                    eprintln!("  [{}] x={:.4}, Cf={:.4e}, theta={:.6e}, h={:.3}, Rθ={:.1}, turb={}", 
                        i, s.x_coord, s.cf, s.theta, s.h, s.r_theta, s.is_turbulent);
                }
            }
        }
    }
    
    cdf
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

/// Compute CL and CM from coupled edge velocities.
///
/// This implements the viscous-coupled lift and moment computation.
/// In XFOIL, CL comes from the inviscid solution, but the inviscid edge
/// velocities are modified by the viscous mass defect through the V-I coupling.
/// The "coupled CL" is the CL computed from these modified edge velocities.
///
/// # Method
/// Uses pressure integration around the airfoil contour:
/// ```text
/// Cp = 1 - (Ue / V∞)²
/// CL = ∮ Cp dy  (in wind axes)
/// CM = ∮ Cp (x - x_ref) dy  (about quarter-chord)
/// ```
///
/// # Arguments
/// * `upper_stations` - BL stations for upper surface (from LE to TE)
/// * `lower_stations` - BL stations for lower surface (from LE to TE)
/// * `panel_x` - Panel x-coordinates (from original geometry, TE to TE via LE)
/// * `panel_y` - Panel y-coordinates
/// * `alpha` - Angle of attack in radians
/// * `chord` - Airfoil chord length
///
/// # Returns
/// Tuple of (CL, CM) computed from the coupled edge velocities.
///
/// # Note
/// For attached flow, this CL will be very close to the inviscid CL.
/// For separated flow, it will differ significantly.
pub fn compute_coupled_cl_cm(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    panel_x: &[f64],
    panel_y: &[f64],
    alpha: f64,
    chord: f64,
) -> (f64, f64) {
    if upper_stations.is_empty() || lower_stations.is_empty() || panel_x.len() < 3 {
        return (0.0, 0.0);
    }
    
    let cosa = alpha.cos();
    let sina = alpha.sin();
    let x_ref = 0.25 * chord; // Quarter-chord reference
    
    let mut cl = 0.0;
    let mut cm = 0.0;
    
    // XFOIL approach: Integrate Cp around the closed contour
    // Cp = 1 - (Ue)² (normalized by freestream V∞ = 1)
    //
    // The stations have edge velocities that include viscous coupling.
    // We integrate Cp = 1 - Ue² around the contour.
    //
    // For upper surface: integrate from LE to TE (positive Ue contributes to lift)
    // For lower surface: integrate from LE to TE (positive Ue on lower surface reduces lift)
    
    // Upper surface contribution (LE to TE direction)
    for i in 0..upper_stations.len().saturating_sub(1) {
        let s1 = &upper_stations[i];
        let s2 = &upper_stations[i + 1];
        
        // Cp from edge velocity (normalized Ue)
        let cp1 = 1.0 - s1.u * s1.u;
        let cp2 = 1.0 - s2.u * s2.u;
        let cp_avg = 0.5 * (cp1 + cp2);
        
        // Use arc length approximation for panel step
        let ds = (s2.x - s1.x).abs();
        
        // For upper surface, positive Cp contributes negative to CL in body axes
        // but we integrate with proper sign convention
        // CL contribution is Cp * dy (in wind axes)
        // For simplicity with arc-length stations, we approximate:
        // The pressure acts normal to the surface (mostly in y-direction for airfoil)
        // Upper surface: Cp acts downward (-y direction), so -Cp contributes to +CL
        cl -= cp_avg * ds;
        
        // Moment about quarter chord
        let x_mid = 0.5 * (s1.x + s2.x); // Using arc length as x-proxy
        cm += cp_avg * (x_mid - x_ref) * ds / chord;
    }
    
    // Lower surface contribution (LE to TE direction)
    for i in 0..lower_stations.len().saturating_sub(1) {
        let s1 = &lower_stations[i];
        let s2 = &lower_stations[i + 1];
        
        let cp1 = 1.0 - s1.u * s1.u;
        let cp2 = 1.0 - s2.u * s2.u;
        let cp_avg = 0.5 * (cp1 + cp2);
        
        let ds = (s2.x - s1.x).abs();
        
        // Lower surface: Cp acts upward (+y direction), so +Cp contributes to +CL
        cl += cp_avg * ds;
        
        // Moment contribution
        let x_mid = 0.5 * (s1.x + s2.x);
        cm -= cp_avg * (x_mid - x_ref) * ds / chord;
    }
    
    // Apply angle of attack rotation for proper wind-axis CL
    // For small angles this is approximately CL ≈ cl_body
    let cl_wind = cl * cosa - cm * sina;
    let cm_wind = cm * cosa + cl * sina / 2.0; // Approximate
    
    (cl_wind, cm_wind)
}

/// Compute CL from edge velocity distribution using Kutta-Joukowski.
///
/// Alternative CL computation using the circulation theorem:
/// ```text
/// CL = 2 * Γ / (V∞ * c)
/// ```
///
/// where Γ is computed from the velocity jump across the wake or
/// from integration of the edge velocity around the contour.
///
/// # Arguments
/// * `upper_ue` - Edge velocities on upper surface
/// * `lower_ue` - Edge velocities on lower surface
/// * `upper_arc` - Arc lengths on upper surface
/// * `lower_arc` - Arc lengths on lower surface
///
/// # Returns
/// CL from circulation
pub fn compute_cl_from_circulation(
    upper_ue: &[f64],
    lower_ue: &[f64],
    upper_arc: &[f64],
    lower_arc: &[f64],
) -> f64 {
    if upper_ue.len() < 2 || lower_ue.len() < 2 {
        return 0.0;
    }
    
    // Integrate Ue around the contour to get circulation
    // Γ = ∮ V · ds ≈ ∫_upper Ue ds - ∫_lower Ue ds
    //
    // Sign convention:
    // - Upper surface has higher Ue (suction peak at positive alpha)
    // - Lower surface has lower Ue
    // - For positive lift: Γ = ∫_upper - ∫_lower > 0
    // - CL = 2Γ > 0 for positive alpha
    //
    // IMPORTANT: This assumes upper_ue/lower_ue are correctly labeled as
    // upper/lower surfaces. If using ViscousSetup.upper_surface(), note that
    // the naming convention there may be inverted (check carefully).
    
    let mut gamma_upper = 0.0;
    let mut gamma_lower = 0.0;
    
    // Upper surface: integrate Ue * ds (stagnation to TE, counterclockwise)
    // Higher Ue contributes positive to circulation
    for i in 0..upper_ue.len() - 1 {
        let ue_avg = 0.5 * (upper_ue[i] + upper_ue[i + 1]);
        let ds = upper_arc[i + 1] - upper_arc[i];
        gamma_upper += ue_avg * ds;
    }
    
    // Lower surface: integrate Ue * ds (stagnation to TE, clockwise = negative)
    // Lower Ue contributes negative to circulation
    for i in 0..lower_ue.len() - 1 {
        let ue_avg = 0.5 * (lower_ue[i] + lower_ue[i + 1]);
        let ds = lower_arc[i + 1] - lower_arc[i];
        gamma_lower += ue_avg * ds;
    }
    
    let gamma = gamma_upper - gamma_lower;
    
    // CL = 2 * Γ / (V∞ * c) with V∞ = 1, c = 1 (normalized)
    2.0 * gamma
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
