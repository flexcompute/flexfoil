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
use rustfoil_bl::{add_event, is_debug_active, CdBreakdownEvent, ClDetailEvent, DebugEvent};

use super::config::ViscousSolverConfig;
use super::state::{CanonicalBlRow, XfoilLikeViscousState};

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

/// Compute XFOIL-style force coefficients from the canonical panel gamma array.
///
/// This matches the inviscid `CLCALC` integration path and returns both lift and
/// pitching moment about quarter chord.
pub fn compute_panel_forces_from_gamma(
    panel_x: &[f64],
    panel_y: &[f64],
    gamma: &[f64],
    alpha_rad: f64,
) -> (f64, f64) {
    let n = panel_x.len();
    if n < 3 || panel_y.len() != n || gamma.len() != n {
        return (0.0, 0.0);
    }

    let cosa = alpha_rad.cos();
    let sina = alpha_rad.sin();

    let x_min = panel_x.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = panel_x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let chord = (x_max - x_min).abs().max(1.0e-12);
    let x_ref = x_min + 0.25 * chord;

    let mut cl = 0.0;
    let mut cm = 0.0;

    for i in 0..n {
        let ip = (i + 1) % n;

        let dx = panel_x[ip] - panel_x[i];
        let dy = panel_y[ip] - panel_y[i];
        let dx_wind = dx * cosa + dy * sina;
        let dy_wind = dy * cosa - dx * sina;

        let cp_i = 1.0 - gamma[i] * gamma[i];
        let cp_ip = 1.0 - gamma[ip] * gamma[ip];
        let cp_avg = 0.5 * (cp_i + cp_ip);
        let dcp = cp_ip - cp_i;

        cl += cp_avg * dx_wind;

        let x_mid = 0.5 * (panel_x[i] + panel_x[ip]);
        let y_mid = 0.5 * (panel_y[i] + panel_y[ip]);
        let ax = (x_mid - x_ref) * cosa + y_mid * sina;
        let ay = y_mid * cosa - (x_mid - x_ref) * sina;

        cm -= cp_avg * (ax * dx_wind / chord + ay * dy_wind / chord);
        cm -= dcp * dx_wind * dx_wind / (12.0 * chord);
        cm -= dcp * dy_wind * dy_wind / (12.0 * chord);
    }

    (cl, cm)
}

/// Compute aerodynamic forces from the canonical viscous state.
pub fn compute_forces_from_canonical_state(
    state: &XfoilLikeViscousState,
    panel_x: &[f64],
    panel_y: &[f64],
    alpha_rad: f64,
    _config: &ViscousSolverConfig,
) -> AeroForces {
    let (cl, cm) = compute_panel_forces_from_gamma(panel_x, panel_y, state.panel_gamma(), alpha_rad);

    let cd_friction_upper = compute_row_friction_drag(&state.upper_rows);
    let cd_friction_lower = compute_row_friction_drag(&state.lower_rows);
    let cd_friction = cd_friction_upper + cd_friction_lower;

    let lower_wake: Vec<&CanonicalBlRow> = state.lower_rows.iter().filter(|row| row.is_wake).collect();
    let te_extrapolated_cd = compute_cd_from_te_rows(&state.upper_rows, &state.lower_rows);
    let cd_total = if !lower_wake.is_empty() {
        let cd_wake = compute_cd_from_wake_rows(&lower_wake);
        if cd_wake.is_finite() && cd_wake > 0.0 {
            cd_wake
        } else {
            te_extrapolated_cd.max(1.0e-6)
        }
    } else {
        te_extrapolated_cd.max(1.0e-6)
    };
    // Match XFOIL CDCALC: the wake/TE total drag can be lower than the
    // integrated skin-friction contribution, which yields a negative form-drag
    // component at low alpha. Do not floor total CD to CDf here.
    let cd = cd_total;
    let cd_pressure = cd - cd_friction;

    if is_debug_active() {
        let upper_te = state.upper_rows.iter().rev().find(|row| !row.is_wake);
        let lower_te = state.lower_rows.iter().rev().find(|row| !row.is_wake);
        let (theta_te_upper, h_te_upper, ue_te_upper) = upper_te
            .map(|row| (row.theta, row.h, row.uedg))
            .unwrap_or((0.0, 0.0, 0.0));
        let (theta_te_lower, h_te_lower, ue_te_lower) = lower_te
            .map(|row| (row.theta, row.h, row.uedg))
            .unwrap_or((0.0, 0.0, 0.0));
        let cd_pressure_upper = if ue_te_upper.abs() > 0.01 {
            squire_young_drag(theta_te_upper, ue_te_upper.abs(), h_te_upper.clamp(1.0, 4.0))
        } else {
            0.0
        };
        let cd_pressure_lower = if ue_te_lower.abs() > 0.01 {
            squire_young_drag(theta_te_lower, ue_te_lower.abs(), h_te_lower.clamp(1.0, 4.0))
        } else {
            0.0
        };

        add_event(DebugEvent::cd_breakdown(
            0,
            CdBreakdownEvent {
                cd_friction,
                cd_pressure,
                cd_total: cd,
                cd_friction_upper,
                cd_friction_lower,
                cd_pressure_upper,
                cd_pressure_lower,
                theta_te_upper,
                theta_te_lower,
                h_te_upper,
                h_te_lower,
                ue_te_upper,
                ue_te_lower,
            },
        ));
        add_event(DebugEvent::cl_detail(
            0,
            ClDetailEvent {
                cl,
                cm,
                cdp: cd_pressure,
                alpha_rad,
                alpha_deg: alpha_rad.to_degrees(),
            },
        ));
    }

    AeroForces {
        cl,
        cd,
        cm,
        cd_pressure,
        cd_friction,
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
    _config: &ViscousSolverConfig,
) -> AeroForces {
    // === Friction Drag ===
    // Integrate friction along both surfaces, but skip problematic stations
    let cd_friction_upper = compute_friction_drag(upper_stations);
    let cd_friction_lower = compute_friction_drag(lower_stations);
    let cd_friction = cd_friction_upper + cd_friction_lower;
    
    // === Total Drag (prefer converged wake, else TE fallback) ===
    // When the lower surface already carries a converged wake, use that directly
    // for the Squire-Young total drag path. Fall back to the older TE-based
    // extrapolation only when no wake stations are available.
    let upper_te = upper_stations.iter().rev().find(|s| !s.is_wake);
    let lower_te = lower_stations.iter().rev().find(|s| !s.is_wake);
    let lower_wake: Vec<BlStation> = lower_stations
        .iter()
        .filter(|s| s.is_wake)
        .cloned()
        .collect();

    let te_extrapolated_cd = match (upper_te, lower_te) {
        (Some(u), Some(l)) => {
            let theta_te = u.theta + l.theta;
            let ue_te = 0.5 * (u.u.abs() + l.u.abs()).max(0.1);
            let ue_far_wake = 0.99;
            let ue_ratio_sq = (ue_te / ue_far_wake).powi(2);
            let theta_far_wake = theta_te * ue_ratio_sq;
            squire_young_drag(theta_far_wake, ue_far_wake, 1.04)
        }
        (Some(u), None) => squire_young_drag(u.theta, 0.99, 1.04),
        (None, Some(l)) => squire_young_drag(l.theta, 0.99, 1.04),
        (None, None) => 0.0,
    };

    let cd_total = if !lower_wake.is_empty() {
        let cd_wake = compute_cd_from_wake(&lower_wake, cd_friction);
        if cd_wake.is_finite() && cd_wake > 0.0 {
            cd_wake
        } else {
            te_extrapolated_cd.max(cd_friction).max(1.0e-6)
        }
    } else {
        te_extrapolated_cd
    };
    
    // CRITICAL FIX: Ensure cd_total >= cd_friction
    // Physically, total drag must be at least friction drag (cd_total = cd_friction + cd_pressure)
    // If Squire-Young gives cd_total < cd_friction, it means theta_te was too small or the
    // extrapolation underestimated the drag. In this case, use cd_friction as a lower bound.
    // This prevents cd_pressure from becoming negative (which would cause cd = 0.0).
    let cd_total = cd_total.max(cd_friction);
    
    // Pressure drag = Total - Friction
    // (This is how XFOIL separates the components)
    // Note: XFOIL's CL_DETAIL reports 'cdp' which is pressure drag only
    let cd_pressure = cd_total - cd_friction;
    
    // Report total drag (cd_total) as the primary CD output
    // XFOIL's CD_BREAKDOWN reports cd_total, and that's what comparisons typically use
    // Pressure drag (cd_pressure) is also stored separately for detailed analysis
    let cd = cd_total;
    
    // Debug drag computation
    if std::env::var("RUSTFOIL_DRAG_DEBUG").is_ok() {
        if let (Some(u), Some(l)) = (upper_te, lower_te) {
            let theta_te = u.theta + l.theta;
            let ue_te = 0.5 * (u.u.abs() + l.u.abs());
            let ue_ratio_sq = (ue_te / 0.99).powi(2);
            let theta_far_wake = theta_te * ue_ratio_sq;
            eprintln!("[DRAG_DEBUG] Upper TE: x={:.4} θ={:.4e} H={:.3} Ue={:.4}", u.x, u.theta, u.h, u.u);
            eprintln!("[DRAG_DEBUG] Lower TE: x={:.4} θ={:.4e} H={:.3} Ue={:.4}", l.x, l.theta, l.h, l.u);
            eprintln!("[DRAG_DEBUG] Combined TE: θ={:.4e} Ue={:.4}", theta_te, ue_te);
            eprintln!("[DRAG_DEBUG] Far wake (θ*Ue²=const): θ={:.4e}", theta_far_wake);
            eprintln!("[DRAG_DEBUG] CD_f={:.5} (upper={:.5} lower={:.5})", cd_friction, cd_friction_upper, cd_friction_lower);
            eprintln!("[DRAG_DEBUG] CD_total={:.5} CD_p={:.5}", cd_total, cd_pressure);
        }
    }
    
    // === CL and CM from coupled edge velocities ===
    // Extract x-coordinates (for XFOIL-style integration) and edge velocities
    // CRITICAL: Use x_coord (actual panel x-coordinate), NOT x (arc length)!
    // Near the curved leading edge, arc length >> x-coordinate, causing CL inflation.
    // CRITICAL: Exclude wake stations! CL is only from airfoil surface pressure.
    let upper_x: Vec<f64> = upper_stations.iter().filter(|s| !s.is_wake).map(|s| s.x_coord).collect();
    let lower_x: Vec<f64> = lower_stations.iter().filter(|s| !s.is_wake).map(|s| s.x_coord).collect();
    let upper_ue: Vec<f64> = upper_stations.iter().filter(|s| !s.is_wake).map(|s| s.u).collect();
    let lower_ue: Vec<f64> = lower_stations.iter().filter(|s| !s.is_wake).map(|s| s.u).collect();
    
    // Use pressure integration method for CL (matches XFOIL's CLCALC)
    let cl = compute_cl_from_circulation(&upper_ue, &lower_ue, &upper_x, &lower_x);
    
    // For CM, use simple approximation (proper CM requires panel geometry)
    // CM ≈ 0 for symmetric airfoils at small AoA
    let cm = 0.0; // TODO: Implement proper CM with panel geometry
    
    // === Emit debug events for force comparison with XFOIL ===
    if is_debug_active() {
        // Extract TE quantities for breakdown
        let (theta_te_upper, h_te_upper, ue_te_upper) = upper_te
            .map(|s| (s.theta, s.h, s.u))
            .unwrap_or((0.0, 0.0, 0.0));
        let (theta_te_lower, h_te_lower, ue_te_lower) = lower_te
            .map(|s| (s.theta, s.h, s.u))
            .unwrap_or((0.0, 0.0, 0.0));
        
        // Compute per-surface pressure drag (Squire-Young at TE)
        let cd_pressure_upper = if ue_te_upper.abs() > 0.01 {
            squire_young_drag(theta_te_upper, ue_te_upper.abs(), h_te_upper.clamp(1.0, 4.0))
        } else {
            0.0
        };
        let cd_pressure_lower = if ue_te_lower.abs() > 0.01 {
            squire_young_drag(theta_te_lower, ue_te_lower.abs(), h_te_lower.clamp(1.0, 4.0))
        } else {
            0.0
        };
        
        // Emit CD breakdown event
        add_event(DebugEvent::cd_breakdown(
            0, // iteration would be passed from caller if needed
            CdBreakdownEvent {
                cd_friction,
                cd_pressure,
                cd_total: cd,
                cd_friction_upper,
                cd_friction_lower,
                cd_pressure_upper,
                cd_pressure_lower,
                theta_te_upper,
                theta_te_lower,
                h_te_upper,
                h_te_lower,
                ue_te_upper,
                ue_te_lower,
            },
        ));
        
        // Emit CL detail event
        add_event(DebugEvent::cl_detail(
            0, // iteration would be passed from caller if needed
            ClDetailEvent {
                cl,
                cm,
                cdp: cd_pressure,
                alpha_rad: 0.0, // Would need to be passed from caller
                alpha_deg: 0.0,
            },
        ));
    }
    
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
            // Skip wake stations - they have no wall so Cf = 0
            // XFOIL integrates friction only from IBL=3 to IBLTE
            if pair[0].is_wake || pair[1].is_wake {
                return 0.0;
            }
            
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
            
            // XFOIL uses TAU = 0.5 * ρ * Ue² * Cf for wall shear stress
            // Friction drag coefficient: CDF = ∫ TAU * dx * 2/(ρ * V∞²)
            // With normalized ρ=1, V∞=1: CDF = ∫ Ue² * Cf * dx
            let ue1 = pair[0].u.abs().max(0.01);
            let ue2 = pair[1].u.abs().max(0.01);
            let ue_sq_avg = 0.5 * (ue1 * ue1 + ue2 * ue2);
            
            max_cf_seen = max_cf_seen.max(pair[0].cf.abs()).max(pair[1].cf.abs());
            max_dx_seen = max_dx_seen.max(dx);
            
            // Correct integration: CDF = ∫ Ue² · Cf · dx
            ue_sq_avg * cf_avg * dx
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

fn compute_row_friction_drag(rows: &[CanonicalBlRow]) -> f64 {
    if rows.len() < 4 {
        return 0.0;
    }

    let max_physical_cf = 0.02;
    let r_theta_min = 80.0;

    rows[3..]
        .windows(2)
        .map(|pair| {
            if pair[0].is_wake || pair[1].is_wake {
                return 0.0;
            }
            if pair[0].theta < 1.0e-10 || pair[1].theta < 1.0e-10 {
                return 0.0;
            }
            if pair[0].r_theta < r_theta_min || pair[1].r_theta < r_theta_min {
                return 0.0;
            }

            let dx = (pair[1].x_coord - pair[0].x_coord).abs();
            let dx = if dx < 1.0e-12 {
                (pair[1].x - pair[0].x).abs()
            } else {
                dx
            };

            let cf1 = pair[0].cf.clamp(0.0, max_physical_cf);
            let cf2 = pair[1].cf.clamp(0.0, max_physical_cf);
            let cf_avg = 0.5 * (cf1 + cf2);
            let ue1 = pair[0].uedg.abs().max(0.01);
            let ue2 = pair[1].uedg.abs().max(0.01);
            let ue_sq_avg = 0.5 * (ue1 * ue1 + ue2 * ue2);

            ue_sq_avg * cf_avg * dx
        })
        .sum()
}

fn compute_cd_from_te_rows(upper_rows: &[CanonicalBlRow], lower_rows: &[CanonicalBlRow]) -> f64 {
    let upper_te = upper_rows.iter().rev().find(|row| !row.is_wake);
    let lower_te = lower_rows.iter().rev().find(|row| !row.is_wake);

    match (upper_te, lower_te) {
        (Some(u), Some(l)) => {
            let theta_te = u.theta + l.theta;
            let ue_te = 0.5 * (u.uedg.abs() + l.uedg.abs()).max(0.1);
            let ue_far_wake = 0.99;
            let theta_far_wake = theta_te * (ue_te / ue_far_wake).powi(2);
            squire_young_drag(theta_far_wake, ue_far_wake, 1.04)
        }
        (Some(u), None) => squire_young_drag(u.theta, 0.99, 1.04),
        (None, Some(l)) => squire_young_drag(l.theta, 0.99, 1.04),
        (None, None) => 0.0,
    }
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
    if !theta.is_finite() || !ue.is_finite() || !h.is_finite() || theta <= 0.0 {
        return 0.0;
    }
    let exponent = (5.0 + h) / 2.0;
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

/// Compute total drag using far-wake Squire-Young formula.
///
/// XFOIL computes total CD using Squire-Young at the far wake where:
/// - Ue ≈ 1.0 (approaches freestream)
/// - Hk ≈ 1.0-1.1 (thin wake limit)
/// - θ comes from momentum conservation through the wake
///
/// The Squire-Young formula is: CD = 2 * θ * Ue^((5+H)/2)
///
/// At the TE, H≈2.5 gives exponent≈3.75, which overestimates drag by ~50%
/// compared to using far-wake values where H≈1.04 gives exponent≈3.02.
///
/// # Arguments
/// * `wake_stations` - Marched wake stations from TE downstream
/// * `cd_friction` - Already-computed friction drag coefficient
///
/// # Returns
/// Total drag coefficient (CD) from Squire-Young at far wake
pub fn compute_cd_from_wake(wake_stations: &[BlStation], cd_friction: f64) -> f64 {
    if wake_stations.is_empty() {
        return cd_friction;
    }

    let far_wake = wake_stations.last().unwrap();
    let theta_wake = far_wake.theta;
    let ue_wake = far_wake.u.abs();
    // XFOIL CDCALC uses SHWAKE = DSTR/THET at the last lower wake row,
    // where DSTR is the total wake displacement thickness.
    let h_wake = if far_wake.theta > 1e-12 {
        (far_wake.delta_star + far_wake.dw.max(0.0)) / far_wake.theta
    } else {
        0.0
    };

    let mut cd_total = squire_young_drag(theta_wake, ue_wake, h_wake);
    cd_total = cd_total.max(cd_friction);
    
    if std::env::var("RUSTFOIL_DRAG_DEBUG").is_ok() {
        eprintln!("[WAKE_CD_DEBUG] Far wake: x={:.4} θ={:.4e} H={:.3} Ue={:.4}", 
            far_wake.x, theta_wake, h_wake, ue_wake);
        eprintln!("[WAKE_CD_DEBUG] CD_total={:.5} (from S-Y at far wake, clamped to >= cd_friction={:.5})", cd_total, cd_friction);
    }
    
    cd_total
}

fn compute_cd_from_wake_rows(wake_rows: &[&CanonicalBlRow]) -> f64 {
    if wake_rows.is_empty() {
        return 0.0;
    }

    let far_wake = *wake_rows.last().expect("non-empty wake rows");

    let theta_wake = far_wake.theta;
    let ue_wake = far_wake.uedg.abs();
    let h_wake = if far_wake.theta > 1.0e-12 {
        (far_wake.dstr + far_wake.dw.max(0.0)) / far_wake.theta
    } else {
        0.0
    };

    squire_young_drag(theta_wake, ue_wake, h_wake)
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
    _panel_y: &[f64],
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

/// Compute CL from edge velocity distribution using pressure integration.
///
/// This matches XFOIL's CLCALC method (xfoil.f line ~1172):
/// ```text
/// CL = ∮ Cp · dx  (counterclockwise contour integral)
/// ```
///
/// where Cp = 1 - Ue² (incompressible pressure coefficient).
///
/// For two separate surfaces stored LE→TE:
/// - Upper surface contributes: -∫ Cp dx (negative sign for counterclockwise)
/// - Lower surface contributes: +∫ Cp dx
///
/// # Arguments
/// * `upper_ue` - Edge velocities on upper surface (LE to TE)
/// * `lower_ue` - Edge velocities on lower surface (LE to TE)
/// * `upper_x` - X-coordinates on upper surface (MUST be actual x, not arc length!)
/// * `lower_x` - X-coordinates on lower surface (MUST be actual x, not arc length!)
///
/// # Returns
/// CL from pressure integration (matches XFOIL's method)
///
/// # Note
/// Using arc length instead of x-coordinates will inflate CL by ~10% due to
/// the curved leading edge where arc >> x.
pub fn compute_cl_from_circulation(
    upper_ue: &[f64],
    lower_ue: &[f64],
    upper_x: &[f64],
    lower_x: &[f64],
) -> f64 {
    if upper_ue.len() < 2 || lower_ue.len() < 2 {
        return 0.0;
    }
    
    // XFOIL's CLCALC integrates Cp around the closed contour counterclockwise:
    //   CL = ∮ Cp · dx
    //
    // For incompressible flow: Cp = 1 - (Ue/V∞)² = 1 - Ue² (normalized V∞=1)
    //
    // With surfaces stored LE→TE:
    // - Upper surface (counterclockwise = TE→LE): reverse sign → -∫ Cp dx
    // - Lower surface (counterclockwise = LE→TE): keep sign → +∫ Cp dx
    //
    // Physical interpretation:
    // - Upper surface suction (Cp < 0) with negative sign → positive lift
    // - Lower surface pressure (Cp > 0 or less negative) → positive lift
    
    // Upper surface contribution (negative sign for counterclockwise direction)
    let cl_upper: f64 = upper_ue
        .windows(2)
        .zip(upper_x.windows(2))
        .map(|(ue_pair, x_pair)| {
            // Average Cp over panel
            let cp1 = 1.0 - ue_pair[0] * ue_pair[0];
            let cp2 = 1.0 - ue_pair[1] * ue_pair[1];
            let cp_avg = 0.5 * (cp1 + cp2);
            // Chord-wise step (dx)
            let dx = x_pair[1] - x_pair[0];
            // Negative sign: counterclockwise on upper surface is TE→LE
            -cp_avg * dx
        })
        .sum();
    
    // Lower surface contribution (positive sign for counterclockwise direction)
    let cl_lower: f64 = lower_ue
        .windows(2)
        .zip(lower_x.windows(2))
        .map(|(ue_pair, x_pair)| {
            // Average Cp over panel
            let cp1 = 1.0 - ue_pair[0] * ue_pair[0];
            let cp2 = 1.0 - ue_pair[1] * ue_pair[1];
            let cp_avg = 0.5 * (cp1 + cp2);
            // Chord-wise step (dx)
            let dx = x_pair[1] - x_pair[0];
            // Positive sign: counterclockwise on lower surface is LE→TE
            cp_avg * dx
        })
        .sum();
    
    cl_upper + cl_lower
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
                s.x_coord = x;
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
    fn test_compute_forces_two_surfaces_prefers_converged_wake() {
        let upper = create_test_stations(8);
        let mut lower = create_test_stations(8);
        for station in lower.iter_mut().skip(6) {
            station.is_wake = true;
            station.is_turbulent = true;
            station.is_laminar = false;
            station.u = 0.98;
            station.theta = 0.0025;
            station.delta_star = 0.0028;
            station.h = station.delta_star / station.theta;
            station.cf = 0.0;
        }

        let config = ViscousSolverConfig::default();
        let forces = compute_forces_two_surfaces(&upper, &lower, &config);
        let lower_wake: Vec<BlStation> = lower.iter().filter(|s| s.is_wake).cloned().collect();
        let expected_cd = compute_cd_from_wake(&lower_wake, forces.cd_friction);

        assert!((forces.cd - expected_cd).abs() < 1e-10);
        assert!(forces.cd_pressure >= 0.0);
    }

    #[test]
    fn test_compute_forces_from_canonical_state_allows_negative_pressure_drag() {
        let mut upper = create_test_stations(8);
        let mut lower = create_test_stations(8);

        for station in upper.iter_mut().chain(lower.iter_mut()) {
            station.cf = 0.01;
            station.u = 1.0;
            station.is_wake = false;
            station.is_laminar = true;
            station.is_turbulent = false;
        }

        for station in lower.iter_mut().skip(6) {
            station.is_wake = true;
            station.is_laminar = false;
            station.is_turbulent = true;
            station.cf = 0.0;
            station.theta = 0.001;
            station.delta_star = 0.00105;
            station.h = station.delta_star / station.theta;
            station.u = 0.99;
        }

        let state = XfoilLikeViscousState::from_station_views(&upper, &lower, 4);
        let forces = compute_forces_from_canonical_state(
            &state,
            &[0.0, 1.0, 1.0, 0.0],
            &[0.0, 0.0, 0.1, 0.1],
            0.0,
            &ViscousSolverConfig::default(),
        );

        assert!(forces.cd_friction > 0.0);
        assert!(forces.cd > 0.0);
        assert!(forces.cd < forces.cd_friction);
        assert!(forces.cd_pressure < 0.0);
    }

    #[test]
    fn test_compute_forces_from_canonical_state_keeps_large_valid_wake_drag() {
        let mut upper = create_test_stations(8);
        let mut lower = create_test_stations(8);

        for station in upper.iter_mut().chain(lower.iter_mut()) {
            station.u = 1.0;
            station.cf = 0.001;
            station.is_wake = false;
        }

        // Keep the airfoil TE estimate modest.
        if let Some(upper_te) = upper.get_mut(7) {
            upper_te.theta = 0.002;
            upper_te.delta_star = 0.004;
            upper_te.h = upper_te.delta_star / upper_te.theta;
        }
        if let Some(lower_te) = lower.get_mut(5) {
            lower_te.theta = 0.0015;
            lower_te.delta_star = 0.003;
            lower_te.h = lower_te.delta_star / lower_te.theta;
        }

        // Make the far wake materially thicker than the TE estimate so the
        // regression would fail if we reintroduce the old clipping heuristic.
        for (idx, station) in lower.iter_mut().enumerate().skip(6) {
            station.is_wake = true;
            station.is_laminar = false;
            station.is_turbulent = true;
            station.cf = 0.0;
            station.theta = if idx == 6 { 0.02 } else { 0.03 };
            station.delta_star = station.theta * 1.15;
            station.h = station.delta_star / station.theta;
            station.u = 0.97;
        }

        let state = XfoilLikeViscousState::from_station_views(&upper, &lower, 4);
        let forces = compute_forces_from_canonical_state(
            &state,
            &[0.0, 1.0, 1.0, 0.0],
            &[0.0, 0.0, 0.1, 0.1],
            0.0,
            &ViscousSolverConfig::default(),
        );

        let lower_wake: Vec<_> = state.lower_rows.iter().filter(|row| row.is_wake).collect();
        let expected_cd = compute_cd_from_wake_rows(&lower_wake);
        let te_cd = compute_cd_from_te_rows(&state.upper_rows, &state.lower_rows);

        assert!(expected_cd > 2.0 * te_cd);
        assert!((forces.cd - expected_cd).abs() < 1e-10);
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
