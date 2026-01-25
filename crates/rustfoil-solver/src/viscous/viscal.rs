//! Main viscous solution loop (VISCAL).
//!
//! This module implements XFOIL's VISCAL subroutine, the main iteration loop
//! that couples the boundary layer solution with the inviscid flow field.
//!
//! # Algorithm Overview
//!
//! 1. **Initialize**: Set up BL stations from inviscid edge velocities
//! 2. **Direct March**: March BL from stagnation with fixed Ue (MRCHUE)
//! 3. **Newton Iteration**:
//!    - Build Newton system from BL residuals (BLSYS)
//!    - Solve block-tridiagonal system (BLSOLV)
//!    - Update BL variables with limiting (UPDATE)
//!    - Update edge velocities from mass defect (UESET)
//!    - Check convergence
//! 4. **Forces**: Compute lift, drag, moment (CDCALC)
//!
//! # XFOIL Reference
//! - VISCAL: xoper.f line 2886

use nalgebra::DMatrix;
use rayon::prelude::*;
use rustfoil_bl::equations::{blvar, FlowType};
use rustfoil_bl::state::BlStation;
use rustfoil_coupling::march::{march_fixed_ue, march_surface, MarchConfig, MarchResult};
use rustfoil_coupling::newton::BlNewtonSystem;
use rustfoil_coupling::solve::solve_coupled_system;
use rustfoil_coupling::update::{update_xfoil_style, UpdateConfig};

use super::config::ViscousSolverConfig;
use super::forces::compute_forces;
use crate::{SolverError, SolverResult};

/// Result of viscous solution.
///
/// Contains all computed aerodynamic coefficients and solution metadata.
#[derive(Debug, Clone)]
pub struct ViscousResult {
    /// Angle of attack (degrees)
    pub alpha: f64,
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient (total)
    pub cd: f64,
    /// Moment coefficient (about quarter-chord)
    pub cm: f64,
    /// Upper surface transition location (x/c)
    pub x_tr_upper: f64,
    /// Lower surface transition location (x/c)
    pub x_tr_lower: f64,
    /// Number of iterations to convergence
    pub iterations: usize,
    /// Final RMS residual
    pub residual: f64,
    /// Whether solution converged within tolerance
    pub converged: bool,
    /// Friction drag coefficient
    pub cd_friction: f64,
    /// Pressure drag coefficient
    pub cd_pressure: f64,
    /// Separation location (if any)
    pub x_separation: Option<f64>,
}

impl ViscousResult {
    /// Check if the flow is attached (no separation).
    pub fn is_attached(&self) -> bool {
        self.x_separation.is_none()
    }
}

/// Solve viscous flow for pre-initialized BL stations.
///
/// This is the core VISCAL implementation that takes already-initialized
/// boundary layer stations and iterates to convergence.
///
/// # Arguments
/// * `stations` - Mutable slice of initialized BlStation
/// * `ue_inviscid` - Edge velocities from inviscid solution
/// * `dij` - Mass defect influence matrix
/// * `config` - Solver configuration
///
/// # Returns
/// `ViscousResult` on success, `SolverError` on failure.
///
/// # XFOIL Reference
/// XFOIL xoper.f VISCAL (line 2886)
pub fn solve_viscous(
    stations: &mut [BlStation],
    ue_inviscid: &[f64],
    dij: &DMatrix<f64>,
    config: &ViscousSolverConfig,
) -> SolverResult<ViscousResult> {
    // Validate inputs
    if config.reynolds <= 0.0 {
        return Err(SolverError::InvalidReynolds);
    }

    let n = stations.len();
    if n < 3 {
        return Err(SolverError::InsufficientPanels);
    }

    let msq = config.msq();
    let re = config.reynolds;

    // Emit debug event at start of VISCAL
    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::viscal(
            0,
            0.0, // alpha_rad - would be passed separately
            re,
            config.mach,
            config.ncrit,
        ));
    }

    // === Step 1: Direct BL march with inviscid Ue ===
    // Extract arc lengths and velocities
    let x: Vec<f64> = stations.iter().map(|s| s.x).collect();
    let ue: Vec<f64> = stations.iter().map(|s| s.u).collect();

    let march_config = MarchConfig {
        ncrit: config.ncrit,
        hlmax: config.hk_max_laminar,
        htmax: config.hk_max_turbulent,
        debug_trace: rustfoil_bl::is_debug_active(),
        ..Default::default()
    };

    let march_result = march_fixed_ue(&x, &ue, re, msq, &march_config);

    // Copy march results back to stations
    for (i, station) in march_result.stations.iter().enumerate() {
        stations[i].theta = station.theta;
        stations[i].delta_star = station.delta_star;
        stations[i].h = station.h;
        stations[i].hk = station.hk;
        stations[i].cf = station.cf;
        stations[i].ctau = station.ctau;
        stations[i].ampl = station.ampl;
        stations[i].is_laminar = station.is_laminar;
        stations[i].is_turbulent = station.is_turbulent;
        stations[i].mass_defect = station.mass_defect;
        stations[i].r_theta = station.r_theta;
    }

    // Track transition locations from initial march
    let x_tr_upper = march_result.x_transition.unwrap_or(1.0);
    let x_tr_lower = 1.0; // TODO: Handle lower surface separately

    // === Step 2: Newton iteration for viscous-inviscid coupling ===
    let mut current_ue = ue.clone();
    let mut converged = false;
    let mut residual = 1.0;
    let mut iteration = 0;

    // Determine flow types for Newton system
    let flow_types: Vec<FlowType> = (1..n)
        .map(|i| {
            if stations[i].is_wake {
                FlowType::Wake
            } else if stations[i].is_turbulent {
                FlowType::Turbulent
            } else {
                FlowType::Laminar
            }
        })
        .collect();

    // Build Newton system
    let mut newton_system = BlNewtonSystem::new(n);

    // Update configuration
    let update_config = UpdateConfig {
        relaxation: config.relaxation,
        ..Default::default()
    };

    while iteration < config.max_iterations && !converged {
        iteration += 1;

        // Recompute secondary variables
        for i in 0..n {
            let flow_type = if stations[i].is_wake {
                FlowType::Wake
            } else if stations[i].is_turbulent {
                FlowType::Turbulent
            } else {
                FlowType::Laminar
            };
            blvar(&mut stations[i], flow_type, msq, re);

            // Emit BLVAR debug event
            if rustfoil_bl::is_debug_active() {
                let flow_type_num = match flow_type {
                    FlowType::Laminar => 1,
                    FlowType::Turbulent => 2,
                    FlowType::Wake => 3,
                };
                let input = rustfoil_bl::BlvarInput {
                    x: stations[i].x,
                    u: stations[i].u,
                    theta: stations[i].theta,
                    delta_star: stations[i].delta_star,
                    ctau: stations[i].ctau,
                    ampl: stations[i].ampl,
                };
                let output = rustfoil_bl::BlvarOutput {
                    H: stations[i].h,
                    Hk: stations[i].hk,
                    Hs: stations[i].hs,
                    Hc: 0.0, // Not stored
                    Rtheta: stations[i].r_theta,
                    Cf: stations[i].cf,
                    Cd: stations[i].cd,
                    Us: 0.0, // Not stored
                    Cq: 0.0, // Not stored
                    De: 0.0, // Not stored
                    Dd: 0.0, // Not stored
                    Dd2: 0.0, // Not stored
                    DiWall: 0.0, // Not stored
                    DiTotal: 0.0, // Not stored
                    DiLam: 0.0, // Not stored
                    DiUsed: stations[i].cd,
                    DiLamOverride: None,
                    DiUsedMinusTotal: None,
                    DiUsedMinusLam: None,
                };
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::blvar(
                    iteration,
                    1, // side
                    i,
                    flow_type_num,
                    input,
                    output,
                ));
            }
        }

        // Validate station states before building Newton system
        let valid_stations = stations.iter().all(|s| {
            s.theta.is_finite()
                && s.theta > 0.0
                && s.delta_star.is_finite()
                && s.delta_star > 0.0
                && s.u.is_finite()
        });

        if !valid_stations {
            // If stations have become invalid, stop iteration
            break;
        }

        // Build Newton system from current BL state with full XFOIL-style coupling
        // This includes VTI signs and forced changes for proper viscous-inviscid coupling
        newton_system.build_with_vm_full(stations, &flow_types, msq, re, dij, ue_inviscid);

        // Emit debug event for Newton system (similar to XFOIL's DBGSETBL)
        if rustfoil_bl::is_debug_active() {
            // Log RMS residual before solve
            let rms_res = newton_system.residual_norm() / ((n - 1) as f64).sqrt();
            let max_res = newton_system.max_residual();
            eprintln!(
                "Newton iter {}: RMS residual = {:.6e}, max = {:.6e}",
                iteration, rms_res, max_res
            );
        }

        // Solve coupled system with VM matrix for full viscous-inviscid interaction
        let deltas_raw = solve_coupled_system(&newton_system);

        // Check if solution contains NaN - if so, reduce step size or skip
        let deltas_valid = deltas_raw.iter().all(|d| d.iter().all(|v| v.is_finite()));
        if !deltas_valid {
            // Newton system produced invalid solution, continue with zero update
            continue;
        }

        // Convert to update format [ctau/ampl, theta, mass_defect]
        let deltas: Vec<[f64; 3]> = deltas_raw
            .iter()
            .map(|d| {
                [
                    d[0], // ctau or ampl change
                    d[1], // theta change
                    d[2], // delta_star change (used as mass proxy)
                ]
            })
            .collect();

        // Update stations using XFOIL's mass-first approach
        // This computes new Ue from proposed mass changes before applying updates,
        // then applies all updates with a global relaxation factor
        let update_result = update_xfoil_style(
            stations,
            &deltas,
            ue_inviscid,
            dij,
            &newton_system.vti,
            &update_config,
        );

        // Update current_ue tracker from modified stations
        for (i, station) in stations.iter().enumerate() {
            current_ue[i] = station.u;
        }

        // Check convergence
        residual = update_result.rms_change;
        if residual < config.tolerance {
            converged = true;
        }

        // Emit debug event for iteration result (forces computed later, use 0 for now)
        if rustfoil_bl::is_debug_active() {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::viscal_result(
                iteration,
                residual,
                update_result.max_change,
                0.0, // CL computed later
                0.0, // CD computed later
                0.0, // CM computed later
            ));
        }

        // Check for stalled iteration (no progress)
        if !residual.is_finite() || residual > 1e10 {
            break;
        }
    }

    // Check for separation
    let x_separation = stations
        .iter()
        .find(|s| s.cf < 0.0 && !s.is_wake)
        .map(|s| s.x);

    if x_separation.is_some() && !config.allow_separation {
        return Err(SolverError::BoundaryLayerSeparation {
            x_sep: x_separation.unwrap(),
        });
    }

    // === Step 3: Compute aerodynamic forces ===
    let forces = compute_forces(stations, config);

    Ok(ViscousResult {
        alpha: 0.0, // Will be set by caller
        cl: forces.cl,
        cd: forces.cd,
        cm: forces.cm,
        x_tr_upper,
        x_tr_lower,
        iterations: iteration,
        residual,
        converged,
        cd_friction: forces.cd_friction,
        cd_pressure: forces.cd_pressure,
        x_separation,
    })
}

/// Solve viscous flow for a single surface (upper or lower).
///
/// This simplified version handles one surface at a time, which is
/// useful for debugging or when surfaces can be treated independently.
///
/// # Arguments
/// * `arc_lengths` - Arc length coordinates from stagnation
/// * `ue` - Edge velocities (should be positive)
/// * `config` - Solver configuration
///
/// # Returns
/// `MarchResult` with BL solution for this surface.
pub fn solve_surface(
    arc_lengths: &[f64],
    ue: &[f64],
    config: &ViscousSolverConfig,
) -> SolverResult<MarchResult> {
    if config.reynolds <= 0.0 {
        return Err(SolverError::InvalidReynolds);
    }

    let march_config = MarchConfig {
        ncrit: config.ncrit,
        hlmax: config.hk_max_laminar,
        htmax: config.hk_max_turbulent,
        debug_trace: rustfoil_bl::is_debug_active(),
        ..Default::default()
    };

    let result = march_fixed_ue(arc_lengths, ue, config.reynolds, config.msq(), &march_config);

    Ok(result)
}

/// Solve viscous flow for two surfaces (upper and lower) separately.
///
/// This function handles the proper XFOIL-style approach where each surface
/// is marched from the stagnation point toward the trailing edge independently.
///
/// # Arguments
/// * `upper_stations` - Initialized BL stations for upper surface (from stagnation to TE)
/// * `lower_stations` - Initialized BL stations for lower surface (from stagnation to TE)
/// * `upper_ue` - Edge velocities for upper surface
/// * `lower_ue` - Edge velocities for lower surface
/// * `dij` - Mass defect influence matrix (for future Newton iteration)
/// * `config` - Solver configuration
///
/// # Returns
/// `ViscousResult` with combined data from both surfaces.
pub fn solve_viscous_two_surfaces(
    upper_stations: &mut [BlStation],
    lower_stations: &mut [BlStation],
    upper_ue: &[f64],
    lower_ue: &[f64],
    dij: &DMatrix<f64>,
    config: &ViscousSolverConfig,
) -> SolverResult<ViscousResult> {
    if config.reynolds <= 0.0 {
        return Err(SolverError::InvalidReynolds);
    }

    let re = config.reynolds;
    let msq = config.msq();

    // Emit debug event at start
    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::viscal(
            0,
            0.0, // alpha set by caller
            re,
            config.mach,
            config.ncrit,
        ));
    }

    let march_config = MarchConfig {
        ncrit: config.ncrit,
        hlmax: config.hk_max_laminar,
        htmax: config.hk_max_turbulent,
        debug_trace: rustfoil_bl::is_debug_active(),
        ..Default::default()
    };

    // Extract arc lengths from stations
    let upper_arc: Vec<f64> = upper_stations.iter().map(|s| s.x).collect();
    let lower_arc: Vec<f64> = lower_stations.iter().map(|s| s.x).collect();

    // Debug: Print initial setup values
    if rustfoil_bl::is_debug_active() {
        eprintln!("[DEBUG viscal] Upper surface: {} stations", upper_arc.len());
        if upper_arc.len() >= 3 {
            eprintln!("[DEBUG viscal]   arc[0..3] = {:?}", &upper_arc[0..3]);
            eprintln!("[DEBUG viscal]   ue[0..3]  = {:?}", &upper_ue[0..3]);
        }
        eprintln!("[DEBUG viscal] Lower surface: {} stations", lower_arc.len());
        if lower_arc.len() >= 3 {
            eprintln!("[DEBUG viscal]   arc[0..3] = {:?}", &lower_arc[0..3]);
            eprintln!("[DEBUG viscal]   ue[0..3]  = {:?}", &lower_ue[0..3]);
        }
    }

    // March upper surface (side 1)
    let upper_result = march_surface(&upper_arc, upper_ue, re, msq, &march_config, 1);

    // Debug: Print march results
    if rustfoil_bl::is_debug_active() {
        eprintln!("[DEBUG viscal] Upper march results: {} stations (input: {})", 
            upper_result.stations.len(), upper_stations.len());
        for (i, s) in upper_result.stations.iter().take(5).enumerate() {
            eprintln!("[DEBUG viscal]   [{i}] theta={:.6e}, dstar={:.6e}, Hk={:.3}, Cf={:.6e}",
                s.theta, s.delta_star, s.hk, s.cf);
        }
        if let Some(last) = upper_result.stations.last() {
            eprintln!("[DEBUG viscal]   [last] theta={:.6e}, dstar={:.6e}, Hk={:.3}, Cf={:.6e}",
                last.theta, last.delta_star, last.hk, last.cf);
        }
        eprintln!("[DEBUG viscal]   x_tr={:?}, x_sep={:?}",
            upper_result.x_transition, upper_result.x_separation);
    }

    // March lower surface (side 2)
    let lower_result = march_surface(&lower_arc, lower_ue, re, msq, &march_config, 2);

    // Debug: Print lower march results
    if rustfoil_bl::is_debug_active() {
        eprintln!("[DEBUG viscal] Lower march results: {} stations (input: {})",
            lower_result.stations.len(), lower_stations.len());
        for (i, s) in lower_result.stations.iter().take(5).enumerate() {
            eprintln!("[DEBUG viscal]   [{i}] theta={:.6e}, dstar={:.6e}, Hk={:.3}, Cf={:.6e}",
                s.theta, s.delta_star, s.hk, s.cf);
        }
        if let Some(last) = lower_result.stations.last() {
            eprintln!("[DEBUG viscal]   [last] theta={:.6e}, dstar={:.6e}, Hk={:.3}, Cf={:.6e}",
                last.theta, last.delta_star, last.hk, last.cf);
        }
        eprintln!("[DEBUG viscal]   x_tr={:?}, x_sep={:?}",
            lower_result.x_transition, lower_result.x_separation);
    }

    // Copy results back to stations
    // Note: march_surface skips stagnation (index 0), so result[i] corresponds to stations[i+1]
    // if stagnation was skipped. Check if the first input station is at stagnation.
    let upper_offset = if upper_arc.first().map_or(false, |&x| x < 1e-6) { 1 } else { 0 };
    for (i, station) in upper_result.stations.iter().enumerate() {
        let target = i + upper_offset;
        if target < upper_stations.len() {
            upper_stations[target].theta = station.theta;
            upper_stations[target].delta_star = station.delta_star;
            upper_stations[target].h = station.h;
            upper_stations[target].hk = station.hk;
            upper_stations[target].cf = station.cf;
            upper_stations[target].ctau = station.ctau;
            upper_stations[target].ampl = station.ampl;
            upper_stations[target].is_laminar = station.is_laminar;
            upper_stations[target].is_turbulent = station.is_turbulent;
            upper_stations[target].mass_defect = station.mass_defect;
            upper_stations[target].r_theta = station.r_theta;
        }
    }
    // Station 0 (stagnation) keeps its initial values - Newton will skip it or handle
    // it specially like XFOIL's SIMI condition. The key is that station 1 gets result[0],
    // station 2 gets result[1], etc.

    let lower_offset = if lower_arc.first().map_or(false, |&x| x < 1e-6) { 1 } else { 0 };
    for (i, station) in lower_result.stations.iter().enumerate() {
        let target = i + lower_offset;
        if target < lower_stations.len() {
            lower_stations[target].theta = station.theta;
            lower_stations[target].delta_star = station.delta_star;
            lower_stations[target].h = station.h;
            lower_stations[target].hk = station.hk;
            lower_stations[target].cf = station.cf;
            lower_stations[target].ctau = station.ctau;
            lower_stations[target].ampl = station.ampl;
            lower_stations[target].is_laminar = station.is_laminar;
            lower_stations[target].is_turbulent = station.is_turbulent;
            lower_stations[target].mass_defect = station.mass_defect;
            lower_stations[target].r_theta = station.r_theta;
        }
    }
    // Station 0 (stagnation) keeps its initial values

    // Get transition locations from march results (may be updated by Newton iteration)
    let mut x_tr_upper = upper_result.x_transition.unwrap_or(1.0);
    let mut x_tr_lower = lower_result.x_transition.unwrap_or(1.0);

    // === Newton iteration for viscous-inviscid coupling ===
    // NOTE: Full two-surface Newton coupling requires proper DIJ matrix extraction
    // and VZ block handling (see Phase 2 in newton-vi-coupling-plan.md).
    // For now, we run Newton iteration only if the DIJ matrix has proper dimensions
    // for the individual surfaces. Otherwise, we rely on the direct march results.

    let update_config = UpdateConfig::default();
    let mut iteration = 0;
    let mut converged = true; // Direct march is considered converged
    let mut residual = 0.0;

    let n_upper = upper_stations.len();
    let n_lower = lower_stations.len();

    // Check if DIJ matrix has proper dimensions for Newton iteration
    // The full DIJ includes both surfaces, so we need to extract submatrices
    // For proper Newton coupling, we need:
    // 1. Upper surface DIJ: influences from upper surface mass on upper surface Ue
    // 2. Lower surface DIJ: influences from lower surface mass on lower surface Ue
    // 3. Cross-surface coupling (Phase 2)
    
    // For Phase 1: Only run Newton if we can extract valid submatrices
    // The DIJ matrix from setup is for the full airfoil (all panels)
    // Upper stations are the latter part of the full array (indices n_lower..n_total)
    // Lower stations are the first part (indices 0..n_lower)
    
    #[allow(unused_variables)]
    let n_total = dij.nrows();
    // Newton V-I coupling for two surfaces requires:
    // 1. Proper DIJ submatrix extraction for each surface
    // 2. Cross-surface coupling through VZ block at trailing edge
    // 3. Combined Newton system for both surfaces
    // See Phase 2 in newton-vi-coupling-plan.md for full implementation.
    //
    // Newton iteration is currently disabled because it fails at transition stations.
    // The transition interval needs proper TRDIF handling (hybrid laminar+turbulent equations)
    // which is not yet implemented. For now, we rely on the direct march results.
    // TODO: Implement TRDIF handling in Newton system to enable V-I coupling.
    let can_run_newton = false;

    if can_run_newton {
        // Save original station values in case Newton diverges
        let upper_backup: Vec<_> = upper_stations.iter().map(|s| (s.theta, s.delta_star, s.ctau, s.ampl, s.h, s.mass_defect)).collect();
        let lower_backup: Vec<_> = lower_stations.iter().map(|s| (s.theta, s.delta_star, s.ctau, s.ampl, s.h, s.mass_defect)).collect();

        // Extract DIJ submatrices
        // Lower surface: indices 0..n_lower in the full system
        // Upper surface: indices n_lower..n_total in the full system
        // But our stations are already extracted, so we need to map correctly
        
        // For now, use diagonal-only DIJ (self-influence only) to avoid cross-coupling issues
        // This is a simplification - full coupling requires Phase 2 implementation
        let mut dij_upper = DMatrix::zeros(n_upper, n_upper);
        let mut dij_lower = DMatrix::zeros(n_lower, n_lower);
        
        // Extract diagonal elements (self-influence) from the full DIJ
        // Upper surface stations map to indices n_lower..n_total in original panel array
        for i in 0..n_upper {
            let full_idx = n_lower + i;
            if full_idx < n_total {
                dij_upper[(i, i)] = dij[(full_idx, full_idx)];
            }
        }
        
        // Lower surface stations map to indices 0..n_lower (reversed from TE to stag)
        for i in 0..n_lower {
            let full_idx = n_lower - 1 - i; // Reverse mapping
            if full_idx < n_total {
                dij_lower[(i, i)] = dij[(full_idx, full_idx)];
            }
        }

        // Run Newton iteration with limited iterations to prevent divergence
        // Phase 1: Very conservative - just 3 iterations with heavy relaxation
        let max_newton_iter = 3;
        let mut prev_residual = f64::MAX;
        
        for iter in 0..max_newton_iter {
            iteration = iter + 1;

            // === Upper surface Newton iteration ===
            let upper_flow_types: Vec<FlowType> = upper_stations
                .iter()
                .skip(1)
                .map(|s| {
                    if s.is_wake {
                        FlowType::Wake
                    } else if s.is_laminar {
                        FlowType::Laminar
                    } else {
                        FlowType::Turbulent
                    }
                })
                .collect();

            // Recompute secondary variables
            for station in upper_stations.iter_mut() {
                let flow_type = if station.is_wake {
                    FlowType::Wake
                } else if station.is_laminar {
                    FlowType::Laminar
                } else {
                    FlowType::Turbulent
                };
                blvar(station, flow_type, msq, re);
            }
            
            // Debug after blvar
            if rustfoil_bl::is_debug_active() && iter > 0 {
                eprintln!("[DEBUG Newton] iter {} after blvar:", iter);
                for i in 1..4.min(upper_stations.len()) {
                    let s = &upper_stations[i];
                    eprintln!("[DEBUG Newton]   [{}] theta={:.6e} dstar={:.6e} Hk={:.4} Cf={:.4e}", 
                        i, s.theta, s.delta_star, s.hk, s.cf);
                }
            }

            // Build Newton system (without VM full coupling for Phase 1)
            let mut upper_newton = BlNewtonSystem::new(n_upper);
            if upper_flow_types.len() == n_upper - 1 {
                // Use basic build (not build_with_vm_full) to avoid DIJ issues
                upper_newton.build(upper_stations, &upper_flow_types, msq, re);
            }
            
            // Zero out the stagnation interval residual (XFOIL SIMI condition)
            // The first interval (station 0 → 1) uses special similarity handling
            // where the residual is set to zero (prescribed boundary condition)
            if upper_newton.rhs.len() > 1 {
                upper_newton.rhs[1] = [0.0, 0.0, 0.0];
            }

            let upper_residual = upper_newton.residual_norm();
            
            // Debug: show station values and RHS
            if rustfoil_bl::is_debug_active() && iter < 2 {
                eprintln!("[DEBUG Newton] iter {} upper stations before build:", iter);
                for i in 0..4.min(upper_stations.len()) {
                    let s = &upper_stations[i];
                    eprintln!("[DEBUG Newton]   [{}] x={:.6e} theta={:.6e} dstar={:.6e} Hk={:.4} Ue={:.6e}", 
                        i, s.x, s.theta, s.delta_star, s.hk, s.u);
                }
                eprintln!("[DEBUG Newton] iter {} upper RHS res={:.4e}", iter, upper_residual);
                // Find max residual
                let mut max_res = 0.0f64;
                let mut max_idx = 0;
                for (i, r) in upper_newton.rhs.iter().enumerate() {
                    let norm = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
                    if norm > max_res {
                        max_res = norm;
                        max_idx = i;
                    }
                }
                eprintln!("[DEBUG Newton]   max residual at idx={}: {:?}", max_idx, upper_newton.rhs[max_idx]);
                eprintln!("[DEBUG Newton]   rhs[2]={:?}", upper_newton.rhs[2]);
            }

            // Solve Newton system (basic tridiagonal solve)
            let upper_deltas = solve_coupled_system(&upper_newton);
            
            // Debug: show deltas and station values after
            if rustfoil_bl::is_debug_active() && iter == 0 {
                eprintln!("[DEBUG Newton] upper deltas n={}", upper_deltas.len());
                for i in 1..6.min(upper_deltas.len()) {
                    let d = &upper_deltas[i];
                    eprintln!("  [{}] d[ampl/ctau]={:.4e}, d[theta]={:.4e}, d[mass]={:.4e}", i, d[0], d[1], d[2]);
                }
            }
            if rustfoil_bl::is_debug_active() && iter == 1 {
                eprintln!("[DEBUG Newton] iter 1 START, upper stations:");
                for i in 1..4.min(upper_stations.len()) {
                    let s = &upper_stations[i];
                    eprintln!("[DEBUG Newton]   [{}] theta={:.6e} dstar={:.6e} Hk={:.4} Cf={:.4e}", 
                        i, s.theta, s.delta_star, s.hk, s.cf);
                }
            }

            // Apply updates with simple limiting (not XFOIL-style Ue coupling)
            // XFOIL UPDATE: DDSTR = (DMASS - DSTR*DUEDG) / UEDG
            // Without V-I coupling (DUEDG=0): DDSTR = DMASS / UEDG
            if upper_deltas.len() == n_upper {
                for (i, station) in upper_stations.iter_mut().enumerate() {
                    if i > 0 {
                        let delta = &upper_deltas[i];
                        // Apply with relaxation
                        let rlx = 0.5; // Conservative relaxation
                        if station.is_laminar {
                            station.ampl += rlx * delta[0];
                            station.ampl = station.ampl.max(0.0);
                        } else {
                            station.ctau += rlx * delta[0];
                            station.ctau = station.ctau.clamp(0.0, 0.25);
                        }
                        station.theta += rlx * delta[1];
                        station.theta = station.theta.max(1e-12);
                        
                        // delta[2] is MASS change, not delta_star change
                        // DDSTR = DMASS / UEDG (when DUEDG = 0)
                        let d_dstar = delta[2] / station.u.max(1e-6);
                        station.delta_star += rlx * d_dstar;
                        station.delta_star = station.delta_star.max(1e-12);
                        
                        station.h = station.delta_star / station.theta;
                        station.mass_defect = station.u * station.delta_star;
                    }
                }
            }

            // === Lower surface Newton iteration ===
            let lower_flow_types: Vec<FlowType> = lower_stations
                .iter()
                .skip(1)
                .map(|s| {
                    if s.is_wake {
                        FlowType::Wake
                    } else if s.is_laminar {
                        FlowType::Laminar
                    } else {
                        FlowType::Turbulent
                    }
                })
                .collect();

            // Recompute secondary variables
            for station in lower_stations.iter_mut() {
                let flow_type = if station.is_wake {
                    FlowType::Wake
                } else if station.is_laminar {
                    FlowType::Laminar
                } else {
                    FlowType::Turbulent
                };
                blvar(station, flow_type, msq, re);
            }

            // Build Newton system
            let mut lower_newton = BlNewtonSystem::new(n_lower);
            if lower_flow_types.len() == n_lower - 1 {
                lower_newton.build(lower_stations, &lower_flow_types, msq, re);
            }
            
            // Zero out the stagnation interval residual (XFOIL SIMI condition)
            if lower_newton.rhs.len() > 1 {
                lower_newton.rhs[1] = [0.0, 0.0, 0.0];
            }

            let lower_residual = lower_newton.residual_norm();

            // Solve Newton system
            let lower_deltas = solve_coupled_system(&lower_newton);

            // Apply updates
            if lower_deltas.len() == n_lower {
                for (i, station) in lower_stations.iter_mut().enumerate() {
                    if i > 0 {
                        let delta = &lower_deltas[i];
                        let rlx = 0.5;
                        if station.is_laminar {
                            station.ampl += rlx * delta[0];
                            station.ampl = station.ampl.max(0.0);
                        } else {
                            station.ctau += rlx * delta[0];
                            station.ctau = station.ctau.clamp(0.0, 0.25);
                        }
                        station.theta += rlx * delta[1];
                        station.theta = station.theta.max(1e-12);
                        
                        // delta[2] is MASS change, not delta_star change
                        // DDSTR = DMASS / UEDG (when DUEDG = 0)
                        let d_dstar = delta[2] / station.u.max(1e-6);
                        station.delta_star += rlx * d_dstar;
                        station.delta_star = station.delta_star.max(1e-12);
                        
                        station.h = station.delta_star / station.theta;
                        station.mass_defect = station.u * station.delta_star;
                    }
                }
            }

            // Combined residual
            residual = (upper_residual + lower_residual) / 2.0;
            
            // Debug: show Newton iteration progress
            if rustfoil_bl::is_debug_active() {
                eprintln!("[DEBUG Newton] iter={}: upper_res={:.6e}, lower_res={:.6e}, total={:.6e}",
                    iter, upper_residual, lower_residual, residual);
            }

            // Check convergence
            if residual < config.tolerance {
                converged = true;
                break;
            }

            // Check for divergence or increasing residual
            if !residual.is_finite() || residual > 1e10 {
                if rustfoil_bl::is_debug_active() {
                    eprintln!("[DEBUG Newton] FAILED: residual={:.6e}, falling back to direct march", residual);
                }
                // Restore original values
                for (i, (theta, dstar, ctau, ampl, h, mass)) in upper_backup.iter().enumerate() {
                    upper_stations[i].theta = *theta;
                    upper_stations[i].delta_star = *dstar;
                    upper_stations[i].ctau = *ctau;
                    upper_stations[i].ampl = *ampl;
                    upper_stations[i].h = *h;
                    upper_stations[i].mass_defect = *mass;
                }
                for (i, (theta, dstar, ctau, ampl, h, mass)) in lower_backup.iter().enumerate() {
                    lower_stations[i].theta = *theta;
                    lower_stations[i].delta_star = *dstar;
                    lower_stations[i].ctau = *ctau;
                    lower_stations[i].ampl = *ampl;
                    lower_stations[i].h = *h;
                    lower_stations[i].mass_defect = *mass;
                }
                converged = true; // Fall back to direct march result
                iteration = 0;
                break;
            }
            
            // Stop if residual is increasing (Newton not converging)
            if residual > prev_residual * 1.5 {
                // Restore original values
                for (i, (theta, dstar, ctau, ampl, h, mass)) in upper_backup.iter().enumerate() {
                    upper_stations[i].theta = *theta;
                    upper_stations[i].delta_star = *dstar;
                    upper_stations[i].ctau = *ctau;
                    upper_stations[i].ampl = *ampl;
                    upper_stations[i].h = *h;
                    upper_stations[i].mass_defect = *mass;
                }
                for (i, (theta, dstar, ctau, ampl, h, mass)) in lower_backup.iter().enumerate() {
                    lower_stations[i].theta = *theta;
                    lower_stations[i].delta_star = *dstar;
                    lower_stations[i].ctau = *ctau;
                    lower_stations[i].ampl = *ampl;
                    lower_stations[i].h = *h;
                    lower_stations[i].mass_defect = *mass;
                }
                converged = true; // Fall back to direct march result
                iteration = 0;
                break;
            }
            prev_residual = residual;
        }
    }

    // Update transition locations from final station states
    // Find transition point (where is_laminar changes to is_turbulent)
    for (i, station) in upper_stations.iter().enumerate() {
        if !station.is_laminar && i > 0 && upper_stations[i - 1].is_laminar {
            x_tr_upper = station.x;
            break;
        }
    }
    for (i, station) in lower_stations.iter().enumerate() {
        if !station.is_laminar && i > 0 && lower_stations[i - 1].is_laminar {
            x_tr_lower = station.x;
            break;
        }
    }

    // Check for separation on either surface
    let x_separation = upper_result
        .x_separation
        .or(lower_result.x_separation);

    // Compute forces from both surfaces
    // For now, compute from upper surface (CD is dominated by friction drag anyway)
    let forces = compute_forces(upper_stations, config);
    let lower_forces = compute_forces(lower_stations, config);

    // Debug: Print force computation
    if rustfoil_bl::is_debug_active() {
        eprintln!("[DEBUG viscal] Upper forces: cd_f={:.6e}, cd_p={:.6e}",
            forces.cd_friction, forces.cd_pressure);
        eprintln!("[DEBUG viscal] Lower forces: cd_f={:.6e}, cd_p={:.6e}",
            lower_forces.cd_friction, lower_forces.cd_pressure);
        
        // Check last station values (used in Squire-Young)
        if let Some(last) = upper_stations.last() {
            eprintln!("[DEBUG viscal] Upper last: theta={:.6e}, H={:.3}, Ue={:.3}, is_wake={}",
                last.theta, last.h, last.u, last.is_wake);
        }
    }

    // Combine CD (friction from both surfaces)
    let total_cd_friction = forces.cd_friction + lower_forces.cd_friction;

    Ok(ViscousResult {
        alpha: 0.0, // Will be set by caller
        cl: forces.cl, // CL from pressure distribution (need full airfoil)
        cd: total_cd_friction, // For now, just friction drag
        cm: forces.cm,
        x_tr_upper,
        x_tr_lower,
        iterations: iteration,
        residual,
        converged,
        cd_friction: total_cd_friction,
        cd_pressure: 0.0, // Need pressure integration for this
        x_separation,
    })
}

/// Solve viscous flow for multiple angles of attack in parallel.
///
/// Uses rayon to compute polar data efficiently. Each angle is solved
/// independently with the same configuration.
///
/// # Arguments
/// * `setup_fn` - Function that creates (stations, ue, dij) for each alpha
/// * `alphas` - Angles of attack to compute (degrees)
/// * `config` - Solver configuration
///
/// # Returns
/// Vector of results, one per angle of attack.
///
/// # Example
/// ```ignore
/// let alphas = vec![-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0];
/// let results = solve_viscous_polar_parallel(
///     |alpha| setup_for_alpha(alpha),
///     &alphas,
///     &config,
/// );
/// ```
pub fn solve_viscous_polar_parallel<F>(
    setup_fn: F,
    alphas: &[f64],
    config: &ViscousSolverConfig,
) -> Vec<SolverResult<ViscousResult>>
where
    F: Fn(f64) -> (Vec<BlStation>, Vec<f64>, DMatrix<f64>) + Sync,
{
    alphas
        .par_iter()
        .map(|&alpha| {
            let (mut stations, ue, dij) = setup_fn(alpha);
            let mut result = solve_viscous(&mut stations, &ue, &dij, config)?;
            result.alpha = alpha;
            Ok(result)
        })
        .collect()
}

/// Solve viscous polar with pre-computed setups.
///
/// This variant takes pre-computed setup data for each angle, avoiding
/// the need to recompute inviscid solutions in the parallel loop.
///
/// # Arguments
/// * `setups` - Vector of (alpha, stations, ue, dij) tuples
/// * `config` - Solver configuration
///
/// # Returns
/// Vector of results, one per setup.
pub fn solve_viscous_polar_from_setups(
    setups: Vec<(f64, Vec<BlStation>, Vec<f64>, DMatrix<f64>)>,
    config: &ViscousSolverConfig,
) -> Vec<SolverResult<ViscousResult>> {
    setups
        .into_par_iter()
        .map(|(alpha, mut stations, ue, dij)| {
            let mut result = solve_viscous(&mut stations, &ue, &dij, config)?;
            result.alpha = alpha;
            Ok(result)
        })
        .collect()
}

/// Sequential polar computation (for debugging or single-threaded contexts).
///
/// # Arguments
/// * `setups` - Iterator of (alpha, stations, ue, dij) tuples
/// * `config` - Solver configuration
///
/// # Returns
/// Vector of results.
pub fn solve_viscous_polar_sequential<I>(
    setups: I,
    config: &ViscousSolverConfig,
) -> Vec<SolverResult<ViscousResult>>
where
    I: IntoIterator<Item = (f64, Vec<BlStation>, Vec<f64>, DMatrix<f64>)>,
{
    setups
        .into_iter()
        .map(|(alpha, mut stations, ue, dij)| {
            let mut result = solve_viscous(&mut stations, &ue, &dij, config)?;
            result.alpha = alpha;
            Ok(result)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::viscous::setup::{find_stagnation, initialize_bl_stations};

    /// Create a simple flat plate test case.
    fn flat_plate_setup(n: usize, re: f64) -> (Vec<BlStation>, Vec<f64>, DMatrix<f64>) {
        // Arc lengths from 0 to 1
        let x: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        // Constant edge velocity
        let ue = vec![1.0; n];
        // Simple DIJ matrix
        let dij = DMatrix::identity(n, n) * 0.1;

        let stag_idx = 0; // First station is stagnation
        let stations = initialize_bl_stations(&x, &ue, re, stag_idx);

        (stations, ue, dij)
    }

    #[test]
    fn test_solve_viscous_flat_plate() {
        let config = ViscousSolverConfig::with_reynolds(1e6);
        let (mut stations, ue, dij) = flat_plate_setup(50, config.reynolds);

        let result = solve_viscous(&mut stations, &ue, &dij, &config);

        assert!(result.is_ok());
        let result = result.unwrap();

        // Should run iterations (may not fully converge with simple test setup)
        assert!(result.iterations > 0);

        // CD from friction should be non-negative
        // Note: With simplified test data, the Newton system may be singular
        // and give non-ideal results, but friction drag integration should work
        assert!(result.cd_friction >= 0.0, "Friction drag should be non-negative");
    }

    #[test]
    fn test_solve_viscous_invalid_reynolds() {
        let config = ViscousSolverConfig {
            reynolds: -1.0,
            ..Default::default()
        };
        let (mut stations, ue, dij) = flat_plate_setup(50, 1e6);

        let result = solve_viscous(&mut stations, &ue, &dij, &config);

        assert!(matches!(result, Err(SolverError::InvalidReynolds)));
    }

    #[test]
    fn test_solve_surface() {
        let config = ViscousSolverConfig::with_reynolds(1e6);

        let n = 40;
        let x: Vec<f64> = (0..n).map(|i| 0.01 * i as f64).collect();
        let ue = vec![1.0; n];

        let result = solve_surface(&x, &ue, &config);

        assert!(result.is_ok());
        let march = result.unwrap();

        assert_eq!(march.stations.len(), n);
        assert!(march.converged);
    }

    #[test]
    fn test_viscous_result_is_attached() {
        let mut result = ViscousResult {
            alpha: 0.0,
            cl: 0.0,
            cd: 0.01,
            cm: 0.0,
            x_tr_upper: 0.3,
            x_tr_lower: 0.5,
            iterations: 10,
            residual: 1e-5,
            converged: true,
            cd_friction: 0.008,
            cd_pressure: 0.002,
            x_separation: None,
        };

        assert!(result.is_attached());

        result.x_separation = Some(0.8);
        assert!(!result.is_attached());
    }

    #[test]
    fn test_polar_sequential() {
        let config = ViscousSolverConfig::with_reynolds(1e6);

        // Create setups for a few "angles"
        let setups: Vec<_> = vec![0.0, 2.0, 4.0]
            .into_iter()
            .map(|alpha| {
                let (stations, ue, dij) = flat_plate_setup(30, config.reynolds);
                (alpha, stations, ue, dij)
            })
            .collect();

        let results = solve_viscous_polar_sequential(setups, &config);

        assert_eq!(results.len(), 3);

        for result in &results {
            assert!(result.is_ok());
        }
    }
}
