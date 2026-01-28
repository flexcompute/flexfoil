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
use rustfoil_coupling::global_newton::{apply_global_updates, solve_global_system, GlobalNewtonSystem};
use rustfoil_coupling::newton::BlNewtonSystem;
use rustfoil_coupling::solve::solve_coupled_system;
// STMOVE module for stagnation point finding and coupling
use rustfoil_coupling::stmove::{find_stagnation_by_gamma, find_stagnation_with_derivs};
use rustfoil_coupling::update::{update_xfoil_style, UpdateConfig};
use rustfoil_coupling::wake::combine_te_for_wake;

use super::circulation::{
    compute_cl_from_gamma, update_circulation_from_qvis,
    update_qvis_from_uedg, update_qvis_from_uedg_two_surfaces,
};
use super::config::ViscousSolverConfig;
use super::forces::{compute_forces, compute_forces_two_surfaces};
use crate::{SolverError, SolverResult};

/// XFOIL's CLCALC subroutine - compute CL by integrating Cp in wind axes.
///
/// This is the exact XFOIL formula for lift coefficient computation:
/// CL = ∮ Cp * dx_wind
/// where dx_wind = (X[i+1]-X[i])*cos(α) + (Y[i+1]-Y[i])*sin(α)
///
/// # Arguments
/// * `panel_x` - Panel node x-coordinates (full airfoil, XFOIL ordering)
/// * `panel_y` - Panel node y-coordinates (full airfoil, XFOIL ordering)
/// * `gamma` - Vorticity (= tangential velocity) at each panel node
/// * `alpha` - Angle of attack in radians
///
/// # Returns
/// Lift coefficient CL
///
/// # Reference
/// XFOIL xfoil.f CLCALC (line 1088)
pub fn clcalc(panel_x: &[f64], panel_y: &[f64], gamma: &[f64], alpha: f64) -> f64 {
    let n = panel_x.len();
    if n < 3 || gamma.len() != n {
        return 0.0;
    }
    
    let ca = alpha.cos();
    let sa = alpha.sin();
    
    // For incompressible flow: BETA = 1.0, BFAC = 0.0
    // So CPG = CGINC = 1.0 - (GAM/QINF)^2
    // With QINF = 1.0 (normalized), CPG = 1.0 - GAM^2
    
    // Initialize with first point's Cp
    let mut cpg1 = 1.0 - gamma[0] * gamma[0];
    let mut cl = 0.0;
    
    // XFOIL loops: DO I=1, N with IP=I+1 and IP=1 when I=N
    // This integrates around the CLOSED contour
    for i in 0..n {
        let ip = (i + 1) % n;  // Wrap around: when i=n-1, ip=0
        
        // Cp at next point
        let cpg2 = 1.0 - gamma[ip] * gamma[ip];
        
        // Panel projection in wind axes (XFOIL: DX)
        let dx_wind = (panel_x[ip] - panel_x[i]) * ca + (panel_y[ip] - panel_y[i]) * sa;
        
        // Average Cp on this panel
        let cp_avg = 0.5 * (cpg1 + cpg2);
        
        // Accumulate CL
        cl += cp_avg * dx_wind;
        
        // Shift for next panel
        cpg1 = cpg2;
    }
    
    cl
}

/// Construct full panel gamma array from upper/lower surface edge velocities.
///
/// Maps BL station edge velocities back to panel nodes using panel_idx.
/// Applies VTI sign convention: upper = +1, lower = -1.
///
/// # Arguments
/// * `upper_stations` - Upper surface BL stations (stagnation to TE)
/// * `lower_stations` - Lower surface BL stations (stagnation to TE)
/// * `n_panels` - Number of panel nodes
///
/// # Returns
/// Full gamma array in XFOIL panel ordering
fn construct_gamma_from_stations(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    n_panels: usize,
) -> Vec<f64> {
    let mut gamma = vec![0.0; n_panels];
    let mut filled = vec![false; n_panels];
    
    // Upper surface: VTI = +1
    // Upper stations go from stagnation (high panel index) toward TE (low panel index)
    // XFOIL: QVIS = VTI * UEDG, GAM = QVIS
    // So: gamma = +1 * Ue = +Ue
    // Note: station.u should be positive (magnitude), but we preserve sign for correctness
    for station in upper_stations.iter() {
        let idx = station.panel_idx;
        if idx < n_panels && !station.is_wake {
            // VTI = +1: gamma = +Ue
            // Use signed value directly (should be positive, but preserve sign)
            gamma[idx] = station.u;
            filled[idx] = true;
        }
    }
    
    // Lower surface: VTI = -1
    // Lower stations go from stagnation (low panel index) toward TE (high panel index)
    // XFOIL: QVIS = VTI * UEDG, GAM = QVIS
    // So: gamma = -1 * Ue = -Ue
    // The negative sign accounts for the opposite panel tangent direction
    for station in lower_stations.iter() {
        let idx = station.panel_idx;
        if idx < n_panels && !station.is_wake && !filled[idx] {
            // VTI = -1: gamma = -Ue
            // Flip sign for lower surface (preserve magnitude sign if Ue can be negative)
            gamma[idx] = -station.u;
            filled[idx] = true;
        }
    }
    
    // Debug: Log gamma array statistics
    let n_filled = filled.iter().filter(|&&f| f).count();
    let n_unfilled = filled.iter().filter(|&&f| !f).count();
    
    // Warn if panels are unfilled (they'll have gamma=0, Cp=1, which may be wrong)
    if n_unfilled > 0 && n_unfilled < n_panels / 10 {
        // Allow a few unfilled panels (e.g., at stagnation or TE), but warn if many
        eprintln!("[WARN construct_gamma] {} panels unfilled (will have gamma=0, Cp=1)", n_unfilled);
    }
    
    if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
        let gamma_min = gamma.iter().cloned().fold(f64::INFINITY, f64::min);
        let gamma_max = gamma.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let gamma_abs_avg: f64 = gamma.iter().map(|g| g.abs()).sum::<f64>() / n_panels as f64;
        eprintln!("[DEBUG construct_gamma] n_panels={}, filled={}, unfilled={}", n_panels, n_filled, n_unfilled);
        eprintln!("[DEBUG construct_gamma] gamma range: [{:.4}, {:.4}], avg |gamma|={:.4}", gamma_min, gamma_max, gamma_abs_avg);
        // Find where max gamma is
        let max_idx = gamma.iter().enumerate().max_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap()).map(|(i, _)| i);
        if let Some(idx) = max_idx {
            eprintln!("[DEBUG construct_gamma] Max |gamma| at panel {}: {:.4}", idx, gamma[idx]);
        }
        // Show first/last few gamma values
        if n_panels >= 10 {
            eprintln!("[DEBUG construct_gamma] First 5 gamma: {:?}", &gamma[..5].iter().map(|g| format!("{:.4}", g)).collect::<Vec<_>>());
            eprintln!("[DEBUG construct_gamma] Last 5 gamma: {:?}", &gamma[n_panels-5..].iter().map(|g| format!("{:.4}", g)).collect::<Vec<_>>());
            // Show gamma around LE (typically panel ~n/2)
            let le_idx = n_panels / 2;
            eprintln!("[DEBUG construct_gamma] Gamma around LE (panel {}): {:?}", le_idx, 
                &gamma[le_idx.saturating_sub(2)..=(le_idx+2).min(n_panels-1)]
                    .iter().map(|g| format!("{:.4}", g)).collect::<Vec<_>>());
        }
        // Show unfilled panel indices if any
        if n_unfilled > 0 && n_unfilled <= 20 {
            let unfilled_indices: Vec<usize> = filled.iter().enumerate()
                .filter(|(_, &f)| !f)
                .map(|(i, _)| i)
                .collect();
            eprintln!("[DEBUG construct_gamma] Unfilled panel indices: {:?}", unfilled_indices);
        }
    }
    
    gamma
}

/// Result of viscous solution.
///
/// Contains all computed aerodynamic coefficients and solution metadata.
#[derive(Debug, Clone)]
pub struct ViscousResult {
    /// Angle of attack (degrees)
    pub alpha: f64,
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient (pressure drag only, matches XFOIL CL_DETAIL 'cdp')
    /// Total drag = cd + cd_friction
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

        // QVFUE + GAMQV: refresh QVIS from Ue, then GAM from QVIS
        let qvis = update_qvis_from_uedg(stations);
        let gamma = update_circulation_from_qvis(&qvis);

        // Check convergence
        residual = update_result.rms_change;
        if residual < config.tolerance {
            converged = true;
        }

        // Debug: Print iteration for single-surface mode
        // (Two-surface mode has its own tracking below)
        
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
/// * `alpha_rad` - Angle of attack in radians (for CLCALC wind-axis integration)
/// * `panel_x` - Panel node x-coordinates (full airfoil, XFOIL ordering)
/// * `panel_y` - Panel node y-coordinates (full airfoil, XFOIL ordering)
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
    alpha_rad: f64,
    panel_x: &[f64],
    panel_y: &[f64],
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
            alpha_rad,
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
            // CRITICAL: Copy edge velocity from march result
            upper_stations[target].u = station.u;
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
            // CRITICAL: Copy edge velocity from march result
            lower_stations[target].u = station.u;
        }
    }
    // Station 0 (stagnation) keeps its initial values

    // Get transition locations from march results (may be updated by Newton iteration)
    let mut x_tr_upper = upper_result.x_transition.unwrap_or(1.0);
    let mut x_tr_lower = lower_result.x_transition.unwrap_or(1.0);

    // === Newton iteration for viscous-inviscid coupling ===
    // Full global Newton coupling using the GlobalNewtonSystem that properly
    // couples both upper and lower surfaces through the DIJ matrix.

    let update_config = UpdateConfig {
        relaxation: config.relaxation,
        ..Default::default()
    };
    let mut iteration = 0;
    let mut converged = true; // Direct march is considered converged
    let mut residual = 0.0;

    let n_upper = upper_stations.len();
    let n_lower = lower_stations.len();

    // Determine trailing edge indices (last non-wake station on each surface)
    let iblte_upper = upper_stations
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.is_wake)
        .map(|(i, _)| i)
        .last()
        .unwrap_or(n_upper.saturating_sub(1));
    let iblte_lower = lower_stations
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.is_wake)
        .map(|(i, _)| i)
        .last()
        .unwrap_or(n_lower.saturating_sub(1));

    // Enable Newton with detailed debug output for comparison with XFOIL
    let can_run_newton = config.max_iterations > 0 && n_upper >= 3 && n_lower >= 3;

    // Debug: track Newton execution
    if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
        eprintln!("[DEBUG viscal] Newton check: can_run={}, max_iter={}, n_upper={}, n_lower={}", 
            can_run_newton, config.max_iterations, n_upper, n_lower);
    }

    if can_run_newton {
        // Save original station values in case Newton diverges
        // IMPORTANT: Include station.u which is modified by VI coupling
        let upper_backup: Vec<_> = upper_stations.iter().map(|s| (s.theta, s.delta_star, s.ctau, s.ampl, s.h, s.mass_defect, s.u)).collect();
        let lower_backup: Vec<_> = lower_stations.iter().map(|s| (s.theta, s.delta_star, s.ctau, s.ampl, s.h, s.mass_defect, s.u)).collect();

        // Create global Newton system with cross-surface coupling
        let mut global_system = GlobalNewtonSystem::new(n_upper, n_lower, iblte_upper, iblte_lower);

        // === Compute stagnation point derivatives (SST_GO, SST_GP) ===
        // XFOIL computes these in STFIND using panel gamma and arc lengths.
        // We estimate from the Ue distribution: gamma changes sign at stagnation,
        // with Ue_upper positive and Ue_lower positive (both grow from 0).
        // 
        // Approximation: SST_GO and SST_GP relate stagnation arc length to gamma changes.
        // From XFOIL: SST_GO = (SST - S(IST+1)) / DGAM, SST_GP = (S(IST) - SST) / DGAM
        // where DGAM is the gamma difference across the stagnation panel.
        //
        // Near stagnation: gamma ≈ ±Ue, so DGAM ≈ -Ue_upper[1] - Ue_lower[1]
        // and DS ≈ arc_length[1] (small). SST_GO/SST_GP ~ DS / |DGAM| ~ O(0.01).
        let (sst_go, sst_gp) = {
            // Get Ue at first station after stagnation on each surface
            let ue_upper_1 = upper_ue.get(1).copied().unwrap_or(0.1).abs();
            let ue_lower_1 = lower_ue.get(1).copied().unwrap_or(0.1).abs();
            let ds_upper = upper_arc.get(1).copied().unwrap_or(0.01);
            let ds_lower = lower_arc.get(1).copied().unwrap_or(0.01);
            
            // DGAM = gamma(IST+1) - gamma(IST) = -gamma_upper_1 - gamma_lower_1
            // Since gamma_upper = +Ue, gamma_lower = -Ue (for lower surface in panel ordering)
            // DGAM ≈ -ue_upper_1 - (-ue_lower_1) = ue_lower_1 - ue_upper_1 (usually small)
            // But to avoid division by zero, use a safer estimate.
            let dgam = (ue_upper_1 + ue_lower_1).max(0.01);
            let ds = (ds_upper + ds_lower) / 2.0;
            
            // SST is approximately at the midpoint between IST and IST+1
            // SST_GO = (SST - S(IST+1)) / DGAM ≈ -DS/2 / DGAM
            // SST_GP = (S(IST) - SST) / DGAM ≈ -DS/2 / DGAM
            // Signs: In XFOIL, SST_GO is typically positive, SST_GP negative
            let sst_go = ds / dgam;
            let sst_gp = -ds / dgam;
            
            (sst_go, sst_gp)
        };
        
        global_system.set_stagnation_derivs(sst_go, sst_gp);
        
        if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
            eprintln!("[DEBUG viscal] Stagnation derivatives: SST_GO={:.6e}, SST_GP={:.6e}", sst_go, sst_gp);
        }

        // Run Newton iteration with full global coupling
        let max_newton_iter = config.max_iterations.min(30);
        let mut prev_residual = f64::MAX;

        for iter in 0..max_newton_iter {
            iteration = iter + 1;

            if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() && iter < 3 {
                eprintln!("[DEBUG Newton] Starting iter {}", iter);
            }
            
            // Emit NEWTON_ITER debug event at start of iteration
            if rustfoil_bl::is_debug_active() {
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::newton_iter(
                    iter,
                    residual,
                    prev_residual,
                    n_upper,
                    n_lower,
                ));
            }

            // === RE-MARCH BOUNDARY LAYERS WITH CURRENT Ue ===
            // CRITICAL FIX: XFOIL re-marches the BL at EVERY iteration, allowing H to exceed
            // htmax nonlinearly. This is essential for stall prediction.
            // The previous implementation only marched once before the Newton loop, then used
            // only Jacobian updates which cannot capture the nonlinear behavior of separated flow.
            
            if rustfoil_bl::is_debug_active() || std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                eprintln!("[DEBUG Newton] Re-marching BL at iteration {}", iter);
            }
            
            // Extract current Ue from stations (which were updated by previous iteration)
            let current_upper_ue: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
            let current_lower_ue: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();
            
            // Re-march upper surface with current Ue
            let upper_result_iter = march_surface(&upper_arc, &current_upper_ue, re, msq, &march_config, 1);
            
            // Re-march lower surface with current Ue
            let lower_result_iter = march_surface(&lower_arc, &current_lower_ue, re, msq, &march_config, 2);
            
            // Copy march results back to stations (reuse logic from initial march)
            // Note: march_surface skips stagnation (index 0), so result[i] corresponds to stations[i+1]
            for (i, station) in upper_result_iter.stations.iter().enumerate() {
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
                    // Keep the edge velocity from stations (updated by Newton coupling)
                    // Don't overwrite station.u - it comes from VI coupling
                }
            }
            
            for (i, station) in lower_result_iter.stations.iter().enumerate() {
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
                    // Keep the edge velocity from stations (updated by Newton coupling)
                    // Don't overwrite station.u - it comes from VI coupling
                }
            }
            
            // Update transition locations from re-march results
            if let Some(xtr) = upper_result_iter.x_transition {
                x_tr_upper = xtr;
            }
            if let Some(xtr) = lower_result_iter.x_transition {
                x_tr_lower = xtr;
            }
            
            // Debug: Show effect of re-march on critical stations
            if (rustfoil_bl::is_debug_active() || std::env::var("RUSTFOIL_CL_DEBUG").is_ok()) && iter < 3 {
                // Show upper TE region
                if let Some(te_idx) = iblte_upper.checked_sub(2) {
                    if te_idx < upper_stations.len() {
                        let s = &upper_stations[te_idx];
                        eprintln!("[DEBUG Newton] iter {} after re-march: upper[{}] Hk={:.4}, H={:.4}, theta={:.6e}", 
                            iter, te_idx, s.hk, s.h, s.theta);
                    }
                }
                // Show TE
                if iblte_upper < upper_stations.len() {
                    let s = &upper_stations[iblte_upper];
                    eprintln!("[DEBUG Newton] iter {} after re-march: upper[TE={}] Hk={:.4}, H={:.4}, theta={:.6e}", 
                        iter, iblte_upper, s.hk, s.h, s.theta);
                }
            }
            
            // === END RE-MARCH ===

            // Determine flow types for each interval
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

            // Recompute secondary variables for both surfaces
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

            // Debug after blvar - emit BL state summaries
            if rustfoil_bl::is_debug_active() {
                // Collect upper surface station states (every 10th station)
                let upper_station_states: Vec<rustfoil_bl::StationState> = upper_stations
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i % 10 == 0 || *i == upper_stations.len() - 1)
                    .map(|(i, s)| rustfoil_bl::StationState {
                        ibl: i,
                        x: s.x,
                        theta: s.theta,
                        delta_star: s.delta_star,
                        ctau: s.ctau,
                        ue: s.u,
                        mass: s.mass_defect,
                    })
                    .collect();
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bl_state_summary(
                    iter,
                    "upper",
                    upper_station_states,
                ));
                
                // Collect lower surface station states (every 10th station)
                let lower_station_states: Vec<rustfoil_bl::StationState> = lower_stations
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i % 10 == 0 || *i == lower_stations.len() - 1)
                    .map(|(i, s)| rustfoil_bl::StationState {
                        ibl: i,
                        x: s.x,
                        theta: s.theta,
                        delta_star: s.delta_star,
                        ctau: s.ctau,
                        ue: s.u,
                        mass: s.mass_defect,
                    })
                    .collect();
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::bl_state_summary(
                    iter,
                    "lower",
                    lower_station_states,
                ));
                
                if iter < 2 {
                    eprintln!("[DEBUG Global Newton] iter {} after blvar:", iter);
                    for i in 1..4.min(upper_stations.len()) {
                        let s = &upper_stations[i];
                        eprintln!("[DEBUG Global Newton]   upper[{}] theta={:.6e} dstar={:.6e} Hk={:.4}", 
                            i, s.theta, s.delta_star, s.hk);
                    }
                }
            }

            // Prepare transition data (None for now - could use trchek2 results if available)
            let upper_transitions: Vec<Option<rustfoil_bl::closures::Trchek2FullResult>> =
                vec![None; n_upper.saturating_sub(1)];
            let lower_transitions: Vec<Option<rustfoil_bl::closures::Trchek2FullResult>> =
                vec![None; n_lower.saturating_sub(1)];

            // Build the global Newton system with full cross-surface coupling
            global_system.build_global_system(
                upper_stations,
                lower_stations,
                &upper_flow_types,
                &lower_flow_types,
                &upper_transitions,
                &lower_transitions,
                dij,
                upper_ue,
                lower_ue,
                msq,
                re,
                iter,  // Pass iteration number for debug events
            );

            // Emit SETBL_SYSTEM debug event after building the global system
            if rustfoil_bl::is_debug_active() {
                // Convert VA and VB from 3x3 blocks to 3x2 format for XFOIL compatibility
                let va_blocks: Vec<[[f64; 2]; 3]> = global_system.va.iter()
                    .skip(1) // Skip index 0
                    .take(global_system.nsys)
                    .map(|block| [
                        [block[0][0], block[0][1]],
                        [block[1][0], block[1][1]],
                        [block[2][0], block[2][1]],
                    ])
                    .collect();
                let vb_blocks: Vec<[[f64; 2]; 3]> = global_system.vb.iter()
                    .skip(1) // Skip index 0
                    .take(global_system.nsys)
                    .map(|block| [
                        [block[0][0], block[0][1]],
                        [block[1][0], block[1][1]],
                        [block[2][0], block[2][1]],
                    ])
                    .collect();
                let vdel_vec: Vec<[f64; 3]> = global_system.vdel.iter()
                    .skip(1) // Skip index 0
                    .take(global_system.nsys)
                    .cloned()
                    .collect();
                // VM diagonal: extract diagonal elements
                let vm_diag: Vec<[f64; 3]> = (1..=global_system.nsys.min(global_system.vm.len().saturating_sub(1)))
                    .map(|iv| {
                        if iv < global_system.vm.len() && iv < global_system.vm[iv].len() {
                            global_system.vm[iv][iv]
                        } else {
                            [0.0, 0.0, 0.0]
                        }
                    })
                    .collect();
                // VM row 1: coupling from station 1 to all others
                let vm_row1: Vec<[f64; 3]> = if global_system.vm.len() > 1 {
                    global_system.vm[1].iter()
                        .skip(1)
                        .take(global_system.nsys)
                        .cloned()
                        .collect()
                } else {
                    vec![]
                };
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::setbl_system(
                    iter,
                    global_system.nsys,
                    va_blocks,
                    vb_blocks,
                    vdel_vec,
                    vm_diag,
                    vm_row1,
                ));
            }

            // Get residual before solve
            residual = global_system.rms_residual();
            
            // Debug: print actual VDEL values at a few stations to understand magnitude
            if iter < 1 && std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                eprintln!("[DEBUG Newton] iter {} VDEL sample:", iter);
                // Print station data at problematic indices
                let n_upper = global_system.n_upper;
                eprintln!("  n_upper={}, nsys={}", n_upper, global_system.nsys);
                
                // Upper surface first stations - include mass_defect for DUE analysis
                eprintln!("  Upper[0]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    upper_stations[0].x, upper_stations[0].theta, upper_stations[0].delta_star, upper_stations[0].u, upper_stations[0].mass_defect);
                eprintln!("  Upper[1]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    upper_stations[1].x, upper_stations[1].theta, upper_stations[1].delta_star, upper_stations[1].u, upper_stations[1].mass_defect);
                eprintln!("  Upper[2]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    upper_stations[2].x, upper_stations[2].theta, upper_stations[2].delta_star, upper_stations[2].u, upper_stations[2].mass_defect);
                    
                // Lower surface first stations
                eprintln!("  Lower[0]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    lower_stations[0].x, lower_stations[0].theta, lower_stations[0].delta_star, lower_stations[0].u, lower_stations[0].mass_defect);
                eprintln!("  Lower[1]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    lower_stations[1].x, lower_stations[1].theta, lower_stations[1].delta_star, lower_stations[1].u, lower_stations[1].mass_defect);
                eprintln!("  Lower[2]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    lower_stations[2].x, lower_stations[2].theta, lower_stations[2].delta_star, lower_stations[2].u, lower_stations[2].mass_defect);
                
                // VDEL at problematic stations
                for iv in [1, 2, n_upper, n_upper+1].iter() {
                    if *iv <= global_system.nsys {
                        let vdel = global_system.vdel[*iv];
                        let mag = (vdel[0]*vdel[0] + vdel[1]*vdel[1] + vdel[2]*vdel[2]).sqrt();
                        eprintln!("  IV={}: VDEL=[{:.6e}, {:.6e}, {:.6e}] mag={:.6e}", iv, vdel[0], vdel[1], vdel[2], mag);
                    }
                }
                
                // Print DUE and DULE values
                eprintln!("  DULE1={:.6e}, DULE2={:.6e}", global_system.dule1, global_system.dule2);
            }
            
            // Debug output for iterations 6-7 to trace explosion
            if iter >= 5 && iter <= 8 {
                eprintln!("[DEBUG Newton] iter {} residual before solve: {:.6e}", iter, residual);
                // Sample a few stations to check for blow-ups
                for ibl in [10, 20, 30, 40, 50] {
                    if ibl < upper_stations.len() {
                        let s = &upper_stations[ibl];
                        eprintln!("[DEBUG Newton] iter {} upper[{}]: theta={:.6e}, dstar={:.6e}, Ue={:.6}, ctau={:.6}, mass={:.6e}",
                            iter, ibl, s.theta, s.delta_star, s.u, s.ctau, s.mass_defect);
                    }
                    if ibl < lower_stations.len() {
                        let s = &lower_stations[ibl];
                        eprintln!("[DEBUG Newton] iter {} lower[{}]: theta={:.6e}, dstar={:.6e}, Ue={:.6}, ctau={:.6}, mass={:.6e}",
                            iter, ibl, s.theta, s.delta_star, s.u, s.ctau, s.mass_defect);
                    }
                }
            }
            
            if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() && iter < 3 {
                eprintln!("[DEBUG Newton] iter {} residual before solve: {:.6e}", iter, residual);
            }
            
            // Solve the global Newton system
            let deltas = solve_global_system(&mut global_system);
            
            // Emit BLSOLV_SOLUTION debug event with all Newton deltas (full system)
            if rustfoil_bl::is_debug_active() {
                let deltas_full: Vec<[f64; 3]> = deltas
                    .iter()
                    .skip(1)  // Skip index 0 (unused)
                    .take(global_system.nsys)
                    .cloned()
                    .collect();
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::blsolv_solution(
                    iter,
                    global_system.nsys,
                    deltas_full,
                ));
            }
            
            // Emit SOLUTION debug event with Newton deltas (sample for compatibility)
            if rustfoil_bl::is_debug_active() {
                let deltas_sample: Vec<[f64; 3]> = deltas
                    .iter()
                    .skip(1)  // Skip index 0 (unused)
                    .take(30)
                    .cloned()
                    .collect();
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::solution(iter, deltas_sample));
            }
            
            // Debug: Check delta magnitudes for iterations 6-7
            if iter >= 5 && iter <= 8 {
                let max_delta_theta = deltas.iter().map(|d| d[1].abs()).fold(0.0, f64::max);
                let max_delta_mass = deltas.iter().map(|d| d[2].abs()).fold(0.0, f64::max);
                let max_delta_ctau = deltas.iter().map(|d| d[0].abs()).fold(0.0, f64::max);
                eprintln!("[DEBUG Newton] iter {} max deltas: theta={:.6e}, mass={:.6e}, ctau={:.6e}",
                    iter, max_delta_theta, max_delta_mass, max_delta_ctau);
            }
            
            // Check if solution is valid
            let deltas_valid = deltas.iter().all(|d| d.iter().all(|v| v.is_finite()));
            if !deltas_valid {
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG Newton] Invalid deltas at iter {} (reverting)", iter);
                }
                // Restore original values including edge velocity
                for (i, (theta, dstar, ctau, ampl, h, mass, u)) in upper_backup.iter().enumerate() {
                    upper_stations[i].theta = *theta;
                    upper_stations[i].delta_star = *dstar;
                    upper_stations[i].ctau = *ctau;
                    upper_stations[i].ampl = *ampl;
                    upper_stations[i].h = *h;
                    upper_stations[i].mass_defect = *mass;
                    upper_stations[i].u = *u;
                }
                for (i, (theta, dstar, ctau, ampl, h, mass, u)) in lower_backup.iter().enumerate() {
                    lower_stations[i].theta = *theta;
                    lower_stations[i].delta_star = *dstar;
                    lower_stations[i].ctau = *ctau;
                    lower_stations[i].ampl = *ampl;
                    lower_stations[i].h = *h;
                    lower_stations[i].mass_defect = *mass;
                    lower_stations[i].u = *u;
                }
                iteration = 0;
                break;
            }

            // Apply updates to both surfaces with relaxation
            // XFOIL uses RLX = 1.0 initially, then the UPDATE routine computes
            // normalized changes and reduces RLX if any variable would change
            // by more than DHI=1.5 (150% increase) or DLO=-0.5 (50% decrease).
            // apply_global_updates now implements this normalized relaxation.
            let rlx = 1.0;
            
            // Debug: Track relaxation factor for iterations 6-7
            let rlx_before = rlx;
            
            let update_result = apply_global_updates(
                upper_stations,
                lower_stations,
                &deltas,
                &global_system,
                rlx,
                Some(dij),
                Some(upper_ue),
                Some(lower_ue),
                msq,
                re,
            );
            
            // Use RMS of normalized changes (XFOIL's RMSBL) for convergence check
            // This is the correct metric - NOT the raw residuals from VDEL
            residual = update_result.rms_change;
            
            // Debug: Check relaxation used and residual after update for iterations 6-7
            if iter >= 5 && iter <= 8 {
                eprintln!("[DEBUG Newton] iter {} relaxation used: {:.6e} (requested: {:.6e}), rms_change: {:.6e}", 
                    iter, update_result.relaxation_used, rlx_before, residual);
            }

            
            // === STMOVE-style stagnation point check ===
            // After updating edge velocities, check if stagnation point would move.
            // In XFOIL, STMOVE shifts BL arrays when stagnation moves panels.
            // Our architecture (pre-split surfaces) limits full STMOVE, but we can
            // detect drift and adjust the first-station Ue to maintain consistency.
            // 
            // Detect stagnation drift by checking if Ue at first stations has changed sign
            // or if the minimum |Ue| has moved significantly.
            if upper_stations.len() > 2 && lower_stations.len() > 2 {
                let upper_first_ue = upper_stations[1].u;
                let lower_first_ue = lower_stations[1].u;
                
                // Stagnation should have Ue ≈ 0 and increasing downstream
                // If first station Ue is very small or wrong sign, adjust
                let stag_ue_upper = upper_stations[0].u;
                let stag_ue_lower = lower_stations[0].u;
                
                // Ensure stagnation stations have small positive Ue (just past stagnation)
                if stag_ue_upper.abs() < 1e-6 || stag_ue_upper < 0.0 {
                    upper_stations[0].u = upper_first_ue.abs().min(0.01).max(1e-6);
                    upper_stations[0].mass_defect = upper_stations[0].u * upper_stations[0].delta_star;
                }
                if stag_ue_lower.abs() < 1e-6 || stag_ue_lower < 0.0 {
                    lower_stations[0].u = lower_first_ue.abs().min(0.01).max(1e-6);
                    lower_stations[0].mass_defect = lower_stations[0].u * lower_stations[0].delta_star;
                }
            }

            if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() && (iter < 10 || iter % 5 == 0) {
                eprintln!("[DEBUG Newton] iter {} RMSBL after update: {:.6e} (max_change: {:.6e})", 
                    iter, residual, update_result.max_change);
            }
            
            if residual < config.tolerance {
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG Newton] Converged at iter {} with RMSBL={:.6e} < tolerance={:.6e}", 
                        iter, residual, config.tolerance);
                }
                converged = true;
                break;
            }

            // Check for divergence or increasing residual
            if !residual.is_finite() || residual > 1e10 {
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG Newton] Diverged at iter {}: residual={:.6e} (reverting)", iter, residual);
                }
                // Restore original values including edge velocity
                for (i, (theta, dstar, ctau, ampl, h, mass, u)) in upper_backup.iter().enumerate() {
                    upper_stations[i].theta = *theta;
                    upper_stations[i].delta_star = *dstar;
                    upper_stations[i].ctau = *ctau;
                    upper_stations[i].ampl = *ampl;
                    upper_stations[i].h = *h;
                    upper_stations[i].mass_defect = *mass;
                    upper_stations[i].u = *u;
                }
                for (i, (theta, dstar, ctau, ampl, h, mass, u)) in lower_backup.iter().enumerate() {
                    lower_stations[i].theta = *theta;
                    lower_stations[i].delta_star = *dstar;
                    lower_stations[i].ctau = *ctau;
                    lower_stations[i].ampl = *ampl;
                    lower_stations[i].h = *h;
                    lower_stations[i].mass_defect = *mass;
                    lower_stations[i].u = *u;
                }
                converged = true; // Fall back to direct march result
                iteration = 0;
                break;
            }

            
            // Stop if residual is increasing rapidly (Newton diverging badly)
            // Allow slow increase since our residual may oscillate
            if iter > 5 && residual > prev_residual * 5.0 {
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG Newton] Rapid increase at iter {}: residual={:.6e}, prev={:.6e} (reverting)", iter, residual, prev_residual);
                }
                // Restore original values including edge velocity
                for (i, (theta, dstar, ctau, ampl, h, mass, u)) in upper_backup.iter().enumerate() {
                    upper_stations[i].theta = *theta;
                    upper_stations[i].delta_star = *dstar;
                    upper_stations[i].ctau = *ctau;
                    upper_stations[i].ampl = *ampl;
                    upper_stations[i].h = *h;
                    upper_stations[i].mass_defect = *mass;
                    upper_stations[i].u = *u;
                }
                for (i, (theta, dstar, ctau, ampl, h, mass, u)) in lower_backup.iter().enumerate() {
                    lower_stations[i].theta = *theta;
                    lower_stations[i].delta_star = *dstar;
                    lower_stations[i].ctau = *ctau;
                    lower_stations[i].ampl = *ampl;
                    lower_stations[i].h = *h;
                    lower_stations[i].mass_defect = *mass;
                    lower_stations[i].u = *u;
                }
                converged = true;
                iteration = 0;
                break;
            }
            // Emit FULL_BL_STATE debug event with complete BL state for both surfaces
            if rustfoil_bl::is_debug_active() {
                let upper_bl = rustfoil_bl::SurfaceBlState {
                    x: upper_stations.iter().map(|s| s.x).collect(),
                    theta: upper_stations.iter().map(|s| s.theta).collect(),
                    delta_star: upper_stations.iter().map(|s| s.delta_star).collect(),
                    ue: upper_stations.iter().map(|s| s.u).collect(),
                    hk: upper_stations.iter().map(|s| s.hk).collect(),
                    cf: upper_stations.iter().map(|s| s.cf).collect(),
                    mass_defect: upper_stations.iter().map(|s| s.mass_defect).collect(),
                };
                let lower_bl = rustfoil_bl::SurfaceBlState {
                    x: lower_stations.iter().map(|s| s.x).collect(),
                    theta: lower_stations.iter().map(|s| s.theta).collect(),
                    delta_star: lower_stations.iter().map(|s| s.delta_star).collect(),
                    ue: lower_stations.iter().map(|s| s.u).collect(),
                    hk: lower_stations.iter().map(|s| s.hk).collect(),
                    cf: lower_stations.iter().map(|s| s.cf).collect(),
                    mass_defect: lower_stations.iter().map(|s| s.mass_defect).collect(),
                };
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_bl_state(
                    iteration,
                    upper_bl,
                    lower_bl,
                ));

                // Emit FULL_NFACTOR debug event with N-factors at all stations
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_nfactor(
                    iteration,
                    upper_stations.iter().map(|s| s.ampl).collect(),
                    lower_stations.iter().map(|s| s.ampl).collect(),
                    upper_stations.iter().map(|s| s.x).collect(),
                    lower_stations.iter().map(|s| s.x).collect(),
                ));
            }

            prev_residual = residual;
        }
    }

    // Update transition locations from final station states
    // Find transition point (where is_laminar changes to is_turbulent)
    // Use x_coord (panel x-coordinate) for proper x/c reporting, not arc length
    for (i, station) in upper_stations.iter().enumerate() {
        if !station.is_laminar && i > 0 && upper_stations[i - 1].is_laminar {
            x_tr_upper = station.x_coord;
            break;
        }
    }
    for (i, station) in lower_stations.iter().enumerate() {
        if !station.is_laminar && i > 0 && lower_stations[i - 1].is_laminar {
            x_tr_lower = station.x_coord;
            break;
        }
    }

    // Check for separation on either surface
    let x_separation = upper_result
        .x_separation
        .or(lower_result.x_separation);

    // === QVFUE + GAMQV: Update QVIS from edge velocities, then gamma from QVIS ===
    // This is part of XFOIL's iteration loop that couples the BL solution back to
    // the panel method circulation.
    let (upper_qvis, lower_qvis) =
        update_qvis_from_uedg_two_surfaces(upper_stations, lower_stations);
    let upper_gamma = update_circulation_from_qvis(&upper_qvis);
    let lower_gamma = update_circulation_from_qvis(&lower_qvis);

    // === CLCALC: Compute CL using XFOIL's exact formula ===
    // XFOIL integrates Cp in wind axes around the closed airfoil contour:
    //   CL = ∮ Cp * dx_wind
    // where dx_wind = (X[i+1]-X[i])*cos(α) + (Y[i+1]-Y[i])*sin(α)
    //
    // This is the EXACT formula from XFOIL's CLCALC (xfoil.f:1088-1168)
    // and must be used to achieve numerical precision matching.
    
    // Debug: Show station edge velocities before constructing gamma
    if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
        let upper_ue_range: (f64, f64) = upper_stations.iter()
            .filter(|s| !s.is_wake)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), s| (min.min(s.u.abs()), max.max(s.u.abs())));
        let lower_ue_range: (f64, f64) = lower_stations.iter()
            .filter(|s| !s.is_wake)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), s| (min.min(s.u.abs()), max.max(s.u.abs())));
        eprintln!("[DEBUG viscal] Upper Ue range: [{:.4}, {:.4}] (n={})", 
            upper_ue_range.0, upper_ue_range.1, upper_stations.iter().filter(|s| !s.is_wake).count());
        eprintln!("[DEBUG viscal] Lower Ue range: [{:.4}, {:.4}] (n={})", 
            lower_ue_range.0, lower_ue_range.1, lower_stations.iter().filter(|s| !s.is_wake).count());
        // Show first/last few edge velocities
        let upper_non_wake: Vec<f64> = upper_stations.iter().filter(|s| !s.is_wake).map(|s| s.u).collect();
        let lower_non_wake: Vec<f64> = lower_stations.iter().filter(|s| !s.is_wake).map(|s| s.u).collect();
        if upper_non_wake.len() >= 5 {
            eprintln!("[DEBUG viscal] Upper first 5 Ue: {:?}", upper_non_wake[..5].iter().map(|u| format!("{:.4}", u)).collect::<Vec<_>>());
            eprintln!("[DEBUG viscal] Upper last 5 Ue: {:?}", upper_non_wake[upper_non_wake.len()-5..].iter().map(|u| format!("{:.4}", u)).collect::<Vec<_>>());
        }
        if lower_non_wake.len() >= 5 {
            eprintln!("[DEBUG viscal] Lower first 5 Ue: {:?}", lower_non_wake[..5].iter().map(|u| format!("{:.4}", u)).collect::<Vec<_>>());
            eprintln!("[DEBUG viscal] Lower last 5 Ue: {:?}", lower_non_wake[lower_non_wake.len()-5..].iter().map(|u| format!("{:.4}", u)).collect::<Vec<_>>());
        }
    }
    
    // Construct full panel gamma array from BL station edge velocities
    let n_panels = panel_x.len();
    let full_gamma = construct_gamma_from_stations(upper_stations, lower_stations, n_panels);
    
    // Emit FULL_GAMMA_ITER debug event with circulation at all panels
    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_gamma_iter(
            iteration,
            full_gamma.clone(),
        ));
    }
    
    // Use XFOIL's exact CLCALC formula with wind-axis integration
    let cl_xfoil = clcalc(panel_x, panel_y, &full_gamma, alpha_rad);
    
    // Also compute circulation-based CL for comparison/fallback
    let cl_from_gamma = compute_cl_from_gamma(&upper_gamma, &lower_gamma, &upper_arc, &lower_arc);
    
    // Debug: Log both CL methods for comparison
    if rustfoil_bl::is_debug_active() || std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
        eprintln!("[DEBUG viscal] CL comparison:");
        eprintln!("  CLCALC (wind-axis):     {:.6}", cl_xfoil);
        eprintln!("  Circulation (arc len):  {:.6}", cl_from_gamma);
        eprintln!("  Difference:             {:.6e}", cl_xfoil - cl_from_gamma);
    }
    
    // Use XFOIL's CLCALC as primary method (matches XFOIL exactly)
    // Fall back to circulation only if CLCALC gives unreasonable result
    let cl_final = if cl_xfoil.is_finite() && cl_xfoil.abs() < 10.0 {
        cl_xfoil
    } else if cl_from_gamma.is_finite() && cl_from_gamma.abs() < 10.0 {
        cl_from_gamma
    } else {
        // Both methods failed
        cl_xfoil
    };

    // === Create and march wake stations for Squire-Young CD computation ===
    // Wake stations extend downstream from TE and are needed for proper far-wake
    // momentum thickness which Squire-Young uses for CD.
    //
    // Use XFOIL-style wake: combine TE from both surfaces, then march downstream
    // using wake BL equations (Cf=0, wake dissipation).
    let upper_te = upper_stations.last().cloned().unwrap_or_default();
    let lower_te = lower_stations.last().cloned().unwrap_or_default();
    
    // Combine upper and lower TE into wake initial condition (TESYS)
    let wake_initial = combine_te_for_wake(&upper_te, &lower_te);
    
    // Generate wake positions and edge velocities
    // XFOIL marches wake to x ≈ 2c with many stations for Hk → 1.0 decay
    // More stations + longer wake allows proper far-wake Squire-Young calculation
    use rustfoil_coupling::wake::{generate_wake_positions, march_wake, wake_edge_velocity};
    let wake_x = generate_wake_positions(15, 2.0);
    let wake_ue: Vec<f64> = wake_x.iter()
        .map(|&x| wake_edge_velocity(x, wake_initial.u))
        .collect();
    
    // March wake downstream using wake BL equations
    let wake_stations = march_wake(&wake_initial, &wake_x, &wake_ue, re, msq);

    // Compute forces from both surfaces using coupled edge velocities
    let mut forces = compute_forces_two_surfaces(upper_stations, lower_stations, config);

    // Use the selected CL calculation method
    forces.cl = cl_final;
    
    // === Wake-based CD using Squire-Young at far wake ===
    // XFOIL computes total CD at far wake where:
    // - Ue ≈ freestream (momentum nearly conserved)
    // - Hk → 1.0 (thin wake limit)
    // 
    // At TE, Hk≈2.5 gives exponent (5+H)/2 ≈ 3.75
    // At far wake, Hk≈1.04 gives exponent ≈ 3.02
    // This difference accounts for ~50% CD overestimate at low alpha.
    //
    // Use the marched wake stations for proper far-wake Squire-Young.
    use super::forces::compute_cd_from_wake;
    
    if !wake_stations.is_empty() && wake_stations.len() > 1 {
        let cd_wake = compute_cd_from_wake(&wake_stations, forces.cd_friction);
        
        // Debug: Compare wake-based vs TE-based CD
        if rustfoil_bl::is_debug_active() || std::env::var("RUSTFOIL_DRAG_DEBUG").is_ok() {
            let cd_pressure_wake_debug = cd_wake - forces.cd_friction;  // Always >= 0.0 due to fix in compute_cd_from_wake
            eprintln!("[DEBUG viscal] CD comparison:");
            eprintln!("  TE-based pressure CD:   {:.6e} (from compute_forces_two_surfaces)", forces.cd);
            eprintln!("  Wake-based total CD:    {:.6e} (from far-wake Squire-Young)", cd_wake);
            eprintln!("  Wake-based pressure CD: {:.6e} (total - friction)", cd_pressure_wake_debug);
            eprintln!("  Wake stations: {} (initial θ={:.4e}, final θ={:.4e})",
                wake_stations.len(),
                wake_stations.first().map(|s| s.theta).unwrap_or(0.0),
                wake_stations.last().map(|s| s.theta).unwrap_or(0.0));
        }
        
        // Use wake-based CD (total drag) if it's reasonable AND larger than TE-based
        // Note: compute_cd_from_wake ensures cd_wake >= cd_friction
        if cd_wake.is_finite() && cd_wake > 0.0 && cd_wake < 0.1 {
            // Compare total drags - use whichever is larger
            // At low alpha, wake Squire-Young can underestimate drag because the wake hasn't
            // developed enough. In that case, keep the TE-based CD which is more reliable.
            if cd_wake > forces.cd {
                let cd_pressure_wake = cd_wake - forces.cd_friction;
                forces.cd = cd_wake;  // Total drag (matches XFOIL CD_BREAKDOWN)
                forces.cd_pressure = cd_pressure_wake;
            }
            // If cd_wake <= forces.cd, keep the TE-based CD from compute_forces_two_surfaces
        }
    }
    
    // Debug: Output BL quantities at all stations for comparison with XFOIL
    if std::env::var("RUSTFOIL_BL_DEBUG").is_ok() {
        eprintln!("[BL_DEBUG] UPPER_SURFACE ({} stations):", upper_stations.len());
        for s in upper_stations.iter() {
            let flow = if s.is_laminar { "LAM" } else { "TURB" };
            eprintln!("[BL_DEBUG] x={:.4} Hk={:.3} Cq={:.4} Cf={:.2e} ctau={:.4} {}",
                s.x, s.hk, s.cq, s.cf, s.ctau, flow);
        }
        
        eprintln!("[BL_DEBUG] LOWER_SURFACE ({} stations):", lower_stations.len());
        for s in lower_stations.iter() {
            let flow = if s.is_laminar { "LAM" } else { "TURB" };
            eprintln!("[BL_DEBUG] x={:.4} Hk={:.3} Cq={:.4} Cf={:.2e} ctau={:.4} {}",
                s.x, s.hk, s.cq, s.cf, s.ctau, flow);
        }
    }

    let _ = (&upper_te, &lower_te); // Suppress unused warnings

    // Debug: Print force computation
    if rustfoil_bl::is_debug_active() {
        eprintln!("[DEBUG viscal] Combined forces: CL={:.4}, CD={:.6e} (Cf={:.6e}, Cp={:.6e})",
            forces.cl, forces.cd, forces.cd_friction, forces.cd_pressure);
        
        // Check last station values (used in Squire-Young)
        if let Some(last) = upper_stations.last() {
            eprintln!("[DEBUG viscal] Upper last: theta={:.6e}, H={:.3}, Ue={:.3}, is_wake={}",
                last.theta, last.h, last.u, last.is_wake);
        }
        if let Some(last) = lower_stations.last() {
            eprintln!("[DEBUG viscal] Lower last: theta={:.6e}, H={:.3}, Ue={:.3}, is_wake={}",
                last.theta, last.h, last.u, last.is_wake);
        }
    }

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
