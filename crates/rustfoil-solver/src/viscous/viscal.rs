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
use rustfoil_coupling::march::{march_fixed_ue, MarchConfig, MarchResult};
use rustfoil_coupling::global_newton::{
    apply_global_updates_from_view, preview_global_update_ue_from_view, solve_global_system, GlobalNewtonSystem,
};
use rustfoil_coupling::newton::BlNewtonSystem;
use rustfoil_coupling::newton_state::{CanonicalNewtonRow, CanonicalNewtonStateView};
use rustfoil_coupling::solve::solve_coupled_system;
// STMOVE module for stagnation point finding and coupling
use rustfoil_coupling::stmove::{
    find_stagnation_with_derivs, StagnationResult,
};
use rustfoil_coupling::update::{update_xfoil_style, UpdateConfig};

use super::circulation::{
    refresh_canonical_panel_arrays, update_circulation_from_qvis, update_qvis_from_uedg,
};
use super::config::{OperatingMode, ViscousSolverConfig};
use super::forces::{compute_forces, compute_forces_from_canonical_state, compute_panel_forces_from_gamma};
use super::setup::{compute_arc_lengths, extract_surface_xfoil};
use super::state::{TransitionalStmoveResult, XfoilLikeViscousState, XfoilSurface};
use crate::{SolverError, SolverResult};

fn surface_state(stations: &[BlStation]) -> rustfoil_bl::SurfaceBlState {
    rustfoil_bl::SurfaceBlState {
        x: stations.iter().map(|station| station.x).collect(),
        theta: stations.iter().map(|station| station.theta).collect(),
        delta_star: stations.iter().map(|station| station.delta_star).collect(),
        ue: stations.iter().map(|station| station.u).collect(),
        hk: stations.iter().map(|station| station.hk).collect(),
        cf: stations.iter().map(|station| station.cf).collect(),
        mass_defect: stations.iter().map(|station| station.mass_defect).collect(),
    }
}

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
    compute_panel_forces_from_gamma(panel_x, panel_y, gamma, alpha).0
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
#[derive(Debug, Clone, Copy)]
enum GammaWriteOrder {
    PreserveFirstWriter,
    LowerOverwritesShared,
}

fn build_panel_gamma_from_stations(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    n_panels: usize,
    order: GammaWriteOrder,
) -> Vec<f64> {
    let mut gamma = vec![0.0; n_panels];
    let mut filled = vec![false; n_panels];
    
    // Upper surface: VTI = +1
    // Upper stations go from stagnation (high panel index) toward TE (low panel index)
    // XFOIL: QVIS = VTI * UEDG, GAM = QVIS
    // So: gamma = +1 * Ue = +Ue
    // Note: station.u should be positive (magnitude), but we preserve sign for correctness
    for station in upper_stations.iter().skip(1) {
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
    for station in lower_stations.iter().skip(1) {
        let idx = station.panel_idx;
        if idx >= n_panels {
            continue;
        }
        let lower_writes = match order {
            GammaWriteOrder::PreserveFirstWriter => !filled[idx],
            GammaWriteOrder::LowerOverwritesShared => true,
        };
        if !station.is_wake && lower_writes {
            // VTI = -1: gamma = -Ue.
            // In exact QVFUE ordering, the lower surface is written after the upper
            // surface and therefore overwrites shared indices.
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

fn build_canonical_state(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    n_panel_nodes: usize,
    stagnation: Option<StagnationResult>,
) -> XfoilLikeViscousState {
    let mut state =
        XfoilLikeViscousState::from_station_views(upper_stations, lower_stations, n_panel_nodes);
    if let Some(stag) = stagnation {
        state.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
    }
    state
}

fn refresh_canonical_state_from_station_views(
    state: &mut XfoilLikeViscousState,
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    n_panel_nodes: usize,
    stagnation: Option<StagnationResult>,
) {
    state.sync_from_station_views(upper_stations, lower_stations, n_panel_nodes);
    if let Some(stag) = stagnation {
        state.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
    }
    refresh_canonical_panel_arrays(state);
}

fn build_canonical_newton_view(
    state: &XfoilLikeViscousState,
    upper_ue_operating: Vec<f64>,
    lower_ue_operating: Vec<f64>,
    ante: f64,
    sst_go: f64,
    sst_gp: f64,
) -> CanonicalNewtonStateView {
    let to_row = |row: &super::state::CanonicalBlRow| CanonicalNewtonRow {
        x: row.x,
        x_coord: row.x_coord,
        panel_idx: row.panel_idx,
        uedg: row.uedg,
        uinv: row.uinv,
        uinv_a: row.uinv_a,
        theta: row.theta,
        dstr: row.dstr,
        ctau_or_ampl: row.ctau_or_ampl,
        ctau: row.ctau,
        ampl: row.ampl,
        mass: row.mass,
        h: row.h,
        hk: row.hk,
        hs: row.hs,
        hc: row.hc,
        r_theta: row.r_theta,
        cf: row.cf,
        cd: row.cd,
        us: row.us,
        cq: row.cq,
        de: row.de,
        dw: row.dw,
        is_laminar: row.is_laminar,
        is_turbulent: row.is_turbulent,
        is_wake: row.is_wake,
        derivs: row.derivs.clone(),
    };

    CanonicalNewtonStateView {
        upper_rows: state.upper_rows.iter().map(to_row).collect(),
        lower_rows: state.lower_rows.iter().map(to_row).collect(),
        upper_flow_types: state.flow_type_view(XfoilSurface::Upper),
        lower_flow_types: state.flow_type_view(XfoilSurface::Lower),
        upper_ue_current: state.uedg_view(XfoilSurface::Upper),
        lower_ue_current: state.uedg_view(XfoilSurface::Lower),
        upper_ue_inviscid: state.uinv_view(XfoilSurface::Upper),
        lower_ue_inviscid: state.uinv_view(XfoilSurface::Lower),
        upper_ue_from_mass: Vec::new(),
        lower_ue_from_mass: Vec::new(),
        upper_ue_operating,
        lower_ue_operating,
        sst_go,
        sst_gp,
        ante,
    }
}

fn transition_arc_to_chord_on_surface(
    stagnation: StagnationResult,
    full_arc: &[f64],
    panel_x: &[f64],
    panel_y: &[f64],
    ue_inviscid_full: &[f64],
    is_upper: bool,
    x_transition: Option<f64>,
) -> Option<f64> {
    let xt = x_transition?;
    if full_arc.is_empty() || panel_x.is_empty() || panel_y.is_empty() || ue_inviscid_full.is_empty() {
        return Some(xt);
    }

    let ist_next = (stagnation.ist + 1).min(ue_inviscid_full.len().saturating_sub(1));
    let ue_stag = if ist_next > stagnation.ist && full_arc[ist_next] != full_arc[stagnation.ist] {
        let frac = (stagnation.sst - full_arc[stagnation.ist])
            / (full_arc[ist_next] - full_arc[stagnation.ist]);
        ue_inviscid_full[stagnation.ist]
            + frac * (ue_inviscid_full[ist_next] - ue_inviscid_full[stagnation.ist])
    } else {
        ue_inviscid_full[stagnation.ist]
    };

    let (surface_arc, surface_x, _, _) = extract_surface_xfoil(
        stagnation.ist,
        stagnation.sst,
        ue_stag,
        full_arc,
        panel_x,
        panel_y,
        ue_inviscid_full,
        is_upper,
    );
    if surface_arc.is_empty() || surface_x.is_empty() {
        return Some(xt);
    }
    if xt <= surface_arc[0] {
        return Some(surface_x[0]);
    }

    for i in 1..surface_arc.len() {
        if xt <= surface_arc[i] {
            let ds = surface_arc[i] - surface_arc[i - 1];
            if ds.abs() <= 1.0e-12 {
                return Some(surface_x[i]);
            }
            let frac = ((xt - surface_arc[i - 1]) / ds).clamp(0.0, 1.0);
            return Some(surface_x[i - 1] + frac * (surface_x[i] - surface_x[i - 1]));
        }
    }

    surface_x.last().copied()
}

fn operating_sensitivity_for_mode(
    mode: OperatingMode,
    upper_alpha: &[f64],
    lower_alpha: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    match mode {
        OperatingMode::PrescribedAlpha => (
            vec![0.0; upper_alpha.len()],
            vec![0.0; lower_alpha.len()],
        ),
        OperatingMode::PrescribedCl { .. } => (upper_alpha.to_vec(), lower_alpha.to_vec()),
    }
}

fn build_panel_values_from_surface_values(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    upper_values: &[f64],
    lower_values: &[f64],
    n_panel_nodes: usize,
) -> Vec<f64> {
    let mut panel_values = vec![0.0; n_panel_nodes];

    for (ibl, station) in upper_stations.iter().enumerate().skip(1) {
        if station.is_wake || station.panel_idx >= n_panel_nodes {
            continue;
        }
        panel_values[station.panel_idx] = upper_values.get(ibl).copied().unwrap_or(station.u);
    }

    for (ibl, station) in lower_stations.iter().enumerate().skip(1) {
        if station.is_wake || station.panel_idx >= n_panel_nodes {
            continue;
        }
        panel_values[station.panel_idx] = -lower_values.get(ibl).copied().unwrap_or(station.u);
    }

    if let Some(last) = panel_values.last_mut() {
        *last = 0.0;
    }

    panel_values
}

fn build_panel_values_from_row_values(
    upper_rows: &[CanonicalNewtonRow],
    lower_rows: &[CanonicalNewtonRow],
    upper_values: &[f64],
    lower_values: &[f64],
    n_panel_nodes: usize,
) -> Vec<f64> {
    let mut panel_values = vec![0.0; n_panel_nodes];

    for (ibl, row) in upper_rows.iter().enumerate().skip(1) {
        if row.is_wake || row.panel_idx >= n_panel_nodes {
            continue;
        }
        panel_values[row.panel_idx] = upper_values.get(ibl).copied().unwrap_or(row.uedg);
    }

    for (ibl, row) in lower_rows.iter().enumerate().skip(1) {
        if row.is_wake || row.panel_idx >= n_panel_nodes {
            continue;
        }
        panel_values[row.panel_idx] = -lower_values.get(ibl).copied().unwrap_or(row.uedg);
    }

    if let Some(last) = panel_values.last_mut() {
        *last = 0.0;
    }

    panel_values
}

fn compute_trial_cl_terms(
    panel_x: &[f64],
    panel_y: &[f64],
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    upper_u_new: &[f64],
    lower_u_new: &[f64],
    upper_u_ac: &[f64],
    lower_u_ac: &[f64],
    alpha_rad: f64,
) -> (f64, f64, f64) {
    let qnew = build_panel_values_from_surface_values(
        upper_stations,
        lower_stations,
        upper_u_new,
        lower_u_new,
        panel_x.len(),
    );
    let qac = build_panel_values_from_surface_values(
        upper_stations,
        lower_stations,
        upper_u_ac,
        lower_u_ac,
        panel_x.len(),
    );
    let cl_new = clcalc(panel_x, panel_y, &qnew, alpha_rad);

    let eps = 1.0e-6;
    let cl_alpha = (clcalc(panel_x, panel_y, &qnew, alpha_rad + eps)
        - clcalc(panel_x, panel_y, &qnew, alpha_rad - eps))
        / (2.0 * eps);

    let cl_ac = if qac.iter().any(|value| value.abs() > 0.0) {
        let qnew_operating: Vec<f64> = qnew
            .iter()
            .zip(qac.iter())
            .map(|(base, operating)| base + eps * operating)
            .collect();
        (clcalc(panel_x, panel_y, &qnew_operating, alpha_rad) - cl_new) / eps
    } else {
        0.0
    };

    (cl_new, cl_alpha, cl_ac)
}

fn compute_trial_cl_terms_from_view(
    panel_x: &[f64],
    panel_y: &[f64],
    state: &CanonicalNewtonStateView,
    upper_u_new: &[f64],
    lower_u_new: &[f64],
    upper_u_ac: &[f64],
    lower_u_ac: &[f64],
    alpha_rad: f64,
) -> (f64, f64, f64) {
    let qnew = build_panel_values_from_row_values(
        &state.upper_rows,
        &state.lower_rows,
        upper_u_new,
        lower_u_new,
        panel_x.len(),
    );
    let qac = build_panel_values_from_row_values(
        &state.upper_rows,
        &state.lower_rows,
        upper_u_ac,
        lower_u_ac,
        panel_x.len(),
    );
    let cl_new = clcalc(panel_x, panel_y, &qnew, alpha_rad);

    let eps = 1.0e-6;
    let cl_alpha = (clcalc(panel_x, panel_y, &qnew, alpha_rad + eps)
        - clcalc(panel_x, panel_y, &qnew, alpha_rad - eps))
        / (2.0 * eps);

    let cl_ac = if qac.iter().any(|value| value.abs() > 0.0) {
        let qnew_operating: Vec<f64> = qnew
            .iter()
            .zip(qac.iter())
            .map(|(base, operating)| base + eps * operating)
            .collect();
        (clcalc(panel_x, panel_y, &qnew_operating, alpha_rad) - cl_new) / eps
    } else {
        0.0
    };

    (cl_new, cl_alpha, cl_ac)
}

fn compute_operating_correction(
    mode: OperatingMode,
    cl_new: f64,
    cl_alpha: f64,
    cl_ac: f64,
) -> f64 {
    match mode {
        OperatingMode::PrescribedAlpha => 0.0,
        OperatingMode::PrescribedCl { target_cl } => {
            let denominator = -(cl_ac + cl_alpha);
            if denominator.abs() < 1.0e-12 {
                0.0
            } else {
                const DALMAX: f64 = 0.5_f64.to_radians();
                const DALMIN: f64 = -0.5_f64.to_radians();
                ((cl_new - target_cl) / denominator).clamp(DALMIN, DALMAX)
            }
        }
    }
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
                    d[2], // mass-defect change (converted to delta_star during update)
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

#[allow(dead_code)]
fn refresh_transition_state_from_march(stations: &mut [BlStation], march_result: &MarchResult) {
    let offset = if stations.first().map_or(false, |s| s.x < 1e-6) {
        1
    } else {
        0
    };

    for (i, station) in march_result.stations.iter().enumerate() {
        let target = i + offset;
        if target < stations.len() {
            stations[target].theta = station.theta;
            stations[target].delta_star = station.delta_star;
            stations[target].h = station.h;
            stations[target].hk = station.hk;
            stations[target].ampl = station.ampl;
            stations[target].is_laminar = station.is_laminar;
            stations[target].is_turbulent = station.is_turbulent;
            stations[target].ctau = station.ctau;
            stations[target].cf = station.cf;
            stations[target].mass_defect = station.mass_defect;
            stations[target].r_theta = station.r_theta;
        }
    }
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
    upper_stations: &mut Vec<BlStation>,
    lower_stations: &mut Vec<BlStation>,
    upper_ue: &[f64],
    lower_ue: &[f64],
    dij: &DMatrix<f64>,
    config: &ViscousSolverConfig,
    alpha_rad: f64,
    panel_x: &[f64],
    panel_y: &[f64],
    ue_inviscid_full: &[f64],
    ue_inviscid_alpha_full: &[f64],
) -> SolverResult<ViscousResult> {
    if config.reynolds <= 0.0 {
        return Err(SolverError::InvalidReynolds);
    }

    let re = config.reynolds;
    let msq = config.msq();
    let operating_mode = config.operating_mode;

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
    let full_arc = compute_arc_lengths(panel_x, panel_y);
    let initial_stagnation = find_stagnation_with_derivs(ue_inviscid_full, &full_arc);
    let mut canonical_state = build_canonical_state(
        upper_stations,
        lower_stations,
        panel_x.len(),
        initial_stagnation,
    );
    canonical_state.set_operating_mode(operating_mode);
    canonical_state.set_panel_inviscid_arrays(ue_inviscid_full, ue_inviscid_alpha_full);
    refresh_canonical_panel_arrays(&mut canonical_state);


    // March upper surface (side 1) through the canonical state.
    let upper_result =
        canonical_state.march_surface(XfoilSurface::Upper, re, msq, &march_config);

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

    // March lower surface (side 2) through the canonical state.
    let lower_result =
        canonical_state.march_surface(XfoilSurface::Lower, re, msq, &march_config);

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

    canonical_state.write_back_station_views(upper_stations, lower_stations);

    canonical_state.set_transition_metadata(
        XfoilSurface::Upper,
        upper_result.transition_index,
        upper_result.x_transition,
    );
    canonical_state.set_transition_metadata(
        XfoilSurface::Lower,
        lower_result.transition_index,
        lower_result.x_transition,
    );
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

    // Enable Newton with detailed debug output for comparison with XFOIL
    let can_run_newton = config.max_iterations > 0 && n_upper >= 3 && n_lower >= 3;

    // Debug: track Newton execution
    if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
        eprintln!("[DEBUG viscal] Newton check: can_run={}, max_iter={}, n_upper={}, n_lower={}", 
            can_run_newton, config.max_iterations, n_upper, n_lower);
    }

    if can_run_newton {
        converged = false;

        // === Compute stagnation point derivatives (SST_GO, SST_GP) ===
        // Build a two-point local STFIND proxy around the stagnation panel:
        //   gamma_upstream   = +Ue on the upper first post-stagnation station
        //   gamma_downstream = -Ue on the lower first post-stagnation station
        // and use the actual first-segment arc lengths on each surface. This follows
        // the STFIND interpolation formulas instead of the previous midpoint heuristic.
        let (sst_go, sst_gp) = {
            let ue_upper_1 = upper_ue.get(1).copied().unwrap_or(0.1).abs();
            let ue_lower_1 = lower_ue.get(1).copied().unwrap_or(0.1).abs();
            let ds_upper = upper_arc.get(1).copied().unwrap_or(0.01).abs().max(1.0e-8);
            let ds_lower = lower_arc.get(1).copied().unwrap_or(0.01).abs().max(1.0e-8);

            let gamma_proxy = [ue_upper_1, -ue_lower_1];
            let s_proxy = [-ds_lower, ds_upper];

            find_stagnation_with_derivs(&gamma_proxy, &s_proxy)
                .map(|result| (result.sst_go, result.sst_gp))
                .unwrap_or((0.0, 0.0))
        };

        if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
            eprintln!("[DEBUG viscal] Stagnation derivatives: SST_GO={:.6e}, SST_GP={:.6e}", sst_go, sst_gp);
        }

        // Run Newton iteration with full global coupling
        let max_newton_iter = config.max_iterations.min(30);
        let mut prev_residual = f64::MAX;

        for iter in 0..max_newton_iter {
            iteration = iter + 1;
            let n_upper = canonical_state.upper_rows.len();
            let n_lower = canonical_state.lower_rows.len();
            let iblte_upper = canonical_state.iblte_upper;
            let iblte_lower = canonical_state.iblte_lower;
            let mut global_system =
                GlobalNewtonSystem::new(n_upper, n_lower, iblte_upper, iblte_lower);

            // ANTE: trailing edge gap thickness (XFOIL's WGAP(1) = DWTE).
            // For sharp TE airfoils this is zero; for blunt TE it's the gap.
            // Computed from the distance between upper and lower TE panel nodes.
            {
                let ist_upper_te = canonical_state.upper_rows
                    .get(iblte_upper)
                    .map(|row| row.panel_idx)
                    .unwrap_or(0);
                let ist_lower_te = canonical_state.lower_rows
                    .get(iblte_lower)
                    .map(|row| row.panel_idx)
                    .unwrap_or(0);
                if ist_upper_te < panel_x.len() && ist_lower_te < panel_x.len() {
                    let dx = panel_x[ist_upper_te] - panel_x[ist_lower_te];
                    let dy = panel_y[ist_upper_te] - panel_y[ist_lower_te];
                    global_system.ante = (dx * dx + dy * dy).sqrt();
                }
            }

            // Compute WGAP per wake station (XFOIL xpanel.f:1586-1598)
            // and assign to station.dw for use in BLPRV/BLDIF/DSLIM.
            {
                let ante = global_system.ante;
                let sharp = ante < 1e-10;
                let telrat = 2.5_f64;
                if !sharp && iblte_lower < canonical_state.lower_rows.len() {
                    let xssi_te = canonical_state.lower_rows[iblte_lower].x;
                    let crosp = 0.0_f64; // sharp TE approximation
                    let dwdxte = crosp / (1.0 - crosp * crosp).sqrt().max(1e-20);
                    let dwdxte = dwdxte.clamp(-3.0 / telrat, 3.0 / telrat);
                    let aa = 3.0 + telrat * dwdxte;
                    let bb = -2.0 - telrat * dwdxte;
                    for iw in 1..=(canonical_state.lower_rows.len() - iblte_lower - 1) {
                        let ibl = iblte_lower + iw;
                        let zn = 1.0 - (canonical_state.lower_rows[ibl].x - xssi_te) / (telrat * ante);
                        canonical_state.lower_rows[ibl].dw = if zn >= 0.0 {
                            ante * (aa + bb * zn) * zn * zn
                        } else {
                            0.0
                        };
                    }
                }
            }
            canonical_state.refresh_panel_arrays_from_rows();
            let upper_ue_alpha = canonical_state.operating_sensitivity_view(XfoilSurface::Upper);
            let lower_ue_alpha = canonical_state.operating_sensitivity_view(XfoilSurface::Lower);
            let (upper_ue_operating, lower_ue_operating) =
                operating_sensitivity_for_mode(operating_mode, &upper_ue_alpha, &lower_ue_alpha);
            let (sst_go_iter, sst_gp_iter) =
                find_stagnation_with_derivs(canonical_state.panel_gamma(), &full_arc)
                    .map(|stag| (stag.sst_go, stag.sst_gp))
                    .unwrap_or((sst_go, sst_gp));
            global_system.set_stagnation_derivs(sst_go_iter, sst_gp_iter);

            // Save the current accepted state so a rejected Newton step rolls
            // back to the previous iterate, matching XFOIL's late-iteration
            // behavior near the stagnation region.
            let canonical_iter_backup = canonical_state.clone();

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

            // === RE-MARCH BOUNDARY LAYERS WITH CURRENT Ue (MRCHDU style) ===
            // XFOIL calls MRCHDU at the top of every SETBL (xbl.f line 98).
            // MRCHDU walks the existing BL arrays in-place, using the current
            // theta/dstar/ctau as initial guesses and refining with a single
            // Newton pass at each station.  This preserves the evolved laminar
            // bubble and N-factor history across Newton iterations, unlike the
            // previous march_surface call which re-initialised from Thwaites.

            canonical_state.march_mixed_du(XfoilSurface::Upper, re, msq, &march_config);
            canonical_state.march_mixed_du(XfoilSurface::Lower, re, msq, &march_config);
            let upper_debug = canonical_state.upper_station_view();
            let lower_debug = canonical_state.lower_station_view();

            // Debug after blvar - emit BL state summaries
            if rustfoil_bl::is_debug_active() {
                // Collect upper surface station states (every 10th station)
                let upper_station_states: Vec<rustfoil_bl::StationState> = upper_debug
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i % 10 == 0 || *i == upper_debug.len() - 1)
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
                let lower_station_states: Vec<rustfoil_bl::StationState> = lower_debug
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i % 10 == 0 || *i == lower_debug.len() - 1)
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
                    for i in 1..4.min(upper_debug.len()) {
                        let s = &upper_debug[i];
                        eprintln!("[DEBUG Global Newton]   upper[{}] theta={:.6e} dstar={:.6e} Hk={:.4}", 
                            i, s.theta, s.delta_star, s.hk);
                    }
                }
            }

            // Build the global Newton system with full cross-surface coupling
            let mut newton_view = build_canonical_newton_view(
                &canonical_state,
                upper_ue_operating.clone(),
                lower_ue_operating.clone(),
                global_system.ante,
                sst_go_iter,
                sst_gp_iter,
            );
            global_system.build_global_system_from_view(
                &newton_view,
                dij,
                config.ncrit,
                msq,
                re,
                iter,
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
                    None,
                    None,
                    None,
                    None,
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
                    upper_debug[0].x, upper_debug[0].theta, upper_debug[0].delta_star, upper_debug[0].u, upper_debug[0].mass_defect);
                eprintln!("  Upper[1]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    upper_debug[1].x, upper_debug[1].theta, upper_debug[1].delta_star, upper_debug[1].u, upper_debug[1].mass_defect);
                eprintln!("  Upper[2]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    upper_debug[2].x, upper_debug[2].theta, upper_debug[2].delta_star, upper_debug[2].u, upper_debug[2].mass_defect);
                    
                // Lower surface first stations
                eprintln!("  Lower[0]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    lower_debug[0].x, lower_debug[0].theta, lower_debug[0].delta_star, lower_debug[0].u, lower_debug[0].mass_defect);
                eprintln!("  Lower[1]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    lower_debug[1].x, lower_debug[1].theta, lower_debug[1].delta_star, lower_debug[1].u, lower_debug[1].mass_defect);
                eprintln!("  Lower[2]: x={:.6e}, theta={:.6e}, dstar={:.6e}, u={:.6}, mass={:.6e}", 
                    lower_debug[2].x, lower_debug[2].theta, lower_debug[2].delta_star, lower_debug[2].u, lower_debug[2].mass_defect);
                
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
                    if ibl < upper_debug.len() {
                        let s = &upper_debug[ibl];
                        eprintln!("[DEBUG Newton] iter {} upper[{}]: theta={:.6e}, dstar={:.6e}, Ue={:.6}, ctau={:.6}, mass={:.6e}",
                            iter, ibl, s.theta, s.delta_star, s.u, s.ctau, s.mass_defect);
                    }
                    if ibl < lower_debug.len() {
                        let s = &lower_debug[ibl];
                        eprintln!("[DEBUG Newton] iter {} lower[{}]: theta={:.6e}, dstar={:.6e}, Ue={:.6}, ctau={:.6}, mass={:.6e}",
                            iter, ibl, s.theta, s.delta_star, s.u, s.ctau, s.mass_defect);
                    }
                }
            }
            
            if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() && iter < 3 {
                eprintln!("[DEBUG Newton] iter {} residual before solve: {:.6e}", iter, residual);
            }

            // Solve the global Newton system
            let solve_result = solve_global_system(&mut global_system);
            let preview = preview_global_update_ue_from_view(
                &newton_view,
                &solve_result.state_deltas,
                Some(&solve_result.operating_deltas),
                &global_system,
                dij,
            );
            let (cl_new, cl_alpha, cl_ac) = compute_trial_cl_terms_from_view(
                panel_x,
                panel_y,
                &newton_view,
                &preview.upper_u_new,
                &preview.lower_u_new,
                &preview.upper_u_ac,
                &preview.lower_u_ac,
                alpha_rad,
            );
            let dac = compute_operating_correction(operating_mode, cl_new, cl_alpha, cl_ac);
            let deltas: Vec<[f64; 3]> = solve_result
                .state_deltas
                .iter()
                .zip(solve_result.operating_deltas.iter())
                .map(|(base, operating)| {
                    [
                        base[0] - dac * operating[0],
                        base[1] - dac * operating[1],
                        base[2] - dac * operating[2],
                    ]
                })
                .collect();
            
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
                canonical_state = canonical_iter_backup;
                canonical_state.write_back_station_views(upper_stations, lower_stations);
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
            
            let update_result = apply_global_updates_from_view(
                &mut newton_view,
                &solve_result.state_deltas,
                Some(&solve_result.operating_deltas),
                &global_system,
                iter + 1,
                rlx,
                dij,
                dac,
                msq,
                re,
            );
            canonical_state.set_operating_scratch(
                &update_result.upper_u_new,
                &update_result.lower_u_new,
                &update_result.upper_u_ac,
                &update_result.lower_u_ac,
                dac,
                update_result.relaxation_used,
            );
            
            // Use RMS of normalized changes (XFOIL's RMSBL) for convergence check
            // This is the correct metric - NOT the raw residuals from VDEL
            residual = update_result.rms_change;
            
            // Debug: Check relaxation used and residual after update for iterations 6-7
            if iter >= 5 && iter <= 8 {
                eprintln!("[DEBUG Newton] iter {} relaxation used: {:.6e} (requested: {:.6e}), rms_change: {:.6e}", 
                    iter, update_result.relaxation_used, rlx_before, residual);
            }

            // Match XFOIL's DBGFULL_BL_STATE timing: emit the updated station
            // state before QVIS/GAM refresh and before STMOVE relocates the split.
            if rustfoil_bl::is_debug_active() {
                let mut debug_state = canonical_state.clone();
                debug_state.overwrite_from_newton_view(&newton_view);
                let upper_post_update = debug_state.upper_station_view();
                let lower_post_update = debug_state.lower_station_view();
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_bl_state(
                    iteration,
                    surface_state(&upper_post_update),
                    surface_state(&lower_post_update),
                ));
            }

            
            // XFOIL updates QVIS/GAM from the current UEDG, then relocates the
            // stagnation split with STMOVE before the next Newton iteration.
            canonical_state.overwrite_from_newton_view(&newton_view);
            let full_gamma_iter = canonical_state.panel_gamma().to_vec();
            let old_ist = canonical_state.ist;
            if let Some(TransitionalStmoveResult { ist: new_ist }) = canonical_state.apply_stmove_like_xfoil(
                ue_inviscid_full,
                panel_x,
                panel_y,
                &full_arc,
                &full_gamma_iter,
                re,
                old_ist,
            ) {
                if new_ist != old_ist {
                    if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                        eprintln!(
                            "[DEBUG Newton] STMOVE shifted stagnation panel: {} -> {}",
                            old_ist, new_ist
                        );
                    }
                    canonical_state.adjust_transition_metadata_for_stmove(old_ist, new_ist);
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
                canonical_state = canonical_iter_backup;
                canonical_state.write_back_station_views(upper_stations, lower_stations);
                converged = false;
                iteration = iter + 1;
                break;
            }

            
            // Stop if residual is increasing rapidly (Newton diverging badly)
            // Allow slow increase since our residual may oscillate
            if iter > 5 && residual > prev_residual * 5.0 {
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG Newton] Rapid increase at iter {}: residual={:.6e}, prev={:.6e} (reverting)", iter, residual, prev_residual);
                }
                canonical_state = canonical_iter_backup;
                canonical_state.write_back_station_views(upper_stations, lower_stations);
                converged = false;
                iteration = iter + 1;
                break;
            }
            if rustfoil_bl::is_debug_active() {
                let upper_iter_end = canonical_state.upper_station_view();
                let lower_iter_end = canonical_state.lower_station_view();
                // Emit FULL_NFACTOR debug event with N-factors at all stations
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_nfactor(
                    iteration,
                    upper_iter_end.iter().map(|s| s.ampl).collect(),
                    lower_iter_end.iter().map(|s| s.ampl).collect(),
                    upper_iter_end.iter().map(|s| s.x).collect(),
                    lower_iter_end.iter().map(|s| s.x).collect(),
                ));
            }

            prev_residual = residual;
        }

    }

    // Update transition locations from the current canonical owner state.
    let final_stagnation = find_stagnation_with_derivs(ue_inviscid_full, &full_arc);
    if let Some(stag) = final_stagnation {
        canonical_state.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
    }

    // Convert arc-length transition locations on the extracted surface geometry,
    // matching XFOIL's surface-based xt -> x/c ownership more closely than BL-row
    // interpolation.
    if let Some(stag) = final_stagnation {
        canonical_state.xtran_upper = transition_arc_to_chord_on_surface(
            stag,
            &full_arc,
            panel_x,
            panel_y,
            ue_inviscid_full,
            true,
            canonical_state.xtran_upper,
        );
        canonical_state.xtran_lower = transition_arc_to_chord_on_surface(
            stag,
            &full_arc,
            panel_x,
            panel_y,
            ue_inviscid_full,
            false,
            canonical_state.xtran_lower,
        );
    }

    // Check for separation on either surface
    let x_separation = upper_result
        .x_separation
        .or(lower_result.x_separation);

    // === QVFUE + GAMQV: Update QVIS from edge velocities, then gamma from QVIS ===
    // This is part of XFOIL's iteration loop that couples the BL solution back to
    // the panel method circulation. The canonical state already owns the accepted
    // BL rows here, so just refresh its panel arrays in place.
    refresh_canonical_panel_arrays(&mut canonical_state);
    canonical_state.write_back_station_views(upper_stations, lower_stations);

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
    
    let full_gamma = canonical_state.panel_gamma().to_vec();
    
    // Emit FULL_GAMMA_ITER debug event with circulation at all panels
    if rustfoil_bl::is_debug_active() {
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_gamma_iter(
            iteration,
            full_gamma.clone(),
        ));
    }
    
    // Use XFOIL's exact CLCALC formula with the canonical GAM array.
    let forces = compute_forces_from_canonical_state(
        &canonical_state,
        panel_x,
        panel_y,
        alpha_rad,
        config,
    );
    
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
        x_tr_upper: canonical_state.xtran_upper.unwrap_or(1.0),
        x_tr_lower: canonical_state.xtran_lower.unwrap_or(1.0),
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
    use crate::viscous::setup::{
        compute_arc_lengths, extract_surface_xfoil, find_stagnation, initialize_bl_stations,
        initialize_surface_stations_with_panel_idx,
    };

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

    #[test]
    fn test_apply_stmove_preserves_lower_wake_on_changed_ist() {
        let panel_x = vec![1.0, 0.6, 0.1, 0.0, 0.3, 0.7, 1.0];
        let panel_y = vec![0.0, 0.08, 0.06, 0.0, -0.05, -0.04, 0.0];
        let full_arc = compute_arc_lengths(&panel_x, &panel_y);
        let old_gamma = vec![0.8, 0.5, 0.2, -0.1, -0.4, -0.7, -0.9];
        let old_ist = 2;
        let old_sst = 2.5;
        let old_ue_stag = 0.0;

        let (upper_arc, upper_x, _, upper_ue) = extract_surface_xfoil(
            old_ist,
            old_sst,
            old_ue_stag,
            &full_arc,
            &panel_x,
            &panel_y,
            &old_gamma,
            true,
        );
        let (lower_arc, lower_x, _, lower_ue) = extract_surface_xfoil(
            old_ist,
            old_sst,
            old_ue_stag,
            &full_arc,
            &panel_x,
            &panel_y,
            &old_gamma,
            false,
        );

        let n_airfoil = panel_x.len();
        let re = 1.0e6;
        let mut upper_stations = initialize_surface_stations_with_panel_idx(
            &upper_arc,
            &upper_ue,
            &upper_x,
            old_ist,
            n_airfoil,
            true,
            re,
        );
        let mut lower_stations = initialize_surface_stations_with_panel_idx(
            &lower_arc,
            &lower_ue,
            &lower_x,
            old_ist,
            n_airfoil,
            false,
            re,
        );
        for wake_idx in 0..2 {
            let mut wake_station = BlStation::new();
            wake_station.x = lower_stations.last().unwrap().x + 0.1 * (wake_idx as f64 + 1.0);
            wake_station.x_coord = 1.0 + 0.1 * wake_idx as f64;
            wake_station.panel_idx = n_airfoil + wake_idx;
            wake_station.u = 0.7 - 0.05 * wake_idx as f64;
            wake_station.mass_defect = wake_station.u * wake_station.delta_star;
            wake_station.is_wake = true;
            lower_stations.push(wake_station);
        }
        let old_lower_wake_count = lower_stations.iter().filter(|station| station.is_wake).count();
        let current_gamma = vec![0.8, 0.2, -0.15, -0.4, -0.6, -0.8, -0.95];

        let mut canonical_state = build_canonical_state(
            &upper_stations,
            &lower_stations,
            panel_x.len(),
            find_stagnation_with_derivs(&old_gamma, &full_arc),
        );
        let stmove = canonical_state
            .apply_stmove_like_xfoil(
                &old_gamma,
                &panel_x,
                &panel_y,
                &full_arc,
                &current_gamma,
                re,
                old_ist,
            )
            .expect("stagnation relocation");
        canonical_state.write_back_station_views(&mut upper_stations, &mut lower_stations);
        let new_ist = stmove.ist;

        assert_eq!(new_ist, 1);
        assert_eq!(
            lower_stations.iter().filter(|station| station.is_wake).count(),
            old_lower_wake_count
        );
        assert!(lower_stations.last().unwrap().is_wake);
    }

    #[test]
    fn test_build_panel_gamma_respects_write_order() {
        let mut upper = vec![BlStation::new(), BlStation::new()];
        let mut lower = vec![BlStation::new(), BlStation::new()];

        upper[1].panel_idx = 3;
        upper[1].u = 0.8;
        lower[1].panel_idx = 3;
        lower[1].u = 0.6;

        let preserve = build_panel_gamma_from_stations(
            &upper,
            &lower,
            8,
            GammaWriteOrder::PreserveFirstWriter,
        );
        let lower_overwrites = build_panel_gamma_from_stations(
            &upper,
            &lower,
            8,
            GammaWriteOrder::LowerOverwritesShared,
        );

        assert!((preserve[3] - 0.8).abs() < 1e-12);
        assert!((lower_overwrites[3] + 0.6).abs() < 1e-12);
    }

    #[test]
    fn test_canonical_state_bridge_round_trips_views() {
        let mut upper = vec![BlStation::new(), BlStation::new()];
        let mut lower = vec![BlStation::new(), BlStation::new(), BlStation::new()];

        upper[0].panel_idx = usize::MAX;
        upper[1].panel_idx = 5;
        upper[1].u = 0.8;

        lower[0].panel_idx = usize::MAX;
        lower[1].panel_idx = 6;
        lower[1].u = 0.7;
        lower[2].panel_idx = 7;
        lower[2].u = 0.6;
        lower[2].is_wake = true;
        lower[2].is_laminar = false;
        lower[2].is_turbulent = true;

        let stagnation = Some(StagnationResult {
            ist: 5,
            sst: 0.12,
            sst_go: -0.03,
            sst_gp: 0.02,
        });
        let state = build_canonical_state(&upper, &lower, 12, stagnation);
        let upper_view = state.upper_station_view();
        let lower_view = state.lower_station_view();

        assert_eq!(state.ist, 5);
        assert_eq!(state.row(2, 1).unwrap().panel_idx, 5);
        assert_eq!(upper_view[1].panel_idx, 5);
        assert_eq!(lower_view[2].panel_idx, 7);
        assert!(lower_view[2].is_wake);
    }
}
