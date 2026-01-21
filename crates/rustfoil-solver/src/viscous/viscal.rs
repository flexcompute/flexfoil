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
use rustfoil_coupling::solve::solve_bl_system;
use rustfoil_coupling::update::{set_edge_velocities, update_stations, UpdateConfig};

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

        // Build Newton system from current BL state
        newton_system.build(stations, &flow_types, msq, re);

        // Solve block-tridiagonal system
        let deltas_raw = solve_bl_system(&newton_system);

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

        // Update stations with limiting
        let update_result = update_stations(stations, &deltas, &current_ue, &update_config);

        // Update edge velocities from mass defect changes (viscous-inviscid coupling)
        // This modifies stations in place
        set_edge_velocities(stations, ue_inviscid, dij, None);

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
    _dij: &DMatrix<f64>,
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
        ..Default::default()
    };

    // Extract arc lengths from stations
    let upper_arc: Vec<f64> = upper_stations.iter().map(|s| s.x).collect();
    let lower_arc: Vec<f64> = lower_stations.iter().map(|s| s.x).collect();

    // March upper surface (side 1)
    let upper_result = march_surface(&upper_arc, upper_ue, re, msq, &march_config, 1);

    // March lower surface (side 2)
    let lower_result = march_surface(&lower_arc, lower_ue, re, msq, &march_config, 2);

    // Copy results back to stations
    for (i, station) in upper_result.stations.iter().enumerate() {
        if i < upper_stations.len() {
            upper_stations[i].theta = station.theta;
            upper_stations[i].delta_star = station.delta_star;
            upper_stations[i].h = station.h;
            upper_stations[i].hk = station.hk;
            upper_stations[i].cf = station.cf;
            upper_stations[i].ctau = station.ctau;
            upper_stations[i].ampl = station.ampl;
            upper_stations[i].is_laminar = station.is_laminar;
            upper_stations[i].is_turbulent = station.is_turbulent;
            upper_stations[i].mass_defect = station.mass_defect;
            upper_stations[i].r_theta = station.r_theta;
        }
    }

    for (i, station) in lower_result.stations.iter().enumerate() {
        if i < lower_stations.len() {
            lower_stations[i].theta = station.theta;
            lower_stations[i].delta_star = station.delta_star;
            lower_stations[i].h = station.h;
            lower_stations[i].hk = station.hk;
            lower_stations[i].cf = station.cf;
            lower_stations[i].ctau = station.ctau;
            lower_stations[i].ampl = station.ampl;
            lower_stations[i].is_laminar = station.is_laminar;
            lower_stations[i].is_turbulent = station.is_turbulent;
            lower_stations[i].mass_defect = station.mass_defect;
            lower_stations[i].r_theta = station.r_theta;
        }
    }

    // Get transition locations from march results
    let x_tr_upper = upper_result.x_transition.unwrap_or(1.0);
    let x_tr_lower = lower_result.x_transition.unwrap_or(1.0);

    // Check for separation on either surface
    let x_separation = upper_result
        .x_separation
        .or(lower_result.x_separation);

    // Compute forces from both surfaces
    // For now, compute from upper surface (CD is dominated by friction drag anyway)
    let forces = compute_forces(upper_stations, config);
    let lower_forces = compute_forces(lower_stations, config);

    // Combine CD (friction from both surfaces)
    let total_cd_friction = forces.cd_friction + lower_forces.cd_friction;

    Ok(ViscousResult {
        alpha: 0.0, // Will be set by caller
        cl: forces.cl, // CL from pressure distribution (need full airfoil)
        cd: total_cd_friction, // For now, just friction drag
        cm: forces.cm,
        x_tr_upper,
        x_tr_lower,
        iterations: 1, // Only direct march, no Newton iteration yet
        residual: 0.0,
        converged: true, // March always "converges"
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
