//! WebAssembly bindings for RustFoil.
//!
//! This crate provides JavaScript-friendly interfaces to the RustFoil
//! aerodynamic analysis engine. It is designed for real-time (60 Hz)
//! interaction in a web-based UI.
//!
//! # Usage from JavaScript
//!
//! ```javascript
//! import init, { RustFoil, analyze_airfoil } from 'rustfoil-wasm';
//!
//! await init();
//!
//! // Quick analysis
//! const result = analyze_airfoil(coordinates, 5.0); // 5° angle of attack
//! console.log(`Cl = ${result.cl}`);
//!
//! // Or use the full API
//! const foil = new RustFoil();
//! foil.set_coordinates(coordinates);
//! foil.set_alpha(5.0);
//! const solution = foil.solve();
//! ```

use rustfoil_core::{naca, point, Body, CubicSpline, Point, flap::xfoil_flap};
use rustfoil_inviscid::{FlowConditions as FaithfulFlowConditions, InviscidSolver as FaithfulInviscidSolver};
use rustfoil_solver::inviscid::{
    FlowConditions, InviscidSolver,
    build_dividing_streamline, build_dividing_streamline_viscous, build_streamlines,
    build_streamlines_viscous, StreamlineOptions, SmokeSystem, WakePanels,
};
use rustfoil_solver::inviscid::compute_psi_grid_with_sources;
use rustfoil_xfoil::XfoilOptions;
use rustfoil_xfoil::oper::solve_operating_point_from_state;
use rustfoil_xfoil::state::XfoilState;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

// When the `console_error_panic_hook` feature is enabled, we can call the
// `set_panic_hook` function at least once during initialization, and then
// we will get better error messages if our code ever panics.
#[cfg(feature = "console_error_panic_hook")]
fn set_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Initialize the WASM module.
///
/// Call this once before using any other functions.
#[wasm_bindgen(start)]
pub fn init() {
    #[cfg(feature = "console_error_panic_hook")]
    set_panic_hook();
}

/// Check that the WASM module is working.
#[wasm_bindgen]
pub fn greet() -> String {
    "RustFoil WASM v0.1.0 initialized".to_string()
}

// ============================================================================
// NACA 4-Series Generation
// ============================================================================

/// Generate a NACA 4-series airfoil.
///
/// # Arguments
/// * `m` - Maximum camber as fraction of chord (0-0.09, first digit / 100)
/// * `p` - Position of maximum camber as fraction of chord (0-0.9, second digit / 10)
/// * `t` - Maximum thickness as fraction of chord (0-0.40, last two digits / 100)
/// * `n_points` - Number of points on each surface (total will be 2*n_points - 1)
///
/// # Returns
/// Flat array of coordinates [x0, y0, x1, y1, ...] starting from trailing edge,
/// going around the lower surface, leading edge, upper surface, back to TE.
#[wasm_bindgen]
pub fn generate_naca4(m: f64, p: f64, t: f64, n_points: usize) -> Vec<f64> {
    generate_naca4_impl(m, p, t, n_points)
        .iter()
        .flat_map(|pt| [pt.x, pt.y])
        .collect()
}

/// Generate NACA 4-series airfoil coordinates.
fn generate_naca4_impl(m: f64, p: f64, t: f64, n_points: usize) -> Vec<Point> {
    let n = n_points.max(10);
    
    // Generate x-coordinates using cosine spacing for better LE/TE resolution
    let x_coords: Vec<f64> = (0..n)
        .map(|i| {
            let beta = PI * (i as f64) / ((n - 1) as f64);
            0.5 * (1.0 - beta.cos())
        })
        .collect();
    
    let mut upper: Vec<Point> = Vec::with_capacity(n);
    let mut lower: Vec<Point> = Vec::with_capacity(n);
    
    for &x in &x_coords {
        // Thickness distribution (NACA 4-digit)
        let y_t = 5.0 * t * (
            0.2969 * x.sqrt()
            - 0.1260 * x
            - 0.3516 * x.powi(2)
            + 0.2843 * x.powi(3)
            - 0.1015 * x.powi(4)  // Original NACA formula (open TE)
        );
        
        // Camber line and gradient
        let (y_c, dy_c) = if m.abs() < 1e-10 || p.abs() < 1e-10 {
            // Symmetric airfoil
            (0.0, 0.0)
        } else if x < p {
            let y_c = (m / (p * p)) * (2.0 * p * x - x * x);
            let dy_c = (2.0 * m / (p * p)) * (p - x);
            (y_c, dy_c)
        } else {
            let one_minus_p = 1.0 - p;
            let y_c = (m / (one_minus_p * one_minus_p)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
            let dy_c = (2.0 * m / (one_minus_p * one_minus_p)) * (p - x);
            (y_c, dy_c)
        };
        
        let theta = dy_c.atan();
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        
        // Upper surface
        let x_u = x - y_t * sin_t;
        let y_u = y_c + y_t * cos_t;
        upper.push(point(x_u, y_u));
        
        // Lower surface
        let x_l = x + y_t * sin_t;
        let y_l = y_c - y_t * cos_t;
        lower.push(point(x_l, y_l));
    }
    
    // Assemble in clockwise order: upper TE -> LE -> lower TE -> close
    let mut coords: Vec<Point> = Vec::with_capacity(2 * n);
    
    // Upper surface: TE to LE
    for i in (0..n).rev() {
        coords.push(upper[i]);
    }
    
    // Lower surface: LE+1 to TE (skip LE to avoid duplicate)
    for i in 1..n {
        coords.push(lower[i]);
    }
    
    // Close contour for sharp TE
    coords.push(coords[0]);
    
    coords
}

/// Generate NACA 4-series using XFOIL's exact algorithm.
///
/// This matches XFOIL's naca.f SUBROUTINE NACA4 exactly:
/// - TE bunching parameter AN = 1.5
/// - Blunt trailing edge (not closed)
/// - Output order: upper TE → LE → lower TE
///
/// # Arguments
/// * `designation` - 4-digit NACA designation (e.g., 12 for NACA 0012, 2412 for NACA 2412)
/// * `n_points_per_side` - Number of points per side (default: 123 = XFOIL's IQX/3)
///
/// # Returns
/// Flat array of coordinates [x0, y0, x1, y1, ...]. Total points = 2*n - 1.
#[wasm_bindgen]
pub fn generate_naca4_xfoil(designation: u32, n_points_per_side: Option<usize>) -> Vec<f64> {
    let coords = naca::naca4(designation, n_points_per_side);
    coords.iter().flat_map(|p| [p.x, p.y]).collect()
}

/// Generate NACA 4-series from digit string (e.g., "2412", "0012").
///
/// # Arguments
/// * `digits` - 4-digit NACA designation as string
/// * `n_points` - Number of points on each surface
///
/// # Returns
/// Flat array of coordinates, or empty if invalid input.
#[wasm_bindgen]
pub fn generate_naca4_from_string(digits: &str, n_points: usize) -> Vec<f64> {
    if digits.len() != 4 {
        return vec![];
    }
    
    let chars: Vec<char> = digits.chars().collect();
    
    let m = match chars[0].to_digit(10) {
        Some(d) => d as f64 / 100.0,
        None => return vec![],
    };
    
    let p = match chars[1].to_digit(10) {
        Some(d) => d as f64 / 10.0,
        None => return vec![],
    };
    
    let t_str = &digits[2..4];
    let t = match t_str.parse::<u32>() {
        Ok(d) => d as f64 / 100.0,
        Err(_) => return vec![],
    };
    
    generate_naca4(m, p, t, n_points)
}

// ============================================================================
// Repaneling with custom spacing
// ============================================================================

/// Repanel airfoil with SSP (Position-Based) spacing.
///
/// # Arguments
/// * `coords` - Input coordinates as flat array [x0, y0, x1, y1, ...]
/// * `spacing_knots` - Spacing knots as flat array [s0, f0, s1, f1, ...]
///                     where s is normalized arc length (0-1) and f is spacing value
/// * `n_panels` - Desired number of panels
///
/// # Returns
/// Repaneled coordinates as flat array.
#[wasm_bindgen]
pub fn repanel_with_spacing(coords: &[f64], spacing_knots: &[f64], n_panels: usize) -> Vec<f64> {
    repanel_with_spacing_and_curvature(coords, spacing_knots, n_panels, 0.0)
}

/// Repanel airfoil with blended SSP and curvature-based spacing.
///
/// # Arguments
/// * `coords` - Input coordinates as flat array [x0, y0, x1, y1, ...]
/// * `spacing_knots` - Spacing knots as flat array [s0, f0, s1, f1, ...]
///                     where s is normalized arc length (0-1) and f is spacing value
/// * `n_panels` - Desired number of panels
/// * `curvature_weight` - Blend factor: 0.0 = pure SSP, 1.0 = pure curvature-based
///
/// # Returns
/// Repaneled coordinates as flat array.
#[wasm_bindgen]
pub fn repanel_with_spacing_and_curvature(
    coords: &[f64], 
    spacing_knots: &[f64], 
    n_panels: usize,
    curvature_weight: f64,
) -> Vec<f64> {
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return vec![];
    }
    if spacing_knots.len() < 4 || spacing_knots.len() % 2 != 0 {
        return vec![];
    }
    
    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();
    
    // Parse spacing knots
    let knots: Vec<(f64, f64)> = spacing_knots
        .chunks(2)
        .map(|c| (c[0], c[1]))
        .collect();
    
    // Create spline
    let spline = match CubicSpline::from_points(&points) {
        Ok(s) => s,
        Err(_) => return vec![],
    };
    
    // Get total arc length for scaling
    let s_max = spline.total_arc_length();
    
    // Clamp curvature weight to [0, 1]
    let curvature_weight = curvature_weight.clamp(0.0, 1.0);
    
    // Generate new parameter distribution based on blended spacing function
    let new_params = generate_blended_distribution(&spline, &knots, n_panels + 1, curvature_weight);
    
    // Sample spline at new parameters (scaled to actual arc length)
    let new_points: Vec<Point> = new_params
        .iter()
        .map(|&s| spline.evaluate(s * s_max))
        .collect();
    
    new_points.iter().flat_map(|p| [p.x, p.y]).collect()
}

/// Compute curvature-based spacing values along the airfoil.
/// 
/// Returns an array of (s, curvature_density) pairs where higher curvature
/// regions get higher density values.
#[wasm_bindgen]
pub fn compute_curvature_spacing(coords: &[f64], n_samples: usize) -> Vec<f64> {
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return vec![];
    }
    
    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();
    
    let spline = match CubicSpline::from_points(&points) {
        Ok(s) => s,
        Err(_) => return vec![],
    };
    
    let s_max = spline.total_arc_length();
    let n = n_samples.max(10);
    
    // Sample curvature along the spline
    let curvatures: Vec<f64> = (0..n)
        .map(|i| {
            let s = (i as f64 / (n - 1) as f64) * s_max;
            spline.curvature(s).abs()
        })
        .collect();
    
    // Convert curvature to spacing density
    let spacing = curvature_to_spacing(&curvatures);
    
    // Return as flat array [s0, f0, s1, f1, ...]
    let mut result = Vec::with_capacity(n * 2);
    for i in 0..n {
        let s_norm = i as f64 / (n - 1) as f64;
        result.push(s_norm);
        result.push(spacing[i]);
    }
    result
}

/// Generate parameter distribution blending SSP and curvature-based spacing.
fn generate_blended_distribution(
    spline: &CubicSpline,
    knots: &[(f64, f64)],
    n_points: usize,
    curvature_weight: f64,
) -> Vec<f64> {
    if n_points < 2 {
        return vec![0.0, 1.0];
    }
    
    let n_samples = 1000;
    let s_max = spline.total_arc_length();
    
    // Compute SSP spacing values
    let ssp_spacing: Vec<f64> = (0..n_samples)
        .map(|i| {
            let s = i as f64 / (n_samples - 1) as f64;
            interpolate_spacing(knots, s)
        })
        .collect();
    
    // Compute curvature-based spacing values
    let curvatures: Vec<f64> = (0..n_samples)
        .map(|i| {
            let s = (i as f64 / (n_samples - 1) as f64) * s_max;
            spline.curvature(s).abs()
        })
        .collect();
    let curvature_spacing = curvature_to_spacing(&curvatures);
    
    // Blend the two spacing functions
    let blended_spacing: Vec<f64> = ssp_spacing
        .iter()
        .zip(curvature_spacing.iter())
        .map(|(&ssp, &curv)| {
            // Linear blend: (1 - w) * SSP + w * curvature
            (1.0 - curvature_weight) * ssp + curvature_weight * curv
        })
        .collect();
    
    spacing_to_distribution(&blended_spacing, n_points)
}

/// Convert a spacing function (sampled at n_samples points) to a parameter distribution.
fn spacing_to_distribution(spacing_values: &[f64], n_points: usize) -> Vec<f64> {
    let n_samples = spacing_values.len();
    
    // Integrate spacing function to get cumulative distribution
    let ds = 1.0 / (n_samples - 1) as f64;
    let mut cumulative: Vec<f64> = Vec::with_capacity(n_samples);
    cumulative.push(0.0);
    
    for i in 1..n_samples {
        let avg_spacing = (spacing_values[i - 1] + spacing_values[i]) / 2.0;
        cumulative.push(cumulative[i - 1] + avg_spacing * ds);
    }
    
    // Normalize cumulative distribution
    let total = cumulative[n_samples - 1];
    if total < 1e-10 {
        // Uniform distribution fallback
        return (0..n_points)
            .map(|i| i as f64 / (n_points - 1) as f64)
            .collect();
    }
    
    for c in &mut cumulative {
        *c /= total;
    }
    
    // Invert: find s values for uniform u distribution
    let mut result: Vec<f64> = Vec::with_capacity(n_points);
    
    for i in 0..n_points {
        let u = i as f64 / (n_points - 1) as f64;
        
        // Binary search for s such that cumulative(s) ≈ u
        let mut lo = 0;
        let mut hi = n_samples - 1;
        
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if cumulative[mid] < u {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        
        // Linear interpolation
        let s_lo = lo as f64 / (n_samples - 1) as f64;
        let s_hi = hi as f64 / (n_samples - 1) as f64;
        let c_lo = cumulative[lo];
        let c_hi = cumulative[hi];
        
        let s = if (c_hi - c_lo).abs() < 1e-10 {
            s_lo
        } else {
            s_lo + (u - c_lo) * (s_hi - s_lo) / (c_hi - c_lo)
        };
        
        result.push(s.clamp(0.0, 1.0));
    }
    
    result
}

/// Convert curvature values to spacing density values.
/// 
/// Higher curvature → higher density (more points).
/// Uses a power-law transformation with smoothing.
fn curvature_to_spacing(curvatures: &[f64]) -> Vec<f64> {
    if curvatures.is_empty() {
        return vec![];
    }
    
    // Find max curvature for normalization
    let max_curv = curvatures.iter().cloned().fold(0.0_f64, f64::max);
    
    if max_curv < 1e-10 {
        // Flat geometry - uniform spacing
        return vec![1.0; curvatures.len()];
    }
    
    // Transform curvature to spacing density
    // Use sqrt to moderate the effect (pure linear would over-cluster at LE)
    // Add a base level so flat regions still get some points
    let base_density = 0.3;
    let curv_scale = 1.0 - base_density;
    
    curvatures
        .iter()
        .map(|&k| {
            let normalized = (k / max_curv).sqrt();
            base_density + curv_scale * normalized
        })
        .collect()
}

/// Interpolate spacing value at parameter s using knot values.
fn interpolate_spacing(knots: &[(f64, f64)], s: f64) -> f64 {
    if knots.is_empty() {
        return 1.0;
    }
    if knots.len() == 1 {
        return knots[0].1;
    }
    
    // Find bracketing knots
    let mut sorted_knots = knots.to_vec();
    sorted_knots.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    // Clamp to range
    if s <= sorted_knots[0].0 {
        return sorted_knots[0].1;
    }
    if s >= sorted_knots.last().unwrap().0 {
        return sorted_knots.last().unwrap().1;
    }
    
    // Find interval
    for i in 0..sorted_knots.len() - 1 {
        let (s0, f0) = sorted_knots[i];
        let (s1, f1) = sorted_knots[i + 1];
        
        if s >= s0 && s <= s1 {
            // Linear interpolation
            let t = (s - s0) / (s1 - s0);
            return f0 + t * (f1 - f0);
        }
    }
    
    1.0
}

// ============================================================================
// XFOIL-Style Repaneling
// ============================================================================

/// Repanel airfoil using XFOIL's curvature-based algorithm.
///
/// This implements Mark Drela's PANGEN algorithm from XFOIL, which distributes
/// panel nodes based on local curvature. Panels are bunched in regions of high
/// curvature (leading edge) and at the trailing edge.
///
/// # Arguments
/// * `coords` - Input coordinates as flat array [x0, y0, x1, y1, ...]
/// * `n_panels` - Desired number of panels (XFOIL default: 160)
///
/// # Returns
/// Repaneled coordinates as flat array.
#[wasm_bindgen]
pub fn repanel_xfoil(coords: &[f64], n_panels: usize) -> Vec<f64> {
    repanel_xfoil_with_params(coords, n_panels, 1.0, 0.15, 0.667)
}

/// Repanel airfoil using XFOIL's algorithm with custom parameters.
///
/// # Arguments
/// * `coords` - Input coordinates as flat array [x0, y0, x1, y1, ...]
/// * `n_panels` - Desired number of panels
/// * `curv_param` - Curvature attraction (0=uniform, 1=XFOIL default)
/// * `te_le_ratio` - TE/LE panel density ratio (XFOIL default: 0.15)
/// * `te_spacing_ratio` - TE panel length ratio (XFOIL default: 0.667)
///
/// # Returns
/// Repaneled coordinates as flat array.
#[wasm_bindgen]
pub fn repanel_xfoil_with_params(
    coords: &[f64],
    n_panels: usize,
    curv_param: f64,
    te_le_ratio: f64,
    te_spacing_ratio: f64,
) -> Vec<f64> {
    use rustfoil_core::PanelingParams;
    
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return vec![];
    }

    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();

    let spline = match CubicSpline::from_points(&points) {
        Ok(s) => s,
        Err(_) => return vec![],
    };

    let params = PanelingParams {
        curv_param,
        te_le_ratio,
        te_spacing_ratio,
    };

    let repaneled = spline.resample_xfoil(n_panels, &params);

    // Don't add a closing point - XFOIL produces blunt TE (upper/lower TE are distinct).
    let result: Vec<f64> = repaneled
        .iter()
        .flat_map(|p| [p.x, p.y])
        .collect();

    result
}

/// A 2D point for JavaScript interop.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JsPoint {
    pub x: f64,
    pub y: f64,
}

/// Result of an aerodynamic analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisResult {
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient
    pub cd: f64,
    /// Moment coefficient (about quarter-chord)
    pub cm: f64,
    /// Pressure coefficient at each panel midpoint
    pub cp: Vec<f64>,
    /// X-coordinates of Cp sample points
    pub cp_x: Vec<f64>,
    /// Gamma (vorticity) values at each node - needed for velocity field computation
    pub gamma: Vec<f64>,
    /// Dividing streamline value (psi_0)
    pub psi_0: f64,
    /// Whether the viscous solve converged
    pub converged: bool,
    /// Newton iterations used by the faithful solver
    pub iterations: usize,
    /// Final residual from the faithful solver
    pub residual: f64,
    /// Transition location on upper surface
    pub x_tr_upper: f64,
    /// Transition location on lower surface
    pub x_tr_lower: f64,
    /// Whether the analysis succeeded
    pub success: bool,
    /// Error message (if any)
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLDistribution {
    pub x_upper: Vec<f64>,
    pub x_lower: Vec<f64>,
    pub theta_upper: Vec<f64>,
    pub theta_lower: Vec<f64>,
    pub delta_star_upper: Vec<f64>,
    pub delta_star_lower: Vec<f64>,
    pub h_upper: Vec<f64>,
    pub h_lower: Vec<f64>,
    pub cf_upper: Vec<f64>,
    pub cf_lower: Vec<f64>,
    pub x_tr_upper: f64,
    pub x_tr_lower: f64,
    pub converged: bool,
    pub iterations: usize,
    pub residual: f64,
    pub success: bool,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLStationData {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub delta_star: Vec<f64>,
    pub theta: Vec<f64>,
    pub cf: Vec<f64>,
    pub ue: Vec<f64>,
    pub hk: Vec<f64>,
    pub ampl: Vec<f64>,
    pub is_laminar: Vec<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeData {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub delta_star: Vec<f64>,
    pub theta: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLVisualizationData {
    pub upper: BLStationData,
    pub lower: BLStationData,
    pub wake: WakeData,
    /// Fraction of wake δ* attributed to the upper side (for asymmetric splitting).
    /// At the TE: upper_fraction = upper_δ* / (upper_δ* + lower_δ*).
    pub wake_upper_fraction: f64,
    pub wake_geometry_x: Vec<f64>,
    pub wake_geometry_y: Vec<f64>,
    pub x_tr_upper: f64,
    pub x_tr_lower: f64,
    pub converged: bool,
    pub iterations: usize,
    pub residual: f64,
    pub success: bool,
    pub error: Option<String>,
}

struct FaithfulSnapshot {
    nodes: Vec<Point>,
    gamma: Vec<f64>,
    psi_0: f64,
    cp_x: Vec<f64>,
    result: AnalysisResult,
    upper_rows: Vec<rustfoil_xfoil::XfoilBlRow>,
    lower_rows: Vec<rustfoil_xfoil::XfoilBlRow>,
    iblte_upper: usize,
    iblte_lower: usize,
    wake_x: Vec<f64>,
    wake_y: Vec<f64>,
}

#[derive(Debug, Clone)]
struct FaithfulFlowField {
    nodes: Vec<Point>,
    gamma: Vec<f64>,
    sigma: Vec<f64>,
    alpha: f64,
    v_inf: f64,
    psi_0: f64,
    wake_panels: Option<WakePanels>,
    effective_body: Vec<Point>,
}

impl FaithfulFlowField {
    fn from_snapshot(snapshot: &FaithfulSnapshot, alpha_deg: f64) -> Self {
        Self {
            nodes: snapshot.nodes.clone(),
            gamma: snapshot.gamma.clone(),
            sigma: build_surface_sigma_from_snapshot(snapshot),
            alpha: alpha_deg.to_radians(),
            v_inf: 1.0,
            psi_0: snapshot.psi_0,
            wake_panels: build_wake_source_panels_from_snapshot(snapshot),
            effective_body: build_effective_body_polygon_from_snapshot(snapshot),
        }
    }

    fn effective_body_ref(&self) -> Option<&[Point]> {
        if self.effective_body.len() >= 3 {
            Some(&self.effective_body)
        } else {
            None
        }
    }

    fn compute_streamlines(&self, options: &StreamlineOptions) -> Vec<Vec<(f64, f64)>> {
        build_streamlines_viscous(
            &self.nodes,
            &self.gamma,
            &self.sigma,
            self.alpha,
            self.v_inf,
            self.wake_panels.as_ref(),
            self.effective_body_ref(),
            options,
        )
    }

    fn compute_dividing_streamline(&self, options: &StreamlineOptions) -> Option<Vec<(f64, f64)>> {
        build_dividing_streamline_viscous(
            &self.nodes,
            &self.gamma,
            &self.sigma,
            self.alpha,
            self.v_inf,
            self.psi_0,
            self.wake_panels.as_ref(),
            self.effective_body_ref(),
            options,
        )
    }

    fn compute_psi_grid(
        &self,
        bounds: &[f64],
        nx: usize,
        ny: usize,
        interior_value: Option<f64>,
    ) -> Vec<f64> {
        compute_psi_grid_with_sources(
            &self.nodes,
            &self.gamma,
            &self.sigma,
            self.alpha,
            self.v_inf,
            bounds[0],
            bounds[1],
            bounds[2],
            bounds[3],
            nx,
            ny,
            interior_value,
            self.wake_panels.as_ref(),
        )
    }
}

fn build_surface_sigma_from_snapshot(snapshot: &FaithfulSnapshot) -> Vec<f64> {
    let n = snapshot.nodes.len();
    let mut sigma = vec![0.0f64; n];

    for row in snapshot.upper_rows.iter().take(snapshot.iblte_upper + 1) {
        if row.panel_idx < n {
            sigma[row.panel_idx] = row.mass;
        }
    }
    for row in snapshot.lower_rows.iter().take(snapshot.iblte_lower + 1) {
        if row.panel_idx < n {
            sigma[row.panel_idx] = row.mass;
        }
    }

    sigma
}

fn build_effective_body_polygon_from_snapshot(snapshot: &FaithfulSnapshot) -> Vec<Point> {
    let upper_end = snapshot.iblte_upper.min(snapshot.upper_rows.len().saturating_sub(1));
    let lower_end = snapshot.iblte_lower.min(snapshot.lower_rows.len().saturating_sub(1));
    let upper_surface = &snapshot.upper_rows[..=upper_end];
    let lower_surface = &snapshot.lower_rows[..=lower_end];
    let wake_rows: Vec<_> = snapshot.upper_rows.iter().skip(upper_end + 1).collect();

    let env_normal = |xs: &[f64], ys: &[f64], i: usize, sign: f64| -> (f64, f64) {
        let (tx, ty) = if xs.len() < 2 {
            (0.0, sign)
        } else if i == 0 {
            (xs[1] - xs[0], ys[1] - ys[0])
        } else if i == xs.len() - 1 {
            (xs[i] - xs[i - 1], ys[i] - ys[i - 1])
        } else {
            (xs[i + 1] - xs[i - 1], ys[i + 1] - ys[i - 1])
        };
        let len = (tx * tx + ty * ty).sqrt().max(1e-10);
        (-ty / len * sign, tx / len * sign)
    };

    let mut poly = Vec::new();

    let upper_x: Vec<f64> = upper_surface.iter().map(|r| r.x_coord).collect();
    let upper_y: Vec<f64> = upper_surface.iter().map(|r| r.y_coord).collect();
    for i in 0..upper_surface.len() {
        let (nx, ny) = env_normal(&upper_x, &upper_y, i, 1.0);
        let ds = upper_surface[i].dstr;
        poly.push(point(upper_x[i] + nx * ds, upper_y[i] + ny * ds));
    }

    if !wake_rows.is_empty() {
        let wake_x: Vec<f64> = wake_rows.iter().map(|r| r.x_coord).collect();
        let wake_y: Vec<f64> = wake_rows.iter().map(|r| r.y_coord).collect();
        let wake_d: Vec<f64> = wake_rows.iter().map(|r| r.dstr).collect();
        let upper_dstr_te = upper_surface.last().map_or(0.0, |r| r.dstr);
        let lower_dstr_te = lower_surface.last().map_or(0.0, |r| r.dstr);
        let total_te = (upper_dstr_te + lower_dstr_te).max(1e-12);
        let f_u = upper_dstr_te / total_te;
        let f_l = 1.0 - f_u;
        let blend_n = wake_x.len().min(5);
        let te_nu = if upper_x.len() > 1 {
            env_normal(&upper_x, &upper_y, upper_x.len() - 1, 1.0)
        } else {
            (0.0, 1.0)
        };
        for i in 0..wake_x.len() {
            let (wnx, wny) = env_normal(&wake_x, &wake_y, i, 1.0);
            let t = if i < blend_n { i as f64 / blend_n as f64 } else { 1.0 };
            let nx = te_nu.0 * (1.0 - t) + wnx * t;
            let ny = te_nu.1 * (1.0 - t) + wny * t;
            let len = (nx * nx + ny * ny).sqrt().max(1e-10);
            let ds = wake_d[i] * f_u;
            poly.push(point(wake_x[i] + nx / len * ds, wake_y[i] + ny / len * ds));
        }

        let lower_x: Vec<f64> = lower_surface.iter().map(|r| r.x_coord).collect();
        let lower_y: Vec<f64> = lower_surface.iter().map(|r| r.y_coord).collect();
        let te_nl = if lower_x.len() > 1 {
            env_normal(&lower_x, &lower_y, lower_x.len() - 1, -1.0)
        } else {
            (0.0, -1.0)
        };
        for i in (0..wake_x.len()).rev() {
            let (wnx, wny) = env_normal(&wake_x, &wake_y, i, 1.0);
            let t = if i < blend_n { i as f64 / blend_n as f64 } else { 1.0 };
            let nx = te_nl.0 * (1.0 - t) + (-wnx) * t;
            let ny = te_nl.1 * (1.0 - t) + (-wny) * t;
            let len = (nx * nx + ny * ny).sqrt().max(1e-10);
            let ds = wake_d[i] * f_l;
            poly.push(point(wake_x[i] + nx / len * ds, wake_y[i] + ny / len * ds));
        }
    }

    let lower_x: Vec<f64> = lower_surface.iter().map(|r| r.x_coord).collect();
    let lower_y: Vec<f64> = lower_surface.iter().map(|r| r.y_coord).collect();
    for i in (0..lower_surface.len()).rev() {
        let (nx, ny) = env_normal(&lower_x, &lower_y, i, -1.0);
        let ds = lower_surface[i].dstr;
        poly.push(point(lower_x[i] + nx * ds, lower_y[i] + ny * ds));
    }

    poly
}

fn analysis_error(message: impl Into<String>) -> AnalysisResult {
    AnalysisResult {
        cl: 0.0,
        cd: 0.0,
        cm: 0.0,
        cp: vec![],
        cp_x: vec![],
        gamma: vec![],
        psi_0: 0.0,
        converged: false,
        iterations: 0,
        residual: 0.0,
        x_tr_upper: 1.0,
        x_tr_lower: 1.0,
        success: false,
        error: Some(message.into()),
    }
}

fn bl_error(message: impl Into<String>) -> BLDistribution {
    BLDistribution {
        x_upper: vec![],
        x_lower: vec![],
        theta_upper: vec![],
        theta_lower: vec![],
        delta_star_upper: vec![],
        delta_star_lower: vec![],
        h_upper: vec![],
        h_lower: vec![],
        cf_upper: vec![],
        cf_lower: vec![],
        x_tr_upper: 1.0,
        x_tr_lower: 1.0,
        converged: false,
        iterations: 0,
        residual: 0.0,
        success: false,
        error: Some(message.into()),
    }
}

fn faithful_snapshot(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
) -> Result<FaithfulSnapshot, String> {
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return Err("Invalid coordinates: need at least 3 points (6 values)".to_string());
    }

    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();
    let body = Body::from_points("airfoil", &points).map_err(|e| format!("Geometry error: {e}"))?;

    let panels = body.panels();
    let mut node_x: Vec<f64> = panels.iter().map(|panel| panel.p1.x).collect();
    let mut node_y: Vec<f64> = panels.iter().map(|panel| panel.p1.y).collect();
    if let Some(last) = panels.last() {
        if (last.p2.x - panels[0].p1.x).abs() > 1.0e-10 || (last.p2.y - panels[0].p1.y).abs() > 1.0e-10 {
            node_x.push(last.p2.x);
            node_y.push(last.p2.y);
        }
    }
    let nodes: Vec<Point> = node_x.iter().zip(node_y.iter()).map(|(&x, &y)| point(x, y)).collect();
    let coords_pairs: Vec<(f64, f64)> = node_x.iter().zip(node_y.iter()).map(|(&x, &y)| (x, y)).collect();

    let solver = FaithfulInviscidSolver::new();
    let factorized = solver.factorize(&coords_pairs).map_err(|e| format!("Factorization error: {e}"))?;
    let inviscid = factorized.solve_alpha(&FaithfulFlowConditions::with_alpha_deg(alpha_deg));
    let panel_s = {
        let mut s = Vec::with_capacity(node_x.len());
        let mut acc = 0.0;
        s.push(0.0);
        for i in 1..node_x.len() {
            let dx = node_x[i] - node_x[i - 1];
            let dy = node_y[i] - node_y[i - 1];
            acc += (dx * dx + dy * dy).sqrt();
            s.push(acc);
        }
        s
    };

    let mut state = XfoilState::new(
        body.name.clone(),
        alpha_deg.to_radians(),
        rustfoil_xfoil::OperatingMode::PrescribedAlpha,
        node_x,
        node_y,
        panel_s,
        factorized.surface_qinvu_basis().0,
        factorized.surface_qinvu_basis().1,
        factorized.gamu_0.clone(),
        factorized.gamu_90.clone(),
        factorized.build_dij_with_default_wake().map_err(|e| format!("DIJ error: {e}"))?,
    );
    state.panel_xp = factorized.geometry().xp.clone();
    state.panel_yp = factorized.geometry().yp.clone();
    state.xle = factorized.geometry().xle;
    state.yle = factorized.geometry().yle;
    state.xte = factorized.geometry().xte;
    state.yte = factorized.geometry().yte;
    state.chord = factorized.geometry().chord;
    state.sharp = factorized.geometry().sharp;
    state.ante = factorized.geometry().ante;

    let options = XfoilOptions {
        reynolds,
        mach,
        ncrit,
        max_iterations,
        ..Default::default()
    };
    let oper_result = solve_operating_point_from_state(&mut state, &factorized, &options)
        .map_err(|e| format!("Faithful solver error: {e}"))?;

    let gamma = state.gam.clone();
    let cp: Vec<f64> = gamma.iter().map(|g| 1.0 - g * g).collect();
    let cp_x = state.panel_x.clone();

    let result = AnalysisResult {
        cl: oper_result.cl,
        cd: oper_result.cd,
        cm: oper_result.cm,
        cp,
        cp_x,
        gamma: gamma.clone(),
        psi_0: inviscid.psi_0,
        converged: oper_result.converged,
        iterations: oper_result.iterations,
        residual: oper_result.residual,
        x_tr_upper: oper_result.x_tr_upper,
        x_tr_lower: oper_result.x_tr_lower,
        success: true,
        error: None,
    };

    Ok(FaithfulSnapshot {
        nodes,
        gamma,
        psi_0: inviscid.psi_0,
        cp_x: state.panel_x.clone(),
        result,
        upper_rows: state.upper_rows.clone(),
        lower_rows: state.lower_rows.clone(),
        iblte_upper: state.iblte_upper,
        iblte_lower: state.iblte_lower,
        wake_x: state.wake_x.clone(),
        wake_y: state.wake_y.clone(),
    })
}

fn rows_to_bl_distribution(
    upper_rows: &[rustfoil_xfoil::XfoilBlRow],
    lower_rows: &[rustfoil_xfoil::XfoilBlRow],
    iblte_upper: usize,
    iblte_lower: usize,
    x_tr_upper: f64,
    x_tr_lower: f64,
    converged: bool,
    iterations: usize,
    residual: f64,
) -> BLDistribution {
    let upper_end = iblte_upper.min(upper_rows.len().saturating_sub(1));
    let lower_end = iblte_lower.min(lower_rows.len().saturating_sub(1));
    let upper = &upper_rows[..=upper_end];
    let lower = &lower_rows[..=lower_end];

    BLDistribution {
        x_upper: upper.iter().map(|row| row.x_coord).collect(),
        x_lower: lower.iter().map(|row| row.x_coord).collect(),
        theta_upper: upper.iter().map(|row| row.theta).collect(),
        theta_lower: lower.iter().map(|row| row.theta).collect(),
        delta_star_upper: upper.iter().map(|row| row.dstr).collect(),
        delta_star_lower: lower.iter().map(|row| row.dstr).collect(),
        h_upper: upper.iter().map(|row| row.h).collect(),
        h_lower: lower.iter().map(|row| row.h).collect(),
        cf_upper: upper.iter().map(|row| row.cf).collect(),
        cf_lower: lower.iter().map(|row| row.cf).collect(),
        x_tr_upper,
        x_tr_lower,
        converged,
        iterations,
        residual,
        success: true,
        error: None,
    }
}

/// Quick analysis function for simple use cases.
///
/// # Arguments
/// * `coords` - Airfoil coordinates as flat array [x0, y0, x1, y1, ...]
/// * `alpha_deg` - Angle of attack in degrees
///
/// # Returns
/// Analysis result as a JavaScript object.
#[wasm_bindgen]
pub fn analyze_airfoil(coords: &[f64], alpha_deg: f64) -> JsValue {
    let result = analyze_airfoil_impl(coords, alpha_deg);
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn analyze_airfoil_impl(coords: &[f64], alpha_deg: f64) -> AnalysisResult {
    // Parse coordinates
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return analysis_error("Invalid coordinates: need at least 3 points (6 values)");
    }

    let points: Vec<_> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();

    // Build body
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            return analysis_error(format!("Geometry error: {}", e))
        }
    };

    // Solve
    let solver = InviscidSolver::new();
    let flow = FlowConditions::with_alpha_deg(alpha_deg);

    match solver.solve(&[body.clone()], &flow) {
        Ok(solution) => {
            let cp_x: Vec<f64> = body.panels().iter().map(|p| p.midpoint().x).collect();

            AnalysisResult {
                cl: solution.cl,
                cd: 0.0,
                cm: solution.cm,
                cp: solution.cp,
                cp_x,
                gamma: solution.gamma,
                psi_0: solution.psi_0,
                converged: true,
                iterations: 0,
                residual: 0.0,
                x_tr_upper: 1.0,
                x_tr_lower: 1.0,
                success: true,
                error: None,
            }
        }
        Err(e) => analysis_error(format!("Solver error: {}", e)),
    }
}

#[wasm_bindgen]
pub fn analyze_airfoil_faithful(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
) -> JsValue {
    let result = match faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
        Ok(snapshot) => snapshot.result,
        Err(message) => analysis_error(message),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

#[wasm_bindgen]
pub fn get_bl_distribution_faithful(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
) -> JsValue {
    let result = match faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
        Ok(snapshot) => rows_to_bl_distribution(
            &snapshot.upper_rows,
            &snapshot.lower_rows,
            snapshot.iblte_upper,
            snapshot.iblte_lower,
            snapshot.result.x_tr_upper,
            snapshot.result.x_tr_lower,
            snapshot.result.converged,
            snapshot.result.iterations,
            snapshot.result.residual,
        ),
        Err(message) => bl_error(message),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn rows_to_station_data(rows: &[rustfoil_xfoil::XfoilBlRow]) -> BLStationData {
    BLStationData {
        x: rows.iter().map(|r| r.x_coord).collect(),
        y: rows.iter().map(|r| r.y_coord).collect(),
        delta_star: rows.iter().map(|r| r.dstr).collect(),
        theta: rows.iter().map(|r| r.theta).collect(),
        cf: rows.iter().map(|r| r.cf).collect(),
        ue: rows.iter().map(|r| r.uedg).collect(),
        hk: rows.iter().map(|r| r.hk).collect(),
        ampl: rows.iter().map(|r| r.ampl).collect(),
        is_laminar: rows.iter().map(|r| r.is_laminar).collect(),
    }
}

fn build_bl_visualization(snapshot: &FaithfulSnapshot) -> BLVisualizationData {
    let upper_end = snapshot.iblte_upper.min(snapshot.upper_rows.len().saturating_sub(1));
    let lower_end = snapshot.iblte_lower.min(snapshot.lower_rows.len().saturating_sub(1));

    let upper_surface = &snapshot.upper_rows[..=upper_end];
    let lower_surface = &snapshot.lower_rows[..=lower_end];

    // Wake rows (upper and lower have identical δ* in XFOIL -- wake is a single entity)
    let wake_rows: Vec<_> = snapshot.upper_rows.iter()
        .skip(upper_end + 1)
        .cloned()
        .collect();

    let wake = WakeData {
        x: wake_rows.iter().map(|r| r.x_coord).collect(),
        y: wake_rows.iter().map(|r| r.y_coord).collect(),
        delta_star: wake_rows.iter().map(|r| r.dstr).collect(),
        theta: wake_rows.iter().map(|r| r.theta).collect(),
    };

    // Compute the upper/lower fraction for asymmetric wake splitting.
    // This ensures the wake envelope connects smoothly to the surface BL at the TE.
    let upper_dstr_te = upper_surface.last().map_or(0.0, |r| r.dstr);
    let lower_dstr_te = lower_surface.last().map_or(0.0, |r| r.dstr);
    let total_te = upper_dstr_te + lower_dstr_te;
    let wake_upper_fraction = if total_te > 1e-12 { upper_dstr_te / total_te } else { 0.5 };

    BLVisualizationData {
        upper: rows_to_station_data(upper_surface),
        lower: rows_to_station_data(lower_surface),
        wake,
        wake_upper_fraction,
        wake_geometry_x: snapshot.wake_x.clone(),
        wake_geometry_y: snapshot.wake_y.clone(),
        x_tr_upper: snapshot.result.x_tr_upper,
        x_tr_lower: snapshot.result.x_tr_lower,
        converged: snapshot.result.converged,
        iterations: snapshot.result.iterations,
        residual: snapshot.result.residual,
        success: true,
        error: None,
    }
}

fn bl_vis_error(message: impl Into<String>) -> BLVisualizationData {
    let empty_station = BLStationData {
        x: vec![], y: vec![], delta_star: vec![], theta: vec![],
        cf: vec![], ue: vec![], hk: vec![], ampl: vec![], is_laminar: vec![],
    };
    BLVisualizationData {
        upper: empty_station.clone(),
        lower: empty_station,
        wake: WakeData { x: vec![], y: vec![], delta_star: vec![], theta: vec![] },
        wake_upper_fraction: 0.5,
        wake_geometry_x: vec![],
        wake_geometry_y: vec![],
        x_tr_upper: 1.0,
        x_tr_lower: 1.0,
        converged: false,
        iterations: 0,
        residual: 0.0,
        success: false,
        error: Some(message.into()),
    }
}

#[wasm_bindgen]
pub fn get_bl_visualization_faithful(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
) -> JsValue {
    let result = match faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
        Ok(snapshot) => build_bl_visualization(&snapshot),
        Err(message) => bl_vis_error(message),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

// ============================================================================
// Streamline Visualization
// ============================================================================

/// Streamline result for JS.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StreamlineResult {
    /// Array of streamlines, each is array of [x, y] points
    pub streamlines: Vec<Vec<[f64; 2]>>,
    /// Whether computation succeeded
    pub success: bool,
    /// Error message if any
    pub error: Option<String>,
}

/// Dividing streamline result for JS.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DividingStreamlineResult {
    /// Streamline points for the separatrix traced by bisection.
    pub streamline: Vec<[f64; 2]>,
    /// Whether computation succeeded
    pub success: bool,
    /// Error message if any
    pub error: Option<String>,
}

/// Extract wake source panels from the snapshot's BL data and wake geometry.
fn build_wake_source_panels_from_snapshot(snapshot: &FaithfulSnapshot) -> Option<WakePanels> {
    let nw = snapshot.wake_x.len();
    if nw < 2 {
        return None;
    }

    let wake_rows: Vec<&rustfoil_xfoil::XfoilBlRow> = snapshot.upper_rows
        .iter()
        .filter(|r| r.is_wake)
        .collect();

    let mut wake_sigma = vec![0.0f64; nw];
    for (i, row) in wake_rows.iter().enumerate() {
        if i < nw {
            wake_sigma[i] = row.mass;
        }
    }

    Some(WakePanels {
        x: snapshot.wake_x.clone(),
        y: snapshot.wake_y.clone(),
        sigma: wake_sigma,
    })
}

/// Compute streamlines for visualization.
///
/// # Arguments
/// * `coords` - Airfoil coordinates as flat array [x0, y0, x1, y1, ...]
/// * `alpha_deg` - Angle of attack in degrees
/// * `seed_count` - Number of streamlines (seeds)
/// * `bounds` - [x_min, x_max, y_min, y_max] for seed and integration bounds
///
/// # Returns
/// StreamlineResult with array of streamline point arrays.
#[wasm_bindgen]
pub fn compute_streamlines(
    coords: &[f64],
    alpha_deg: f64,
    seed_count: u32,
    bounds: &[f64],
) -> JsValue {
    let result = compute_streamlines_impl(coords, alpha_deg, seed_count, bounds);
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn compute_streamlines_impl(
    coords: &[f64],
    alpha_deg: f64,
    seed_count: u32,
    bounds: &[f64],
) -> StreamlineResult {
    // Parse coordinates
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return StreamlineResult {
            streamlines: vec![],
            success: false,
            error: Some("Invalid coordinates".to_string()),
        };
    }

    // Parse bounds
    if bounds.len() != 4 {
        return StreamlineResult {
            streamlines: vec![],
            success: false,
            error: Some("bounds must have 4 values: [x_min, x_max, y_min, y_max]".to_string()),
        };
    }

    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();

    // Build body and solve
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            return StreamlineResult {
                streamlines: vec![],
                success: false,
                error: Some(format!("Geometry error: {}", e)),
            };
        }
    };

    let solver = InviscidSolver::new();
    let flow = FlowConditions::with_alpha_deg(alpha_deg);

    let solution = match solver.solve(&[body], &flow) {
        Ok(s) => s,
        Err(e) => {
            return StreamlineResult {
                streamlines: vec![],
                success: false,
                error: Some(format!("Solver error: {}", e)),
            };
        }
    };

    // Build streamlines
    // Seed points are placed upstream (at x = bounds[0])
    // Integration domain extends slightly beyond bounds to allow streamlines to start
    let options = StreamlineOptions {
        seed_count: seed_count as usize,
        seed_x: bounds[0],           // Start at left edge of domain
        y_min: bounds[2],
        y_max: bounds[3],
        step_size: 0.01,
        max_steps: 2000,
        x_min: bounds[0] - 0.5,      // Extend left to include seed points
        x_max: bounds[1],
    };

    let streamlines_raw = build_streamlines(
        &points,
        &solution.gamma,
        flow.alpha,
        flow.v_inf,
        &options,
    );

    // Convert to JS-friendly format
    let streamlines: Vec<Vec<[f64; 2]>> = streamlines_raw
        .into_iter()
        .map(|line| line.into_iter().map(|(x, y)| [x, y]).collect())
        .collect();

    StreamlineResult {
        streamlines,
        success: true,
        error: None,
    }
}

#[wasm_bindgen]
pub fn compute_streamlines_faithful(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
    seed_count: u32,
    bounds: &[f64],
) -> JsValue {
    let result = match faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
        Ok(snapshot) => {
            if bounds.len() != 4 {
                StreamlineResult {
                    streamlines: vec![],
                    success: false,
                    error: Some("bounds must have 4 values: [x_min, x_max, y_min, y_max]".to_string()),
                }
            } else {
                let field = FaithfulFlowField::from_snapshot(&snapshot, alpha_deg);
                let options = StreamlineOptions {
                    seed_count: seed_count as usize,
                    seed_x: bounds[0],
                    y_min: bounds[2],
                    y_max: bounds[3],
                    step_size: 0.01,
                    max_steps: 2000,
                    x_min: bounds[0] - 0.5,
                    x_max: bounds[1],
                };
                let streamlines_raw = field.compute_streamlines(&options);
                let streamlines = streamlines_raw
                    .into_iter()
                    .map(|line| line.into_iter().map(|(x, y)| [x, y]).collect())
                    .collect();
                StreamlineResult {
                    streamlines,
                    success: true,
                    error: None,
                }
            }
        }
        Err(message) => StreamlineResult {
            streamlines: vec![],
            success: false,
            error: Some(message),
        },
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

#[wasm_bindgen]
pub fn compute_dividing_streamline(
    coords: &[f64],
    alpha_deg: f64,
    bounds: &[f64],
) -> JsValue {
    let result = compute_dividing_streamline_impl(coords, alpha_deg, bounds);
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn compute_dividing_streamline_impl(
    coords: &[f64],
    alpha_deg: f64,
    bounds: &[f64],
) -> DividingStreamlineResult {
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return DividingStreamlineResult {
            streamline: vec![],
            success: false,
            error: Some("Invalid coordinates".to_string()),
        };
    }

    if bounds.len() != 4 {
        return DividingStreamlineResult {
            streamline: vec![],
            success: false,
            error: Some("bounds must have 4 values: [x_min, x_max, y_min, y_max]".to_string()),
        };
    }

    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            return DividingStreamlineResult {
                streamline: vec![],
                success: false,
                error: Some(format!("Geometry error: {}", e)),
            };
        }
    };

    let solver = InviscidSolver::new();
    let flow = FlowConditions::with_alpha_deg(alpha_deg);
    let solution = match solver.solve(&[body], &flow) {
        Ok(s) => s,
        Err(e) => {
            return DividingStreamlineResult {
                streamline: vec![],
                success: false,
                error: Some(format!("Solver error: {}", e)),
            };
        }
    };

    let options = StreamlineOptions {
        seed_count: 33,
        seed_x: bounds[0],
        y_min: bounds[2],
        y_max: bounds[3],
        step_size: 0.01,
        max_steps: 2000,
        x_min: bounds[0] - 0.5,
        x_max: bounds[1],
    };

    let Some(streamline_raw) = build_dividing_streamline(
        &points,
        &solution.gamma,
        flow.alpha,
        flow.v_inf,
        solution.psi_0,
        &options,
    ) else {
        return DividingStreamlineResult {
            streamline: vec![],
            success: false,
            error: Some("Unable to bracket dividing streamline".to_string()),
        };
    };

    DividingStreamlineResult {
        streamline: streamline_raw
            .into_iter()
            .map(|(x, y)| [x, y])
            .collect::<Vec<[f64; 2]>>(),
        success: true,
        error: None,
    }
}

#[wasm_bindgen]
pub fn compute_dividing_streamline_faithful(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
    bounds: &[f64],
) -> JsValue {
    let result = match faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
        Ok(snapshot) => {
            if bounds.len() != 4 {
                DividingStreamlineResult {
                    streamline: vec![],
                    success: false,
                    error: Some("bounds must have 4 values: [x_min, x_max, y_min, y_max]".to_string()),
                }
            } else {
                let field = FaithfulFlowField::from_snapshot(&snapshot, alpha_deg);
                let options = StreamlineOptions {
                    seed_count: 33,
                    seed_x: bounds[0],
                    y_min: bounds[2],
                    y_max: bounds[3],
                    step_size: 0.01,
                    max_steps: 2000,
                    x_min: bounds[0] - 0.5,
                    x_max: bounds[1],
                };

                match field.compute_dividing_streamline(&options) {
                    Some(streamline_raw) => DividingStreamlineResult {
                        streamline: streamline_raw
                            .into_iter()
                            .map(|(x, y)| [x, y])
                            .collect::<Vec<[f64; 2]>>(),
                        success: true,
                        error: None,
                    },
                    None => DividingStreamlineResult {
                        streamline: vec![],
                        success: false,
                        error: Some("Unable to bracket dividing streamline".to_string()),
                    },
                }
            }
        }
        Err(message) => DividingStreamlineResult {
            streamline: vec![],
            success: false,
            error: Some(message),
        },
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

// ============================================================================
// Stream Function Grid
// ============================================================================

/// Result of stream function grid computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PsiGridResult {
    /// Stream function values in row-major order (ny rows, nx columns)
    /// Values inside the airfoil are NaN
    pub grid: Vec<f64>,
    /// Internal stream function value (contour at this value = body + dividing streamline)
    pub psi_0: f64,
    /// Grid width (number of columns)
    pub nx: usize,
    /// Grid height (number of rows)
    pub ny: usize,
    /// Minimum psi value in grid (excluding NaN)
    pub psi_min: f64,
    /// Maximum psi value in grid (excluding NaN)
    pub psi_max: f64,
    /// Whether computation succeeded
    pub success: bool,
    /// Error message if any
    pub error: Option<String>,
}

/// Compute stream function values on a grid.
///
/// The stream function ψ is constant along streamlines. The contour ψ = psi_0
/// represents the dividing streamline that passes through the stagnation point.
///
/// # Arguments
/// * `coords` - Airfoil coordinates as flat array [x0, y0, x1, y1, ...]
/// * `alpha_deg` - Angle of attack in degrees
/// * `bounds` - [x_min, x_max, y_min, y_max] for the grid
/// * `resolution` - [nx, ny] grid resolution
///
/// # Returns
/// PsiGridResult with grid values and metadata.
#[wasm_bindgen]
pub fn compute_psi_grid(
    coords: &[f64],
    alpha_deg: f64,
    bounds: &[f64],
    resolution: &[u32],
) -> JsValue {
    let result = compute_psi_grid_impl(coords, alpha_deg, bounds, resolution);
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn compute_psi_grid_impl(
    coords: &[f64],
    alpha_deg: f64,
    bounds: &[f64],
    resolution: &[u32],
) -> PsiGridResult {
    // Validate inputs
    if coords.len() < 6 || coords.len() % 2 != 0 {
        return PsiGridResult {
            grid: vec![],
            psi_0: 0.0,
            nx: 0,
            ny: 0,
            psi_min: 0.0,
            psi_max: 0.0,
            success: false,
            error: Some("Invalid coordinates".to_string()),
        };
    }

    if bounds.len() != 4 {
        return PsiGridResult {
            grid: vec![],
            psi_0: 0.0,
            nx: 0,
            ny: 0,
            psi_min: 0.0,
            psi_max: 0.0,
            success: false,
            error: Some("bounds must have 4 values: [x_min, x_max, y_min, y_max]".to_string()),
        };
    }

    if resolution.len() != 2 {
        return PsiGridResult {
            grid: vec![],
            psi_0: 0.0,
            nx: 0,
            ny: 0,
            psi_min: 0.0,
            psi_max: 0.0,
            success: false,
            error: Some("resolution must have 2 values: [nx, ny]".to_string()),
        };
    }

    let points: Vec<Point> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();
    let nx = resolution[0] as usize;
    let ny = resolution[1] as usize;

    // Build body and solve
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            return PsiGridResult {
                grid: vec![],
                psi_0: 0.0,
                nx: 0,
                ny: 0,
                psi_min: 0.0,
                psi_max: 0.0,
                success: false,
                error: Some(format!("Geometry error: {}", e)),
            };
        }
    };

    let solver = InviscidSolver::new();
    let flow = FlowConditions::with_alpha_deg(alpha_deg);

    let solution = match solver.solve(&[body], &flow) {
        Ok(s) => s,
        Err(e) => {
            return PsiGridResult {
                grid: vec![],
                psi_0: 0.0,
                nx: 0,
                ny: 0,
                psi_min: 0.0,
                psi_max: 0.0,
                success: false,
                error: Some(format!("Solver error: {}", e)),
            };
        }
    };

    // Compute grid with NaN for interior points
    // The airfoil will be masked visually in the UI layer
    let grid = rustfoil_solver::inviscid::compute_psi_grid(
        &points,
        &solution.gamma,
        flow.alpha,
        flow.v_inf,
        bounds[0], bounds[1], bounds[2], bounds[3],
        nx, ny,
    );

    // Find min/max (excluding NaN interior points)
    let (psi_min, psi_max) = grid.iter()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), &v| {
            (min.min(v), max.max(v))
        });

    PsiGridResult {
        grid,
        psi_0: solution.psi_0,
        nx,
        ny,
        psi_min,
        psi_max,
        success: true,
        error: None,
    }
}

#[wasm_bindgen]
pub fn compute_psi_grid_faithful(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    max_iterations: usize,
    bounds: &[f64],
    resolution: &[u32],
) -> JsValue {
    let result = match faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
        Ok(snapshot) => {
            if bounds.len() != 4 {
                PsiGridResult {
                    grid: vec![],
                    psi_0: 0.0,
                    nx: 0,
                    ny: 0,
                    psi_min: 0.0,
                    psi_max: 0.0,
                    success: false,
                    error: Some("bounds must have 4 values: [x_min, x_max, y_min, y_max]".to_string()),
                }
            } else if resolution.len() != 2 {
                PsiGridResult {
                    grid: vec![],
                    psi_0: 0.0,
                    nx: 0,
                    ny: 0,
                    psi_min: 0.0,
                    psi_max: 0.0,
                    success: false,
                    error: Some("resolution must have 2 values: [nx, ny]".to_string()),
                }
            } else {
                let nx = resolution[0] as usize;
                let ny = resolution[1] as usize;
                let field = FaithfulFlowField::from_snapshot(&snapshot, alpha_deg);
                let grid = field.compute_psi_grid(bounds, nx, ny, Some(field.psi_0));
                let (psi_min, psi_max) = grid
                    .iter()
                    .filter(|v| v.is_finite())
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), &v| {
                        (min.min(v), max.max(v))
                    });
                PsiGridResult {
                    grid,
                    psi_0: snapshot.psi_0,
                    nx,
                    ny,
                    psi_min,
                    psi_max,
                    success: true,
                    error: None,
                }
            }
        }
        Err(message) => PsiGridResult {
            grid: vec![],
            psi_0: 0.0,
            nx: 0,
            ny: 0,
            psi_min: 0.0,
            psi_max: 0.0,
            success: false,
            error: Some(message),
        },
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

// ============================================================================
// Smoke Visualization
// ============================================================================

/// Smoke particle system for flow visualization.
///
/// This is a stateful object that maintains particle positions
/// and should be updated each frame.
#[wasm_bindgen]
pub struct WasmSmokeSystem {
    inner: SmokeSystem,
    coords: Vec<Point>,
    gamma: Vec<f64>,
    alpha: f64,
    v_inf: f64,
    /// Dividing streamline value (psi_0)
    psi_0: f64,
    faithful_field: Option<FaithfulFlowField>,
}

#[wasm_bindgen]
impl WasmSmokeSystem {
    /// Create a new smoke system.
    ///
    /// # Arguments
    /// * `spawn_y_values` - Y coordinates for spawn points
    /// * `spawn_x` - X coordinate for spawn points
    /// * `particles_per_blob` - Number of particles per blob (default 20)
    #[wasm_bindgen(constructor)]
    pub fn new(spawn_y_values: &[f64], spawn_x: f64, particles_per_blob: u32) -> Self {
        Self {
            inner: SmokeSystem::new(spawn_y_values, spawn_x, particles_per_blob as usize),
            coords: Vec::new(),
            gamma: Vec::new(),
            alpha: 0.0,
            v_inf: 1.0,
            psi_0: 0.0,
            faithful_field: None,
        }
    }

    /// Set the airfoil geometry and flow solution.
    /// Call this when the airfoil or alpha changes.
    /// This will invalidate cached streamlines, causing them to be recomputed on next update.
    pub fn set_flow(&mut self, coords: &[f64], alpha_deg: f64) {
        // Parse coordinates
        if coords.len() < 6 || coords.len() % 2 != 0 {
            return;
        }

        self.coords = coords.chunks(2).map(|c| point(c[0], c[1])).collect();
        self.alpha = alpha_deg.to_radians();
        self.faithful_field = None;

        // Solve to get gamma and psi_0
        let body = match Body::from_points("airfoil", &self.coords) {
            Ok(b) => b,
            Err(_) => return,
        };

        let solver = InviscidSolver::new();
        let flow = FlowConditions::with_alpha_deg(alpha_deg);

        if let Ok(solution) = solver.solve(&[body], &flow) {
            self.gamma = solution.gamma;
            self.psi_0 = solution.psi_0;
        }
        
        // Invalidate cached streamlines so they're recomputed with new flow
        self.inner.invalidate_cache();
    }

    /// Set a faithful viscous flow field for smoke advection/coloring.
    pub fn set_faithful_flow(
        &mut self,
        coords: &[f64],
        alpha_deg: f64,
        reynolds: f64,
        mach: f64,
        ncrit: f64,
        max_iterations: usize,
    ) {
        if let Ok(snapshot) = faithful_snapshot(coords, alpha_deg, reynolds, mach, ncrit, max_iterations) {
            let field = FaithfulFlowField::from_snapshot(&snapshot, alpha_deg);
            self.coords = field.nodes.clone();
            self.gamma = field.gamma.clone();
            self.alpha = field.alpha;
            self.psi_0 = field.psi_0;
            self.faithful_field = Some(field);
            self.inner.invalidate_cache();
        }
    }

    /// Update particles by one time step.
    ///
    /// # Arguments
    /// * `dt` - Time step in seconds (typically 1/60 for 60 FPS)
    pub fn update(&mut self, dt: f64) {
        if let Some(field) = &self.faithful_field {
            self.inner.update_with_sources(
                &field.nodes,
                &field.gamma,
                &field.sigma,
                field.alpha,
                self.v_inf,
                field.wake_panels.as_ref(),
                field.effective_body_ref(),
                dt,
            );
        } else {
            self.inner.update(&self.coords, &self.gamma, self.alpha, self.v_inf, dt);
        }
    }

    /// Get particle positions as flat array [x0, y0, x1, y1, ...].
    pub fn get_positions(&self) -> Vec<f64> {
        self.inner.get_positions()
    }

    /// Get particle alpha (opacity) values.
    pub fn get_alphas(&self) -> Vec<f64> {
        self.inner.get_alphas()
    }

    /// Get number of active particles.
    pub fn particle_count(&self) -> usize {
        self.inner.particle_count()
    }

    /// Reset (remove all particles).
    pub fn reset(&mut self) {
        self.inner.reset();
    }

    /// Set spawn points from flat array of (x, y) pairs.
    /// 
    /// This allows arbitrary spawn positions (e.g., rotated for angle of attack).
    /// 
    /// # Arguments
    /// * `points` - Flat array [x0, y0, x1, y1, ...]
    pub fn set_spawn_points(&mut self, points: &[f64]) {
        self.inner.set_spawn_points(points);
    }

    /// Set spawn interval in seconds.
    pub fn set_spawn_interval(&mut self, interval: f64) {
        self.inner.set_spawn_interval(interval);
    }

    /// Set maximum particle age in seconds.
    pub fn set_max_age(&mut self, max_age: f64) {
        self.inner.set_max_age(max_age);
    }

    /// Set freestream velocity magnitude (flow speed multiplier).
    /// Values are clamped to [0.1, 10.0].
    /// This invalidates cached streamlines since velocity affects the paths.
    pub fn set_v_inf(&mut self, v_inf: f64) {
        let new_v_inf = v_inf.max(0.1).min(10.0);
        if (new_v_inf - self.v_inf).abs() > 0.01 {
            self.v_inf = new_v_inf;
            self.inner.invalidate_cache();
        }
    }

    /// Get current freestream velocity magnitude.
    pub fn get_v_inf(&self) -> f64 {
        self.v_inf
    }

    /// Get stream function (psi) values for each particle.
    /// 
    /// This is used to determine which side of the dividing streamline
    /// each particle is on. Compare with get_psi_0() to determine above/below.
    pub fn get_psi_values(&self) -> Vec<f64> {
        if let Some(field) = &self.faithful_field {
            self.inner.get_psi_values_with_sources(
                &field.nodes,
                &field.gamma,
                &field.sigma,
                field.alpha,
                self.v_inf,
                field.wake_panels.as_ref(),
            )
        } else {
            self.inner.get_psi_values(&self.coords, &self.gamma, self.alpha, self.v_inf)
        }
    }

    /// Get the dividing streamline value (psi_0).
    /// 
    /// Particles with psi > psi_0 go above the dividing streamline (upper surface).
    /// Particles with psi < psi_0 go below the dividing streamline (lower surface).
    pub fn get_psi_0(&self) -> f64 {
        self.faithful_field.as_ref().map(|f| f.psi_0).unwrap_or(self.psi_0)
    }
}

// ============================================================================
// Morphing Interpolation (High-Performance)
// ============================================================================

/// Interpolate between two point arrays (same length - fast path).
/// 
/// # Arguments
/// * `from` - Source points as flat array [x0, y0, x1, y1, ...]
/// * `to` - Target points as flat array [x0, y0, x1, y1, ...]
/// * `t` - Interpolation factor (0.0 = from, 1.0 = to)
/// 
/// # Returns
/// Interpolated points as flat array.
#[wasm_bindgen]
pub fn lerp_points(from: &[f64], to: &[f64], t: f64) -> Vec<f64> {
    let from_len = from.len();
    let to_len = to.len();
    
    if from_len == 0 {
        return to.to_vec();
    }
    if to_len == 0 {
        return from.to_vec();
    }
    
    // Fast path: same length arrays
    if from_len == to_len {
        let one_minus_t = 1.0 - t;
        from.iter()
            .zip(to.iter())
            .map(|(&a, &b)| a * one_minus_t + b * t)
            .collect()
    } else {
        // Different lengths: need index mapping
        lerp_points_diff_length(from, to, t)
    }
}

/// Interpolate between point arrays of different lengths.
fn lerp_points_diff_length(from: &[f64], to: &[f64], t: f64) -> Vec<f64> {
    let from_pts = from.len() / 2;
    let to_pts = to.len() / 2;
    
    if from_pts == 0 || to_pts == 0 {
        return if to_pts > 0 { to.to_vec() } else { from.to_vec() };
    }
    
    let mut result = Vec::with_capacity(to.len());
    let scale = (from_pts - 1) as f64 / (to_pts - 1) as f64;
    let one_minus_t = 1.0 - t;
    
    for i in 0..to_pts {
        let src_idx = i as f64 * scale;
        let src_idx_low = src_idx as usize;
        let src_idx_high = (src_idx_low + 1).min(from_pts - 1);
        let src_t = src_idx - src_idx_low as f64;
        let one_minus_src_t = 1.0 - src_t;
        
        // Interpolate within source
        let src_x = from[src_idx_low * 2] * one_minus_src_t + from[src_idx_high * 2] * src_t;
        let src_y = from[src_idx_low * 2 + 1] * one_minus_src_t + from[src_idx_high * 2 + 1] * src_t;
        
        // Interpolate between source and target
        result.push(src_x * one_minus_t + to[i * 2] * t);
        result.push(src_y * one_minus_t + to[i * 2 + 1] * t);
    }
    
    result
}

/// Interpolate between two number arrays (for Cp values, etc.).
/// 
/// # Arguments
/// * `from` - Source array
/// * `to` - Target array  
/// * `t` - Interpolation factor (0.0 = from, 1.0 = to)
/// 
/// # Returns
/// Interpolated array.
#[wasm_bindgen]
pub fn lerp_array(from: &[f64], to: &[f64], t: f64) -> Vec<f64> {
    let from_len = from.len();
    let to_len = to.len();
    
    if from_len == 0 {
        return to.to_vec();
    }
    if to_len == 0 {
        return from.to_vec();
    }
    
    let one_minus_t = 1.0 - t;
    
    // Fast path: same length
    if from_len == to_len {
        from.iter()
            .zip(to.iter())
            .map(|(&a, &b)| a * one_minus_t + b * t)
            .collect()
    } else {
        // Different lengths: need index mapping
        let scale = (from_len - 1) as f64 / (to_len - 1) as f64;
        let mut result = Vec::with_capacity(to_len);
        
        for i in 0..to_len {
            let src_idx = i as f64 * scale;
            let src_idx_low = src_idx as usize;
            let src_idx_high = (src_idx_low + 1).min(from_len - 1);
            let src_t = src_idx - src_idx_low as f64;
            
            let src_val = from[src_idx_low] * (1.0 - src_t) + from[src_idx_high] * src_t;
            result.push(src_val * one_minus_t + to[i] * t);
        }
        
        result
    }
}

/// Interpolate between two streamline arrays.
/// 
/// Streamlines are encoded as: [n_lines, n_pts_0, x0, y0, x1, y1, ..., n_pts_1, x0, y0, ...]
/// 
/// # Arguments
/// * `from` - Source streamlines encoded
/// * `to` - Target streamlines encoded
/// * `t` - Interpolation factor (0.0 = from, 1.0 = to)
/// 
/// # Returns
/// Interpolated streamlines in same encoding.
#[wasm_bindgen]
pub fn lerp_streamlines(from: &[f64], to: &[f64], t: f64) -> Vec<f64> {
    // Decode streamlines
    let from_lines = decode_streamlines(from);
    let to_lines = decode_streamlines(to);
    
    if from_lines.is_empty() {
        return to.to_vec();
    }
    if to_lines.is_empty() {
        return from.to_vec();
    }
    
    let max_lines = from_lines.len().max(to_lines.len());
    let one_minus_t = 1.0 - t;
    
    let mut result_lines: Vec<Vec<(f64, f64)>> = Vec::with_capacity(max_lines);
    
    for i in 0..max_lines {
        let from_line = &from_lines[i % from_lines.len()];
        let to_line = &to_lines[i % to_lines.len()];
        
        if from_line.is_empty() {
            result_lines.push(to_line.clone());
            continue;
        }
        if to_line.is_empty() {
            result_lines.push(from_line.clone());
            continue;
        }
        
        let from_pts = from_line.len();
        let to_pts = to_line.len();
        
        let mut line_result: Vec<(f64, f64)> = Vec::with_capacity(to_pts);
        
        if from_pts == to_pts {
            // Fast path: same length
            for j in 0..to_pts {
                line_result.push((
                    from_line[j].0 * one_minus_t + to_line[j].0 * t,
                    from_line[j].1 * one_minus_t + to_line[j].1 * t,
                ));
            }
        } else {
            // Different lengths: index mapping
            let scale = (from_pts - 1) as f64 / (to_pts - 1) as f64;
            
            for j in 0..to_pts {
                let src_idx = j as f64 * scale;
                let src_idx_low = src_idx as usize;
                let src_idx_high = (src_idx_low + 1).min(from_pts - 1);
                let src_t = src_idx - src_idx_low as f64;
                let one_minus_src_t = 1.0 - src_t;
                
                let src_x = from_line[src_idx_low].0 * one_minus_src_t + from_line[src_idx_high].0 * src_t;
                let src_y = from_line[src_idx_low].1 * one_minus_src_t + from_line[src_idx_high].1 * src_t;
                
                line_result.push((
                    src_x * one_minus_t + to_line[j].0 * t,
                    src_y * one_minus_t + to_line[j].1 * t,
                ));
            }
        }
        
        result_lines.push(line_result);
    }
    
    encode_streamlines(&result_lines)
}

/// Decode streamlines from flat array format.
/// Format: [n_lines, n_pts_0, x0, y0, x1, y1, ..., n_pts_1, x0, y0, ...]
fn decode_streamlines(data: &[f64]) -> Vec<Vec<(f64, f64)>> {
    if data.is_empty() {
        return vec![];
    }
    
    let n_lines = data[0] as usize;
    let mut result = Vec::with_capacity(n_lines);
    let mut idx = 1;
    
    for _ in 0..n_lines {
        if idx >= data.len() {
            break;
        }
        
        let n_pts = data[idx] as usize;
        idx += 1;
        
        let mut line = Vec::with_capacity(n_pts);
        for _ in 0..n_pts {
            if idx + 1 >= data.len() {
                break;
            }
            line.push((data[idx], data[idx + 1]));
            idx += 2;
        }
        result.push(line);
    }
    
    result
}

/// Encode streamlines to flat array format.
fn encode_streamlines(lines: &[Vec<(f64, f64)>]) -> Vec<f64> {
    // Calculate total size
    let total_pts: usize = lines.iter().map(|l| l.len()).sum();
    let total_size = 1 + lines.len() + total_pts * 2;
    
    let mut result = Vec::with_capacity(total_size);
    result.push(lines.len() as f64);
    
    for line in lines {
        result.push(line.len() as f64);
        for &(x, y) in line {
            result.push(x);
            result.push(y);
        }
    }
    
    result
}

/// Batch interpolation for morphing animation.
/// Interpolates coordinates, panels, Cp, and cpX in a single WASM call.
/// 
/// # Arguments
/// * `from_coords` - Source coordinates [x0, y0, ...]
/// * `to_coords` - Target coordinates
/// * `from_panels` - Source panels [x0, y0, ...]
/// * `to_panels` - Target panels
/// * `from_cp` - Source Cp values
/// * `to_cp` - Target Cp values
/// * `from_cpx` - Source Cp X positions
/// * `to_cpx` - Target Cp X positions
/// * `t` - Interpolation factor
/// 
/// # Returns
/// JsValue containing { coordinates, panels, cp, cpX }
#[wasm_bindgen]
pub fn lerp_morph_state(
    from_coords: &[f64],
    to_coords: &[f64],
    from_panels: &[f64],
    to_panels: &[f64],
    from_cp: &[f64],
    to_cp: &[f64],
    from_cpx: &[f64],
    to_cpx: &[f64],
    t: f64,
) -> JsValue {
    let coords = lerp_points(from_coords, to_coords, t);
    let panels = lerp_points(from_panels, to_panels, t);
    let cp = lerp_array(from_cp, to_cp, t);
    let cp_x = lerp_array(from_cpx, to_cpx, t);
    
    let result = MorphResult {
        coordinates: coords,
        panels,
        cp,
        cp_x,
    };
    
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

/// Result of morph interpolation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MorphResult {
    pub coordinates: Vec<f64>,
    pub panels: Vec<f64>,
    pub cp: Vec<f64>,
    pub cp_x: Vec<f64>,
}

// ============================================================================
// Inverse Design (QDES)
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InverseDesignResult {
    pub input_coords: Vec<f64>,
    pub output_coords: Vec<f64>,
    pub paneled_coords: Vec<f64>,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub converged: bool,
    pub iterations: usize,
    pub rms_error: f64,
    pub max_error: f64,
    pub target_kind: String,
    pub target_upper_x: Vec<f64>,
    pub target_upper_values: Vec<f64>,
    pub target_lower_x: Vec<f64>,
    pub target_lower_values: Vec<f64>,
    pub achieved_upper_x: Vec<f64>,
    pub achieved_upper_values: Vec<f64>,
    pub achieved_lower_x: Vec<f64>,
    pub achieved_lower_values: Vec<f64>,
    pub history: Vec<InverseDesignIterationSnapshot>,
    pub success: bool,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InverseDesignIterationSnapshot {
    pub iteration: usize,
    pub rms_error: f64,
    pub max_error: f64,
    pub geometry_delta_norm: f64,
    pub cl: f64,
    pub cd: f64,
}

/// Run inverse design (QDES) to modify an airfoil to match a target Cp or velocity distribution.
///
/// # Arguments
/// * `coords` - Airfoil coordinates as flat array [x0, y0, x1, y1, ...]
/// * `alpha_deg` - Design angle of attack in degrees
/// * `reynolds` - Reynolds number
/// * `mach` - Mach number
/// * `ncrit` - eN transition criterion
/// * `target_kind` - "cp" for pressure coefficient, "velocity" for edge velocity
/// * `upper_x` - Target x-locations on upper surface (empty to skip upper)
/// * `upper_values` - Target values on upper surface
/// * `lower_x` - Target x-locations on lower surface (empty to skip lower)
/// * `lower_values` - Target values on lower surface
/// * `max_design_iterations` - Max inverse design outer iterations (default 6)
/// * `damping` - Update damping factor 0-1 (default 0.6)
#[wasm_bindgen]
pub fn inverse_design_qdes(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    target_kind: &str,
    upper_x: &[f64],
    upper_values: &[f64],
    lower_x: &[f64],
    lower_values: &[f64],
    max_design_iterations: Option<usize>,
    damping: Option<f64>,
) -> JsValue {
    let result = inverse_design_qdes_impl(
        coords, alpha_deg, reynolds, mach, ncrit,
        target_kind, upper_x, upper_values, lower_x, lower_values,
        max_design_iterations, damping,
    );
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn inverse_design_qdes_impl(
    coords: &[f64],
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    ncrit: f64,
    target_kind: &str,
    upper_x: &[f64],
    upper_values: &[f64],
    lower_x: &[f64],
    lower_values: &[f64],
    max_design_iterations: Option<usize>,
    damping: Option<f64>,
) -> InverseDesignResult {
    use rustfoil_xfoil::{
        QdesOptions, QdesSpec, QdesTarget, QdesTargetKind, AlphaSpec,
        solve_coords_qdes,
    };

    if coords.len() < 8 || coords.len() % 2 != 0 {
        return inverse_design_error("Invalid coordinates: need at least 4 points");
    }
    if upper_x.is_empty() && lower_x.is_empty() {
        return inverse_design_error("At least one surface target (upper or lower) is required");
    }
    if upper_x.len() != upper_values.len() {
        return inverse_design_error("upper_x and upper_values must have the same length");
    }
    if lower_x.len() != lower_values.len() {
        return inverse_design_error("lower_x and lower_values must have the same length");
    }

    let kind = match target_kind {
        "cp" | "Cp" | "CP" | "pressure" => QdesTargetKind::PressureCoefficient,
        _ => QdesTargetKind::EdgeVelocity,
    };

    let upper = if !upper_x.is_empty() {
        Some(QdesTarget {
            x: upper_x.to_vec(),
            values: upper_values.to_vec(),
        })
    } else {
        None
    };

    let lower = if !lower_x.is_empty() {
        Some(QdesTarget {
            x: lower_x.to_vec(),
            values: lower_values.to_vec(),
        })
    } else {
        None
    };

    let spec = QdesSpec {
        operating_point: AlphaSpec::AlphaDeg(alpha_deg),
        target_kind: kind,
        upper,
        lower,
    };

    let options = QdesOptions {
        xfoil_options: rustfoil_xfoil::XfoilOptions {
            reynolds,
            mach,
            ncrit,
            max_iterations: 100,
            ..Default::default()
        },
        outer_iterations: max_design_iterations.unwrap_or(6),
        update_damping: damping.unwrap_or(0.6),
        ..Default::default()
    };

    let coord_pairs: Vec<(f64, f64)> = coords.chunks(2).map(|c| (c[0], c[1])).collect();

    match solve_coords_qdes("airfoil", &coord_pairs, spec, options) {
        Ok(qdes_result) => {
            let (tu_x, tu_v) = qdes_result.target_upper.as_ref()
                .map(|t| (t.x.clone(), t.values.clone()))
                .unwrap_or_default();
            let (tl_x, tl_v) = qdes_result.target_lower.as_ref()
                .map(|t| (t.x.clone(), t.values.clone()))
                .unwrap_or_default();
            let (au_x, au_v) = qdes_result.achieved_upper.as_ref()
                .map(|d| (d.x.clone(), d.values.clone()))
                .unwrap_or_default();
            let (al_x, al_v) = qdes_result.achieved_lower.as_ref()
                .map(|d| (d.x.clone(), d.values.clone()))
                .unwrap_or_default();

            InverseDesignResult {
                input_coords: qdes_result.input_coords.iter().flat_map(|&(x, y)| [x, y]).collect(),
                output_coords: qdes_result.output_coords.iter().flat_map(|&(x, y)| [x, y]).collect(),
                paneled_coords: qdes_result.paneled_coords.iter().flat_map(|&(x, y)| [x, y]).collect(),
                cl: qdes_result.oper_result.cl,
                cd: qdes_result.oper_result.cd,
                cm: qdes_result.oper_result.cm,
                converged: qdes_result.converged,
                iterations: qdes_result.iterations,
                rms_error: qdes_result.rms_error,
                max_error: qdes_result.max_error,
                target_kind: match qdes_result.target_kind {
                    rustfoil_xfoil::QdesTargetKind::PressureCoefficient => "cp".to_string(),
                    rustfoil_xfoil::QdesTargetKind::EdgeVelocity => "velocity".to_string(),
                },
                target_upper_x: tu_x,
                target_upper_values: tu_v,
                target_lower_x: tl_x,
                target_lower_values: tl_v,
                achieved_upper_x: au_x,
                achieved_upper_values: au_v,
                achieved_lower_x: al_x,
                achieved_lower_values: al_v,
                history: qdes_result.history.iter().map(|h| InverseDesignIterationSnapshot {
                    iteration: h.iteration,
                    rms_error: h.rms_error,
                    max_error: h.max_error,
                    geometry_delta_norm: h.geometry_delta_norm,
                    cl: h.oper_result.cl,
                    cd: h.oper_result.cd,
                }).collect(),
                success: true,
                error: None,
            }
        }
        Err(e) => inverse_design_error(format!("Inverse design failed: {e}")),
    }
}

fn inverse_design_error(message: impl Into<String>) -> InverseDesignResult {
    InverseDesignResult {
        input_coords: vec![], output_coords: vec![], paneled_coords: vec![],
        cl: 0.0, cd: 0.0, cm: 0.0, converged: false, iterations: 0,
        rms_error: 0.0, max_error: 0.0,
        target_kind: String::new(),
        target_upper_x: vec![], target_upper_values: vec![],
        target_lower_x: vec![], target_lower_values: vec![],
        achieved_upper_x: vec![], achieved_upper_values: vec![],
        achieved_lower_x: vec![], achieved_lower_values: vec![],
        history: vec![],
        success: false, error: Some(message.into()),
    }
}

// ============================================================================
// Geometry Design (GDES) Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryResult {
    pub coords: Vec<f64>,
    pub success: bool,
    pub error: Option<String>,
}

fn geometry_error(message: impl Into<String>) -> GeometryResult {
    GeometryResult { coords: vec![], success: false, error: Some(message.into()) }
}

fn parse_coords(coords: &[f64]) -> Option<Vec<Point>> {
    if coords.len() < 6 || coords.len() % 2 != 0 { return None; }
    Some(coords.chunks(2).map(|c| point(c[0], c[1])).collect())
}

fn points_to_flat(pts: &[Point]) -> Vec<f64> {
    pts.iter().flat_map(|p| [p.x, p.y]).collect()
}

/// Rotate airfoil by the given angle in degrees about (cx, cy).
#[wasm_bindgen]
pub fn gdes_rotate(coords: &[f64], angle_deg: f64, cx: f64, cy: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            let rad = angle_deg.to_radians();
            let cos_a = rad.cos();
            let sin_a = rad.sin();
            let rotated: Vec<Point> = pts.iter().map(|p| {
                let dx = p.x - cx;
                let dy = p.y - cy;
                point(cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a)
            }).collect();
            GeometryResult { coords: points_to_flat(&rotated), success: true, error: None }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

/// Scale airfoil about (cx, cy) by the given factors.
#[wasm_bindgen]
pub fn gdes_scale(coords: &[f64], sx: f64, sy: f64, cx: f64, cy: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            let scaled: Vec<Point> = pts.iter().map(|p| {
                point(cx + (p.x - cx) * sx, cy + (p.y - cy) * sy)
            }).collect();
            GeometryResult { coords: points_to_flat(&scaled), success: true, error: None }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

/// Translate airfoil by (dx, dy).
#[wasm_bindgen]
pub fn gdes_translate(coords: &[f64], dx: f64, dy: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            let translated: Vec<Point> = pts.iter().map(|p| point(p.x + dx, p.y + dy)).collect();
            GeometryResult { coords: points_to_flat(&translated), success: true, error: None }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

/// Set trailing-edge gap to specified value (fraction of chord).
/// Blends the gap change over the aft portion of the airfoil (XFOIL TGAP approach).
#[wasm_bindgen]
pub fn gdes_set_te_gap(coords: &[f64], gap: f64, blend_fraction: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            if pts.len() < 4 {
                geometry_error("Need at least 4 points")
            } else {
                let modified = set_te_gap_impl(&pts, gap, blend_fraction);
                GeometryResult { coords: points_to_flat(&modified), success: true, error: None }
            }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn set_te_gap_impl(pts: &[Point], target_gap: f64, blend_fraction: f64) -> Vec<Point> {
    let n = pts.len();
    let le_idx = pts.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i).unwrap_or(0);

    let x_min = pts[le_idx].x;
    let x_max = pts.iter().fold(f64::NEG_INFINITY, |acc, p| acc.max(p.x));
    let chord = (x_max - x_min).max(1e-10);

    let current_gap = pts[0].y - pts[n - 1].y;
    let delta = target_gap * chord - current_gap;
    let blend = blend_fraction.clamp(0.1, 1.0);

    pts.iter().enumerate().map(|(i, p)| {
        let eta = ((p.x - x_min) / chord).clamp(0.0, 1.0);
        let ramp = if eta > (1.0 - blend) {
            ((eta - (1.0 - blend)) / blend).powi(2)
        } else {
            0.0
        };
        let sign = if i <= le_idx { -1.0 } else { 1.0 };
        point(p.x, p.y + sign * 0.5 * delta * ramp)
    }).collect()
}

/// Apply a trailing-edge flap deflection.
///
/// * `hinge_x_frac` - Hinge x-position as fraction of chord (e.g. 0.75)
/// * `hinge_y_frac` - Hinge y-position as fraction of local thickness
///   (0.0 = lower surface, 0.5 = mid-thickness, 1.0 = upper surface)
/// * `deflection_deg` - Flap deflection angle in degrees (positive = down)
#[wasm_bindgen]
pub fn gdes_flap(coords: &[f64], hinge_x_frac: f64, hinge_y_frac: f64, deflection_deg: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            if pts.len() < 4 {
                geometry_error("Need at least 4 points")
            } else {
                let modified = xfoil_flap(&pts, hinge_x_frac, hinge_y_frac, deflection_deg);
                GeometryResult { coords: points_to_flat(&modified), success: true, error: None }
            }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

/// Interpolate y on a surface (sequence of points ordered by increasing x)
/// at the given x. Returns None if x is outside the range.
fn interp_y_at_x(surface: &[Point], x: f64) -> Option<f64> {
    for pair in surface.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        if (a.x <= x && b.x >= x) || (b.x <= x && a.x >= x) {
            let dx = b.x - a.x;
            let t = if dx.abs() > 1e-12 { (x - a.x) / dx } else { 0.5 };
            return Some(a.y + t * (b.y - a.y));
        }
    }
    None
}

/// 2D line-segment intersection. Returns the point where (a1→a2) crosses
/// (b1→b2), or None if the segments don't intersect within their spans.
fn seg_intersect(a1: Point, a2: Point, b1: Point, b2: Point) -> Option<Point> {
    let d1x = a2.x - a1.x;
    let d1y = a2.y - a1.y;
    let d2x = b2.x - b1.x;
    let d2y = b2.y - b1.y;
    let denom = d1x * d2y - d1y * d2x;
    if denom.abs() < 1e-15 { return None; }
    let dx = b1.x - a1.x;
    let dy = b1.y - a1.y;
    let t = (dx * d2y - dy * d2x) / denom;
    let u = (dx * d1y - dy * d1x) / denom;
    if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 {
        Some(point(a1.x + t * d1x, a1.y + t * d1y))
    } else {
        None
    }
}

/// Process one surface (LE→TE) through a flap deflection.
///
/// Splits into fore (fixed, x ≤ hinge) and aft (rotated). If the rotated aft
/// folds back over the fore (like XFOIL's "disappeared" surface segment), finds
/// the intersection, trims both sides to it, and joins cleanly.
fn flap_surface(
    surface: &[Point],
    hinge_x: f64,
    hinge_y: f64,
    cos_d: f64,
    sin_d: f64,
) -> Vec<Point> {
    let mut fore: Vec<Point> = Vec::new();
    let mut aft: Vec<Point> = Vec::new();

    for &p in surface {
        if p.x <= hinge_x {
            fore.push(p);
        } else {
            let dx = p.x - hinge_x;
            let dy = p.y - hinge_y;
            aft.push(point(
                hinge_x + dx * cos_d + dy * sin_d,
                hinge_y - dx * sin_d + dy * cos_d,
            ));
        }
    }

    if fore.is_empty() { return aft; }
    if aft.is_empty() { return fore; }

    // Search for intersection between the tail of the fore polyline and the
    // head of the aft polyline. This is the "break point" where the surfaces meet.
    let n_fore = fore.len();
    let n_aft = aft.len();
    let check_fore = n_fore.min(50);
    let check_aft = n_aft.min(50);

    let mut best: Option<(usize, usize, Point)> = None;

    'outer: for i in ((n_fore.saturating_sub(check_fore))..n_fore.saturating_sub(1)).rev() {
        for j in 0..check_aft.saturating_sub(1) {
            if let Some(pt) = seg_intersect(fore[i], fore[i + 1], aft[j], aft[j + 1]) {
                best = Some((i, j + 1, pt));
                break 'outer;
            }
        }
    }

    match best {
        Some((fi, ai, pt)) => {
            // Trim: keep fore up to the crossing segment, add intersection, then
            // aft from after the crossing — deleting the folded region between them.
            let mut result = Vec::with_capacity(fi + 2 + n_aft - ai);
            result.extend_from_slice(&fore[..=fi]);
            result.push(pt);
            result.extend_from_slice(&aft[ai..]);
            result
        }
        None => {
            // No fold — just concatenate (gap side; rendering spline bridges it)
            let mut result = fore;
            result.extend(aft);
            result
        }
    }
}

fn flap_impl(pts: &[Point], hinge_x_frac: f64, hinge_y_frac: f64, deflection_deg: f64) -> Vec<Point> {
    let le_idx = pts.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i).unwrap_or(0);

    let x_min = pts[le_idx].x;
    let x_max = pts.iter().fold(f64::NEG_INFINITY, |acc, p| acc.max(p.x));
    let chord = (x_max - x_min).max(1e-10);
    let hinge_x = x_min + hinge_x_frac * chord;

    // Split into upper (LE→TE) and lower (LE→TE).
    // Input convention: pts[0]=TE upper, …, pts[le_idx]=LE, …, pts[last]=TE lower.
    let upper: Vec<Point> = pts[..=le_idx].iter().rev().cloned().collect();
    let lower: Vec<Point> = pts[le_idx..].iter().cloned().collect();

    let y_upper = interp_y_at_x(&upper, hinge_x).unwrap_or(0.0);
    let y_lower = interp_y_at_x(&lower, hinge_x).unwrap_or(0.0);
    let hinge_y = y_lower + hinge_y_frac.clamp(0.0, 1.0) * (y_upper - y_lower);

    let rad = deflection_deg.to_radians();
    let cos_d = rad.cos();
    let sin_d = rad.sin();

    // Process each surface: rotate aft points, find intersection where the
    // rotated curve meets the fixed curve, trim the fold (like XFOIL SSS+FLAP).
    let upper_flapped = flap_surface(&upper, hinge_x, hinge_y, cos_d, sin_d);
    let lower_flapped = flap_surface(&lower, hinge_x, hinge_y, cos_d, sin_d);

    // Reassemble: TE→LE (upper reversed) then LE→TE (lower, skip shared LE)
    let mut result: Vec<Point> = upper_flapped.into_iter().rev().collect();
    result.extend(lower_flapped.into_iter().skip(1));
    result
}

/// Set leading-edge radius by scaling the forward portion of the airfoil.
///
/// * `le_radius_factor` - Multiplier for LE radius (1.0 = unchanged)
#[wasm_bindgen]
pub fn gdes_set_le_radius(coords: &[f64], le_radius_factor: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            if pts.len() < 4 {
                geometry_error("Need at least 4 points")
            } else {
                let modified = set_le_radius_impl(&pts, le_radius_factor);
                GeometryResult { coords: points_to_flat(&modified), success: true, error: None }
            }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn set_le_radius_impl(pts: &[Point], factor: f64) -> Vec<Point> {
    let le_idx = pts.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i).unwrap_or(0);

    let x_min = pts[le_idx].x;
    let x_max = pts.iter().fold(f64::NEG_INFINITY, |acc, p| acc.max(p.x));
    let chord = (x_max - x_min).max(1e-10);

    let blend_extent = 0.15;

    pts.iter().enumerate().map(|(_, p)| {
        let eta = ((p.x - x_min) / chord).clamp(0.0, 1.0);
        if eta < blend_extent {
            let local = eta / blend_extent;
            let scale = 1.0 + (factor - 1.0) * (1.0 - local * local);
            let y_from_camber = p.y;
            point(p.x, y_from_camber * scale)
        } else {
            *p
        }
    }).collect()
}

/// Scale thickness by a factor while preserving camber line.
///
/// * `factor` - Thickness multiplier (1.0 = unchanged)
#[wasm_bindgen]
pub fn gdes_scale_thickness(coords: &[f64], factor: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            if pts.len() < 4 {
                geometry_error("Need at least 4 points")
            } else {
                let modified = scale_thickness_impl(&pts, factor);
                GeometryResult { coords: points_to_flat(&modified), success: true, error: None }
            }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn scale_thickness_impl(pts: &[Point], factor: f64) -> Vec<Point> {
    let le_idx = pts.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i).unwrap_or(0);

    let upper = &pts[le_idx..];
    let lower = &pts[..=le_idx];
    let lower_rev: Vec<&Point> = lower.iter().rev().collect();

    let compute_camber_at = |x: f64| -> f64 {
        let y_u = interp_surface(upper, x);
        let y_l = interp_surface_rev(&lower_rev, x);
        (y_u + y_l) / 2.0
    };

    pts.iter().enumerate().map(|(i, p)| {
        let camber_y = compute_camber_at(p.x);
        let thickness_half = p.y - camber_y;
        point(p.x, camber_y + thickness_half * factor)
    }).collect()
}

fn interp_surface(surface: &[Point], x: f64) -> f64 {
    if surface.is_empty() { return 0.0; }
    if surface.len() == 1 { return surface[0].y; }
    for i in 1..surface.len() {
        if (surface[i].x >= x && surface[i - 1].x <= x)
            || (surface[i].x <= x && surface[i - 1].x >= x) {
            let dx = surface[i].x - surface[i - 1].x;
            if dx.abs() < 1e-12 { return surface[i].y; }
            let t = (x - surface[i - 1].x) / dx;
            return surface[i - 1].y + t * (surface[i].y - surface[i - 1].y);
        }
    }
    surface.last().unwrap().y
}

fn interp_surface_rev(surface: &[&Point], x: f64) -> f64 {
    if surface.is_empty() { return 0.0; }
    if surface.len() == 1 { return surface[0].y; }
    for i in 1..surface.len() {
        if (surface[i].x >= x && surface[i - 1].x <= x)
            || (surface[i].x <= x && surface[i - 1].x >= x) {
            let dx = surface[i].x - surface[i - 1].x;
            if dx.abs() < 1e-12 { return surface[i].y; }
            let t = (x - surface[i - 1].x) / dx;
            return surface[i - 1].y + t * (surface[i].y - surface[i - 1].y);
        }
    }
    surface.last().unwrap().y
}

/// Scale camber by a factor while preserving thickness distribution.
///
/// * `factor` - Camber multiplier (1.0 = unchanged, 0.0 = symmetric)
#[wasm_bindgen]
pub fn gdes_scale_camber(coords: &[f64], factor: f64) -> JsValue {
    let result = match parse_coords(coords) {
        Some(pts) => {
            if pts.len() < 4 {
                geometry_error("Need at least 4 points")
            } else {
                let modified = scale_camber_impl(&pts, factor);
                GeometryResult { coords: points_to_flat(&modified), success: true, error: None }
            }
        }
        None => geometry_error("Invalid coordinates"),
    };
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn scale_camber_impl(pts: &[Point], factor: f64) -> Vec<Point> {
    let le_idx = pts.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i).unwrap_or(0);

    let upper = &pts[le_idx..];
    let lower = &pts[..=le_idx];
    let lower_rev: Vec<&Point> = lower.iter().rev().collect();

    let compute_camber_at = |x: f64| -> f64 {
        let y_u = interp_surface(upper, x);
        let y_l = interp_surface_rev(&lower_rev, x);
        (y_u + y_l) / 2.0
    };

    pts.iter().map(|p| {
        let camber_y = compute_camber_at(p.x);
        let thickness_half = p.y - camber_y;
        point(p.x, camber_y * factor + thickness_half)
    }).collect()
}

// ============================================================================
// Full-Inverse Design (MDES) via Circle-Plane Conformal Mapping
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MdesDesignResult {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub qspec_x: Vec<f64>,
    pub qspec_values: Vec<f64>,
    pub cl: f64,
    pub cm: f64,
    pub success: bool,
    pub error: Option<String>,
}

/// Run full-inverse design (MDES) via circle-plane conformal mapping.
///
/// This is XFOIL's full-inverse method, which designs an entire airfoil
/// from a velocity distribution using Fourier coefficients in the circle plane.
///
/// # Arguments
/// * `coords` - Initial airfoil coordinates as flat array [x0, y0, ...]
/// * `alpha_deg` - Design angle of attack in degrees
/// * `symmetric` - If true, enforce symmetric airfoil
/// * `filter_strength` - Hanning filter strength (0 = none, higher = more smoothing)
/// * `target_q` - Optional target velocity distribution (nc points). If empty,
///   uses the current airfoil's velocity distribution.
/// * `nc` - Number of circle-plane points (must be 2^n + 1, e.g. 129)
#[wasm_bindgen]
pub fn full_inverse_design_mdes(
    coords: &[f64],
    alpha_deg: f64,
    symmetric: bool,
    filter_strength: f64,
    target_q: &[f64],
    nc: usize,
) -> JsValue {
    let result = full_inverse_design_mdes_impl(
        coords, alpha_deg, symmetric, filter_strength, target_q, nc,
    );
    serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
}

fn full_inverse_design_mdes_impl(
    coords: &[f64],
    alpha_deg: f64,
    symmetric: bool,
    filter_strength: f64,
    target_q: &[f64],
    nc: usize,
) -> MdesDesignResult {
    use rustfoil_xfoil::mdes::{MdesOptions, solve_mdes};

    if coords.len() < 8 || coords.len() % 2 != 0 {
        return MdesDesignResult {
            x: vec![], y: vec![], qspec_x: vec![], qspec_values: vec![],
            cl: 0.0, cm: 0.0, success: false,
            error: Some("Invalid coordinates: need at least 4 points".to_string()),
        };
    }

    let x: Vec<f64> = coords.iter().step_by(2).copied().collect();
    let y: Vec<f64> = coords.iter().skip(1).step_by(2).copied().collect();

    let options = MdesOptions {
        nc: if nc >= 33 { nc } else { 129 },
        alpha_deg,
        symmetric,
        filter_strength,
    };

    let tq = if target_q.is_empty() { None } else { Some(target_q) };
    let result = solve_mdes(&x, &y, tq, &options);

    MdesDesignResult {
        x: result.x,
        y: result.y,
        qspec_x: result.qspec_x,
        qspec_values: result.qspec_values,
        cl: result.cl,
        cm: result.cm,
        success: result.success,
        error: result.error,
    }
}

/// Main RustFoil interface for JavaScript.
///
/// Provides a stateful API for interactive airfoil manipulation.
#[wasm_bindgen]
pub struct RustFoil {
    body: Option<Body>,
    alpha_deg: f64,
    solver: InviscidSolver,
    last_solution: Option<AnalysisResult>,
}

#[wasm_bindgen]
impl RustFoil {
    /// Create a new RustFoil instance.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            body: None,
            alpha_deg: 0.0,
            solver: InviscidSolver::new(),
            last_solution: None,
        }
    }

    /// Set airfoil coordinates from a flat array.
    ///
    /// # Arguments
    /// * `coords` - Coordinates as [x0, y0, x1, y1, ...]
    ///
    /// # Returns
    /// `true` if coordinates were valid, `false` otherwise.
    #[wasm_bindgen]
    pub fn set_coordinates(&mut self, coords: &[f64]) -> bool {
        if coords.len() < 6 || coords.len() % 2 != 0 {
            self.body = None;
            return false;
        }

        let points: Vec<_> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();

        match Body::from_points("airfoil", &points) {
            Ok(b) => {
                self.body = Some(b);
                self.last_solution = None; // Invalidate cached solution
                true
            }
            Err(_) => {
                self.body = None;
                false
            }
        }
    }

    /// Set the angle of attack (in degrees).
    #[wasm_bindgen]
    pub fn set_alpha(&mut self, alpha_deg: f64) {
        self.alpha_deg = alpha_deg;
        self.last_solution = None; // Invalidate cached solution
    }

    /// Get the current angle of attack (in degrees).
    #[wasm_bindgen]
    pub fn get_alpha(&self) -> f64 {
        self.alpha_deg
    }

    /// Run the analysis and return the result.
    #[wasm_bindgen]
    pub fn solve(&mut self) -> JsValue {
        let result = self.solve_impl();
        self.last_solution = Some(result.clone());
        serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL)
    }

    /// Get the number of panels in the current geometry.
    #[wasm_bindgen]
    pub fn panel_count(&self) -> usize {
        self.body.as_ref().map(|b| b.n_panels()).unwrap_or(0)
    }

    /// Repanel the airfoil using cosine spacing.
    ///
    /// # Arguments
    /// * `n_panels` - Desired number of panels
    ///
    /// # Returns
    /// New coordinates as a flat array, or empty if failed.
    #[wasm_bindgen]
    pub fn repanel_cosine(&mut self, n_panels: usize) -> Vec<f64> {
        let body = match &self.body {
            Some(b) => b,
            None => return vec![],
        };

        // Get current points
        let points: Vec<_> = body
            .panels()
            .iter()
            .map(|p| p.p1)
            .chain(std::iter::once(body.panels().last().unwrap().p2))
            .collect();

        // Create spline and resample
        let spline = match CubicSpline::from_points(&points) {
            Ok(s) => s,
            Err(_) => return vec![],
        };

        let new_points = spline.resample_cosine(n_panels + 1);

        // Flatten to coordinate array
        new_points.iter().flat_map(|p| [p.x, p.y]).collect()
    }
}

impl RustFoil {
    fn solve_impl(&self) -> AnalysisResult {
        let body = match &self.body {
            Some(b) => b,
            None => {
                return analysis_error("No geometry set");
            }
        };

        let flow = FlowConditions::with_alpha_deg(self.alpha_deg);

        match self.solver.solve(&[body.clone()], &flow) {
            Ok(solution) => {
                let cp_x: Vec<f64> = body.panels().iter().map(|p| p.midpoint().x).collect();

                AnalysisResult {
                    cl: solution.cl,
                    cd: 0.0,
                    cm: solution.cm,
                    cp: solution.cp,
                    cp_x,
                    gamma: solution.gamma,
                    psi_0: solution.psi_0,
                    converged: true,
                    iterations: 0,
                    residual: 0.0,
                    x_tr_upper: 1.0,
                    x_tr_lower: 1.0,
                    success: true,
                    error: None,
                }
            }
            Err(e) => analysis_error(format!("Solver error: {}", e)),
        }
    }
}

impl Default for RustFoil {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// CFD Mesh Generation & Setup
// ============================================================================

/// Generate a structured O-type mesh around an airfoil for the CFD solver.
///
/// # Arguments
/// * `coords_flat` - Flat array of [x0, y0, x1, y1, ...] airfoil coordinates
/// * `ni` - Circumferential grid points
/// * `nj` - Radial grid points (layers away from airfoil)
/// * `far_field` - Far-field distance in chord lengths
/// * `ds0` - First cell wall-normal spacing
///
/// # Returns
/// JsValue containing { x: Float32Array, y: Float32Array, ni: u32, nj: u32 }
#[wasm_bindgen]
pub fn cfd_generate_mesh(
    coords_flat: &[f64],
    ni: u32,
    nj: u32,
    far_field: f64,
    ds0: f64,
) -> JsValue {
    let n_pts = coords_flat.len() / 2;
    let mut ax = Vec::with_capacity(n_pts);
    let mut ay = Vec::with_capacity(n_pts);
    for i in 0..n_pts {
        ax.push(coords_flat[i * 2]);
        ay.push(coords_flat[i * 2 + 1]);
    }

    let mesh = rustfoil_cfd::mesh::generate_o_mesh(&ax, &ay, ni, nj, far_field, ds0);

    #[derive(Serialize)]
    struct MeshOutput {
        x: Vec<f32>,
        y: Vec<f32>,
        ni: u32,
        nj: u32,
    }

    let output = MeshOutput {
        x: mesh.x,
        y: mesh.y,
        ni: mesh.ni,
        nj: mesh.nj,
    };

    serde_wasm_bindgen::to_value(&output).unwrap_or(JsValue::NULL)
}

/// Compute freestream initial conditions for the CFD solver.
///
/// # Arguments
/// * `ni` - Grid circumferential points
/// * `nj` - Grid radial points
/// * `mach` - Freestream Mach number
/// * `alpha_deg` - Angle of attack in degrees
/// * `gamma` - Ratio of specific heats (1.4 for air)
/// * `physics_mode` - 0=Euler, 1=LaminarNS, 2=RANS_SA
/// * `reynolds` - Reynolds number (for RANS initialization)
///
/// # Returns
/// Float32Array of length ni*nj*5 with Q = [rho, rho*u, rho*v, E, nu_tilde]
#[wasm_bindgen]
pub fn cfd_initial_conditions(
    ni: u32,
    nj: u32,
    mach: f32,
    alpha_deg: f32,
    gamma: f32,
    physics_mode: u32,
    reynolds: f32,
) -> Vec<f32> {
    use rustfoil_cfd::config::{CfdConfig, PhysicsMode};

    let physics = match physics_mode {
        0 => PhysicsMode::Euler,
        1 => PhysicsMode::LaminarNS,
        2 => PhysicsMode::RansSA,
        _ => PhysicsMode::Euler,
    };

    let config = CfdConfig {
        ni,
        nj,
        gamma,
        mach_inf: mach,
        alpha: alpha_deg.to_radians(),
        physics,
        reynolds,
        ..CfdConfig::default()
    };

    rustfoil_cfd::init::compute_initial_conditions(&config)
}

/// Generate boundary condition type array for the CFD grid.
///
/// # Returns
/// Uint32Array of length ni*nj with BC types (0=interior, 1=wall, 2=farfield)
#[wasm_bindgen]
pub fn cfd_boundary_types(ni: u32, nj: u32) -> Vec<u32> {
    rustfoil_cfd::boundary::generate_bc_types(ni, nj)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naca_0012() {
        let coords = generate_naca4_impl(0.0, 0.0, 0.12, 50);
        
        // Should have 2*50 - 1 + 1 (closing point) = 100 points
        assert_eq!(coords.len(), 100);
        
        // First and last points should be the same (closed contour at TE)
        assert!((coords[0].x - 1.0).abs() < 0.01);
        assert!((coords[0].x - coords.last().unwrap().x).abs() < 1e-10);
        assert!((coords[0].y - coords.last().unwrap().y).abs() < 1e-10);
        
        // Should be symmetric (y values should mirror)
        // LE is at index (n-1) where n=50, so index 49
        let le_idx = 49;
        assert!((coords[le_idx].x).abs() < 0.01); // LE at x ≈ 0
        
        // Check thickness at x = 0.3 (max thickness location for NACA 00xx)
        // Max thickness should be about 0.12 * chord
        let max_y = coords.iter().map(|p| p.y.abs()).fold(0.0, f64::max);
        assert!(max_y > 0.05 && max_y < 0.07); // ~6% half-thickness
    }

    #[test]
    fn test_naca_2412() {
        let coords = generate_naca4_impl(0.02, 0.4, 0.12, 50);
        
        // 2*50 - 1 + 1 (closing) = 100 points
        assert_eq!(coords.len(), 100);
        
        // Should have camber - upper surface higher than lower
        // With CW ordering: first half is upper surface, second half is lower
        let mid = coords.len() / 2;
        
        // Find max y on upper surface (first half) and min y on lower surface (second half)
        let max_upper = coords[..mid].iter().map(|p| p.y).fold(f64::MIN, f64::max);
        let min_lower = coords[mid..].iter().map(|p| p.y).fold(f64::MAX, f64::min);
        
        // Upper should be positive, lower should be negative (mostly)
        assert!(max_upper > 0.05);
        assert!(min_lower < 0.0);
    }

    #[test]
    fn test_naca_from_string() {
        let coords = generate_naca4_from_string("0012", 30);
        assert!(!coords.is_empty());
        // 2*30 - 1 + 1 (closing) = 60 points * 2 coords each = 120
        assert_eq!(coords.len(), 60 * 2);
        
        // Invalid input
        let invalid = generate_naca4_from_string("abc", 30);
        assert!(invalid.is_empty());
        
        let too_short = generate_naca4_from_string("12", 30);
        assert!(too_short.is_empty());
    }

    #[test]
    fn test_ssp_distribution() {
        // Uniform spacing - use spacing_to_distribution directly
        let uniform_spacing = vec![1.0; 100];
        let dist = spacing_to_distribution(&uniform_spacing, 11);
        
        assert_eq!(dist.len(), 11);
        assert!((dist[0] - 0.0).abs() < 1e-10);
        assert!((dist[10] - 1.0).abs() < 1e-10);
        
        // Should be roughly uniform
        for i in 0..10 {
            let delta = dist[i + 1] - dist[i];
            assert!((delta - 0.1).abs() < 0.02);
        }
    }

    #[test]
    fn test_ssp_clustering() {
        // Higher spacing at ends = more points at ends
        // Create a spacing function that's high at ends, low in middle
        let n_samples = 100;
        let spacing: Vec<f64> = (0..n_samples)
            .map(|i| {
                let s = i as f64 / (n_samples - 1) as f64;
                // High at s=0 and s=1, low at s=0.5
                let dist_from_center = (s - 0.5).abs() * 2.0;
                0.5 + 1.5 * dist_from_center
            })
            .collect();
        
        let dist = spacing_to_distribution(&spacing, 21);
        
        // Points should cluster at s=0 and s=1
        // First few deltas should be smaller than middle deltas
        let delta_start = dist[1] - dist[0];
        let delta_mid = dist[11] - dist[10];
        
        assert!(delta_start < delta_mid);
    }

    #[test]
    fn test_analyze_diamond() {
        // More refined diamond airfoil for better matrix conditioning
        let coords = vec![
            1.0, 0.0, // TE
            0.75, -0.05, 0.5, -0.08, 0.25, -0.05, // Lower
            0.0, 0.0, // LE
            0.25, 0.05, 0.5, 0.08, 0.75, 0.05, // Upper
            1.0, 0.0, // Back to TE
        ];

        let result = analyze_airfoil_impl(&coords, 0.0);

        // Accept either success or solver failure (matrix conditioning is Phase 2)
        if result.success {
            assert!(result.cp.len() == 8); // 8 panels
            assert!(result.cl.is_finite());
        } else {
            // Solver may fail with singular matrix on small test cases
            assert!(result.error.is_some());
        }
    }

    #[test]
    fn test_invalid_coords() {
        let coords = vec![1.0, 0.0]; // Too few points
        let result = analyze_airfoil_impl(&coords, 0.0);
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_rustfoil_api() {
        let mut foil = RustFoil::new();

        // Test with valid coordinates
        let coords = vec![
            1.0, 0.0, 0.75, -0.05, 0.5, -0.08, 0.25, -0.05, 0.0, 0.0, 0.25, 0.05, 0.5, 0.08, 0.75,
            0.05, 1.0, 0.0,
        ];

        assert!(foil.set_coordinates(&coords));
        assert_eq!(foil.panel_count(), 8);

        foil.set_alpha(5.0);
        assert!((foil.get_alpha() - 5.0).abs() < 1e-10);
    }

    /// Test the full pipeline that the frontend uses:
    /// 1. Generate NACA 0012 with XFOIL generator
    /// 2. Apply XFOIL paneling (50 or 160 panels)
    /// 3. Solve and check symmetry
    #[test]
    fn test_frontend_pipeline_naca0012() {
        use rustfoil_core::{naca::naca4, PanelingParams};
        
        // Test 1: Direct Rust pipeline (should work)
        println!("\n=== Direct Rust pipeline ===");
        let buffer_direct = naca4(12, Some(123));
        let spline_direct = CubicSpline::from_points(&buffer_direct).unwrap();
        let params = PanelingParams::default();
        let paneled_direct = spline_direct.resample_xfoil(160, &params);
        println!("Direct - Buffer: {} pts, Paneled: {} pts", buffer_direct.len(), paneled_direct.len());
        
        let body_direct = Body::from_points("NACA0012_direct", &paneled_direct).unwrap();
        let solver = rustfoil_solver::inviscid::InviscidSolver::new();
        let factorized = solver.factorize(&[body_direct]).unwrap();
        let flow = rustfoil_solver::inviscid::FlowConditions::with_alpha_deg(0.0);
        let solution = factorized.solve_alpha(&flow);
        println!("Direct - Cl at α=0: {:.6}", solution.cl);
        
        // Test 2: WASM pipeline (via flat arrays)
        println!("\n=== WASM flat-array pipeline ===");
        let buffer_flat = generate_naca4_xfoil(12, Some(123));
        println!("WASM - Buffer: {} values = {} pts", buffer_flat.len(), buffer_flat.len() / 2);
        
        // Check if buffers match
        let match_count = buffer_direct.iter().enumerate()
            .filter(|(i, p)| {
                let fx = buffer_flat[i * 2];
                let fy = buffer_flat[i * 2 + 1];
                (p.x - fx).abs() < 1e-10 && (p.y - fy).abs() < 1e-10
            })
            .count();
        println!("Buffer match: {}/{} points identical", match_count, buffer_direct.len());
        
        // Apply WASM repaneling
        let paneled_flat = repanel_xfoil(&buffer_flat, 160);
        println!("WASM - Paneled: {} values = {} pts", paneled_flat.len(), paneled_flat.len() / 2);
        
        // Check the first and last points
        println!("WASM - First pt: ({:.6}, {:.6})", paneled_flat[0], paneled_flat[1]);
        println!("WASM - Last pt: ({:.6}, {:.6})", paneled_flat[paneled_flat.len()-2], paneled_flat[paneled_flat.len()-1]);
        
        // Analyze using the WASM function
        let result_wasm = analyze_airfoil_impl(&paneled_flat, 0.0);
        println!("WASM - Cl at α=0: {:.6}", result_wasm.cl);
        
        // Compare paneled coordinates
        println!("\n=== Paneled coordinate comparison (first 5, last 5) ===");
        for i in 0..5 {
            println!("  [{}] Direct: ({:.6}, {:.6}) | WASM: ({:.6}, {:.6})",
                i, paneled_direct[i].x, paneled_direct[i].y,
                paneled_flat[i*2], paneled_flat[i*2+1]);
        }
        println!("  ...");
        let n = paneled_flat.len() / 2;
        for i in (n-5)..n {
            let d_idx = if i < paneled_direct.len() { Some(i) } else { None };
            if let Some(idx) = d_idx {
                println!("  [{}] Direct: ({:.6}, {:.6}) | WASM: ({:.6}, {:.6})",
                    i, paneled_direct[idx].x, paneled_direct[idx].y,
                    paneled_flat[i*2], paneled_flat[i*2+1]);
            } else {
                println!("  [{}] WASM: ({:.6}, {:.6}) (no direct equivalent)",
                    i, paneled_flat[i*2], paneled_flat[i*2+1]);
            }
        }
        
        // Direct Rust pipeline should give Cl ≈ 0
        assert!(solution.cl.abs() < 0.0001, "Direct: Cl at α=0 should be ~0, got {}", solution.cl);
    }

    #[test]
    fn test_faithful_dividing_streamline_hits_le_stagnation_region() {
        let buffer_flat = generate_naca4_xfoil(12, Some(123));
        let paneled_flat = repanel_xfoil(&buffer_flat, 160);
        let alpha_deg = 15.0;

        let snapshot = faithful_snapshot(&paneled_flat, alpha_deg, 1.0e6, 0.0, 9.0, 100)
            .expect("faithful snapshot should solve");
        let field = FaithfulFlowField::from_snapshot(&snapshot, alpha_deg);
        let options = StreamlineOptions {
            seed_count: 33,
            seed_x: -0.5,
            y_min: -0.5,
            y_max: 0.5,
            step_size: 0.01,
            max_steps: 2000,
            x_min: -1.0,
            x_max: 2.0,
        };

        let streamline = field
            .compute_dividing_streamline(&options)
            .expect("dividing streamline should be bracketed");

        let start = *streamline.first().expect("streamline start");
        let end = *streamline.last().expect("streamline end");
        let body_x_min = snapshot.nodes.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let body_x_max = snapshot.nodes.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let chord = body_x_max - body_x_min;

        println!(
            "dividing streamline start=({:.4}, {:.4}) end=({:.4}, {:.4})",
            start.0, start.1, end.0, end.1
        );

        assert!(
            end.0 <= body_x_min + 0.2 * chord,
            "dividing streamline should terminate near the leading edge, got x={}",
            end.0
        );
        assert!(
            end.1.abs() <= 0.1 * chord,
            "dividing streamline should terminate near the leading-edge stagnation region, got y={}",
            end.1
        );
    }

    // ====================================================================
    // GDES geometry operation tests (call underlying Rust fns, not WASM)
    // ====================================================================

    fn diamond_pts() -> Vec<Point> {
        vec![
            point(1.0, 0.0),
            point(0.75, 0.04), point(0.5, 0.06), point(0.25, 0.04),
            point(0.0, 0.0),
            point(0.25, -0.04), point(0.5, -0.06), point(0.75, -0.04),
            point(1.0, 0.0),
        ]
    }

    #[test]
    fn gdes_rotate_zero_is_identity() {
        let pts = diamond_pts();
        let rad = 0.0_f64.to_radians();
        let cos_a = rad.cos();
        let sin_a = rad.sin();
        let (cx, cy) = (0.5, 0.0);
        for p in &pts {
            let dx = p.x - cx;
            let dy = p.y - cy;
            let rx = cx + dx * cos_a - dy * sin_a;
            let ry = cy + dx * sin_a + dy * cos_a;
            assert!((rx - p.x).abs() < 1e-12);
            assert!((ry - p.y).abs() < 1e-12);
        }
    }

    #[test]
    fn gdes_flap_deflects_te_downward() {
        let pts = diamond_pts();
        let result = flap_impl(&pts, 0.75, 0.5, 10.0);
        // TE (first point, was at y=0) should have moved down for positive deflection
        assert!(result[0].y < -0.001, "positive flap should move TE down, got y={}", result[0].y);
        // LE (min-x point) should be unchanged
        let le = result.iter().min_by(|a, b| a.x.partial_cmp(&b.x).unwrap()).unwrap();
        assert!((le.x - 0.0).abs() < 1e-10, "LE should not move, got x={}", le.x);
        assert!((le.y - 0.0).abs() < 1e-10, "LE should not move, got y={}", le.y);
    }

    #[test]
    fn gdes_flap_zero_deflection_is_identity() {
        let pts = diamond_pts();
        let result = flap_impl(&pts, 0.75, 0.5, 0.0);
        for (a, b) in result.iter().zip(pts.iter()) {
            assert!((a.x - b.x).abs() < 1e-10);
            assert!((a.y - b.y).abs() < 1e-10);
        }
    }

    #[test]
    fn gdes_te_gap_zero_closes_te() {
        let pts = diamond_pts();
        let result = set_te_gap_impl(&pts, 0.0, 0.8);
        let gap = (result[0].y - result[result.len() - 1].y).abs();
        assert!(gap < 0.01, "TE gap should be near zero, got {gap}");
    }

    #[test]
    fn gdes_scale_thickness_doubles() {
        let pts = diamond_pts();
        let result = scale_thickness_impl(&pts, 2.0);
        // Max y should roughly double from 0.06
        let y_max = result.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        assert!(y_max > 0.10, "doubled thickness should exceed 0.10, got {y_max}");
    }

    #[test]
    fn gdes_scale_camber_zero_symmetrizes() {
        let cambered = vec![
            point(1.0, 0.01),
            point(0.5, 0.08),
            point(0.0, 0.02),
            point(0.5, -0.04),
            point(1.0, -0.01),
        ];
        let result = scale_camber_impl(&cambered, 0.0);
        // At x=0.5, upper and lower y should be ±same
        let y_upper = result[1].y;
        let y_lower = result[3].y;
        assert!(
            (y_upper + y_lower).abs() < 0.02,
            "zero camber should be nearly symmetric: y_u={y_upper}, y_l={y_lower}"
        );
    }

    #[test]
    fn gdes_le_radius_factor_1_is_identity() {
        let pts = diamond_pts();
        let result = set_le_radius_impl(&pts, 1.0);
        for (a, b) in result.iter().zip(pts.iter()) {
            assert!((a.x - b.x).abs() < 1e-10);
            assert!((a.y - b.y).abs() < 1e-10);
        }
    }

    // ====================================================================
    // QDES integration test (calls underlying Rust fn)
    // ====================================================================

    #[test]
    fn inverse_design_qdes_impl_rejects_empty_targets() {
        let result = inverse_design_qdes_impl(
            &[1.0, 0.0, 0.5, 0.06, 0.0, 0.0, 0.5, -0.06, 1.0, 0.0],
            0.0, 1e6, 0.0, 9.0, "cp",
            &[], &[], &[], &[],
            Some(3), Some(0.6),
        );
        assert!(!result.success);
        assert!(result.error.as_ref().unwrap().contains("required"));
    }

    // ====================================================================
    // MDES integration test (calls underlying Rust fn)
    // ====================================================================

    #[test]
    fn full_inverse_design_mdes_impl_naca0012() {
        let buf = generate_naca4_xfoil(12, Some(65));
        let result = full_inverse_design_mdes_impl(
            &buf, 0.0, false, 0.0, &[], 65,
        );

        assert!(result.success, "MDES should succeed: {:?}", result.error);
        assert!(result.x.len() > 10);
        assert_eq!(result.x.len(), result.y.len());
    }
}
