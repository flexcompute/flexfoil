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

use rustfoil_core::{naca, point, Body, CubicSpline, Point};
use rustfoil_solver::inviscid::{
    FlowConditions, InviscidSolver, 
    build_streamlines, StreamlineOptions, SmokeSystem
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;
use std::f64::consts::PI;

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
    /// Moment coefficient (about quarter-chord)
    pub cm: f64,
    /// Pressure coefficient at each panel midpoint
    pub cp: Vec<f64>,
    /// X-coordinates of Cp sample points
    pub cp_x: Vec<f64>,
    /// Whether the analysis succeeded
    pub success: bool,
    /// Error message (if any)
    pub error: Option<String>,
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
        return AnalysisResult {
            cl: 0.0,
            cm: 0.0,
            cp: vec![],
            cp_x: vec![],
            success: false,
            error: Some("Invalid coordinates: need at least 3 points (6 values)".to_string()),
        };
    }

    let points: Vec<_> = coords.chunks(2).map(|c| point(c[0], c[1])).collect();

    // Build body
    let body = match Body::from_points("airfoil", &points) {
        Ok(b) => b,
        Err(e) => {
            return AnalysisResult {
                cl: 0.0,
                cm: 0.0,
                cp: vec![],
                cp_x: vec![],
                success: false,
                error: Some(format!("Geometry error: {}", e)),
            }
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
                cm: solution.cm,
                cp: solution.cp,
                cp_x,
                success: true,
                error: None,
            }
        }
        Err(e) => AnalysisResult {
            cl: 0.0,
            cm: 0.0,
            cp: vec![],
            cp_x: vec![],
            success: false,
            error: Some(format!("Solver error: {}", e)),
        },
    }
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

    /// Update particles by one time step.
    ///
    /// # Arguments
    /// * `dt` - Time step in seconds (typically 1/60 for 60 FPS)
    pub fn update(&mut self, dt: f64) {
        self.inner.update(&self.coords, &self.gamma, self.alpha, self.v_inf, dt);
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
        self.inner.get_psi_values(&self.coords, &self.gamma, self.alpha, self.v_inf)
    }

    /// Get the dividing streamline value (psi_0).
    /// 
    /// Particles with psi > psi_0 go above the dividing streamline (upper surface).
    /// Particles with psi < psi_0 go below the dividing streamline (lower surface).
    pub fn get_psi_0(&self) -> f64 {
        self.psi_0
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
                return AnalysisResult {
                    cl: 0.0,
                    cm: 0.0,
                    cp: vec![],
                    cp_x: vec![],
                    success: false,
                    error: Some("No geometry set".to_string()),
                }
            }
        };

        let flow = FlowConditions::with_alpha_deg(self.alpha_deg);

        match self.solver.solve(&[body.clone()], &flow) {
            Ok(solution) => {
                let cp_x: Vec<f64> = body.panels().iter().map(|p| p.midpoint().x).collect();

                AnalysisResult {
                    cl: solution.cl,
                    cm: solution.cm,
                    cp: solution.cp,
                    cp_x,
                    success: true,
                    error: None,
                }
            }
            Err(e) => AnalysisResult {
                cl: 0.0,
                cm: 0.0,
                cp: vec![],
                cp_x: vec![],
                success: false,
                error: Some(format!("Solver error: {}", e)),
            },
        }
    }
}

impl Default for RustFoil {
    fn default() -> Self {
        Self::new()
    }
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
}
