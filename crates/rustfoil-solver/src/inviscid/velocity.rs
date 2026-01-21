//! Velocity field evaluation, stream function computation, and streamline integration.
//!
//! Provides functions to evaluate the velocity field and stream function from a panel solution
//! and integrate streamlines using RK4.

use rustfoil_core::Point;
use std::f64::consts::PI;

/// 1/(4π) - used for stream function influence coefficients (matches XFOIL's QOPI)
const QOPI: f64 = 0.25 / PI;

/// Evaluate velocity at a point (x, y) given the panel solution.
///
/// Uses K&P VOR2DL linear vorticity panel method (Eq 11.99-11.100).
/// Gamma values are at nodes, varying linearly across each panel.
pub fn velocity_at(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
) -> (f64, f64) {
    let n = nodes.len();
    if n < 2 || gamma.len() != n {
        return (v_inf * alpha.cos(), v_inf * alpha.sin());
    }

    // Freestream velocity
    let mut u = v_inf * alpha.cos();
    let mut v = v_inf * alpha.sin();

    let two_pi = 2.0 * PI;

    // Add contribution from each panel
    // Panel j goes from node j to node (j+1) % n
    for j in 0..n {
        let jp = (j + 1) % n;

        let x1 = nodes[j].x;
        let y1 = nodes[j].y;
        let x2 = nodes[jp].x;
        let y2 = nodes[jp].y;

        let dx = x2 - x1;
        let dy = y2 - y1;
        let panel_len = (dx * dx + dy * dy).sqrt();

        if panel_len < 1e-12 {
            continue;
        }

        // Panel angle and trig values
        let theta = dy.atan2(dx);
        let cos_t = theta.cos();
        let sin_t = theta.sin();

        // Transform field point to panel-local coordinates
        // xp: along panel (0 at start, panel_len at end)
        // yp: perpendicular (positive to LEFT of panel direction)
        let xp = (x - x1) * cos_t + (y - y1) * sin_t;
        let yp = -(x - x1) * sin_t + (y - y1) * cos_t;
        let xp2 = xp - panel_len; // xp relative to panel end

        // Squared distances from field point to panel endpoints
        let r1_sq = xp * xp + yp * yp;
        let r2_sq = xp2 * xp2 + yp * yp;

        // Skip if too close to panel endpoints (singularity)
        if r1_sq < 1e-10 || r2_sq < 1e-10 {
            continue;
        }

        // Angles from field point to panel endpoints
        let theta1 = yp.atan2(xp);
        let theta2 = yp.atan2(xp2);
        let beta = theta2 - theta1;

        // Log term: ln(r₂/r₁) = 0.5 * ln(r₂²/r₁²)
        let logterm = 0.5 * (r2_sq / r1_sq).ln();

        let inv_2pi_l = 1.0 / (two_pi * panel_len);

        // K&P VOR2DL (Eq 11.99-11.100): influence from γⱼ = 1, γⱼ₊₁ = 0
        let u1_local = -(yp * logterm + xp * beta - panel_len * beta) * inv_2pi_l;
        let w1_local = -((panel_len - yp * beta) + xp * logterm - panel_len * logterm) * inv_2pi_l;

        // K&P VOR2DL (Eq 11.99-11.100): influence from γⱼ = 0, γⱼ₊₁ = 1
        let u2_local = (yp * logterm + xp * beta) * inv_2pi_l;
        let w2_local = ((panel_len - yp * beta) + xp * logterm) * inv_2pi_l;

        // Transform velocities back to global coordinates
        let u1 = u1_local * cos_t - w1_local * sin_t;
        let v1 = u1_local * sin_t + w1_local * cos_t;
        let u2 = u2_local * cos_t - w2_local * sin_t;
        let v2 = u2_local * sin_t + w2_local * cos_t;

        // Add weighted by gamma at each node
        u += gamma[j] * u1 + gamma[jp] * u2;
        v += gamma[j] * v1 + gamma[jp] * v2;
    }

    (u, v)
}

/// Evaluate the stream function at a point (x, y) given the panel solution.
///
/// Uses XFOIL's PSILIN formulation with linear vorticity panels.
/// The stream function ψ is constant along streamlines; the body surface is at ψ = ψ₀.
///
/// # Arguments
/// * `x`, `y` - Field point coordinates
/// * `nodes` - Airfoil node positions (closed contour)
/// * `gamma` - Vorticity at each node
/// * `alpha` - Angle of attack (radians)
/// * `v_inf` - Freestream velocity magnitude
///
/// # Returns
/// Stream function value at (x, y). Returns NaN if inside the airfoil.
pub fn psi_at(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
) -> f64 {
    let n = nodes.len();
    if n < 2 || gamma.len() != n {
        // Just freestream
        return v_inf * (alpha.cos() * y - alpha.sin() * x);
    }

    // Return NaN for points inside the airfoil
    if is_inside_airfoil(x, y, nodes) {
        return f64::NAN;
    }

    let mut psi = 0.0;

    // Add contribution from each vortex panel
    // Panel j goes from node j to node (j+1) % n
    for jo in 0..n {
        let jp = (jo + 1) % n;

        // Panel endpoints
        let x_jo = nodes[jo].x;
        let y_jo = nodes[jo].y;
        let x_jp = nodes[jp].x;
        let y_jp = nodes[jp].y;

        // Panel geometry
        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        if ds_sq < 1e-24 {
            continue;
        }

        let ds = ds_sq.sqrt();
        let sx = dx / ds; // Panel tangent x
        let sy = dy / ds; // Panel tangent y

        // Vector from panel endpoints to field point
        let rx1 = x - x_jo;
        let ry1 = y - y_jo;

        // Transform to panel-local coordinates (XFOIL convention)
        // X1: tangential distance from field point to node JO
        // X2: tangential distance from field point to node JP
        // YY: normal distance from field point to panel
        let x1 = sx * rx1 + sy * ry1;
        let x2 = x1 - ds; // = sx*(x-x_jp) + sy*(y-y_jp)
        let yy = sx * ry1 - sy * rx1;

        // Squared distances to panel endpoints
        let rs1 = x1 * x1 + yy * yy;
        let rs2 = x2 * x2 + yy * yy;

        // XFOIL's reflection correction for atan2 branch cuts
        // When YY < 0, reflect arguments and add π offset
        let sgn = if yy >= 0.0 { 1.0 } else { -1.0 };
        let pi_offset = (0.5 - 0.5 * sgn) * PI;

        // Log and arctangent terms, with singularity handling
        let g1 = if rs1 > 1e-20 { rs1.ln() } else { 0.0 };
        let g2 = if rs2 > 1e-20 { rs2.ln() } else { 0.0 };
        let t1 = (sgn * x1).atan2(sgn * yy) + pi_offset;
        let t2 = (sgn * x2).atan2(sgn * yy) + pi_offset;

        // XFOIL's PSIS/PSID formulation (xpanel.f lines 350-351)
        // PSIS = coefficient for (γ_JP + γ_JO) / 2
        // PSID = coefficient for (γ_JP - γ_JO) / 2
        let dxinv = 1.0 / (x1 - x2);
        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

        // Sum and difference of gamma at endpoints
        let gsum = gamma[jp] + gamma[jo];
        let gdif = gamma[jp] - gamma[jo];

        // Accumulate stream function contribution
        psi += QOPI * (psis * gsum + psid * gdif);
    }

    // Freestream contribution: ψ∞ = V∞(cos(α)·y - sin(α)·x)
    psi += v_inf * (alpha.cos() * y - alpha.sin() * x);

    psi
}

/// Compute stream function values on a rectangular grid.
///
/// # Arguments
/// * `nodes` - Airfoil node positions
/// * `gamma` - Vorticity at each node
/// * `alpha` - Angle of attack (radians)
/// * `v_inf` - Freestream velocity magnitude
/// * `x_min`, `x_max`, `y_min`, `y_max` - Grid bounds
/// * `nx`, `ny` - Grid resolution
///
/// # Returns
/// Row-major array of stream function values: psi[iy * nx + ix].
/// Points inside the airfoil have value NaN.
pub fn compute_psi_grid(
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    nx: usize,
    ny: usize,
) -> Vec<f64> {
    compute_psi_grid_with_interior(nodes, gamma, alpha, v_inf, x_min, x_max, y_min, y_max, nx, ny, None)
}

/// Compute stream function values on a rectangular grid, with optional interior value.
///
/// # Arguments
/// * `nodes` - Airfoil node positions
/// * `gamma` - Vorticity at each node
/// * `alpha` - Angle of attack (radians)
/// * `v_inf` - Freestream velocity magnitude
/// * `x_min`, `x_max`, `y_min`, `y_max` - Grid bounds
/// * `nx`, `ny` - Grid resolution
/// * `interior_value` - Value to use inside the airfoil (typically ψ₀). If None, uses NaN.
///
/// # Returns
/// Row-major array of stream function values: psi[iy * nx + ix].
pub fn compute_psi_grid_with_interior(
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    nx: usize,
    ny: usize,
    interior_value: Option<f64>,
) -> Vec<f64> {
    let mut grid = vec![0.0; nx * ny];

    let dx = if nx > 1 {
        (x_max - x_min) / (nx - 1) as f64
    } else {
        0.0
    };
    let dy = if ny > 1 {
        (y_max - y_min) / (ny - 1) as f64
    } else {
        0.0
    };

    for iy in 0..ny {
        let y = y_min + iy as f64 * dy;
        for ix in 0..nx {
            let x = x_min + ix as f64 * dx;
            let psi = psi_at(x, y, nodes, gamma, alpha, v_inf);
            // If inside airfoil (NaN) and we have an interior value, use it
            grid[iy * nx + ix] = if psi.is_nan() {
                interior_value.unwrap_or(f64::NAN)
            } else {
                psi
            };
        }
    }

    grid
}

/// Check if a point is inside the airfoil.
/// Uses ray casting algorithm.
pub fn is_inside_airfoil(x: f64, y: f64, nodes: &[Point]) -> bool {
    let n = nodes.len();
    if n < 3 {
        return false;
    }
    
    let mut inside = false;
    let mut j = n - 1;
    
    for i in 0..n {
        let xi = nodes[i].x;
        let yi = nodes[i].y;
        let xj = nodes[j].x;
        let yj = nodes[j].y;
        
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    
    inside
}

/// RK4 integration step for streamline tracing.
fn rk4_step<F>(field: &F, x: f64, y: f64, dt: f64) -> Option<(f64, f64)>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    let (k1x, k1y) = field(x, y);
    let speed1 = (k1x * k1x + k1y * k1y).sqrt();
    
    if !speed1.is_finite() || speed1 < 1e-8 {
        return None;
    }
    
    let (k2x, k2y) = field(x + 0.5 * dt * k1x, y + 0.5 * dt * k1y);
    let (k3x, k3y) = field(x + 0.5 * dt * k2x, y + 0.5 * dt * k2y);
    let (k4x, k4y) = field(x + dt * k3x, y + dt * k3y);
    
    let new_x = x + (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    let new_y = y + (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
    
    Some((new_x, new_y))
}

/// RK4 step with adaptive arc-length control.
fn rk4_arc_step<F>(
    field: &F,
    x: f64,
    y: f64,
    ds: f64,
    dt_min: f64,
    dt_max: f64,
) -> Option<(f64, f64)>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    let (u, v) = field(x, y);
    let speed = (u * u + v * v).sqrt();
    
    if !speed.is_finite() || speed < 1e-8 {
        return None;
    }
    
    // Choose dt so distance traveled ≈ ds
    let dt = (ds / speed).clamp(dt_min, dt_max);
    rk4_step(field, x, y, dt)
}

/// Integrate a single streamline from a seed point.
fn integrate_streamline<F>(
    field: &F,
    x0: f64,
    y0: f64,
    step_size: f64,
    max_steps: usize,
    nodes: &[Point],
    bounds: (f64, f64, f64, f64),
) -> Vec<(f64, f64)>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    let (x_min, x_max, y_min, y_max) = bounds;
    let mut points = vec![(x0, y0)];
    let mut x = x0;
    let mut y = y0;
    
    for _ in 0..max_steps {
        // Check if inside airfoil
        if is_inside_airfoil(x, y, nodes) {
            break;
        }
        
        // Check bounds
        if x < x_min || x > x_max || y < y_min || y > y_max {
            break;
        }
        
        match rk4_arc_step(field, x, y, step_size, 1e-4, 0.02) {
            Some((new_x, new_y)) => {
                x = new_x;
                y = new_y;
                points.push((x, y));
            }
            None => break,
        }
    }
    
    points
}

/// Streamline generation options.
#[derive(Debug, Clone)]
pub struct StreamlineOptions {
    pub seed_count: usize,
    pub seed_x: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub step_size: f64,
    pub max_steps: usize,
    pub x_min: f64,
    pub x_max: f64,
}

impl Default for StreamlineOptions {
    fn default() -> Self {
        Self {
            seed_count: 25,
            seed_x: -0.5,
            y_min: -0.4,
            y_max: 0.4,
            step_size: 0.01,
            max_steps: 2000,
            x_min: -1.0,
            x_max: 2.0,
        }
    }
}

/// Build streamlines from seed points.
pub fn build_streamlines(
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    options: &StreamlineOptions,
) -> Vec<Vec<(f64, f64)>> {
    let mut streamlines = Vec::with_capacity(options.seed_count);
    
    // Create velocity field closure
    let field = |x: f64, y: f64| velocity_at(x, y, nodes, gamma, alpha, v_inf);
    
    let bounds = (options.x_min, options.x_max, options.y_min, options.y_max);
    
    // Generate seed points
    for i in 0..options.seed_count {
        let t = i as f64 / (options.seed_count - 1).max(1) as f64;
        let y = options.y_min + t * (options.y_max - options.y_min);
        let x = options.seed_x;
        
        // Skip seeds that start inside the airfoil
        if is_inside_airfoil(x, y, nodes) {
            continue;
        }
        
        let streamline = integrate_streamline(
            &field,
            x,
            y,
            options.step_size,
            options.max_steps,
            nodes,
            bounds,
        );
        
        if streamline.len() >= 2 {
            streamlines.push(streamline);
        }
    }
    
    streamlines
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;
    
    fn make_circle(n: usize, radius: f64) -> Vec<Point> {
        (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                point(radius * theta.cos(), radius * theta.sin())
            })
            .collect()
    }
    
    #[test]
    fn test_is_inside_airfoil() {
        let circle = make_circle(32, 0.5);
        
        assert!(is_inside_airfoil(0.0, 0.0, &circle));
        assert!(is_inside_airfoil(0.2, 0.1, &circle));
        assert!(!is_inside_airfoil(1.0, 0.0, &circle));
        assert!(!is_inside_airfoil(0.0, 1.0, &circle));
    }
    
    #[test]
    fn test_velocity_freestream() {
        let nodes: Vec<Point> = vec![];
        let gamma: Vec<f64> = vec![];
        
        let (u, v) = velocity_at(0.0, 0.0, &nodes, &gamma, 0.0, 1.0);
        assert!((u - 1.0).abs() < 1e-10);
        assert!(v.abs() < 1e-10);
        
        let alpha = 5.0_f64.to_radians();
        let (u, v) = velocity_at(0.0, 0.0, &nodes, &gamma, alpha, 1.0);
        assert!((u - alpha.cos()).abs() < 1e-10);
        assert!((v - alpha.sin()).abs() < 1e-10);
    }
    
    #[test]
    fn test_psi_freestream() {
        // Test freestream-only stream function
        let nodes: Vec<Point> = vec![];
        let gamma: Vec<f64> = vec![];
        
        // For α=0, ψ = V∞ * y
        let psi = psi_at(0.5, 0.3, &nodes, &gamma, 0.0, 1.0);
        assert!((psi - 0.3).abs() < 1e-10, "ψ = y for freestream at α=0");
        
        // For α=90°, ψ = -V∞ * x
        let alpha = std::f64::consts::FRAC_PI_2;
        let psi = psi_at(0.5, 0.3, &nodes, &gamma, alpha, 1.0);
        assert!((psi - (-0.5)).abs() < 1e-10, "ψ = -x for freestream at α=90°");
    }
    
    #[test]
    fn test_psi_inside_airfoil() {
        let circle = make_circle(32, 0.5);
        let gamma = vec![0.0; 32];
        
        // Inside should return NaN
        let psi = psi_at(0.0, 0.0, &circle, &gamma, 0.0, 1.0);
        assert!(psi.is_nan(), "ψ should be NaN inside airfoil");
        
        // Outside should return finite value
        let psi = psi_at(1.0, 0.0, &circle, &gamma, 0.0, 1.0);
        assert!(psi.is_finite(), "ψ should be finite outside airfoil");
    }
    
    #[test]
    fn test_psi_grid() {
        let circle = make_circle(32, 0.3);
        let gamma = vec![0.0; 32];
        
        let grid = compute_psi_grid(
            &circle, &gamma, 0.0, 1.0,
            -1.0, 2.0, -1.0, 1.0,
            10, 8
        );
        
        assert_eq!(grid.len(), 80);  // 10 * 8
        
        // Check that some interior points are NaN
        let nan_count = grid.iter().filter(|&&v| v.is_nan()).count();
        assert!(nan_count > 0, "Should have some NaN values inside airfoil");
        
        // Check that most exterior points are finite
        let finite_count = grid.iter().filter(|&&v| v.is_finite()).count();
        assert!(finite_count > 50, "Should have many finite values outside airfoil");
    }
}
