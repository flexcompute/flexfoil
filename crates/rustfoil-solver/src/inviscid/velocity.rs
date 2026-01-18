//! Velocity field evaluation and streamline integration.
//!
//! Provides functions to evaluate the velocity field from a panel solution
//! and integrate streamlines using RK4.

use rustfoil_core::Point;
use std::f64::consts::PI;

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
}
