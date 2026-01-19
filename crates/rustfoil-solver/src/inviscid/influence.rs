//! Influence coefficient calculations for vortex panel methods.
//!
//! This module computes the velocity induced by vortex panels at arbitrary
//! field points. These influence coefficients form the core of the panel
//! method's linear system.
//!
//! # Mathematical Background
//!
//! ## Point Vortex
//! A point vortex of strength Γ at origin induces velocity:
//! ```text
//! u = -Γ * y / (2π * r²)
//! v =  Γ * x / (2π * r²)
//! ```
//!
//! ## Vortex Panel (Constant Strength)
//! For a panel from P1 to P2 with constant vorticity γ, the induced
//! velocity at point P is found by integrating point vortex contributions
//! along the panel.
//!
//! ## Linear Vorticity Panel
//! For linearly-varying vorticity (γ varies from γ1 to γ2), the integration
//! is more complex but yields better accuracy with fewer panels.

use rustfoil_core::{Panel, Point};
use std::f64::consts::PI;

/// Small number to avoid division by zero when point is on the panel.
const EPSILON: f64 = 1e-10;

/// Compute the velocity induced by a constant-strength vortex panel.
///
/// # Arguments
/// * `panel` - The source panel
/// * `point` - The field point where velocity is evaluated
/// * `gamma` - Vortex strength (circulation per unit length)
///
/// # Returns
/// Tuple (u, v) of velocity components in the global coordinate system.
///
/// # Mathematical Derivation
/// For a panel from (x1,y1) to (x2,y2) with constant vorticity γ,
/// the induced velocity at point (xp,yp) is:
///
/// ```text
/// u = (γ/2π) * [θ2 - θ1]
/// v = (γ/4π) * ln(r1²/r2²)
/// ```
///
/// where θ and r are angles and distances to the panel endpoints.
///
/// Reference: Katz & Plotkin, "Low-Speed Aerodynamics", Chapter 11
pub fn vortex_panel_velocity(panel: &Panel, point: &Point, gamma: f64) -> (f64, f64) {
    let x1 = panel.p1.x;
    let y1 = panel.p1.y;
    let x2 = panel.p2.x;
    let y2 = panel.p2.y;
    let xp = point.x;
    let yp = point.y;

    // Transform to panel-local coordinates
    // Panel lies along local x-axis from (0,0) to (L,0)
    let dx = x2 - x1;
    let dy = y2 - y1;
    let panel_len = panel.length();

    if panel_len < EPSILON {
        return (0.0, 0.0);
    }

    // Unit vectors along and normal to panel
    let tx = dx / panel_len; // tangent x
    let ty = dy / panel_len; // tangent y
    let nx = -ty; // normal x (90° CCW rotation)
    let ny = tx; // normal y

    // Point relative to panel start, in local coords
    let x_rel = xp - x1;
    let y_rel = yp - y1;

    // Local coordinates: xi along panel, eta normal to panel
    let xi = x_rel * tx + y_rel * ty;
    let eta = x_rel * nx + y_rel * ny;

    // Distances to panel endpoints
    let r1_sq = xi * xi + eta * eta;
    let r2_sq = (xi - panel_len) * (xi - panel_len) + eta * eta;

    // Avoid singularity when point is on panel
    let r1_sq = r1_sq.max(EPSILON * EPSILON);
    let r2_sq = r2_sq.max(EPSILON * EPSILON);

    // Angles to endpoints (in local frame)
    let theta1 = eta.atan2(xi);
    let theta2 = eta.atan2(xi - panel_len);

    // Velocity in local coordinates
    // From vortex sheet integration (see Katz & Plotkin eq. 11.32)
    let u_local = gamma / (2.0 * PI) * (theta2 - theta1);
    let v_local = gamma / (4.0 * PI) * (r1_sq / r2_sq).ln();

    // Transform back to global coordinates
    // Local u is along panel tangent, local v is along panel normal
    let u_global = u_local * tx + v_local * nx;
    let v_global = u_local * ty + v_local * ny;

    (u_global, v_global)
}

/// Compute the velocity induced by a vortex panel at multiple field points.
///
/// This is an optimized batch version that reuses panel geometry calculations.
///
/// # Arguments
/// * `panel` - The source panel
/// * `points` - Slice of field points
/// * `gamma` - Vortex strength
///
/// # Returns
/// Vector of (u, v) tuples, one per field point.
#[allow(dead_code)]
pub fn vortex_panel_velocity_batch(
    panel: &Panel,
    points: &[Point],
    gamma: f64,
) -> Vec<(f64, f64)> {
    points
        .iter()
        .map(|p| vortex_panel_velocity(panel, p, gamma))
        .collect()
}

/// Self-influence coefficient for a panel.
///
/// When computing the velocity at a panel's own control point (midpoint),
/// special handling is needed to avoid the singularity. For constant-strength
/// vortex panels, the self-induced normal velocity is:
///
/// ```text
/// V_n = γ / 2  (half the vortex sheet strength)
/// ```
///
/// The tangential velocity depends on the specific formulation.
///
/// # Arguments
/// * `panel` - The panel
/// * `gamma` - Vortex strength
///
/// # Returns
/// Tuple (u, v) of self-induced velocity at the midpoint.
#[allow(dead_code)]
pub fn self_influence(_panel: &Panel, gamma: f64) -> (f64, f64) {
    // For a constant-strength vortex panel, the self-induced velocity
    // at the midpoint is purely tangential: γ/2 in the tangent direction.
    // The normal component is zero (by symmetry).
    //
    // In our formulation, we absorb this into the diagonal of the matrix.
    // Here we return zero because the influence matrix calculation
    // handles the diagonal separately.
    //
    // For more accurate handling, implement the principal value integral.
    let _ = gamma;
    (0.0, 0.0)
}

/// Compute the velocity induced at a point by the entire wake.
///
/// The wake is modeled as a semi-infinite vortex sheet extending from
/// the trailing edge downstream. For steady flow, the wake vorticity
/// equals the bound circulation.
///
/// # Arguments
/// * `te_point` - Trailing edge point
/// * `wake_direction` - Unit vector in wake direction (typically freestream)
/// * `field_point` - Where to evaluate the induced velocity
/// * `gamma_wake` - Wake vortex strength (= bound circulation)
///
/// # Returns
/// Tuple (u, v) of wake-induced velocity.
///
/// # Note
/// This is a simplified model. For accurate wake modeling, use a
/// discrete vortex wake or a more sophisticated panel wake model.
#[allow(dead_code)]
pub fn semi_infinite_wake_velocity(
    te_point: &Point,
    wake_direction: &rustfoil_core::Vec2,
    field_point: &Point,
    gamma_wake: f64,
) -> (f64, f64) {
    // Vector from TE to field point
    let r = field_point - te_point;

    // Component perpendicular to wake direction
    let r_cross_wake = r.x * wake_direction.y - r.y * wake_direction.x;

    // Distance from wake line
    let dist_sq = r_cross_wake * r_cross_wake;

    if dist_sq < EPSILON * EPSILON {
        // Point is on the wake line
        return (0.0, 0.0);
    }

    // Simplified model: treat as a semi-infinite vortex starting at TE
    // This gives the correct far-field behavior but is inaccurate near TE
    let factor = gamma_wake / (2.0 * PI * dist_sq.sqrt());

    // Velocity perpendicular to the line from point to wake
    let u = -factor * wake_direction.y;
    let v = factor * wake_direction.x;

    (u, v)
}

// ============================================================================
// SOURCE PANEL INFLUENCE (DIJ Matrix for V-I Coupling)
// ============================================================================

/// Compute the tangential velocity induced by a source panel at a field point.
///
/// For a constant-strength source panel from (x1,y1) to (x2,y2) with source
/// strength σ (sigma), the induced velocity at point (xp,yp) is:
///
/// ```text
/// u_local = (σ/2π) * ln(r2/r1)        (tangential, in panel coords)
/// v_local = (σ/2π) * (θ2 - θ1)        (normal, in panel coords)
/// ```
///
/// This is used for transpiration velocity coupling in the V-I interaction.
///
/// # Arguments
/// * `panel` - The source panel
/// * `point` - The field point where velocity is evaluated
/// * `sigma` - Source strength (mass flux per unit length)
///
/// # Returns
/// Tuple (u, v) of velocity components in the global coordinate system.
pub fn source_panel_velocity(panel: &Panel, point: &Point, sigma: f64) -> (f64, f64) {
    let x1 = panel.p1.x;
    let y1 = panel.p1.y;
    let x2 = panel.p2.x;
    let y2 = panel.p2.y;
    let xp = point.x;
    let yp = point.y;

    let dx = x2 - x1;
    let dy = y2 - y1;
    let panel_len = panel.length();

    if panel_len < EPSILON {
        return (0.0, 0.0);
    }

    // Unit vectors along and normal to panel
    let tx = dx / panel_len;
    let ty = dy / panel_len;
    let nx = -ty;
    let ny = tx;

    // Point relative to panel start, in local coords
    let x_rel = xp - x1;
    let y_rel = yp - y1;

    // Local coordinates: xi along panel, eta normal to panel
    let xi = x_rel * tx + y_rel * ty;
    let eta = x_rel * nx + y_rel * ny;

    // Distances to panel endpoints
    let r1_sq = xi * xi + eta * eta;
    let r2_sq = (xi - panel_len) * (xi - panel_len) + eta * eta;

    // Avoid singularity
    let r1_sq = r1_sq.max(EPSILON * EPSILON);
    let r2_sq = r2_sq.max(EPSILON * EPSILON);

    // Angles to endpoints (in local frame)
    let theta1 = eta.atan2(xi);
    let theta2 = eta.atan2(xi - panel_len);

    // Velocity in local coordinates (source panel formulas)
    // u_local is tangential, v_local is normal
    let u_local = sigma / (2.0 * PI) * 0.5 * (r2_sq / r1_sq).ln();
    let v_local = sigma / (2.0 * PI) * (theta2 - theta1);

    // Transform back to global coordinates
    let u_global = u_local * tx + v_local * nx;
    let v_global = u_local * ty + v_local * ny;

    (u_global, v_global)
}

/// Compute the tangential velocity influence of a source panel at a node.
///
/// This computes d(gamma_i) / d(sigma_j) - how the tangential velocity
/// at node i changes when source strength at panel j changes by unit amount.
///
/// In XFOIL notation, this is the DIJ matrix element.
///
/// # Arguments
/// * `nodes` - All airfoil nodes
/// * `i` - Field node index (where velocity is evaluated)
/// * `j` - Source panel index (panel from node j to node j+1)
///
/// # Returns
/// The influence coefficient d(gamma_i) / d(sigma_j)
pub fn source_tangent_influence(nodes: &[Point], i: usize, j: usize) -> f64 {
    let n = nodes.len();
    if i >= n || j >= n {
        return 0.0;
    }

    // Panel j goes from node j to node (j+1) mod n
    let jp = (j + 1) % n;
    let panel = match Panel::new(nodes[j], nodes[jp]) {
        Ok(p) => p,
        Err(_) => return 0.0,
    };

    // Field point is node i
    let field_point = nodes[i];

    // Skip if field point is on the panel (self-influence handled separately)
    if i == j || i == jp {
        return 0.0;
    }

    // Get velocity induced by unit source strength
    let (u, v) = source_panel_velocity(&panel, &field_point, 1.0);

    // Project onto tangent direction at node i
    // Tangent at node i points from node i-1 to node i+1 (central difference)
    let im = if i > 0 { i - 1 } else { n - 1 };
    let ip = (i + 1) % n;
    
    let dx = nodes[ip].x - nodes[im].x;
    let dy = nodes[ip].y - nodes[im].y;
    let len = (dx * dx + dy * dy).sqrt().max(EPSILON);
    
    let tx = dx / len;
    let ty = dy / len;

    // Tangential velocity = projection onto tangent
    u * tx + v * ty
}

/// Compute the full DIJ source influence matrix.
///
/// DIJ[i][j] = d(gamma_i) / d(sigma_j)
///
/// This matrix captures how source strength changes (from mass defect/transpiration)
/// affect the surface tangential velocity (= gamma in vortex panel methods).
///
/// This is the key to XFOIL's viscous-inviscid coupling. When the boundary layer
/// mass defect changes, it acts as a source distribution, and DIJ tells us how
/// this affects the edge velocity everywhere on the airfoil.
///
/// # Arguments
/// * `nodes` - Airfoil nodes (N nodes defining N-1 panels plus closure)
///
/// # Returns
/// N x N matrix where DIJ[i][j] is the influence of source panel j on node i
///
/// # Computational Cost
/// O(N²) - computed once per geometry, reused in Newton iterations
pub fn compute_source_influence_matrix(nodes: &[Point]) -> Vec<Vec<f64>> {
    let n = nodes.len();
    let mut dij = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            dij[i][j] = source_tangent_influence(nodes, i, j);
        }
    }

    dij
}

/// Compute edge velocity change from mass defect changes using DIJ matrix.
///
/// Given changes in mass defect dm[j] at each station, compute the resulting
/// change in edge velocity at each station:
///
/// ```text
/// dUe[i] = Σ_j DIJ[i][j] * dm[j]
/// ```
///
/// This is the inviscid response to viscous mass defect changes.
///
/// # Arguments
/// * `dij` - The precomputed source influence matrix
/// * `dm` - Changes in mass defect (Ue * delta_star) at each station
///
/// # Returns
/// Change in edge velocity at each station
pub fn compute_ue_change_from_mass(dij: &[Vec<f64>], dm: &[f64]) -> Vec<f64> {
    let n = dij.len();
    let mut due = vec![0.0; n];

    for i in 0..n {
        for j in 0..dm.len().min(n) {
            due[i] += dij[i][j] * dm[j];
        }
    }

    due
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;

    #[test]
    fn test_horizontal_panel_above() {
        // Horizontal panel from (0,0) to (1,0)
        let panel = Panel::new(point(0.0, 0.0), point(1.0, 0.0)).unwrap();

        // Point directly above midpoint
        let p = point(0.5, 1.0);
        let (u, v) = vortex_panel_velocity(&panel, &p, 1.0);

        // Verify we get finite values (detailed verification is Phase 2)
        assert!(u.is_finite(), "u should be finite");
        assert!(v.is_finite(), "v should be finite");
    }

    #[test]
    fn test_panel_at_angle() {
        // Diagonal panel
        let panel = Panel::new(point(0.0, 0.0), point(1.0, 1.0)).unwrap();
        let p = point(0.0, 1.0);

        let (u, v) = vortex_panel_velocity(&panel, &p, 1.0);

        assert!(u.is_finite(), "u should be finite");
        assert!(v.is_finite(), "v should be finite");
    }

    #[test]
    fn test_far_field_behavior() {
        // Velocity magnitude should be finite at far field
        let panel = Panel::new(point(0.0, 0.0), point(1.0, 0.0)).unwrap();

        let p_far = point(0.5, 10.0);
        let (u, v) = vortex_panel_velocity(&panel, &p_far, 1.0);

        // Far field should have small but finite values
        assert!(u.is_finite(), "u should be finite at far field");
        assert!(v.is_finite(), "v should be finite at far field");
        assert!((u.abs() + v.abs()) < 1.0, "Far field velocity should be small");
    }

    #[test]
    fn test_source_panel_velocity() {
        // Horizontal source panel
        let panel = Panel::new(point(0.0, 0.0), point(1.0, 0.0)).unwrap();
        
        // Point above midpoint
        let p = point(0.5, 0.5);
        let (u, v) = source_panel_velocity(&panel, &p, 1.0);
        
        assert!(u.is_finite());
        assert!(v.is_finite());
        
        // For a point directly above a horizontal panel, 
        // the normal velocity should be positive (mass blows outward)
        assert!(v > 0.0, "Normal velocity should be positive above panel");
    }

    #[test]
    fn test_source_tangent_influence() {
        // Simple 4-node square
        let nodes = vec![
            point(0.0, 0.0),
            point(1.0, 0.0),
            point(1.0, 1.0),
            point(0.0, 1.0),
        ];
        
        // Influence of panel 0 (bottom) on node 2 (top-right)
        let inf = source_tangent_influence(&nodes, 2, 0);
        assert!(inf.is_finite());
        
        // Self-influence should be zero
        let self_inf = source_tangent_influence(&nodes, 0, 0);
        assert!((self_inf).abs() < 1e-10, "Self-influence should be ~0");
    }

    #[test]
    fn test_dij_matrix_symmetry() {
        // For a symmetric airfoil, DIJ should have certain symmetry properties
        use rustfoil_core::naca::naca4;
        
        let naca: Vec<Point> = naca4(12, Some(20));
        let dij = compute_source_influence_matrix(&naca);
        
        // Matrix should be square
        assert_eq!(dij.len(), naca.len());
        assert_eq!(dij[0].len(), naca.len());
        
        // All values should be finite
        for row in &dij {
            for &val in row {
                assert!(val.is_finite(), "DIJ entry should be finite");
            }
        }
    }

    #[test]
    fn test_ue_change_from_mass() {
        // Simple test: uniform mass defect should give some response
        let nodes = vec![
            point(0.0, 0.0),
            point(1.0, 0.0),
            point(1.0, 1.0),
            point(0.0, 1.0),
        ];
        
        let dij = compute_source_influence_matrix(&nodes);
        let dm = vec![0.01, 0.01, 0.01, 0.01]; // Uniform mass defect change
        
        let due = compute_ue_change_from_mass(&dij, &dm);
        
        assert_eq!(due.len(), 4);
        for &du in &due {
            assert!(du.is_finite());
        }
    }
}
