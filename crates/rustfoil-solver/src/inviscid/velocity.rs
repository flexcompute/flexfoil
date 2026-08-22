//! Velocity field evaluation, stream function computation, and streamline integration.
//!
//! Provides functions to evaluate the velocity field and stream function from a panel solution
//! and integrate streamlines using RK4.

use rustfoil_core::{Layout, Point};
use std::f64::consts::PI;

/// 1/(4π) - used for stream function influence coefficients (matches XFOIL's QOPI)
const QOPI: f64 = 0.25 / PI;

/// Connectivity for a node array that holds one closed contour.
///
/// Every panel loop in this module closes the panel starting at node `i` onto
/// [`Layout::next_node(i)`](Layout::next_node), which wraps inside the owning
/// element: an element's last node connects back to its own first node and
/// never to the next element's. With a single element that is exactly
/// `(i + 1) % n`, the closure the single-airfoil entry points have always used,
/// so those results are unchanged.
fn single_element_layout(n: usize) -> Layout {
    if n == 0 {
        // No nodes, so no element to describe.
        return Layout::from_node_counts(&[]).expect("an empty layout has no element to reject");
    }
    // `from_node_counts` rejects only an element with no nodes.
    Layout::from_node_counts(&[n]).expect("one element of at least one node")
}

/// Wake panel geometry and source strengths for viscous streamline computation.
///
/// Wake panels are an open polyline extending from the trailing edge downstream.
/// Each node carries a source strength (mass defect = Ue * delta_star) that
/// displaces streamlines outward to account for the viscous wake thickness.
#[derive(Debug, Clone)]
pub struct WakePanels {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub sigma: Vec<f64>,
}

/// Evaluate velocity at a point (x, y) given the panel solution.
///
/// Uses K&P VOR2DL linear vorticity panel method (Eq 11.99-11.100).
/// Gamma values are at nodes, varying linearly across each panel.
///
/// `nodes` is one closed contour. Use [`velocity_at_multi`] for a node array
/// that concatenates several elements.
pub fn velocity_at(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
) -> (f64, f64) {
    let layout = single_element_layout(nodes.len());
    velocity_at_multi(x, y, nodes, gamma, alpha, v_inf, &layout)
}

/// Evaluate velocity at a point (x, y) for a node array of several elements.
///
/// As [`velocity_at`], but `layout` supplies the element connectivity, so each
/// element's trailing-edge panel closes onto that element's own first node
/// instead of onto the next element's — no panel is created across the gap
/// between two elements.
///
/// `layout` must describe `nodes`: if [`Layout::total_nodes`] disagrees with
/// `nodes.len()`, only the freestream is returned.
pub fn velocity_at_multi(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    layout: &Layout,
) -> (f64, f64) {
    let n = nodes.len();
    if n < 2 || gamma.len() != n || layout.total_nodes() != n {
        return (v_inf * alpha.cos(), v_inf * alpha.sin());
    }

    // Freestream velocity
    let mut u = v_inf * alpha.cos();
    let mut v = v_inf * alpha.sin();

    let two_pi = 2.0 * PI;

    // Add contribution from each panel
    // Panel j goes from node j to the next node in the same element
    for j in 0..n {
        let jp = layout.next_node(j);

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

/// Evaluate velocity including source panel contributions from viscous coupling.
///
/// Extends `velocity_at` with:
/// - `sigma`: source strength (mass defect = Ue * δ*) at each airfoil node.
///   The source velocity influence is the 90° rotation of the vortex influence.
/// - `wake_panels`: optional wake geometry with source strengths. Wake panels
///   are an open polyline (not closed), so panel j goes from node j to j+1.
///
/// `nodes` is one closed contour. Use [`velocity_at_with_sources_multi`] for a
/// node array that concatenates several elements.
pub fn velocity_at_with_sources(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    wake_panels: Option<&WakePanels>,
) -> (f64, f64) {
    let layout = single_element_layout(nodes.len());
    velocity_at_with_sources_multi(
        x,
        y,
        nodes,
        gamma,
        sigma,
        alpha,
        v_inf,
        wake_panels,
        &layout,
    )
}

/// Evaluate velocity including source panels, for a node array of several
/// elements.
///
/// As [`velocity_at_with_sources`], but `layout` supplies the element
/// connectivity, so no panel is created across the gap between two elements.
///
/// `layout` must describe `nodes`: if [`Layout::total_nodes`] disagrees with
/// `nodes.len()`, only the freestream is returned.
#[allow(clippy::too_many_arguments)]
pub fn velocity_at_with_sources_multi(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    wake_panels: Option<&WakePanels>,
    layout: &Layout,
) -> (f64, f64) {
    let n = nodes.len();
    if n < 2 || gamma.len() != n || sigma.len() != n || layout.total_nodes() != n {
        return (v_inf * alpha.cos(), v_inf * alpha.sin());
    }

    let mut u = v_inf * alpha.cos();
    let mut v = v_inf * alpha.sin();

    let two_pi = 2.0 * PI;

    for j in 0..n {
        let jp = layout.next_node(j);

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

        let theta = dy.atan2(dx);
        let cos_t = theta.cos();
        let sin_t = theta.sin();

        let xp = (x - x1) * cos_t + (y - y1) * sin_t;
        let yp = -(x - x1) * sin_t + (y - y1) * cos_t;
        let xp2 = xp - panel_len;

        let r1_sq = xp * xp + yp * yp;
        let r2_sq = xp2 * xp2 + yp * yp;

        if r1_sq < 1e-10 || r2_sq < 1e-10 {
            continue;
        }

        let theta1 = yp.atan2(xp);
        let theta2 = yp.atan2(xp2);
        let beta = theta2 - theta1;
        let logterm = 0.5 * (r2_sq / r1_sq).ln();
        let inv_2pi_l = 1.0 / (two_pi * panel_len);

        // Vortex influence (K&P VOR2DL)
        let u1_local = -(yp * logterm + xp * beta - panel_len * beta) * inv_2pi_l;
        let w1_local = -((panel_len - yp * beta) + xp * logterm - panel_len * logterm) * inv_2pi_l;
        let u2_local = (yp * logterm + xp * beta) * inv_2pi_l;
        let w2_local = ((panel_len - yp * beta) + xp * logterm) * inv_2pi_l;

        let uv1 = u1_local * cos_t - w1_local * sin_t;
        let vv1 = u1_local * sin_t + w1_local * cos_t;
        let uv2 = u2_local * cos_t - w2_local * sin_t;
        let vv2 = u2_local * sin_t + w2_local * cos_t;

        u += gamma[j] * uv1 + gamma[jp] * uv2;
        v += gamma[j] * vv1 + gamma[jp] * vv2;

        // Source influence: 90° rotation of vortex influence
        let us1_local = w1_local;
        let ws1_local = -u1_local;
        let us2_local = w2_local;
        let ws2_local = -u2_local;

        let us1 = us1_local * cos_t - ws1_local * sin_t;
        let vs1 = us1_local * sin_t + ws1_local * cos_t;
        let us2 = us2_local * cos_t - ws2_local * sin_t;
        let vs2 = us2_local * sin_t + ws2_local * cos_t;

        u += sigma[j] * us1 + sigma[jp] * us2;
        v += sigma[j] * vs1 + sigma[jp] * vs2;
    }

    // Wake panels: source-only (open polyline, not closed)
    if let Some(wake) = wake_panels {
        let nw = wake.x.len();
        if nw >= 2 && wake.sigma.len() == nw {
            for j in 0..(nw - 1) {
                let jp = j + 1;

                let x1 = wake.x[j];
                let y1 = wake.y[j];
                let x2 = wake.x[jp];
                let y2 = wake.y[jp];

                let dx = x2 - x1;
                let dy = y2 - y1;
                let panel_len = (dx * dx + dy * dy).sqrt();

                if panel_len < 1e-12 {
                    continue;
                }

                let theta = dy.atan2(dx);
                let cos_t = theta.cos();
                let sin_t = theta.sin();

                let xp = (x - x1) * cos_t + (y - y1) * sin_t;
                let yp = -(x - x1) * sin_t + (y - y1) * cos_t;
                let xp2 = xp - panel_len;

                let r1_sq = xp * xp + yp * yp;
                let r2_sq = xp2 * xp2 + yp * yp;

                if r1_sq < 1e-10 || r2_sq < 1e-10 {
                    continue;
                }

                let theta1 = yp.atan2(xp);
                let theta2 = yp.atan2(xp2);
                let beta = theta2 - theta1;
                let logterm = 0.5 * (r2_sq / r1_sq).ln();
                let inv_2pi_l = 1.0 / (two_pi * panel_len);

                let u1_local = -(yp * logterm + xp * beta - panel_len * beta) * inv_2pi_l;
                let w1_local = -((panel_len - yp * beta) + xp * logterm - panel_len * logterm) * inv_2pi_l;
                let u2_local = (yp * logterm + xp * beta) * inv_2pi_l;
                let w2_local = ((panel_len - yp * beta) + xp * logterm) * inv_2pi_l;

                let us1_local = w1_local;
                let ws1_local = -u1_local;
                let us2_local = w2_local;
                let ws2_local = -u2_local;

                let us1 = us1_local * cos_t - ws1_local * sin_t;
                let vs1 = us1_local * sin_t + ws1_local * cos_t;
                let us2 = us2_local * cos_t - ws2_local * sin_t;
                let vs2 = us2_local * sin_t + ws2_local * cos_t;

                u += wake.sigma[j] * us1 + wake.sigma[jp] * us2;
                v += wake.sigma[j] * vs1 + wake.sigma[jp] * vs2;
            }
        }
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
///
/// `nodes` is one closed contour. Use [`psi_at_multi`] for a node array that
/// concatenates several elements.
pub fn psi_at(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
) -> f64 {
    let layout = single_element_layout(nodes.len());
    psi_at_multi(x, y, nodes, gamma, alpha, v_inf, &layout)
}

/// Evaluate the stream function for a node array of several elements.
///
/// As [`psi_at`], but `layout` supplies the element connectivity, so no panel
/// is created across the gap between two elements, and the interior test is
/// [`is_inside_any_element`] rather than a single-contour one.
///
/// `layout` must describe `nodes`: if [`Layout::total_nodes`] disagrees with
/// `nodes.len()`, only the freestream stream function is returned.
pub fn psi_at_multi(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    layout: &Layout,
) -> f64 {
    let n = nodes.len();
    if n < 2 || gamma.len() != n || layout.total_nodes() != n {
        // Just freestream
        return v_inf * (alpha.cos() * y - alpha.sin() * x);
    }

    // Return NaN for points inside any element
    if is_inside_any_element(x, y, nodes, layout) {
        return f64::NAN;
    }

    let mut psi = 0.0;

    // Add contribution from each vortex panel
    // Panel jo goes from node jo to the next node in the same element
    for jo in 0..n {
        let jp = layout.next_node(jo);

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

/// Evaluate the stream function including source panel contributions.
///
/// Extends `psi_at` with source panels from viscous coupling.
/// For each panel, the source stream function for constant strength sigma is:
///   psi_source = sigma/(2*pi) * [x1*atan2(yy,x1) - x2*atan2(yy,x2) + 0.5*yy*(g1-g2)]
/// where (x1, x2, yy) are panel-local coordinates and g1, g2 are ln(r^2) terms.
///
/// `nodes` is one closed contour. Use [`psi_at_with_sources_multi`] for a node
/// array that concatenates several elements.
pub fn psi_at_with_sources(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    wake_panels: Option<&WakePanels>,
) -> f64 {
    let layout = single_element_layout(nodes.len());
    psi_at_with_sources_multi(
        x,
        y,
        nodes,
        gamma,
        sigma,
        alpha,
        v_inf,
        wake_panels,
        &layout,
    )
}

/// Evaluate the stream function including source panels, for a node array of
/// several elements.
///
/// As [`psi_at_with_sources`], but `layout` supplies the element connectivity,
/// so no panel is created across the gap between two elements, and the interior
/// test is [`is_inside_any_element`].
///
/// `layout` must describe `nodes`: if [`Layout::total_nodes`] disagrees with
/// `nodes.len()`, only the freestream stream function is returned.
#[allow(clippy::too_many_arguments)]
pub fn psi_at_with_sources_multi(
    x: f64,
    y: f64,
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    wake_panels: Option<&WakePanels>,
    layout: &Layout,
) -> f64 {
    let n = nodes.len();
    if n < 2 || gamma.len() != n || sigma.len() != n || layout.total_nodes() != n {
        return v_inf * (alpha.cos() * y - alpha.sin() * x);
    }

    if is_inside_any_element(x, y, nodes, layout) {
        return f64::NAN;
    }

    let mut psi = 0.0;
    let two_pi = 2.0 * PI;

    for jo in 0..n {
        let jp = layout.next_node(jo);

        let x_jo = nodes[jo].x;
        let y_jo = nodes[jo].y;
        let x_jp = nodes[jp].x;
        let y_jp = nodes[jp].y;

        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        if ds_sq < 1e-24 {
            continue;
        }

        let ds = ds_sq.sqrt();
        let sx = dx / ds;
        let sy = dy / ds;

        let rx1 = x - x_jo;
        let ry1 = y - y_jo;

        let x1 = sx * rx1 + sy * ry1;
        let x2 = x1 - ds;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = x1 * x1 + yy * yy;
        let rs2 = x2 * x2 + yy * yy;

        let sgn = if yy >= 0.0 { 1.0 } else { -1.0 };
        let pi_offset = (0.5 - 0.5 * sgn) * PI;

        let g1 = if rs1 > 1e-20 { rs1.ln() } else { 0.0 };
        let g2 = if rs2 > 1e-20 { rs2.ln() } else { 0.0 };
        let t1 = (sgn * x1).atan2(sgn * yy) + pi_offset;
        let t2 = (sgn * x2).atan2(sgn * yy) + pi_offset;

        // Vortex contribution (unchanged from psi_at)
        let dxinv = 1.0 / (x1 - x2);
        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

        let gsum = gamma[jp] + gamma[jo];
        let gdif = gamma[jp] - gamma[jo];
        psi += QOPI * (psis * gsum + psid * gdif);

        // Source contribution: constant-strength approximation per panel
        let sigma_avg = 0.5 * (sigma[jo] + sigma[jp]);
        if sigma_avg.abs() > 1e-20 {
            let at1 = yy.atan2(x1);
            let at2 = yy.atan2(x2);
            psi += sigma_avg / two_pi * (x1 * at1 - x2 * at2 + 0.5 * yy * (g1 - g2));
        }
    }

    // Wake panels: source-only (open polyline)
    if let Some(wake) = wake_panels {
        let nw = wake.x.len();
        if nw >= 2 && wake.sigma.len() == nw {
            for j in 0..(nw - 1) {
                let jp = j + 1;

                let dx = wake.x[jp] - wake.x[j];
                let dy = wake.y[jp] - wake.y[j];
                let ds_sq = dx * dx + dy * dy;

                if ds_sq < 1e-24 {
                    continue;
                }

                let ds = ds_sq.sqrt();
                let sx = dx / ds;
                let sy = dy / ds;

                let rx1 = x - wake.x[j];
                let ry1 = y - wake.y[j];

                let x1 = sx * rx1 + sy * ry1;
                let x2 = x1 - ds;
                let yy = sx * ry1 - sy * rx1;

                let rs1 = x1 * x1 + yy * yy;
                let rs2 = x2 * x2 + yy * yy;

                let g1 = if rs1 > 1e-20 { rs1.ln() } else { 0.0 };
                let g2 = if rs2 > 1e-20 { rs2.ln() } else { 0.0 };

                let sigma_avg = 0.5 * (wake.sigma[j] + wake.sigma[jp]);
                if sigma_avg.abs() > 1e-20 {
                    let at1 = yy.atan2(x1);
                    let at2 = yy.atan2(x2);
                    psi += sigma_avg / two_pi * (x1 * at1 - x2 * at2 + 0.5 * yy * (g1 - g2));
                }
            }
        }
    }

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

    // One contour, so the connectivity is the same at every grid point.
    let layout = single_element_layout(nodes.len());

    for iy in 0..ny {
        let y = y_min + iy as f64 * dy;
        for ix in 0..nx {
            let x = x_min + ix as f64 * dx;
            let psi = psi_at_multi(x, y, nodes, gamma, alpha, v_inf, &layout);
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

/// Compute stream function grid including source panel contributions.
pub fn compute_psi_grid_with_sources(
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    nx: usize,
    ny: usize,
    interior_value: Option<f64>,
    wake_panels: Option<&WakePanels>,
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

    // One contour, so the connectivity is the same at every grid point.
    let layout = single_element_layout(nodes.len());

    for iy in 0..ny {
        let y = y_min + iy as f64 * dy;
        for ix in 0..nx {
            let x = x_min + ix as f64 * dx;
            let psi = psi_at_with_sources_multi(
                x, y, nodes, gamma, sigma, alpha, v_inf, wake_panels, &layout,
            );
            grid[iy * nx + ix] = if psi.is_nan() {
                interior_value.unwrap_or(f64::NAN)
            } else {
                psi
            };
        }
    }

    grid
}

/// Check if a point lies inside a polygon using ray casting.
///
/// The polygon is closed implicitly: the edge from the last vertex back to the
/// first is included.
pub fn is_inside_polygon(x: f64, y: f64, polygon: &[Point]) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }

    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let xi = polygon[i].x;
        let yi = polygon[i].y;
        let xj = polygon[j].x;
        let yj = polygon[j].y;

        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }

    inside
}

/// Check if a point is inside the airfoil.
/// Uses ray casting algorithm.
///
/// `nodes` is one closed contour. Use [`is_inside_any_element`] for a node array
/// that concatenates several elements: casting one ray over all of them would
/// treat the whole configuration as a single polygon, whose implied edges bridge
/// the gaps between elements.
pub fn is_inside_airfoil(x: f64, y: f64, nodes: &[Point]) -> bool {
    is_inside_polygon(x, y, nodes)
}

/// Check if a point is inside **any** element of a configuration.
///
/// Each element's nodes are ray cast as their own closed contour, so a point in
/// the slot between two elements is outside both, and a point inside a slat or a
/// flap is reported as inside even though it is well clear of the main element.
///
/// The elements are tested in configuration order and the first hit wins, so the
/// cost is at worst one ray cast over `nodes` — the same work the single-contour
/// test does.
///
/// `layout` must describe `nodes`; nodes beyond `nodes.len()` are ignored.
pub fn is_inside_any_element(x: f64, y: f64, nodes: &[Point], layout: &Layout) -> bool {
    layout.spans().iter().any(|span| {
        let start = span.start().min(nodes.len());
        let end = span.end().min(nodes.len());
        is_inside_polygon(x, y, &nodes[start..end])
    })
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
///
/// `layout` must describe `nodes`; a streamline terminates on entering **any**
/// element, so it cannot be drawn through a slat or a flap.
#[allow(clippy::too_many_arguments)]
fn integrate_streamline<F>(
    field: &F,
    x0: f64,
    y0: f64,
    step_size: f64,
    max_steps: usize,
    nodes: &[Point],
    layout: &Layout,
    bounds: (f64, f64, f64, f64),
    effective_body: Option<&[Point]>,
) -> Vec<(f64, f64)>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    let (x_min, x_max, y_min, y_max) = bounds;
    let mut points = vec![(x0, y0)];
    let mut x = x0;
    let mut y = y0;

    for _ in 0..max_steps {
        // Check if inside any element
        if is_inside_any_element(x, y, nodes, layout) {
            break;
        }
        if effective_body.is_some_and(|poly| is_inside_polygon(x, y, poly)) {
            break;
        }
        
        // Check bounds
        if x < x_min || x > x_max || y < y_min || y > y_max {
            break;
        }
        
        match rk4_arc_step(field, x, y, step_size, 1e-4, 0.02) {
            Some((new_x, new_y)) => {
                if effective_body.is_some_and(|poly| is_inside_polygon(new_x, new_y, poly)) {
                    break;
                }
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
/// 
/// Seeds are placed along boundaries where flow enters the domain:
/// - Left boundary: always (flow generally goes left-to-right)
/// - Top boundary: when alpha > 0 (flow enters from top-left)
/// - Bottom boundary: when alpha < 0 (flow enters from bottom-left)
pub fn build_streamlines(
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    options: &StreamlineOptions,
) -> Vec<Vec<(f64, f64)>> {
    let layout = single_element_layout(nodes.len());

    // Create velocity field closure
    let field = |x: f64, y: f64| velocity_at_multi(x, y, nodes, gamma, alpha, v_inf, &layout);

    let bounds = (options.x_min, options.x_max, options.y_min, options.y_max);

    // Calculate how many seeds to allocate to each boundary based on alpha
    // At alpha=0, all seeds go to left boundary
    // As |alpha| increases, more seeds go to top/bottom boundaries
    let alpha_factor = alpha.abs().sin().min(0.8); // Cap at 80% to keep some left-edge seeds
    let left_count = ((1.0 - alpha_factor) * options.seed_count as f64).round() as usize;
    let edge_count = options.seed_count.saturating_sub(left_count);

    let mut streamlines = Vec::with_capacity(options.seed_count + edge_count);

    // Helper to add a streamline from a seed point
    let mut add_streamline = |x: f64, y: f64| {
        if is_inside_any_element(x, y, nodes, &layout) {
            return;
        }

        let streamline = integrate_streamline(
            &field,
            x,
            y,
            options.step_size,
            options.max_steps,
            nodes,
            &layout,
            bounds,
            None,
        );

        if streamline.len() >= 2 {
            streamlines.push(streamline);
        }
    };
    
    // 1. Seeds along left boundary (primary - always active)
    let left_seed_count = left_count.max(5); // At least 5 seeds on left
    for i in 0..left_seed_count {
        let t = i as f64 / (left_seed_count - 1).max(1) as f64;
        let y = options.y_min + t * (options.y_max - options.y_min);
        let x = options.seed_x;
        add_streamline(x, y);
    }
    
    // 2. Seeds along top/bottom boundaries based on flow direction
    if edge_count > 0 {
        // Freestream direction components
        let v_y = alpha.sin();
        
        if v_y > 0.01 {
            // Alpha > 0: flow enters from top boundary (y = y_max)
            // Seed along top edge, from left to somewhere past the airfoil
            let x_start = options.x_min;
            let x_end = 0.5; // Cover area around the airfoil
            for i in 0..edge_count {
                let t = i as f64 / edge_count.max(1) as f64;
                let x = x_start + t * (x_end - x_start);
                let y = options.y_max - 0.01; // Slightly inside the boundary
                add_streamline(x, y);
            }
        } else if v_y < -0.01 {
            // Alpha < 0: flow enters from bottom boundary (y = y_min)
            // Seed along bottom edge
            let x_start = options.x_min;
            let x_end = 0.5;
            for i in 0..edge_count {
                let t = i as f64 / edge_count.max(1) as f64;
                let x = x_start + t * (x_end - x_start);
                let y = options.y_min + 0.01; // Slightly inside the boundary
                add_streamline(x, y);
            }
        }
    }
    
    streamlines
}

/// `layout` must describe `nodes`; the body the streamline is classified against
/// is the whole configuration, and the closest-approach walk closes each
/// element's contour on itself.
fn build_dividing_streamline_internal<F>(
    field: &F,
    nodes: &[Point],
    layout: &Layout,
    effective_body: Option<&[Point]>,
    options: &StreamlineOptions,
) -> Option<Vec<(f64, f64)>>
where
    F: Fn(f64, f64) -> (f64, f64),
{
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum StreamlineSide {
        Above,
        Below,
    }

    #[derive(Debug, Clone)]
    struct TraceResult {
        seed_y: f64,
        side: Option<StreamlineSide>,
        closest_distance: f64,
        streamline: Vec<(f64, f64)>,
    }

    fn point_segment_distance_sq(px: f64, py: f64, a: Point, b: Point) -> f64 {
        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let len_sq = dx * dx + dy * dy;
        if len_sq <= 1e-16 {
            let ex = px - a.x;
            let ey = py - a.y;
            return ex * ex + ey * ey;
        }
        let t = (((px - a.x) * dx + (py - a.y) * dy) / len_sq).clamp(0.0, 1.0);
        let qx = a.x + t * dx;
        let qy = a.y + t * dy;
        let ex = px - qx;
        let ey = py - qy;
        ex * ex + ey * ey
    }

    fn classify_streamline(
        streamline: &[(f64, f64)],
        nodes: &[Point],
        layout: &Layout,
    ) -> Option<(Option<StreamlineSide>, f64)> {
        if streamline.len() < 2 || nodes.len() < 2 || layout.total_nodes() != nodes.len() {
            return None;
        }

        let body_x_min = nodes.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let body_x_max = nodes.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let body_y_center = nodes.iter().map(|p| p.y).sum::<f64>() / nodes.len() as f64;
        let chord = (body_x_max - body_x_min).abs().max(1e-6);
        let x_pad = 0.05 * chord;
        let x_ref = body_x_min + 0.65 * chord;

        let mut best: Option<(f64, f64)> = None;
        let consider_point = |best: &mut Option<(f64, f64)>, x: f64, y: f64| {
            let mut min_dist_sq = f64::INFINITY;
            for i in 0..nodes.len() {
                let a = nodes[i];
                // Closes on this element's own first node, so no segment spans
                // the gap between two elements.
                let b = nodes[layout.next_node(i)];
                min_dist_sq = min_dist_sq.min(point_segment_distance_sq(x, y, a, b));
            }
            let signed_offset = y - body_y_center;
            if best.as_ref().is_none_or(|(best_dist_sq, _)| min_dist_sq < *best_dist_sq) {
                *best = Some((min_dist_sq, signed_offset));
            }
        };

        for &(x, y) in streamline {
            if x >= body_x_min - x_pad && x <= body_x_max + x_pad {
                consider_point(&mut best, x, y);
            }
        }
        if best.is_none() {
            for &(x, y) in streamline {
                consider_point(&mut best, x, y);
            }
        }

        let (dist_sq, signed_offset) = best?;
        let closest_distance = dist_sq.sqrt();
        let end = *streamline.last()?;
        if end.0 <= body_x_min + 0.12 * chord && closest_distance <= 0.04 * chord {
            let end_offset = end.1 - body_y_center;
            if end_offset.abs() <= 0.02 * chord {
                return Some((None, closest_distance));
            }
            let side = if end_offset >= 0.0 {
                Some(StreamlineSide::Above)
            } else {
                Some(StreamlineSide::Below)
            };
            return Some((side, closest_distance));
        }

        let mut station_sample: Option<(f64, f64)> = None;
        for &(x, y) in streamline {
            let dx = (x - x_ref).abs();
            if station_sample
                .as_ref()
                .is_none_or(|(best_dx, _)| dx < *best_dx)
            {
                station_sample = Some((dx, y));
            }
        }
        let station_y = station_sample.map(|(_, y)| y).unwrap_or(body_y_center + signed_offset);
        let side = if station_y >= body_y_center {
            Some(StreamlineSide::Above)
        } else {
            Some(StreamlineSide::Below)
        };
        Some((side, closest_distance))
    }

    let bounds = (options.x_min, options.x_max, options.y_min, options.y_max);
    let sample_count = options.seed_count.max(25).min(129);
    let y_span = (options.y_max - options.y_min).abs();
    let seed_tol = (1e-4 * y_span).max(1e-5);
    let mut previous: Option<TraceResult> = None;

    let trace_seed = |seed_y: f64| -> Option<TraceResult> {
        if is_inside_any_element(options.seed_x, seed_y, nodes, layout)
            || effective_body.is_some_and(|poly| is_inside_polygon(options.seed_x, seed_y, poly))
        {
            return None;
        }

        let streamline = integrate_streamline(
            field,
            options.seed_x,
            seed_y,
            options.step_size,
            options.max_steps,
            nodes,
            layout,
            bounds,
            effective_body,
        );
        if streamline.len() < 2 {
            return None;
        }

        let (side, closest_distance) = classify_streamline(&streamline, nodes, layout)?;
        Some(TraceResult {
            seed_y,
            side,
            closest_distance,
            streamline,
        })
    };

    let mut bracket: Option<(TraceResult, TraceResult)> = None;

    for i in 0..sample_count {
        let t = i as f64 / (sample_count - 1).max(1) as f64;
        let seed_y = options.y_min + t * (options.y_max - options.y_min);
        let Some(trace) = trace_seed(seed_y) else {
            continue;
        };
        if trace.side.is_none() {
            return Some(trace.streamline);
        }

        if let Some(prev) = &previous {
            if prev.side.is_some() && trace.side.is_some() && trace.side != prev.side {
                bracket = Some((prev.clone(), trace));
                break;
            }
        }

        previous = Some(trace);
    }

    let Some((mut lo, mut hi)) = bracket else {
        return None;
    };
    if lo.seed_y > hi.seed_y {
        std::mem::swap(&mut lo, &mut hi);
    }
    let mut best = if lo.closest_distance <= hi.closest_distance {
        lo.clone()
    } else {
        hi.clone()
    };

    for _ in 0..32 {
        if (hi.seed_y - lo.seed_y).abs() <= seed_tol {
            break;
        }

        let y_mid = 0.5 * (lo.seed_y + hi.seed_y);
        let Some(mid) = trace_seed(y_mid) else {
            break;
        };
        if mid.side.is_none() {
            return Some(mid.streamline);
        }
        if mid.closest_distance < best.closest_distance {
            best = mid.clone();
        }

        if mid.side == lo.side {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    Some(best.streamline)
}

/// Build the streamline whose stream-function value brackets `psi_0`.
pub fn build_dividing_streamline(
    nodes: &[Point],
    gamma: &[f64],
    alpha: f64,
    v_inf: f64,
    _psi_0: f64,
    options: &StreamlineOptions,
) -> Option<Vec<(f64, f64)>> {
    let layout = single_element_layout(nodes.len());
    let field = |x: f64, y: f64| velocity_at_multi(x, y, nodes, gamma, alpha, v_inf, &layout);
    build_dividing_streamline_internal(&field, nodes, &layout, None, options)
}

/// Build streamlines using the viscous velocity field (vortex + source panels).
pub fn build_streamlines_viscous(
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    wake_panels: Option<&WakePanels>,
    effective_body: Option<&[Point]>,
    options: &StreamlineOptions,
) -> Vec<Vec<(f64, f64)>> {
    let layout = single_element_layout(nodes.len());
    let field = |x: f64, y: f64| {
        velocity_at_with_sources_multi(
            x, y, nodes, gamma, sigma, alpha, v_inf, wake_panels, &layout,
        )
    };

    let bounds = (options.x_min, options.x_max, options.y_min, options.y_max);

    let alpha_factor = alpha.abs().sin().min(0.8);
    let left_count = ((1.0 - alpha_factor) * options.seed_count as f64).round() as usize;
    let edge_count = options.seed_count.saturating_sub(left_count);

    let mut streamlines = Vec::with_capacity(options.seed_count + edge_count);

    let mut add_streamline = |x: f64, y: f64| {
        if is_inside_any_element(x, y, nodes, &layout) {
            return;
        }
        let streamline = integrate_streamline(
            &field,
            x,
            y,
            options.step_size,
            options.max_steps,
            nodes,
            &layout,
            bounds,
            effective_body,
        );
        if streamline.len() >= 2 {
            streamlines.push(streamline);
        }
    };

    let left_seed_count = left_count.max(5);
    for i in 0..left_seed_count {
        let t = i as f64 / (left_seed_count - 1).max(1) as f64;
        let y = options.y_min + t * (options.y_max - options.y_min);
        add_streamline(options.seed_x, y);
    }

    if edge_count > 0 {
        let v_y = alpha.sin();
        if v_y > 0.01 {
            for i in 0..edge_count {
                let t = i as f64 / edge_count.max(1) as f64;
                let x = options.x_min + t * (0.5 - options.x_min);
                add_streamline(x, options.y_max - 0.01);
            }
        } else if v_y < -0.01 {
            for i in 0..edge_count {
                let t = i as f64 / edge_count.max(1) as f64;
                let x = options.x_min + t * (0.5 - options.x_min);
                add_streamline(x, options.y_min + 0.01);
            }
        }
    }

    streamlines
}

/// Build the viscous dividing streamline by bracketing `psi_0` on the inflow edge
/// and bisecting between streamlines that lie above and below the separatrix.
pub fn build_dividing_streamline_viscous(
    nodes: &[Point],
    gamma: &[f64],
    sigma: &[f64],
    alpha: f64,
    v_inf: f64,
    _psi_0: f64,
    wake_panels: Option<&WakePanels>,
    effective_body: Option<&[Point]>,
    options: &StreamlineOptions,
) -> Option<Vec<(f64, f64)>> {
    let layout = single_element_layout(nodes.len());
    let field = |x: f64, y: f64| {
        velocity_at_with_sources_multi(
            x, y, nodes, gamma, sigma, alpha, v_inf, wake_panels, &layout,
        )
    };
    build_dividing_streamline_internal(&field, nodes, &layout, effective_body, options)
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

    #[test]
    fn test_build_dividing_streamline_brackets_circle_stagnation() {
        let circle = make_circle(96, 0.5);
        let gamma = vec![0.0; circle.len()];
        let options = StreamlineOptions {
            seed_count: 25,
            seed_x: -1.0,
            y_min: -1.0,
            y_max: 1.0,
            step_size: 0.01,
            max_steps: 400,
            x_min: -1.2,
            x_max: 1.5,
        };

        let streamline = build_dividing_streamline(&circle, &gamma, 0.0, 1.0, 0.0, &options)
            .expect("expected dividing streamline");

        let seed = streamline.first().copied().expect("seed point");
        let last = streamline.last().copied().expect("terminal point");

        assert!(seed.1.abs() < 1e-3, "seed should converge to y=0, got {}", seed.1);
        assert!(
            (last.0 + 0.5).abs() < 0.05,
            "streamline should end near stagnation x=-0.5, got {}",
            last.0
        );
        assert!(last.1.abs() < 0.05, "streamline should remain near y=0, got {}", last.1);
    }

    #[test]
    fn test_build_dividing_streamline_viscous_matches_zero_source_case() {
        let circle = make_circle(96, 0.5);
        let gamma = vec![0.0; circle.len()];
        let sigma = vec![0.0; circle.len()];
        let options = StreamlineOptions {
            seed_count: 25,
            seed_x: -1.0,
            y_min: -1.0,
            y_max: 1.0,
            step_size: 0.01,
            max_steps: 400,
            x_min: -1.2,
            x_max: 1.5,
        };

        let streamline = build_dividing_streamline_viscous(
            &circle,
            &gamma,
            &sigma,
            0.0,
            1.0,
            0.0,
            None,
            None,
            &options,
        )
        .expect("expected viscous dividing streamline");

        let seed = streamline.first().copied().expect("seed point");
        assert!(seed.1.abs() < 1e-3, "seed should converge to y=0, got {}", seed.1);
    }

    // --- element-aware connectivity --------------------------------------

    /// A circle of `n` nodes centred on `(cx, cy)`.
    fn make_circle_at(n: usize, radius: f64, cx: f64, cy: f64) -> Vec<Point> {
        (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                point(cx + radius * theta.cos(), cy + radius * theta.sin())
            })
            .collect()
    }

    /// Two well-separated circles, concatenated as one node array, with the
    /// layout that describes them.
    fn two_bodies() -> (Vec<Point>, Layout) {
        let mut nodes = make_circle_at(48, 0.3, 0.0, 0.0);
        nodes.extend(make_circle_at(32, 0.2, 1.5, 0.4));
        let layout = Layout::from_node_counts(&[48, 32]).unwrap();
        assert_eq!(layout.total_nodes(), nodes.len());
        (nodes, layout)
    }

    #[test]
    fn single_element_layout_reproduces_the_modulo_closure() {
        // The backward-compatibility guarantee the single-airfoil entry points
        // rest on: with one element the successor is the closed-contour wrap.
        for n in [2usize, 3, 61, 160] {
            let layout = single_element_layout(n);
            assert_eq!(layout.n_elements(), 1);
            assert_eq!(layout.total_nodes(), n);
            for i in 0..n {
                assert_eq!(layout.next_node(i), (i + 1) % n, "n = {n}, i = {i}");
            }
        }
        // An empty node array has no element at all, rather than an empty one.
        assert_eq!(single_element_layout(0).n_elements(), 0);
        assert_eq!(single_element_layout(0).total_nodes(), 0);
    }

    #[test]
    fn multi_element_velocity_is_the_sum_of_the_separate_bodies() {
        // The property a modulo-n closure breaks: with two elements in one node
        // array, the induced velocity must be what the two bodies induce
        // separately. A panel spanning the gap between them would add a
        // contribution that neither body has.
        let first = make_circle_at(48, 0.3, 0.0, 0.0);
        let second = make_circle_at(32, 0.2, 1.5, 0.4);
        let (nodes, layout) = two_bodies();

        let gamma_first: Vec<f64> = (0..first.len()).map(|i| 0.1 + 0.01 * i as f64).collect();
        let gamma_second: Vec<f64> = (0..second.len()).map(|i| -0.2 + 0.02 * i as f64).collect();
        let mut gamma = gamma_first.clone();
        gamma.extend_from_slice(&gamma_second);

        let alpha = 4.0_f64.to_radians();
        let v_inf = 1.0;
        let (u_free, v_free) = (v_inf * alpha.cos(), v_inf * alpha.sin());

        for &(x, y) in &[(-0.8, 0.0), (0.7, 0.2), (1.5, -0.6), (2.4, 0.9)] {
            let (u, v) = velocity_at_multi(x, y, &nodes, &gamma, alpha, v_inf, &layout);

            let (u1, v1) = velocity_at(x, y, &first, &gamma_first, alpha, v_inf);
            let (u2, v2) = velocity_at(x, y, &second, &gamma_second, alpha, v_inf);
            // Each single-body call carries the freestream, so remove one copy.
            let u_expected = u1 + u2 - u_free;
            let v_expected = v1 + v2 - v_free;

            assert!(
                (u - u_expected).abs() < 1e-12 && (v - v_expected).abs() < 1e-12,
                "at ({x}, {y}): got ({u}, {v}), expected ({u_expected}, {v_expected})"
            );
        }
    }

    #[test]
    fn multi_element_velocity_differs_from_treating_the_array_as_one_contour() {
        // And the two really are different, so the closure rule is doing work:
        // one contour over the same nodes adds two panels across the gap.
        let (nodes, layout) = two_bodies();
        let gamma: Vec<f64> = (0..nodes.len()).map(|i| 0.1 + 0.01 * i as f64).collect();

        let (u_multi, _) = velocity_at_multi(0.7, 0.2, &nodes, &gamma, 0.0, 1.0, &layout);
        let (u_single, _) = velocity_at(0.7, 0.2, &nodes, &gamma, 0.0, 1.0);
        assert!(
            (u_multi - u_single).abs() > 1e-6,
            "expected the gap-spanning panels to change the answer"
        );
    }

    #[test]
    fn inside_test_covers_every_element_and_not_the_gap() {
        let (nodes, layout) = two_bodies();

        // Inside the first element, inside the second, and in the gap between.
        assert!(is_inside_any_element(0.0, 0.0, &nodes, &layout));
        assert!(is_inside_any_element(1.5, 0.4, &nodes, &layout));
        assert!(!is_inside_any_element(0.75, 0.2, &nodes, &layout));
        assert!(!is_inside_any_element(-2.0, 0.0, &nodes, &layout));

        // The single-contour test cannot answer this geometry: sweeping a grid,
        // the two disagree somewhere, which is why the element-aware form
        // exists.
        let mut disagreements = 0usize;
        for iy in 0..41 {
            let y = -1.0 + 0.05 * iy as f64;
            for ix in 0..61 {
                let x = -1.0 + 0.05 * ix as f64;
                if is_inside_airfoil(x, y, &nodes) != is_inside_any_element(x, y, &nodes, &layout) {
                    disagreements += 1;
                }
            }
        }
        assert!(
            disagreements > 0,
            "the single-contour and per-element tests should not agree everywhere"
        );
    }

    #[test]
    fn inside_test_matches_the_single_contour_test_for_one_element() {
        // One element: the two must agree at every point, which is what keeps
        // single-airfoil ψ grids and streamline termination unchanged.
        let circle = make_circle(64, 0.4);
        let layout = single_element_layout(circle.len());
        for iy in 0..21 {
            let y = -0.6 + 0.06 * iy as f64;
            for ix in 0..21 {
                let x = -0.6 + 0.06 * ix as f64;
                assert_eq!(
                    is_inside_any_element(x, y, &circle, &layout),
                    is_inside_airfoil(x, y, &circle),
                    "at ({x}, {y})"
                );
            }
        }
    }

    #[test]
    fn psi_is_nan_inside_either_element() {
        let (nodes, layout) = two_bodies();
        let gamma = vec![0.0; nodes.len()];

        assert!(psi_at_multi(0.0, 0.0, &nodes, &gamma, 0.0, 1.0, &layout).is_nan());
        assert!(psi_at_multi(1.5, 0.4, &nodes, &gamma, 0.0, 1.0, &layout).is_nan());
        assert!(psi_at_multi(0.75, 0.2, &nodes, &gamma, 0.0, 1.0, &layout).is_finite());
    }

    #[test]
    fn a_streamline_stops_at_the_second_element() {
        // Streamline termination is the visible half of the inside test: a line
        // aimed at the downstream element must stop on it rather than being
        // drawn straight through.
        let (nodes, layout) = two_bodies();
        let field = |_x: f64, _y: f64| (1.0, 0.0);

        let streamline = integrate_streamline(
            &field,
            0.6,
            0.4,
            0.01,
            2000,
            &nodes,
            &layout,
            (-2.0, 3.0, -2.0, 2.0),
            None,
        );

        let last = *streamline.last().expect("at least the seed point");
        assert!(
            last.0 < 1.5,
            "should stop on the second element's upstream side, got x = {}",
            last.0
        );
        // Integration stops on the step after entering a body, so the final
        // point may sit just inside it; every earlier point is outside.
        for &(x, y) in &streamline[..streamline.len() - 1] {
            assert!(
                !is_inside_any_element(x, y, &nodes, &layout),
                "streamline passed through an element at ({x}, {y})"
            );
        }
    }

    #[test]
    fn a_mismatched_layout_falls_back_to_the_freestream() {
        // The layout has to describe the nodes; a disagreement is reported as
        // no panels rather than by indexing into the wrong element.
        let circle = make_circle(32, 0.4);
        let gamma = vec![0.5; circle.len()];
        let wrong = Layout::from_node_counts(&[10, 10]).unwrap();

        let (u, v) = velocity_at_multi(2.0, 0.0, &circle, &gamma, 0.0, 1.0, &wrong);
        assert_eq!((u, v), (1.0, 0.0));
        assert_eq!(psi_at_multi(2.0, 0.5, &circle, &gamma, 0.0, 1.0, &wrong), 0.5);
    }
}
