/**
 * Velocity Field Compute Shader
 * 
 * Computes velocity (u, v), speed |V|, and stream function ψ at each grid point.
 * Port of velocity_at() and psi_at() from rustfoil-solver/src/inviscid/velocity.rs
 * 
 * Output texture RGBA32Float:
 *   R = u (x-velocity)
 *   G = v (y-velocity)  
 *   B = |V| (speed)
 *   A = ψ (stream function)
 */

const PI: f32 = 3.14159265358979323846;
const TWO_PI: f32 = 6.28318530717958647692;
const QOPI: f32 = 0.07957747154594767; // 1/(4π)

// Uniforms
struct Params {
    // Grid bounds
    x_min: f32,
    x_max: f32,
    y_min: f32,
    y_max: f32,
    // Flow conditions
    alpha: f32,      // Angle of attack (radians)
    v_inf: f32,      // Freestream velocity
    // Grid size
    nx: u32,
    ny: u32,
    // Panel count
    n_panels: u32,
    _pad: u32,
}

// Panel node: (x, y) position
struct Node {
    x: f32,
    y: f32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> nodes: array<Node>;
@group(0) @binding(2) var<storage, read> gamma: array<f32>;
@group(0) @binding(3) var output_texture: texture_storage_2d<rgba32float, write>;

/**
 * Evaluate velocity at point (x, y) using K&P VOR2DL linear vorticity panel method.
 */
fn velocity_at(x: f32, y: f32) -> vec2f {
    let n = params.n_panels;
    if (n < 2u) {
        return vec2f(params.v_inf * cos(params.alpha), params.v_inf * sin(params.alpha));
    }

    // Freestream velocity
    var u = params.v_inf * cos(params.alpha);
    var v = params.v_inf * sin(params.alpha);

    // Add contribution from each panel
    for (var j = 0u; j < n; j++) {
        let jp = (j + 1u) % n;

        let x1 = nodes[j].x;
        let y1 = nodes[j].y;
        let x2 = nodes[jp].x;
        let y2 = nodes[jp].y;

        let dx = x2 - x1;
        let dy = y2 - y1;
        let panel_len = sqrt(dx * dx + dy * dy);

        if (panel_len < 1e-12) {
            continue;
        }

        // Panel angle
        let theta = atan2(dy, dx);
        let cos_t = cos(theta);
        let sin_t = sin(theta);

        // Transform to panel-local coordinates
        let xp = (x - x1) * cos_t + (y - y1) * sin_t;
        let yp = -(x - x1) * sin_t + (y - y1) * cos_t;
        let xp2 = xp - panel_len;

        // Squared distances
        let r1_sq = xp * xp + yp * yp;
        let r2_sq = xp2 * xp2 + yp * yp;

        // Skip if too close
        if (r1_sq < 1e-10 || r2_sq < 1e-10) {
            continue;
        }

        // Angles from field point to panel endpoints
        let theta1 = atan2(yp, xp);
        let theta2 = atan2(yp, xp2);
        let beta = theta2 - theta1;

        // Log term
        let logterm = 0.5 * log(r2_sq / r1_sq);

        let inv_2pi_l = 1.0 / (TWO_PI * panel_len);

        // K&P VOR2DL influence coefficients
        let u1_local = -(yp * logterm + xp * beta - panel_len * beta) * inv_2pi_l;
        let w1_local = -((panel_len - yp * beta) + xp * logterm - panel_len * logterm) * inv_2pi_l;

        let u2_local = (yp * logterm + xp * beta) * inv_2pi_l;
        let w2_local = ((panel_len - yp * beta) + xp * logterm) * inv_2pi_l;

        // Transform back to global
        let u1 = u1_local * cos_t - w1_local * sin_t;
        let v1 = u1_local * sin_t + w1_local * cos_t;
        let u2 = u2_local * cos_t - w2_local * sin_t;
        let v2 = u2_local * sin_t + w2_local * cos_t;

        // Add weighted by gamma
        u += gamma[j] * u1 + gamma[jp] * u2;
        v += gamma[j] * v1 + gamma[jp] * v2;
    }

    return vec2f(u, v);
}

/**
 * Check if point is inside airfoil using ray casting.
 */
fn is_inside_airfoil(x: f32, y: f32) -> bool {
    let n = params.n_panels;
    if (n < 3u) {
        return false;
    }

    var inside = false;
    var j = n - 1u;

    for (var i = 0u; i < n; i++) {
        let xi = nodes[i].x;
        let yi = nodes[i].y;
        let xj = nodes[j].x;
        let yj = nodes[j].y;

        if (((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            inside = !inside;
        }
        j = i;
    }

    return inside;
}

/**
 * Evaluate stream function at point (x, y) using XFOIL's PSILIN formulation.
 */
fn psi_at(x: f32, y: f32) -> f32 {
    let n = params.n_panels;
    if (n < 2u) {
        return params.v_inf * (cos(params.alpha) * y - sin(params.alpha) * x);
    }

    // Return NaN equivalent (large value) for interior points
    if (is_inside_airfoil(x, y)) {
        return 1e10;  // Will be masked in rendering
    }

    var psi: f32 = 0.0;

    for (var jo = 0u; jo < n; jo++) {
        let jp = (jo + 1u) % n;

        let x_jo = nodes[jo].x;
        let y_jo = nodes[jo].y;
        let x_jp = nodes[jp].x;
        let y_jp = nodes[jp].y;

        let dx = x_jp - x_jo;
        let dy = y_jp - y_jo;
        let ds_sq = dx * dx + dy * dy;

        if (ds_sq < 1e-24) {
            continue;
        }

        let ds = sqrt(ds_sq);
        let sx = dx / ds;
        let sy = dy / ds;

        let rx1 = x - x_jo;
        let ry1 = y - y_jo;

        // Transform to panel-local coordinates
        let x1 = sx * rx1 + sy * ry1;
        let x2 = x1 - ds;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = x1 * x1 + yy * yy;
        let rs2 = x2 * x2 + yy * yy;

        // XFOIL's reflection correction
        var sgn: f32 = 1.0;
        var pi_offset: f32 = 0.0;
        if (yy < 0.0) {
            sgn = -1.0;
            pi_offset = PI;
        }

        // Log and arctangent terms
        var g1: f32 = 0.0;
        var g2: f32 = 0.0;
        if (rs1 > 1e-20) { g1 = log(rs1); }
        if (rs2 > 1e-20) { g2 = log(rs2); }
        let t1 = atan2(sgn * x1, sgn * yy) + pi_offset;
        let t2 = atan2(sgn * x2, sgn * yy) + pi_offset;

        // PSIS/PSID formulation
        let dxinv = 1.0 / (x1 - x2);
        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

        let gsum = gamma[jp] + gamma[jo];
        let gdif = gamma[jp] - gamma[jo];

        psi += QOPI * (psis * gsum + psid * gdif);
    }

    // Freestream contribution
    psi += params.v_inf * (cos(params.alpha) * y - sin(params.alpha) * x);

    return psi;
}

@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) id: vec3u) {
    let ix = id.x;
    let iy = id.y;

    if (ix >= params.nx || iy >= params.ny) {
        return;
    }

    // Compute world position
    let x = params.x_min + f32(ix) * (params.x_max - params.x_min) / f32(params.nx - 1u);
    let y = params.y_min + f32(iy) * (params.y_max - params.y_min) / f32(params.ny - 1u);

    // Compute velocity
    let vel = velocity_at(x, y);
    let speed = length(vel);

    // Compute stream function
    let psi = psi_at(x, y);

    // Write to texture
    textureStore(output_texture, vec2i(i32(ix), i32(iy)), vec4f(vel.x, vel.y, speed, psi));
}
