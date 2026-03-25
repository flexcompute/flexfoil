// Roe approximate Riemann solver.
// Computes numerical flux at each cell face from reconstructed left/right primitive states.
//
// Two entry points: roe_flux_xi (xi-direction faces) and roe_flux_eta (eta-direction faces).

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> metrics: array<f32>;
@group(0) @binding(2) var<storage, read> left_xi: array<f32>;
@group(0) @binding(3) var<storage, read> right_xi: array<f32>;
@group(0) @binding(4) var<storage, read> left_eta: array<f32>;
@group(0) @binding(5) var<storage, read> right_eta: array<f32>;
@group(0) @binding(6) var<storage, read_write> flux_xi: array<f32>;
@group(0) @binding(7) var<storage, read_write> flux_eta: array<f32>;

const EPSROE: f32 = 0.1;

// Compute Roe flux in a given direction (nx, ny = metric-weighted face normal).
// W_L, W_R are primitive states [rho, u, v, p, nu_tilde].
// Returns flux vector [F0, F1, F2, F3, F4].
fn roe_flux_1d(
    rhoL: f32, uL: f32, vL: f32, pL: f32, nuL: f32,
    rhoR: f32, uR: f32, vR: f32, pR: f32, nuR: f32,
    nx: f32, ny: f32, gamma: f32
) -> array<f32, 5> {
    // Roe-averaged quantities
    let sqrtL = sqrt(max(rhoL, 1e-10));
    let sqrtR = sqrt(max(rhoR, 1e-10));
    let inv_sum = 1.0 / (sqrtL + sqrtR);

    let u_roe = (sqrtL * uL + sqrtR * uR) * inv_sum;
    let v_roe = (sqrtL * vL + sqrtR * vR) * inv_sum;

    let hL = (gamma * pL / ((gamma - 1.0) * rhoL)) + 0.5 * (uL * uL + vL * vL);
    let hR = (gamma * pR / ((gamma - 1.0) * rhoR)) + 0.5 * (uR * uR + vR * vR);
    let h_roe = (sqrtL * hL + sqrtR * hR) * inv_sum;

    let q2_roe = u_roe * u_roe + v_roe * v_roe;
    let c2_roe = max((gamma - 1.0) * (h_roe - 0.5 * q2_roe), 1e-10);
    let c_roe = sqrt(c2_roe);

    // Contravariant velocity
    let face_mag = sqrt(nx * nx + ny * ny);
    let face_inv = select(1.0 / face_mag, 0.0, face_mag < 1e-15);
    let nx_hat = nx * face_inv;
    let ny_hat = ny * face_inv;

    let V_roe = u_roe * nx_hat + v_roe * ny_hat;

    // Eigenvalues with entropy fix
    let lam1 = abs(V_roe - c_roe * face_mag);
    let lam2 = abs(V_roe);
    let lam3 = abs(V_roe + c_roe * face_mag);

    let eps = EPSROE * c_roe * face_mag;
    let l1 = select(lam1, (lam1 * lam1 + eps * eps) / (2.0 * eps), lam1 < eps);
    let l2 = select(lam2, (lam2 * lam2 + eps * eps) / (2.0 * eps), lam2 < eps);
    let l3 = select(lam3, (lam3 * lam3 + eps * eps) / (2.0 * eps), lam3 < eps);

    // Jump in primitive variables
    let drho = rhoR - rhoL;
    let du = uR - uL;
    let dv = vR - vL;
    let dp = pR - pL;
    let dnu = nuR - nuL;

    let rho_roe = sqrtL * sqrtR;
    let dVn = du * nx_hat + dv * ny_hat;

    // Wave strengths
    let w1 = (dp - rho_roe * c_roe * face_mag * dVn) / (2.0 * c2_roe);
    let w2 = drho - dp / c2_roe;
    let w3 = (dp + rho_roe * c_roe * face_mag * dVn) / (2.0 * c2_roe);

    // Physical flux from left and right states
    let VnL = uL * nx + vL * ny;
    let VnR = uR * nx + vR * ny;

    let eL = pL / (gamma - 1.0) + 0.5 * rhoL * (uL * uL + vL * vL);
    let eR = pR / (gamma - 1.0) + 0.5 * rhoR * (uR * uR + vR * vR);

    // Left flux
    var fL: array<f32, 5>;
    fL[0] = rhoL * VnL;
    fL[1] = rhoL * uL * VnL + pL * nx;
    fL[2] = rhoL * vL * VnL + pL * ny;
    fL[3] = (eL + pL) * VnL;
    fL[4] = rhoL * nuL * VnL;

    // Right flux
    var fR: array<f32, 5>;
    fR[0] = rhoR * VnR;
    fR[1] = rhoR * uR * VnR + pR * nx;
    fR[2] = rhoR * vR * VnR + pR * ny;
    fR[3] = (eR + pR) * VnR;
    fR[4] = rhoR * nuR * VnR;

    // Roe dissipation terms
    // Wave 1 (u - c)
    var d: array<f32, 5>;
    d[0] = l1 * w1;
    d[1] = l1 * w1 * (u_roe - c_roe * nx_hat);
    d[2] = l1 * w1 * (v_roe - c_roe * ny_hat);
    d[3] = l1 * w1 * (h_roe - c_roe * face_mag * V_roe);
    d[4] = 0.0;

    // Wave 2 (u) - entropy wave + shear
    d[0] += l2 * w2;
    d[1] += l2 * (w2 * u_roe + rho_roe * (du - dVn * nx_hat));
    d[2] += l2 * (w2 * v_roe + rho_roe * (dv - dVn * ny_hat));
    d[3] += l2 * (w2 * 0.5 * q2_roe + rho_roe * (u_roe * du + v_roe * dv - V_roe * dVn));
    d[4] += l2 * (drho * nuR + rho_roe * dnu); // SA variable: passive advection

    // Wave 3 (u + c)
    d[0] += l3 * w3;
    d[1] += l3 * w3 * (u_roe + c_roe * nx_hat);
    d[2] += l3 * w3 * (v_roe + c_roe * ny_hat);
    d[3] += l3 * w3 * (h_roe + c_roe * face_mag * V_roe);

    // Roe flux = 0.5 * (F_L + F_R - dissipation)
    var flux: array<f32, 5>;
    for (var k: u32 = 0u; k < 5u; k = k + 1u) {
        flux[k] = 0.5 * (fL[k] + fR[k] - d[k]);
    }

    return flux;
}

@compute @workgroup_size(16, 16, 1)
fn roe_flux_xi(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni || j >= nj) {
        return;
    }

    // Get left/right primitive states at face (i+1/2, j)
    let face = (j * ni + i) * NVAR;
    let rhoL = left_xi[face + 0u]; let uL = left_xi[face + 1u];
    let vL = left_xi[face + 2u]; let pL = left_xi[face + 3u]; let nuL = left_xi[face + 4u];
    let rhoR = right_xi[face + 0u]; let uR = right_xi[face + 1u];
    let vR = right_xi[face + 2u]; let pR = right_xi[face + 3u]; let nuR = right_xi[face + 4u];

    // Face normal in xi-direction: (xi_x, xi_y) from metrics at cell (i, j)
    let m_base = cell_idx(i, j, ni) * 5u;
    let nx = metrics[m_base + 0u]; // xi_x
    let ny = metrics[m_base + 1u]; // xi_y

    let flux = roe_flux_1d(rhoL, uL, vL, pL, nuL, rhoR, uR, vR, pR, nuR, nx, ny, params.gamma);

    for (var k: u32 = 0u; k < 5u; k = k + 1u) {
        flux_xi[face + k] = flux[k];
    }
}

@compute @workgroup_size(16, 16, 1)
fn roe_flux_eta(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni || j >= nj - 1u) {
        return;
    }

    let face = (j * ni + i) * NVAR;
    let rhoL = left_eta[face + 0u]; let uL = left_eta[face + 1u];
    let vL = left_eta[face + 2u]; let pL = left_eta[face + 3u]; let nuL = left_eta[face + 4u];
    let rhoR = right_eta[face + 0u]; let uR = right_eta[face + 1u];
    let vR = right_eta[face + 2u]; let pR = right_eta[face + 3u]; let nuR = right_eta[face + 4u];

    // Face normal in eta-direction: (eta_x, eta_y) from metrics at cell (i, j)
    let m_base = cell_idx(i, j, ni) * 5u;
    let nx = metrics[m_base + 2u]; // eta_x
    let ny = metrics[m_base + 3u]; // eta_y

    let flux = roe_flux_1d(rhoL, uL, vL, pL, nuL, rhoR, uR, vR, pR, nuR, nx, ny, params.gamma);

    for (var k: u32 = 0u; k < 5u; k = k + 1u) {
        flux_eta[face + k] = flux[k];
    }
}
