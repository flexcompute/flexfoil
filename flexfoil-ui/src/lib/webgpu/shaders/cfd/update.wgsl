// Update solution: Q_new = Q_old + dt * residual
// For explicit forward Euler time stepping.

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> residual: array<f32>;
@group(0) @binding(2) var<storage, read_write> Q: array<f32>;
@group(0) @binding(3) var<storage, read> bc_type: array<u32>;

@compute @workgroup_size(16, 16, 1)
fn update_solution(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni || j >= nj) {
        return;
    }

    // Only update interior cells (boundaries handled by BC shader)
    let cell = cell_idx(i, j, ni);
    let bc = bc_type[cell];
    if (bc != BC_INTERIOR) {
        return;
    }

    let dt = params.dt;
    for (var v: u32 = 0u; v < NVAR; v = v + 1u) {
        let idx = q_idx(i, j, v, ni);
        Q[idx] = Q[idx] + dt * residual[idx];
    }

    // Enforce positivity of density and pressure
    let base = q_idx(i, j, 0u, ni);
    Q[base + 0u] = max(Q[base + 0u], 1e-10); // rho > 0

    // Recompute pressure to ensure positivity
    let rho = Q[base + 0u];
    let u = Q[base + 1u] / rho;
    let v2 = Q[base + 2u] / rho;
    let E = Q[base + 3u];
    let p = (params.gamma - 1.0) * (E - 0.5 * rho * (u * u + v2 * v2));
    if (p < 1e-10) {
        // Reset energy to maintain minimum pressure
        Q[base + 3u] = 1e-10 / (params.gamma - 1.0) + 0.5 * rho * (u * u + v2 * v2);
    }
}
