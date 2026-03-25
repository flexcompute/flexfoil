// Convert conservative variables Q to primitive variables W.
// Q = [rho, rho*u, rho*v, E, nu_tilde]
// W = [rho, u, v, p, nu_tilde]

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> Q: array<f32>;
@group(0) @binding(2) var<storage, read_write> W: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn conservative_to_primitive(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni || j >= nj) {
        return;
    }

    let base = q_idx(i, j, 0u, ni);
    let rho = max(Q[base + 0u], 1e-10);
    let rhou = Q[base + 1u];
    let rhov = Q[base + 2u];
    let E = Q[base + 3u];
    let nu_t = Q[base + 4u];

    let inv_rho = 1.0 / rho;
    let u = rhou * inv_rho;
    let v = rhov * inv_rho;
    let p = max((params.gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v)), 1e-10);

    W[base + 0u] = rho;
    W[base + 1u] = u;
    W[base + 2u] = v;
    W[base + 3u] = p;
    W[base + 4u] = nu_t;
}
