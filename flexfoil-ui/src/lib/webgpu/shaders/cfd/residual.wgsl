// Compute the RHS residual: R(Q) = -(1/J) * (dF/dxi + dG/deta)
// For explicit time stepping: Q_new = Q_old + dt * R(Q)

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> metrics: array<f32>;
@group(0) @binding(2) var<storage, read> flux_xi: array<f32>;
@group(0) @binding(3) var<storage, read> flux_eta: array<f32>;
@group(0) @binding(4) var<storage, read_write> residual: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn compute_residual(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni || j >= nj) {
        return;
    }

    // Get Jacobian at this cell
    let m_base = cell_idx(i, j, ni) * 5u;
    let J = metrics[m_base + 4u];
    let inv_J = select(1.0 / J, 0.0, abs(J) < 1e-15);

    for (var v: u32 = 0u; v < NVAR; v = v + 1u) {
        // xi-direction flux difference: F(i+1/2) - F(i-1/2)
        // Face (i) stores flux at face i+1/2
        let face_ip = (j * ni + i) * NVAR + v;
        let im = wrap_i(i32(i) - 1, ni);
        let face_im = (j * ni + im) * NVAR + v;

        let dF_dxi = flux_xi[face_ip] - flux_xi[face_im];

        // eta-direction flux difference: G(j+1/2) - G(j-1/2)
        var dG_deta: f32 = 0.0;
        if (j > 0u && j < nj - 1u) {
            let face_jp = (j * ni + i) * NVAR + v;
            let face_jm = ((j - 1u) * ni + i) * NVAR + v;
            dG_deta = flux_eta[face_jp] - flux_eta[face_jm];
        } else if (j == 0u) {
            // At wall: only outward flux
            let face_jp = (0u * ni + i) * NVAR + v;
            dG_deta = flux_eta[face_jp]; // G(1/2) - 0 (wall flux handled by BC)
        }

        // Residual: R = -(1/J) * (dF/dxi + dG/deta)
        let r_idx = q_idx(i, j, v, ni);
        residual[r_idx] = -inv_J * (dF_dxi + dG_deta);
    }
}
