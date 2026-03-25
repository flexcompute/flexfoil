// 2nd-order MUSCL reconstruction with minmod limiter.
// Reconstructs left/right states at cell faces for the Roe solver.
//
// Two entry points: reconstruct_xi (along i-direction) and reconstruct_eta (along j-direction).
// Each workgroup processes one grid line.

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> W: array<f32>;
@group(0) @binding(2) var<storage, read_write> left_xi: array<f32>;
@group(0) @binding(3) var<storage, read_write> right_xi: array<f32>;
@group(0) @binding(4) var<storage, read_write> left_eta: array<f32>;
@group(0) @binding(5) var<storage, read_write> right_eta: array<f32>;

fn minmod(a: f32, b: f32) -> f32 {
    if (a * b <= 0.0) {
        return 0.0;
    }
    return select(b, a, abs(a) < abs(b));
}

// Get primitive variable v at wrapped cell (i, j)
fn get_w(i: i32, j: i32, v: u32, ni: u32, nj: u32) -> f32 {
    let ii = wrap_i(i, ni);
    let jj = clamp(u32(j), 0u, nj - 1u);
    return W[q_idx(ii, jj, v, ni)];
}

// Reconstruct along xi (i-direction) for a fixed j-line.
// Face (i+1/2, j) has left state from cell (i) and right state from cell (i+1).
@compute @workgroup_size(1, 1, 1)
fn reconstruct_xi(@builtin(global_invocation_id) gid: vec3<u32>) {
    let j = gid.x;
    let ni = params.ni;
    let nj = params.nj;

    if (j >= nj) {
        return;
    }

    // Process each face along this j-line
    for (var i: u32 = 0u; i < ni; i = i + 1u) {
        for (var v: u32 = 0u; v < NVAR; v = v + 1u) {
            // Get stencil values
            let wim1 = get_w(i32(i) - 1, i32(j), v, ni, nj);
            let wi   = get_w(i32(i),     i32(j), v, ni, nj);
            let wip1 = get_w(i32(i) + 1, i32(j), v, ni, nj);
            let wip2 = get_w(i32(i) + 2, i32(j), v, ni, nj);

            // Slopes
            let dL = wi - wim1;
            let dR = wip1 - wi;
            let dL2 = wip1 - wi;
            let dR2 = wip2 - wip1;

            // MUSCL reconstruction with minmod
            let slope_L = minmod(dL, dR);
            let slope_R = minmod(dL2, dR2);

            // Left state at face i+1/2 (extrapolated from cell i)
            let wL = wi + 0.5 * slope_L;
            // Right state at face i+1/2 (extrapolated from cell i+1)
            let wR = wip1 - 0.5 * slope_R;

            // Store: face index = j * ni + i (face between cell i and i+1)
            let face_idx = (j * ni + i) * NVAR + v;
            left_xi[face_idx] = wL;
            right_xi[face_idx] = wR;
        }
    }
}

// Reconstruct along eta (j-direction) for a fixed i-column.
// Face (i, j+1/2) has left state from cell (i,j) and right state from cell (i,j+1).
@compute @workgroup_size(1, 1, 1)
fn reconstruct_eta(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni) {
        return;
    }

    // Process each face along this i-column
    for (var j: u32 = 0u; j < nj - 1u; j = j + 1u) {
        for (var v: u32 = 0u; v < NVAR; v = v + 1u) {
            let wjm1 = get_w(i32(i), i32(j) - 1, v, ni, nj);
            let wj   = get_w(i32(i), i32(j),     v, ni, nj);
            let wjp1 = get_w(i32(i), i32(j) + 1, v, ni, nj);
            let wjp2 = get_w(i32(i), i32(j) + 2, v, ni, nj);

            let dL = wj - wjm1;
            let dR = wjp1 - wj;
            let dL2 = wjp1 - wj;
            let dR2 = wjp2 - wjp1;

            let slope_L = minmod(dL, dR);
            let slope_R = minmod(dL2, dR2);

            let wL = wj + 0.5 * slope_L;
            let wR = wjp1 - 0.5 * slope_R;

            // Face index: j * ni + i (face between cell j and j+1 at column i)
            let face_idx = (j * ni + i) * NVAR + v;
            left_eta[face_idx] = wL;
            right_eta[face_idx] = wR;
        }
    }
}
