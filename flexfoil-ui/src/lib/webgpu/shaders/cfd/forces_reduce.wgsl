// Parallel reduction to compute aerodynamic force coefficients.
// Integrates pressure and (optionally) shear stress along the wall (j=0).
//
// Outputs: [Cl, Cd, Cm] in the force_accum buffer.
// Also computes L2 residual norm in residual_norm buffer.

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> Q: array<f32>;
@group(0) @binding(2) var<storage, read> metrics: array<f32>;
@group(0) @binding(3) var<storage, read> mesh_x: array<f32>;
@group(0) @binding(4) var<storage, read> mesh_y: array<f32>;
@group(0) @binding(5) var<storage, read> residual: array<f32>;
@group(0) @binding(6) var<storage, read_write> force_accum: array<f32>;   // [Cl, Cd, Cm]
@group(0) @binding(7) var<storage, read_write> residual_norm: array<f32>; // [L2_rho, L2_rhou, L2_rhov, L2_E]

var<workgroup> shared_cl: array<f32, 256>;
var<workgroup> shared_cd: array<f32, 256>;
var<workgroup> shared_cm: array<f32, 256>;
var<workgroup> shared_res: array<f32, 256>;

@compute @workgroup_size(256, 1, 1)
fn reduce_wall_forces(@builtin(local_invocation_id) lid: vec3<u32>,
                       @builtin(workgroup_id) wid: vec3<u32>,
                       @builtin(num_workgroups) nwg: vec3<u32>) {
    let tid = lid.x;
    let ni = params.ni;
    let nj = params.nj;
    let gamma = params.gamma;
    let mach = params.mach_inf;
    let alpha = params.alpha;

    // Touch all bindings so auto-layout includes them
    let _m_touch = metrics[0];
    _ = _m_touch;

    // Reference quantities for force coefficients
    let p_inf = 1.0 / gamma;
    let q_inf = 0.5 * gamma * p_inf * mach * mach; // Dynamic pressure

    // Each thread processes multiple wall cells
    var local_cl: f32 = 0.0;
    var local_cd: f32 = 0.0;
    var local_cm: f32 = 0.0;
    var local_res_sq: f32 = 0.0;

    var idx = tid;
    while (idx < ni) {
        // Surface pressure at wall (j=0)
        let base = q_idx(idx, 0u, 0u, ni);
        let rho = Q[base + 0u];
        let rhou = Q[base + 1u];
        let rhov = Q[base + 2u];
        let E = Q[base + 3u];
        let u = rhou / max(rho, 1e-10);
        let v = rhov / max(rho, 1e-10);
        let p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
        let cp = (p - p_inf) / q_inf;

        // Surface element: ds = |dx/di| at wall
        let ip = select(idx + 1u, 0u, idx + 1u >= ni);
        let im = select(idx - 1u, ni - 1u, idx == 0u);
        let dx = 0.5 * (mesh_x[cell_idx(ip, 0u, ni)] - mesh_x[cell_idx(im, 0u, ni)]);
        let dy = 0.5 * (mesh_y[cell_idx(ip, 0u, ni)] - mesh_y[cell_idx(im, 0u, ni)]);

        // Outward normal (rotate tangent 90 degrees: tangent (dx,dy) -> normal (-dy, dx))
        // Convention: positive Cl is lift (upward)
        let fx = -p * (-dy); // Force in x per unit span (pressure * normal_x * ds)
        let fy = -p * dx;    // Force in y

        // Rotate to wind axes
        let ca = cos(alpha);
        let sa = sin(alpha);
        let drag_contrib = fx * ca + fy * sa;
        let lift_contrib = -fx * sa + fy * ca;

        // Moment about (0.25, 0) - quarter chord
        let xc = mesh_x[cell_idx(idx, 0u, ni)] - 0.25;
        let yc = mesh_y[cell_idx(idx, 0u, ni)];
        let moment_contrib = xc * fy - yc * fx; // Pitching moment

        local_cd += drag_contrib / q_inf;
        local_cl += lift_contrib / q_inf;
        local_cm += moment_contrib / q_inf;

        // Accumulate residual L2 norm (all interior cells in this thread's range)
        for (var j: u32 = 1u; j < nj - 1u; j = j + 1u) {
            let r0 = residual[q_idx(idx, j, 0u, ni)];
            local_res_sq += r0 * r0;
        }

        idx += 256u;
    }

    shared_cl[tid] = local_cl;
    shared_cd[tid] = local_cd;
    shared_cm[tid] = local_cm;
    shared_res[tid] = local_res_sq;
    workgroupBarrier();

    // Parallel reduction
    var stride: u32 = 128u;
    while (stride > 0u) {
        if (tid < stride) {
            shared_cl[tid] += shared_cl[tid + stride];
            shared_cd[tid] += shared_cd[tid + stride];
            shared_cm[tid] += shared_cm[tid + stride];
            shared_res[tid] += shared_res[tid + stride];
        }
        workgroupBarrier();
        stride = stride >> 1u;
    }

    if (tid == 0u) {
        force_accum[0] = shared_cl[0];
        force_accum[1] = shared_cd[0];
        force_accum[2] = shared_cm[0];
        residual_norm[0] = sqrt(shared_res[0] / f32(ni * (nj - 2u)));
    }
}
