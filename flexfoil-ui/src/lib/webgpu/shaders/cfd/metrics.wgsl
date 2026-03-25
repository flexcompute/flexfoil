// Compute grid metrics from mesh node coordinates.
// Runs once at mesh setup (not per timestep).
//
// Computes: xi_x, xi_y, eta_x, eta_y, J (Jacobian)
// Using central differences on the (x,y) coordinates.

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read> mesh_x: array<f32>;
@group(0) @binding(2) var<storage, read> mesh_y: array<f32>;
@group(0) @binding(3) var<storage, read_write> metrics: array<f32>;
// metrics layout: [xi_x, xi_y, eta_x, eta_y, J] per cell = 5 floats

@compute @workgroup_size(16, 16, 1)
fn compute_metrics(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni || j >= nj) {
        return;
    }

    // Compute x_xi, y_xi (derivatives w.r.t. xi = i-direction)
    let ip = wrap_i(i32(i) + 1, ni);
    let im = wrap_i(i32(i) - 1, ni);

    let x_xi = 0.5 * (mesh_x[cell_idx(ip, j, ni)] - mesh_x[cell_idx(im, j, ni)]);
    let y_xi = 0.5 * (mesh_y[cell_idx(ip, j, ni)] - mesh_y[cell_idx(im, j, ni)]);

    // Compute x_eta, y_eta (derivatives w.r.t. eta = j-direction)
    var x_eta: f32;
    var y_eta: f32;
    if (j == 0u) {
        x_eta = mesh_x[cell_idx(i, 1u, ni)] - mesh_x[cell_idx(i, 0u, ni)];
        y_eta = mesh_y[cell_idx(i, 1u, ni)] - mesh_y[cell_idx(i, 0u, ni)];
    } else if (j == nj - 1u) {
        x_eta = mesh_x[cell_idx(i, nj - 1u, ni)] - mesh_x[cell_idx(i, nj - 2u, ni)];
        y_eta = mesh_y[cell_idx(i, nj - 1u, ni)] - mesh_y[cell_idx(i, nj - 2u, ni)];
    } else {
        x_eta = 0.5 * (mesh_x[cell_idx(i, j + 1u, ni)] - mesh_x[cell_idx(i, j - 1u, ni)]);
        y_eta = 0.5 * (mesh_y[cell_idx(i, j + 1u, ni)] - mesh_y[cell_idx(i, j - 1u, ni)]);
    }

    // Jacobian of the transformation: J = x_xi * y_eta - x_eta * y_xi
    let J = x_xi * y_eta - x_eta * y_xi;
    let inv_J = select(1.0 / J, 1e10, abs(J) < 1e-15);

    // Metric terms: inverse mapping
    // xi_x = y_eta / J, xi_y = -x_eta / J
    // eta_x = -y_xi / J, eta_y = x_xi / J
    let base = cell_idx(i, j, ni) * 5u;
    metrics[base + 0u] = y_eta * inv_J;      // xi_x
    metrics[base + 1u] = -x_eta * inv_J;     // xi_y
    metrics[base + 2u] = -y_xi * inv_J;      // eta_x
    metrics[base + 3u] = x_xi * inv_J;       // eta_y
    metrics[base + 4u] = J;                   // Jacobian
}
