// Boundary condition application for the CFD solver.
//
// Handles:
// - Wall BC (j=0): slip wall (Euler) or no-slip (NS/RANS)
// - Far-field BC (j=nj-1): characteristic-based
// - Periodicity in i-direction is implicit in the O-grid topology

@group(0) @binding(0) var<uniform> params: CfdParams;
@group(0) @binding(1) var<storage, read_write> Q: array<f32>;
@group(0) @binding(2) var<storage, read> metrics: array<f32>;
@group(0) @binding(3) var<storage, read> bc_type: array<u32>;
@group(0) @binding(4) var<storage, read> mesh_x: array<f32>;
@group(0) @binding(5) var<storage, read> mesh_y: array<f32>;

@compute @workgroup_size(256, 1, 1)
fn apply_boundary_conditions(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let ni = params.ni;
    let nj = params.nj;

    if (i >= ni) {
        return;
    }

    let gamma = params.gamma;

    // === Wall boundary (j=0) ===
    {
        let j = 0u;
        let cell = cell_idx(i, j, ni);

        // Get interior cell values (j=1)
        let base_int = q_idx(i, 1u, 0u, ni);
        let rho_int = Q[base_int + 0u];
        let rhou_int = Q[base_int + 1u];
        let rhov_int = Q[base_int + 2u];
        let E_int = Q[base_int + 3u];
        let nu_int = Q[base_int + 4u];

        // Surface normal (eta-direction at wall, pointing outward from body)
        let m_base = cell * 5u;
        let eta_x = metrics[m_base + 2u];
        let eta_y = metrics[m_base + 3u];
        let mag = sqrt(eta_x * eta_x + eta_y * eta_y);
        let nx = select(eta_x / mag, 0.0, mag < 1e-15);
        let ny = select(eta_y / mag, 0.0, mag < 1e-15);

        let u_int = rhou_int / rho_int;
        let v_int = rhov_int / rho_int;

        // Reflect normal velocity component
        let Vn = u_int * nx + v_int * ny;

        var u_wall: f32;
        var v_wall: f32;

        if (params.physics_mode == 0u) {
            // Euler: slip wall (zero normal velocity, keep tangential)
            u_wall = u_int - 2.0 * Vn * nx;
            v_wall = v_int - 2.0 * Vn * ny;
        } else {
            // NS/RANS: no-slip wall (zero velocity)
            u_wall = -u_int; // Ghost cell for zero at wall
            v_wall = -v_int;
        }

        let base_wall = q_idx(i, 0u, 0u, ni);
        Q[base_wall + 0u] = rho_int;  // Extrapolate density
        Q[base_wall + 1u] = rho_int * u_wall;
        Q[base_wall + 2u] = rho_int * v_wall;
        // Extrapolate pressure, set energy accordingly
        let p_int = (gamma - 1.0) * (E_int - 0.5 * rho_int * (u_int * u_int + v_int * v_int));
        Q[base_wall + 3u] = p_int / (gamma - 1.0) + 0.5 * rho_int * (u_wall * u_wall + v_wall * v_wall);
        // SA: zero at wall
        Q[base_wall + 4u] = select(nu_int, 0.0, params.physics_mode == 2u);
    }

    // === Far-field boundary (j=nj-1) ===
    {
        let j = nj - 1u;

        // Freestream state
        let mach = params.mach_inf;
        let alpha = params.alpha;
        let rho_inf: f32 = 1.0;
        let u_inf = mach * cos(alpha);
        let v_inf = mach * sin(alpha);
        let p_inf = 1.0 / gamma;
        let E_inf = p_inf / (gamma - 1.0) + 0.5 * rho_inf * (u_inf * u_inf + v_inf * v_inf);

        // Simple characteristic-based: use freestream for inflow, extrapolate for outflow
        let base_int = q_idx(i, j - 1u, 0u, ni);
        let rho_int = Q[base_int + 0u];
        let u_int = Q[base_int + 1u] / rho_int;
        let v_int = Q[base_int + 2u] / rho_int;

        // Radial direction at far-field
        let m_base = cell_idx(i, j, ni) * 5u;
        let eta_x = metrics[m_base + 2u];
        let eta_y = metrics[m_base + 3u];
        let mag = sqrt(eta_x * eta_x + eta_y * eta_y);
        let nx = select(eta_x / mag, 0.0, mag < 1e-15);
        let ny = select(eta_y / mag, 0.0, mag < 1e-15);

        // Check if inflow or outflow based on interior velocity
        let Vn = u_int * nx + v_int * ny;

        let base_ff = q_idx(i, j, 0u, ni);
        if (Vn >= 0.0) {
            // Outflow: extrapolate from interior
            Q[base_ff + 0u] = Q[base_int + 0u];
            Q[base_ff + 1u] = Q[base_int + 1u];
            Q[base_ff + 2u] = Q[base_int + 2u];
            Q[base_ff + 3u] = Q[base_int + 3u];
            Q[base_ff + 4u] = Q[base_int + 4u];
        } else {
            // Inflow: set to freestream
            Q[base_ff + 0u] = rho_inf;
            Q[base_ff + 1u] = rho_inf * u_inf;
            Q[base_ff + 2u] = rho_inf * v_inf;
            Q[base_ff + 3u] = E_inf;
            Q[base_ff + 4u] = 0.0;
        }
    }
}
