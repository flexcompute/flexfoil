//! Initial condition computation for the CFD solver.
//!
//! Generates freestream uniform initial conditions for the conservative
//! variables Q = [rho, rho*u, rho*v, E, nu_tilde].

use crate::config::CfdConfig;

/// Compute freestream initial conditions for the entire grid.
///
/// Returns a flat f32 array of length ni*nj*5 containing:
/// Q = [rho, rho*u, rho*v, E, nu_tilde] at each cell.
///
/// Uses non-dimensionalization where:
/// - rho_inf = 1.0
/// - |V_inf| = M_inf * a_inf = M_inf (with a_inf = 1.0)
/// - p_inf = 1.0 / gamma
pub fn compute_initial_conditions(config: &CfdConfig) -> Vec<f32> {
    let ni = config.ni as usize;
    let nj = config.nj as usize;
    let gamma = config.gamma;
    let mach = config.mach_inf;
    let alpha = config.alpha; // radians

    // Freestream state (non-dimensional)
    let rho_inf: f32 = 1.0;
    let u_inf: f32 = mach * alpha.cos();
    let v_inf: f32 = mach * alpha.sin();
    let p_inf: f32 = 1.0 / gamma;
    let e_inf: f32 = p_inf / (gamma - 1.0) + 0.5 * rho_inf * (u_inf * u_inf + v_inf * v_inf);

    // SA turbulent viscosity initial value
    let nu_tilde_inf: f32 = if config.physics as u32 == 2 {
        // RANS: initialize to ~3x molecular viscosity
        3.0 * mach / config.reynolds
    } else {
        0.0
    };

    let n_total = ni * nj * 5;
    let mut q = vec![0.0f32; n_total];

    for j in 0..nj {
        for i in 0..ni {
            let base = (j * ni + i) * 5;
            q[base] = rho_inf;
            q[base + 1] = rho_inf * u_inf;
            q[base + 2] = rho_inf * v_inf;
            q[base + 3] = e_inf;
            q[base + 4] = nu_tilde_inf;
        }
    }

    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initial_conditions_euler() {
        let config = CfdConfig {
            ni: 4,
            nj: 4,
            mach_inf: 0.5,
            alpha: 0.0,
            gamma: 1.4,
            ..CfdConfig::default()
        };

        let q = compute_initial_conditions(&config);
        assert_eq!(q.len(), 4 * 4 * 5);

        // Check first cell
        let rho = q[0];
        let rhou = q[1];
        let rhov = q[2];
        let e = q[3];
        let nu = q[4];

        assert!((rho - 1.0).abs() < 1e-6);
        assert!((rhou - 0.5).abs() < 1e-6); // rho * M * cos(0)
        assert!(rhov.abs() < 1e-6); // rho * M * sin(0)
        assert!(nu.abs() < 1e-6); // Euler: no SA
        // E = p/(gamma-1) + 0.5*rho*V^2 = (1/1.4)/0.4 + 0.5*0.25 = 1.7857 + 0.125 = 1.9107
        assert!((e - 1.9107).abs() < 0.01);
    }
}
