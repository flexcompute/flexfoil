//! Head's entrainment method for turbulent boundary layers.
//!
//! Head's method solves the momentum and entrainment equations
//! for turbulent boundary layers:
//!
//! Momentum: dθ/ds + (H+2)(θ/Ue)(dUe/ds) = Cf/2
//! Entrainment: d(H1*θ)/ds = Ce
//!
//! where H1 = (δ - δ*)/θ is the entrainment shape factor.
//!
//! Reference: Head, M.R. (1958) "Entrainment in the Turbulent Boundary Layer"

use super::state::BLState;
use super::closure;

/// Solve Head's entrainment equations for one step.
///
/// # Arguments
/// * `state_prev` - BL state at previous station
/// * `s` - Current arc-length
/// * `ue` - Edge velocity at current station
/// * `due_ds` - Velocity gradient dUe/ds
/// * `reynolds` - Chord Reynolds number
///
/// # Returns
/// Updated BL state at current station.
pub fn head_solve(
    state_prev: &BLState,
    s: f64,
    ue: f64,
    due_ds: f64,
    reynolds: f64,
) -> BLState {
    let mut state = state_prev.clone();
    state.s = s;
    state.ue = ue;
    state.is_turbulent = true;

    let ds = (s - state_prev.s).abs().max(1e-10);
    let ue_avg = 0.5 * (ue + state_prev.ue).max(1e-10);

    // Previous values
    let theta_prev = state_prev.theta.max(1e-10);
    let h_prev = state_prev.h.max(1.1);

    // Compute skin friction using Ludwieg-Tillmann formula
    let re_theta = ue_avg * theta_prev * reynolds;
    let cf = closure::ludwieg_tillmann_cf(re_theta, h_prev);

    // Head's shape factor H1 = (delta - delta*) / theta
    let h1_prev = closure::head_h1(h_prev);

    // Entrainment coefficient Ce = F(H1)
    let ce = closure::head_entrainment(h1_prev);

    // Momentum equation: dθ/ds = Cf/2 - (H+2)(θ/Ue)(dUe/ds)
    let dtheta_ds = cf / 2.0 - (h_prev + 2.0) * theta_prev * due_ds / ue_avg;

    // Update momentum thickness
    state.theta = (theta_prev + dtheta_ds * ds).max(1e-10);

    // Entrainment equation: d(H1*θ)/ds = Ce
    // Expanding: H1 * dθ/ds + θ * dH1/ds = Ce
    // => dH1/ds = (Ce - H1 * dθ/ds) / θ
    let dh1_ds = (ce - h1_prev * dtheta_ds) / theta_prev;
    let h1_new = (h1_prev + dh1_ds * ds).clamp(2.0, 20.0);

    // Convert H1 back to H
    state.h = closure::head_h_from_h1(h1_new);
    state.h = state.h.clamp(1.0, 10.0); // Prevent unrealistic values

    // Update displacement thickness
    state.delta_star = state.h * state.theta;

    // Update skin friction with new values
    let re_theta_new = ue * state.theta * reynolds;
    state.cf = closure::ludwieg_tillmann_cf(re_theta_new, state.h);

    // Ctau for turbulence models (simplified)
    state.ctau = 0.024 * state.cf;

    state
}

/// Green's lag-entrainment method (more advanced alternative).
///
/// This adds a lag equation for the entrainment coefficient to better
/// capture history effects.
#[allow(dead_code)]
pub fn green_lag_solve(
    state_prev: &BLState,
    s: f64,
    ue: f64,
    due_ds: f64,
    reynolds: f64,
    _ctau_eq: f64, // Equilibrium Ctau
) -> BLState {
    // For now, fall back to Head's method
    // Full Green's lag implementation would solve:
    // d(Ctau)/ds = (Ctau_eq - Ctau) / L
    // where L is a lag length scale
    
    head_solve(state_prev, s, ue, due_ds, reynolds)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_turbulent_flat_plate() {
        // Turbulent flat plate: θ/x ≈ 0.036 * Rex^(-1/5)
        // Start with an initial turbulent state
        let mut state = BLState::default();
        state.theta = 0.001;
        state.h = 1.4;
        state.delta_star = 1.4 * 0.001;
        state.is_turbulent = true;
        state.s = 0.1;
        state.ue = 1.0;

        let reynolds = 1e6;

        // March 10 steps
        for i in 1..10 {
            let s = 0.1 + 0.01 * i as f64;
            state = head_solve(&state, s, 1.0, 0.0, reynolds);
        }

        // Theta should have grown
        assert!(state.theta > 0.001);
        
        // H should be in turbulent range (1.3 - 2.0)
        assert!(state.h > 1.2 && state.h < 2.5);
        
        // Cf should be positive and small
        assert!(state.cf > 0.0 && state.cf < 0.01);
    }
}
