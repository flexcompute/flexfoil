//! XFOIL-style Cτ lag-dissipation turbulent boundary layer model.
//!
//! This implements Mark Drela's two-equation turbulent BL formulation
//! used in XFOIL and ISES codes. The primary variables are:
//! - θ (momentum thickness)
//! - Cτ (shear stress coefficient with lag)
//!
//! The key advantage over Head's method is the lag equation which
//! captures history effects better, especially in separating flows.
//!
//! Reference: Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds
//! Number Airfoils", MIT, 1989.

use super::state::BLState;

/// XFOIL turbulent boundary layer constants.
///
/// These are the default values from XFOIL's BLPINI subroutine.
#[derive(Debug, Clone, Copy)]
pub struct XfoilConstants {
    /// Shear coefficient lag constant (SCCON)
    pub sccon: f64,
    /// G-β equilibrium constant A (GACON)
    pub gacon: f64,
    /// G-β equilibrium constant B (GBCON)
    pub gbcon: f64,
    /// G-β wall term constant (GCCON)
    pub gccon: f64,
    /// Wall/wake dissipation length ratio (DLCON)
    pub dlcon: f64,
    /// Cf scaling factor (CFFAC)
    pub cffac: f64,
}

impl Default for XfoilConstants {
    fn default() -> Self {
        // Values from XFOIL's BLPINI
        Self {
            sccon: 5.6,
            gacon: 6.70,
            gbcon: 0.75,
            gccon: 18.0,
            dlcon: 0.9,
            cffac: 1.0,
        }
    }
}

/// Compute turbulent skin friction using Coles' wall law (XFOIL CFT).
///
/// This is the Cf correlation used in XFOIL, based on the Coles
/// velocity defect law with compressibility corrections.
///
/// # Arguments
/// * `hk` - Kinematic shape factor H
/// * `re_theta` - Reynolds number based on θ
/// * `msq` - Mach number squared (0 for incompressible)
/// * `cffac` - Cf scaling factor (default 1.0)
///
/// # Returns
/// (Cf, dCf/dH, dCf/dReθ, dCf/dM²)
pub fn coles_cf(hk: f64, re_theta: f64, msq: f64, cffac: f64) -> (f64, f64, f64, f64) {
    const GAM: f64 = 1.4;
    const LN10: f64 = 2.302585093;  // ln(10)
    let gm1 = GAM - 1.0;
    
    // Compressibility factor
    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();
    
    // Log Reynolds number: GRT = ln(Reθ/fc), then GRT/ln(10) = log10(Reθ/fc)
    let grt = (re_theta / fc).max(10.0).ln().max(3.0);
    
    // Exponent
    let gex = -1.74 - 0.31 * hk;
    
    // Exponential term (from -1.33*H)
    let arg = (-1.33 * hk).max(-20.0);
    
    // Tanh term for low H behavior (kicks in when H > 3.5)
    let thk = (4.0 - hk / 0.875).tanh();
    
    // Base Cf: XFOIL formula
    // CFO = CFFAC * 0.3 * exp(-1.33*H) * (log10(Reθ))^(-1.74 - 0.31*H)
    let log10_ret = grt / LN10;
    let cfo = cffac * 0.3 * arg.exp() * log10_ret.powf(gex);
    
    // Total Cf with low-H correction and compressibility
    let cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
    
    // Derivatives (simplified)
    let cf_hk = (-1.33 * cfo - 0.31 * log10_ret.ln() * cfo
                 - 1.1e-4 * (1.0 - thk * thk) / 0.875) / fc;
    let cf_rt = gex * cfo / (fc * grt * re_theta);
    let cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / (fc * fc)) 
                 - 0.25 * gm1 * cf / (fc * fc);
    
    (cf.max(1e-6), cf_hk, cf_rt, cf_msq)
}

/// Compute equilibrium Cτ from the G-β locus.
///
/// The G-β relationship defines the equilibrium shear stress
/// as a function of the pressure gradient parameter β.
///
/// # Arguments
/// * `h` - Shape factor H
/// * `hs` - Kinematic energy shape factor H*
/// * `re_theta` - Reynolds number based on θ
/// * `cf` - Skin friction coefficient
/// * `beta` - Pressure gradient parameter
/// * `consts` - XFOIL constants
///
/// # Returns
/// Equilibrium Cτ value
pub fn ctau_equilibrium(
    h: f64,
    hs: f64,
    re_theta: f64,
    cf: f64,
    beta: f64,
    consts: &XfoilConstants,
) -> f64 {
    // G-β equilibrium locus (XFOIL formulation)
    // G = GACON * sqrt(1 + GBCON * β) + GCCON / (H * Reθ * sqrt(Cf/2))
    
    let cf_half = (cf / 2.0).max(1e-10);
    let wall_term = consts.gccon / (h * re_theta * cf_half.sqrt()).max(1e-10);
    
    let g = if beta >= 0.0 {
        consts.gacon * (1.0 + consts.gbcon * beta).sqrt() + wall_term
    } else {
        // Favorable pressure gradient
        consts.gacon * (1.0 + consts.gbcon * beta.abs()).sqrt().recip() + wall_term
    };
    
    // Equilibrium Cτ from G
    // Cτ_eq = (Hs * Cf/2) * (G / 6.7)²
    let g_ratio = g / 6.7;
    let ctau_eq = hs * cf_half * g_ratio * g_ratio;
    
    ctau_eq.max(1e-10)
}

/// Compute the kinematic energy shape factor H* (Hs).
///
/// This relates to the energy thickness δ**.
///
/// # Arguments
/// * `hk` - Kinematic shape factor H
/// * `re_theta` - Reynolds number based on θ
///
/// # Returns
/// H* value
pub fn energy_shape_factor(hk: f64, re_theta: f64) -> f64 {
    // XFOIL correlation for H*
    // From Drela's thesis
    
    let hk_clamped = hk.clamp(1.0, 10.0);
    
    if hk_clamped < 4.0 {
        // Attached flow correlation
        let hk_ref = hk_clamped - 1.0;
        1.505 + 4.0 / (re_theta.max(100.0) + 1000.0) 
            + 0.165 * hk_ref.powi(2) - 0.035 * hk_ref.powi(3)
    } else {
        // Separated flow correlation  
        let hs4 = 1.505 + 4.0 / (re_theta.max(100.0) + 1000.0) + 0.165 * 9.0 - 0.035 * 27.0;
        hs4 + 0.04 * (hk_clamped - 4.0)
    }
}

/// Compute the lag length scale L.
///
/// This determines how quickly Cτ relaxes to equilibrium.
///
/// # Arguments
/// * `theta` - Momentum thickness
/// * `hs` - Energy shape factor H*
/// * `consts` - XFOIL constants
///
/// # Returns
/// Lag length scale
pub fn lag_length(theta: f64, hs: f64, consts: &XfoilConstants) -> f64 {
    // L = SCCON * θ * Hs
    consts.sccon * theta * hs
}

/// Solve XFOIL-style Cτ turbulent equations for one step.
///
/// Marches the boundary layer using the two-equation Cτ lag formulation.
///
/// # Arguments
/// * `state_prev` - BL state at previous station
/// * `s` - Current arc-length
/// * `ue` - Edge velocity at current station  
/// * `due_ds` - Velocity gradient dUe/ds
/// * `reynolds` - Chord Reynolds number
/// * `consts` - XFOIL constants
///
/// # Returns
/// Updated BL state at current station.
pub fn xfoil_turb_solve(
    state_prev: &BLState,
    s: f64,
    ue: f64,
    due_ds: f64,
    reynolds: f64,
    consts: &XfoilConstants,
) -> BLState {
    let mut state = state_prev.clone();
    state.s = s;
    state.ue = ue;
    state.is_turbulent = true;
    
    let ds = (s - state_prev.s).abs().max(1e-10);
    let ue_avg = 0.5 * (ue + state_prev.ue).max(1e-10);
    
    // Previous values
    let theta_prev = state_prev.theta.max(1e-10);
    let h_prev = state_prev.h.clamp(1.0, 10.0);
    let ctau_prev = state_prev.ctau.max(1e-6);
    
    // Compute Re_theta
    let re_theta = (ue_avg * theta_prev * reynolds).max(100.0);
    
    // Compute Cf using Coles correlation (incompressible)
    let (cf, _cf_h, _cf_rt, _) = coles_cf(h_prev, re_theta, 0.0, consts.cffac);
    let cf = cf.max(1e-6);
    
    // Energy shape factor
    let hs = energy_shape_factor(h_prev, re_theta);
    
    // Pressure gradient parameter β
    // β = (θ/Ue)(dUe/ds) × 2/(Cf × Hs)
    let beta = theta_prev * due_ds / ue_avg * 2.0 / (cf * hs).max(1e-10);
    
    // Equilibrium Cτ
    let ctau_eq = ctau_equilibrium(h_prev, hs, re_theta, cf, beta, consts);
    
    // Lag length
    let lag_l = lag_length(theta_prev, hs, consts);
    
    // === Solve the two equations ===
    
    // 1. Momentum equation: dθ/ds = Cf/2 - (H+2)(θ/Ue)(dUe/ds)
    let dtheta_ds = cf / 2.0 - (h_prev + 2.0) * theta_prev * due_ds / ue_avg;
    state.theta = (theta_prev + dtheta_ds * ds).max(1e-10);
    
    // 2. Shear lag equation: dCτ/ds = (Cτ_eq - Cτ) / L
    let dctau_ds = (ctau_eq - ctau_prev) / lag_l.max(1e-10);
    state.ctau = (ctau_prev + dctau_ds * ds).max(1e-6);
    
    // Derive H from Cτ using inverse of equilibrium relation
    // For simplicity, use an iterative approach or correlation
    // Here we use a simplified relation: Cτ ≈ Hs × Cf/2 at equilibrium
    // So Cf_new ≈ 2 × Cτ / Hs
    let cf_from_ctau = 2.0 * state.ctau / hs.max(0.1);
    
    // Invert Coles correlation to get H from Cf
    // Simplified: use Newton iteration or table lookup
    // For now, use a linearized update
    let h_change = (cf - cf_from_ctau) / (cf + 1e-10) * 0.5;
    state.h = (h_prev + h_change).clamp(1.0, 10.0);
    
    // Update displacement thickness
    state.delta_star = state.h * state.theta;
    
    // Update Cf with new H
    let re_theta_new = (ue * state.theta * reynolds).max(100.0);
    let (cf_new, _, _, _) = coles_cf(state.h, re_theta_new, 0.0, consts.cffac);
    state.cf = cf_new.max(1e-6);
    
    state
}

/// Initialize turbulent boundary layer at transition.
///
/// Sets initial conditions for turbulent marching after transition.
///
/// # Arguments
/// * `laminar_state` - State at end of laminar region
/// * `re_theta` - Reynolds number at transition
///
/// # Returns
/// Initial turbulent state
pub fn initialize_turbulent(laminar_state: &BLState, re_theta: f64) -> BLState {
    let mut state = laminar_state.clone();
    state.is_turbulent = true;
    
    // Initial Cτ from equilibrium at zero pressure gradient
    let h = state.h.clamp(1.0, 10.0);
    let hs = energy_shape_factor(h, re_theta);
    let (cf, _, _, _) = coles_cf(h, re_theta, 0.0, 1.0);
    
    // For β=0, G ≈ GACON, so Cτ_eq ≈ Hs × Cf/2 × (GACON/6.7)²
    let consts = XfoilConstants::default();
    state.ctau = hs * cf / 2.0 * (consts.gacon / 6.7).powi(2);
    state.ctau = state.ctau.max(0.03); // XFOIL default initial Cτ
    
    state.cf = cf;
    
    state
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_coles_cf() {
        // At H=1.4, Re_theta=1000, incompressible
        let (cf, _, _, _) = coles_cf(1.4, 1000.0, 0.0, 1.0);
        assert!(cf > 0.002 && cf < 0.006, "Cf={} should be ~0.003-0.004", cf);
    }
    
    #[test]
    fn test_energy_shape_factor() {
        // For attached flow, Hs ≈ 1.5-1.6
        let hs = energy_shape_factor(1.4, 1000.0);
        assert!(hs > 1.4 && hs < 1.8, "Hs={} should be ~1.5-1.6", hs);
    }
    
    #[test]
    fn test_ctau_equilibrium() {
        let consts = XfoilConstants::default();
        let ctau_eq = ctau_equilibrium(1.4, 1.5, 1000.0, 0.004, 0.0, &consts);
        assert!(ctau_eq > 0.001 && ctau_eq < 0.01, "Ctau_eq={} should be ~0.003", ctau_eq);
    }
    
    #[test]
    fn test_xfoil_constants_default() {
        let c = XfoilConstants::default();
        assert_eq!(c.sccon, 5.6);
        assert_eq!(c.gacon, 6.70);
    }
    
    #[test]
    fn test_turbulent_march() {
        let consts = XfoilConstants::default();
        
        // Start with typical post-transition state
        let mut state = BLState::default();
        state.theta = 0.001;
        state.h = 1.4;
        state.delta_star = 0.0014;
        state.ctau = 0.03;
        state.is_turbulent = true;
        state.s = 0.1;
        state.ue = 1.0;
        
        let reynolds = 1e6;
        
        // March 10 steps along flat plate
        for i in 1..10 {
            let s = 0.1 + 0.01 * i as f64;
            state = xfoil_turb_solve(&state, s, 1.0, 0.0, reynolds, &consts);
        }
        
        // Theta should have grown
        assert!(state.theta > 0.001);
        
        // H should be in reasonable range
        assert!(state.h > 1.2 && state.h < 2.0);
        
        // Cf should be positive and small
        assert!(state.cf > 0.001 && state.cf < 0.01);
    }
}
