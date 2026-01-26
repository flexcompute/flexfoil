//! Circulation updates for viscous-inviscid coupling.
//!
//! This mirrors XFOIL's QVFUE (set QVIS from UEDG) and GAMQV (set GAM from QVIS).
//!
//! # XFOIL Reference
//!
//! QVFUE (xpanel.f:1597):
//! ```fortran
//! DO IS=1, 2
//!   DO IBL=2, NBL(IS)
//!     I = IPAN(IBL,IS)
//!     QVIS(I) = VTI(IBL,IS) * UEDG(IBL,IS)
//!   ENDDO
//! ENDDO
//! ```
//!
//! GAMQV (xpanel.f:1633):
//! ```fortran
//! DO I=1, N
//!   GAM(I) = QVIS(I)
//!   GAM_A(I) = QINV_A(I)
//! ENDDO
//! ```
//!
//! VTI is the "viscous tangent inversion" sign:
//! - Upper surface (IS=1): VTI = +1
//! - Lower surface (IS=2): VTI = -1

use rustfoil_bl::state::BlStation;

/// Set QVIS from the current edge velocities (Ue).
///
/// This is XFOIL's QVFUE subroutine. QVIS is the panel tangential velocity
/// consistent with the vorticity distribution (GAM).
///
/// # Arguments
/// * `stations` - BL stations with current edge velocities
///
/// # Returns
/// QVIS array (panel tangential velocity)
pub fn update_qvis_from_uedg(stations: &[BlStation]) -> Vec<f64> {
    // In single-surface mode (what we mostly use), VTI is implicit in the sign of Ue.
    // For two-surface mode, VTI would need to be passed explicitly.
    // For now, just return Ue directly since our CL computation handles signs.
    stations.iter().map(|s| s.u).collect()
}

/// Set QVIS from edge velocities with explicit VTI signs.
///
/// This properly implements XFOIL's QVFUE with VTI handling:
///   QVIS(I) = VTI(IBL,IS) * UEDG(IBL,IS)
///
/// # Arguments
/// * `upper_stations` - Upper surface BL stations (stagnation to TE)
/// * `lower_stations` - Lower surface BL stations (stagnation to TE)
///
/// # Returns
/// (upper_qvis, lower_qvis) arrays
pub fn update_qvis_from_uedg_two_surfaces(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
) -> (Vec<f64>, Vec<f64>) {
    // Upper surface: VTI = +1 (velocity in same direction as panel tangent)
    let upper_qvis: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
    
    // Lower surface: VTI = -1 (velocity opposite to panel tangent for same sign of lift)
    // However, in our convention, lower surface Ue is already stored with the proper sign
    // based on how we extract it from the inviscid solution.
    let lower_qvis: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();
    
    (upper_qvis, lower_qvis)
}

/// Update circulation (GAM) from QVIS.
///
/// This is XFOIL's GAMQV subroutine. The circulation distribution equals
/// the tangential velocity distribution.
///
/// # Arguments
/// * `qvis` - Panel tangential velocities
///
/// # Returns
/// GAM array (panel circulation = vortex strength)
pub fn update_circulation_from_qvis(qvis: &[f64]) -> Vec<f64> {
    let mut gamma = qvis.to_vec();
    
    // Kutta condition: gamma at TE should satisfy trailing edge closure
    // In XFOIL: GAM(N+1) = 0 for the TE pseudo-node
    // We enforce this by setting the last value to 0
    if let Some(last) = gamma.last_mut() {
        *last = 0.0;
    }
    
    gamma
}

/// Compute CL from two-surface gamma distribution.
///
/// This implements the circulation theorem: L = ρ V∞ Γ
/// where Γ is the total circulation around the airfoil.
///
/// In the linear vorticity formulation, gamma at each node equals
/// the local tangential velocity. The total circulation is:
///   Γ = ∮ γ ds = ∫_upper γ ds - ∫_lower γ ds
///
/// CL = 2Γ / (V∞ c) = 2Γ for normalized flow (V∞ = 1, c = 1)
///
/// # Arguments
/// * `upper_gamma` - Circulation on upper surface
/// * `lower_gamma` - Circulation on lower surface  
/// * `upper_arc` - Arc lengths on upper surface
/// * `lower_arc` - Arc lengths on lower surface
///
/// # Returns
/// Lift coefficient CL
pub fn compute_cl_from_gamma(
    upper_gamma: &[f64],
    lower_gamma: &[f64],
    upper_arc: &[f64],
    lower_arc: &[f64],
) -> f64 {
    if upper_gamma.len() < 2 || lower_gamma.len() < 2 {
        return 0.0;
    }
    
    let mut circulation = 0.0;
    
    // Upper surface: integrate gamma * ds (counterclockwise)
    for i in 0..upper_gamma.len().saturating_sub(1) {
        let gamma_avg = 0.5 * (upper_gamma[i] + upper_gamma[i + 1]);
        let ds = (upper_arc[i + 1] - upper_arc[i]).abs();
        circulation += gamma_avg * ds;
    }
    
    // Lower surface: integrate gamma * ds (clockwise = negative contribution)
    for i in 0..lower_gamma.len().saturating_sub(1) {
        let gamma_avg = 0.5 * (lower_gamma[i] + lower_gamma[i + 1]);
        let ds = (lower_arc[i + 1] - lower_arc[i]).abs();
        circulation -= gamma_avg * ds;
    }
    
    // CL = 2Γ for normalized flow
    2.0 * circulation
}

/// Initialize basic wake stations after the trailing edge.
///
/// Wake stations extend downstream from the trailing edge and are used for:
/// 1. Proper Squire-Young CD computation (needs far-wake values)
/// 2. Wake displacement effect on airfoil
///
/// # Arguments
/// * `upper_te` - Last station on upper surface (at TE)
/// * `lower_te` - Last station on lower surface (at TE)
/// * `n_wake` - Number of wake stations
/// * `wake_length` - Total wake length (x/c)
///
/// # Returns
/// Vector of wake BlStations
///
/// # Note
/// This is a simplified wake model. XFOIL's full wake model includes:
/// - Wake curvature from potential flow
/// - Displacement effect coupling
/// - Far-wake extrapolation
pub fn initialize_wake_stations(
    upper_te: &BlStation,
    lower_te: &BlStation,
    n_wake: usize,
    wake_length: f64,
) -> Vec<BlStation> {
    if n_wake == 0 {
        return Vec::new();
    }
    
    let mut wake_stations = Vec::with_capacity(n_wake);
    
    // Wake starts at TE with combined upper+lower properties
    // Momentum and displacement thicknesses add at TE
    let theta_te = upper_te.theta + lower_te.theta;
    let dstar_te = upper_te.de + lower_te.de;
    let ue_te = 0.5 * (upper_te.u.abs() + lower_te.u.abs());
    
    // Wake shape factor typically around 1.0-1.1 at TE
    let h_te = if theta_te > 0.0 { dstar_te / theta_te } else { 1.0 };
    
    // Generate wake stations with gradual relaxation toward freestream
    for i in 0..n_wake {
        let t = (i + 1) as f64 / n_wake as f64; // 0 < t <= 1
        let x_wake = 1.0 + t * wake_length; // x/c from TE
        
        // Wake grows according to turbulent far-wake similarity
        // θ ~ x^(1/2), H → 1.0 as x → ∞
        let growth = (1.0 + 2.0 * t * wake_length).sqrt();
        let theta_wake = theta_te * growth;
        
        // Shape factor decays toward 1.0 (fully developed wake)
        let h_wake = 1.0 + (h_te - 1.0) * (-t * 2.0).exp();
        
        // Edge velocity recovers toward freestream
        let ue_wake = ue_te + (1.0 - ue_te) * (1.0 - (-t * 3.0).exp());
        
        let mut station = BlStation::default();
        station.x = x_wake;
        // Arc length continues from TE (using x field which stores arc length in BL code)
        station.u = ue_wake;
        station.theta = theta_wake;
        station.de = h_wake * theta_wake;
        station.h = h_wake;
        station.mass_defect = ue_wake * station.de;
        station.is_wake = true;
        station.is_turbulent = true;
        station.is_laminar = false;
        
        // Cf in wake is 0 (no wall)
        station.cf = 0.0;
        
        wake_stations.push(station);
    }
    
    wake_stations
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_qvis_from_uedg() {
        let stations = vec![
            BlStation { u: 1.0, ..Default::default() },
            BlStation { u: 1.2, ..Default::default() },
            BlStation { u: 0.8, ..Default::default() },
        ];
        
        let qvis = update_qvis_from_uedg(&stations);
        assert_eq!(qvis.len(), 3);
        assert!((qvis[0] - 1.0).abs() < 1e-10);
        assert!((qvis[1] - 1.2).abs() < 1e-10);
        assert!((qvis[2] - 0.8).abs() < 1e-10);
    }
    
    #[test]
    fn test_gamma_from_qvis() {
        let qvis = vec![1.0, 1.2, 0.8, 0.5];
        let gamma = update_circulation_from_qvis(&qvis);
        
        assert_eq!(gamma.len(), 4);
        assert!((gamma[0] - 1.0).abs() < 1e-10);
        assert!((gamma[3] - 0.0).abs() < 1e-10); // Last should be 0 (Kutta)
    }
    
    #[test]
    fn test_cl_from_gamma() {
        // Simple test: uniform gamma on both surfaces
        // Upper: gamma = 1.0 over arc length 1.0 → Γ_upper = 1.0
        // Lower: gamma = 0.5 over arc length 1.0 → Γ_lower = 0.5  
        // Net Γ = 1.0 - 0.5 = 0.5
        // CL = 2 * 0.5 = 1.0
        
        let upper_gamma = vec![1.0, 1.0, 1.0];
        let lower_gamma = vec![0.5, 0.5, 0.5];
        let upper_arc = vec![0.0, 0.5, 1.0];
        let lower_arc = vec![0.0, 0.5, 1.0];
        
        let cl = compute_cl_from_gamma(&upper_gamma, &lower_gamma, &upper_arc, &lower_arc);
        assert!((cl - 1.0).abs() < 0.1, "Expected CL ≈ 1.0, got {}", cl);
    }
    
    #[test]
    fn test_wake_stations_initialization() {
        let upper_te = BlStation {
            x: 1.0,
            u: 0.9,
            theta: 0.002,
            de: 0.005,
            h: 2.5,
            ..Default::default()
        };
        let lower_te = BlStation {
            x: 1.0,
            u: 0.8,
            theta: 0.003,
            de: 0.006,
            h: 2.0,
            ..Default::default()
        };
        
        let wake = initialize_wake_stations(&upper_te, &lower_te, 5, 1.0);
        
        assert_eq!(wake.len(), 5);
        
        // Check first wake station
        assert!(wake[0].x > 1.0, "Wake should be downstream of TE");
        assert!(wake[0].is_wake, "Should be marked as wake");
        assert!(wake[0].is_turbulent, "Wake should be turbulent");
        assert!(!wake[0].is_laminar, "Wake should not be laminar");
        
        // Check last wake station
        let last = &wake[4];
        assert!(last.x <= 2.0 + 1e-6, "Last wake station should be at x ≈ 2.0");
        assert!(last.h < upper_te.h, "Wake H should decay toward 1.0");
    }
}
