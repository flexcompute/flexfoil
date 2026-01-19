//! Boundary layer state variables.
//!
//! Defines the primary variables used in integral boundary layer methods:
//! - θ (theta): momentum thickness
//! - δ* (delta_star): displacement thickness  
//! - H: shape factor = δ*/θ
//! - Cf: skin friction coefficient
//! - Ctau: maximum shear stress coefficient (for turbulent closure)

// =============================================================================
// XFOIL Separation Thresholds (from xbl.f lines 556-557)
// =============================================================================

/// Maximum kinematic shape factor for laminar boundary layers.
/// When Hk exceeds this value, the flow is considered separated and
/// the solver switches to inverse mode.
/// From XFOIL: HLMAX = 3.8
pub const HK_MAX_LAMINAR: f64 = 3.8;

/// Maximum kinematic shape factor for turbulent boundary layers.
/// When Hk exceeds this value, the flow is considered separated and
/// the solver switches to inverse mode.
/// From XFOIL: HTMAX = 2.5
pub const HK_MAX_TURBULENT: f64 = 2.5;

/// Which surface of the airfoil.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Surface {
    Upper,
    Lower,
}

/// Boundary layer state at a single station.
#[derive(Debug, Clone)]
pub struct BLState {
    /// Arc-length from stagnation point
    pub s: f64,
    /// X-coordinate
    pub x: f64,
    /// Edge velocity (from inviscid solution)
    pub ue: f64,
    /// Momentum thickness θ
    pub theta: f64,
    /// Displacement thickness δ*
    pub delta_star: f64,
    /// Shape factor H = δ*/θ
    pub h: f64,
    /// Skin friction coefficient Cf
    pub cf: f64,
    /// Maximum shear stress coefficient (turbulent)
    pub ctau: f64,
    /// Whether this station is turbulent
    pub is_turbulent: bool,
    /// Amplification factor (for transition prediction)
    pub n_amp: f64,
}

impl Default for BLState {
    fn default() -> Self {
        Self {
            s: 0.0,
            x: 0.0,
            ue: 1.0,
            theta: 1e-6,
            delta_star: 2.6e-6,
            h: 2.6,
            cf: 0.0,
            ctau: 0.0,
            is_turbulent: false,
            n_amp: 0.0,
        }
    }
}

impl BLState {
    /// Create initial state at the stagnation point.
    ///
    /// Uses Hiemenz stagnation point flow solution as starting condition.
    pub fn stagnation_point(reynolds: f64) -> Self {
        // Hiemenz flow: theta/x = 0.292 / sqrt(Re_x * dU/dx)
        // At stagnation, we use a small initial theta
        let theta_init = 1e-4 / reynolds.sqrt();
        let h_init = 2.59; // Blasius value
        
        Self {
            s: 0.0,
            x: 0.0,
            ue: 0.0,
            theta: theta_init,
            delta_star: h_init * theta_init,
            h: h_init,
            cf: 0.0,
            ctau: 0.0,
            is_turbulent: false,
            n_amp: 0.0,
        }
    }

    /// Momentum thickness Reynolds number Re_θ = Ue * θ / ν
    pub fn re_theta(&self, reynolds: f64) -> f64 {
        self.ue * self.theta * reynolds
    }

    /// Displacement thickness Reynolds number Re_δ* = Ue * δ* / ν
    pub fn re_delta_star(&self, reynolds: f64) -> f64 {
        self.ue * self.delta_star * reynolds
    }

    /// Kinematic shape factor Hk (used in XFOIL)
    /// For incompressible flow, Hk = H
    pub fn hk(&self) -> f64 {
        self.h
    }

    /// Energy thickness shape factor Hs = θ*/θ
    /// Approximation for turbulent flow
    pub fn hs(&self) -> f64 {
        if self.is_turbulent {
            // Turbulent correlation (simplified)
            let hk = self.hk();
            if hk < 4.0 {
                1.505 + 4.0 / (hk - 1.0).max(0.01) + 0.165 - 1.6 / hk.max(1.01)
            } else {
                1.505 + 4.0 / (hk - 1.0) + (hk - 4.0).powi(2) / (6.0 * hk)
            }
        } else {
            // Laminar: from Falkner-Skan
            let hk = self.hk();
            0.0111 * (hk - 1.0).powi(2) / (hk - 1.0 + 0.0278) + 1.528 - 0.0002 * (hk.powi(2) - 1.0)
        }
    }
    
    /// Check if this station is in separated/inverse mode.
    /// 
    /// Separation is detected when the kinematic shape factor Hk exceeds
    /// the threshold for the current flow regime:
    /// - Laminar: Hk >= 3.8 (HK_MAX_LAMINAR)
    /// - Turbulent: Hk >= 2.5 (HK_MAX_TURBULENT)
    /// 
    /// In separated regions, the BL solver should switch to inverse mode,
    /// where Hk is prescribed and Ue is solved for.
    pub fn is_separated(&self) -> bool {
        let hk = self.hk();
        let hk_max = if self.is_turbulent {
            HK_MAX_TURBULENT
        } else {
            HK_MAX_LAMINAR
        };
        hk >= hk_max
    }
    
    /// Get the maximum shape factor threshold for current flow state.
    pub fn hk_max(&self) -> f64 {
        if self.is_turbulent {
            HK_MAX_TURBULENT
        } else {
            HK_MAX_LAMINAR
        }
    }
}

/// A boundary layer station with surface and index info.
#[derive(Debug, Clone)]
pub struct BLStation {
    /// Which surface
    pub surface: Surface,
    /// Index in the original point array
    pub idx: usize,
    /// BL state at this station
    pub state: BLState,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stagnation_point() {
        let state = BLState::stagnation_point(1e6);
        assert!(state.theta > 0.0);
        assert!(state.theta < 1e-3);
        assert!((state.h - 2.59).abs() < 0.1);
        assert!(!state.is_turbulent);
    }

    #[test]
    fn test_re_theta() {
        let mut state = BLState::default();
        state.ue = 1.0;
        state.theta = 0.001;
        let re_theta = state.re_theta(1e6);
        assert!((re_theta - 1000.0).abs() < 0.1);
    }
    
    #[test]
    fn test_is_separated_laminar() {
        let mut state = BLState::default();
        state.is_turbulent = false;
        
        // Below threshold
        state.h = 3.0;
        assert!(!state.is_separated());
        
        // At threshold
        state.h = HK_MAX_LAMINAR;
        assert!(state.is_separated());
        
        // Above threshold
        state.h = 4.5;
        assert!(state.is_separated());
    }
    
    #[test]
    fn test_is_separated_turbulent() {
        let mut state = BLState::default();
        state.is_turbulent = true;
        
        // Below threshold
        state.h = 2.0;
        assert!(!state.is_separated());
        
        // At threshold
        state.h = HK_MAX_TURBULENT;
        assert!(state.is_separated());
        
        // Above threshold
        state.h = 3.0;
        assert!(state.is_separated());
    }
}
