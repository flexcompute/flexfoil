//! Solution computation and forces (XFOIL's SPECAL, QISET, CLCALC).
//!
//! This module handles:
//! - Combining base solutions for arbitrary angle of attack
//! - Computing pressure coefficients
//! - Integrating forces (lift and moment)
//!
//! # XFOIL Reference
//!
//! - `xoper.f`: SPECAL, SPECCL, QISET, CLCALC, CPCALC

use std::f64::consts::PI;

/// Flow conditions for the analysis.
#[derive(Debug, Clone, Copy)]
pub struct FlowConditions {
    /// Angle of attack in radians
    pub alpha: f64,
    /// Freestream velocity magnitude (typically normalized to 1.0)
    pub v_inf: f64,
    /// Mach number (0.0 for incompressible)
    pub mach: f64,
}

impl Default for FlowConditions {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            v_inf: 1.0,
            mach: 0.0,
        }
    }
}

impl FlowConditions {
    /// Create flow conditions with the given angle of attack (in degrees).
    pub fn with_alpha_deg(alpha_deg: f64) -> Self {
        Self {
            alpha: alpha_deg.to_radians(),
            v_inf: 1.0,
            mach: 0.0,
        }
    }

    /// Create flow conditions with the given angle of attack (in radians).
    pub fn with_alpha_rad(alpha_rad: f64) -> Self {
        Self {
            alpha: alpha_rad,
            v_inf: 1.0,
            mach: 0.0,
        }
    }

    /// Set the freestream velocity.
    pub fn with_velocity(mut self, v_inf: f64) -> Self {
        self.v_inf = v_inf;
        self
    }

    /// Set the Mach number.
    pub fn with_mach(mut self, mach: f64) -> Self {
        self.mach = mach;
        self
    }

    /// Get angle of attack in degrees.
    pub fn alpha_deg(&self) -> f64 {
        self.alpha.to_degrees()
    }

    /// Cosine of angle of attack.
    pub fn cosa(&self) -> f64 {
        self.alpha.cos()
    }

    /// Sine of angle of attack.
    pub fn sina(&self) -> f64 {
        self.alpha.sin()
    }
}

/// Solution from the inviscid panel method.
#[derive(Debug, Clone)]
pub struct InviscidSolution {
    /// Vorticity values at each node (γᵢ = surface tangential velocity)
    pub gamma: Vec<f64>,
    /// Vorticity derivative w.r.t. alpha (∂γ/∂α)
    pub gamma_a: Vec<f64>,
    /// Pressure coefficient at each node
    pub cp: Vec<f64>,
    /// Lift coefficient (per unit span)
    pub cl: f64,
    /// Moment coefficient about quarter-chord
    pub cm: f64,
    /// Internal stream function value
    pub psi_0: f64,
    /// Number of nodes
    pub n: usize,
}

impl InviscidSolution {
    /// Get the edge velocity at a node (same as gamma for inviscid).
    pub fn edge_velocity(&self, i: usize) -> f64 {
        self.gamma[i]
    }

    /// Get the tangential velocity (same as gamma).
    pub fn qinv(&self, i: usize) -> f64 {
        self.gamma[i]
    }

    /// Find the stagnation point (where gamma crosses zero).
    ///
    /// Returns (ist, sst_fraction) where:
    /// - ist: panel index containing stagnation
    /// - sst_fraction: fractional position within panel [0, 1]
    pub fn find_stagnation(&self) -> Option<(usize, f64)> {
        let n = self.n;

        // Look for sign change in gamma
        for i in 0..n - 1 {
            if self.gamma[i] >= 0.0 && self.gamma[i + 1] < 0.0 {
                // Stagnation between i and i+1
                let dgam = self.gamma[i + 1] - self.gamma[i];
                let frac = if dgam.abs() > 1e-12 {
                    -self.gamma[i] / dgam
                } else {
                    0.5
                };
                return Some((i, frac.clamp(0.0, 1.0)));
            }
        }

        None
    }

    /// Get maximum velocity magnitude.
    pub fn max_velocity(&self) -> f64 {
        self.gamma
            .iter()
            .map(|g| g.abs())
            .fold(0.0, f64::max)
    }

    /// Get minimum Cp (suction peak).
    pub fn min_cp(&self) -> f64 {
        self.cp.iter().cloned().fold(f64::INFINITY, f64::min)
    }

    /// Check if flow is valid (no NaN/Inf).
    pub fn is_valid(&self) -> bool {
        self.gamma.iter().all(|g| g.is_finite())
            && self.cp.iter().all(|c| c.is_finite())
            && self.cl.is_finite()
            && self.cm.is_finite()
    }
}

/// Compute pressure coefficient from velocity.
///
/// # Arguments
///
/// * `gamma` - Surface tangential velocity (= vortex strength)
/// * `v_inf` - Freestream velocity
/// * `mach` - Mach number (0 for incompressible)
///
/// # Returns
///
/// Pressure coefficient Cp
pub fn compute_cp(gamma: f64, v_inf: f64, mach: f64) -> f64 {
    let q_ratio = gamma / v_inf;
    
    if mach < 0.01 {
        // Incompressible: Cp = 1 - (V/V∞)²
        1.0 - q_ratio * q_ratio
    } else {
        // Compressible (Karman-Tsien or Prandtl-Glauert)
        // For now, use incompressible + Prandtl-Glauert correction
        let beta = (1.0 - mach * mach).sqrt();
        (1.0 - q_ratio * q_ratio) / beta
    }
}

/// Compute lift curve slope estimate.
///
/// For thin airfoils, the theoretical value is 2π/rad.
pub fn thin_airfoil_cl_alpha() -> f64 {
    2.0 * PI
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flow_conditions_default() {
        let flow = FlowConditions::default();
        assert!((flow.alpha - 0.0).abs() < 1e-10);
        assert!((flow.v_inf - 1.0).abs() < 1e-10);
        assert!((flow.mach - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_flow_conditions_degrees() {
        let flow = FlowConditions::with_alpha_deg(4.0);
        assert!((flow.alpha_deg() - 4.0).abs() < 1e-10);
        assert!((flow.alpha - 4.0_f64.to_radians()).abs() < 1e-10);
    }

    #[test]
    fn test_flow_conditions_builder() {
        let flow = FlowConditions::with_alpha_deg(5.0)
            .with_velocity(2.0)
            .with_mach(0.3);
        
        assert!((flow.alpha_deg() - 5.0).abs() < 1e-10);
        assert!((flow.v_inf - 2.0).abs() < 1e-10);
        assert!((flow.mach - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_cp_incompressible() {
        // At stagnation: V=0, Cp=1
        assert!((compute_cp(0.0, 1.0, 0.0) - 1.0).abs() < 1e-10);
        
        // At freestream velocity: V=V∞, Cp=0
        assert!((compute_cp(1.0, 1.0, 0.0) - 0.0).abs() < 1e-10);
        
        // Accelerated flow: V=1.5*V∞, Cp=-1.25
        assert!((compute_cp(1.5, 1.0, 0.0) - (-1.25)).abs() < 1e-10);
    }

    #[test]
    fn test_thin_airfoil_theory() {
        let cl_alpha = thin_airfoil_cl_alpha();
        assert!((cl_alpha - 2.0 * PI).abs() < 1e-10);
    }
}
