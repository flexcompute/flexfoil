//! Viscous solver configuration.
//!
//! This module defines the configuration parameters for the viscous solver,
//! including flow conditions (Reynolds, Mach), transition criteria (Ncrit),
//! and iteration control parameters.
//!
//! # XFOIL Equivalents
//! - Reynolds → REINF1
//! - Mach → MINF1
//! - Ncrit → ACRIT
//! - max_iterations → VARONE iteration limit
//! - tolerance → RMSBL convergence threshold

/// Configuration for the viscous solver.
///
/// Controls the physical flow conditions and numerical parameters
/// for the viscous-inviscid coupled solution.
#[derive(Debug, Clone)]
pub struct ViscousSolverConfig {
    /// Reynolds number (based on chord and freestream velocity).
    /// Typical values: 1e5 - 1e7 for airfoils
    pub reynolds: f64,

    /// Freestream Mach number.
    /// For incompressible flow, use 0.0.
    /// Compressibility corrections applied for M > 0.
    pub mach: f64,

    /// Critical N-factor for transition prediction (e^N method).
    /// - 9.0: Clean flight conditions, low turbulence
    /// - 3.0-5.0: Noisy wind tunnels, high freestream turbulence
    /// - 0.0: Forces immediate transition
    pub ncrit: f64,

    /// Maximum number of global Newton iterations.
    /// XFOIL typically converges in 10-30 iterations.
    pub max_iterations: usize,

    /// Convergence tolerance on RMS residual.
    /// Solution is converged when RMS(ΔBL/BL) < tolerance.
    pub tolerance: f64,

    /// Under-relaxation factor (0-1).
    /// - 1.0: Full Newton update (fastest but may diverge)
    /// - 0.5-0.8: Safer for difficult cases
    pub relaxation: f64,

    /// Maximum shape factor Hk for laminar flow before inverse mode.
    /// XFOIL default: 3.8
    pub hk_max_laminar: f64,

    /// Maximum shape factor Hk for turbulent flow before inverse mode.
    /// XFOIL default: 2.5
    pub hk_max_turbulent: f64,

    /// Whether to allow separation bubbles.
    /// If false, solver returns error on separation.
    pub allow_separation: bool,
}

impl Default for ViscousSolverConfig {
    fn default() -> Self {
        Self {
            reynolds: 1e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 50,
            tolerance: 1e-4,
            relaxation: 1.0,
            hk_max_laminar: 3.8,
            hk_max_turbulent: 2.5,
            allow_separation: true,
        }
    }
}

impl ViscousSolverConfig {
    /// Create configuration with specified Reynolds number.
    ///
    /// Uses default values for all other parameters.
    ///
    /// # Example
    /// ```
    /// use rustfoil_solver::viscous::ViscousSolverConfig;
    /// let config = ViscousSolverConfig::with_reynolds(5e6);
    /// assert_eq!(config.reynolds, 5e6);
    /// ```
    pub fn with_reynolds(re: f64) -> Self {
        Self {
            reynolds: re,
            ..Default::default()
        }
    }

    /// Create configuration with Reynolds and Mach numbers.
    ///
    /// # Example
    /// ```
    /// use rustfoil_solver::viscous::ViscousSolverConfig;
    /// let config = ViscousSolverConfig::with_re_mach(1e6, 0.3);
    /// assert_eq!(config.mach, 0.3);
    /// ```
    pub fn with_re_mach(re: f64, mach: f64) -> Self {
        Self {
            reynolds: re,
            mach,
            ..Default::default()
        }
    }

    /// Create configuration for high-turbulence environment.
    ///
    /// Sets Ncrit = 3.0 for early transition typical of
    /// noisy wind tunnels or high freestream turbulence.
    pub fn high_turbulence(re: f64) -> Self {
        Self {
            reynolds: re,
            ncrit: 3.0,
            ..Default::default()
        }
    }

    /// Create configuration for low-turbulence (flight-like) conditions.
    ///
    /// Sets Ncrit = 9.0 for late transition typical of
    /// clean flight conditions or quiet wind tunnels.
    pub fn low_turbulence(re: f64) -> Self {
        Self {
            reynolds: re,
            ncrit: 9.0,
            ..Default::default()
        }
    }

    /// Set under-relaxation for difficult convergence cases.
    ///
    /// Returns self for builder-style chaining.
    pub fn with_relaxation(mut self, relax: f64) -> Self {
        self.relaxation = relax.clamp(0.1, 1.0);
        self
    }

    /// Set maximum iterations.
    ///
    /// Returns self for builder-style chaining.
    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iterations = max_iter;
        self
    }

    /// Set convergence tolerance.
    ///
    /// Returns self for builder-style chaining.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Compute Mach number squared (cached for closure evaluations).
    #[inline]
    pub fn msq(&self) -> f64 {
        self.mach * self.mach
    }

    /// Validate configuration parameters.
    ///
    /// Returns an error if any parameters are invalid.
    pub fn validate(&self) -> Result<(), &'static str> {
        if self.reynolds <= 0.0 {
            return Err("Reynolds number must be positive");
        }
        if self.mach < 0.0 || self.mach >= 1.0 {
            return Err("Mach number must be in range [0, 1)");
        }
        if self.ncrit < 0.0 {
            return Err("Ncrit must be non-negative");
        }
        if self.tolerance <= 0.0 {
            return Err("Tolerance must be positive");
        }
        if self.relaxation <= 0.0 || self.relaxation > 1.0 {
            return Err("Relaxation must be in range (0, 1]");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = ViscousSolverConfig::default();
        assert_eq!(config.reynolds, 1e6);
        assert_eq!(config.mach, 0.0);
        assert_eq!(config.ncrit, 9.0);
        assert_eq!(config.max_iterations, 50);
        assert!((config.tolerance - 1e-4).abs() < 1e-10);
        assert_eq!(config.relaxation, 1.0);
    }

    #[test]
    fn test_with_reynolds() {
        let config = ViscousSolverConfig::with_reynolds(5e6);
        assert_eq!(config.reynolds, 5e6);
        assert_eq!(config.ncrit, 9.0); // Should keep default
    }

    #[test]
    fn test_with_re_mach() {
        let config = ViscousSolverConfig::with_re_mach(2e6, 0.5);
        assert_eq!(config.reynolds, 2e6);
        assert_eq!(config.mach, 0.5);
    }

    #[test]
    fn test_turbulence_presets() {
        let high = ViscousSolverConfig::high_turbulence(1e6);
        assert_eq!(high.ncrit, 3.0);

        let low = ViscousSolverConfig::low_turbulence(1e6);
        assert_eq!(low.ncrit, 9.0);
    }

    #[test]
    fn test_builder_pattern() {
        let config = ViscousSolverConfig::with_reynolds(1e6)
            .with_relaxation(0.8)
            .with_max_iterations(100)
            .with_tolerance(1e-5);

        assert_eq!(config.reynolds, 1e6);
        assert_eq!(config.relaxation, 0.8);
        assert_eq!(config.max_iterations, 100);
        assert!((config.tolerance - 1e-5).abs() < 1e-10);
    }

    #[test]
    fn test_relaxation_clamping() {
        let config = ViscousSolverConfig::default().with_relaxation(2.0);
        assert_eq!(config.relaxation, 1.0); // Clamped to max

        let config = ViscousSolverConfig::default().with_relaxation(0.01);
        assert_eq!(config.relaxation, 0.1); // Clamped to min
    }

    #[test]
    fn test_msq() {
        let config = ViscousSolverConfig::with_re_mach(1e6, 0.5);
        assert!((config.msq() - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_validation() {
        // Valid config
        assert!(ViscousSolverConfig::default().validate().is_ok());

        // Invalid Reynolds
        let mut config = ViscousSolverConfig::default();
        config.reynolds = -1.0;
        assert!(config.validate().is_err());

        // Invalid Mach
        config = ViscousSolverConfig::default();
        config.mach = 1.0;
        assert!(config.validate().is_err());

        // Invalid Ncrit
        config = ViscousSolverConfig::default();
        config.ncrit = -1.0;
        assert!(config.validate().is_err());

        // Invalid tolerance
        config = ViscousSolverConfig::default();
        config.tolerance = 0.0;
        assert!(config.validate().is_err());

        // Invalid relaxation
        config = ViscousSolverConfig::default();
        config.relaxation = 0.0;
        assert!(config.validate().is_err());
    }
}
