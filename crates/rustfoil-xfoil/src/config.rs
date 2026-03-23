#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReType {
    /// Type 1: constant Re (default). `Re_eff = Re`.
    Type1,
    /// Type 2: fixed `Re * sqrt(CL)`. `Re_eff = Re / sqrt(|CL|)`.
    Type2,
    /// Type 3: fixed `Re * CL`. `Re_eff = Re / |CL|`.
    Type3,
}

impl Default for ReType {
    fn default() -> Self {
        Self::Type1
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OperatingMode {
    PrescribedAlpha,
    PrescribedCl { target_cl: f64 },
}

impl Default for OperatingMode {
    fn default() -> Self {
        Self::PrescribedAlpha
    }
}

#[derive(Debug, Clone)]
pub struct XfoilOptions {
    pub reynolds: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub max_iterations: usize,
    pub tolerance: f64,
    pub wake_length_chords: f64,
    pub operating_mode: OperatingMode,
    pub re_type: ReType,
}

impl Default for XfoilOptions {
    fn default() -> Self {
        Self {
            reynolds: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 100,
            tolerance: 1.0e-4,
            wake_length_chords: 1.0,
            operating_mode: OperatingMode::PrescribedAlpha,
            re_type: ReType::Type1,
        }
    }
}

impl XfoilOptions {
    pub fn with_target_cl(mut self, target_cl: f64) -> Self {
        self.operating_mode = OperatingMode::PrescribedCl { target_cl };
        self
    }

    pub fn with_prescribed_alpha(mut self) -> Self {
        self.operating_mode = OperatingMode::PrescribedAlpha;
        self
    }

    pub fn effective_reynolds(&self, cl: f64) -> f64 {
        if cl.abs() < 1e-6 {
            return self.reynolds;
        }
        match self.re_type {
            ReType::Type1 => self.reynolds,
            ReType::Type2 => self.reynolds / cl.abs().sqrt(),
            ReType::Type3 => self.reynolds / cl.abs(),
        }
    }

    pub fn validate(&self) -> Result<(), &'static str> {
        if self.reynolds <= 0.0 {
            return Err("Reynolds number must be positive");
        }
        if !(0.0..1.0).contains(&self.mach) {
            return Err("Mach number must be in [0, 1)");
        }
        if self.ncrit < 0.0 {
            return Err("Ncrit must be non-negative");
        }
        if self.max_iterations == 0 {
            return Err("max_iterations must be non-zero");
        }
        if self.tolerance <= 0.0 {
            return Err("tolerance must be positive");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(re_type: ReType, reynolds: f64) -> XfoilOptions {
        XfoilOptions {
            reynolds,
            re_type,
            ..Default::default()
        }
    }

    #[test]
    fn type1_returns_nominal() {
        let o = opts(ReType::Type1, 1e6);
        assert_eq!(o.effective_reynolds(0.5), 1e6);
        assert_eq!(o.effective_reynolds(-1.2), 1e6);
    }

    #[test]
    fn type2_scales_by_sqrt_cl() {
        let o = opts(ReType::Type2, 1e6);
        let re = o.effective_reynolds(0.25);
        // Re / sqrt(0.25) = 1e6 / 0.5 = 2e6
        assert!((re - 2e6).abs() < 1.0);
    }

    #[test]
    fn type3_scales_by_cl() {
        let o = opts(ReType::Type3, 1e6);
        let re = o.effective_reynolds(0.5);
        // Re / 0.5 = 2e6
        assert!((re - 2e6).abs() < 1.0);
    }

    #[test]
    fn negative_cl_uses_abs() {
        let o = opts(ReType::Type2, 1e6);
        let re_pos = o.effective_reynolds(0.5);
        let re_neg = o.effective_reynolds(-0.5);
        assert!((re_pos - re_neg).abs() < 1.0);
    }

    #[test]
    fn near_zero_cl_returns_nominal() {
        let o = opts(ReType::Type2, 1e6);
        assert_eq!(o.effective_reynolds(0.0), 1e6);
        assert_eq!(o.effective_reynolds(1e-8), 1e6);

        let o3 = opts(ReType::Type3, 1e6);
        assert_eq!(o3.effective_reynolds(0.0), 1e6);
    }

    #[test]
    fn default_is_type1() {
        assert_eq!(ReType::default(), ReType::Type1);
        let o = XfoilOptions::default();
        assert_eq!(o.re_type, ReType::Type1);
    }
}
