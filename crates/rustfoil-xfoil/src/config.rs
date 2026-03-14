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
