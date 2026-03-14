//! Boundary layer state at a single station
//!
//! Matches XFOIL's XBL.INC common block structure.
//!
//! # XFOIL Reference
//! The XBL.INC common block contains 73+ variables per BL station:
//! - Primary variables: X, U, T (theta), D (delta_star), S (ctau), AMPL
//! - Secondary variables: H, HK, HS, HC, RT, CF, DI, US, CQ, DE
//! - Mode flags: SIMI, TRAN, TURB, WAKE
//! - 50+ partial derivatives for Jacobian construction

/// Partial derivatives of secondary variables w.r.t. primary variables
///
/// These are needed to construct the Newton system Jacobian.
/// Notation: `foo_bar` means ∂foo/∂bar
#[derive(Debug, Clone, Default)]
pub struct BlDerivatives {
    // === H derivatives ===
    /// ∂H/∂θ
    pub h_theta: f64,
    /// ∂H/∂δ*
    pub h_delta_star: f64,

    // === Hk derivatives ===
    /// ∂Hk/∂H
    pub hk_h: f64,
    /// ∂Hk/∂M²
    pub hk_msq: f64,

    // === Hs derivatives ===
    /// ∂Hs/∂Hk
    pub hs_hk: f64,
    /// ∂Hs/∂Rθ
    pub hs_rt: f64,
    /// ∂Hs/∂M²
    pub hs_msq: f64,

    // === Cf derivatives ===
    /// ∂Cf/∂Hk
    pub cf_hk: f64,
    /// ∂Cf/∂Rθ
    pub cf_rt: f64,
    /// ∂Cf/∂M²
    pub cf_msq: f64,

    // === CD (dissipation) derivatives ===
    /// ∂CD/∂Hk
    pub cd_hk: f64,
    /// ∂CD/∂Rθ
    pub cd_rt: f64,
    /// ∂CD/∂θ
    pub cd_t: f64,
    /// ∂CD/∂δ*
    pub cd_d: f64,
    /// ∂CD/∂Ue
    pub cd_u: f64,
    /// ∂CD/∂S (ctau)
    pub cd_s: f64,

    // === Hc (density shape factor) derivatives ===
    /// ∂Hc/∂Hk
    pub hc_hk: f64,
    /// ∂Hc/∂M²
    pub hc_msq: f64,

    // === Us (slip velocity) derivatives ===
    /// ∂Us/∂θ
    pub us_t: f64,
    /// ∂Us/∂δ*
    pub us_d: f64,
    /// ∂Us/∂Ue
    pub us_u: f64,

    // === CQ (equilibrium shear) derivatives ===
    /// ∂CQ/∂θ
    pub cq_t: f64,
    /// ∂CQ/∂δ*
    pub cq_d: f64,
    /// ∂CQ/∂Ue
    pub cq_u: f64,

    // === DE (energy thickness) derivatives ===
    /// ∂DE/∂θ
    pub de_t: f64,
    /// ∂DE/∂δ*
    pub de_d: f64,

    // === Edge velocity derivative ===
    /// dUe/dx (velocity gradient)
    pub u_x: f64,
}

/// Complete boundary layer state at one station
///
/// This struct holds all the state variables needed to describe the boundary
/// layer at a single streamwise location. It includes:
/// - Primary variables that are the unknowns in the Newton iteration
/// - Secondary variables computed from closures (HKIN, HS, CF, etc.)
/// - Mode flags indicating laminar/turbulent/wake state
/// - Partial derivatives for Jacobian construction
///
/// # XFOIL Correspondence
/// Maps to XFOIL's V_VAR1/V_VAR2 common blocks in XBL.INC:
/// - `x` → X1/X2 (arc length)
/// - `u` → U1/U2 (edge velocity Ue)
/// - `theta` → T1/T2 (momentum thickness θ)
/// - `delta_star` → D1/D2 (displacement thickness δ*)
/// - `ctau` → S1/S2 (shear stress coefficient √Cτ)
/// - `ampl` → AMPL1/AMPL2 (amplification factor N)
#[derive(Debug, Clone)]
pub struct BlStation {
    // === Primary Variables (Newton unknowns) ===
    /// Arc length position along surface (used in BL equations)
    pub x: f64,
    /// Panel x-coordinate for reporting x/c transition location (not arc length)
    pub x_coord: f64,
    /// Panel index in the global DIJ matrix (for VI coupling)
    /// Maps this BL station to the corresponding panel in the inviscid solution.
    pub panel_idx: usize,
    /// Edge velocity Ue (normalized by freestream)
    pub u: f64,
    /// Momentum thickness θ
    pub theta: f64,
    /// Displacement thickness δ*
    pub delta_star: f64,
    /// Shear stress coefficient √Cτ (turbulent only)
    pub ctau: f64,
    /// Amplification factor N (laminar only, for transition prediction)
    pub ampl: f64,

    // === Secondary Variables (computed by closures) ===
    /// Shape factor H = δ*/θ
    pub h: f64,
    /// Kinematic shape factor Hk (compressibility-corrected)
    pub hk: f64,
    /// Energy shape factor Hs
    pub hs: f64,
    /// Density shape factor Hc (compressible flows)
    pub hc: f64,
    /// Reynolds number Rθ = Ue·θ/ν
    pub r_theta: f64,
    /// Skin friction coefficient Cf
    pub cf: f64,
    /// Dissipation coefficient CD
    pub cd: f64,
    /// Normalized slip velocity Us (XFOIL US)
    pub us: f64,
    /// Equilibrium shear stress coefficient CQ (√Cτ_eq)
    pub cq: f64,
    /// Energy thickness DE from Green's correlation
    pub de: f64,
    /// Mass defect Ue·δ* (used in viscous-inviscid coupling)
    pub mass_defect: f64,

    /// Wake gap thickness DW (XFOIL's DSWAKI from WGAP array).
    /// For non-wake stations this is 0. BLPRV sets D2 = DSI - DSWAKI, DW2 = DSWAKI.
    pub dw: f64,

    // === Mode Flags ===
    /// True if boundary layer is laminar
    pub is_laminar: bool,
    /// True if in wake region (behind trailing edge)
    pub is_wake: bool,
    /// True if boundary layer is turbulent
    pub is_turbulent: bool,

    // === Partial Derivatives ===
    /// All partial derivatives for Jacobian construction
    pub derivs: BlDerivatives,
}

impl BlStation {
    /// Create a new station with default values suitable for initialization
    ///
    /// Default values represent a typical attached laminar boundary layer:
    /// - H = 2.0 (typical laminar shape factor)
    /// - Hk = 2.0 (incompressible limit)
    /// - Hs = 1.5 (typical energy shape factor)
    /// - Rθ = 1000 (moderate Reynolds number)
    /// - Cf = 0.003 (typical laminar skin friction)
    pub fn new() -> Self {
        Self {
            x: 0.0,
            x_coord: 0.0,
            panel_idx: 0,
            u: 1.0,
            theta: 0.001,
            delta_star: 0.002,
            ctau: 0.0,
            ampl: 0.0,
            h: 2.0,
            hk: 2.0,
            hs: 1.5,
            hc: 0.0,
            r_theta: 1000.0,
            cf: 0.003,
            cd: 0.001,
            us: 0.5,
            cq: 0.03,
            de: 0.006,
            mass_defect: 0.002,
            dw: 0.0,
            is_laminar: true,
            is_wake: false,
            is_turbulent: false,
            derivs: BlDerivatives::default(),
        }
    }

    /// Total displacement thickness used for wake mass defect.
    pub fn total_delta_star(&self) -> f64 {
        if self.is_wake {
            self.delta_star + self.dw.max(0.0)
        } else {
            self.delta_star
        }
    }

    /// Recompute the stored mass defect from the current station state.
    pub fn refresh_mass_defect(&mut self) {
        self.mass_defect = self.u * self.total_delta_star();
    }

    /// Initialize for stagnation point using Thwaites' method (XFOIL-compatible)
    ///
    /// Uses Thwaites' formula for stagnation point initialization, matching XFOIL's
    /// approach in xbl.f lines 590-598. This is critical for correct early boundary
    /// layer thickness.
    ///
    /// # Arguments
    /// * `x` - Arc length from stagnation (must be > 0, typically ~0.001)
    /// * `ue` - Edge velocity at the station (typically small but non-zero)
    /// * `re` - Reynolds number (Ue·L/ν where L is reference length)
    ///
    /// # Thwaites Formula (XFOIL xbl.f:597)
    /// With BULE = 1.0 (stagnation flow):
    /// - TSQ = 0.45 / (UCON * (5.0*BULE+1.0) * REYBL) * XSI^(1.0-BULE)
    /// - Where UCON = Ue/x and REYBL = REINF for incompressible flow
    /// - θ = √(TSQ)
    /// - δ* = 2.2·θ
    /// - H = δ*/θ = 2.2
    ///
    /// # Note on REYBL
    /// XFOIL computes REYBL = REINF * SQRT(HERAT³) * (1+HVRAT)/(HERAT+HVRAT),
    /// which reduces to REINF for incompressible flow, so the stagnation seed uses
    /// the same Reynolds scaling as the marched boundary-layer stations.
    ///
    /// # Example
    /// ```
    /// use rustfoil_bl::state::BlStation;
    /// let station = BlStation::stagnation(0.001104, 0.060676, 3e6);
    /// // theta should be ~3.69e-5 (matches XFOIL)
    /// ```
    pub fn stagnation(x: f64, ue: f64, re: f64) -> Self {
        let mut station = Self::new();
        station.x = x.max(1e-12); // Ensure x > 0
        station.u = ue.abs().max(1e-6); // Ensure Ue > 0
        // x_coord will be set by the caller (typically near LE, so ~0)
        station.x_coord = 0.0;

        // Thwaites' formula for stagnation point (XFOIL xbl.f:597)
        // BULE = 1.0 (stagnation flow: Ue = K*x)
        // TSQ = 0.45 / (UCON * (5.0*BULE+1.0) * REYBL) * XSI^(1.0-BULE)
        // With BULE=1.0: TSQ = 0.45 / ((Ue/x) * 6 * REYBL) * x^0 = 0.45*x / (6*Ue*REYBL)
        // 
        let bule = 1.0;
        let ucon = station.u / station.x.powf(bule);
        let reybl = re;
        let tsq = 0.45 / (ucon * (5.0 * bule + 1.0) * reybl) * station.x.powf(1.0 - bule);
        
        station.theta = tsq.sqrt().max(1e-12);

        // Displacement thickness: DSI = 2.2*THI (XFOIL xbl.f:578)
        station.delta_star = 2.2 * station.theta;

        // Shape factor from displacement/momentum ratio
        station.h = station.delta_star / station.theta;

        // Kinematic shape factor equals H at incompressible (M=0) stagnation
        station.hk = station.h;

        // Mass defect
        station.refresh_mass_defect();

        // Reynolds number based on momentum thickness
        station.r_theta = station.u * station.theta * re;

        // Stagnation point is always laminar
        station.is_laminar = true;
        station.is_turbulent = false;
        station.is_wake = false;

        station
    }
}

impl Default for BlStation {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blstation_new_defaults() {
        let station = BlStation::new();

        // Primary variables have sensible defaults
        assert_eq!(station.x, 0.0);
        assert_eq!(station.u, 1.0);
        assert!(station.theta > 0.0, "theta should be positive");
        assert!(station.delta_star > 0.0, "delta_star should be positive");

        // Secondary variables are physically reasonable
        assert!(
            station.h > 1.0 && station.h < 5.0,
            "H should be between 1 and 5 for physical BL"
        );
        assert!(
            station.hk > 1.0 && station.hk < 5.0,
            "Hk should be between 1 and 5"
        );
        assert!(station.hs > 1.0, "Hs should be > 1");
        assert!(station.r_theta > 0.0, "Rtheta should be positive");
        assert!(station.cf > 0.0, "Cf should be positive");
        assert!(station.cd > 0.0, "CD should be positive");

        // Default is laminar flow
        assert!(station.is_laminar, "Default should be laminar");
        assert!(!station.is_turbulent, "Default should not be turbulent");
        assert!(!station.is_wake, "Default should not be wake");
    }

    #[test]
    fn test_blstation_stagnation() {
        let x = 0.001;
        let ue = 0.1;
        let re = 1e6;
        let station = BlStation::stagnation(x, ue, re);

        // Thwaites shape factor for stagnation (BULE=1.0): H = 2.2
        assert!(
            (station.h - 2.2).abs() < 1e-10,
            "Thwaites H should be 2.2, got {}",
            station.h
        );

        // Verify theta scaling: theta = sqrt(0.45 * x / (6 * Ue * Re))
        let expected_theta = (0.45 * x / (6.0 * ue * re)).sqrt();
        assert!(
            (station.theta - expected_theta).abs() < 1e-15,
            "theta should match Thwaites scaling"
        );

        // Verify delta_star = 2.2 * theta
        let expected_delta_star = 2.2 * station.theta;
        assert!(
            (station.delta_star - expected_delta_star).abs() < 1e-15,
            "delta_star = 2.2 * theta"
        );

        // Stagnation is always laminar
        assert!(station.is_laminar);
        assert!(!station.is_turbulent);
    }

    #[test]
    fn test_blstation_stagnation_physics() {
        // Test at different Reynolds numbers to verify scaling
        let x = 0.001;
        let ue = 0.5;

        let station_low_re = BlStation::stagnation(x, ue, 1e5);
        let station_high_re = BlStation::stagnation(x, ue, 1e7);

        // Higher Re should give thinner BL (Thwaites: theta ~ 1/sqrt(Re))
        assert!(
            station_high_re.theta < station_low_re.theta,
            "Higher Re should give smaller theta"
        );

        // With Thwaites: theta = sqrt(0.45*x/(6*Ue*Re)), so ratio = sqrt(Re_low/Re_high) = sqrt(1e5/1e7) = 0.1
        let theta_ratio = station_high_re.theta / station_low_re.theta;
        let expected_ratio = (1e5_f64 / 1e7_f64).sqrt();
        assert!(
            (theta_ratio - expected_ratio).abs() < 1e-10,
            "theta should scale as 1/sqrt(Re)"
        );

        // Shape factor should be same regardless of Re
        assert!(
            (station_low_re.h - station_high_re.h).abs() < 1e-10,
            "H is independent of Re"
        );
    }

    #[test]
    fn test_blderivatives_default() {
        let derivs = BlDerivatives::default();

        // All derivatives should initialize to zero
        assert_eq!(derivs.h_theta, 0.0);
        assert_eq!(derivs.h_delta_star, 0.0);
        assert_eq!(derivs.hk_h, 0.0);
        assert_eq!(derivs.hk_msq, 0.0);
        assert_eq!(derivs.hs_hk, 0.0);
        assert_eq!(derivs.hs_rt, 0.0);
        assert_eq!(derivs.hs_msq, 0.0);
        assert_eq!(derivs.cf_hk, 0.0);
        assert_eq!(derivs.cf_rt, 0.0);
        assert_eq!(derivs.cf_msq, 0.0);
        assert_eq!(derivs.cd_hk, 0.0);
        assert_eq!(derivs.cd_rt, 0.0);
        assert_eq!(derivs.cd_t, 0.0);
        assert_eq!(derivs.cd_d, 0.0);
        assert_eq!(derivs.cd_u, 0.0);
        assert_eq!(derivs.cd_s, 0.0);
        assert_eq!(derivs.us_t, 0.0);
        assert_eq!(derivs.us_d, 0.0);
        assert_eq!(derivs.us_u, 0.0);
        assert_eq!(derivs.cq_t, 0.0);
        assert_eq!(derivs.cq_d, 0.0);
        assert_eq!(derivs.cq_u, 0.0);
        assert_eq!(derivs.de_t, 0.0);
        assert_eq!(derivs.de_d, 0.0);
        assert_eq!(derivs.u_x, 0.0);
    }

    #[test]
    fn test_blstation_mode_flags_exclusive() {
        // Test that mode flags are self-consistent in default state
        let station = BlStation::new();

        // Laminar and turbulent should be mutually exclusive
        assert!(
            !(station.is_laminar && station.is_turbulent),
            "Cannot be both laminar and turbulent"
        );

        // Default is laminar, not wake
        assert!(station.is_laminar);
        assert!(!station.is_wake);
    }

    #[test]
    fn test_blstation_clone() {
        let mut station = BlStation::new();
        station.x = 0.5;
        station.theta = 0.002;
        station.derivs.hk_h = 1.0;

        let cloned = station.clone();

        assert_eq!(cloned.x, 0.5);
        assert_eq!(cloned.theta, 0.002);
        assert_eq!(cloned.derivs.hk_h, 1.0);

        // Verify it's a true clone (modifying original doesn't affect clone)
        // This is implicitly tested by the Clone derive working correctly
    }

    #[test]
    fn test_blstation_default_trait() {
        let station: BlStation = Default::default();

        // Should be identical to new()
        let new_station = BlStation::new();

        assert_eq!(station.x, new_station.x);
        assert_eq!(station.u, new_station.u);
        assert_eq!(station.h, new_station.h);
        assert_eq!(station.is_laminar, new_station.is_laminar);
    }

    #[test]
    fn test_blstation_mass_defect_consistency() {
        let station = BlStation::new();

        // Mass defect should equal u * delta_star
        let expected_mass_defect = station.u * station.delta_star;
        assert!(
            (station.mass_defect - expected_mass_defect).abs() < 1e-15,
            "mass_defect should equal u * delta_star"
        );
    }
}
