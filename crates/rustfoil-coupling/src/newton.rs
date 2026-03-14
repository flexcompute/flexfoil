//! Newton system construction for viscous-inviscid coupling
//!
//! This module implements the Newton system builder for the coupled
//! viscous-inviscid iteration, porting XFOIL's BLSYS subroutine.
//!
//! # System Structure
//!
//! The BL Newton system has the form:
//! ```text
//! |       ||dA1|     |       ||dA2|       |     |
//! |  VS1  ||dT1|  +  |  VS2  ||dT2|   =   |VSREZ|
//! |       ||dD1|     |       ||dD2|       |     |
//!          |dU1|              |dU2|
//!          |dX1|              |dX2|
//!
//!   3x5    5x1         3x5    5x1          3x1
//! ```
//!
//! This forms a block-tridiagonal system where each station has 3 BL equations
//! (amplification/shear-lag, momentum, shape parameter) with coupling to
//! edge velocity through the VS terms.
//!
//! # XFOIL Reference
//! - BLSYS: xblsys.f line 583
//! - SETBL: xbl.f line 21

use nalgebra::{DMatrix, DVector};
use rustfoil_bl::closures::Trchek2FullResult;
use rustfoil_bl::equations::{bldif, trdif_full, FlowType};
use rustfoil_bl::state::BlStation;

/// Block-tridiagonal system for BL Newton iteration
///
/// This structure stores the linearized BL equations in block-tridiagonal form.
/// Each interval between stations i-1 and i contributes:
/// - VA[i]: Diagonal block (derivatives w.r.t. station i variables)
/// - VB[i]: Lower diagonal block (derivatives w.r.t. station i-1 variables)
/// - RHS[i]: Residual vector
/// - VS[i]: Edge velocity coupling (for viscous-inviscid interaction)
///
/// The 3 equations at each station are:
/// 1. Amplification (laminar) or shear-lag (turbulent/wake)
/// 2. Momentum integral
/// 3. Shape parameter (kinetic energy)
///
/// # XFOIL Reference
/// XFOIL xblsys.f BLSYS (line 583)
#[derive(Debug, Clone)]
pub struct BlNewtonSystem {
    /// Number of stations
    pub n: usize,
    /// Diagonal blocks VA[i] (3x3) - derivatives w.r.t. downstream station
    /// Index 0 is unused (boundary condition), indices 1..n are valid
    pub va: Vec<[[f64; 3]; 3]>,
    /// Lower diagonal blocks VB[i] (3x3) - derivatives w.r.t. upstream station
    /// Index 0 is unused, indices 1..n are valid
    pub vb: Vec<[[f64; 3]; 3]>,
    /// Right-hand side residuals (3 per interval)
    /// Index 0 is unused, indices 1..n are valid
    pub rhs: Vec<[f64; 3]>,
    /// Coupling to inviscid (dUe influence)
    /// These are the VS2[:,3] entries from bldif - edge velocity sensitivity
    pub vs: Vec<[f64; 3]>,
    /// VM matrix for mass defect coupling (global viscous-inviscid interaction)
    /// VM[i][j][k] = sensitivity of equation k at station i to mass defect at station j
    /// This is the key matrix for full Newton coupling as in XFOIL's SETBL
    /// Dimensions: n x n x 3
    pub vm: Vec<Vec<[f64; 3]>>,
    /// VS columns for delta_star derivatives (column 2 = index 2)
    /// Used in VM construction: VS1[:,2] and VS2[:,2]
    pub vs_delta: Vec<[f64; 3]>,
    /// VS columns for Ue derivatives (column 3 = index 3)
    /// Used in VM construction: VS1[:,3] and VS2[:,3]
    pub vs_ue: Vec<[f64; 3]>,
    /// Surface direction signs (VTI in XFOIL)
    /// +1.0 for upper surface stations, -1.0 for lower surface stations
    /// Used in VM matrix construction: U_M = -VTI[i]*VTI[j]*DIJ[i,j]
    /// For single-surface mode, all values are +1.0
    pub vti: Vec<f64>,
    /// VS1 columns for delta_star derivatives (stored for forced changes)
    /// VS1[:,2] from bldif - upstream station delta_star sensitivity
    pub vs1_delta: Vec<[f64; 3]>,
    /// VS1 columns for Ue derivatives (stored for forced changes)
    /// VS1[:,3] from bldif - upstream station Ue sensitivity
    pub vs1_ue: Vec<[f64; 3]>,
}

impl BlNewtonSystem {
    /// Create a new Newton system for n stations
    ///
    /// # Arguments
    /// * `n` - Number of BL stations (including stagnation point)
    ///
    /// # Returns
    /// A new BlNewtonSystem with all blocks initialized to zero
    pub fn new(n: usize) -> Self {
        Self {
            n,
            va: vec![[[0.0; 3]; 3]; n],
            vb: vec![[[0.0; 3]; 3]; n],
            rhs: vec![[0.0; 3]; n],
            vs: vec![[0.0; 3]; n],
            vm: vec![vec![[0.0; 3]; n]; n],
            vs_delta: vec![[0.0; 3]; n],
            vs_ue: vec![[0.0; 3]; n],
            vti: vec![1.0; n], // Default to +1 for single surface
            vs1_delta: vec![[0.0; 3]; n],
            vs1_ue: vec![[0.0; 3]; n],
        }
    }

    /// Build the Newton system from BL stations
    ///
    /// This method constructs the block-tridiagonal Newton system by calling
    /// `bldif()` for each interval between adjacent stations. The resulting
    /// Jacobian blocks and residuals are extracted and stored.
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval (Laminar, Turbulent, or Wake)
    /// * `msq` - Mach number squared (M²)
    /// * `re` - Reference Reynolds number
    ///
    /// # Note
    /// The `stations` array should have secondary variables already computed
    /// via `blvar()`. This method extracts the 3x3 core blocks from the full
    /// 3x5 Jacobian matrices returned by `bldif()`.
    ///
    /// # Reference
    /// XFOIL xblsys.f BLSYS (line 583)
    pub fn build(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
    ) {
        assert!(
            stations.len() >= 2,
            "Need at least 2 stations to build Newton system"
        );
        assert_eq!(
            stations.len(),
            self.n,
            "Station count must match system size"
        );
        assert_eq!(
            flow_types.len(),
            self.n - 1,
            "Need n-1 flow types for n stations"
        );

        // For each interval (between stations i-1 and i)
        for i in 1..self.n {
            let s1 = &stations[i - 1]; // Upstream station
            let s2 = &stations[i]; // Downstream station
            let flow_type = flow_types[i - 1];

            // Compute residuals and Jacobian for this interval
            let (residuals, jacobian) = bldif(s1, s2, flow_type, msq, re);

            // Store residuals (negated to match Newton iteration convention)
            // Note: bldif already negates residuals, so we use them directly
            self.rhs[i] = [residuals.res_third, residuals.res_mom, residuals.res_shape];

            // Extract 3x3 blocks from the 3x5 Jacobian matrices
            // The first 3 columns correspond to [ampl/ctau, theta, delta_star]
            //
            // VA[i] comes from vs2 (downstream station derivatives)
            // VB[i] comes from vs1 (upstream station derivatives)

            for eq in 0..3 {
                for var in 0..3 {
                    self.va[i][eq][var] = jacobian.vs2[eq][var];
                    self.vb[i][eq][var] = jacobian.vs1[eq][var];
                }
                // Store edge velocity coupling (column 3 = u variable)
                self.vs[i][eq] = jacobian.vs2[eq][3];
                // Store VS2 delta_star derivatives (column 2)
                self.vs_delta[i][eq] = jacobian.vs2[eq][2];
                // Store VS2 Ue derivatives (column 3)
                self.vs_ue[i][eq] = jacobian.vs2[eq][3];
                // Store VS1 delta_star derivatives (column 2) for forced changes
                self.vs1_delta[i][eq] = jacobian.vs1[eq][2];
                // Store VS1 Ue derivatives (column 3) for forced changes
                self.vs1_ue[i][eq] = jacobian.vs1[eq][3];
            }
        }
    }

    /// Build the Newton system with TRDIF support for transition intervals
    ///
    /// This method extends `build()` to handle transition intervals using XFOIL's
    /// TRDIF formulation. When a transition result is provided for an interval,
    /// the system uses `trdif_full()` to properly compute the hybrid laminar-turbulent
    /// residuals and Jacobians with full chain rule through the transition point.
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval (Laminar, Turbulent, or Wake)
    /// * `msq` - Mach number squared (M²)
    /// * `re` - Reference Reynolds number
    /// * `transitions` - Optional transition results for each interval. If `Some`, 
    ///   the interval uses TRDIF; if `None`, regular BLDIF is used.
    ///
    /// # Note
    /// The `transitions` slice should have length `n-1` (same as `flow_types`).
    /// Each entry can be `None` (use regular BLDIF) or `Some(Trchek2FullResult)`
    /// (use TRDIF with full Jacobian transformation).
    ///
    /// # Reference
    /// XFOIL xblsys.f TRDIF (lines 1195-1549)
    pub fn build_with_transition(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
        transitions: &[Option<Trchek2FullResult>],
    ) {
        assert!(
            stations.len() >= 2,
            "Need at least 2 stations to build Newton system"
        );
        assert_eq!(
            stations.len(),
            self.n,
            "Station count must match system size"
        );
        assert_eq!(
            flow_types.len(),
            self.n - 1,
            "Need n-1 flow types for n stations"
        );
        assert_eq!(
            transitions.len(),
            self.n - 1,
            "Need n-1 transition entries for n stations"
        );

        // For each interval (between stations i-1 and i)
        for i in 1..self.n {
            let s1 = &stations[i - 1]; // Upstream station
            let s2 = &stations[i]; // Downstream station
            let flow_type = flow_types[i - 1];
            let transition = &transitions[i - 1];

            // Compute residuals and Jacobian for this interval
            // Use TRDIF if transition info is available and indicates transition occurred
            let (residuals, jacobian) = if let Some(tr) = transition {
                if tr.transition {
                    // Use full TRDIF with chain rule through XT derivatives
                    trdif_full(s1, s2, tr, 9.0, msq, re)
                } else {
                    // No transition occurred in this interval, use regular BLDIF
                    bldif(s1, s2, flow_type, msq, re)
                }
            } else {
                // No transition info, use regular BLDIF
                bldif(s1, s2, flow_type, msq, re)
            };

            // Store residuals (negated to match Newton iteration convention)
            // Note: bldif already negates residuals, so we use them directly
            self.rhs[i] = [residuals.res_third, residuals.res_mom, residuals.res_shape];

            // Extract 3x3 blocks from the 3x5 Jacobian matrices
            // The first 3 columns correspond to [ampl/ctau, theta, delta_star]
            //
            // VA[i] comes from vs2 (downstream station derivatives)
            // VB[i] comes from vs1 (upstream station derivatives)

            for eq in 0..3 {
                for var in 0..3 {
                    self.va[i][eq][var] = jacobian.vs2[eq][var];
                    self.vb[i][eq][var] = jacobian.vs1[eq][var];
                }
                // Store edge velocity coupling (column 3 = u variable)
                self.vs[i][eq] = jacobian.vs2[eq][3];
                // Store VS2 delta_star derivatives (column 2)
                self.vs_delta[i][eq] = jacobian.vs2[eq][2];
                // Store VS2 Ue derivatives (column 3)
                self.vs_ue[i][eq] = jacobian.vs2[eq][3];
                // Store VS1 delta_star derivatives (column 2) for forced changes
                self.vs1_delta[i][eq] = jacobian.vs1[eq][2];
                // Store VS1 Ue derivatives (column 3) for forced changes
                self.vs1_ue[i][eq] = jacobian.vs1[eq][3];
            }
        }
    }

    /// Build the Newton system automatically detecting transition intervals
    ///
    /// This method automatically detects transition intervals by checking if
    /// the upstream station is laminar and downstream is turbulent. When detected,
    /// it computes the transition location and uses TRDIF.
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval
    /// * `msq` - Mach number squared
    /// * `re` - Reference Reynolds number
    /// * `ncrit` - Critical N-factor for transition
    ///
    /// # Note
    /// This is a convenience method that detects transitions automatically.
    /// For more control, use `build_with_transition()` with explicit transition data.
    pub fn build_auto_transition(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
        ncrit: f64,
    ) {
        use rustfoil_bl::closures::trchek2_full;

        assert!(
            stations.len() >= 2,
            "Need at least 2 stations to build Newton system"
        );
        assert_eq!(
            stations.len(),
            self.n,
            "Station count must match system size"
        );
        assert_eq!(
            flow_types.len(),
            self.n - 1,
            "Need n-1 flow types for n stations"
        );

        // For each interval (between stations i-1 and i)
        for i in 1..self.n {
            let s1 = &stations[i - 1]; // Upstream station
            let s2 = &stations[i]; // Downstream station
            let flow_type = flow_types[i - 1];

            // Check if this is a transition interval
            let is_transition_interval = s1.is_laminar && s2.is_turbulent;

            // Compute residuals and Jacobian for this interval
            let (residuals, jacobian) = if is_transition_interval && ncrit > 0.0 {
                // Compute transition location using TRCHEK2
                let tr = trchek2_full(
                    s1.x,
                    s2.x,
                    s1.theta,
                    s2.theta,
                    s1.delta_star,
                    s2.delta_star,
                    s1.u,
                    s2.u,
                    s1.hk,
                    s2.hk,
                    s1.r_theta,
                    s2.r_theta,
                    s1.ampl,
                    ncrit,
                    None,
                    msq,
                    re,
                );

                if tr.transition {
                    // Use full TRDIF with chain rule
                    trdif_full(s1, s2, &tr, ncrit, msq, re)
                } else {
                    // No transition, use regular BLDIF
                    bldif(s1, s2, flow_type, msq, re)
                }
            } else {
                // Not a transition interval, use regular BLDIF
                bldif(s1, s2, flow_type, msq, re)
            };

            // Store residuals
            self.rhs[i] = [residuals.res_third, residuals.res_mom, residuals.res_shape];

            // Extract 3x3 blocks from the 3x5 Jacobian matrices
            for eq in 0..3 {
                for var in 0..3 {
                    self.va[i][eq][var] = jacobian.vs2[eq][var];
                    self.vb[i][eq][var] = jacobian.vs1[eq][var];
                }
                self.vs[i][eq] = jacobian.vs2[eq][3];
                self.vs_delta[i][eq] = jacobian.vs2[eq][2];
                self.vs_ue[i][eq] = jacobian.vs2[eq][3];
                self.vs1_delta[i][eq] = jacobian.vs1[eq][2];
                self.vs1_ue[i][eq] = jacobian.vs1[eq][3];
            }
        }
    }

    /// Build the Newton system with VM matrix for mass defect coupling
    ///
    /// This extended build method also constructs the VM matrix that couples
    /// mass defect changes at all stations to the BL equations. This is required
    /// for the full viscous-inviscid coupling as in XFOIL's SETBL.
    ///
    /// The VM matrix entries are computed as:
    /// ```text
    /// VM[i][j][k] = VS2[k,2]*D2_M[j] + VS2[k,3]*U2_M[j]
    ///             + VS1[k,2]*D1_M[j] + VS1[k,3]*U1_M[j]
    /// ```
    /// where:
    /// - D_M[j] = derivative of delta_star w.r.t. mass defect at station j
    /// - U_M[j] = derivative of Ue w.r.t. mass defect at station j (from DIJ)
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval
    /// * `msq` - Mach number squared
    /// * `re` - Reference Reynolds number
    /// * `dij` - DIJ matrix for Ue-mass defect coupling
    ///
    /// # Reference
    /// XFOIL xbl.f SETBL (lines 327-393)
    pub fn build_with_vm(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
        dij: &DMatrix<f64>,
    ) {
        // First build the basic system
        self.build(stations, flow_types, msq, re);

        // Clear VM matrix
        for i in 0..self.n {
            for j in 0..self.n {
                self.vm[i][j] = [0.0, 0.0, 0.0];
            }
        }

        // Build VM matrix for viscous-inviscid coupling
        // This implements XFOIL's SETBL VM construction
        //
        // For each equation station i (interval between i-1 and i):
        // For each influence station j:
        //   VM[i][j][k] = VS2[k,2]*D2_M[j] + VS2[k,3]*U2_M[j]
        //               + VS1[k,2]*D1_M[j] + VS1[k,3]*U1_M[j]
        //
        // where:
        //   D_M[j] = 1/Ue for j==i, 0 otherwise (delta_star sensitivity)
        //   U_M[j] = -DIJ[i,j] (Ue sensitivity from DIJ matrix)

        for i in 1..self.n {
            let s1 = &stations[i - 1];
            let s2 = &stations[i];
            let flow_type = flow_types[i - 1];

            // Recompute Jacobian (or use cached values)
            let (_, jacobian) = bldif(s1, s2, flow_type, msq, re);

            // VS1 and VS2 columns for delta_star and Ue
            let vs1_d: [f64; 3] = [jacobian.vs1[0][2], jacobian.vs1[1][2], jacobian.vs1[2][2]];
            let vs1_u: [f64; 3] = [jacobian.vs1[0][3], jacobian.vs1[1][3], jacobian.vs1[2][3]];
            let vs2_d: [f64; 3] = [jacobian.vs2[0][2], jacobian.vs2[1][2], jacobian.vs2[2][2]];
            let vs2_u: [f64; 3] = [jacobian.vs2[0][3], jacobian.vs2[1][3], jacobian.vs2[2][3]];

            // Derivatives of delta_star w.r.t. Ue (from mass defect relation)
            // delta_star = mass_defect / Ue
            // d(delta_star)/d(Ue) = -delta_star / Ue
            let d2_u2 = if s2.u.abs() > 1e-10 {
                -s2.delta_star / s2.u
            } else {
                0.0
            };
            let d1_u1 = if s1.u.abs() > 1e-10 {
                -s1.delta_star / s1.u
            } else {
                0.0
            };

            // For each influence station j
            for j in 0..self.n {
                // U_M[j] = derivative of Ue at current station w.r.t. mass defect at station j
                // This comes from the DIJ matrix
                let u2_m_j = if i < dij.nrows() && j < dij.ncols() {
                    -dij[(i, j)]
                } else {
                    0.0
                };
                let u1_m_j = if i > 0 && (i - 1) < dij.nrows() && j < dij.ncols() {
                    -dij[(i - 1, j)]
                } else {
                    0.0
                };

                // D_M[j] = derivative of delta_star w.r.t. mass defect at station j
                // d(delta_star)/d(mass) = 1/Ue at the same station
                // Plus the coupling through Ue change
                let d2_m_j = if j == i && s2.u.abs() > 1e-10 {
                    1.0 / s2.u + d2_u2 * u2_m_j
                } else {
                    d2_u2 * u2_m_j
                };
                let d1_m_j = if j == i - 1 && s1.u.abs() > 1e-10 {
                    1.0 / s1.u + d1_u1 * u1_m_j
                } else {
                    d1_u1 * u1_m_j
                };

                // Build VM entry for each equation
                for k in 0..3 {
                    self.vm[i][j][k] = vs2_d[k] * d2_m_j
                        + vs2_u[k] * u2_m_j
                        + vs1_d[k] * d1_m_j
                        + vs1_u[k] * u1_m_j;
                }
            }
        }
    }

    /// Build the Newton system with full XFOIL-style viscous-inviscid coupling
    ///
    /// This method implements XFOIL's SETBL algorithm including:
    /// 1. Building the basic block-tridiagonal system via bldif()
    /// 2. Building the VM matrix with VTI signs for proper Ue coupling
    /// 3. Computing and adding forced changes (DUE, DDS) to residuals
    ///
    /// The forced changes represent the mismatch between the current edge
    /// velocities and what UESET would compute from the current mass defect.
    /// This is crucial for Newton convergence when Ue has been modified
    /// outside the linearization.
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval
    /// * `msq` - Mach number squared
    /// * `re` - Reference Reynolds number
    /// * `dij` - DIJ matrix for Ue-mass defect coupling
    /// * `ue_inviscid` - Inviscid edge velocities (UINV in XFOIL)
    ///
    /// # Reference
    /// XFOIL xbl.f SETBL (lines 95-110 for USAV, 215-217 for forced changes)
    pub fn build_with_vm_full(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
        dij: &DMatrix<f64>,
        ue_inviscid: &[f64],
    ) {
        // Step 1: Save current Ue values (USAV in XFOIL)
        let ue_current: Vec<f64> = stations.iter().map(|s| s.u).collect();

        // Step 2: Compute what UESET would give (Ue from mass defect via DIJ)
        // XFOIL: UEDG(IBL,IS) = UINV(IBL,IS) + sum_j(-VTI*VTI*DIJ * MASS)
        let ue_from_mass = self.compute_ue_from_mass(stations, ue_inviscid, dij);

        // Step 3: Compute forced changes (mismatch between current Ue and UESET result)
        // XFOIL: DUE2 = UEDG(IBL,IS) - USAV(IBL,IS)
        // Note: In XFOIL, after calling UESET, USAV has new values and UEDG has old
        // So DUE = old_Ue - new_Ue = current - from_mass
        let due: Vec<f64> = (0..stations.len())
            .map(|i| ue_current[i] - ue_from_mass[i])
            .collect();

        // Step 4: Build basic system (calls bldif for each interval)
        // This populates VA, VB, RHS, and stores VS columns
        self.build(stations, flow_types, msq, re);

        // Step 5: Build VM matrix with VTI signs
        self.build_vm_with_vti(stations, flow_types, msq, re, dij);

        // Step 6: Add forced changes to residuals
        // XFOIL: VDEL(k,1,IV) = VSREZ(k) + VS1*DDS1 + VS2*DDS2 + VS1*DUE1 + VS2*DUE2
        self.add_forced_changes(stations, &due);
    }

    /// Compute edge velocities from mass defect using UESET formula
    ///
    /// XFOIL formula: Ue[i] = Ue_inv[i] + sum_j(-VTI[i]*VTI[j]*DIJ[i,j] * mass[j])
    ///
    /// # Arguments
    /// * `stations` - BL stations with mass_defect values
    /// * `ue_inviscid` - Inviscid edge velocities
    /// * `dij` - DIJ influence matrix
    ///
    /// # Returns
    /// Vector of computed edge velocities
    fn compute_ue_from_mass(
        &self,
        stations: &[BlStation],
        ue_inviscid: &[f64],
        dij: &DMatrix<f64>,
    ) -> Vec<f64> {
        let n = stations.len();
        let mut ue = vec![0.0; n];

        // Debug: print DIJ dimensions and sample panel indices (only if debug active)
        if rustfoil_bl::is_debug_active() {
            eprintln!("[DEBUG UESET] n_stations={}, DIJ dims={}x{}", 
                n, dij.nrows(), dij.ncols());
            eprintln!("[DEBUG UESET] panel_idx[0..5]={:?}", 
                stations.iter().take(5).map(|s| s.panel_idx).collect::<Vec<_>>());
        }

        for i in 0..n {
            let panel_i = stations[i].panel_idx;
            let mut dui = 0.0;
            for j in 0..n {
                let panel_j = stations[j].panel_idx;
                if panel_i < dij.nrows() && panel_j < dij.ncols() {
                    // UE_M = -VTI[i] * VTI[j] * DIJ[panel_i, panel_j]
                    let ue_m = -self.vti[i] * self.vti[j] * dij[(panel_i, panel_j)];
                    dui += ue_m * stations[j].mass_defect;
                }
            }
            ue[i] = if i < ue_inviscid.len() {
                ue_inviscid[i] + dui
            } else {
                dui
            };
        }

        ue
    }

    /// Build VM matrix with VTI signs for proper surface direction handling
    ///
    /// This implements XFOIL's VM construction with the correct VTI sign convention:
    /// U_M[j] = -VTI[i]*VTI[j]*DIJ[i,j]
    ///
    /// For single-surface mode (all VTI = +1), this reduces to U_M = -DIJ.
    /// For two-surface mode, the signs handle upper/lower velocity directions.
    fn build_vm_with_vti(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
        dij: &DMatrix<f64>,
    ) {
        // Clear VM matrix
        for i in 0..self.n {
            for j in 0..self.n {
                self.vm[i][j] = [0.0, 0.0, 0.0];
            }
        }

        for i in 1..self.n {
            let s1 = &stations[i - 1];
            let s2 = &stations[i];
            let flow_type = flow_types[i - 1];

            // Get Jacobian (recompute or use cached values from build())
            let (_, jacobian) = bldif(s1, s2, flow_type, msq, re);

            // VS1 and VS2 columns for delta_star and Ue
            let vs1_d: [f64; 3] = [jacobian.vs1[0][2], jacobian.vs1[1][2], jacobian.vs1[2][2]];
            let vs1_u: [f64; 3] = [jacobian.vs1[0][3], jacobian.vs1[1][3], jacobian.vs1[2][3]];
            let vs2_d: [f64; 3] = [jacobian.vs2[0][2], jacobian.vs2[1][2], jacobian.vs2[2][2]];
            let vs2_u: [f64; 3] = [jacobian.vs2[0][3], jacobian.vs2[1][3], jacobian.vs2[2][3]];

            // Derivatives of delta_star w.r.t. Ue (from mass defect relation)
            // delta_star = mass_defect / Ue
            // d(delta_star)/d(Ue) = -delta_star / Ue
            let d2_u2 = if s2.u.abs() > 1e-10 {
                -s2.delta_star / s2.u
            } else {
                0.0
            };
            let d1_u1 = if s1.u.abs() > 1e-10 {
                -s1.delta_star / s1.u
            } else {
                0.0
            };

            // Panel index for this station (maps BL station to DIJ matrix index)
            let panel_i = s2.panel_idx;
            let panel_im1 = s1.panel_idx;
            
            // For each influence station j
            for j in 0..self.n {
                // U_M[j] with VTI signs (XFOIL formula)
                // U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
                // I = panel index of current station, J = panel index of influence station
                let panel_j = stations[j].panel_idx;
                let u2_m_j = if panel_i < dij.nrows() && panel_j < dij.ncols() {
                    -self.vti[i] * self.vti[j] * dij[(panel_i, panel_j)]
                } else {
                    0.0
                };
                // U1_M for the upstream station (station i-1)
                let u1_m_j = if i > 1 && panel_im1 < dij.nrows() && panel_j < dij.ncols() {
                    -self.vti[i - 1] * self.vti[j] * dij[(panel_im1, panel_j)]
                } else {
                    0.0
                };

                // D_M[j] = derivative of delta_star w.r.t. mass defect at station j
                // d(delta_star)/d(mass) = 1/Ue at the same station
                // Plus the coupling through Ue change: d(ds)/d(Ue) * d(Ue)/d(mass)
                let d2_m_j = if j == i && s2.u.abs() > 1e-10 {
                    1.0 / s2.u + d2_u2 * u2_m_j
                } else {
                    d2_u2 * u2_m_j
                };
                let d1_m_j = if j == i - 1 && i > 1 && s1.u.abs() > 1e-10 {
                    1.0 / s1.u + d1_u1 * u1_m_j
                } else {
                    d1_u1 * u1_m_j
                };

                // Build VM entry for each equation
                // VM[k,j,i] = VS2[k,2]*D2_M[j] + VS2[k,3]*U2_M[j]
                //           + VS1[k,2]*D1_M[j] + VS1[k,3]*U1_M[j]
                for k in 0..3 {
                    self.vm[i][j][k] = vs2_d[k] * d2_m_j
                        + vs2_u[k] * u2_m_j
                        + vs1_d[k] * d1_m_j
                        + vs1_u[k] * u1_m_j;
                }
            }
        }
    }

    /// Add forced changes to residuals
    ///
    /// The forced changes account for the mismatch between current Ue values
    /// and what UESET would compute. This is essential for Newton convergence.
    ///
    /// XFOIL formula (xbl.f lines 351-355):
    /// ```text
    /// VDEL(k,1,IV) = VSREZ(k)
    ///              + VS1(k,3)*DUE1 + VS1(k,2)*DDS1
    ///              + VS2(k,3)*DUE2 + VS2(k,2)*DDS2
    /// ```
    ///
    /// where DDS = d(delta_star)/d(Ue) * DUE = -delta_star/Ue * DUE
    fn add_forced_changes(&mut self, stations: &[BlStation], due: &[f64]) {
        for i in 1..self.n {
            let s1 = &stations[i - 1];
            let s2 = &stations[i];

            // Compute d(delta_star)/d(Ue) = -delta_star / Ue
            let d2_u2 = if s2.u.abs() > 1e-10 {
                -s2.delta_star / s2.u
            } else {
                0.0
            };
            let d1_u1 = if i > 1 && s1.u.abs() > 1e-10 {
                -s1.delta_star / s1.u
            } else {
                0.0
            };

            // Forced changes in delta_star due to Ue mismatch
            let dds2 = d2_u2 * due[i];
            let dds1 = if i > 1 { d1_u1 * due[i - 1] } else { 0.0 };
            let due1 = if i > 1 { due[i - 1] } else { 0.0 };
            let due2 = due[i];

            // Add to residuals: VDEL += VS1*DDS1 + VS2*DDS2 + VS1*DUE1 + VS2*DUE2
            // Using stored VS columns
            for k in 0..3 {
                self.rhs[i][k] += self.vs1_delta[i][k] * dds1
                    + self.vs_delta[i][k] * dds2
                    + self.vs1_ue[i][k] * due1
                    + self.vs_ue[i][k] * due2;
            }
        }
    }

    /// Get maximum absolute residual (for convergence checking)
    ///
    /// This returns the infinity norm of the residual vector, which is
    /// commonly used to check Newton iteration convergence.
    ///
    /// # Returns
    /// The maximum absolute value of any residual component
    pub fn max_residual(&self) -> f64 {
        self.rhs
            .iter()
            .skip(1) // Skip index 0 (unused)
            .flat_map(|r| r.iter())
            .map(|&v| v.abs())
            .fold(0.0, f64::max)
    }

    /// Get the L2 norm of residuals (alternative convergence metric)
    ///
    /// # Returns
    /// The Euclidean norm of the residual vector
    pub fn residual_norm(&self) -> f64 {
        self.rhs
            .iter()
            .skip(1)
            .flat_map(|r| r.iter())
            .map(|&v| v * v)
            .sum::<f64>()
            .sqrt()
    }

    /// Check if the system has converged
    ///
    /// # Arguments
    /// * `tolerance` - Maximum allowed residual
    ///
    /// # Returns
    /// True if max_residual() < tolerance
    pub fn is_converged(&self, tolerance: f64) -> bool {
        self.max_residual() < tolerance
    }
}

/// Full coupled Newton system including inviscid interaction
///
/// This structure combines the BL block-tridiagonal system with the
/// inviscid coupling through the DIJ matrix. The DIJ matrix relates
/// changes in mass defect to changes in edge velocity.
///
/// # System Structure
/// The full coupled system solves for:
/// - BL primary variables at each station (ampl/ctau, theta, delta_star)
/// - Edge velocity updates through viscous-inviscid interaction
///
/// # XFOIL Reference
/// XFOIL xbl.f SETBL (line 21)
#[derive(Debug, Clone)]
pub struct CoupledNewtonSystem {
    /// BL system (block tridiagonal)
    pub bl: BlNewtonSystem,
    /// DIJ matrix for mass defect coupling
    /// ΔUe_i = Σ_j DIJ[i,j] × Δ(Ue×δ*)_j
    pub dij: DMatrix<f64>,
    /// Global system matrix (assembled from BL and coupling)
    pub matrix: DMatrix<f64>,
    /// Global RHS vector
    pub rhs: DVector<f64>,
}

impl CoupledNewtonSystem {
    /// Create a new coupled system
    ///
    /// # Arguments
    /// * `n_stations` - Number of BL stations
    /// * `dij` - Pre-computed DIJ influence matrix
    ///
    /// # Returns
    /// A new CoupledNewtonSystem with the BL system and DIJ matrix initialized
    pub fn new(n_stations: usize, dij: DMatrix<f64>) -> Self {
        // Each station has 3 BL variables
        let n_vars = 3 * n_stations;

        Self {
            bl: BlNewtonSystem::new(n_stations),
            dij,
            matrix: DMatrix::zeros(n_vars, n_vars),
            rhs: DVector::zeros(n_vars),
        }
    }

    /// Build the BL portion of the coupled system
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval
    /// * `msq` - Mach number squared
    /// * `re` - Reference Reynolds number
    pub fn build_bl_system(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
    ) {
        self.bl.build(stations, flow_types, msq, re);
    }

    /// Assemble the full coupled system matrix
    ///
    /// This combines the BL block-tridiagonal system with the inviscid
    /// coupling through the DIJ matrix. The assembled system can then
    /// be solved for the Newton update.
    ///
    /// # Note
    /// The BL system must be built first via `build_bl_system()`.
    ///
    /// # Reference
    /// XFOIL xbl.f SETBL (line 21)
    pub fn assemble(&mut self) {
        let n = self.bl.n;
        let n_vars = 3 * n;

        // Reset matrix and RHS
        self.matrix = DMatrix::zeros(n_vars, n_vars);
        self.rhs = DVector::zeros(n_vars);

        // Assemble block-tridiagonal BL system into global matrix
        for i in 1..n {
            let row_base = 3 * i;
            let col_base_diag = 3 * i; // Diagonal block (station i)
            let col_base_lower = 3 * (i - 1); // Lower diagonal (station i-1)

            // Insert VA[i] (diagonal block)
            for eq in 0..3 {
                for var in 0..3 {
                    self.matrix[(row_base + eq, col_base_diag + var)] = self.bl.va[i][eq][var];
                }
            }

            // Insert VB[i] (lower diagonal block)
            for eq in 0..3 {
                for var in 0..3 {
                    self.matrix[(row_base + eq, col_base_lower + var)] = self.bl.vb[i][eq][var];
                }
            }

            // Copy RHS
            for eq in 0..3 {
                self.rhs[row_base + eq] = self.bl.rhs[i][eq];
            }
        }

        // Add viscous-inviscid coupling through DIJ
        // The coupling adds terms: dR/dUe * dUe/d(mass_defect) * d(mass_defect)/d(primary_vars)
        // This is a simplified version - full implementation would include
        // the complete SETBL coupling logic
        //
        // TODO: Implement full SETBL coupling with DIJ matrix
    }

    /// Get maximum residual from the BL system
    pub fn max_residual(&self) -> f64 {
        self.bl.max_residual()
    }

    /// Check convergence of the coupled system
    pub fn is_converged(&self, tolerance: f64) -> bool {
        self.bl.is_converged(tolerance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_bl::equations::{blvar, FlowType};

    // =========================================================================
    // BlNewtonSystem Creation Tests
    // =========================================================================

    #[test]
    fn test_system_creation() {
        let system = BlNewtonSystem::new(50);

        assert_eq!(system.n, 50);
        assert_eq!(system.va.len(), 50);
        assert_eq!(system.vb.len(), 50);
        assert_eq!(system.rhs.len(), 50);
        assert_eq!(system.vs.len(), 50);

        // All values should be initialized to zero
        for i in 0..50 {
            for eq in 0..3 {
                for var in 0..3 {
                    assert_eq!(system.va[i][eq][var], 0.0);
                    assert_eq!(system.vb[i][eq][var], 0.0);
                }
                assert_eq!(system.rhs[i][eq], 0.0);
                assert_eq!(system.vs[i][eq], 0.0);
            }
        }
    }

    #[test]
    fn test_system_small() {
        // Minimum size system (2 stations)
        let system = BlNewtonSystem::new(2);

        assert_eq!(system.n, 2);
        assert_eq!(system.va.len(), 2);
    }

    // =========================================================================
    // Build Tests - Laminar Flow
    // =========================================================================

    #[test]
    fn test_build_laminar() {
        // Create a simple laminar BL with 5 stations
        let n = 5;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        // Create stations with increasing x and growing BL thickness
        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 + 0.1 * i as f64;
            station.u = 1.0 - 0.01 * i as f64; // Slight deceleration
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);
            station.ampl = 0.5 * i as f64; // Growing amplification

            // Compute secondary variables
            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        // Build the Newton system
        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Check that residuals are finite
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Residual [{},{}] should be finite, got {}",
                    i,
                    eq,
                    system.rhs[i][eq]
                );
            }
        }

        // Check that Jacobian blocks are non-trivial
        let mut has_nonzero_va = false;
        let mut has_nonzero_vb = false;
        for i in 1..n {
            for eq in 0..3 {
                for var in 0..3 {
                    if system.va[i][eq][var].abs() > 1e-15 {
                        has_nonzero_va = true;
                    }
                    if system.vb[i][eq][var].abs() > 1e-15 {
                        has_nonzero_vb = true;
                    }
                }
            }
        }
        assert!(has_nonzero_va, "VA blocks should have non-zero entries");
        assert!(has_nonzero_vb, "VB blocks should have non-zero entries");
    }

    #[test]
    fn test_build_turbulent() {
        // Create a turbulent BL
        let n = 5;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.3 + 0.1 * i as f64;
            station.u = 0.9 - 0.02 * i as f64;
            station.theta = 0.003 * (1.0 + 0.15 * i as f64);
            station.delta_star = 0.005 * (1.0 + 0.15 * i as f64);
            station.ctau = 0.1 + 0.01 * i as f64; // Shear stress coefficient
            station.is_laminar = false;
            station.is_turbulent = true;

            blvar(&mut station, FlowType::Turbulent, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Turbulent; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Check residuals are finite
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Turbulent residual [{},{}] should be finite",
                    i,
                    eq
                );
            }
        }
    }

    #[test]
    fn test_build_wake() {
        // Create a wake region
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 1.0 + 0.1 * i as f64;
            station.u = 0.5 + 0.02 * i as f64; // Recovering velocity
            station.theta = 0.01 * (1.0 + 0.05 * i as f64);
            station.delta_star = 0.02 * (1.0 + 0.05 * i as f64);
            station.ctau = 0.05 - 0.005 * i as f64; // Decaying shear
            station.is_wake = true;
            station.is_laminar = false;

            blvar(&mut station, FlowType::Wake, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Wake; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Verify residuals are finite for wake
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Wake residual [{},{}] should be finite",
                    i,
                    eq
                );
            }
        }
    }

    // =========================================================================
    // Jacobian Structure Tests
    // =========================================================================

    #[test]
    fn test_jacobian_block_structure_laminar() {
        // For laminar flow, the amplification equation should have
        // vs1[0][0] = -1 and vs2[0][0] = +1 (from bldif)
        let n = 3;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 * (i + 1) as f64;
            station.u = 1.0;
            station.theta = 0.001;
            station.delta_star = 0.0026;
            station.ampl = i as f64;

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Check structure of first equation (amplification)
        // VA[i][0][0] should be close to +1 (dres/dampl2)
        // VB[i][0][0] should be close to -1 (dres/dampl1)
        for i in 1..n {
            assert!(
                (system.va[i][0][0] - 1.0).abs() < 0.1,
                "VA[{}][0][0] should be ~1 for amplification, got {}",
                i,
                system.va[i][0][0]
            );
            assert!(
                (system.vb[i][0][0] - (-1.0)).abs() < 0.1,
                "VB[{}][0][0] should be ~-1 for amplification, got {}",
                i,
                system.vb[i][0][0]
            );
        }
    }

    // =========================================================================
    // max_residual Tests
    // =========================================================================

    #[test]
    fn test_max_residual_zero() {
        let system = BlNewtonSystem::new(5);

        // All zeros should give zero max residual
        assert_eq!(system.max_residual(), 0.0);
    }

    #[test]
    fn test_max_residual_finds_max() {
        let mut system = BlNewtonSystem::new(5);

        // Set some residuals
        system.rhs[1] = [0.1, -0.2, 0.05];
        system.rhs[2] = [0.01, 0.5, -0.3];
        system.rhs[3] = [-0.15, 0.02, 0.8]; // Maximum is 0.8

        let max_res = system.max_residual();
        assert!(
            (max_res - 0.8).abs() < 1e-10,
            "max_residual should be 0.8, got {}",
            max_res
        );
    }

    #[test]
    fn test_max_residual_skips_zero_index() {
        let mut system = BlNewtonSystem::new(3);

        // Put large value at index 0 (should be skipped)
        system.rhs[0] = [100.0, 100.0, 100.0];
        system.rhs[1] = [0.1, 0.2, 0.3];
        system.rhs[2] = [0.05, 0.15, 0.25];

        let max_res = system.max_residual();
        assert!(
            (max_res - 0.3).abs() < 1e-10,
            "max_residual should skip index 0, got {}",
            max_res
        );
    }

    #[test]
    fn test_residual_norm() {
        let mut system = BlNewtonSystem::new(3);

        system.rhs[1] = [3.0, 0.0, 0.0];
        system.rhs[2] = [0.0, 4.0, 0.0];

        // L2 norm should be sqrt(9 + 16) = 5
        let norm = system.residual_norm();
        assert!(
            (norm - 5.0).abs() < 1e-10,
            "residual_norm should be 5.0, got {}",
            norm
        );
    }

    #[test]
    fn test_is_converged() {
        let mut system = BlNewtonSystem::new(3);

        system.rhs[1] = [1e-6, 1e-7, 1e-8];
        system.rhs[2] = [1e-7, 1e-6, 1e-7];

        assert!(system.is_converged(1e-5), "Should converge with tolerance 1e-5");
        assert!(!system.is_converged(1e-7), "Should not converge with tolerance 1e-7");
    }

    // =========================================================================
    // Edge Velocity Coupling Tests
    // =========================================================================

    #[test]
    fn test_vs_coupling_nonzero() {
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 * (i + 1) as f64;
            station.u = 1.0 - 0.05 * i as f64; // Decelerating flow
            station.theta = 0.001 * (1.0 + 0.3 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.3 * i as f64);

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // VS coupling should have non-zero entries (edge velocity sensitivity)
        let mut has_nonzero_vs = false;
        for i in 1..n {
            for eq in 0..3 {
                if system.vs[i][eq].abs() > 1e-15 {
                    has_nonzero_vs = true;
                }
            }
        }

        assert!(
            has_nonzero_vs,
            "VS coupling terms should have non-zero entries for viscous-inviscid interaction"
        );
    }

    // =========================================================================
    // CoupledNewtonSystem Tests
    // =========================================================================

    #[test]
    fn test_coupled_system_creation() {
        use crate::dij::build_dij_matrix;

        let n = 5;
        let x: Vec<f64> = (0..n).map(|i| 0.1 * (i + 1) as f64).collect();
        let y: Vec<f64> = (0..n).map(|i| 0.01 * i as f64).collect();

        let dij = build_dij_matrix(&x, &y);
        let coupled = CoupledNewtonSystem::new(n, dij);

        assert_eq!(coupled.bl.n, n);
        assert_eq!(coupled.dij.nrows(), n);
        assert_eq!(coupled.dij.ncols(), n);
        assert_eq!(coupled.matrix.nrows(), 3 * n);
        assert_eq!(coupled.matrix.ncols(), 3 * n);
        assert_eq!(coupled.rhs.len(), 3 * n);
    }

    #[test]
    fn test_coupled_system_build_and_assemble() {
        use crate::dij::build_dij_matrix;

        let n = 4;
        let x: Vec<f64> = (0..n).map(|i| 0.1 * (i + 1) as f64).collect();
        let y: Vec<f64> = (0..n).map(|i| 0.01 * i as f64).collect();

        let dij = build_dij_matrix(&x, &y);
        let mut coupled = CoupledNewtonSystem::new(n, dij);

        // Create stations
        let re = 1e6;
        let msq = 0.0;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = x[i];
            station.u = 1.0;
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let flow_types = vec![FlowType::Laminar; n - 1];
        coupled.build_bl_system(&stations, &flow_types, msq, re);
        coupled.assemble();

        // Check that matrix has non-zero entries
        let mut has_nonzero = false;
        for i in 0..coupled.matrix.nrows() {
            for j in 0..coupled.matrix.ncols() {
                if coupled.matrix[(i, j)].abs() > 1e-15 {
                    has_nonzero = true;
                    break;
                }
            }
            if has_nonzero {
                break;
            }
        }

        assert!(has_nonzero, "Assembled matrix should have non-zero entries");
    }

    #[test]
    fn test_coupled_system_max_residual() {
        use crate::dij::build_dij_matrix;

        let n = 3;
        let x = vec![0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015];

        let dij = build_dij_matrix(&x, &y);
        let mut coupled = CoupledNewtonSystem::new(n, dij);

        // Manually set BL residuals
        coupled.bl.rhs[1] = [0.1, 0.2, 0.3];
        coupled.bl.rhs[2] = [0.05, 0.4, 0.15];

        let max_res = coupled.max_residual();
        assert!(
            (max_res - 0.4).abs() < 1e-10,
            "Coupled max_residual should be 0.4, got {}",
            max_res
        );
    }

    // =========================================================================
    // TRDIF (Transition) Tests
    // =========================================================================

    #[test]
    fn test_build_with_transition_no_transition() {
        // Test build_with_transition when no transition occurs (all None)
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 + 0.1 * i as f64;
            station.u = 1.0 - 0.01 * i as f64;
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);
            station.ampl = 0.5 * i as f64;

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        let transitions: Vec<Option<Trchek2FullResult>> = vec![None; n - 1];
        
        system.build_with_transition(&stations, &flow_types, msq, re, &transitions);

        // Should produce same results as regular build
        let mut system_regular = BlNewtonSystem::new(n);
        system_regular.build(&stations, &flow_types, msq, re);

        // Compare residuals
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    (system.rhs[i][eq] - system_regular.rhs[i][eq]).abs() < 1e-12,
                    "Residuals should match regular build when no transition"
                );
            }
        }
    }

    #[test]
    fn test_build_with_transition_interval() {
        // Test build_with_transition with a transition interval
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        // Create stations: first 2 laminar, last 2 turbulent
        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 + 0.1 * i as f64;
            station.u = 1.0 - 0.01 * i as f64;
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);
            
            if i < 2 {
                // Laminar stations
                station.ampl = 2.0 * i as f64; // Growing amplification
                station.is_laminar = true;
                station.is_turbulent = false;
                blvar(&mut station, FlowType::Laminar, msq, re);
            } else {
                // Turbulent stations
                station.ctau = 0.05;
                station.ampl = 0.0;
                station.is_laminar = false;
                station.is_turbulent = true;
                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
            stations.push(station);
        }

        // Create transition result for interval 1 (between stations 1 and 2)
        // This is where transition occurs (laminar -> turbulent)
        let mut tr = Trchek2FullResult::default();
        tr.transition = true;
        tr.xt = 0.25; // Transition at x=0.25 (between 0.2 and 0.3)
        tr.ampl2 = 9.0; // N-crit
        tr.wf1 = 0.5;
        tr.wf2 = 0.5;
        // Set some non-zero derivatives for testing
        tr.tt_t1 = 0.5;
        tr.tt_t2 = 0.5;
        tr.dt_d1 = 0.5;
        tr.dt_d2 = 0.5;
        tr.ut_u1 = 0.5;
        tr.ut_u2 = 0.5;

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar, FlowType::Turbulent, FlowType::Turbulent];
        let transitions: Vec<Option<Trchek2FullResult>> = vec![None, Some(tr), None];
        
        system.build_with_transition(&stations, &flow_types, msq, re, &transitions);

        // Check that residuals are finite
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Residual [{},{}] should be finite with transition, got {}",
                    i, eq, system.rhs[i][eq]
                );
            }
        }

        // The transition interval (i=2) should have different Jacobian structure
        // due to TRDIF combining laminar and turbulent parts
        let mut has_nonzero_va = false;
        let mut has_nonzero_vb = false;
        for eq in 0..3 {
            for var in 0..3 {
                if system.va[2][eq][var].abs() > 1e-15 {
                    has_nonzero_va = true;
                }
                if system.vb[2][eq][var].abs() > 1e-15 {
                    has_nonzero_vb = true;
                }
            }
        }
        assert!(has_nonzero_va, "VA at transition interval should have non-zero entries");
        assert!(has_nonzero_vb, "VB at transition interval should have non-zero entries");
    }

    #[test]
    fn test_build_auto_transition() {
        // Test automatic transition detection
        let n = 5;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;
        let ncrit = 9.0;

        // Create stations with laminar -> turbulent transition
        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.05 + 0.1 * i as f64;
            station.u = 1.0 - 0.02 * i as f64;
            station.theta = 0.0005 * (1.0 + 0.3 * i as f64);
            station.delta_star = 0.0013 * (1.0 + 0.3 * i as f64);
            
            if i < 3 {
                // Laminar stations
                station.ampl = 3.0 * i as f64;
                station.is_laminar = true;
                station.is_turbulent = false;
                blvar(&mut station, FlowType::Laminar, msq, re);
            } else {
                // Turbulent stations
                station.ctau = 0.04;
                station.ampl = 0.0;
                station.is_laminar = false;
                station.is_turbulent = true;
                blvar(&mut station, FlowType::Turbulent, msq, re);
            }
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![
            FlowType::Laminar, 
            FlowType::Laminar, 
            FlowType::Turbulent, // Transition interval
            FlowType::Turbulent,
        ];
        
        system.build_auto_transition(&stations, &flow_types, msq, re, ncrit);

        // Check that residuals are finite
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Residual [{},{}] should be finite with auto-transition, got {}",
                    i, eq, system.rhs[i][eq]
                );
            }
        }

        // Check that Jacobian blocks are non-trivial
        let mut total_nonzero = 0;
        for i in 1..n {
            for eq in 0..3 {
                for var in 0..3 {
                    if system.va[i][eq][var].abs() > 1e-15 {
                        total_nonzero += 1;
                    }
                }
            }
        }
        assert!(total_nonzero > 0, "VA blocks should have non-zero entries with auto-transition");
    }

    #[test]
    fn test_build_auto_transition_ncrit_zero_disables() {
        // When ncrit=0, auto transition should behave like regular build
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 + 0.1 * i as f64;
            station.u = 1.0;
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);
            station.ampl = i as f64;

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let mut system_auto = BlNewtonSystem::new(n);
        let mut system_regular = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        
        system_auto.build_auto_transition(&stations, &flow_types, msq, re, 0.0);
        system_regular.build(&stations, &flow_types, msq, re);

        // Results should match when ncrit=0 (TRDIF disabled)
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    (system_auto.rhs[i][eq] - system_regular.rhs[i][eq]).abs() < 1e-12,
                    "Residuals should match regular build when ncrit=0"
                );
            }
        }
    }
}
