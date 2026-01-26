//! Global Newton system for viscous-inviscid coupling
//!
//! This module implements XFOIL's full global Newton system that couples
//! upper and lower boundary layer surfaces through the mass defect influence
//! matrix (DIJ). Unlike the separate-surface approach, this properly handles
//! cross-surface coupling where mass defect changes on one surface affect
//! edge velocities on the other.
//!
//! # System Structure
//!
//! The global system combines upper and lower surface BL equations:
//! ```text
//! |  A  |  |  .  |  |  .  |    d       R
//! |  B  A  |  .  |  |  .  |    d       R
//! |  |  B  A  .  |  |  .  |    d       R
//! |  .  .  .  .  |  |  .  |    d   =   R
//! |  |  |  |  B  A  |  .  |    d       R
//! |  |  Z  |  |  B  A  .  |    d       R  <- VZ at upper TE
//! |  .  .  .  .  .  .  .  |    d       R
//! |  |  |  |  |  |  |  B  A    d       R
//!
//! A, B, Z = 3x3 blocks
//! |       = 3x1 VM columns (mass influence)
//! d       = 3x1 unknowns (dCtau/dAmpl, dTheta, dMass)
//! R       = 3x1 residuals
//! ```
//!
//! # Key Features
//!
//! - Full VM matrix with cross-surface coupling through DIJ
//! - VZ block at trailing edge coupling upper and lower surfaces
//! - VTI sign convention: +1 upper, -1 lower
//! - VACCEL threshold for sparse elimination
//!
//! # XFOIL Reference
//! - SETBL: xbl.f line 21 (builds global system)
//! - BLSOLV: xsolve.f line 283 (solves the system)

use nalgebra::DMatrix;
use rustfoil_bl::closures::Trchek2FullResult;
use rustfoil_bl::equations::{bldif, bldif_full_simi, blvar, trdif_full, FlowType};
use rustfoil_bl::state::BlStation;

/// Global Newton system for coupled upper/lower surface BL iteration
///
/// This structure holds the complete Newton system for full viscous-inviscid
/// coupling. The system combines both surfaces into a single matrix solve,
/// allowing proper coupling of mass defect effects across surfaces.
///
/// # Station Ordering
/// - Stations 0..n_upper: Upper surface (stagnation to TE to wake)
/// - Stations n_upper..nsys: Lower surface (stagnation to TE to wake)
///
/// # Index Convention
/// Global index IV maps to local surface coordinates via `from_global()`.
/// The VTI array contains surface direction signs (+1 upper, -1 lower).
#[derive(Debug, Clone)]
pub struct GlobalNewtonSystem {
    /// Total number of system stations (upper + lower, minus shared stagnation)
    pub nsys: usize,
    /// Number of stations on upper surface (including stagnation and wake)
    pub n_upper: usize,
    /// Number of stations on lower surface (including stagnation and wake)
    pub n_lower: usize,
    /// BL station index of trailing edge on upper surface
    pub iblte_upper: usize,
    /// BL station index of trailing edge on lower surface
    pub iblte_lower: usize,

    // === Block matrices (indexed by global IV) ===
    /// Diagonal blocks VA[IV] (3x3) - derivatives w.r.t. downstream station
    pub va: Vec<[[f64; 3]; 3]>,
    /// Lower diagonal blocks VB[IV] (3x3) - derivatives w.r.t. upstream station
    pub vb: Vec<[[f64; 3]; 3]>,
    /// TE coupling block VZ (3x2) - couples upper TE to lower wake start
    pub vz: [[f64; 3]; 2],
    /// Mass influence matrix VM[IV][JV][k] - sensitivity of eq k at IV to mass at JV
    /// Dimensions: nsys x nsys x 3
    pub vm: Vec<Vec<[f64; 3]>>,
    /// RHS residual / solution vector (3 per station, 2 for RHS/solution)
    /// VDEL[IV] = [res_third, res_mom, res_shape]
    pub vdel: Vec<[f64; 3]>,

    // === Index mappings ===
    /// Surface direction signs: +1.0 for upper, -1.0 for lower
    pub vti: Vec<f64>,
    /// Panel index for each global station (maps to DIJ matrix)
    pub panel_idx: Vec<usize>,
    /// Flag for stations at or after TE (where wake equations apply)
    pub is_wake: Vec<bool>,

    /// VS1 columns for delta_star derivatives (stored for forced changes)
    pub vs1_delta: Vec<[f64; 3]>,
    /// VS1 columns for Ue derivatives (stored for forced changes)
    pub vs1_ue: Vec<[f64; 3]>,
    /// VS2 columns for delta_star derivatives
    pub vs2_delta: Vec<[f64; 3]>,
    /// VS2 columns for Ue derivatives
    pub vs2_ue: Vec<[f64; 3]>,
}

impl GlobalNewtonSystem {
    /// Create a new global Newton system
    ///
    /// # Arguments
    /// * `n_upper` - Number of stations on upper surface
    /// * `n_lower` - Number of stations on lower surface
    /// * `iblte_upper` - TE station index on upper surface
    /// * `iblte_lower` - TE station index on lower surface
    ///
    /// # Note
    /// The total system size is n_upper + n_lower - 2 because:
    /// - Stagnation point is shared (don't double-count)
    /// - First interval on each surface is boundary condition (SIMI)
    pub fn new(n_upper: usize, n_lower: usize, iblte_upper: usize, iblte_lower: usize) -> Self {
        // Total system size: upper + lower, with adjustments
        // Following XFOIL: NSYS = NBL(1) + NBL(2) - 2
        // This accounts for shared stagnation and boundary conditions
        let nsys = n_upper + n_lower - 2;

        Self {
            nsys,
            n_upper,
            n_lower,
            iblte_upper,
            iblte_lower,
            va: vec![[[0.0; 3]; 3]; nsys + 1],
            vb: vec![[[0.0; 3]; 3]; nsys + 1],
            vz: [[0.0; 3]; 2],
            vm: vec![vec![[0.0; 3]; nsys + 1]; nsys + 1],
            vdel: vec![[0.0; 3]; nsys + 1],
            vti: vec![1.0; nsys + 1],
            panel_idx: vec![0; nsys + 1],
            is_wake: vec![false; nsys + 1],
            vs1_delta: vec![[0.0; 3]; nsys + 1],
            vs1_ue: vec![[0.0; 3]; nsys + 1],
            vs2_delta: vec![[0.0; 3]; nsys + 1],
            vs2_ue: vec![[0.0; 3]; nsys + 1],
        }
    }

    /// Map (surface, local_idx) to global system index
    ///
    /// # Arguments
    /// * `surface` - 0 for upper, 1 for lower
    /// * `local_ibl` - Local BL station index on that surface
    ///
    /// # Returns
    /// Global system index IV
    ///
    /// # Station Mapping
    /// - Upper surface: IV = local_ibl (1-indexed, so IBL=1 → IV=1)
    /// - Lower surface: IV = n_upper - 1 + local_ibl
    pub fn to_global(&self, surface: usize, local_ibl: usize) -> usize {
        if surface == 0 {
            // Upper surface: direct mapping
            local_ibl
        } else {
            // Lower surface: offset by upper count, minus shared stagnation
            self.n_upper - 1 + local_ibl
        }
    }

    /// Map global index to (surface, local_idx)
    ///
    /// # Arguments
    /// * `iv` - Global system index
    ///
    /// # Returns
    /// (surface, local_ibl) where surface is 0=upper, 1=lower
    pub fn from_global(&self, iv: usize) -> (usize, usize) {
        if iv < self.n_upper {
            // Upper surface
            (0, iv)
        } else {
            // Lower surface
            (1, iv - (self.n_upper - 1))
        }
    }

    /// Build the global Newton system from BL stations
    ///
    /// This implements XFOIL's SETBL algorithm including:
    /// 1. Building VA, VB blocks via bldif() for each interval
    /// 2. Building the full VM matrix with cross-surface coupling
    /// 3. Computing forced changes and adding to residuals
    /// 4. Setting up VZ block at trailing edge
    ///
    /// # Arguments
    /// * `upper_stations` - Upper surface BL stations (pre-computed secondary vars)
    /// * `lower_stations` - Lower surface BL stations
    /// * `upper_flow_types` - Flow types for upper surface intervals
    /// * `lower_flow_types` - Flow types for lower surface intervals
    /// * `upper_transitions` - Optional transition results for upper intervals
    /// * `lower_transitions` - Optional transition results for lower intervals
    /// * `dij` - Global DIJ matrix for mass defect coupling
    /// * `ue_inviscid_upper` - Inviscid edge velocities (upper)
    /// * `ue_inviscid_lower` - Inviscid edge velocities (lower)
    /// * `msq` - Mach number squared
    /// * `re` - Reynolds number
    ///
    /// # XFOIL Reference
    /// XFOIL xbl.f SETBL (lines 21-500)
    #[allow(clippy::too_many_arguments)]
    pub fn build_global_system(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        upper_flow_types: &[FlowType],
        lower_flow_types: &[FlowType],
        upper_transitions: &[Option<Trchek2FullResult>],
        lower_transitions: &[Option<Trchek2FullResult>],
        dij: &DMatrix<f64>,
        ue_inviscid_upper: &[f64],
        ue_inviscid_lower: &[f64],
        msq: f64,
        re: f64,
    ) {
        // Reset arrays
        for iv in 0..=self.nsys {
            self.va[iv] = [[0.0; 3]; 3];
            self.vb[iv] = [[0.0; 3]; 3];
            self.vdel[iv] = [0.0; 3];
            self.vs1_delta[iv] = [0.0; 3];
            self.vs1_ue[iv] = [0.0; 3];
            self.vs2_delta[iv] = [0.0; 3];
            self.vs2_ue[iv] = [0.0; 3];
            for jv in 0..=self.nsys {
                self.vm[iv][jv] = [0.0; 3];
            }
        }
        self.vz = [[0.0; 3]; 2];

        // === Setup index mappings and VTI signs ===
        self.setup_index_mappings(upper_stations, lower_stations);

        // === Step 1: Save current Ue values (USAV) for forced changes ===
        let ue_current_upper: Vec<f64> = upper_stations.iter().map(|s| s.u).collect();
        let ue_current_lower: Vec<f64> = lower_stations.iter().map(|s| s.u).collect();

        // === Step 2: Compute Ue from mass defect (UESET) ===
        // IMPORTANT: Sum over BOTH surfaces to get proper mass coupling
        let ue_from_mass_upper =
            self.compute_ue_from_mass_both(upper_stations, lower_stations, ue_inviscid_upper, dij, 0);
        let ue_from_mass_lower =
            self.compute_ue_from_mass_both(upper_stations, lower_stations, ue_inviscid_lower, dij, 1);

        // === Step 3: Build equations for upper surface ===
        self.build_surface_equations(
            upper_stations,
            upper_flow_types,
            upper_transitions,
            &ue_current_upper,
            &ue_from_mass_upper,
            msq,
            re,
            0, // surface = upper
        );

        // === Step 4: Build equations for lower surface ===
        self.build_surface_equations(
            lower_stations,
            lower_flow_types,
            lower_transitions,
            &ue_current_lower,
            &ue_from_mass_lower,
            msq,
            re,
            1, // surface = lower
        );

        // === Step 5: Build VM matrix with cross-surface coupling ===
        self.build_vm_global(upper_stations, lower_stations, dij, msq, re);

        // === Step 6: Build VZ block at trailing edge ===
        self.build_vz_block(upper_stations, lower_stations);

        // === Step 7: Add forced changes to residuals ===
        // XFOIL formula (xbl.f lines 383-387):
        // VDEL(k,1,IV) = VSREZ(k) + (VS1(k,4)*DUE1 + VS1(k,3)*DDS1)
        //                        + (VS2(k,4)*DUE2 + VS2(k,3)*DDS2) + ...
        // 
        // DUE = current_Ue - ue_from_mass (change needed to match mass coupling)
        // DDS = D_U * DUE (delta_star change induced by Ue change)
        self.add_forced_changes(
            upper_stations,
            lower_stations,
            &ue_current_upper,
            &ue_current_lower,
            &ue_from_mass_upper,
            &ue_from_mass_lower,
        );
        
        // === Step 8: Negate residuals for Newton sign convention ===
        // Our bldif computes res = -F (negative of equation value)
        // XFOIL's VSREZ = -REZC, -REZT, -REZH (also negative)
        // But the solver expects vdel = -F for Newton: J*delta = -F, x += delta
        // So we need to ensure the sign is correct. Actually XFOIL does NOT negate
        // after storing in VSREZ, so we shouldn't either.
    }

    /// Emit debug event for the built Newton system
    ///
    /// This dumps VA, VB, VDEL, and VM samples matching XFOIL's DBGSETBLSYSTEM.
    /// Should be called after `build_global_system` completes.
    ///
    /// # Arguments
    /// * `iteration` - Viscous iteration number
    pub fn emit_setbl_debug(&self, iteration: usize) {
        if !rustfoil_bl::is_debug_active() {
            return;
        }

        let nout = self.nsys.min(20);

        // Extract VA blocks (3 equations x 2 variables)
        // XFOIL VA is (3,2,IVX) - we have (3,3) but only use first 2 columns
        let mut va_blocks: Vec<[[f64; 2]; 3]> = Vec::with_capacity(nout);
        for iv in 1..=nout {
            let mut block = [[0.0; 2]; 3];
            for k in 0..3 {
                block[k][0] = self.va[iv][k][0];
                block[k][1] = self.va[iv][k][1];
            }
            va_blocks.push(block);
        }

        // Extract VB blocks
        let mut vb_blocks: Vec<[[f64; 2]; 3]> = Vec::with_capacity(nout);
        for iv in 1..=nout {
            let mut block = [[0.0; 2]; 3];
            for k in 0..3 {
                block[k][0] = self.vb[iv][k][0];
                block[k][1] = self.vb[iv][k][1];
            }
            vb_blocks.push(block);
        }

        // Extract VDEL (RHS)
        let mut vdel_rhs: Vec<[f64; 3]> = Vec::with_capacity(nout);
        for iv in 1..=nout {
            vdel_rhs.push(self.vdel[iv]);
        }

        // Extract VM diagonal
        let mut vm_diag: Vec<[f64; 3]> = Vec::with_capacity(nout);
        for iv in 1..=nout {
            vm_diag.push(self.vm[iv][iv]);
        }

        // Extract VM row 1
        let mut vm_row1: Vec<[f64; 3]> = Vec::with_capacity(nout);
        for jv in 1..=nout {
            vm_row1.push(self.vm[1][jv]);
        }

        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::setbl_system(
            iteration,
            self.nsys,
            va_blocks,
            vb_blocks,
            vdel_rhs,
            vm_diag,
            vm_row1,
        ));
    }

    /// Setup index mappings and VTI signs
    fn setup_index_mappings(&mut self, upper_stations: &[BlStation], lower_stations: &[BlStation]) {
        // Upper surface: IV = 1..n_upper, VTI = +1
        for ibl in 1..self.n_upper {
            let iv = self.to_global(0, ibl);
            if iv <= self.nsys {
                self.vti[iv] = 1.0;
                self.panel_idx[iv] = upper_stations[ibl].panel_idx;
                self.is_wake[iv] = ibl > self.iblte_upper;
            }
        }

        // Lower surface: IV = n_upper..nsys, VTI = -1
        for ibl in 1..self.n_lower {
            let iv = self.to_global(1, ibl);
            if iv <= self.nsys {
                self.vti[iv] = -1.0;
                self.panel_idx[iv] = lower_stations[ibl].panel_idx;
                self.is_wake[iv] = ibl > self.iblte_lower;
            }
        }
    }

    /// Compute Ue from mass defect using UESET formula
    ///
    /// XFOIL formula: Ue[i] = Ue_inv[i] + sum_j(-VTI[i]*VTI[j]*DIJ[i,j] * mass[j])
    /// NOTE: Sum is over BOTH surfaces, not just the current one!
    fn compute_ue_from_mass_both(
        &self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        ue_inviscid: &[f64],
        dij: &DMatrix<f64>,
        surface: usize,
    ) -> Vec<f64> {
        let stations = if surface == 0 { upper_stations } else { lower_stations };
        let n = stations.len();
        let mut ue = vec![0.0; n];

        for i in 0..n {
            let panel_i = stations[i].panel_idx;
            let vti_i = if surface == 0 { 1.0 } else { -1.0 };

            let mut dui = 0.0;

            // Sum over upper surface (VTI = +1)
            // XFOIL starts from JBL=2 (first BL station, skipping stagnation)
            for j in 1..upper_stations.len() {
                let panel_j = upper_stations[j].panel_idx;
                let vti_j = 1.0;  // Upper surface

                if panel_i < dij.nrows() && panel_j < dij.ncols() {
                    let ue_m = -vti_i * vti_j * dij[(panel_i, panel_j)];
                    dui += ue_m * upper_stations[j].mass_defect;
                }
            }

            // Sum over lower surface (VTI = -1)
            // XFOIL starts from JBL=2 (first BL station, skipping stagnation)
            for j in 1..lower_stations.len() {
                let panel_j = lower_stations[j].panel_idx;
                let vti_j = -1.0;  // Lower surface

                if panel_i < dij.nrows() && panel_j < dij.ncols() {
                    let ue_m = -vti_i * vti_j * dij[(panel_i, panel_j)];
                    dui += ue_m * lower_stations[j].mass_defect;
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

    /// Build BL equations for one surface
    #[allow(clippy::too_many_arguments)]
    fn build_surface_equations(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        transitions: &[Option<Trchek2FullResult>],
        _ue_current: &[f64],
        _ue_from_mass: &[f64],
        msq: f64,
        re: f64,
        surface: usize, // 0 = upper, 1 = lower
    ) {
        let n = stations.len();

        // For each interval (between stations i-1 and i)
        for i in 1..n {
            let iv = self.to_global(surface, i);
            if iv > self.nsys {
                continue;
            }

            let s2 = &stations[i];
            let flow_type = flow_types.get(i - 1).copied().unwrap_or(FlowType::Laminar);
            let transition = transitions.get(i - 1).and_then(|t| t.as_ref());
            
            // Compute residuals and Jacobian
            // For first interval (i==1), use similarity-specific bldif with XFOIL's
            // fixed log values (XLOG=1, ULOG=1, TLOG=0, HLOG=0) and DDLOG=0
            // which zeros out the log-derivative Jacobian terms
            let (residuals, jacobian) = if i == 1 {
                bldif_full_simi(s2, flow_type, msq, re)
            } else if let Some(tr) = transition {
                if tr.transition {
                    trdif_full(&stations[i - 1], s2, tr, msq, re)
                } else {
                    bldif(&stations[i - 1], s2, flow_type, msq, re)
                }
            } else {
                bldif(&stations[i - 1], s2, flow_type, msq, re)
            };


            // Store residuals
            self.vdel[iv] = [residuals.res_third, residuals.res_mom, residuals.res_shape];
            

            // Extract VA and VB blocks from Jacobian
            // IMPORTANT: XFOIL's VA only has 2 columns: [ampl/ctau, theta]
            // The delta_star/mass derivatives go ENTIRELY in VM, not VA!
            // 
            // XFOIL (xbl.f lines 340-341, 370-371, 400-401):
            //   VA(k,1,IV) = VS2(k,1)  -- ∂Fk/∂(ampl or ctau)
            //   VA(k,2,IV) = VS2(k,2)  -- ∂Fk/∂theta
            //
            // For similarity station (i==1), XFOIL combines: VS2 = VS1 + VS2, VS1 = 0
            if i == 1 {
                // Similarity station: XFOIL combines VS2 = VS1 + VS2, VS1 = 0
                // BUT when s1=s2 (both point to same station), vs1 and vs2 from bldif
                // are IDENTICAL (both compute derivatives w.r.t. same station).
                // So we should NOT add them - just use vs2 directly!
                //
                // In XFOIL's BLDIF, VS1 and VS2 are computed with conceptually different
                // "upstream" vs "downstream" formulas even when using same state variables.
                // For SIMI, the final combined VS2 is what we want.
                //
                // Since our bldif_full_simi returns vs1=vs2 (same derivatives), we just
                // use vs2 without adding vs1 to avoid 2x factor.
                for eq in 0..3 {
                    // Only columns 0,1 go in VA (ampl/ctau and theta)
                    // For SIMI, DON'T add vs1+vs2 since they're identical
                    self.va[iv][eq][0] = jacobian.vs2[eq][0];
                    self.va[iv][eq][1] = jacobian.vs2[eq][1];
                    self.va[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    self.vb[iv][eq][0] = 0.0;
                    self.vb[iv][eq][1] = 0.0;
                    self.vb[iv][eq][2] = 0.0;
                    
                    // Store delta_star and Ue columns for VM construction
                    // Same logic - just use vs2, not vs1+vs2
                    self.vs1_delta[iv][eq] = 0.0;
                    self.vs1_ue[iv][eq] = 0.0;
                    self.vs2_delta[iv][eq] = jacobian.vs2[eq][2];
                    self.vs2_ue[iv][eq] = jacobian.vs2[eq][3];
                }
            } else {
                for eq in 0..3 {
                    // Only columns 0,1 go in VA and VB (ampl/ctau and theta)
                    self.va[iv][eq][0] = jacobian.vs2[eq][0];
                    self.va[iv][eq][1] = jacobian.vs2[eq][1];
                    self.va[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    self.vb[iv][eq][0] = jacobian.vs1[eq][0];
                    self.vb[iv][eq][1] = jacobian.vs1[eq][1];
                    self.vb[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    
                    // Store delta_star and Ue columns for VM construction
                    self.vs1_delta[iv][eq] = jacobian.vs1[eq][2];
                    self.vs1_ue[iv][eq] = jacobian.vs1[eq][3];
                    self.vs2_delta[iv][eq] = jacobian.vs2[eq][2];
                    self.vs2_ue[iv][eq] = jacobian.vs2[eq][3];
                }
            }
            
        }

        // NOTE: XFOIL does NOT zero the SIMI residuals - they're computed normally
        // with fixed log values (XLOG=1, ULOG=BULE, TLOG=0, HLOG=0).
        // The similarity condition is handled by combining Jacobians (VS2 = VS1 + VS2)
        // which is done in the build_surface_equations loop for i==1.
    }

    /// Build global VM matrix with cross-surface coupling
    ///
    /// VM[IV][JV][k] = sensitivity of equation k at station IV to mass defect at station JV
    ///
    /// The key insight is that mass changes on one surface affect Ue on BOTH surfaces
    /// through the DIJ matrix, so we need full cross-surface coupling.
    fn build_vm_global(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        dij: &DMatrix<f64>,
        _msq: f64,
        _re: f64,
    ) {
        // For each equation station IV
        for iv in 1..=self.nsys {
            let (surface_i, ibl_i) = self.from_global(iv);
            let stations_i = if surface_i == 0 {
                upper_stations
            } else {
                lower_stations
            };

            if ibl_i == 0 || ibl_i >= stations_i.len() {
                continue;
            }

            let s1 = &stations_i[ibl_i - 1];
            let s2 = &stations_i[ibl_i];
            let panel_i = s2.panel_idx;
            let panel_im1 = s1.panel_idx;

            // Derivatives of delta_star w.r.t. Ue
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

            // For each influence station JV (on BOTH surfaces)
            for jv in 1..=self.nsys {
                let (surface_j, ibl_j) = self.from_global(jv);
                let stations_j = if surface_j == 0 {
                    upper_stations
                } else {
                    lower_stations
                };

                if ibl_j >= stations_j.len() {
                    continue;
                }

                let panel_j = stations_j[ibl_j].panel_idx;

                // U_M with VTI signs (XFOIL formula)
                // U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
                let u2_m_j = if panel_i < dij.nrows() && panel_j < dij.ncols() {
                    -self.vti[iv] * self.vti[jv] * dij[(panel_i, panel_j)]
                } else {
                    0.0
                };

                let u1_m_j = if ibl_i > 1 && panel_im1 < dij.nrows() && panel_j < dij.ncols() {
                    -self.vti[iv - 1].max(self.vti[iv]) * self.vti[jv] * dij[(panel_im1, panel_j)]
                } else {
                    0.0
                };

                // D_M[j] = derivative of delta_star w.r.t. mass defect at station j
                // Self-term: d(delta_star)/d(mass) = 1/Ue at same station
                let same_station_i = surface_i == surface_j && ibl_j == ibl_i;
                let same_station_im1 = surface_i == surface_j && ibl_j == ibl_i - 1;

                let d2_m_j = if same_station_i && s2.u.abs() > 1e-10 {
                    1.0 / s2.u + d2_u2 * u2_m_j
                } else {
                    d2_u2 * u2_m_j
                };

                let d1_m_j = if same_station_im1 && ibl_i > 1 && s1.u.abs() > 1e-10 {
                    1.0 / s1.u + d1_u1 * u1_m_j
                } else {
                    d1_u1 * u1_m_j
                };

                // Build VM entry for each equation
                for k in 0..3 {
                    self.vm[iv][jv][k] = self.vs2_delta[iv][k] * d2_m_j
                        + self.vs2_ue[iv][k] * u2_m_j
                        + self.vs1_delta[iv][k] * d1_m_j
                        + self.vs1_ue[iv][k] * u1_m_j;
                }
            }
        }
    }

    /// Build VZ block at trailing edge
    ///
    /// The VZ block couples the upper TE to the lower wake start.
    /// This handles the discontinuity where upper and lower surfaces meet.
    fn build_vz_block(&mut self, upper_stations: &[BlStation], lower_stations: &[BlStation]) {
        // VZ block couples upper TE (IV = iblte_upper) to lower wake (JV = n_upper + iblte_lower)
        let iv_upper_te = self.to_global(0, self.iblte_upper);
        let jv_lower_wake = self.to_global(1, self.iblte_lower + 1);

        if iv_upper_te > self.nsys || jv_lower_wake > self.nsys {
            return;
        }

        // Get stations at TE
        if self.iblte_upper < upper_stations.len() && self.iblte_lower + 1 < lower_stations.len() {
            let _s_upper_te = &upper_stations[self.iblte_upper];
            let _s_lower_wake = &lower_stations[self.iblte_lower + 1];

            // VZ stores coupling coefficients
            // VZ[k,0] = theta coupling
            // VZ[k,1] = ctau coupling
            // In XFOIL: VZ(k,1) = VS1(k,2), VZ(k,2) = VS1(k,1)
            for k in 0..3 {
                self.vz[k][0] = self.vs1_delta[jv_lower_wake][k];
                self.vz[k][1] = self.vs1_ue[jv_lower_wake][k];
            }
        }
    }

    /// Add forced changes to residuals
    ///
    /// The forced changes account for the mismatch between current Ue values
    /// and what UESET would compute from current mass defect.
    fn add_forced_changes(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        ue_current_upper: &[f64],
        ue_current_lower: &[f64],
        ue_from_mass_upper: &[f64],
        ue_from_mass_lower: &[f64],
    ) {
        // Compute DUE for upper surface
        let due_upper: Vec<f64> = (0..upper_stations.len())
            .map(|i| {
                let curr = ue_current_upper.get(i).copied().unwrap_or(0.0);
                let from_mass = ue_from_mass_upper.get(i).copied().unwrap_or(0.0);
                curr - from_mass
            })
            .collect();

        // Compute DUE for lower surface
        let due_lower: Vec<f64> = (0..lower_stations.len())
            .map(|i| {
                let curr = ue_current_lower.get(i).copied().unwrap_or(0.0);
                let from_mass = ue_from_mass_lower.get(i).copied().unwrap_or(0.0);
                curr - from_mass
            })
            .collect();

        // Add forced changes to upper surface residuals
        // Start from IBL=1 (which is IV=1 in global indexing, XFOIL's IBL=2)
        for ibl in 1..upper_stations.len() {
            let iv = self.to_global(0, ibl);
            if iv > self.nsys {
                continue;
            }

            let s1 = &upper_stations[ibl - 1];
            let s2 = &upper_stations[ibl];

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

            let dds2 = d2_u2 * due_upper.get(ibl).copied().unwrap_or(0.0);
            let dds1 = d1_u1 * due_upper.get(ibl - 1).copied().unwrap_or(0.0);
            let due2 = due_upper.get(ibl).copied().unwrap_or(0.0);
            let due1 = due_upper.get(ibl - 1).copied().unwrap_or(0.0);

            for k in 0..3 {
                self.vdel[iv][k] += self.vs1_delta[iv][k] * dds1
                    + self.vs2_delta[iv][k] * dds2
                    + self.vs1_ue[iv][k] * due1
                    + self.vs2_ue[iv][k] * due2;
            }
        }

        // Add forced changes to lower surface residuals
        // Start from IBL=1 (which is IV=1+n_upper in global indexing, XFOIL's IBL=2)
        for ibl in 1..lower_stations.len() {
            let iv = self.to_global(1, ibl);
            if iv > self.nsys {
                continue;
            }

            let s1 = &lower_stations[ibl - 1];
            let s2 = &lower_stations[ibl];

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

            let dds2 = d2_u2 * due_lower.get(ibl).copied().unwrap_or(0.0);
            let dds1 = d1_u1 * due_lower.get(ibl - 1).copied().unwrap_or(0.0);
            let due2 = due_lower.get(ibl).copied().unwrap_or(0.0);
            let due1 = due_lower.get(ibl - 1).copied().unwrap_or(0.0);

            for k in 0..3 {
                self.vdel[iv][k] += self.vs1_delta[iv][k] * dds1
                    + self.vs2_delta[iv][k] * dds2
                    + self.vs1_ue[iv][k] * due1
                    + self.vs2_ue[iv][k] * due2;
            }
        }
    }

    /// Get maximum absolute residual
    pub fn max_residual(&self) -> f64 {
        self.vdel
            .iter()
            .skip(1)
            .take(self.nsys)
            .flat_map(|r| r.iter())
            .map(|&v| v.abs())
            .fold(0.0, f64::max)
    }

    /// Get RMS residual
    pub fn rms_residual(&self) -> f64 {
        let sum_sq: f64 = self
            .vdel
            .iter()
            .skip(1)
            .take(self.nsys)
            .flat_map(|r| r.iter())
            .map(|&v| v * v)
            .sum();
        (sum_sq / (3.0 * self.nsys as f64)).sqrt()
    }

    /// Check convergence
    pub fn is_converged(&self, tolerance: f64) -> bool {
        self.max_residual() < tolerance
    }
}

/// Solve the global Newton system using XFOIL's BLSOLV algorithm
///
/// This implements the full XFOIL BLSOLV algorithm with:
/// - Forward sweep: eliminate VB blocks and lower VM columns
/// - VZ block handling at trailing edge
/// - Back substitution with VM columns
///
/// # Arguments
/// * `system` - The global Newton system to solve (modified in place)
///
/// # Returns
/// Solution vector indexed by global IV, containing [dAmpl/dCtau, dTheta, dMass]
///
/// # XFOIL Reference
/// XFOIL xsolve.f BLSOLV (lines 283-485)
pub fn solve_global_system(system: &mut GlobalNewtonSystem) -> Vec<[f64; 3]> {
    let nsys = system.nsys;
    if nsys < 2 {
        return vec![[0.0; 3]; nsys + 1];
    }


    // VACCEL threshold for sparse VM elimination
    let vaccel = 0.005;

    // Working copies
    let mut vdel = system.vdel.clone();
    let mut va_mod = system.va.clone();
    let mut vm_mod = system.vm.clone();

    // Global index of upper TE (where VZ block applies)
    let ivte = system.to_global(0, system.iblte_upper);
    // Global index of lower wake start
    let ivz_target = system.to_global(1, system.iblte_lower + 1);

    // === Forward Sweep ===
    for iv in 1..=nsys {
        let ivp = iv + 1;

        // === Invert VA block (rows 1-2) ===

        // Normalize first row by VA[0][0]
        let pivot = va_mod[iv][0][0];
        if pivot.abs() < 1e-30 {
            continue;
        }
        let pivot_inv = 1.0 / pivot;
        va_mod[iv][0][1] *= pivot_inv;
        for l in iv..=nsys {
            vm_mod[iv][l][0] *= pivot_inv;
        }
        vdel[iv][0] *= pivot_inv;

        // Eliminate lower first column in VA block (rows 2, 3)
        for k in 1..3 {
            let vtmp = va_mod[iv][k][0];
            va_mod[iv][k][1] -= vtmp * va_mod[iv][0][1];
            for l in iv..=nsys {
                vm_mod[iv][l][k] -= vtmp * vm_mod[iv][l][0];
            }
            vdel[iv][k] -= vtmp * vdel[iv][0];
        }

        // Normalize second row by VA[1][1]
        let pivot = va_mod[iv][1][1];
        if pivot.abs() < 1e-30 {
            continue;
        }
        let pivot_inv = 1.0 / pivot;
        for l in iv..=nsys {
            vm_mod[iv][l][1] *= pivot_inv;
        }
        vdel[iv][1] *= pivot_inv;

        // Eliminate lower second column (row 3 only)
        let vtmp = va_mod[iv][2][1];
        for l in iv..=nsys {
            vm_mod[iv][l][2] -= vtmp * vm_mod[iv][l][1];
        }
        vdel[iv][2] -= vtmp * vdel[iv][1];

        // Normalize third row by VM(3,IV,IV) - mass coupling diagonal
        let pivot = vm_mod[iv][iv][2];
        if pivot.abs() < 1e-30 {
            continue;
        }
        let pivot_inv = 1.0 / pivot;
        for l in ivp..=nsys {
            vm_mod[iv][l][2] *= pivot_inv;
        }
        vdel[iv][2] *= pivot_inv;

        // Back-substitute rows 1-2 from row 3
        let vtmp1 = vm_mod[iv][iv][0];
        let vtmp2 = vm_mod[iv][iv][1];
        for l in ivp..=nsys {
            vm_mod[iv][l][0] -= vtmp1 * vm_mod[iv][l][2];
            vm_mod[iv][l][1] -= vtmp2 * vm_mod[iv][l][2];
        }
        vdel[iv][0] -= vtmp1 * vdel[iv][2];
        vdel[iv][1] -= vtmp2 * vdel[iv][2];

        // Back-substitute row 1 from row 2
        let vtmp = va_mod[iv][0][1];
        for l in ivp..=nsys {
            vm_mod[iv][l][0] -= vtmp * vm_mod[iv][l][1];
        }
        vdel[iv][0] -= vtmp * vdel[iv][1];

        if iv >= nsys {
            continue;
        }

        // === Handle VZ block at upper TE ===
        if iv == ivte && ivz_target <= nsys {
            // Add VZ contribution to lower wake equations
            for k in 0..3 {
                let vtmp1 = system.vz[k][0];
                let vtmp2 = system.vz[k][1];
                for l in ivp..=nsys {
                    vm_mod[ivz_target][l][k] -= vtmp1 * vm_mod[iv][l][0] + vtmp2 * vm_mod[iv][l][1];
                }
                vdel[ivz_target][k] -= vtmp1 * vdel[iv][0] + vtmp2 * vdel[iv][1];
            }
        }

        // === Eliminate VB[IV+1] block ===
        for k in 0..3 {
            let vtmp1 = system.vb[ivp][k][0];
            let vtmp2 = system.vb[ivp][k][1];
            let vtmp3 = vm_mod[ivp][iv][k];
            for l in ivp..=nsys {
                vm_mod[ivp][l][k] -= vtmp1 * vm_mod[iv][l][0]
                    + vtmp2 * vm_mod[iv][l][1]
                    + vtmp3 * vm_mod[iv][l][2];
            }
            vdel[ivp][k] -= vtmp1 * vdel[iv][0] + vtmp2 * vdel[iv][1] + vtmp3 * vdel[iv][2];
        }

        if ivp >= nsys {
            continue;
        }

        // === Eliminate lower VM column (sparse, with VACCEL threshold) ===
        for kv in (iv + 2)..=nsys {
            let vtmp1 = vm_mod[kv][iv][0];
            let vtmp2 = vm_mod[kv][iv][1];
            let vtmp3 = vm_mod[kv][iv][2];

            if vtmp1.abs() > vaccel {
                for l in ivp..=nsys {
                    vm_mod[kv][l][0] -= vtmp1 * vm_mod[iv][l][2];
                }
                vdel[kv][0] -= vtmp1 * vdel[iv][2];
            }
            if vtmp2.abs() > vaccel {
                for l in ivp..=nsys {
                    vm_mod[kv][l][1] -= vtmp2 * vm_mod[iv][l][2];
                }
                vdel[kv][1] -= vtmp2 * vdel[iv][2];
            }
            if vtmp3.abs() > vaccel {
                for l in ivp..=nsys {
                    vm_mod[kv][l][2] -= vtmp3 * vm_mod[iv][l][2];
                }
                vdel[kv][2] -= vtmp3 * vdel[iv][2];
            }
        }
    }

    // === Back Substitution ===
    // Eliminate upper VM columns using row 3 (mass) solution
    for iv in (2..=nsys).rev() {
        let vtmp = vdel[iv][2];

        for kv in (1..iv).rev() {
            vdel[kv][0] -= vm_mod[kv][iv][0] * vtmp;
            vdel[kv][1] -= vm_mod[kv][iv][1] * vtmp;
            vdel[kv][2] -= vm_mod[kv][iv][2] * vtmp;
        }
    }

    vdel
}

/// Emit debug event for the BLSOLV solution
///
/// This dumps the solution deltas [dCtau/dAmpl, dTheta, dMass] matching XFOIL's DBGBLSOLVSOLUTION.
/// Should be called after `solve_global_system` returns.
///
/// # Arguments
/// * `iteration` - Viscous iteration number
/// * `deltas` - Solution from `solve_global_system`
pub fn emit_blsolv_solution_debug(iteration: usize, deltas: &[[f64; 3]]) {
    if !rustfoil_bl::is_debug_active() {
        return;
    }

    let nsys = deltas.len().saturating_sub(1); // deltas[0] is unused
    let nout = nsys.min(30);

    // Extract solution deltas (skip index 0)
    let mut solution_deltas: Vec<[f64; 3]> = Vec::with_capacity(nout);
    for iv in 1..=nout {
        if iv < deltas.len() {
            solution_deltas.push(deltas[iv]);
        }
    }

    rustfoil_bl::add_event(rustfoil_bl::DebugEvent::blsolv_solution(
        iteration,
        nsys,
        solution_deltas,
    ));
}

/// Apply global Newton updates to both surfaces
///
/// # Arguments
/// * `upper_stations` - Upper surface stations to update
/// * `lower_stations` - Lower surface stations to update
/// * `deltas` - Solution from `solve_global_system`
/// * `system` - The global Newton system (for index mapping)
/// * `relaxation` - Relaxation factor (0 < rlx <= 1)
/// * `dij` - Optional DIJ influence matrix for VI coupling
/// * `upper_ue_inv` - Optional inviscid edge velocities (upper surface)
/// * `lower_ue_inv` - Optional inviscid edge velocities (lower surface)
/// * `msq` - Mach number squared for blvar
/// * `re` - Reynolds number for blvar
pub fn apply_global_updates(
    upper_stations: &mut [BlStation],
    lower_stations: &mut [BlStation],
    deltas: &[[f64; 3]],
    system: &GlobalNewtonSystem,
    relaxation: f64,
    dij: Option<&nalgebra::DMatrix<f64>>,
    upper_ue_inv: Option<&[f64]>,
    lower_ue_inv: Option<&[f64]>,
    msq: f64,
    re: f64,
) {
    // === XFOIL-style normalized relaxation (xbl.f lines 1527-1597) ===
    // Compute global relaxation factor based on normalized changes
    let dhi = 1.5;   // Max 150% relative increase
    let dlo = -0.5;  // Max 50% relative decrease (CRITICAL for ctau stability)
    
    let mut rlx = relaxation.clamp(0.0, 1.0);
    
    // First pass: compute required relaxation from all stations
    // Upper surface
    for ibl in 1..upper_stations.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = &deltas[iv];
        let station = &upper_stations[ibl];
        
        // Ctau/ampl normalized change (XFOIL lines 1545-1546)
        let dn1 = if station.is_laminar {
            delta[0] / 10.0  // Laminar: fixed scale
        } else if station.ctau.abs() > 1e-10 {
            delta[0] / station.ctau  // Turbulent: normalize by current value
        } else {
            delta[0] / 0.03  // Fallback for near-zero ctau
        };
        
        // Theta normalized change
        let dn2 = if station.theta.abs() > 1e-12 {
            delta[1] / station.theta
        } else {
            0.0
        };
        
        // Apply relaxation limits (XFOIL lines 1563-1576)
        let rdn1 = rlx * dn1;
        if rdn1 > dhi && dn1.abs() > 1e-12 { rlx = dhi / dn1; }
        if rdn1 < dlo && dn1.abs() > 1e-12 { rlx = dlo / dn1; }
        
        let rdn2 = rlx * dn2;
        if rdn2 > dhi && dn2.abs() > 1e-12 { rlx = dhi / dn2; }
        if rdn2 < dlo && dn2.abs() > 1e-12 { rlx = dlo / dn2; }
    }
    
    // Lower surface
    for ibl in 1..lower_stations.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = &deltas[iv];
        let station = &lower_stations[ibl];
        
        // Ctau/ampl normalized change
        let dn1 = if station.is_laminar {
            delta[0] / 10.0
        } else if station.ctau.abs() > 1e-10 {
            delta[0] / station.ctau
        } else {
            delta[0] / 0.03
        };
        
        // Theta normalized change
        let dn2 = if station.theta.abs() > 1e-12 {
            delta[1] / station.theta
        } else {
            0.0
        };
        
        // Apply relaxation limits
        let rdn1 = rlx * dn1;
        if rdn1 > dhi && dn1.abs() > 1e-12 { rlx = dhi / dn1; }
        if rdn1 < dlo && dn1.abs() > 1e-12 { rlx = dlo / dn1; }
        
        let rdn2 = rlx * dn2;
        if rdn2 > dhi && dn2.abs() > 1e-12 { rlx = dhi / dn2; }
        if rdn2 < dlo && dn2.abs() > 1e-12 { rlx = dlo / dn2; }
    }
    
    // Ensure relaxation stays positive
    rlx = rlx.max(0.01);

    // === Step 1: Compute new edge velocities via DIJ coupling ===
    // XFOIL's UPDATE: UNEW(IBL,IS) = UINV(IBL,IS) + sum_j(DIJ * proposed_mass)
    let (upper_new_ue, lower_new_ue) = if let (Some(dij), Some(ue_inv_u), Some(ue_inv_l)) = 
        (dij, upper_ue_inv, lower_ue_inv) 
    {
        compute_new_ue_via_dij(
            upper_stations, lower_stations, deltas, system, dij, ue_inv_u, ue_inv_l
        )
    } else {
        // No DIJ coupling - keep current Ue values
        (
            upper_stations.iter().map(|s| s.u).collect::<Vec<_>>(),
            lower_stations.iter().map(|s| s.u).collect::<Vec<_>>()
        )
    };

    // === Step 2: Update upper surface BL variables ===
    for ibl in 1..upper_stations.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }

        let delta = &deltas[iv];
        let station = &mut upper_stations[ibl];

        // XFOIL UPDATE (xbl.f:1622-1625): x = x + RLX * delta
        // The normalized relaxation above ensures relative changes are bounded
        if station.is_laminar {
            station.ampl += rlx * delta[0];
            station.ampl = station.ampl.max(0.0);
        } else {
            // Store old ctau for additional safeguard
            let ctau_old = station.ctau;
            station.ctau += rlx * delta[0];
            
            // XFOIL post-update safeguard (xbl.f line 1639-1640):
            // Only upper bound clamp since normalized relaxation prevents large decreases
            // But add a floor based on old value as extra protection near transition
            let ctau_min = (ctau_old * 0.5).max(1e-6);  // Never drop below 50% or 1e-6
            station.ctau = station.ctau.clamp(ctau_min, 0.25);
        }

        // Store old theta for safeguard check
        let theta_old = station.theta;
        station.theta += rlx * delta[1];
        
        // Safeguard: don't let theta decrease by more than 50% per iteration
        // (This should already be handled by normalized relaxation, but add as backup)
        let theta_min = theta_old * 0.5;
        station.theta = station.theta.max(theta_min.max(1e-8));

        // delta[2] is mass change; convert to delta_star using XFOIL formula:
        // DDSTR = (DMASS - DSTR*DUEDG) / UEDG
        // where DUEDG = UNEW - UEDG is the *proposed* Ue change
        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u);
        let due = new_u - station.u;
        
        // XFOIL formula: delta_dstar = (delta_mass - dstar * delta_ue) / ue
        let delta_mass = delta[2];
        let d_dstar = (delta_mass - station.delta_star * due) / station.u.max(1e-6);
        
        // Apply Ue update with relaxation (XFOIL: UEDG = UEDG + RLX*DUEDG)
        station.u += rlx * due;
        // Safeguard: Ue must remain positive for attached flow
        station.u = station.u.max(0.01);
        
        // Apply delta_star update (XFOIL: DSTR = DSTR + RLX*DDSTR)
        station.delta_star += rlx * d_dstar;
        // Ensure H is physical: H >= 1.0 (displacement thickness >= momentum thickness)
        station.delta_star = station.delta_star.max(station.theta);

        station.h = (station.delta_star / station.theta).clamp(1.0, 20.0);
        station.mass_defect = station.u * station.delta_star;
        
        // Recompute secondary variables (Hk, Cf, Cd, etc.) for next iteration
        let flow_type = if station.is_laminar { FlowType::Laminar } else { FlowType::Turbulent };
        blvar(station, flow_type, msq, re);
    }

    // === Step 3: Update lower surface BL variables ===
    for ibl in 1..lower_stations.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }

        let delta = &deltas[iv];
        let station = &mut lower_stations[ibl];

        // XFOIL UPDATE (xbl.f:1622-1625): x = x + RLX * delta
        // The normalized relaxation above ensures relative changes are bounded
        if station.is_laminar {
            station.ampl += rlx * delta[0];
            station.ampl = station.ampl.max(0.0);
        } else {
            // Store old ctau for additional safeguard
            let ctau_old = station.ctau;
            station.ctau += rlx * delta[0];
            
            // XFOIL post-update safeguard (xbl.f line 1639-1640):
            // Only upper bound clamp since normalized relaxation prevents large decreases
            // But add a floor based on old value as extra protection near transition
            let ctau_min = (ctau_old * 0.5).max(1e-6);  // Never drop below 50% or 1e-6
            station.ctau = station.ctau.clamp(ctau_min, 0.25);
        }

        // Store old theta for safeguard check
        let theta_old = station.theta;
        station.theta += rlx * delta[1];
        
        // Safeguard: don't let theta decrease by more than 50% per iteration
        // (This should already be handled by normalized relaxation, but add as backup)
        let theta_min = theta_old * 0.5;
        station.theta = station.theta.max(theta_min.max(1e-8));

        // delta[2] is mass change; convert to delta_star using XFOIL formula:
        // DDSTR = (DMASS - DSTR*DUEDG) / UEDG
        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(station.u);
        let due = new_u - station.u;
        
        // XFOIL formula: delta_dstar = (delta_mass - dstar * delta_ue) / ue
        let delta_mass = delta[2];
        let d_dstar = (delta_mass - station.delta_star * due) / station.u.max(1e-6);
        
        // Apply Ue update with relaxation (XFOIL: UEDG = UEDG + RLX*DUEDG)
        station.u += rlx * due;
        // Safeguard: Ue must remain positive for attached flow
        station.u = station.u.max(0.01);

        // Apply delta_star update (XFOIL: DSTR = DSTR + RLX*DDSTR)
        station.delta_star += rlx * d_dstar;
        // Ensure H is physical: H >= 1.0 (displacement thickness >= momentum thickness)
        station.delta_star = station.delta_star.max(station.theta);

        station.h = (station.delta_star / station.theta).clamp(1.0, 20.0);
        station.mass_defect = station.u * station.delta_star;
        
        // Recompute secondary variables (Hk, Cf, Cd, etc.) for next iteration
        let flow_type = if station.is_laminar { FlowType::Laminar } else { FlowType::Turbulent };
        blvar(station, flow_type, msq, re);
    }
}

/// Compute new edge velocities via DIJ influence matrix.
///
/// XFOIL's UPDATE formula:
/// UNEW(i) = UINV(i) + sum_j(-VTI_i * VTI_j * DIJ(i,j) * (MASS(j) + VDEL(j)))
fn compute_new_ue_via_dij(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    deltas: &[[f64; 3]],
    system: &GlobalNewtonSystem,
    dij: &nalgebra::DMatrix<f64>,
    upper_ue_inv: &[f64],
    lower_ue_inv: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n_upper = upper_stations.len();
    let n_lower = lower_stations.len();
    
    let mut upper_new_ue = vec![0.0; n_upper];
    let mut lower_new_ue = vec![0.0; n_lower];
    
    
    // VTI signs: Upper surface = +1, Lower surface = -1
    let vti_upper = 1.0_f64;
    let vti_lower = -1.0_f64;
    
    // Compute new Ue for upper surface stations (IS=1)
    for i in 0..n_upper {
        let panel_i = upper_stations[i].panel_idx;
        if panel_i >= dij.nrows() {
            upper_new_ue[i] = upper_stations[i].u;
            continue;
        }
        
        // Start with inviscid Ue
        let uinv_i = upper_ue_inv.get(i).copied().unwrap_or(upper_stations[i].u);
        let mut dui = 0.0_f64;
        
        // XFOIL's UPDATE formula:
        // UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
        // DUI = sum_j(UE_M * (MASS(j) + VDEL_mass(j)))
        
        // Upper surface contribution (JS=1, VTI_j = +1)
        // UE_M = -VTI_upper * VTI_upper * DIJ = -1 * 1 * DIJ = -DIJ
        for j in 0..n_upper {
            let panel_j = upper_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = upper_stations[j].mass_defect + delta_mass;
            
            // UE_M = -VTI_i * VTI_j * DIJ(i,j)
            let ue_m = -vti_upper * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        
        // Lower surface contribution (JS=2, VTI_j = -1)
        // UE_M = -VTI_upper * VTI_lower * DIJ = -1 * (-1) * DIJ = +DIJ
        for j in 0..n_lower {
            let panel_j = lower_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = lower_stations[j].mass_defect + delta_mass;
            
            let ue_m = -vti_upper * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        
        let ue_new = uinv_i + dui;
        
        // Safeguard: limit change per iteration
        let due_unclamped = ue_new - uinv_i;
        let max_due = 0.05 * uinv_i.abs().max(0.1);
        let due_clamped = due_unclamped.clamp(-max_due, max_due);
        upper_new_ue[i] = (uinv_i + due_clamped).clamp(0.01, 5.0);
    }
    
    // Compute new Ue for lower surface stations (IS=2)
    for i in 0..n_lower {
        let panel_i = lower_stations[i].panel_idx;
        if panel_i >= dij.nrows() {
            lower_new_ue[i] = lower_stations[i].u;
            continue;
        }
        
        let uinv_i = lower_ue_inv.get(i).copied().unwrap_or(lower_stations[i].u);
        let mut dui = 0.0_f64;
        
        // Upper surface contribution (JS=1, VTI_j = +1)
        // UE_M = -VTI_lower * VTI_upper * DIJ = -(-1) * 1 * DIJ = +DIJ
        for j in 0..n_upper {
            let panel_j = upper_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = upper_stations[j].mass_defect + delta_mass;
            
            let ue_m = -vti_lower * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        
        // Lower surface contribution (JS=2, VTI_j = -1)
        // UE_M = -VTI_lower * VTI_lower * DIJ = -(-1)*(-1) * DIJ = -DIJ
        for j in 0..n_lower {
            let panel_j = lower_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = lower_stations[j].mass_defect + delta_mass;
            
            let ue_m = -vti_lower * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        
        let ue_new = uinv_i + dui;
        
        // Safeguard: limit change per iteration
        let due_unclamped = ue_new - uinv_i;
        let max_due = 0.05 * uinv_i.abs().max(0.1);
        let due_clamped = due_unclamped.clamp(-max_due, max_due);
        lower_new_ue[i] = (uinv_i + due_clamped).clamp(0.01, 5.0);
    }
    
    (upper_new_ue, lower_new_ue)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_global_system_creation() {
        let system = GlobalNewtonSystem::new(80, 80, 75, 75);

        // NSYS = 80 + 80 - 2 = 158
        assert_eq!(system.nsys, 158);
        assert_eq!(system.n_upper, 80);
        assert_eq!(system.n_lower, 80);
        assert_eq!(system.iblte_upper, 75);
        assert_eq!(system.iblte_lower, 75);
    }

    #[test]
    fn test_index_mapping() {
        let system = GlobalNewtonSystem::new(80, 80, 75, 75);

        // Upper surface: direct mapping
        assert_eq!(system.to_global(0, 1), 1);
        assert_eq!(system.to_global(0, 50), 50);
        assert_eq!(system.to_global(0, 79), 79);

        // Lower surface: offset by n_upper - 1
        assert_eq!(system.to_global(1, 1), 80);
        assert_eq!(system.to_global(1, 50), 129);

        // Reverse mapping
        assert_eq!(system.from_global(1), (0, 1));
        assert_eq!(system.from_global(50), (0, 50));
        assert_eq!(system.from_global(80), (1, 1));
        assert_eq!(system.from_global(129), (1, 50));
    }

    #[test]
    fn test_vti_signs() {
        use rustfoil_bl::equations::blvar;

        let n = 10;
        let mut upper_stations: Vec<BlStation> = (0..n)
            .map(|i| {
                let mut s = BlStation::new();
                s.x = 0.1 * i as f64;
                s.panel_idx = i;
                blvar(&mut s, FlowType::Laminar, 0.0, 1e6);
                s
            })
            .collect();
        let mut lower_stations = upper_stations.clone();
        for (i, s) in lower_stations.iter_mut().enumerate() {
            s.panel_idx = n + i;
        }

        let mut system = GlobalNewtonSystem::new(n, n, n - 2, n - 2);
        system.setup_index_mappings(&upper_stations, &lower_stations);

        // Upper surface should have VTI = +1
        for ibl in 1..n {
            let iv = system.to_global(0, ibl);
            if iv <= system.nsys {
                assert_eq!(
                    system.vti[iv], 1.0,
                    "Upper surface VTI should be +1 at IV={}",
                    iv
                );
            }
        }

        // Lower surface should have VTI = -1
        for ibl in 1..n {
            let iv = system.to_global(1, ibl);
            if iv <= system.nsys {
                assert_eq!(
                    system.vti[iv], -1.0,
                    "Lower surface VTI should be -1 at IV={}",
                    iv
                );
            }
        }
    }
}
