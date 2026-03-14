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
use rustfoil_bl::closures::hkin::hkin;
use rustfoil_bl::closures::trchek2_full;
use rustfoil_bl::equations::{bldif_full_simi, bldif_ncrit, blvar, trdif_full, FlowType};
use rustfoil_bl::state::BlStation;

use crate::newton_state::{CanonicalNewtonRow, CanonicalNewtonStateView};

#[derive(Debug, Clone, Copy)]
struct TeWakeTerms {
    tte: f64,
    dte: f64,
    cte: f64,
    dte_mte1: f64,
    dte_ute1: f64,
    dte_mte2: f64,
    dte_ute2: f64,
    cte_cte1: f64,
    cte_cte2: f64,
    cte_tte1: f64,
    cte_tte2: f64,
}

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
    /// Dimensions: 3 equations (rows) x 2 coupling terms (columns)
    pub vz: [[f64; 2]; 3],
    /// Mass influence matrix VM[IV][JV][k] - sensitivity of eq k at IV to mass at JV
    /// Dimensions: nsys x nsys x 3
    pub vm: Vec<Vec<[f64; 3]>>,
    /// Primary RHS residual / solution vector.
    /// VDEL[IV] = [res_third, res_mom, res_shape]
    pub vdel: Vec<[f64; 3]>,
    /// Operating-variable sensitivity RHS, analogous to XFOIL's second VDEL column.
    pub vdel_operating: Vec<[f64; 3]>,

    // === Index mappings ===
    /// Surface direction signs: +1.0 for upper, -1.0 for lower
    pub vti: Vec<f64>,
    /// Panel index for each global station (maps to DIJ matrix)
    pub panel_idx: Vec<usize>,
    /// Flag for stations at or after TE (where wake equations apply)
    pub is_wake: Vec<bool>,
    /// Approximate airfoil arc length used for XFOIL-style VACCEL scaling
    pub airfoil_arc_length: f64,

    /// VS1 columns for delta_star derivatives (stored for forced changes)
    pub vs1_delta: Vec<[f64; 3]>,
    /// VS1 columns for Ue derivatives (stored for forced changes)
    pub vs1_ue: Vec<[f64; 3]>,
    /// VS2 columns for delta_star derivatives
    pub vs2_delta: Vec<[f64; 3]>,
    /// VS2 columns for Ue derivatives
    pub vs2_ue: Vec<[f64; 3]>,
    
    // === Stagnation point coupling terms (XFOIL XI_ULE) ===
    /// VS1 column 5 (X derivatives) for stagnation coupling
    pub vs1_x: Vec<[f64; 3]>,
    /// VS2 column 5 (X derivatives) for stagnation coupling
    pub vs2_x: Vec<[f64; 3]>,
    /// SST_GO: dSST/dGAM(IST) - stagnation arc length derivative
    pub sst_go: f64,
    /// SST_GP: dSST/dGAM(IST+1) - stagnation arc length derivative
    pub sst_gp: f64,
    /// DULE1: Forced change in leading edge Ue (upper surface)
    pub dule1: f64,
    /// DULE2: Forced change in leading edge Ue (lower surface)
    pub dule2: f64,
    /// Current Newton iteration number (for debug output)
    pub current_iteration: usize,
    /// ANTE: base thickness contribution at blunt trailing edge (XFOIL WGAP(1))
    pub ante: f64,
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
            vz: [[0.0; 2]; 3],
            vm: vec![vec![[0.0; 3]; nsys + 1]; nsys + 1],
            vdel: vec![[0.0; 3]; nsys + 1],
            vdel_operating: vec![[0.0; 3]; nsys + 1],
            vti: vec![1.0; nsys + 1],
            panel_idx: vec![0; nsys + 1],
            is_wake: vec![false; nsys + 1],
            airfoil_arc_length: 1.0,
            vs1_delta: vec![[0.0; 3]; nsys + 1],
            vs1_ue: vec![[0.0; 3]; nsys + 1],
            vs2_delta: vec![[0.0; 3]; nsys + 1],
            vs2_ue: vec![[0.0; 3]; nsys + 1],
            // Stagnation coupling - initialized to zero, computed in build
            vs1_x: vec![[0.0; 3]; nsys + 1],
            vs2_x: vec![[0.0; 3]; nsys + 1],
            sst_go: 0.0,
            sst_gp: 0.0,
            dule1: 0.0,
            dule2: 0.0,
            // Current iteration - set in build_global_system
            current_iteration: 0,
            ante: 0.0,
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

    /// Set the stagnation point derivatives for XI_ULE coupling.
    ///
    /// These values come from `find_stagnation_with_derivs` and are used
    /// to couple the stagnation point location to the Newton system.
    ///
    /// # Arguments
    /// * `sst_go` - dSST/dGAM(IST) - sensitivity to upstream gamma
    /// * `sst_gp` - dSST/dGAM(IST+1) - sensitivity to downstream gamma
    ///
    /// # XFOIL Reference
    /// xpanel.f STFIND (lines 1438-1439)
    pub fn set_stagnation_derivs(&mut self, sst_go: f64, sst_gp: f64) {
        self.sst_go = if sst_go.is_finite() { sst_go } else { 0.0 };
        self.sst_gp = if sst_gp.is_finite() { sst_gp } else { 0.0 };
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
    /// * `dij` - Global DIJ matrix for mass defect coupling
    /// * `ue_inviscid_upper` - Inviscid edge velocities (upper)
    /// * `ue_inviscid_lower` - Inviscid edge velocities (lower)
    /// * `ue_from_mass_upper` - Canonical UESET result for upper stations
    /// * `ue_from_mass_lower` - Canonical UESET result for lower stations
    /// * `msq` - Mach number squared
    /// * `re` - Reynolds number
    /// * `iteration` - Newton iteration number (for debug output)
    ///
    /// # XFOIL Reference
    /// XFOIL xbl.f SETBL (lines 21-500)
    #[allow(clippy::too_many_arguments)]
    pub fn build_global_system(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        ue_current_upper: &[f64],
        ue_current_lower: &[f64],
        upper_flow_types: &[FlowType],
        lower_flow_types: &[FlowType],
        dij: &DMatrix<f64>,
        ue_inviscid_upper: &[f64],
        ue_inviscid_lower: &[f64],
        ue_from_mass_upper: &[f64],
        ue_from_mass_lower: &[f64],
        ue_operating_upper: &[f64],
        ue_operating_lower: &[f64],
        ncrit: f64,
        msq: f64,
        re: f64,
        iteration: usize,
    ) {
        // Store iteration for debug events
        self.current_iteration = iteration;
        // Reset arrays
        for iv in 0..=self.nsys {
            self.va[iv] = [[0.0; 3]; 3];
            self.vb[iv] = [[0.0; 3]; 3];
            self.vdel[iv] = [0.0; 3];
            self.vdel_operating[iv] = [0.0; 3];
            self.vs1_delta[iv] = [0.0; 3];
            self.vs1_ue[iv] = [0.0; 3];
            self.vs2_delta[iv] = [0.0; 3];
            self.vs2_ue[iv] = [0.0; 3];
            self.vs1_x[iv] = [0.0; 3];
            self.vs2_x[iv] = [0.0; 3];
            for jv in 0..=self.nsys {
                self.vm[iv][jv] = [0.0; 3];
            }
        }
        self.vz = [[0.0; 2]; 3];
        // Note: sst_go, sst_gp, dule1, dule2 are computed in add_forced_changes

        // === Setup index mappings and VTI signs ===
        self.setup_index_mappings(upper_stations, lower_stations);
        self.airfoil_arc_length = upper_stations
            .get(self.iblte_upper)
            .map(|s| s.x)
            .unwrap_or(0.0)
            + lower_stations
                .get(self.iblte_lower)
                .map(|s| s.x)
                .unwrap_or(0.0);

        // === Step 1: Save current Ue values (USAV) for forced changes ===
        // Preserve the signed row UEDG values here so the stagnation coupling
        // path can see the same near-zero sign flips as XFOIL.
        let ue_current_upper: Vec<f64> = if ue_current_upper.len() == upper_stations.len() {
            ue_current_upper.to_vec()
        } else {
            upper_stations.iter().map(|s| s.u).collect()
        };
        let ue_current_lower: Vec<f64> = if ue_current_lower.len() == lower_stations.len() {
            ue_current_lower.to_vec()
        } else {
            lower_stations.iter().map(|s| s.u).collect()
        };

        // === Step 2: Compute Ue from mass defect (UESET) ===
        // IMPORTANT: Sum over BOTH surfaces to get proper mass coupling
        let ue_from_mass_upper = if ue_from_mass_upper.len() == upper_stations.len() {
            ue_from_mass_upper.to_vec()
        } else {
            self.compute_ue_from_mass_both(
                upper_stations,
                lower_stations,
                ue_inviscid_upper,
                dij,
                0,
            )
        };
        let ue_from_mass_lower = if ue_from_mass_lower.len() == lower_stations.len() {
            ue_from_mass_lower.to_vec()
        } else {
            self.compute_ue_from_mass_both(
                upper_stations,
                lower_stations,
                ue_inviscid_lower,
                dij,
                1,
            )
        };

        // === Step 3: Build equations for upper surface ===
        self.build_surface_equations(
            upper_stations,
            upper_flow_types,
            &ue_current_upper,
            &ue_from_mass_upper,
            ncrit,
            msq,
            re,
            0, // surface = upper
        );

        // === Step 4: Build equations for lower surface ===
        self.build_surface_equations(
            lower_stations,
            lower_flow_types,
            &ue_current_lower,
            &ue_from_mass_lower,
            ncrit,
            msq,
            re,
            1, // surface = lower
        );

        // === Step 4.5: Replace first wake interval with TESYS coupling ===
        self.apply_te_wake_interface(upper_stations, lower_stations);

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
            ue_inviscid_upper,
            ue_inviscid_lower,
        );
        self.build_operating_rhs(
            upper_stations,
            lower_stations,
            ue_operating_upper,
            ue_operating_lower,
        );
        
        // Debug: print final VDEL after forced changes
        if std::env::var("RUSTFOIL_NEWTON_DEBUG").is_ok() {
            eprintln!("\n[NEWTON DEBUG] After forced changes (FINAL VDEL):");
            for ibl in 1..=5.min(upper_stations.len().saturating_sub(1)) {
                let iv = self.to_global(0, ibl);
                if iv <= self.nsys {
                    eprintln!("  upper IBL={} IV={}: VDEL=[{:>12.6e}, {:>12.6e}, {:>12.6e}]",
                        ibl, iv, self.vdel[iv][0], self.vdel[iv][1], self.vdel[iv][2]);
                }
            }
            for ibl in 1..=3.min(lower_stations.len().saturating_sub(1)) {
                let iv = self.to_global(1, ibl);
                if iv <= self.nsys {
                    eprintln!("  lower IBL={} IV={}: VDEL=[{:>12.6e}, {:>12.6e}, {:>12.6e}]",
                        ibl, iv, self.vdel[iv][0], self.vdel[iv][1], self.vdel[iv][2]);
                }
            }
        }
        
        // === Step 8: Zero SIMI station residuals ===
        // IMPORTANT: XFOIL does NOT zero SIMI residuals. The large residuals at SIMI
        // stations (-42.30, 31.47 for NACA 0012 α=4°) are the forced change terms
        // that drive convergence. We previously zeroed these, which may have been
        // causing divergence issues.
        //
        // DISABLED: Let SIMI residuals contribute to the Newton system
        // let iv_upper_simi = 1;
        // let iv_lower_simi = self.n_upper;
        // 
        // if iv_upper_simi <= self.nsys {
        //     self.vdel[iv_upper_simi][1] = 0.0; // Momentum
        //     self.vdel[iv_upper_simi][2] = 0.0; // Shape
        // }
        // if iv_lower_simi <= self.nsys {
        //     self.vdel[iv_lower_simi][1] = 0.0; // Momentum
        //     self.vdel[iv_lower_simi][2] = 0.0; // Shape
        // }
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

        // Emit the full assembled system so late-iteration lower-side rows can be
        // compared directly against XFOIL's `FULL_ITER` dump.
        let nout = self.nsys;

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

        if std::env::var("RUSTFOIL_SETBL_VDEL_DEBUG").is_ok() {
            for iv in [27usize, 28, 29, 30, 159, 160, 161, 162, 163] {
                if iv <= self.nsys {
                    let vals = self.vdel[iv];
                    eprintln!(
                        "[RUST SETBL VDEL] iv={} vals=[{:.12e}, {:.12e}, {:.12e}]",
                        iv, vals[0], vals[1], vals[2]
                    );
                }
            }
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

        let focused_vm_row = |row: usize| -> Option<Vec<[f64; 3]>> {
            if row > self.nsys {
                None
            } else {
                let mut out = Vec::with_capacity(nout);
                for jv in 1..=nout {
                    out.push(self.vm[row][jv]);
                }
                Some(out)
            }
        };
        let vm_row24 = focused_vm_row(24);
        let vm_row25 = focused_vm_row(25);
        let vm_row26 = focused_vm_row(26);
        let vm_full = if std::env::var("RUSTFOIL_SETBL_FULL_VM").is_ok() {
            Some(
                (1..=nout)
                    .map(|iv| (1..=nout).map(|jv| self.vm[iv][jv]).collect())
                    .collect(),
            )
        } else {
            None
        };

        if std::env::var("RUSTFOIL_VM_FOCUS_DEBUG").is_ok() {
            for focus_iv in [24usize, 25, 26] {
                if focus_iv <= self.nsys {
                    let end = nout.min(self.nsys);
                    eprintln!("[VM FOCUS] row={focus_iv}");
                    for jv in 1..=end {
                        let entry = self.vm[focus_iv][jv];
                        if entry[0].abs() > 1.0e-6 || entry[1].abs() > 1.0e-6 || entry[2].abs() > 1.0e-6 {
                            eprintln!(
                                "[VM FOCUS ENTRY] row={} col={} vals=[{:.8e}, {:.8e}, {:.8e}]",
                                focus_iv, jv, entry[0], entry[1], entry[2]
                            );
                        }
                    }
                }
            }
        }

        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::setbl_system(
            iteration,
            self.nsys,
            va_blocks,
            vb_blocks,
            vdel_rhs,
            vm_diag,
            vm_row1,
            vm_row24,
            vm_row25,
            vm_row26,
            vm_full,
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
            let mut dui_upper = 0.0;
            let mut dui_lower = 0.0;

            // Sum over upper surface (VTI = +1)
            // XFOIL starts from JBL=2 (first BL station, skipping stagnation)
            for j in 1..upper_stations.len() {
                let panel_j = upper_stations[j].panel_idx;
                let vti_j = 1.0;  // Upper surface

                if panel_i < dij.nrows() && panel_j < dij.ncols() {
                    let ue_m = -vti_i * vti_j * dij[(panel_i, panel_j)];
                    dui += ue_m * upper_stations[j].mass_defect;
                    dui_upper += ue_m * upper_stations[j].mass_defect;
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
                    dui_lower += ue_m * lower_stations[j].mass_defect;
                }
            }

            // Debug: print DUI breakdown at first few stations
            if std::env::var("RUSTFOIL_DUI_DEBUG").is_ok() && i <= 3 {
                let surface_name = if surface == 0 { "upper" } else { "lower" };
                let ue_inv_i = ue_inviscid.get(i).copied().unwrap_or(0.0);
                
                // Also print DIJ diagonal for this panel
                let dij_diag = if panel_i < dij.nrows() { dij[(panel_i, panel_i)] } else { 0.0 };
                
                eprintln!("[DUI DEBUG] {}[{}] panel={}: ue_inv={:.6}, dui_upper={:.6e}, dui_lower={:.6e}, dui_total={:.6e}, dij_diag={:.6e}",
                    surface_name, i, panel_i, ue_inv_i, dui_upper, dui_lower, dui, dij_diag);
                    
                // For station 1, trace cumulative DUI sum
                if i == 1 && surface == 0 {
                    eprintln!("[DUI TRACE] Breaking down upper[1] DIJ sum:");
                    let mut cum_upper = 0.0;
                    for jj in 1..upper_stations.len().min(20) {
                        let panel_jj = upper_stations[jj].panel_idx;
                        let vti_jj = 1.0;
                        if panel_i < dij.nrows() && panel_jj < dij.ncols() {
                            let dij_val = dij[(panel_i, panel_jj)];
                            let ue_m_jj = -vti_i * vti_jj * dij_val;
                            let mass_jj = upper_stations[jj].mass_defect;
                            let contrib = ue_m_jj * mass_jj;
                            cum_upper += contrib;
                            if jj <= 10 {
                                eprintln!("  j={:2} panel={:3}: DIJ={:>10.4e}, mass={:>10.4e}, contrib={:>10.4e}, cum={:>10.4e}",
                                    jj, panel_jj, dij_val, mass_jj, contrib, cum_upper);
                            }
                        }
                    }
                    eprintln!("  ... cumulative after j=20: {:.4e}", cum_upper);
                    eprintln!("  TOTAL dui_upper: {:.4e}", dui_upper);
                    
                    // Also trace lower surface contribution
                    eprintln!("[DUI TRACE] Breaking down lower contribution to upper[1]:");
                    let mut cum_lower = 0.0;
                    for jj in 1..lower_stations.len().min(10) {
                        let panel_jj = lower_stations[jj].panel_idx;
                        let vti_jj = -1.0;  // Lower surface VTI
                        if panel_i < dij.nrows() && panel_jj < dij.ncols() {
                            let dij_val = dij[(panel_i, panel_jj)];
                            let ue_m_jj = -vti_i * vti_jj * dij_val;
                            let mass_jj = lower_stations[jj].mass_defect;
                            let contrib = ue_m_jj * mass_jj;
                            cum_lower += contrib;
                            eprintln!("  j={:2} panel={:3}: DIJ={:>10.4e}, mass={:>10.4e}, contrib={:>10.4e}, cum={:>10.4e}",
                                jj, panel_jj, dij_val, mass_jj, contrib, cum_lower);
                        }
                    }
                    eprintln!("  TOTAL dui_lower: {:.4e}", dui_lower);
                }
            }

            if std::env::var("RUSTFOIL_UESET_FIRST_DEBUG").is_ok() && surface == 0 && i == 1 {
                let ue_inv_i = ue_inviscid.get(i).copied().unwrap_or(0.0);
                eprintln!(
                    "[UESET FIRST] panel_i={} ue_inv={:.8e} dui_upper={:.8e} dui_lower={:.8e} ue_after={:.8e} mass={:.8e}",
                    panel_i,
                    ue_inv_i,
                    dui_upper,
                    dui_lower,
                    ue_inv_i + dui,
                    upper_stations.get(i).map(|s| s.mass_defect).unwrap_or(0.0)
                );
                for jj in 1..6.min(upper_stations.len()) {
                    let panel_jj = upper_stations[jj].panel_idx;
                    if panel_i < dij.nrows() && panel_jj < dij.ncols() {
                        let dij_val = dij[(panel_i, panel_jj)];
                        let contrib = -dij_val * upper_stations[jj].mass_defect;
                        eprintln!(
                            "[UESET FIRST] upper j={} panel={} dij={:.8e} mass={:.8e} contrib={:.8e}",
                            jj, panel_jj, dij_val, upper_stations[jj].mass_defect, contrib
                        );
                    }
                }
                for jj in 1..6.min(lower_stations.len()) {
                    let panel_jj = lower_stations[jj].panel_idx;
                    if panel_i < dij.nrows() && panel_jj < dij.ncols() {
                        let dij_val = dij[(panel_i, panel_jj)];
                        let contrib = dij_val * lower_stations[jj].mass_defect;
                        eprintln!(
                            "[UESET FIRST] lower j={} panel={} dij={:.8e} mass={:.8e} contrib={:.8e}",
                            jj, panel_jj, dij_val, lower_stations[jj].mass_defect, contrib
                        );
                    }
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
        _ue_current: &[f64],
        _ue_from_mass: &[f64],
        ncrit: f64,
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
            // Note: transitions parameter is no longer used - XFOIL doesn't use TRDIF in global Newton
            // Compute residuals and Jacobian.
            // XFOIL BLSYS uses:
            // - SIMI-special BLDIF at the first station
            // - TRDIF on the transition interval
            // - regular BLDIF everywhere else
            let s1 = &stations[i - 1];
            let is_transition_interval =
                !s1.is_turbulent && !s1.is_wake && s2.is_turbulent && !s2.is_wake;
            let x_forced = if surface == 0 {
                stations.get(self.iblte_upper).map(|s| s.x)
            } else {
                stations.get(self.iblte_lower).map(|s| s.x)
            };
            let live_transition = if is_transition_interval {
                Some(trchek2_full(
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
                    x_forced,
                    msq,
                    re,
                ))
            } else {
                None
            };
            if std::env::var("RUSTFOIL_TRANSITION_BRANCH_DEBUG").is_ok()
                && ((surface == 0 && (30..=31).contains(&i))
                    || (surface == 1 && (61..=62).contains(&i)))
            {
                eprintln!(
                    "[TRANSITION BRANCH] {}[{i}] flow={flow_type:?} has_transition={} transition_flag={} live_transition={} xt={:.8e} ampl2={:.8e}",
                    if surface == 0 { "upper" } else { "lower" },
                    is_transition_interval,
                    live_transition.as_ref().map(|tr| tr.transition).unwrap_or(false),
                    live_transition.as_ref().is_some_and(|tr| tr.transition),
                    live_transition.as_ref().map(|tr| tr.xt).unwrap_or(0.0),
                    live_transition.as_ref().map(|tr| tr.ampl2).unwrap_or(0.0),
                );
                if let Some(tr) = &live_transition {
                    eprintln!(
                        "[TRANSITION DERIVS] {}[{i}] xt_x1={:.12e} xt_x2={:.12e} tt_x2={:.12e} dt_x2={:.12e} ut_x2={:.12e}",
                        if surface == 0 { "upper" } else { "lower" },
                        tr.xt_x1,
                        tr.xt_x2,
                        tr.tt_x2,
                        tr.dt_x2,
                        tr.ut_x2,
                    );
                }
            }
            let (residuals, jacobian) = if i == 1 {
                bldif_full_simi(s2, flow_type, msq, re)
            } else if let Some(tr) = live_transition.as_ref().filter(|tr| tr.transition) {
                trdif_full(s1, s2, tr, ncrit, msq, re)
            } else {
                bldif_ncrit(s1, s2, flow_type, msq, re, ncrit)
            };
            if std::env::var("RUSTFOIL_TRANSITION_BLOCK_DEBUG").is_ok()
                && ((surface == 0 && (i == 30 || i == 31))
                    || (surface == 1 && (i == 61 || i == 62)))
            {
                eprintln!(
                    "[RUST TRANSITION REGION BLOCK] {}[{i}] residuals=[{:.12e}, {:.12e}, {:.12e}] vs1_row1=[{:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}] vs1_row2=[{:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}] vs1_row3=[{:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}] vs2_row1=[{:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}] vs2_row2=[{:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}] vs2_row3=[{:.12e}, {:.12e}, {:.12e}, {:.12e}, {:.12e}]",
                    if surface == 0 { "upper" } else { "lower" },
                    residuals.res_third,
                    residuals.res_mom,
                    residuals.res_shape,
                    jacobian.vs1[0][0], jacobian.vs1[0][1], jacobian.vs1[0][2], jacobian.vs1[0][3], jacobian.vs1[0][4],
                    jacobian.vs1[1][0], jacobian.vs1[1][1], jacobian.vs1[1][2], jacobian.vs1[1][3], jacobian.vs1[1][4],
                    jacobian.vs1[2][0], jacobian.vs1[2][1], jacobian.vs1[2][2], jacobian.vs1[2][3], jacobian.vs1[2][4],
                    jacobian.vs2[0][0], jacobian.vs2[0][1], jacobian.vs2[0][2], jacobian.vs2[0][3], jacobian.vs2[0][4],
                    jacobian.vs2[1][0], jacobian.vs2[1][1], jacobian.vs2[1][2], jacobian.vs2[1][3], jacobian.vs2[1][4],
                    jacobian.vs2[2][0], jacobian.vs2[2][1], jacobian.vs2[2][2], jacobian.vs2[2][3], jacobian.vs2[2][4],
                );
            }

        if std::env::var("RUSTFOIL_BLDIF_DEBUG").is_ok()
            && ((surface == 0 && (i <= 2 || (20..=35).contains(&i) || (99..=105).contains(&i)))
                || (surface == 1 && (56..=63).contains(&i)))
        {
            eprintln!(
                "[BLDIF DEBUG] {}[{i}] flow={flow_type:?} x1={:.6e} x2={:.6e} ampl1={:.6e} ampl2={:.6e} ctau1={:.6e} ctau2={:.6e} vrez=[{:.6e}, {:.6e}, {:.6e}] vs2_row1=[{:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}]",
                if surface == 0 { "upper" } else { "lower" },
                stations[i - 1].x,
                s2.x,
                stations[i - 1].ampl,
                s2.ampl,
                stations[i - 1].ctau,
                s2.ctau,
                residuals.res_third,
                residuals.res_mom,
                residuals.res_shape,
                jacobian.vs2[0][0],
                jacobian.vs2[0][1],
                jacobian.vs2[0][2],
                jacobian.vs2[0][3],
                jacobian.vs2[0][4],
            );
                            if (24..=29).contains(&i)
                                || (surface == 0 && (99..=105).contains(&i))
                                || (surface == 1 && (56..=63).contains(&i))
            {
                for eq in 0..3 {
                    eprintln!(
                                        "[BLDIF DEBUG FULL] {}[{i}] row{} vs1=[{:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}] vs2=[{:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}]",
                                        if surface == 0 { "upper" } else { "lower" },
                        eq + 1,
                        jacobian.vs1[eq][0],
                        jacobian.vs1[eq][1],
                        jacobian.vs1[eq][2],
                        jacobian.vs1[eq][3],
                        jacobian.vs1[eq][4],
                        jacobian.vs2[eq][0],
                        jacobian.vs2[eq][1],
                        jacobian.vs2[eq][2],
                        jacobian.vs2[eq][3],
                        jacobian.vs2[eq][4],
                    );
                }
            }
        }


            // Store residuals
            self.vdel[iv] = [residuals.res_third, residuals.res_mom, residuals.res_shape];
            
            // Debug print moved AFTER VA/VB assignment to show stored values
            
            // Debug: print raw bldif residuals at first few stations (before forced changes)
            if std::env::var("RUSTFOIL_RAW_RES_DEBUG").is_ok() && i <= 3 {
                let surface_name = if surface == 0 { "upper" } else { "lower" };
                eprintln!("[RAW RES] {}[{}] IV={}: raw=[{:.6e}, {:.6e}, {:.6e}]", 
                    surface_name, i, iv, residuals.res_third, residuals.res_mom, residuals.res_shape);
            }
            
            // Debug: print Newton system at i=2 (XFOIL IBL=3) for comparison
            if std::env::var("RUSTFOIL_NEWTON_CMP").is_ok() && i == 2 && surface == 0 {
                eprintln!("\n=== RustFoil Newton System at i=2 (XFOIL IBL=3) ===");
                eprintln!("VS1 (upstream Jacobian, 3x5):");
                for eq in 0..3 {
                    eprintln!("  [{:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}]",
                        jacobian.vs1[eq][0], jacobian.vs1[eq][1], jacobian.vs1[eq][2], 
                        jacobian.vs1[eq][3], jacobian.vs1[eq][4]);
                }
                eprintln!("VS2 (downstream Jacobian, 3x5):");
                for eq in 0..3 {
                    eprintln!("  [{:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}]",
                        jacobian.vs2[eq][0], jacobian.vs2[eq][1], jacobian.vs2[eq][2], 
                        jacobian.vs2[eq][3], jacobian.vs2[eq][4]);
                }
                eprintln!("VDEL (raw residuals):");
                eprintln!("  [{:.6e}, {:.6e}, {:.6e}]", 
                    residuals.res_third, residuals.res_mom, residuals.res_shape);
            }
            
            // Debug: log residuals at transition stations
            if std::env::var("RUSTFOIL_TRANS_DEBUG").is_ok() {
                // Check if this is near the transition station
                let surface_name = if surface == 0 { "upper" } else { "lower" };
                let trans_idx = stations.iter().position(|s| !s.is_laminar);
                if let Some(trans) = trans_idx {
                    // Log stations from trans-2 to trans+5
                    if i >= trans.saturating_sub(2) && i <= trans + 5 {
                        let is_trans_station = i == trans;
                        let trans_mark = if is_trans_station { " <-- TRANSITION" } else { "" };
                        eprintln!("[VSREZ DEBUG] {surface_name}[{i}]: vsrez[0]={:.4e}, vsrez[1]={:.4e}, vsrez[2]={:.4e}, ctau={:.4e}, cq={:.4e}{trans_mark}",
                            residuals.res_third, residuals.res_mom, residuals.res_shape,
                            stations[i].ctau, stations[i].cq);
                    }
                }
            }
            
            // Emit JACOBIAN debug event for every 10th station
            if rustfoil_bl::is_debug_active() && (i % 10 == 0 || i <= 3) {
                // Convert 4x5 jacobian to 3x5 for debug output (skip 4th equation)
                let vs1_3x5: [[f64; 5]; 3] = [
                    jacobian.vs1[0],
                    jacobian.vs1[1],
                    jacobian.vs1[2],
                ];
                let vs2_3x5: [[f64; 5]; 3] = [
                    jacobian.vs2[0],
                    jacobian.vs2[1],
                    jacobian.vs2[2],
                ];
                let vsrez_3: [f64; 3] = [residuals.res_third, residuals.res_mom, residuals.res_shape];
                
                let flow_type_str = match flow_type {
                    FlowType::Laminar => "laminar",
                    FlowType::Turbulent => "turbulent",
                    FlowType::Wake => "wake",
                };
                
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::jacobian(
                    self.current_iteration,
                    i,
                    surface,
                    vs1_3x5,
                    vs2_3x5,
                    vsrez_3,
                    flow_type_str,
                ));
            }

            // Extract VA and VB blocks from Jacobian
            // IMPORTANT: XFOIL's VA only has 2 columns: [ampl/ctau, theta]
            // The delta_star/mass derivatives go ENTIRELY in VM, not VA!
            // 
            // XFOIL (xbl.f lines 340-341, 370-371, 400-401):
            //   VA(k,1,IV) = VS2(k,1)  -- ∂Fk/∂(ampl or ctau)
            //   VA(k,2,IV) = VS2(k,2)  -- ∂Fk/∂theta
            //
            // For similarity station (i==1), XFOIL combines: VS2 = VS1 + VS2, VS1 = 0
            // This DOUBLES the jacobian because at SIMI, VS1 and VS2 are identical
            // for columns 1-4 (theta, dstr, ue, x). Column 0 (ampl) is special.
            let finite_or_zero = |value: f64| if value.is_finite() { value } else { 0.0 };

            if i == 1 {
                // `bldif_full_simi()` already applies XFOIL's SIMI combine
                // (VS2 = VS1 + VS2, VS1 = 0), so store the returned blocks
                // directly without combining them a second time here.
                for eq in 0..3 {
                    // Only columns 0,1 go in VA (ampl/ctau and theta)
                    self.va[iv][eq][0] = finite_or_zero(jacobian.vs2[eq][0]);
                    self.va[iv][eq][1] = finite_or_zero(jacobian.vs2[eq][1]);
                    self.va[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    self.vb[iv][eq][0] = 0.0;
                    self.vb[iv][eq][1] = 0.0;
                    self.vb[iv][eq][2] = 0.0;
                    
                    // Store the already-combined delta_star and Ue columns for VM construction.
                    self.vs1_delta[iv][eq] = 0.0;
                    self.vs1_ue[iv][eq] = 0.0;
                    self.vs2_delta[iv][eq] = finite_or_zero(jacobian.vs2[eq][2]);
                    self.vs2_ue[iv][eq] = finite_or_zero(jacobian.vs2[eq][3]);
                    self.vs1_x[iv][eq] = 0.0;
                    self.vs2_x[iv][eq] =
                        finite_or_zero(jacobian.vs2[eq][4] + jacobian.vsx[eq]);
                }
            } else {
                for eq in 0..3 {
                    // Only columns 0,1 go in VA and VB (ampl/ctau and theta)
                    self.va[iv][eq][0] = finite_or_zero(jacobian.vs2[eq][0]);
                    self.va[iv][eq][1] = finite_or_zero(jacobian.vs2[eq][1]);
                    self.va[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    self.vb[iv][eq][0] = finite_or_zero(jacobian.vs1[eq][0]);
                    self.vb[iv][eq][1] = finite_or_zero(jacobian.vs1[eq][1]);
                    self.vb[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    
                    // Store delta_star and Ue columns for VM construction
                    self.vs1_delta[iv][eq] = finite_or_zero(jacobian.vs1[eq][2]);
                    self.vs1_ue[iv][eq] = finite_or_zero(jacobian.vs1[eq][3]);
                    self.vs2_delta[iv][eq] = finite_or_zero(jacobian.vs2[eq][2]);
                    self.vs2_ue[iv][eq] = finite_or_zero(jacobian.vs2[eq][3]);
                    
                    // Store X derivatives (column 4, 0-indexed) for stagnation coupling
                    // XFOIL: VS1(k,5), VS2(k,5) - derivatives w.r.t. arc length
                    self.vs1_x[iv][eq] = finite_or_zero(jacobian.vs1[eq][4]);
                    self.vs2_x[iv][eq] =
                        finite_or_zero(jacobian.vs2[eq][4] + jacobian.vsx[eq]);
                }
            }
            
            // Debug: print Newton matrix comparison data for first few stations
            // Placed AFTER VA/VB assignment to show actual stored values
            if std::env::var("RUSTFOIL_NEWTON_DEBUG").is_ok() && i <= 5 {
                let surface_name = if surface == 0 { "upper" } else { "lower" };
                eprintln!("\n[NEWTON DEBUG] {} IBL={} IV={}", surface_name, i, iv);
                eprintln!("  VA (stored in system):");
                eprintln!("    [0]: [{:>12.6e}, {:>12.6e}]", self.va[iv][0][0], self.va[iv][0][1]);
                eprintln!("    [1]: [{:>12.6e}, {:>12.6e}]", self.va[iv][1][0], self.va[iv][1][1]);
                eprintln!("    [2]: [{:>12.6e}, {:>12.6e}]", self.va[iv][2][0], self.va[iv][2][1]);
                eprintln!("  VB (stored in system):");
                eprintln!("    [0]: [{:>12.6e}, {:>12.6e}]", self.vb[iv][0][0], self.vb[iv][0][1]);
                eprintln!("    [1]: [{:>12.6e}, {:>12.6e}]", self.vb[iv][1][0], self.vb[iv][1][1]);
                eprintln!("    [2]: [{:>12.6e}, {:>12.6e}]", self.vb[iv][2][0], self.vb[iv][2][1]);
                eprintln!("  VSREZ (raw residuals):");
                eprintln!("    [0] res_ampl : {:>12.6e}", residuals.res_third);
                eprintln!("    [1] res_mom  : {:>12.6e}", residuals.res_mom);
                eprintln!("    [2] res_shape: {:>12.6e}", residuals.res_shape);
                if i == 1 {
                    eprintln!("  [SIMI combining: vs1[1][1]={:.6e} + vs2[1][1]={:.6e} = {:.6e}]",
                        jacobian.vs1[1][1], jacobian.vs2[1][1], jacobian.vs1[1][1] + jacobian.vs2[1][1]);
                }
            }
            
        }

        // NOTE: XFOIL does NOT zero the SIMI residuals - they're computed normally
        // with fixed log values (XLOG=1, ULOG=BULE, TLOG=0, HLOG=0).
        // The similarity condition is handled by combining Jacobians (VS2 = VS1 + VS2)
        // which is done in the build_surface_equations loop for i==1.
    }

    pub fn build_global_system_from_view(
        &mut self,
        state: &CanonicalNewtonStateView,
        dij: &DMatrix<f64>,
        ncrit: f64,
        msq: f64,
        re: f64,
        iteration: usize,
    ) {
        self.set_stagnation_derivs(state.sst_go, state.sst_gp);
        self.ante = state.ante;
        let upper_stations = state.upper_stations();
        let lower_stations = state.lower_stations();
        self.build_global_system(
            &upper_stations,
            &lower_stations,
            &state.upper_ue_current,
            &state.lower_ue_current,
            &state.upper_flow_types,
            &state.lower_flow_types,
            dij,
            &state.upper_ue_inviscid,
            &state.lower_ue_inviscid,
            &state.upper_ue_from_mass,
            &state.lower_ue_from_mass,
            &state.upper_ue_operating,
            &state.lower_ue_operating,
            ncrit,
            msq,
            re,
            iteration,
        );
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
            // XFOIL: DSI = MDI/UEI, D2_U2 = -DSI/UEI = -MASS/UEDG^2
            let dsi2 = if s2.u.abs() > 1e-20 { s2.mass_defect / s2.u } else { s2.delta_star };
            let d2_u2 = if s2.u.abs() > 1e-20 {
                -dsi2 / s2.u
            } else {
                0.0
            };
            let dsi1 = if s1.u.abs() > 1e-20 { s1.mass_defect / s1.u } else { s1.delta_star };
            let d1_u1 = if s1.u.abs() > 1e-20 {
                -dsi1 / s1.u
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

                let first_lower_wake = surface_i == 1 && ibl_i == self.iblte_lower + 1;
                let upstream_vti = if surface_i == 0 { 1.0 } else { -1.0 };
                let u1_m_j = if first_lower_wake {
                    0.0
                } else if ibl_i > 1 && panel_im1 < dij.nrows() && panel_j < dij.ncols() {
                    -upstream_vti * self.vti[jv] * dij[(panel_im1, panel_j)]
                } else {
                    0.0
                };

                // D_M[j] = derivative of delta_star w.r.t. mass defect at station j
                // Self-term: d(delta_star)/d(mass) = 1/Ue at same station
                let same_station_i = surface_i == surface_j && ibl_j == ibl_i;
                let same_station_im1 = surface_i == surface_j && ibl_j == ibl_i - 1;

                let d2_m_j = if same_station_i && s2.u.abs() > 1e-20 {
                    1.0 / s2.u + d2_u2 * u2_m_j
                } else {
                    d2_u2 * u2_m_j
                };

                let d1_m_j = if first_lower_wake {
                    if let Some(te) = self.te_wake_terms(upper_stations, lower_stations) {
                        let upper_te = &upper_stations[self.iblte_upper];
                        let lower_te = &lower_stations[self.iblte_lower];
                        let ivte_upper = self.to_global(0, self.iblte_upper);
                        let ivte_lower = self.to_global(1, self.iblte_lower);

                        let ute1_m_j =
                            if upper_te.panel_idx < dij.nrows() && panel_j < dij.ncols() {
                                -1.0 * self.vti[jv] * dij[(upper_te.panel_idx, panel_j)]
                            } else {
                                0.0
                            };
                        let ute2_m_j =
                            if lower_te.panel_idx < dij.nrows() && panel_j < dij.ncols() {
                                -(-1.0) * self.vti[jv] * dij[(lower_te.panel_idx, panel_j)]
                            } else {
                                0.0
                            };

                        let mut d1 =
                            te.dte_ute1 * ute1_m_j + te.dte_ute2 * ute2_m_j;
                        if jv == ivte_upper {
                            d1 += te.dte_mte1;
                        }
                        if jv == ivte_lower {
                            d1 += te.dte_mte2;
                        }
                        d1
                    } else {
                        0.0
                    }
                } else if same_station_im1 && ibl_i > 1 && s1.u.abs() > 1e-20 {
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
                    
                    // === Stagnation point coupling (XFOIL XI_ULE terms for VM) ===
                    // XFOIL: + (VS1(k,5) + VS2(k,5) + VSX(k)) * (XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
                    // ULE1_M(JV) = -VTI(2,1)*VTI(JBL,JS)*DIJ(ILE1,J)  
                    // ULE2_M(JV) = -VTI(2,2)*VTI(JBL,JS)*DIJ(ILE2,J)
                    // where ILE1/ILE2 are panel indices of first station after stagnation
                    
                    // Get panel indices for leading edge stations
                    let ile1 = if upper_stations.len() > 1 {
                        upper_stations[1].panel_idx
                    } else {
                        0
                    };
                    let ile2 = if lower_stations.len() > 1 {
                        lower_stations[1].panel_idx
                    } else {
                        0
                    };
                    
                    // Compute ULE_M derivatives
                    // VTI(2,1) = +1.0 (first station after stagnation on upper)
                    // VTI(2,2) = -1.0 (first station after stagnation on lower)
                    let vti_le1 = 1.0;  // Upper surface
                    let vti_le2 = -1.0; // Lower surface
                    
                    let ule1_m_j = if ile1 < dij.nrows() && panel_j < dij.ncols() {
                        -vti_le1 * self.vti[jv] * dij[(ile1, panel_j)]
                    } else {
                        0.0
                    };
                    
                    let ule2_m_j = if ile2 < dij.nrows() && panel_j < dij.ncols() {
                        -vti_le2 * self.vti[jv] * dij[(ile2, panel_j)]
                    } else {
                        0.0
                    };
                    
                    // XI_ULE depends on which surface we're on
                    // Upper surface (IS=1): XI_ULE1 = SST_GO, XI_ULE2 = -SST_GP
                    // Lower surface (IS=2): XI_ULE1 = -SST_GO, XI_ULE2 = SST_GP
                    let (xi_ule1, xi_ule2) = if surface_i == 0 {
                        (self.sst_go, -self.sst_gp)
                    } else {
                        (-self.sst_go, self.sst_gp)
                    };
                    
                    // Add stagnation coupling term
                    let vs_x = self.vs1_x[iv][k] + self.vs2_x[iv][k]; // VSX typically 0
                    self.vm[iv][jv][k] += vs_x * (xi_ule1 * ule1_m_j + xi_ule2 * ule2_m_j);
                }
            }
            
            // Emit VM_BLOCK debug event for every 10th station
            if rustfoil_bl::is_debug_active() && iv % 10 == 0 {
                // Collect near-diagonal VM entries (|IV - JV| <= 5)
                let mut vm_near: Vec<(usize, [f64; 3])> = Vec::new();
                let jv_start = if iv > 5 { iv - 5 } else { 1 };
                let jv_end = (iv + 5).min(self.nsys);
                for jv in jv_start..=jv_end {
                    if jv != iv {
                        vm_near.push((jv, self.vm[iv][jv]));
                    }
                }
                
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::vm_block(
                    self.current_iteration,
                    iv,
                    self.vm[iv][iv],
                    vm_near,
                ));
            }
        }
    }

    /// Build VZ block at trailing edge
    ///
    /// The VZ block couples the upper TE to the lower wake start.
    /// This handles the discontinuity where upper and lower surfaces meet.
    fn build_vz_block(&mut self, upper_stations: &[BlStation], lower_stations: &[BlStation]) {
        let iv_upper_te = self.to_global(0, self.iblte_upper);
        let jv_lower_wake = self.to_global(1, self.iblte_lower + 1);
        if iv_upper_te > self.nsys || jv_lower_wake > self.nsys {
            return;
        }

        if let Some(te) = self.te_wake_terms(upper_stations, lower_stations) {
            self.vz = [[0.0; 2]; 3];
            self.vz[0][0] = -te.cte_cte1;
            self.vz[0][1] = -te.cte_tte1;
            self.vz[1][1] = -1.0;
        }
    }

    /// Add forced changes to residuals
    ///
    /// The forced changes account for the mismatch between current Ue values
    /// and what UESET would compute from current mass defect.
    ///
    /// This also includes the stagnation point coupling (XI_ULE terms) from XFOIL,
    /// which accounts for how leading edge Ue changes affect the entire BL solution.
    fn add_forced_changes(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        ue_current_upper: &[f64],
        ue_current_lower: &[f64],
        ue_from_mass_upper: &[f64],
        ue_from_mass_lower: &[f64],
        ue_inviscid_upper: &[f64],
        ue_inviscid_lower: &[f64],
    ) {
        // Compute DUE for upper surface
        // XFOIL swaps after UESET: UEDG = march_Ue, USAV = mass_Ue
        // Then: DUE = UEDG - USAV = current - from_mass
        let due_upper: Vec<f64> = (0..upper_stations.len())
            .map(|i| {
                let curr = ue_current_upper.get(i).copied().unwrap_or(0.0);
                let from_mass = ue_from_mass_upper.get(i).copied().unwrap_or(0.0);
                curr - from_mass  // XFOIL: DUE = UEDG - USAV = current - from_mass
            })
            .collect();

        // Compute DUE for lower surface
        let due_lower: Vec<f64> = (0..lower_stations.len())
            .map(|i| {
                let curr = ue_current_lower.get(i).copied().unwrap_or(0.0);
                let from_mass = ue_from_mass_lower.get(i).copied().unwrap_or(0.0);
                curr - from_mass  // XFOIL: DUE = UEDG - USAV = current - from_mass
            })
            .collect();

        // === Stagnation point coupling (XFOIL XI_ULE terms) ===
        // DULE1, DULE2 are the forced changes in leading edge Ue (station 1 = first after stagnation)
        // XFOIL: DULE1 = UEDG(2,1) - USAV(2,1)  (IBL=2 is first station after stagnation)
        self.dule1 = due_upper.get(1).copied().unwrap_or(0.0);
        self.dule2 = due_lower.get(1).copied().unwrap_or(0.0);
        
        // Debug: print DUE for first few stations to understand the mismatch
        if std::env::var("RUSTFOIL_DUE_DEBUG").is_ok() {
            eprintln!("[DUE DEBUG] Upper surface DUE at first stations:");
            for i in 0..5.min(due_upper.len()) {
                let curr = ue_current_upper.get(i).copied().unwrap_or(0.0);
                let from_mass = ue_from_mass_upper.get(i).copied().unwrap_or(0.0);
                eprintln!("  Upper[{}]: ue_curr={:.6}, ue_mass={:.6}, DUE={:.6e}", 
                    i, curr, from_mass, due_upper[i]);
            }
            eprintln!("[DUE DEBUG] Lower surface DUE at first stations:");
            for i in 0..5.min(due_lower.len()) {
                let curr = ue_current_lower.get(i).copied().unwrap_or(0.0);
                let from_mass = ue_from_mass_lower.get(i).copied().unwrap_or(0.0);
                eprintln!("  Lower[{}]: ue_curr={:.6}, ue_mass={:.6}, DUE={:.6e}", 
                    i, curr, from_mass, due_lower[i]);
            }
        }

        // Debug: print VDEL before forced changes for first few stations
        if std::env::var("RUSTFOIL_NEWTON_DEBUG").is_ok() {
            eprintln!("\n[NEWTON DEBUG] Before forced changes:");
            for ibl in 1..=5.min(upper_stations.len().saturating_sub(1)) {
                let iv = self.to_global(0, ibl);
                if iv <= self.nsys {
                    eprintln!("  upper IBL={} IV={}: VDEL=[{:>12.6e}, {:>12.6e}, {:>12.6e}]",
                        ibl, iv, self.vdel[iv][0], self.vdel[iv][1], self.vdel[iv][2]);
                }
            }
        }

        // Add forced changes to upper surface residuals
        // Start from IBL=1 (which is IV=1 in global indexing, XFOIL's IBL=2)
        for ibl in 1..upper_stations.len() {
            let iv = self.to_global(0, ibl);
            if iv > self.nsys {
                continue;
            }

            let s1 = &upper_stations[ibl - 1];
            let s2 = &upper_stations[ibl];

            // d(delta_star)/d(Ue) using DSI = MASS/UEDG (XFOIL xbl.f:194,204-205)
            let dsi2 = if s2.u.abs() > 1e-20 { s2.mass_defect / s2.u } else { s2.delta_star };
            let d2_u2 = if s2.u.abs() > 1e-20 {
                -dsi2 / s2.u
            } else {
                0.0
            };
            let dsi1 = if s1.u.abs() > 1e-20 { s1.mass_defect / s1.u } else { s1.delta_star };
            let d1_u1 = if s1.u.abs() > 1e-20 {
                -dsi1 / s1.u
            } else {
                0.0
            };

            let dds2 = d2_u2 * due_upper.get(ibl).copied().unwrap_or(0.0);
            let due2 = due_upper.get(ibl).copied().unwrap_or(0.0);
            let (dds1, due1) = if ibl == 1 {
                (0.0, 0.0)
            } else {
                (
                    d1_u1 * due_upper.get(ibl - 1).copied().unwrap_or(0.0),
                    due_upper.get(ibl - 1).copied().unwrap_or(0.0),
                )
            };

            if std::env::var("RUSTFOIL_FORCED_DEBUG").is_ok()
                && ((20..=35).contains(&ibl) || (99..=105).contains(&ibl))
            {
                let fc0 = self.vs1_delta[iv][0] * dds1
                    + self.vs2_delta[iv][0] * dds2
                    + self.vs1_ue[iv][0] * due1
                    + self.vs2_ue[iv][0] * due2;
                let fc1 = self.vs1_delta[iv][1] * dds1
                    + self.vs2_delta[iv][1] * dds2
                    + self.vs1_ue[iv][1] * due1
                    + self.vs2_ue[iv][1] * due2;
                let fc2 = self.vs1_delta[iv][2] * dds1
                    + self.vs2_delta[iv][2] * dds2
                    + self.vs1_ue[iv][2] * due1
                    + self.vs2_ue[iv][2] * due2;
                eprintln!(
                    "[FORCED DEBUG] upper[{ibl}] due1={:.12e} due2={:.12e} dds1={:.12e} dds2={:.12e} fc=[{:.12e}, {:.12e}, {:.12e}]",
                    due1, due2, dds1, dds2, fc0, fc1, fc2
                );
                eprintln!(
                    "[FORCED DEBUG UE] upper[{ibl}] ue_prev={:.12e} ue_prev_mass={:.12e} ue_curr={:.12e} ue_curr_mass={:.12e}",
                    ue_current_upper.get(ibl - 1).copied().unwrap_or(0.0),
                    ue_from_mass_upper.get(ibl - 1).copied().unwrap_or(0.0),
                    ue_current_upper.get(ibl).copied().unwrap_or(0.0),
                    ue_from_mass_upper.get(ibl).copied().unwrap_or(0.0),
                );
                let xi_ule1 = self.sst_go;
                let xi_ule2 = -self.sst_gp;
                let xi_term = xi_ule1 * self.dule1 + xi_ule2 * self.dule2;
                eprintln!(
                    "[FORCED DEBUG XI] upper[{ibl}] vs_x=[{:.12e}, {:.12e}, {:.12e}] xi_term={:.12e}",
                    self.vs1_x[iv][0] + self.vs2_x[iv][0],
                    self.vs1_x[iv][1] + self.vs2_x[iv][1],
                    self.vs1_x[iv][2] + self.vs2_x[iv][2],
                    xi_term,
                );
            }

            // Emit forced changes debug event at key stations and leading edge
            if rustfoil_bl::is_debug_active()
                && (ibl <= 3 || ibl % 10 == 0 || (99..=105).contains(&ibl))
            {
                let ue_from_mass_i = ue_from_mass_upper.get(ibl).copied().unwrap_or(0.0);
                let ue_inviscid_i = ue_inviscid_upper.get(ibl).copied().unwrap_or(0.0);
                let ue_current_i = ue_current_upper.get(ibl).copied().unwrap_or(0.0);
                let xi_ule1 = self.sst_go;
                let xi_ule2 = -self.sst_gp;
                let xi_term = xi_ule1 * self.dule1 + xi_ule2 * self.dule2;
                
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::forced_changes(
                    self.current_iteration,
                    ibl,
                    0,  // side = upper
                    upper_stations.get(ibl).map(|s| s.panel_idx).unwrap_or(usize::MAX),
                    due1,
                    due2,
                    dds1,
                    dds2,
                    ue_from_mass_i,
                    ue_inviscid_i,
                    ue_current_i,
                    self.dule1,
                    self.dule2,
                    self.sst_go,
                    self.sst_gp,
                    xi_term,
                ));
            }

            // Debug: print forced change contributions at station 2 (first iteration only)
            if std::env::var("RUSTFOIL_FC_DEBUG").is_ok() && ibl <= 2 {
                eprintln!("[FC DEBUG] upper[{ibl}]: due1={:.6e}, due2={:.6e}, dds1={:.6e}, dds2={:.6e}",
                    due1, due2, dds1, dds2);
                eprintln!("[FC DEBUG]   vs1_ue=[{:.6e}, {:.6e}, {:.6e}]",
                    self.vs1_ue[iv][0], self.vs1_ue[iv][1], self.vs1_ue[iv][2]);
                eprintln!("[FC DEBUG]   vs2_ue=[{:.6e}, {:.6e}, {:.6e}]",
                    self.vs2_ue[iv][0], self.vs2_ue[iv][1], self.vs2_ue[iv][2]);
                let fc0 = self.vs1_delta[iv][0] * dds1 + self.vs2_delta[iv][0] * dds2
                        + self.vs1_ue[iv][0] * due1 + self.vs2_ue[iv][0] * due2;
                let fc1 = self.vs1_delta[iv][1] * dds1 + self.vs2_delta[iv][1] * dds2
                        + self.vs1_ue[iv][1] * due1 + self.vs2_ue[iv][1] * due2;
                let fc2 = self.vs1_delta[iv][2] * dds1 + self.vs2_delta[iv][2] * dds2
                        + self.vs1_ue[iv][2] * due1 + self.vs2_ue[iv][2] * due2;
                eprintln!("[FC DEBUG]   forced_change=[{:.6e}, {:.6e}, {:.6e}]", fc0, fc1, fc2);
            }
            
            for k in 0..3 {
                self.vdel[iv][k] += self.vs1_delta[iv][k] * dds1
                    + self.vs2_delta[iv][k] * dds2
                    + self.vs1_ue[iv][k] * due1
                    + self.vs2_ue[iv][k] * due2;
                
                // === Stagnation point coupling (XFOIL XI_ULE terms) ===
                // XFOIL: + (VS1(k,5) + VS2(k,5) + VSX(k)) * (XI_ULE1*DULE1 + XI_ULE2*DULE2)
                // For upper surface (IS=1): XI_ULE1 = SST_GO, XI_ULE2 = -SST_GP
                // VSX is usually 0 for turbulent flow, we skip it for now
                let xi_ule1 = self.sst_go;
                let xi_ule2 = -self.sst_gp;
                let vs_x = self.vs1_x[iv][k] + self.vs2_x[iv][k]; // VSX typically 0
                self.vdel[iv][k] += vs_x * (xi_ule1 * self.dule1 + xi_ule2 * self.dule2);
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

            let dsi2 = if s2.u.abs() > 1e-20 { s2.mass_defect / s2.u } else { s2.delta_star };
            let d2_u2 = if s2.u.abs() > 1e-20 {
                -dsi2 / s2.u
            } else {
                0.0
            };
            let dsi1 = if s1.u.abs() > 1e-20 { s1.mass_defect / s1.u } else { s1.delta_star };
            let d1_u1 = if s1.u.abs() > 1e-20 {
                -dsi1 / s1.u
            } else {
                0.0
            };

            let dds2 = d2_u2 * due_lower.get(ibl).copied().unwrap_or(0.0);
            let (dds1, due1) = if ibl == self.iblte_lower + 1 {
                if let Some(te) = self.te_wake_terms(upper_stations, lower_stations) {
                    let due_te_upper = due_upper.get(self.iblte_upper).copied().unwrap_or(0.0);
                    let due_te_lower = due_lower.get(self.iblte_lower).copied().unwrap_or(0.0);
                    (
                        te.dte_ute1 * due_te_upper + te.dte_ute2 * due_te_lower,
                        0.0,
                    )
                } else {
                    (d1_u1 * due_lower.get(ibl - 1).copied().unwrap_or(0.0), 0.0)
                }
            } else {
                (
                    d1_u1 * due_lower.get(ibl - 1).copied().unwrap_or(0.0),
                    due_lower.get(ibl - 1).copied().unwrap_or(0.0),
                )
            };
            let due2 = due_lower.get(ibl).copied().unwrap_or(0.0);

            // Emit forced changes debug event at key stations and leading edge
            if rustfoil_bl::is_debug_active() && (ibl <= 3 || ibl % 10 == 0) {
                let ue_from_mass_i = ue_from_mass_lower.get(ibl).copied().unwrap_or(0.0);
                let ue_inviscid_i = ue_inviscid_lower.get(ibl).copied().unwrap_or(0.0);
                let ue_current_i = ue_current_lower.get(ibl).copied().unwrap_or(0.0);
                let xi_ule1 = -self.sst_go;
                let xi_ule2 = self.sst_gp;
                let xi_term = xi_ule1 * self.dule1 + xi_ule2 * self.dule2;
                
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::forced_changes(
                    self.current_iteration,
                    ibl,
                    1,  // side = lower
                    lower_stations.get(ibl).map(|s| s.panel_idx).unwrap_or(usize::MAX),
                    due1,
                    due2,
                    dds1,
                    dds2,
                    ue_from_mass_i,
                    ue_inviscid_i,
                    ue_current_i,
                    self.dule1,
                    self.dule2,
                    self.sst_go,
                    self.sst_gp,
                    xi_term,
                ));
            }

            if std::env::var("RUSTFOIL_FORCED_DEBUG").is_ok() && (56..=63).contains(&ibl) {
                let fc0 = self.vs1_delta[iv][0] * dds1
                    + self.vs2_delta[iv][0] * dds2
                    + self.vs1_ue[iv][0] * due1
                    + self.vs2_ue[iv][0] * due2;
                let fc1 = self.vs1_delta[iv][1] * dds1
                    + self.vs2_delta[iv][1] * dds2
                    + self.vs1_ue[iv][1] * due1
                    + self.vs2_ue[iv][1] * due2;
                let fc2 = self.vs1_delta[iv][2] * dds1
                    + self.vs2_delta[iv][2] * dds2
                    + self.vs1_ue[iv][2] * due1
                    + self.vs2_ue[iv][2] * due2;
                eprintln!(
                    "[FORCED DEBUG LOWER] lower[{ibl}] due1={:.12e} due2={:.12e} dds1={:.12e} dds2={:.12e} fc=[{:.12e}, {:.12e}, {:.12e}]",
                    due1, due2, dds1, dds2, fc0, fc1, fc2
                );
                eprintln!(
                    "[FORCED DEBUG LOWER UE] lower[{ibl}] ue_prev={:.12e} ue_prev_mass={:.12e} ue_curr={:.12e} ue_curr_mass={:.12e}",
                    ue_current_lower.get(ibl - 1).copied().unwrap_or(0.0),
                    ue_from_mass_lower.get(ibl - 1).copied().unwrap_or(0.0),
                    ue_current_lower.get(ibl).copied().unwrap_or(0.0),
                    ue_from_mass_lower.get(ibl).copied().unwrap_or(0.0),
                );
                let xi_ule1 = -self.sst_go;
                let xi_ule2 = self.sst_gp;
                let xi_term = xi_ule1 * self.dule1 + xi_ule2 * self.dule2;
                eprintln!(
                    "[FORCED DEBUG LOWER XI] lower[{ibl}] vs_x=[{:.12e}, {:.12e}, {:.12e}] xi_term={:.12e}",
                    self.vs1_x[iv][0] + self.vs2_x[iv][0],
                    self.vs1_x[iv][1] + self.vs2_x[iv][1],
                    self.vs1_x[iv][2] + self.vs2_x[iv][2],
                    xi_term,
                );
            }

            for k in 0..3 {
                self.vdel[iv][k] += self.vs1_delta[iv][k] * dds1
                    + self.vs2_delta[iv][k] * dds2
                    + self.vs1_ue[iv][k] * due1
                    + self.vs2_ue[iv][k] * due2;
                
                // === Stagnation point coupling (XFOIL XI_ULE terms) ===
                // XFOIL: + (VS1(k,5) + VS2(k,5) + VSX(k)) * (XI_ULE1*DULE1 + XI_ULE2*DULE2)
                // For lower surface (IS=2): XI_ULE1 = -SST_GO, XI_ULE2 = SST_GP
                // VSX is usually 0 for turbulent flow, we skip it for now
                let xi_ule1 = -self.sst_go;
                let xi_ule2 = self.sst_gp;
                let vs_x = self.vs1_x[iv][k] + self.vs2_x[iv][k]; // VSX typically 0
                self.vdel[iv][k] += vs_x * (xi_ule1 * self.dule1 + xi_ule2 * self.dule2);
            }
        }
    }

    fn build_operating_rhs(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        ue_operating_upper: &[f64],
        ue_operating_lower: &[f64],
    ) {
        let dule1 = ue_operating_upper.get(1).copied().unwrap_or(0.0);
        let dule2 = ue_operating_lower.get(1).copied().unwrap_or(0.0);

        for ibl in 1..upper_stations.len() {
            let iv = self.to_global(0, ibl);
            if iv > self.nsys {
                continue;
            }

            let s1 = &upper_stations[ibl - 1];
            let s2 = &upper_stations[ibl];
            let dsi2 = if s2.u.abs() > 1e-20 { s2.mass_defect / s2.u } else { s2.delta_star };
            let d2_u2 = if s2.u.abs() > 1e-20 { -dsi2 / s2.u } else { 0.0 };
            let dsi1 = if s1.u.abs() > 1e-20 { s1.mass_defect / s1.u } else { s1.delta_star };
            let d1_u1 = if s1.u.abs() > 1e-20 { -dsi1 / s1.u } else { 0.0 };

            let due2 = ue_operating_upper.get(ibl).copied().unwrap_or(0.0);
            let due1 = ue_operating_upper.get(ibl - 1).copied().unwrap_or(0.0);
            let dds2 = d2_u2 * due2;
            let dds1 = d1_u1 * due1;

            for k in 0..3 {
                self.vdel_operating[iv][k] = self.vs1_delta[iv][k] * dds1
                    + self.vs2_delta[iv][k] * dds2
                    + self.vs1_ue[iv][k] * due1
                    + self.vs2_ue[iv][k] * due2;

                let xi_ule1 = self.sst_go;
                let xi_ule2 = -self.sst_gp;
                let vs_x = self.vs1_x[iv][k] + self.vs2_x[iv][k];
                self.vdel_operating[iv][k] += vs_x * (xi_ule1 * dule1 + xi_ule2 * dule2);
            }
        }

        for ibl in 1..lower_stations.len() {
            let iv = self.to_global(1, ibl);
            if iv > self.nsys {
                continue;
            }

            let s1 = &lower_stations[ibl - 1];
            let s2 = &lower_stations[ibl];
            let dsi2 = if s2.u.abs() > 1e-20 { s2.mass_defect / s2.u } else { s2.delta_star };
            let d2_u2 = if s2.u.abs() > 1e-20 { -dsi2 / s2.u } else { 0.0 };
            let dsi1 = if s1.u.abs() > 1e-20 { s1.mass_defect / s1.u } else { s1.delta_star };
            let d1_u1 = if s1.u.abs() > 1e-20 { -dsi1 / s1.u } else { 0.0 };

            let due2 = ue_operating_lower.get(ibl).copied().unwrap_or(0.0);
            let (dds1, due1) = if ibl == self.iblte_lower + 1 {
                if let Some(te) = self.te_wake_terms(upper_stations, lower_stations) {
                    let due_te_upper = ue_operating_upper.get(self.iblte_upper).copied().unwrap_or(0.0);
                    let due_te_lower = ue_operating_lower.get(self.iblte_lower).copied().unwrap_or(0.0);
                    (
                        te.dte_ute1 * due_te_upper + te.dte_ute2 * due_te_lower,
                        0.0,
                    )
                } else {
                    let due_prev = ue_operating_lower.get(ibl - 1).copied().unwrap_or(0.0);
                    (d1_u1 * due_prev, due_prev)
                }
            } else {
                let due_prev = ue_operating_lower.get(ibl - 1).copied().unwrap_or(0.0);
                (d1_u1 * due_prev, due_prev)
            };
            let dds2 = d2_u2 * due2;

            for k in 0..3 {
                self.vdel_operating[iv][k] = self.vs1_delta[iv][k] * dds1
                    + self.vs2_delta[iv][k] * dds2
                    + self.vs1_ue[iv][k] * due1
                    + self.vs2_ue[iv][k] * due2;

                let xi_ule1 = -self.sst_go;
                let xi_ule2 = self.sst_gp;
                let vs_x = self.vs1_x[iv][k] + self.vs2_x[iv][k];
                self.vdel_operating[iv][k] += vs_x * (xi_ule1 * dule1 + xi_ule2 * dule2);
            }
        }
    }

    fn te_wake_terms(
        &self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
    ) -> Option<TeWakeTerms> {
        if self.iblte_upper >= upper_stations.len()
            || self.iblte_lower >= lower_stations.len()
            || self.iblte_lower + 1 >= lower_stations.len()
        {
            return None;
        }

        let upper_te = &upper_stations[self.iblte_upper];
        let lower_te = &lower_stations[self.iblte_lower];
        let tte = (upper_te.theta + lower_te.theta).max(1.0e-12);
        // XFOIL TESYS (xblsys.f line 713): DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
        let dte = upper_te.delta_star + lower_te.delta_star + self.ante;
        let cte = (upper_te.ctau * upper_te.theta + lower_te.ctau * lower_te.theta) / tte;

        Some(TeWakeTerms {
            tte,
            dte,
            cte,
            dte_mte1: 1.0 / upper_te.u.max(1.0e-12),
            dte_ute1: -upper_te.delta_star / upper_te.u.max(1.0e-12),
            dte_mte2: 1.0 / lower_te.u.max(1.0e-12),
            dte_ute2: -lower_te.delta_star / lower_te.u.max(1.0e-12),
            cte_cte1: upper_te.theta / tte,
            cte_cte2: lower_te.theta / tte,
            cte_tte1: (upper_te.ctau - cte) / tte,
            cte_tte2: (lower_te.ctau - cte) / tte,
        })
    }

    fn apply_te_wake_interface(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
    ) {
        let Some(te) = self.te_wake_terms(upper_stations, lower_stations) else {
            return;
        };
        let wake_ibl = self.iblte_lower + 1;
        if wake_ibl >= lower_stations.len() {
            return;
        }

        let iv = self.to_global(1, wake_ibl);
        if iv > self.nsys {
            return;
        }

        let wake_station = &lower_stations[wake_ibl];

        self.vdel[iv] = [
            te.cte - wake_station.ctau,
            te.tte - wake_station.theta,
            te.dte - (wake_station.delta_star + wake_station.dw),
        ];

        if std::env::var("RUSTFOIL_TESYS_DEBUG").is_ok() {
            eprintln!(
                "[TESYS DEBUG] iv={} wake_ibl={} upper_te(ctau={:.12e}, theta={:.12e}, dstar={:.12e}, u={:.12e}) lower_te(ctau={:.12e}, theta={:.12e}, dstar={:.12e}, u={:.12e}) wake(ctau={:.12e}, theta={:.12e}, dstar={:.12e}, dw={:.12e}, u={:.12e}) cte={:.12e} tte={:.12e} dte={:.12e} base_vdel=[{:.12e}, {:.12e}, {:.12e}]",
                iv,
                wake_ibl,
                upper_stations[self.iblte_upper].ctau,
                upper_stations[self.iblte_upper].theta,
                upper_stations[self.iblte_upper].delta_star,
                upper_stations[self.iblte_upper].u,
                lower_stations[self.iblte_lower].ctau,
                lower_stations[self.iblte_lower].theta,
                lower_stations[self.iblte_lower].delta_star,
                lower_stations[self.iblte_lower].u,
                wake_station.ctau,
                wake_station.theta,
                wake_station.delta_star,
                wake_station.dw,
                wake_station.u,
                te.cte,
                te.tte,
                te.dte,
                self.vdel[iv][0],
                self.vdel[iv][1],
                self.vdel[iv][2],
            );
        }

        self.va[iv] = [[0.0; 3]; 3];
        self.vb[iv] = [[0.0; 3]; 3];

        self.va[iv][0][0] = 1.0;
        self.va[iv][1][1] = 1.0;

        self.vb[iv][0][0] = -te.cte_cte2;
        self.vb[iv][0][1] = -te.cte_tte2;
        self.vb[iv][1][1] = -1.0;

        self.vs1_delta[iv] = [0.0, 0.0, -1.0];
        self.vs1_ue[iv] = [0.0; 3];
        self.vs2_delta[iv] = [0.0, 0.0, 1.0];
        self.vs2_ue[iv] = [0.0; 3];
        self.vs1_x[iv] = [0.0; 3];
        self.vs2_x[iv] = [0.0; 3];
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
/// Solution vectors for the primary state correction and the operating-variable
/// sensitivity column.
///
/// # XFOIL Reference
/// XFOIL xsolve.f BLSOLV (lines 283-485)
#[derive(Debug, Clone)]
pub struct GlobalSolveResult {
    pub state_deltas: Vec<[f64; 3]>,
    pub operating_deltas: Vec<[f64; 3]>,
}

pub fn solve_global_system(system: &mut GlobalNewtonSystem) -> GlobalSolveResult {
    let nsys = system.nsys;
    if nsys < 2 {
        return GlobalSolveResult {
            state_deltas: vec![[0.0; 3]; nsys + 1],
            operating_deltas: vec![[0.0; 3]; nsys + 1],
        };
    }


    // XFOIL VACCEL = 0.01 (xfoil.f line 701)
    let vaccel1 = 0.01;
    let vaccel23 = 0.01 * 2.0 / system.airfoil_arc_length.max(1.0);

    // Working copies
    let mut vdel = system.vdel.clone();
    let mut vdel_operating = system.vdel_operating.clone();
    let mut va_mod = system.va.clone();
    let mut vm_mod = system.vm.clone();

    // Global index of upper TE (where VZ block applies)
    let ivte = system.to_global(0, system.iblte_upper);
    // Global index of lower wake start
    let ivz_target = system.to_global(1, system.iblte_lower + 1);

    // === Forward Sweep ===
    for iv in 1..=nsys {
        let ivp = iv + 1;

        // === Invert VA block ===
        // XFOIL divides unconditionally (xsolve.f line 327). For robustness,
        // we use a safe inverse that returns 0 for near-zero pivots rather than
        // skipping the entire station (which would corrupt all downstream rows).
        let safe_inv = |x: f64| -> f64 {
            if x.abs() < 1e-25 { 0.0 } else { 1.0 / x }
        };

        // Normalize first row by VA[0][0]
        let pivot_inv = safe_inv(va_mod[iv][0][0]);
        va_mod[iv][0][1] *= pivot_inv;
        for l in iv..=nsys {
            vm_mod[iv][l][0] *= pivot_inv;
        }
        vdel[iv][0] *= pivot_inv;
        vdel_operating[iv][0] *= pivot_inv;

        // Eliminate lower first column in VA block (rows 2, 3)
        for k in 1..3 {
            let vtmp = va_mod[iv][k][0];
            va_mod[iv][k][1] -= vtmp * va_mod[iv][0][1];
            for l in iv..=nsys {
                vm_mod[iv][l][k] -= vtmp * vm_mod[iv][l][0];
            }
            vdel[iv][k] -= vtmp * vdel[iv][0];
            vdel_operating[iv][k] -= vtmp * vdel_operating[iv][0];
        }

        // Normalize second row by VA[1][1]
        let pivot_inv = safe_inv(va_mod[iv][1][1]);
        for l in iv..=nsys {
            vm_mod[iv][l][1] *= pivot_inv;
        }
        vdel[iv][1] *= pivot_inv;
        vdel_operating[iv][1] *= pivot_inv;

        // Eliminate lower second column (row 3 only)
        let vtmp = va_mod[iv][2][1];
        for l in iv..=nsys {
            vm_mod[iv][l][2] -= vtmp * vm_mod[iv][l][1];
        }
        vdel[iv][2] -= vtmp * vdel[iv][1];
        vdel_operating[iv][2] -= vtmp * vdel_operating[iv][1];

        // Normalize third row by VM(3,IV,IV)
        let pivot_inv = safe_inv(vm_mod[iv][iv][2]);
        for l in ivp..=nsys {
            vm_mod[iv][l][2] *= pivot_inv;
        }
        vdel[iv][2] *= pivot_inv;
        vdel_operating[iv][2] *= pivot_inv;

        // Back-substitute rows 1-2 from row 3
        let vtmp1 = vm_mod[iv][iv][0];
        let vtmp2 = vm_mod[iv][iv][1];
        for l in ivp..=nsys {
            vm_mod[iv][l][0] -= vtmp1 * vm_mod[iv][l][2];
            vm_mod[iv][l][1] -= vtmp2 * vm_mod[iv][l][2];
        }
        vdel[iv][0] -= vtmp1 * vdel[iv][2];
        vdel[iv][1] -= vtmp2 * vdel[iv][2];
        vdel_operating[iv][0] -= vtmp1 * vdel_operating[iv][2];
        vdel_operating[iv][1] -= vtmp2 * vdel_operating[iv][2];

        // Back-substitute row 1 from row 2
        let vtmp = va_mod[iv][0][1];
        for l in ivp..=nsys {
            vm_mod[iv][l][0] -= vtmp * vm_mod[iv][l][1];
        }
        vdel[iv][0] -= vtmp * vdel[iv][1];
        vdel_operating[iv][0] -= vtmp * vdel_operating[iv][1];

        if iv >= nsys {
            continue;
        }

        // === Eliminate VB[IV+1] block (XFOIL does VB first, then VZ) ===
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
            vdel_operating[ivp][k] -= vtmp1 * vdel_operating[iv][0]
                + vtmp2 * vdel_operating[iv][1]
                + vtmp3 * vdel_operating[iv][2];
        }

        // === Handle VZ block at upper TE (after VB, matching XFOIL xsolve.f line 418) ===
        if iv == ivte && ivz_target <= nsys {
            for k in 0..3 {
                let vtmp1 = system.vz[k][0];
                let vtmp2 = system.vz[k][1];
                for l in ivp..=nsys {
                    vm_mod[ivz_target][l][k] -= vtmp1 * vm_mod[iv][l][0] + vtmp2 * vm_mod[iv][l][1];
                }
                vdel[ivz_target][k] -= vtmp1 * vdel[iv][0] + vtmp2 * vdel[iv][1];
                vdel_operating[ivz_target][k] -=
                    vtmp1 * vdel_operating[iv][0] + vtmp2 * vdel_operating[iv][1];
            }
        }

        if ivp >= nsys {
            continue;
        }

        // === Eliminate lower VM column (sparse, with VACCEL threshold) ===
        for kv in (iv + 2)..=nsys {
            let vtmp1 = vm_mod[kv][iv][0];
            let vtmp2 = vm_mod[kv][iv][1];
            let vtmp3 = vm_mod[kv][iv][2];

            if vtmp1.abs() > vaccel1 {
                for l in ivp..=nsys {
                    vm_mod[kv][l][0] -= vtmp1 * vm_mod[iv][l][2];
                }
                vdel[kv][0] -= vtmp1 * vdel[iv][2];
                vdel_operating[kv][0] -= vtmp1 * vdel_operating[iv][2];
            }
            if vtmp2.abs() > vaccel23 {
                for l in ivp..=nsys {
                    vm_mod[kv][l][1] -= vtmp2 * vm_mod[iv][l][2];
                }
                vdel[kv][1] -= vtmp2 * vdel[iv][2];
                vdel_operating[kv][1] -= vtmp2 * vdel_operating[iv][2];
            }
            if vtmp3.abs() > vaccel23 {
                for l in ivp..=nsys {
                    vm_mod[kv][l][2] -= vtmp3 * vm_mod[iv][l][2];
                }
                vdel[kv][2] -= vtmp3 * vdel[iv][2];
                vdel_operating[kv][2] -= vtmp3 * vdel_operating[iv][2];
            }
        }

        if std::env::var("RUSTFOIL_BLSOLV_STEP_DEBUG").is_ok()
            && ((20..=26).contains(&iv) || (94..=100).contains(&iv) || (156..=160).contains(&iv))
        {
            eprintln!("[RUST BLSOLV STEP] iv={iv}");
            for row in [24usize, 25, 26, 98, 99, 159, 160] {
                if row <= nsys {
                    eprintln!(
                        "[RUST BLSOLV ROW] row={} vals=[{:.8e}, {:.8e}, {:.8e}]",
                        row, vdel[row][0], vdel[row][1], vdel[row][2]
                    );
                }
            }
        }
    }

    // === Back Substitution ===
    // Eliminate upper VM columns using row 3 (mass) solution
    for iv in (2..=nsys).rev() {
        let vtmp = vdel[iv][2];

        if std::env::var("RUSTFOIL_BLSOLV_STEP_DEBUG").is_ok() && iv > 26 {
            eprintln!(
                "[RUST BLSOLV BACK MASS] iv={} coef={:.8e} mass={:.8e} contrib={:.8e}",
                iv,
                vm_mod[26][iv][2],
                vtmp,
                -vm_mod[26][iv][2] * vtmp
            );
            eprintln!(
                "[RUST BLSOLV BACK COEF26] iv={} vals=[{:.8e}, {:.8e}, {:.8e}]",
                iv,
                vm_mod[26][iv][0],
                vm_mod[26][iv][1],
                vm_mod[26][iv][2],
            );
        }

        if std::env::var("RUSTFOIL_BLSOLV_STEP_DEBUG").is_ok() && iv == 26 {
            eprintln!(
                "[RUST BLSOLV BACK COEF] row=24 col=26 vals=[{:.8e}, {:.8e}, {:.8e}] mass={:.8e}",
                vm_mod[24][26][0],
                vm_mod[24][26][1],
                vm_mod[24][26][2],
                vtmp
            );
            eprintln!(
                "[RUST BLSOLV BACK COEF] row=25 col=26 vals=[{:.8e}, {:.8e}, {:.8e}] mass={:.8e}",
                vm_mod[25][26][0],
                vm_mod[25][26][1],
                vm_mod[25][26][2],
                vtmp
            );
        }

        for kv in (1..iv).rev() {
            vdel[kv][0] -= vm_mod[kv][iv][0] * vtmp;
            vdel[kv][1] -= vm_mod[kv][iv][1] * vtmp;
            vdel[kv][2] -= vm_mod[kv][iv][2] * vtmp;
            vdel_operating[kv][0] -= vm_mod[kv][iv][0] * vdel_operating[iv][2];
            vdel_operating[kv][1] -= vm_mod[kv][iv][1] * vdel_operating[iv][2];
            vdel_operating[kv][2] -= vm_mod[kv][iv][2] * vdel_operating[iv][2];
        }

        if std::env::var("RUSTFOIL_BLSOLV_STEP_DEBUG").is_ok()
            && ((20..=26).contains(&iv) || (94..=100).contains(&iv) || (156..=160).contains(&iv))
        {
            eprintln!("[RUST BLSOLV BACK] iv={iv}");
            for row in [24usize, 25, 26, 98, 99, 159, 160] {
                if row <= nsys {
                    eprintln!(
                        "[RUST BLSOLV BACK ROW] row={} vals=[{:.8e}, {:.8e}, {:.8e}]",
                        row, vdel[row][0], vdel[row][1], vdel[row][2]
                    );
                }
            }
        }
    }

    if std::env::var("RUSTFOIL_FINAL_DELTA_DEBUG").is_ok() {
        let start = 156.min(nsys);
        for iv in start..=nsys {
            eprintln!(
                "[RUST FINAL DELTA] iv={} vals=[{:.12e}, {:.12e}, {:.12e}]",
                iv, vdel[iv][0], vdel[iv][1], vdel[iv][2]
            );
        }
    }

    GlobalSolveResult {
        state_deltas: vdel,
        operating_deltas: vdel_operating,
    }
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
    let nout = nsys;

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

/// Result of applying global updates
///
/// Contains the actual relaxation factor used and the RMS of normalized changes
/// (XFOIL's RMSBL metric for convergence).
#[derive(Debug, Clone)]
pub struct GlobalUpdateResult {
    /// Actual relaxation factor used (may be reduced from requested)
    pub relaxation_used: f64,
    /// RMS of normalized changes (XFOIL's RMSBL)
    /// This is sqrt(sum(DN1² + DN2² + DN3² + DN4²) / (4*N))
    pub rms_change: f64,
    /// Maximum normalized change across all stations
    pub max_change: f64,
    /// Trial Ue distribution before under-relaxation.
    pub upper_u_new: Vec<f64>,
    pub lower_u_new: Vec<f64>,
    /// Operating-column Ue sensitivity used in the DAC correction.
    pub upper_u_ac: Vec<f64>,
    pub lower_u_ac: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct GlobalOperatingPreview {
    pub upper_u_new: Vec<f64>,
    pub lower_u_new: Vec<f64>,
    pub upper_u_ac: Vec<f64>,
    pub lower_u_ac: Vec<f64>,
}

fn apply_dslim_like_xfoil(station: &mut BlStation, msq: f64) {
    let hklim = if station.is_wake { 1.00005 } else { 1.02 };
    let theta = station.theta.max(1.0e-12);
    let dsw = station.delta_star;
    let h = dsw / theta;
    let hkin_result = hkin(h, msq);

    if hkin_result.hk_h.abs() > 1.0e-12 {
        let dh = (hklim - hkin_result.hk).max(0.0) / hkin_result.hk_h;
        station.delta_star = dsw + dh * theta;
    }
}

fn total_dstr_for_update(station: &BlStation) -> f64 {
    if station.is_wake {
        station.delta_star + station.dw
    } else {
        station.delta_star
    }
}

fn total_dstr_for_row_update(row: &CanonicalNewtonRow) -> f64 {
    if row.is_wake {
        row.dstr + row.dw
    } else {
        row.dstr
    }
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
///
/// # Returns
/// `GlobalUpdateResult` containing the actual relaxation used and RMS of normalized changes
pub fn apply_global_updates(
    upper_stations: &mut [BlStation],
    lower_stations: &mut [BlStation],
    deltas: &[[f64; 3]],
    operating_deltas: Option<&[[f64; 3]]>,
    system: &GlobalNewtonSystem,
    iter: usize,
    relaxation: f64,
    dij: Option<&nalgebra::DMatrix<f64>>,
    upper_ue_inv: Option<&[f64]>,
    lower_ue_inv: Option<&[f64]>,
    upper_ue_operating: Option<&[f64]>,
    lower_ue_operating: Option<&[f64]>,
    dac: f64,
    msq: f64,
    re: f64,
) -> GlobalUpdateResult {
    // === XFOIL-style normalized relaxation (xbl.f lines 1527-1597) ===
    // Compute global relaxation factor based on normalized changes
    let dhi: f64 = 1.5;
    let dlo: f64 = -0.5;

    let mut rlx = relaxation.clamp(0.0, 1.0);
    
    // XFOIL forms UNEW before computing DN1..DN4 for the relaxation limit.
    let preview = preview_global_update_ue(
        upper_stations,
        lower_stations,
        deltas,
        operating_deltas,
        system,
        dij,
        upper_ue_inv,
        lower_ue_inv,
        upper_ue_operating,
        lower_ue_operating,
    );
    let upper_new_ue = preview.upper_u_new;
    let lower_new_ue = preview.lower_u_new;
    let upper_operating_ue = preview.upper_u_ac;
    let lower_operating_ue = preview.lower_u_ac;

    let corrected_delta = |iv: usize| -> [f64; 3] {
        if let Some(operating) = operating_deltas.filter(|_| iv < deltas.len()) {
            let op = operating.get(iv).copied().unwrap_or([0.0; 3]);
            [
                deltas[iv][0] - dac * op[0],
                deltas[iv][1] - dac * op[1],
                deltas[iv][2] - dac * op[2],
            ]
        } else if iv < deltas.len() {
            deltas[iv]
        } else {
            [0.0; 3]
        }
    };

    let should_log_update_newton = |surface: usize, ibl: usize| -> bool {
        if std::env::var("RUSTFOIL_UPDATE_NEWTON_DEBUG").is_err() {
            return false;
        }

        let side = surface + 1;
        if (6..=10).contains(&iter) {
            ((side == 2 && (54..=60).contains(&ibl))
                || (side == 2 && (2..=8).contains(&ibl))
                || (side == 1 && (18..=24).contains(&ibl))
                || (side == 1 && (26..=32).contains(&ibl))
                || (side == 1 && (98..=106).contains(&ibl)))
        } else {
            ibl % 20 == 0 && (iter <= 5 || (14..=20).contains(&iter))
        }
    };

    // First pass: compute required relaxation and RMSBL from DN1..DN4 exactly as XFOIL does.
    // XFOIL computes RMSBL and relaxation BEFORE applying any updates (P7).
    // XFOIL does NOT have 'if new_rlx < rlx' guards — it overwrites RLX directly (P3).
    let mut rms_sum = 0.0;
    let mut max_dn = 0.0_f64;
    let mut n_stations = 0usize;

    // Upper surface
    for ibl in 1..upper_stations.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = corrected_delta(iv);
        let station = &upper_stations[ibl];
        
        let dn1 = if station.is_laminar {
            delta[0] / 10.0
        } else if station.ctau.abs() > 1e-10 {
            delta[0] / station.ctau
        } else {
            delta[0] / 0.03
        };
        
        let dn2 = if station.theta.abs() > 1e-12 {
            delta[1] / station.theta
        } else {
            0.0
        };
        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u)
            + dac * upper_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - station.u;
        let total_dstr = total_dstr_for_update(station);
        let ddstr = (delta[2] - total_dstr * due) / station.u.max(1e-12);
        let dn3 = if total_dstr.abs() > 1e-12 {
            ddstr / total_dstr
        } else {
            0.0
        };
        let dn4 = due.abs() / 0.25;
        
        // Accumulate RMSBL (XFOIL xbl.f:1592)
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
        max_dn = max_dn.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());
        n_stations += 1;

        // Apply relaxation limits — XFOIL overwrites RLX directly, no min guard (P3)
        let rdn1 = rlx * dn1;
        if rdn1 > dhi { rlx = dhi / dn1; }
        if rdn1 < dlo { rlx = dlo / dn1; }
        
        let rdn2 = rlx * dn2;
        if rdn2 > dhi { rlx = dhi / dn2; }
        if rdn2 < dlo { rlx = dlo / dn2; }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi { rlx = dhi / dn3; }
        if rdn3 < dlo { rlx = dlo / dn3; }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi { rlx = dhi / dn4; }

        if rustfoil_bl::is_debug_active() && should_log_update_newton(0, ibl) {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::update_newton(
                iter,
                1,
                ibl,
                rlx,
                dn1,
                dn2,
                dn3,
                dn4,
                delta[0],
                delta[1],
                ddstr,
                due,
            ));
        }
    }
    
    // Lower surface
    for ibl in 1..lower_stations.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = corrected_delta(iv);
        let station = &lower_stations[ibl];
        
        let dn1 = if station.is_laminar {
            delta[0] / 10.0
        } else if station.ctau.abs() > 1e-10 {
            delta[0] / station.ctau
        } else {
            delta[0] / 0.03
        };
        
        let dn2 = if station.theta.abs() > 1e-12 {
            delta[1] / station.theta
        } else {
            0.0
        };
        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(station.u)
            + dac * lower_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - station.u;
        let total_dstr = total_dstr_for_update(station);
        let ddstr = (delta[2] - total_dstr * due) / station.u.max(1e-12);
        let dn3 = if total_dstr.abs() > 1e-12 {
            ddstr / total_dstr
        } else {
            0.0
        };
        let dn4 = due.abs() / 0.25;
        
        // Accumulate RMSBL
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
        max_dn = max_dn.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());
        n_stations += 1;

        // Apply relaxation limits — XFOIL overwrites RLX directly, no min guard (P3)
        let rdn1 = rlx * dn1;
        if rdn1 > dhi { rlx = dhi / dn1; }
        if rdn1 < dlo { rlx = dlo / dn1; }
        
        let rdn2 = rlx * dn2;
        if rdn2 > dhi { rlx = dhi / dn2; }
        if rdn2 < dlo { rlx = dlo / dn2; }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi { rlx = dhi / dn3; }
        if rdn3 < dlo { rlx = dlo / dn3; }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi { rlx = dhi / dn4; }

        if rustfoil_bl::is_debug_active() && should_log_update_newton(1, ibl) {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::update_newton(
                iter,
                2,
                ibl,
                rlx,
                dn1,
                dn2,
                dn3,
                dn4,
                delta[0],
                delta[1],
                ddstr,
                due,
            ));
        }
    }
    
    // XFOIL has no relaxation floor — just clamp to [0, 1]
    rlx = rlx.max(0.0).min(1.0);
    
    // Store computed relaxation before applying
    let rlx_computed = rlx;
    
    // Collect sample update data for debug output
    let mut sample_updates: Vec<rustfoil_bl::SampleUpdateData> = Vec::new();
    let mut limiting_station: Option<usize> = None;
    let mut limiting_reason: Option<String> = None;
    
    // Track which station limited relaxation most
    if rlx < relaxation * 0.99 {
        let mut best_limit = 0.0_f64;

        let mut scan_surface = |
            surface: usize,
            side_name: &str,
            stations: &[BlStation],
            new_ue: &[f64],
            operating_ue: &[f64],
        | {
            for ibl in 1..stations.len() {
                let iv = system.to_global(surface, ibl);
                if iv >= deltas.len() {
                    continue;
                }

                let station = &stations[ibl];
                let delta = corrected_delta(iv);
                let new_u = new_ue.get(ibl).copied().unwrap_or(station.u)
                    + dac * operating_ue.get(ibl).copied().unwrap_or(0.0);
                let due = new_u - station.u;
                let total_dstr = total_dstr_for_update(station);
                let d_dstar = (delta[2] - total_dstr * due) / station.u.max(1.0e-12);

                let dn1 = if station.is_laminar {
                    delta[0] / 10.0
                } else if station.ctau.abs() > 1e-10 {
                    delta[0] / station.ctau
                } else {
                    delta[0] / 0.03
                };
                let dn2 = if station.theta.abs() > 1e-12 {
                    delta[1] / station.theta
                } else {
                    0.0
                };
                let dn3 = if total_dstr.abs() > 1.0e-12 {
                    d_dstar / total_dstr
                } else {
                    0.0
                };
                let dn4 = due.abs() / 0.25;

                let candidates = [
                    ("ctau", "dn1", dn1),
                    ("theta", "dn2", dn2),
                    ("dstar", "dn3", dn3),
                    ("Ue", "dn4", dn4),
                ];

                for (quantity, dn_name, dn_value) in candidates {
                    let scaled = (rlx * dn_value).abs();
                    if scaled >= dhi && scaled > best_limit {
                        best_limit = scaled;
                        limiting_station = Some(ibl);
                        limiting_reason = Some(format!(
                            "{} {} {}={:.3}",
                            side_name, quantity, dn_name, dn_value
                        ));
                    }
                }
            }
        };

        scan_surface(0, "upper", &upper_stations, &upper_new_ue, &upper_operating_ue);
        scan_surface(1, "lower", &lower_stations, &lower_new_ue, &lower_operating_ue);
        if best_limit == 0.0 {
            limiting_station = None;
            limiting_reason = None;
        }
    }

    // === Step 2: Update upper surface BL variables ===
    for ibl in 1..upper_stations.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }

        let delta = corrected_delta(iv);
        let station = &mut upper_stations[ibl];

        // XFOIL: CTAU = CTAU + RLX*DCTAU, then MIN(CTAU, 0.25) for turbulent only
        if station.is_laminar {
            station.ampl += rlx * delta[0];
        } else {
            station.ctau += rlx * delta[0];
            station.ctau = station.ctau.min(0.25);
        }

        station.theta += rlx * delta[1];

        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u)
            + dac * upper_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - station.u;
        
        // XFOIL: DDSTR = (DMASS - DSTR*DUEDG) / UEDG
        let delta_mass = delta[2];
        let total_dstr = total_dstr_for_update(station);
        let d_dstar = (delta_mass - total_dstr * due) / station.u;
        
        // XFOIL: UEDG = UEDG + RLX*DUEDG
        station.u += rlx * due;
        
        // XFOIL evolves total DSTR, which includes the wake gap for wake rows.
        let mut total_dstr_after = total_dstr + rlx * d_dstar;
        station.delta_star = if station.is_wake {
            total_dstr_after - station.dw
        } else {
            total_dstr_after
        };

        // P5: DSLIM with WGAP — subtract wake gap before DSLIM, add back after
        let dswaki = station.dw;
        station.delta_star -= dswaki;
        apply_dslim_like_xfoil(station, msq);
        station.delta_star += dswaki;
        total_dstr_after = total_dstr_for_update(station);

        // XFOIL: MASS = DSTR * UEDG (nonlinear update, total DSTR for wake rows)
        station.h = if station.theta.abs() > 1e-20 { station.delta_star / station.theta } else { 2.0 };
        station.mass_defect = station.u * total_dstr_after;
    }

    // === Step 3: Update lower surface BL variables ===
    for ibl in 1..lower_stations.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }

        let delta = corrected_delta(iv);
        let station = &mut lower_stations[ibl];

        if station.is_laminar {
            station.ampl += rlx * delta[0];
        } else {
            station.ctau += rlx * delta[0];
            station.ctau = station.ctau.min(0.25);
        }

        station.theta += rlx * delta[1];

        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(station.u)
            + dac * lower_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - station.u;
        
        let delta_mass = delta[2];
        let total_dstr = total_dstr_for_update(station);
        let d_dstar = (delta_mass - total_dstr * due) / station.u;
        
        station.u += rlx * due;
        
        let mut total_dstr_after = total_dstr + rlx * d_dstar;
        station.delta_star = if station.is_wake {
            total_dstr_after - station.dw
        } else {
            total_dstr_after
        };

        // P5: DSLIM with WGAP
        let dswaki = station.dw;
        station.delta_star -= dswaki;
        apply_dslim_like_xfoil(station, msq);
        station.delta_star += dswaki;
        total_dstr_after = total_dstr_for_update(station);

        station.h = if station.theta.abs() > 1e-20 { station.delta_star / station.theta } else { 2.0 };
        station.mass_defect = station.u * total_dstr_after;
    }
    
    // Emit UPDATE_SUMMARY debug event with relaxation details
    if rustfoil_bl::is_debug_active() {
        // Collect sample update data for every 10th station (upper surface)
        // Skip station 0 (stagnation point) - start from ibl=1
        for ibl in (1..upper_stations.len()).filter(|i| i % 10 == 0) {
            let iv = system.to_global(0, ibl);
            if iv >= deltas.len() { continue; }
            
            let station = &upper_stations[ibl];
            let delta = corrected_delta(iv);
            let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u)
                + dac * upper_operating_ue.get(ibl).copied().unwrap_or(0.0);
            let due = new_u - station.u;
            
            let dn1 = if station.is_laminar {
                delta[0] / 10.0
            } else if station.ctau.abs() > 1e-10 {
                delta[0] / station.ctau
            } else {
                delta[0] / 0.03
            };
            let dn2 = if station.theta.abs() > 1e-12 { delta[1] / station.theta } else { 0.0 };
            let total_dstr = total_dstr_for_update(station);
            let d_dstar = (delta[2] - total_dstr * due) / station.u.max(1e-6);
            let dn3 = if total_dstr.abs() > 1e-12 { d_dstar / total_dstr } else { 0.0 };
            // XFOIL: DN4 = ABS(DUEDG)/0.25 (fixed constant)
            let dn4 = due.abs() / 0.25;
            
            sample_updates.push(rustfoil_bl::SampleUpdateData {
                ibl,
                dn1,
                dn2,
                dn3,
                dn4,
                d_ctau: delta[0],
                d_theta: delta[1],
                d_dstar,
                d_ue: due,
            });
        }
        
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::update_summary(
            iter,
            rlx_computed,
            rlx,
            limiting_station,
            limiting_reason,
            sample_updates,
        ));
    }
    
    // === Compute final RMS of normalized changes (XFOIL's RMSBL) ===
    // RMSBL = sqrt(sum(DN1² + DN2² + DN3² + DN4²) / (4*N))
    let rms_change = if n_stations > 0 {
        (rms_sum / (4.0 * n_stations as f64)).sqrt()
    } else {
        0.0
    };
    
    GlobalUpdateResult {
        relaxation_used: rlx,
        rms_change,
        max_change: max_dn,
        upper_u_new: upper_new_ue,
        lower_u_new: lower_new_ue,
        upper_u_ac: upper_operating_ue,
        lower_u_ac: lower_operating_ue,
    }
}

pub fn apply_global_updates_from_view(
    state: &mut CanonicalNewtonStateView,
    deltas: &[[f64; 3]],
    operating_deltas: Option<&[[f64; 3]]>,
    system: &GlobalNewtonSystem,
    iter: usize,
    relaxation: f64,
    dij: &nalgebra::DMatrix<f64>,
    dac: f64,
    msq: f64,
    re: f64,
) -> GlobalUpdateResult {
    let dhi: f64 = 1.5;
    let dlo: f64 = -0.5;
    let mut rlx = relaxation.clamp(0.0, 1.0);

    let preview = preview_global_update_ue_from_view(
        state,
        deltas,
        operating_deltas,
        system,
        dij,
    );
    let upper_new_ue = preview.upper_u_new;
    let lower_new_ue = preview.lower_u_new;
    let upper_operating_ue = preview.upper_u_ac;
    let lower_operating_ue = preview.lower_u_ac;

    let corrected_delta = |iv: usize| -> [f64; 3] {
        if let Some(operating) = operating_deltas.filter(|_| iv < deltas.len()) {
            let op = operating.get(iv).copied().unwrap_or([0.0; 3]);
            [
                deltas[iv][0] - dac * op[0],
                deltas[iv][1] - dac * op[1],
                deltas[iv][2] - dac * op[2],
            ]
        } else if iv < deltas.len() {
            deltas[iv]
        } else {
            [0.0; 3]
        }
    };

    let mut rms_sum = 0.0;
    let mut max_dn = 0.0_f64;
    let mut n_stations = 0usize;

    for ibl in 1..state.upper_rows.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = corrected_delta(iv);
        let row = &state.upper_rows[ibl];
        let dn1 = if row.is_laminar {
            delta[0] / 10.0
        } else if row.ctau.abs() > 1e-10 {
            delta[0] / row.ctau
        } else {
            delta[0] / 0.03
        };
        let dn2 = if row.theta.abs() > 1e-12 { delta[1] / row.theta } else { 0.0 };
        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(row.uedg)
            + dac * upper_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - row.uedg;
        let total_dstr = total_dstr_for_row_update(row);
        let ddstr = (delta[2] - total_dstr * due) / row.uedg.max(1e-12);
        let dn3 = if total_dstr.abs() > 1e-12 { ddstr / total_dstr } else { 0.0 };
        let dn4 = due.abs() / 0.25;
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
        max_dn = max_dn.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());
        n_stations += 1;
        let rdn1 = rlx * dn1;
        if rdn1 > dhi { rlx = dhi / dn1; }
        if rdn1 < dlo { rlx = dlo / dn1; }
        let rdn2 = rlx * dn2;
        if rdn2 > dhi { rlx = dhi / dn2; }
        if rdn2 < dlo { rlx = dlo / dn2; }
        let rdn3 = rlx * dn3;
        if rdn3 > dhi { rlx = dhi / dn3; }
        if rdn3 < dlo { rlx = dlo / dn3; }
        let rdn4 = rlx * dn4;
        if rdn4 > dhi { rlx = dhi / dn4; }
    }

    for ibl in 1..state.lower_rows.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = corrected_delta(iv);
        let row = &state.lower_rows[ibl];
        let dn1 = if row.is_laminar {
            delta[0] / 10.0
        } else if row.ctau.abs() > 1e-10 {
            delta[0] / row.ctau
        } else {
            delta[0] / 0.03
        };
        let dn2 = if row.theta.abs() > 1e-12 { delta[1] / row.theta } else { 0.0 };
        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(row.uedg)
            + dac * lower_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - row.uedg;
        let total_dstr = total_dstr_for_row_update(row);
        let ddstr = (delta[2] - total_dstr * due) / row.uedg.max(1e-12);
        let dn3 = if total_dstr.abs() > 1e-12 { ddstr / total_dstr } else { 0.0 };
        let dn4 = due.abs() / 0.25;
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
        max_dn = max_dn.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());
        n_stations += 1;
        let rdn1 = rlx * dn1;
        if rdn1 > dhi { rlx = dhi / dn1; }
        if rdn1 < dlo { rlx = dlo / dn1; }
        let rdn2 = rlx * dn2;
        if rdn2 > dhi { rlx = dhi / dn2; }
        if rdn2 < dlo { rlx = dlo / dn2; }
        let rdn3 = rlx * dn3;
        if rdn3 > dhi { rlx = dhi / dn3; }
        if rdn3 < dlo { rlx = dlo / dn3; }
        let rdn4 = rlx * dn4;
        if rdn4 > dhi { rlx = dhi / dn4; }
    }

    rlx = rlx.max(0.0).min(1.0);
    let rlx_computed = rlx;

    for ibl in 1..state.upper_rows.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = corrected_delta(iv);
        let flow_type = state
            .upper_flow_types
            .get(ibl.saturating_sub(1))
            .copied()
            .unwrap_or(FlowType::Turbulent);
        let row = &mut state.upper_rows[ibl];
        if row.is_laminar {
            row.ampl += rlx * delta[0];
        } else {
            row.ctau += rlx * delta[0];
            row.ctau = row.ctau.min(0.25);
        }
        row.theta += rlx * delta[1];
        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(row.uedg)
            + dac * upper_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - row.uedg;
        let delta_mass = delta[2];
        let total_dstr = total_dstr_for_row_update(row);
        let d_dstar = (delta_mass - total_dstr * due) / row.uedg.max(1e-12);
        row.uedg += rlx * due;
        let total_dstr_after = total_dstr + rlx * d_dstar;
        row.dstr = if row.is_wake { total_dstr_after - row.dw } else { total_dstr_after };
        let mut station = row.as_station(flow_type);
        station.u = row.uedg;
        station.theta = row.theta;
        station.delta_star = row.dstr;
        station.ctau = row.ctau;
        station.ampl = row.ampl;
        station.dw = row.dw;
        apply_dslim_like_xfoil(&mut station, msq);
        blvar(&mut station, flow_type, msq, re);
        station.mass_defect = station.u * total_dstr_for_update(&station);
        row.overwrite_from_station(&station);
    }

    for ibl in 1..state.lower_rows.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = corrected_delta(iv);
        let flow_type = state
            .lower_flow_types
            .get(ibl.saturating_sub(1))
            .copied()
            .unwrap_or(FlowType::Turbulent);
        let row = &mut state.lower_rows[ibl];
        if row.is_laminar {
            row.ampl += rlx * delta[0];
        } else {
            row.ctau += rlx * delta[0];
            row.ctau = row.ctau.min(0.25);
        }
        row.theta += rlx * delta[1];
        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(row.uedg)
            + dac * lower_operating_ue.get(ibl).copied().unwrap_or(0.0);
        let due = new_u - row.uedg;
        let delta_mass = delta[2];
        let total_dstr = total_dstr_for_row_update(row);
        let d_dstar = (delta_mass - total_dstr * due) / row.uedg.max(1e-12);
        row.uedg += rlx * due;
        let total_dstr_after = total_dstr + rlx * d_dstar;
        row.dstr = if row.is_wake { total_dstr_after - row.dw } else { total_dstr_after };
        let mut station = row.as_station(flow_type);
        station.u = row.uedg;
        station.theta = row.theta;
        station.delta_star = row.dstr;
        station.ctau = row.ctau;
        station.ampl = row.ampl;
        station.dw = row.dw;
        apply_dslim_like_xfoil(&mut station, msq);
        blvar(&mut station, flow_type, msq, re);
        station.mass_defect = station.u * total_dstr_for_update(&station);
        row.overwrite_from_station(&station);
    }

    state.upper_ue_current = state.upper_rows.iter().map(|row| row.uedg).collect();
    state.lower_ue_current = state.lower_rows.iter().map(|row| row.uedg).collect();

    let rms_change = if n_stations > 0 {
        (rms_sum / (4.0 * n_stations as f64)).sqrt()
    } else {
        0.0
    };

    GlobalUpdateResult {
        relaxation_used: rlx,
        rms_change,
        max_change: max_dn,
        upper_u_new: upper_new_ue,
        lower_u_new: lower_new_ue,
        upper_u_ac: upper_operating_ue,
        lower_u_ac: lower_operating_ue,
    }
}

pub fn preview_global_update_ue(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    deltas: &[[f64; 3]],
    operating_deltas: Option<&[[f64; 3]]>,
    system: &GlobalNewtonSystem,
    dij: Option<&nalgebra::DMatrix<f64>>,
    upper_ue_inv: Option<&[f64]>,
    lower_ue_inv: Option<&[f64]>,
    upper_ue_operating: Option<&[f64]>,
    lower_ue_operating: Option<&[f64]>,
) -> GlobalOperatingPreview {
    let (upper_u_new, lower_u_new) = if let (Some(dij), Some(ue_inv_u), Some(ue_inv_l)) =
        (dij, upper_ue_inv, lower_ue_inv)
    {
        compute_new_ue_via_dij(
            upper_stations,
            lower_stations,
            deltas,
            system,
            dij,
            ue_inv_u,
            ue_inv_l,
        )
    } else {
        (
            upper_stations.iter().map(|s| s.u).collect::<Vec<_>>(),
            lower_stations.iter().map(|s| s.u).collect::<Vec<_>>(),
        )
    };

    let (upper_u_ac, lower_u_ac) = if let (
        Some(dij),
        Some(operating_deltas),
        Some(ue_op_u),
        Some(ue_op_l),
    ) = (dij, operating_deltas, upper_ue_operating, lower_ue_operating)
    {
        compute_operating_ue_via_dij(
            upper_stations,
            lower_stations,
            operating_deltas,
            system,
            dij,
            ue_op_u,
            ue_op_l,
        )
    } else {
        (
            vec![0.0; upper_stations.len()],
            vec![0.0; lower_stations.len()],
        )
    };

    GlobalOperatingPreview {
        upper_u_new,
        lower_u_new,
        upper_u_ac,
        lower_u_ac,
    }
}

pub fn preview_global_update_ue_from_view(
    state: &CanonicalNewtonStateView,
    deltas: &[[f64; 3]],
    operating_deltas: Option<&[[f64; 3]]>,
    system: &GlobalNewtonSystem,
    dij: &nalgebra::DMatrix<f64>,
) -> GlobalOperatingPreview {
    let (upper_u_new, lower_u_new) = compute_new_ue_via_dij_from_rows(
        &state.upper_rows,
        &state.lower_rows,
        deltas,
        system,
        dij,
        &state.upper_ue_inviscid,
        &state.lower_ue_inviscid,
    );
    let (upper_u_ac, lower_u_ac) = if let Some(operating_deltas) = operating_deltas {
        compute_operating_ue_via_dij_from_rows(
            &state.upper_rows,
            &state.lower_rows,
            operating_deltas,
            system,
            dij,
            &state.upper_ue_operating,
            &state.lower_ue_operating,
        )
    } else {
        (
            vec![0.0; state.upper_rows.len()],
            vec![0.0; state.lower_rows.len()],
        )
    };
    GlobalOperatingPreview {
        upper_u_new,
        lower_u_new,
        upper_u_ac,
        lower_u_ac,
    }
}

fn compute_new_ue_via_dij_from_rows(
    upper_rows: &[CanonicalNewtonRow],
    lower_rows: &[CanonicalNewtonRow],
    deltas: &[[f64; 3]],
    system: &GlobalNewtonSystem,
    dij: &nalgebra::DMatrix<f64>,
    upper_ue_inv: &[f64],
    lower_ue_inv: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n_upper = upper_rows.len();
    let n_lower = lower_rows.len();
    let mut upper_new_ue = vec![0.0; n_upper];
    let mut lower_new_ue = vec![0.0; n_lower];
    let vti_upper = 1.0_f64;
    let vti_lower = -1.0_f64;
    upper_new_ue[0] = upper_rows[0].uedg;
    for i in 1..n_upper {
        let panel_i = upper_rows[i].panel_idx;
        if panel_i >= dij.nrows() {
            upper_new_ue[i] = upper_rows[i].uedg;
            continue;
        }
        let uinv_i = upper_ue_inv.get(i).copied().unwrap_or(upper_rows[i].uedg);
        let mut dui = 0.0;
        for j in 1..n_upper {
            let panel_j = upper_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = upper_rows[j].mass + delta_mass;
            let ue_m = -vti_upper * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        for j in 1..n_lower {
            let panel_j = lower_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = lower_rows[j].mass + delta_mass;
            let ue_m = -vti_upper * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        upper_new_ue[i] = uinv_i + dui;
    }
    lower_new_ue[0] = lower_rows[0].uedg;
    for i in 1..n_lower {
        let panel_i = lower_rows[i].panel_idx;
        if panel_i >= dij.nrows() {
            lower_new_ue[i] = lower_rows[i].uedg;
            continue;
        }
        let uinv_i = lower_ue_inv.get(i).copied().unwrap_or(lower_rows[i].uedg);
        let mut dui = 0.0;
        for j in 1..n_upper {
            let panel_j = upper_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = upper_rows[j].mass + delta_mass;
            let ue_m = -vti_lower * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        for j in 1..n_lower {
            let panel_j = lower_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < deltas.len() { deltas[iv][2] } else { 0.0 };
            let proposed_mass = lower_rows[j].mass + delta_mass;
            let ue_m = -vti_lower * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * proposed_mass;
        }
        lower_new_ue[i] = uinv_i + dui;
    }
    (upper_new_ue, lower_new_ue)
}

fn compute_operating_ue_via_dij_from_rows(
    upper_rows: &[CanonicalNewtonRow],
    lower_rows: &[CanonicalNewtonRow],
    operating_deltas: &[[f64; 3]],
    system: &GlobalNewtonSystem,
    dij: &nalgebra::DMatrix<f64>,
    upper_ue_operating: &[f64],
    lower_ue_operating: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n_upper = upper_rows.len();
    let n_lower = lower_rows.len();
    let mut upper_u_ac = vec![0.0; n_upper];
    let mut lower_u_ac = vec![0.0; n_lower];
    let vti_upper = 1.0_f64;
    let vti_lower = -1.0_f64;

    for i in 1..n_upper {
        let panel_i = upper_rows[i].panel_idx;
        if panel_i >= dij.nrows() {
            continue;
        }
        let mut dui = upper_ue_operating.get(i).copied().unwrap_or(0.0);
        for j in 1..n_upper {
            let panel_j = upper_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < operating_deltas.len() { operating_deltas[iv][2] } else { 0.0 };
            let ue_m = -vti_upper * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        for j in 1..n_lower {
            let panel_j = lower_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < operating_deltas.len() { operating_deltas[iv][2] } else { 0.0 };
            let ue_m = -vti_upper * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        upper_u_ac[i] = dui;
    }

    for i in 1..n_lower {
        let panel_i = lower_rows[i].panel_idx;
        if panel_i >= dij.nrows() {
            continue;
        }
        let mut dui = lower_ue_operating.get(i).copied().unwrap_or(0.0);
        for j in 1..n_upper {
            let panel_j = upper_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < operating_deltas.len() { operating_deltas[iv][2] } else { 0.0 };
            let ue_m = -vti_lower * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        for j in 1..n_lower {
            let panel_j = lower_rows[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < operating_deltas.len() { operating_deltas[iv][2] } else { 0.0 };
            let ue_m = -vti_lower * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        lower_u_ac[i] = dui;
    }

    (upper_u_ac, lower_u_ac)
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
    
    // Skip station 0 (stagnation point) - XFOIL starts from ibl=2
    // Keep original Ue at stagnation point (Ue≈0.001)
    upper_new_ue[0] = upper_stations[0].u;
    
    // Compute new Ue for upper surface stations (IS=1), starting from station 1
    for i in 1..n_upper {
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
        // XFOIL starts from JBL=2 (skips stagnation), so we start from j=1
        for j in 1..n_upper {
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
        // XFOIL starts from JBL=2 (skips stagnation), so we start from j=1
        for j in 1..n_lower {
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
        
        // Debug: trace dui computation at sample stations
        if std::env::var("RUSTFOIL_LE_UPDATE_DEBUG").is_ok()
            && (i <= 3 || i == 20 || i == 40)
        {
            eprintln!("[DEBUG compute_new_ue] upper i={}: uinv={:.6}, dui={:.6e}, ue_new={:.6}, current_u={:.6}",
                i, uinv_i, dui, ue_new, upper_stations[i].u);
        }
        
        // XFOIL's UPDATE uses the raw inviscid-plus-source result here and lets
        // STMOVE apply the tiny post-update UEPS floor near stagnation.
        upper_new_ue[i] = ue_new;
    }
    
    // Skip station 0 (stagnation point) - XFOIL starts from ibl=2
    // Keep original Ue at stagnation point (Ue≈0.001)
    lower_new_ue[0] = lower_stations[0].u;
    
    // Compute new Ue for lower surface stations (IS=2), starting from station 1
    for i in 1..n_lower {
        let panel_i = lower_stations[i].panel_idx;
        if panel_i >= dij.nrows() {
            lower_new_ue[i] = lower_stations[i].u;
            continue;
        }
        
        let uinv_i = lower_ue_inv.get(i).copied().unwrap_or(lower_stations[i].u);
        let mut dui = 0.0_f64;
        
        // Upper surface contribution (JS=1, VTI_j = +1)
        // UE_M = -VTI_lower * VTI_upper * DIJ = -(-1) * 1 * DIJ = +DIJ
        // XFOIL starts from JBL=2 (skips stagnation), so we start from j=1
        for j in 1..n_upper {
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
        // XFOIL starts from JBL=2 (skips stagnation), so we start from j=1
        for j in 1..n_lower {
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

        if std::env::var("RUSTFOIL_LE_UPDATE_DEBUG").is_ok() && i <= 3 {
            eprintln!(
                "[DEBUG compute_new_ue] lower i={}: uinv={:.6}, dui={:.6e}, ue_new={:.6}, current_u={:.6}",
                i, uinv_i, dui, ue_new, lower_stations[i].u
            );
        }
        
        // XFOIL's UPDATE uses the raw inviscid-plus-source result here and lets
        // STMOVE apply the tiny post-update UEPS floor near stagnation.
        lower_new_ue[i] = ue_new;
    }
    
    (upper_new_ue, lower_new_ue)
}

fn compute_operating_ue_via_dij(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    operating_deltas: &[[f64; 3]],
    system: &GlobalNewtonSystem,
    dij: &nalgebra::DMatrix<f64>,
    upper_ue_operating: &[f64],
    lower_ue_operating: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n_upper = upper_stations.len();
    let n_lower = lower_stations.len();
    let mut upper_u_ac = vec![0.0; n_upper];
    let mut lower_u_ac = vec![0.0; n_lower];

    upper_u_ac[0] = 0.0;
    lower_u_ac[0] = 0.0;

    let vti_upper = 1.0_f64;
    let vti_lower = -1.0_f64;

    for i in 1..n_upper {
        let panel_i = upper_stations[i].panel_idx;
        if panel_i >= dij.nrows() {
            continue;
        }
        let mut dui = upper_ue_operating.get(i).copied().unwrap_or(0.0);
        for j in 1..n_upper {
            let panel_j = upper_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < operating_deltas.len() {
                -operating_deltas[iv][2]
            } else {
                0.0
            };
            let ue_m = -vti_upper * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        for j in 1..n_lower {
            let panel_j = lower_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < operating_deltas.len() {
                -operating_deltas[iv][2]
            } else {
                0.0
            };
            let ue_m = -vti_upper * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        upper_u_ac[i] = dui;
    }

    for i in 1..n_lower {
        let panel_i = lower_stations[i].panel_idx;
        if panel_i >= dij.nrows() {
            continue;
        }
        let mut dui = lower_ue_operating.get(i).copied().unwrap_or(0.0);
        for j in 1..n_upper {
            let panel_j = upper_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(0, j);
            let delta_mass = if iv < operating_deltas.len() {
                -operating_deltas[iv][2]
            } else {
                0.0
            };
            let ue_m = -vti_lower * vti_upper * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        for j in 1..n_lower {
            let panel_j = lower_stations[j].panel_idx;
            if panel_j >= dij.ncols() {
                continue;
            }
            let iv = system.to_global(1, j);
            let delta_mass = if iv < operating_deltas.len() {
                -operating_deltas[iv][2]
            } else {
                0.0
            };
            let ue_m = -vti_lower * vti_lower * dij[(panel_i, panel_j)];
            dui += ue_m * delta_mass;
        }
        lower_u_ac[i] = dui;
    }

    (upper_u_ac, lower_u_ac)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    use rustfoil_bl::equations::FlowType;

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

    fn make_station(panel_idx: usize, x: f64, u: f64, theta: f64, delta_star: f64, ctau: f64) -> BlStation {
        let mut s = BlStation::new();
        s.panel_idx = panel_idx;
        s.x = x;
        s.u = u;
        s.theta = theta;
        s.delta_star = delta_star;
        s.ctau = ctau;
        s.mass_defect = u * delta_star;
        s.is_laminar = false;
        s.is_turbulent = true;
        s
    }

    #[test]
    fn test_te_wake_interface_uses_combined_te_state() {
        let mut upper = vec![
            make_station(0, 0.0, 0.01, 0.001, 0.002, 0.02),
            make_station(1, 0.4, 0.8, 0.002, 0.003, 0.03),
            make_station(2, 1.0, 0.9, 0.004, 0.006, 0.05),
        ];
        let mut lower = vec![
            make_station(3, 0.0, 0.01, 0.001, 0.002, 0.02),
            make_station(4, 0.4, 0.7, 0.0025, 0.0035, 0.025),
            make_station(5, 1.0, 0.85, 0.003, 0.005, 0.04),
            make_station(6, 1.1, 0.82, 0.010, 0.012, 0.01),
        ];
        lower[3].is_wake = true;
        upper[0].is_laminar = true;
        lower[0].is_laminar = true;

        let mut system = GlobalNewtonSystem::new(upper.len(), lower.len(), 2, 2);
        system.apply_te_wake_interface(&upper, &lower);
        system.build_vz_block(&upper, &lower);

        let wake_iv = system.to_global(1, 3);
        let tte = upper[2].theta + lower[2].theta;
        let dte = upper[2].delta_star + lower[2].delta_star;
        let cte = (upper[2].ctau * upper[2].theta + lower[2].ctau * lower[2].theta) / tte;

        assert!((system.vdel[wake_iv][0] - (cte - lower[3].ctau)).abs() < 1e-12);
        assert!((system.vdel[wake_iv][1] - (tte - lower[3].theta)).abs() < 1e-12);
        assert!((system.vdel[wake_iv][2] - (dte - lower[3].delta_star)).abs() < 1e-12);

        assert!((system.vz[0][0] + upper[2].theta / tte).abs() < 1e-12);
        assert!((system.vz[1][1] + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_second_rhs_survives_global_solve() {
        let mut system = GlobalNewtonSystem::new(2, 2, 1, 1);

        for iv in 1..=system.nsys {
            system.va[iv][0][0] = 1.0;
            system.va[iv][1][1] = 1.0;
            system.vm[iv][iv][2] = 1.0;
        }

        system.vdel[1] = [1.0, 2.0, 3.0];
        system.vdel[2] = [4.0, 5.0, 6.0];
        system.vdel_operating[1] = [0.5, 0.25, 0.125];
        system.vdel_operating[2] = [1.5, 1.25, 1.125];

        let result = solve_global_system(&mut system);

        assert_eq!(result.state_deltas[1], [1.0, 2.0, 3.0]);
        assert_eq!(result.state_deltas[2], [4.0, 5.0, 6.0]);
        assert_eq!(result.operating_deltas[1], [0.5, 0.25, 0.125]);
        assert_eq!(result.operating_deltas[2], [1.5, 1.25, 1.125]);
    }

    #[test]
    fn test_operating_rhs_uses_surface_alpha_sensitivity() {
        let mut upper = vec![
            make_station(0, 0.0, 0.01, 0.001, 0.002, 0.02),
            make_station(1, 0.2, 0.8, 0.002, 0.003, 0.03),
            make_station(2, 0.6, 0.85, 0.0025, 0.0036, 0.035),
        ];
        let mut lower = vec![
            make_station(3, 0.0, 0.01, 0.001, 0.002, 0.02),
            make_station(4, 0.2, 0.75, 0.0021, 0.0032, 0.03),
            make_station(5, 0.6, 0.8, 0.0023, 0.0034, 0.032),
        ];
        upper[0].is_laminar = true;
        lower[0].is_laminar = true;

        let mut system = GlobalNewtonSystem::new(upper.len(), lower.len(), 2, 2);
        system.setup_index_mappings(&upper, &lower);

        let ue_current_upper: Vec<f64> = upper.iter().map(|s| s.u).collect();
        let ue_current_lower: Vec<f64> = lower.iter().map(|s| s.u).collect();
        let ue_from_mass_upper = ue_current_upper.clone();
        let ue_from_mass_lower = ue_current_lower.clone();
        let flow_types = vec![FlowType::Turbulent; 2];
        let transitions = vec![None, None];

        system.build_surface_equations(
            &upper,
            &flow_types,
            &transitions,
            &ue_current_upper,
            &ue_from_mass_upper,
            9.0,
            0.0,
            1.0e6,
            0,
        );
        system.build_surface_equations(
            &lower,
            &flow_types,
            &transitions,
            &ue_current_lower,
            &ue_from_mass_lower,
            9.0,
            0.0,
            1.0e6,
            1,
        );
        system.set_stagnation_derivs(0.0, 0.0);
        system.build_operating_rhs(&upper, &lower, &[0.0, 0.1, 0.2], &[0.0, -0.1, -0.2]);

        let upper_iv = system.to_global(0, 1);
        let lower_iv = system.to_global(1, 1);
        assert!(system.vdel_operating[upper_iv].iter().any(|v| v.abs() > 0.0));
        assert!(system.vdel_operating[lower_iv].iter().any(|v| v.abs() > 0.0));
    }

    #[test]
    fn test_apply_global_updates_uses_operating_correction() {
        let mut upper = vec![BlStation::new(), BlStation::new()];
        let mut lower = vec![BlStation::new(), BlStation::new()];

        upper[1].u = 1.0;
        upper[1].theta = 1.0;
        upper[1].delta_star = 1.0;
        upper[1].ampl = 0.2;
        upper[1].mass_defect = 1.0;
        upper[1].is_laminar = true;
        upper[1].is_turbulent = false;
        upper[1].panel_idx = 1;

        lower[1].u = 1.0;
        lower[1].theta = 1.0;
        lower[1].delta_star = 1.0;
        lower[1].ampl = 0.1;
        lower[1].mass_defect = 1.0;
        lower[1].is_laminar = true;
        lower[1].is_turbulent = false;
        lower[1].panel_idx = 2;

        let system = GlobalNewtonSystem::new(upper.len(), lower.len(), 1, 1);
        let mut base = vec![[0.0; 3]; system.nsys + 1];
        let mut operating = vec![[0.0; 3]; system.nsys + 1];
        let upper_iv = system.to_global(0, 1);

        base[upper_iv] = [1.0, 1.0, 1.0];
        operating[upper_iv] = [0.25, 0.25, 0.25];

        let result = apply_global_updates(
            &mut upper,
            &mut lower,
            &base,
            Some(&operating),
            &system,
            1,
            1.0,
            None,
            None,
            None,
            None,
            None,
            2.0,
            0.0,
            1.0e6,
        );

        assert!((result.relaxation_used - 1.0).abs() < 1.0e-12);
        assert!((upper[1].ampl - 0.7).abs() < 1.0e-12);
        assert!((upper[1].theta - 1.5).abs() < 1.0e-12);
        assert!(upper[1].delta_star > 1.0);
        assert!((upper[1].u - 1.0).abs() < 1.0e-12);
    }
}
