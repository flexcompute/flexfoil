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
use rustfoil_bl::equations::{bldif, bldif_full_simi, blvar, FlowType};
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
    /// Dimensions: 3 equations (rows) x 2 coupling terms (columns)
    pub vz: [[f64; 2]; 3],
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
            vti: vec![1.0; nsys + 1],
            panel_idx: vec![0; nsys + 1],
            is_wake: vec![false; nsys + 1],
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
        self.sst_go = sst_go;
        self.sst_gp = sst_gp;
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
    /// * `iteration` - Newton iteration number (for debug output)
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
        iteration: usize,
    ) {
        // Store iteration for debug events
        self.current_iteration = iteration;
        // Reset arrays
        for iv in 0..=self.nsys {
            self.va[iv] = [[0.0; 3]; 3];
            self.vb[iv] = [[0.0; 3]; 3];
            self.vdel[iv] = [0.0; 3];
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
            ue_inviscid_upper,
            ue_inviscid_lower,
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
            // Note: transitions parameter is no longer used - XFOIL doesn't use TRDIF in global Newton
            let _ = transitions;  // Silence unused warning
            
            // Compute residuals and Jacobian
            // For first interval (i==1), use similarity-specific bldif with XFOIL's
            // fixed log values (XLOG=1, ULOG=1, TLOG=0, HLOG=0) and DDLOG=0
            // which zeros out the log-derivative Jacobian terms
            // 
            // IMPORTANT: XFOIL uses regular BLSYS/BLDIF for ALL stations in global Newton,
            // including the transition station. TRDIF is only used during MRCHUE (march phase)
            // for local Newton iterations. Using TRDIF in global Newton causes the shear-lag
            // residual to be 10 orders of magnitude too large!
            //
            // NOTE: For SIMI (i==1), XFOIL's march already converged the similarity equations.
            // Our march uses universal correction factors that give good approximations but
            // don't exactly satisfy bldif_full_simi. To avoid huge SIMI residuals disrupting
            // Newton convergence, we set SIMI residuals to zero (treating it as a boundary
            // condition). The Jacobians are still computed to provide proper coupling.
            let (mut residuals, jacobian) = if i == 1 {
                bldif_full_simi(s2, flow_type, msq, re)
            } else {
                bldif(&stations[i - 1], s2, flow_type, msq, re)
            };
            
            // Zero out SIMI residuals (boundary condition approach)
            // The march has already set similarity values; we just need to couple them
            if i == 1 {
                residuals.res_mom = 0.0;
                residuals.res_shape = 0.0;
                // Keep res_third (ampl) as-is since it's already -ampl ≈ 0 for laminar
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
            if rustfoil_bl::is_debug_active() && i % 10 == 0 {
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
            if i == 1 {
                // Similarity station: XFOIL combines VS2 = VS1 + VS2, VS1 = 0
                // This effectively doubles the jacobian for columns 1-4.
                // The combining is done in XFOIL's BLDIF after computing both VS matrices.
                //
                // CRITICAL FIX: We MUST add vs1+vs2 to match XFOIL's combined jacobian.
                // This is NOT a 2x error - XFOIL intentionally doubles the sensitivity
                // at the SIMI station because both "upstream" and "downstream" affect
                // the same physical location.
                for eq in 0..3 {
                    // Only columns 0,1 go in VA (ampl/ctau and theta)
                    // For SIMI, add vs1+vs2 to match XFOIL's combining
                    self.va[iv][eq][0] = jacobian.vs1[eq][0] + jacobian.vs2[eq][0];
                    self.va[iv][eq][1] = jacobian.vs1[eq][1] + jacobian.vs2[eq][1];
                    self.va[iv][eq][2] = 0.0; // Not used - mass goes in VM
                    self.vb[iv][eq][0] = 0.0;
                    self.vb[iv][eq][1] = 0.0;
                    self.vb[iv][eq][2] = 0.0;
                    
                    // Store delta_star and Ue columns for VM construction
                    // These also need combining: vs2_combined = vs1 + vs2
                    self.vs1_delta[iv][eq] = 0.0;
                    self.vs1_ue[iv][eq] = 0.0;
                    self.vs2_delta[iv][eq] = jacobian.vs1[eq][2] + jacobian.vs2[eq][2];
                    self.vs2_ue[iv][eq] = jacobian.vs1[eq][3] + jacobian.vs2[eq][3];
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
                    
                    // Store X derivatives (column 4, 0-indexed) for stagnation coupling
                    // XFOIL: VS1(k,5), VS2(k,5) - derivatives w.r.t. arc length
                    self.vs1_x[iv][eq] = jacobian.vs1[eq][4];
                    self.vs2_x[iv][eq] = jacobian.vs2[eq][4];
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
                    -self.vti[iv - 1] * self.vti[jv] * dij[(panel_im1, panel_j)]
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

            // Emit forced changes debug event at key stations (every 10th)
            if rustfoil_bl::is_debug_active() && ibl % 10 == 0 {
                let ue_from_mass_i = ue_from_mass_upper.get(ibl).copied().unwrap_or(0.0);
                let ue_inviscid_i = ue_inviscid_upper.get(ibl).copied().unwrap_or(0.0);
                let ue_current_i = ue_current_upper.get(ibl).copied().unwrap_or(0.0);
                
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::forced_changes(
                    self.current_iteration,
                    ibl,
                    0,  // side = upper
                    due1,
                    due2,
                    dds1,
                    dds2,
                    ue_from_mass_i,
                    ue_inviscid_i,
                    ue_current_i,
                ));
            }

            // Debug: print forced change contributions at station 2 (first iteration only)
            if std::env::var("RUSTFOIL_FC_DEBUG").is_ok() && ibl == 2 {
                eprintln!("[FC DEBUG] upper[2]: due1={:.6e}, due2={:.6e}, dds1={:.6e}, dds2={:.6e}",
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

            // Emit forced changes debug event at key stations (every 10th)
            if rustfoil_bl::is_debug_active() && ibl % 10 == 0 {
                let ue_from_mass_i = ue_from_mass_lower.get(ibl).copied().unwrap_or(0.0);
                let ue_inviscid_i = ue_inviscid_lower.get(ibl).copied().unwrap_or(0.0);
                let ue_current_i = ue_current_lower.get(ibl).copied().unwrap_or(0.0);
                
                rustfoil_bl::add_event(rustfoil_bl::DebugEvent::forced_changes(
                    self.current_iteration,
                    ibl,
                    1,  // side = lower
                    due1,
                    due2,
                    dds1,
                    dds2,
                    ue_from_mass_i,
                    ue_inviscid_i,
                    ue_current_i,
                ));
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
    system: &GlobalNewtonSystem,
    relaxation: f64,
    dij: Option<&nalgebra::DMatrix<f64>>,
    upper_ue_inv: Option<&[f64]>,
    lower_ue_inv: Option<&[f64]>,
    msq: f64,
    re: f64,
) -> GlobalUpdateResult {
    // === XFOIL-style normalized relaxation (xbl.f lines 1527-1597) ===
    // Compute global relaxation factor based on normalized changes
    let dhi = 1.5;   // Max 150% relative increase
    let dlo = -0.5;  // Max 50% relative decrease (CRITICAL for ctau stability)
    
    let mut rlx = relaxation.clamp(0.0, 1.0);
    
    // CRITICAL FIX: We need to compute new Ue FIRST to check Ue changes,
    // but we'll do a preliminary check on mass changes using current Ue
    // to get an initial relaxation estimate, then refine it after computing new Ue.
    
    // First pass: compute required relaxation from theta and ctau/ampl changes
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
        // CRITICAL FIX: Check absolute value of normalized change to catch both positive and negative large changes
        let rdn1 = rlx * dn1;
        if rdn1.abs() > dhi && dn1.abs() > 1e-12 {
            let new_rlx = dhi / dn1.abs();
            if new_rlx < rlx {
                rlx = new_rlx;
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG apply_global_updates] Reduced rlx to {:.6e} due to ctau/ampl change at lower[{}]: dn1={:.6e}, delta[0]={:.6e}, ctau={:.6e}",
                        rlx, ibl, dn1, delta[0], station.ctau);
                }
            }
        }
        if rdn1.abs() > dlo.abs() && dn1.abs() > 1e-12 && dn1 < 0.0 {
            let new_rlx = dlo.abs() / dn1.abs();
            if new_rlx < rlx {
                rlx = new_rlx;
                if std::env::var("RUSTFOIL_CL_DEBUG").is_ok() {
                    eprintln!("[DEBUG apply_global_updates] Reduced rlx to {:.6e} due to negative ctau/ampl change at lower[{}]: dn1={:.6e}",
                        rlx, ibl, dn1);
                }
            }
        }
        
        let rdn2 = rlx * dn2;
        if rdn2.abs() > dhi && dn2.abs() > 1e-12 {
            let new_rlx = dhi / dn2.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rdn2.abs() > dlo.abs() && dn2.abs() > 1e-12 && dn2 < 0.0 {
            let new_rlx = dlo.abs() / dn2.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
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
    
    // === Step 1.5: Refine relaxation based on Ue and delta_star changes ===
    // CRITICAL FIX: Check normalized Ue and delta_star changes to prevent explosion
    // Upper surface
    for ibl in 1..upper_stations.len() {
        let iv = system.to_global(0, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = &deltas[iv];
        let station = &upper_stations[ibl];
        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u);
        let due = new_u - station.u;
        
        // Normalized Ue change - XFOIL uses fixed constant 0.25, NOT current Ue!
        // xbl.f line 1572: DN4 = ABS(DUEDG)/0.25
        // This prevents huge normalized changes at stations near stagnation (small Ue)
        let dn4 = due.abs() / 0.25;
        
        // Normalized delta_star change
        // Compute d_dstar using current Ue (will be refined after relaxation)
        let d_dstar = (delta[2] - station.delta_star * due) / station.u.max(1e-6);
        let dn3 = if station.delta_star.abs() > 1e-12 {
            d_dstar / station.delta_star
        } else {
            d_dstar / 1e-6  // Fallback for very small delta_star
        };
        
        // Apply relaxation limits to Ue changes
        let rdn4 = rlx * dn4;
        let rlx_before_ue = rlx;
        if rdn4.abs() > dhi && dn4.abs() > 1e-12 {
            let new_rlx = dhi / dn4.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rdn4.abs() > dlo.abs() && dn4.abs() > 1e-12 && dn4 < 0.0 {
            let new_rlx = dlo.abs() / dn4.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rlx < rlx_before_ue && rustfoil_bl::is_debug_active() {
            eprintln!("[DEBUG apply_global_updates] Reduced rlx from {:.6e} to {:.6e} due to Ue change at upper[{}]: dn4={:.6e}",
                rlx_before_ue, rlx, ibl, dn4);
        }
        
        // Apply relaxation limits to delta_star changes
        let rdn3 = rlx * dn3;
        let rlx_before_dstar = rlx;
        if rdn3.abs() > dhi && dn3.abs() > 1e-12 {
            let new_rlx = dhi / dn3.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rdn3.abs() > dlo.abs() && dn3.abs() > 1e-12 && dn3 < 0.0 {
            let new_rlx = dlo.abs() / dn3.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rlx < rlx_before_dstar && rustfoil_bl::is_debug_active() {
            eprintln!("[DEBUG apply_global_updates] Reduced rlx from {:.6e} to {:.6e} due to delta_star change at upper[{}]: dn3={:.6e}",
                rlx_before_dstar, rlx, ibl, dn3);
        }
    }
    
    // Lower surface
    for ibl in 1..lower_stations.len() {
        let iv = system.to_global(1, ibl);
        if iv >= deltas.len() {
            continue;
        }
        let delta = &deltas[iv];
        let station = &lower_stations[ibl];
        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(station.u);
        let due = new_u - station.u;
        
        // Normalized Ue change - XFOIL uses fixed constant 0.25
        // xbl.f line 1572: DN4 = ABS(DUEDG)/0.25
        let dn4 = due.abs() / 0.25;
        
        // Normalized delta_star change
        let d_dstar = (delta[2] - station.delta_star * due) / station.u.max(1e-6);
        let dn3 = if station.delta_star.abs() > 1e-12 {
            d_dstar / station.delta_star
        } else {
            d_dstar / 1e-6
        };
        
        // Apply relaxation limits
        let rdn4 = rlx * dn4;
        let rlx_before_ue = rlx;
        if rdn4.abs() > dhi && dn4.abs() > 1e-12 {
            let new_rlx = dhi / dn4.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rdn4.abs() > dlo.abs() && dn4.abs() > 1e-12 && dn4 < 0.0 {
            let new_rlx = dlo.abs() / dn4.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rlx < rlx_before_ue && rustfoil_bl::is_debug_active() {
            eprintln!("[DEBUG apply_global_updates] Reduced rlx from {:.6e} to {:.6e} due to Ue change at lower[{}]: dn4={:.6e}",
                rlx_before_ue, rlx, ibl, dn4);
        }
        
        let rdn3 = rlx * dn3;
        let rlx_before_dstar = rlx;
        if rdn3.abs() > dhi && dn3.abs() > 1e-12 {
            let new_rlx = dhi / dn3.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rdn3.abs() > dlo.abs() && dn3.abs() > 1e-12 && dn3 < 0.0 {
            let new_rlx = dlo.abs() / dn3.abs();
            if new_rlx < rlx { rlx = new_rlx; }
        }
        if rlx < rlx_before_dstar && rustfoil_bl::is_debug_active() {
            eprintln!("[DEBUG apply_global_updates] Reduced rlx from {:.6e} to {:.6e} due to delta_star change at lower[{}]: dn3={:.6e}",
                rlx_before_dstar, rlx, ibl, dn3);
        }
    }
    
    // Ensure relaxation stays positive and reasonable
    rlx = rlx.max(0.01).min(1.0);
    
    // Store computed relaxation before applying
    let rlx_computed = relaxation;
    
    // === Accumulate RMS of normalized changes (XFOIL's RMSBL) ===
    // RMSBL = sqrt(sum(DN1² + DN2² + DN3² + DN4²) / (4*N))
    // where DN1-DN4 are normalized changes in ctau/ampl, theta, delta*, and Ue
    let mut rms_sum = 0.0;
    let mut max_dn = 0.0_f64;
    let mut n_stations = 0usize;
    
    // Collect sample update data for debug output
    let mut sample_updates: Vec<rustfoil_bl::SampleUpdateData> = Vec::new();
    let mut limiting_station: Option<usize> = None;
    let mut limiting_reason: Option<String> = None;
    
    // Track which station limited relaxation most
    if rlx < relaxation * 0.99 {
        // Find the limiting station by checking normalized changes
        for ibl in 1..upper_stations.len() {
            let iv = system.to_global(0, ibl);
            if iv >= deltas.len() { continue; }
            let station = &upper_stations[ibl];
            let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u);
            let due = new_u - station.u;
            
            let dn1 = if station.is_laminar {
                deltas[iv][0] / 10.0
            } else if station.ctau.abs() > 1e-10 {
                deltas[iv][0] / station.ctau
            } else {
                deltas[iv][0] / 0.03
            };
            let dn2 = if station.theta.abs() > 1e-12 { deltas[iv][1] / station.theta } else { 0.0 };
            // XFOIL: DN4 = ABS(DUEDG)/0.25 (fixed constant, not relative to Ue)
            let dn4 = due.abs() / 0.25;
            
            if (rlx * dn1).abs() >= dhi || (rlx * dn2).abs() >= dhi || (rlx * dn4).abs() >= dhi {
                limiting_station = Some(ibl);
                if (rlx * dn1).abs() >= dhi {
                    limiting_reason = Some(format!("upper ctau dn1={:.3}", dn1));
                } else if (rlx * dn2).abs() >= dhi {
                    limiting_reason = Some(format!("upper theta dn2={:.3}", dn2));
                } else {
                    limiting_reason = Some(format!("upper Ue dn4={:.3}", dn4));
                }
                break;
            }
        }
    }

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
        
        // Safeguard: limit theta changes to prevent blow-up or collapse
        // - Don't let theta decrease by more than 50% per iteration
        // - Don't let theta increase by more than 2x per iteration (matches direct march)
        // - Absolute maximum: 0.1 (physically reasonable for TE stations)
        // CRITICAL FIX: Without upper bound, Newton can produce huge theta values
        // (e.g., theta=0.398 at α=-4° causing CD=0.59 instead of 0.0008)
        let theta_min = theta_old * 0.5;
        let theta_max = theta_old * 2.0;  // Max 2x increase per iteration
        let theta_abs_max = 0.1;  // Absolute maximum (matches march.rs line 543)
        station.theta = station.theta
            .max(theta_min.max(1e-8))
            .min(theta_max.min(theta_abs_max));

        // delta[2] is mass change; convert to delta_star using XFOIL formula:
        // DDSTR = (DMASS - DSTR*DUEDG) / UEDG
        // where DUEDG = UNEW - UEDG is the *proposed* Ue change
        let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u);
        let due = new_u - station.u;
        
        // Debug: trace Ue updates at sample stations
        if rustfoil_bl::is_debug_active() && (ibl == 20 || ibl == 40 || ibl == 60) {
            eprintln!("[DEBUG Ue UPDATE] upper ibl={}: new_u={:.6}, current_u={:.6}, due={:.6e}, delta_mass={:.6e}, rlx={:.3}",
                ibl, new_u, station.u, due, delta[2], rlx);
        }
        
        // XFOIL formula: delta_dstar = (delta_mass - dstar * delta_ue) / ue
        let delta_mass = delta[2];
        let d_dstar = (delta_mass - station.delta_star * due) / station.u.max(1e-6);
        
        // Apply Ue update with relaxation (XFOIL: UEDG = UEDG + RLX*DUEDG)
        let u_old = station.u;
        station.u += rlx * due;
        // Safeguard: Ue must remain positive for attached flow
        station.u = station.u.max(0.01);
        
        // Debug: show Ue change
        if rustfoil_bl::is_debug_active() && (ibl == 20 || ibl == 40 || ibl == 60) {
            eprintln!("[DEBUG Ue UPDATE] upper ibl={}: Ue changed {:.6} -> {:.6} (delta={:.6e})",
                ibl, u_old, station.u, station.u - u_old);
        }
        
        // Apply delta_star update (XFOIL: DSTR = DSTR + RLX*DDSTR)
        station.delta_star += rlx * d_dstar;
        // Ensure H is physical: H >= 1.0 (displacement thickness >= momentum thickness)
        station.delta_star = station.delta_star.max(station.theta);

        station.h = (station.delta_star / station.theta).clamp(1.0, 20.0);
        station.mass_defect = station.u * station.delta_star;
        
        // Recompute secondary variables (Hk, Cf, Cd, etc.) for next iteration
        let flow_type = if station.is_laminar { FlowType::Laminar } else { FlowType::Turbulent };
        blvar(station, flow_type, msq, re);
        
        // === Accumulate normalized changes for RMSBL (XFOIL xbl.f lines 1545-1572) ===
        // CRITICAL: RMSBL uses PROPOSED changes (no rlx), not APPLIED changes
        // XFOIL computes DN1-DN4 from raw deltas, then applies relaxation separately
        // DN1 = normalized ctau/ampl change
        // DN2 = normalized theta change  
        // DN3 = normalized delta_star change
        // DN4 = normalized Ue change (XFOIL uses fixed 0.25 constant)
        let dn1 = if upper_stations[ibl].is_laminar {
            delta[0] / 10.0  // Laminar: fixed scale (NO rlx!)
        } else if upper_stations[ibl].ctau.abs() > 1e-10 {
            delta[0] / upper_stations[ibl].ctau
        } else {
            delta[0] / 0.03
        };
        
        let dn2 = if upper_stations[ibl].theta.abs() > 1e-12 {
            delta[1] / upper_stations[ibl].theta
        } else {
            0.0
        };
        
        // For DN3 and DN4, use the PROPOSED changes (no rlx)
        let dn3 = if upper_stations[ibl].delta_star.abs() > 1e-12 {
            d_dstar / upper_stations[ibl].delta_star
        } else {
            d_dstar / 1e-6
        };
        
        // XFOIL: DN4 = ABS(DUEDG)/0.25 (fixed constant, not relative to Ue)
        let dn4 = due.abs() / 0.25;
        
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
        max_dn = max_dn.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());
        n_stations += 1;
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
        
        // Safeguard: limit theta changes to prevent blow-up or collapse
        // - Don't let theta decrease by more than 50% per iteration
        // - Don't let theta increase by more than 2x per iteration (matches direct march)
        // - Absolute maximum: 0.1 (physically reasonable for TE stations)
        // CRITICAL FIX: Without upper bound, Newton can produce huge theta values
        // (e.g., theta=0.398 at α=-4° causing CD=0.59 instead of 0.0008)
        let theta_min = theta_old * 0.5;
        let theta_max = theta_old * 2.0;  // Max 2x increase per iteration
        let theta_abs_max = 0.1;  // Absolute maximum (matches march.rs line 543)
        station.theta = station.theta
            .max(theta_min.max(1e-8))
            .min(theta_max.min(theta_abs_max));

        // delta[2] is mass change; convert to delta_star using XFOIL formula:
        // DDSTR = (DMASS - DSTR*DUEDG) / UEDG
        let new_u = lower_new_ue.get(ibl).copied().unwrap_or(station.u);
        let due = new_u - station.u;
        
        // Debug: trace Ue updates at sample stations
        if rustfoil_bl::is_debug_active() && (ibl == 20 || ibl == 40) {
            eprintln!("[DEBUG Ue UPDATE] lower ibl={}: new_u={:.6}, current_u={:.6}, due={:.6e}, delta_mass={:.6e}, rlx={:.3}",
                ibl, new_u, station.u, due, delta[2], rlx);
        }
        
        // XFOIL formula: delta_dstar = (delta_mass - dstar * delta_ue) / ue
        let delta_mass = delta[2];
        let d_dstar = (delta_mass - station.delta_star * due) / station.u.max(1e-6);
        
        // Apply Ue update with relaxation (XFOIL: UEDG = UEDG + RLX*DUEDG)
        let u_old = station.u;
        station.u += rlx * due;
        // Safeguard: Ue must remain positive for attached flow
        station.u = station.u.max(0.01);
        
        // Debug: show Ue change
        if rustfoil_bl::is_debug_active() && (ibl == 20 || ibl == 40) {
            eprintln!("[DEBUG Ue UPDATE] lower ibl={}: Ue changed {:.6} -> {:.6} (delta={:.6e})",
                ibl, u_old, station.u, station.u - u_old);
        }

        // Apply delta_star update (XFOIL: DSTR = DSTR + RLX*DDSTR)
        station.delta_star += rlx * d_dstar;
        // Ensure H is physical: H >= 1.0 (displacement thickness >= momentum thickness)
        station.delta_star = station.delta_star.max(station.theta);

        station.h = (station.delta_star / station.theta).clamp(1.0, 20.0);
        station.mass_defect = station.u * station.delta_star;
        
        // Recompute secondary variables (Hk, Cf, Cd, etc.) for next iteration
        let flow_type = if station.is_laminar { FlowType::Laminar } else { FlowType::Turbulent };
        blvar(station, flow_type, msq, re);
        
        // === Accumulate normalized changes for RMSBL (XFOIL xbl.f lines 1545-1572) ===
        // CRITICAL: RMSBL uses PROPOSED changes (no rlx), not APPLIED changes
        // XFOIL: IF(IBL.GE.ITRAN(IS)) DN1 = DCTAU / CTAU(IBL,IS)
        let dn1 = if lower_stations[ibl].is_laminar {
            delta[0] / 10.0  // NO rlx!
        } else if lower_stations[ibl].ctau.abs() > 1e-10 {
            delta[0] / lower_stations[ibl].ctau
        } else {
            delta[0] / 0.03
        };
        
        let dn2 = if lower_stations[ibl].theta.abs() > 1e-12 {
            delta[1] / lower_stations[ibl].theta
        } else {
            0.0
        };
        
        // Use PROPOSED changes (no rlx) for DN3 and DN4
        let dn3 = if lower_stations[ibl].delta_star.abs() > 1e-12 {
            d_dstar / lower_stations[ibl].delta_star
        } else {
            d_dstar / 1e-6
        };
        
        // XFOIL: DN4 = ABS(DUEDG)/0.25 (fixed constant)
        let dn4 = due.abs() / 0.25;
        
        // Debug transition stations (60-70 on lower surface)
        if ibl >= 60 && ibl <= 70 && std::env::var("RUSTFOIL_TRANS_DEBUG").is_ok() {
            // Find first turbulent station index for reference
            let trans_idx = lower_stations.iter()
                .position(|s| !s.is_laminar)
                .unwrap_or(9999);
            let is_trans = ibl == trans_idx;
            let normalizer = if lower_stations[ibl].is_laminar {
                "10.0".to_string()
            } else if lower_stations[ibl].ctau.abs() > 1e-10 {
                format!("ctau={:.6e}", lower_stations[ibl].ctau)
            } else {
                "0.03".to_string()
            };
            eprintln!("[TRANS DEBUG] lower[{}]: is_laminar={}, ctau={:.6e}, theta={:.6e}, dstr={:.6e}, Ue={:.4}{}",
                ibl, lower_stations[ibl].is_laminar, lower_stations[ibl].ctau, 
                lower_stations[ibl].theta, lower_stations[ibl].delta_star, lower_stations[ibl].u,
                if is_trans { " <-- TRANSITION" } else { "" });
            eprintln!("              delta[0]={:.6e}, delta[1]={:.6e}, delta[2]={:.6e}",
                delta[0], delta[1], delta[2]);
            eprintln!("              dn1={:.4} (div by {}), dn2={:.4}, dn3={:.4}, dn4={:.4}",
                dn1, normalizer, dn2, dn3, dn4);
            // Show the raw calculation for transition station
            if is_trans {
                eprintln!("              [TRANS] dn1 = delta[0]/ctau = {:.6e}/{:.6e} = {:.4}",
                    delta[0], lower_stations[ibl].ctau, delta[0] / lower_stations[ibl].ctau.max(1e-10));
            }
        }
        
        rms_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
        max_dn = max_dn.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());
        n_stations += 1;
    }
    
    // Emit UPDATE_SUMMARY debug event with relaxation details
    if rustfoil_bl::is_debug_active() {
        // Collect sample update data for every 10th station (upper surface)
        // Skip station 0 (stagnation point) - start from ibl=1
        for ibl in (1..upper_stations.len()).filter(|i| i % 10 == 0) {
            let iv = system.to_global(0, ibl);
            if iv >= deltas.len() { continue; }
            
            let station = &upper_stations[ibl];
            let new_u = upper_new_ue.get(ibl).copied().unwrap_or(station.u);
            let due = new_u - station.u;
            
            let dn1 = if station.is_laminar {
                deltas[iv][0] / 10.0
            } else if station.ctau.abs() > 1e-10 {
                deltas[iv][0] / station.ctau
            } else {
                deltas[iv][0] / 0.03
            };
            let dn2 = if station.theta.abs() > 1e-12 { deltas[iv][1] / station.theta } else { 0.0 };
            let d_dstar = (deltas[iv][2] - station.delta_star * due) / station.u.max(1e-6);
            let dn3 = if station.delta_star.abs() > 1e-12 { d_dstar / station.delta_star } else { 0.0 };
            // XFOIL: DN4 = ABS(DUEDG)/0.25 (fixed constant)
            let dn4 = due.abs() / 0.25;
            
            sample_updates.push(rustfoil_bl::SampleUpdateData {
                ibl,
                dn1,
                dn2,
                dn3,
                dn4,
                d_ctau: deltas[iv][0],
                d_theta: deltas[iv][1],
                d_dstar,
                d_ue: due,
            });
        }
        
        rustfoil_bl::add_event(rustfoil_bl::DebugEvent::update_summary(
            0,  // iteration will be set by caller context
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
        if rustfoil_bl::is_debug_active() && (i == 20 || i == 40) {
            eprintln!("[DEBUG compute_new_ue] upper i={}: uinv={:.6}, dui={:.6e}, ue_new={:.6}, current_u={:.6}",
                i, uinv_i, dui, ue_new, upper_stations[i].u);
        }
        
        // XFOIL's UPDATE: UNEW is computed directly from mass defect, no artificial clamping.
        // The relaxation in apply_global_updates handles stability.
        // Only clamp to physical bounds (Ue must be positive for attached flow).
        upper_new_ue[i] = ue_new.clamp(0.01, 5.0);
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
        
        // XFOIL's UPDATE: UNEW is computed directly from mass defect, no artificial clamping.
        // The relaxation in apply_global_updates handles stability.
        // Only clamp to physical bounds (Ue must be positive for attached flow).
        lower_new_ue[i] = ue_new.clamp(0.01, 5.0);
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
