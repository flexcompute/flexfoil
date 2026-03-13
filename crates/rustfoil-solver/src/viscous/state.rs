use rustfoil_bl::equations::FlowType;
use rustfoil_bl::state::{BlDerivatives, BlStation};
use rustfoil_coupling::march::{march_mixed_du, march_surface, MarchConfig, MarchResult};
use rustfoil_coupling::newton_state::CanonicalNewtonStateView;
use rustfoil_coupling::stmove::{
    adjust_transition_for_stmove, find_stagnation_with_derivs, StagnationResult,
};

use super::config::OperatingMode;
use super::setup::extract_surface_xfoil;

/// XFOIL surface selector using XFOIL's `IS` numbering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XfoilSurface {
    Upper,
    Lower,
}

impl XfoilSurface {
    pub const fn from_xfoil_is(is: usize) -> Option<Self> {
        match is {
            1 => Some(Self::Upper),
            2 => Some(Self::Lower),
            _ => None,
        }
    }

    pub const fn xfoil_is(self) -> usize {
        match self {
            Self::Upper => 1,
            Self::Lower => 2,
        }
    }

    pub const fn vti(self) -> f64 {
        match self {
            Self::Upper => 1.0,
            Self::Lower => -1.0,
        }
    }
}

/// XFOIL-style `(IBL, IS)` index.
///
/// `ibl` uses XFOIL's 1-based station numbering, while `surface` carries the
/// corresponding `IS` selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XfoilBlIndex {
    pub ibl: usize,
    pub surface: XfoilSurface,
}

impl XfoilBlIndex {
    pub fn from_xfoil(ibl: usize, is: usize) -> Option<Self> {
        if ibl == 0 {
            return None;
        }

        Some(Self {
            ibl,
            surface: XfoilSurface::from_xfoil_is(is)?,
        })
    }
}

fn flow_type_for_station_index(ibl: usize, iblte: usize, itran: Option<usize>) -> FlowType {
    let itran = itran.unwrap_or(iblte + 1);
    if ibl > iblte {
        FlowType::Wake
    } else if ibl + 1 >= itran {
        FlowType::Turbulent
    } else {
        FlowType::Laminar
    }
}

fn flow_type_for_row(
    row: &CanonicalBlRow,
    ibl: usize,
    iblte: usize,
    itran: Option<usize>,
) -> FlowType {
    if row.is_wake {
        FlowType::Wake
    } else {
        flow_type_for_station_index(ibl, iblte, itran)
    }
}

/// Canonical per-station row owned by the transitional XFOIL-style state.
///
/// This is intentionally close to the existing `BlStation` shape so the current
/// solver can move data in and out with minimal disruption during Phase 1.
#[derive(Debug, Clone, Default)]
pub struct CanonicalBlRow {
    pub x: f64,
    pub x_coord: f64,
    pub y_coord: f64,
    pub panel_idx: usize,
    pub uedg: f64,
    pub uinv: f64,
    pub uinv_a: f64,
    pub theta: f64,
    pub dstr: f64,
    pub ctau_or_ampl: f64,
    pub ctau: f64,
    pub ampl: f64,
    pub mass: f64,
    pub h: f64,
    pub hk: f64,
    pub hs: f64,
    pub hc: f64,
    pub r_theta: f64,
    pub cf: f64,
    pub cd: f64,
    pub us: f64,
    pub cq: f64,
    pub de: f64,
    pub dw: f64,
    pub is_laminar: bool,
    pub is_turbulent: bool,
    pub is_wake: bool,
    pub derivs: BlDerivatives,
}

impl CanonicalBlRow {
    pub fn from_station(station: &BlStation) -> Self {
        let mut row = Self {
            x: station.x,
            x_coord: station.x_coord,
            y_coord: 0.0,
            panel_idx: station.panel_idx,
            uedg: station.u,
            uinv: station.u,
            uinv_a: 0.0,
            theta: station.theta,
            dstr: station.delta_star,
            ctau_or_ampl: 0.0,
            ctau: station.ctau,
            ampl: station.ampl,
            mass: station.mass_defect,
            h: station.h,
            hk: station.hk,
            hs: station.hs,
            hc: station.hc,
            r_theta: station.r_theta,
            cf: station.cf,
            cd: station.cd,
            us: station.us,
            cq: station.cq,
            de: station.de,
            dw: station.dw,
            is_laminar: station.is_laminar,
            is_turbulent: station.is_turbulent,
            is_wake: station.is_wake,
            derivs: station.derivs.clone(),
        };
        row.sync_ctau_or_ampl();
        row
    }

    pub fn as_station_with_flow_type(&self, flow_type: FlowType) -> BlStation {
        let is_wake = matches!(flow_type, FlowType::Wake);
        let is_turbulent = matches!(flow_type, FlowType::Turbulent | FlowType::Wake);
        let is_laminar = matches!(flow_type, FlowType::Laminar);
        BlStation {
            x: self.x,
            x_coord: self.x_coord,
            panel_idx: self.panel_idx,
            u: self.uedg,
            theta: self.theta,
            delta_star: self.dstr,
            ctau: self.ctau,
            ampl: self.ampl,
            h: self.h,
            hk: self.hk,
            hs: self.hs,
            hc: self.hc,
            r_theta: self.r_theta,
            cf: self.cf,
            cd: self.cd,
            us: self.us,
            cq: self.cq,
            de: self.de,
            mass_defect: self.mass,
            dw: self.dw,
            is_laminar,
            is_wake,
            is_turbulent,
            derivs: self.derivs.clone(),
        }
    }

    pub fn overwrite_from_station(&mut self, station: &BlStation) {
        *self = Self::from_station(station);
    }

    fn sync_ctau_or_ampl(&mut self) {
        self.ctau_or_ampl = if self.is_laminar && !self.is_turbulent {
            self.ampl
        } else {
            self.ctau
        };
    }
}

/// Canonical XFOIL-style viscous state owner.
///
/// Phase 1 keeps this state additive and transitional: existing solver paths can
/// still operate on `Vec<BlStation>` while this object accumulates the eventual
/// panel/BL/wake arrays and metadata under one owner.
#[derive(Debug, Clone, Default)]
pub struct XfoilLikeViscousState {
    pub ist: usize,
    pub sst: f64,
    pub sst_go: f64,
    pub sst_gp: f64,
    pub nbl_upper: usize,
    pub nbl_lower: usize,
    pub iblte_upper: usize,
    pub iblte_lower: usize,
    pub ipan_upper: Vec<usize>,
    pub ipan_lower: Vec<usize>,
    pub vti_upper: Vec<f64>,
    pub vti_lower: Vec<f64>,
    pub isys_upper: Vec<Option<usize>>,
    pub isys_lower: Vec<Option<usize>>,
    pub qinv: Vec<f64>,
    pub qinv_a: Vec<f64>,
    pub uinv: Vec<f64>,
    pub qvis: Vec<f64>,
    pub gam: Vec<f64>,
    pub gam_a: Vec<f64>,
    pub upper_rows: Vec<CanonicalBlRow>,
    pub lower_rows: Vec<CanonicalBlRow>,
    pub wake_x: Vec<f64>,
    pub wake_y: Vec<f64>,
    pub wake_s: Vec<f64>,
    pub wake_panel_indices: Vec<usize>,
    pub wake_uedg: Vec<f64>,
    pub wake_mass: Vec<f64>,
    pub itran_upper: Option<usize>,
    pub itran_lower: Option<usize>,
    pub xtran_upper: Option<f64>,
    pub xtran_lower: Option<f64>,
    pub operating_mode: OperatingMode,
    pub u_new: Vec<f64>,
    pub u_ac: Vec<f64>,
    pub dac: f64,
    pub rlx: f64,
}

#[derive(Debug, Clone)]
pub struct TransitionalStmoveResult {
    pub ist: usize,
}

impl XfoilLikeViscousState {
    pub fn new(n_panel_nodes: usize) -> Self {
        Self {
            qinv: vec![0.0; n_panel_nodes],
            qinv_a: vec![0.0; n_panel_nodes],
            uinv: vec![0.0; n_panel_nodes],
            qvis: vec![0.0; n_panel_nodes],
            gam: vec![0.0; n_panel_nodes],
            gam_a: vec![0.0; n_panel_nodes],
            u_new: vec![0.0; n_panel_nodes],
            u_ac: vec![0.0; n_panel_nodes],
            rlx: 1.0,
            ..Default::default()
        }
    }

    pub fn from_station_views(
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        n_panel_nodes: usize,
    ) -> Self {
        let mut state = Self::new(n_panel_nodes);
        state.sync_from_station_views(upper_stations, lower_stations, n_panel_nodes);
        state
    }

    pub fn set_stagnation_metadata(&mut self, ist: usize, sst: f64, sst_go: f64, sst_gp: f64) {
        self.ist = ist;
        self.sst = sst;
        self.sst_go = sst_go;
        self.sst_gp = sst_gp;
    }

    pub fn set_transition_metadata(
        &mut self,
        surface: XfoilSurface,
        transition_index: Option<usize>,
        x_transition: Option<f64>,
    ) {
        match surface {
            XfoilSurface::Upper => {
                self.itran_upper = transition_index;
                self.xtran_upper = x_transition;
            }
            XfoilSurface::Lower => {
                self.itran_lower = transition_index;
                self.xtran_lower = x_transition;
            }
        }
    }

    pub fn update_transition_metadata(
        &mut self,
        surface: XfoilSurface,
        transition_index: Option<usize>,
        x_transition: Option<f64>,
    ) {
        match surface {
            XfoilSurface::Upper => {
                if let Some(idx) = transition_index {
                    self.itran_upper = Some(idx);
                }
                if let Some(x) = x_transition {
                    self.xtran_upper = Some(x);
                }
            }
            XfoilSurface::Lower => {
                if let Some(idx) = transition_index {
                    self.itran_lower = Some(idx);
                }
                if let Some(x) = x_transition {
                    self.xtran_lower = Some(x);
                }
            }
        }
    }

    pub fn adjust_transition_metadata_for_stmove(&mut self, old_ist: usize, new_ist: usize) {
        self.itran_upper = self
            .itran_upper
            .map(|itran| adjust_transition_for_stmove(itran, old_ist, new_ist, true));
        self.itran_lower = self
            .itran_lower
            .map(|itran| adjust_transition_for_stmove(itran, old_ist, new_ist, false));
    }

    pub fn set_panel_inviscid_arrays(&mut self, qinv: &[f64], qinv_a: &[f64]) {
        resize_if_needed(&mut self.qinv, qinv.len());
        resize_if_needed(&mut self.qinv_a, qinv_a.len());
        resize_if_needed(&mut self.uinv, qinv.len());
        self.qinv.clone_from_slice(qinv);
        self.qinv_a.clone_from_slice(qinv_a);
        self.uinv.clone_from_slice(qinv);
        self.sync_row_inviscid_from_panel_arrays();
    }

    pub fn set_operating_mode(&mut self, operating_mode: OperatingMode) {
        self.operating_mode = operating_mode;
    }

    pub fn set_operating_scratch(
        &mut self,
        upper_u_new: &[f64],
        lower_u_new: &[f64],
        upper_u_ac: &[f64],
        lower_u_ac: &[f64],
        dac: f64,
        rlx: f64,
    ) {
        let total = upper_u_new.len() + lower_u_new.len();
        resize_if_needed(&mut self.u_new, total);
        resize_if_needed(&mut self.u_ac, total);

        self.u_new.clear();
        self.u_new.extend_from_slice(upper_u_new);
        self.u_new.extend_from_slice(lower_u_new);

        self.u_ac.clear();
        self.u_ac.extend_from_slice(upper_u_ac);
        self.u_ac.extend_from_slice(lower_u_ac);

        self.dac = dac;
        self.rlx = rlx;
    }

    pub fn set_wake_geometry(&mut self, wake_x: &[f64], wake_y: &[f64], wake_s: &[f64]) {
        self.wake_x = wake_x.to_vec();
        self.wake_y = wake_y.to_vec();
        self.wake_s = wake_s.to_vec();
    }

    pub fn n_panel_nodes(&self) -> usize {
        self.qvis
            .len()
            .max(self.gam.len())
            .max(self.qinv.len())
            .max(self.qinv_a.len())
    }

    pub fn nbl(&self, surface: XfoilSurface) -> usize {
        match surface {
            XfoilSurface::Upper => self.nbl_upper,
            XfoilSurface::Lower => self.nbl_lower,
        }
    }

    pub fn iblte(&self, surface: XfoilSurface) -> usize {
        match surface {
            XfoilSurface::Upper => self.iblte_upper,
            XfoilSurface::Lower => self.iblte_lower,
        }
    }

    pub fn xfoil_index(&self, ibl: usize, is: usize) -> Option<(XfoilSurface, usize)> {
        let index = XfoilBlIndex::from_xfoil(ibl, is)?;
        let row_idx = index.ibl - 1;
        let len = self.nbl(index.surface);
        (row_idx < len).then_some((index.surface, row_idx))
    }

    pub fn row(&self, ibl: usize, is: usize) -> Option<&CanonicalBlRow> {
        let (surface, row_idx) = self.xfoil_index(ibl, is)?;
        self.surface_rows(surface).get(row_idx)
    }

    pub fn row_mut(&mut self, ibl: usize, is: usize) -> Option<&mut CanonicalBlRow> {
        let (surface, row_idx) = self.xfoil_index(ibl, is)?;
        self.surface_rows_mut(surface).get_mut(row_idx)
    }

    /// Derive an upper-surface `BlStation` view from canonical rows.
    pub fn upper_station_view(&self) -> Vec<BlStation> {
        self.station_view(XfoilSurface::Upper)
    }

    /// Derive a lower-surface `BlStation` view from canonical rows.
    pub fn lower_station_view(&self) -> Vec<BlStation> {
        self.station_view(XfoilSurface::Lower)
    }

    /// Derive a `BlStation` view for one surface.
    pub fn station_view(&self, surface: XfoilSurface) -> Vec<BlStation> {
        let iblte = self.iblte(surface);
        let itran = match surface {
            XfoilSurface::Upper => self.itran_upper,
            XfoilSurface::Lower => self.itran_lower,
        };
        self.surface_rows(surface)
            .iter()
            .enumerate()
            .map(|(ibl, row)| row.as_station_with_flow_type(flow_type_for_row(row, ibl, iblte, itran)))
            .collect()
    }

    pub fn flow_type_view(&self, surface: XfoilSurface) -> Vec<FlowType> {
        let iblte = self.iblte(surface);
        let itran = match surface {
            XfoilSurface::Upper => self.itran_upper,
            XfoilSurface::Lower => self.itran_lower,
        };
        self.surface_rows(surface)
            .iter()
            .enumerate()
            .skip(1)
            .map(|(ibl, row)| flow_type_for_row(row, ibl, iblte, itran))
            .collect()
    }

    /// Derive per-surface `Ue` values from canonical rows.
    pub fn uedg_view(&self, surface: XfoilSurface) -> Vec<f64> {
        self.surface_rows(surface)
            .iter()
            .map(|row| row.uedg)
            .collect()
    }

    pub fn uinv_view(&self, surface: XfoilSurface) -> Vec<f64> {
        self.surface_rows(surface)
            .iter()
            .map(|row| row.uinv)
            .collect()
    }

    pub fn operating_sensitivity_view(&self, surface: XfoilSurface) -> Vec<f64> {
        self.surface_rows(surface)
            .iter()
            .map(|row| row.uinv_a)
            .collect()
    }

    pub fn upper_uedg_view(&self) -> Vec<f64> {
        self.uedg_view(XfoilSurface::Upper)
    }

    pub fn lower_uedg_view(&self) -> Vec<f64> {
        self.uedg_view(XfoilSurface::Lower)
    }

    /// Run the direct march using canonical rows as the source of truth, with a
    /// temporary `BlStation` wrapper only at the coupling API boundary.
    pub fn march_surface(
        &mut self,
        surface: XfoilSurface,
        re: f64,
        msq: f64,
        config: &MarchConfig,
    ) -> MarchResult {
        let stations = self.station_view(surface);
        let x: Vec<f64> = stations.iter().map(|station| station.x).collect();
        let ue: Vec<f64> = stations.iter().map(|station| station.u).collect();
        let result = march_surface(&x, &ue, re, msq, config, surface.xfoil_is());
        self.apply_march_result_to_surface(surface, &result);
        self.set_transition_metadata(
            surface,
            result.transition_index,
            result.x_transition,
        );
        self.refresh_panel_arrays_from_rows();
        result
    }

    /// Run MRCHDU-style mixed marching with the canonical state owning the
    /// accepted station rows.
    pub fn march_mixed_du(
        &mut self,
        surface: XfoilSurface,
        re: f64,
        msq: f64,
        config: &MarchConfig,
    ) -> MarchResult {
        let mut stations = self.station_view(surface);
        let result = march_mixed_du(&mut stations, re, msq, config, surface.xfoil_is());
        self.overwrite_surface_rows_from_stations(surface, &stations);
        self.update_transition_metadata(
            surface,
            result.transition_index,
            result.x_transition,
        );
        self.refresh_panel_arrays_from_rows();
        result
    }

    /// Overwrite the canonical rows from split-station solver views.
    pub fn sync_from_station_views(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
        n_panel_nodes: usize,
    ) {
        self.upper_rows = upper_stations.iter().map(CanonicalBlRow::from_station).collect();
        self.lower_rows = lower_stations.iter().map(CanonicalBlRow::from_station).collect();
        self.sync_row_inviscid_from_panel_arrays();
        self.sync_surface_metadata();
        self.sync_panel_arrays(n_panel_nodes);
        self.sync_wake_arrays();
    }

    /// Write the canonical rows back out into split-station representation.
    pub fn write_back_station_views(
        &self,
        upper_stations: &mut Vec<BlStation>,
        lower_stations: &mut Vec<BlStation>,
    ) {
        *upper_stations = self.upper_station_view();
        *lower_stations = self.lower_station_view();
    }

    pub fn overwrite_from_newton_view(&mut self, view: &CanonicalNewtonStateView) {
        self.upper_rows = view
            .upper_rows
            .iter()
            .map(|row| CanonicalBlRow {
                x: row.x,
                x_coord: row.x_coord,
                y_coord: 0.0,
                panel_idx: row.panel_idx,
                uedg: row.uedg,
                uinv: row.uinv,
                uinv_a: row.uinv_a,
                theta: row.theta,
                dstr: row.dstr,
                ctau_or_ampl: row.ctau_or_ampl,
                ctau: row.ctau,
                ampl: row.ampl,
                mass: row.mass,
                h: row.h,
                hk: row.hk,
                hs: row.hs,
                hc: row.hc,
                r_theta: row.r_theta,
                cf: row.cf,
                cd: row.cd,
                us: row.us,
                cq: row.cq,
                de: row.de,
                dw: row.dw,
                is_laminar: row.is_laminar,
                is_turbulent: row.is_turbulent,
                is_wake: row.is_wake,
                derivs: row.derivs.clone(),
            })
            .collect();
        self.lower_rows = view
            .lower_rows
            .iter()
            .map(|row| CanonicalBlRow {
                x: row.x,
                x_coord: row.x_coord,
                y_coord: 0.0,
                panel_idx: row.panel_idx,
                uedg: row.uedg,
                uinv: row.uinv,
                uinv_a: row.uinv_a,
                theta: row.theta,
                dstr: row.dstr,
                ctau_or_ampl: row.ctau_or_ampl,
                ctau: row.ctau,
                ampl: row.ampl,
                mass: row.mass,
                h: row.h,
                hk: row.hk,
                hs: row.hs,
                hc: row.hc,
                r_theta: row.r_theta,
                cf: row.cf,
                cd: row.cd,
                us: row.us,
                cq: row.cq,
                de: row.de,
                dw: row.dw,
                is_laminar: row.is_laminar,
                is_turbulent: row.is_turbulent,
                is_wake: row.is_wake,
                derivs: row.derivs.clone(),
            })
            .collect();
        self.refresh_panel_arrays_from_rows();
    }

    /// Phase 3 canonical owner path: refresh `QVIS/GAM/GAM_A` from the currently
    /// owned BL rows while preserving the panel-space inviscid arrays.
    pub fn refresh_panel_arrays_from_rows(&mut self) {
        let n_panel_nodes = self.n_panel_nodes();
        self.sync_row_inviscid_from_panel_arrays();
        self.sync_surface_metadata();
        self.sync_panel_arrays(n_panel_nodes);
        self.sync_wake_arrays();
    }

    pub fn finalize_viscal_state(
        &mut self,
        ue_inviscid_full: &[f64],
        full_arc: &[f64],
        panel_x: &[f64],
        panel_y: &[f64],
    ) {
        let final_stagnation = find_stagnation_with_derivs(ue_inviscid_full, full_arc);
        if let Some(stag) = final_stagnation {
            self.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
            self.xtran_upper = transition_arc_to_chord_on_surface(
                stag,
                full_arc,
                panel_x,
                panel_y,
                ue_inviscid_full,
                true,
                self.xtran_upper,
            );
            self.xtran_lower = transition_arc_to_chord_on_surface(
                stag,
                full_arc,
                panel_x,
                panel_y,
                ue_inviscid_full,
                false,
                self.xtran_lower,
            );
        }
        self.refresh_panel_arrays_from_rows();
    }

    pub fn panel_qvis(&self) -> &[f64] {
        &self.qvis
    }

    pub fn panel_gamma(&self) -> &[f64] {
        &self.gam
    }

    pub fn panel_gamma_alpha(&self) -> &[f64] {
        &self.gam_a
    }

    /// Mutate the canonical state in-place like XFOIL's `STMOVE`, while still
    /// returning split inviscid views for the not-yet-migrated Newton assembly
    /// code.
    pub fn apply_stmove_like_xfoil(
        &mut self,
        ue_inviscid_full: &[f64],
        panel_x: &[f64],
        panel_y: &[f64],
        full_arc: &[f64],
        current_gamma: &[f64],
        re: f64,
        old_ist: usize,
    ) -> Option<TransitionalStmoveResult> {
        let stag = find_stagnation_with_derivs(current_gamma, full_arc)?;
        let ue_stag = interpolate_stagnation_velocity(ue_inviscid_full, full_arc, stag.ist, stag.sst);

        let upper_geom = build_surface_geometry(
            stag.ist,
            stag.sst,
            ue_stag,
            full_arc,
            panel_x,
            panel_y,
            ue_inviscid_full,
            true,
            re,
        );
        let lower_geom = build_surface_geometry(
            stag.ist,
            stag.sst,
            ue_stag,
            full_arc,
            panel_x,
            panel_y,
            ue_inviscid_full,
            false,
            re,
        );

        let old_lower_wake_start = self
            .lower_rows
            .iter()
            .position(|row| row.is_wake)
            .unwrap_or(self.lower_rows.len());
        let old_lower_airfoil_te_x = self
            .lower_rows
            .get(old_lower_wake_start.saturating_sub(1))
            .map(|row| row.x)
            .or_else(|| self.lower_rows.last().map(|row| row.x))
            .unwrap_or(0.0);
        let new_lower_airfoil_te_x = lower_geom
            .rows
            .last()
            .map(|row| row.x)
            .unwrap_or(old_lower_airfoil_te_x);
        let wake_x_shift = new_lower_airfoil_te_x - old_lower_airfoil_te_x;

        let old_upper = self.upper_rows.clone();
        let old_lower = self.lower_rows.clone();
        let old_lower_wake = old_lower
            .get(old_lower_wake_start..)
            .unwrap_or(&[])
            .to_vec();

        if stag.ist == old_ist {
            apply_geometry_in_place(&mut self.upper_rows, &upper_geom.rows);
            apply_geometry_in_place_with_shifted_wake(
                &mut self.lower_rows,
                &lower_geom.rows,
                wake_x_shift,
            );
            self.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
            self.refresh_panel_arrays_from_rows();
            return Some(TransitionalStmoveResult {
                ist: stag.ist,
            });
        }

        let idif = stag.ist.abs_diff(old_ist);
        let mut new_upper_rows = upper_geom.rows;
        let mut new_lower_rows = lower_geom.rows;
        for wake_row in &old_lower_wake {
            let mut shifted = wake_row.clone();
            shifted.x += wake_x_shift;
            shifted.is_wake = true;
            new_lower_rows.push(shifted);
        }

        if stag.ist > old_ist {
            for i in (1 + idif)..new_upper_rows.len() {
                let src = i - idif;
                if src < old_upper.len() {
                    copy_row_state(&mut new_upper_rows[i], &old_upper[src]);
                }
            }
            if idif + 1 < new_upper_rows.len() {
                let ref_idx = idif + 1;
                let ref_row = new_upper_rows[ref_idx].clone();
                let dudx = ref_row.uedg / ref_row.x.max(1.0e-12);
                for i in (1..=idif).rev() {
                    copy_row_state(&mut new_upper_rows[i], &ref_row);
                    new_upper_rows[i].uedg = (dudx * new_upper_rows[i].x).max(1.0e-7);
                    new_upper_rows[i].mass = new_upper_rows[i].uedg * new_upper_rows[i].dstr;
                }
            }

            for (i, dst) in new_lower_rows.iter_mut().enumerate().skip(1) {
                let src = i + idif;
                if src < old_lower.len() {
                    copy_row_state(dst, &old_lower[src]);
                }
            }
        } else {
            for i in (1 + idif)..new_lower_rows.len() {
                let src = i - idif;
                if src < old_lower.len() {
                    copy_row_state(&mut new_lower_rows[i], &old_lower[src]);
                }
            }
            if idif + 1 < new_lower_rows.len() {
                let ref_idx = idif + 1;
                let ref_row = new_lower_rows[ref_idx].clone();
                let dudx = ref_row.uedg / ref_row.x.max(1.0e-12);
                for i in (1..=idif).rev() {
                    copy_row_state(&mut new_lower_rows[i], &ref_row);
                    new_lower_rows[i].uedg = (dudx * new_lower_rows[i].x).max(1.0e-7);
                    new_lower_rows[i].mass = new_lower_rows[i].uedg * new_lower_rows[i].dstr;
                }
            }

            for (i, dst) in new_upper_rows.iter_mut().enumerate().skip(1) {
                let src = i + idif;
                if src < old_upper.len() {
                    copy_row_state(dst, &old_upper[src]);
                }
            }
        }

        clamp_minimum_uedg(&mut new_upper_rows);
        clamp_minimum_uedg(&mut new_lower_rows);

        self.upper_rows = new_upper_rows;
        self.lower_rows = new_lower_rows;
        self.set_stagnation_metadata(stag.ist, stag.sst, stag.sst_go, stag.sst_gp);
        self.refresh_panel_arrays_from_rows();
        Some(TransitionalStmoveResult {
            ist: stag.ist,
        })
    }

    fn surface_rows(&self, surface: XfoilSurface) -> &[CanonicalBlRow] {
        match surface {
            XfoilSurface::Upper => &self.upper_rows,
            XfoilSurface::Lower => &self.lower_rows,
        }
    }

    fn surface_rows_mut(&mut self, surface: XfoilSurface) -> &mut Vec<CanonicalBlRow> {
        match surface {
            XfoilSurface::Upper => &mut self.upper_rows,
            XfoilSurface::Lower => &mut self.lower_rows,
        }
    }

    fn overwrite_surface_rows_from_stations(
        &mut self,
        surface: XfoilSurface,
        stations: &[BlStation],
    ) {
        let rows = self.surface_rows_mut(surface);
        if rows.len() != stations.len() {
            *rows = stations.iter().map(CanonicalBlRow::from_station).collect();
            return;
        }
        for (row, station) in rows.iter_mut().zip(stations.iter()) {
            row.overwrite_from_station(station);
        }
    }

    fn apply_march_result_to_surface(&mut self, surface: XfoilSurface, result: &MarchResult) {
        let rows = self.surface_rows_mut(surface);
        if rows.is_empty() || result.stations.is_empty() {
            return;
        }

        let offset = if rows.first().map_or(false, |row| row.x < 1e-6) {
            1
        } else {
            0
        };

        for (i, station) in result.stations.iter().enumerate() {
            let target = i + offset;
            if target < rows.len() {
                update_row_from_march_station(&mut rows[target], station);
            }
        }
    }

    fn sync_surface_metadata(&mut self) {
        self.nbl_upper = self.upper_rows.len();
        self.nbl_lower = self.lower_rows.len();
        self.iblte_upper = trailing_edge_row_index(&self.upper_rows);
        self.iblte_lower = trailing_edge_row_index(&self.lower_rows);
        self.ipan_upper = self.upper_rows.iter().map(|row| row.panel_idx).collect();
        self.ipan_lower = self.lower_rows.iter().map(|row| row.panel_idx).collect();
        self.vti_upper = vec![XfoilSurface::Upper.vti(); self.upper_rows.len()];
        self.vti_lower = vec![XfoilSurface::Lower.vti(); self.lower_rows.len()];
        self.isys_upper = build_isys_placeholder(&self.upper_rows);
        self.isys_lower = build_isys_placeholder(&self.lower_rows);
    }

    fn sync_panel_arrays(&mut self, n_panel_nodes: usize) {
        resize_if_needed(&mut self.qinv, n_panel_nodes);
        resize_if_needed(&mut self.qinv_a, n_panel_nodes);
        resize_if_needed(&mut self.uinv, n_panel_nodes);
        resize_if_needed(&mut self.qvis, n_panel_nodes);
        resize_if_needed(&mut self.gam, n_panel_nodes);
        resize_if_needed(&mut self.gam_a, n_panel_nodes);
        resize_if_needed(&mut self.u_new, n_panel_nodes);
        resize_if_needed(&mut self.u_ac, n_panel_nodes);

        self.qvis.fill(0.0);
        self.gam.fill(0.0);
        self.gam_a.clone_from(&self.qinv_a);

        for (surface, rows) in [
            (XfoilSurface::Upper, &self.upper_rows),
            (XfoilSurface::Lower, &self.lower_rows),
        ] {
            let vti = surface.vti();
            for row in rows {
                if row.is_wake || row.panel_idx >= n_panel_nodes || row.panel_idx == usize::MAX {
                    continue;
                }
                self.qvis[row.panel_idx] = vti * row.uedg;
                self.gam[row.panel_idx] = self.qvis[row.panel_idx];
            }
        }

        if let Some(last) = self.gam.last_mut() {
            *last = 0.0;
        }
    }

    fn sync_wake_arrays(&mut self) {
        self.wake_x.clear();
        self.wake_y.clear();
        self.wake_s.clear();
        self.wake_panel_indices.clear();
        self.wake_uedg.clear();
        self.wake_mass.clear();

        for row in self.lower_rows.iter().filter(|row| row.is_wake) {
            self.wake_x.push(row.x_coord);
            self.wake_y.push(row.y_coord);
            self.wake_s.push(row.x);
            self.wake_panel_indices.push(row.panel_idx);
            self.wake_uedg.push(row.uedg);
            self.wake_mass.push(row.mass);
        }
    }

    fn sync_row_inviscid_from_panel_arrays(&mut self) {
        for row in &mut self.upper_rows {
            if row.panel_idx == usize::MAX {
                row.uinv = 0.0;
                row.uinv_a = 0.0;
            } else if row.is_wake {
                row.uinv = row.uedg;
                row.uinv_a = 0.0;
            } else {
                row.uinv = self.qinv.get(row.panel_idx).copied().unwrap_or(row.uedg).abs();
                row.uinv_a = self.qinv_a.get(row.panel_idx).copied().unwrap_or(0.0);
            }
        }
        for row in &mut self.lower_rows {
            if row.panel_idx == usize::MAX {
                row.uinv = 0.0;
                row.uinv_a = 0.0;
            } else if row.is_wake {
                row.uinv = row.uedg;
                row.uinv_a = 0.0;
            } else {
                row.uinv = self.qinv.get(row.panel_idx).copied().unwrap_or(row.uedg).abs();
                row.uinv_a = -self.qinv_a.get(row.panel_idx).copied().unwrap_or(0.0);
            }
        }
    }
}

fn resize_if_needed(values: &mut Vec<f64>, size: usize) {
    if values.len() != size {
        values.resize(size, 0.0);
    }
}

fn update_row_from_march_station(row: &mut CanonicalBlRow, station: &BlStation) {
    row.uedg = station.u;
    row.uinv = station.u;
    row.theta = station.theta;
    row.dstr = station.delta_star;
    row.ctau = station.ctau;
    row.ampl = station.ampl;
    row.mass = station.mass_defect;
    row.h = station.h;
    row.hk = station.hk;
    row.hs = station.hs;
    row.hc = station.hc;
    row.r_theta = station.r_theta;
    row.cf = station.cf;
    row.cd = station.cd;
    row.us = station.us;
    row.cq = station.cq;
    row.de = station.de;
    row.is_laminar = station.is_laminar;
    row.is_turbulent = station.is_turbulent;
    row.derivs = station.derivs.clone();
    row.sync_ctau_or_ampl();
}

#[derive(Debug, Clone)]
struct SurfaceGeometrySeed {
    rows: Vec<CanonicalBlRow>,
}

fn build_surface_geometry(
    ist: usize,
    sst: f64,
    ue_stag: f64,
    full_arc: &[f64],
    panel_x: &[f64],
    panel_y: &[f64],
    ue_source: &[f64],
    is_upper: bool,
    re: f64,
) -> SurfaceGeometrySeed {
    let (arc, x, y, ue_inv) = extract_surface_xfoil(
        ist,
        sst,
        ue_stag,
        full_arc,
        panel_x,
        panel_y,
        ue_source,
        is_upper,
    );

    let n = arc.len();
    let mut rows = Vec::with_capacity(n);
    let seed_x = arc.get(1).copied().unwrap_or(1.0e-6).abs().max(1.0e-12);
    let seed_ue = ue_inv.get(1).copied().unwrap_or(0.01).abs().max(0.01);
    let mut stag_row = CanonicalBlRow::from_station(&BlStation::stagnation(seed_x, seed_ue, re));
    stag_row.x = arc[0];
    stag_row.x_coord = x[0];
    stag_row.y_coord = y[0];
    stag_row.panel_idx = usize::MAX;
    stag_row.uedg = 0.0;
    stag_row.uinv = 0.0;
    stag_row.mass = 0.0;
    stag_row.is_wake = false;
    rows.push(stag_row);

    for i in 1..n {
        let mut row = CanonicalBlRow::default();
        row.x = arc[i];
        row.x_coord = x[i];
        row.y_coord = y[i];
        row.panel_idx = surface_panel_idx(ist, panel_x.len(), is_upper, i);
        row.uedg = ue_inv[i].abs();
        row.uinv = ue_inv[i].abs();
        row.theta = 0.001;
        row.dstr = 0.002;
        row.mass = row.uedg * row.dstr;
        row.h = 2.0;
        row.hk = 2.0;
        row.hs = 1.5;
        row.cf = 0.003;
        row.cd = 0.001;
        row.us = 0.5;
        row.cq = 0.03;
        row.de = 0.006;
        row.is_laminar = true;
        row.is_turbulent = false;
        row.is_wake = false;
        row.r_theta = re * row.uedg * row.theta;
        rows.push(row);
    }

    SurfaceGeometrySeed { rows }
}

fn surface_panel_idx(
    stagnation_idx: usize,
    n_airfoil_panels: usize,
    is_upper: bool,
    station_idx: usize,
) -> usize {
    if is_upper {
        let panel = stagnation_idx as i64 - (station_idx as i64 - 1);
        if panel >= 0 {
            panel as usize
        } else {
            n_airfoil_panels + (-panel - 1) as usize
        }
    } else {
        stagnation_idx + station_idx
    }
}

fn copy_row_state(dst: &mut CanonicalBlRow, src: &CanonicalBlRow) {
    dst.theta = src.theta;
    dst.dstr = src.dstr;
    dst.ctau_or_ampl = src.ctau_or_ampl;
    dst.ctau = src.ctau;
    dst.ampl = src.ampl;
    dst.uedg = src.uedg;
    dst.uinv = src.uinv;
    dst.uinv_a = src.uinv_a;
    dst.h = src.h;
    dst.hk = src.hk;
    dst.hs = src.hs;
    dst.hc = src.hc;
    dst.r_theta = src.r_theta;
    dst.cf = src.cf;
    dst.cd = src.cd;
    dst.us = src.us;
    dst.cq = src.cq;
    dst.de = src.de;
    dst.mass = src.mass;
    dst.dw = src.dw;
    dst.is_laminar = src.is_laminar;
    dst.is_turbulent = src.is_turbulent;
    dst.is_wake = src.is_wake;
    dst.derivs = src.derivs.clone();
}

fn apply_geometry_in_place(dst: &mut [CanonicalBlRow], geom_rows: &[CanonicalBlRow]) {
    for (row, geom) in dst.iter_mut().zip(geom_rows.iter()) {
        row.x = geom.x;
        row.x_coord = geom.x_coord;
        row.y_coord = geom.y_coord;
        row.panel_idx = geom.panel_idx;
        row.is_wake = geom.is_wake;
    }
}

fn apply_geometry_in_place_with_shifted_wake(
    dst: &mut [CanonicalBlRow],
    airfoil_geom_rows: &[CanonicalBlRow],
    wake_x_shift: f64,
) {
    apply_geometry_in_place(dst, airfoil_geom_rows);
    for row in dst.iter_mut().skip(airfoil_geom_rows.len()) {
        row.x += wake_x_shift;
        row.is_wake = true;
    }
}

fn clamp_minimum_uedg(rows: &mut [CanonicalBlRow]) {
    for row in rows.iter_mut().skip(1) {
        if row.uedg <= 1.0e-7 {
            row.uedg = 1.0e-7;
            row.mass = row.dstr * row.uedg;
        }
    }
}

fn interpolate_stagnation_velocity(
    ue_source: &[f64],
    full_arc: &[f64],
    ist: usize,
    sst: f64,
) -> f64 {
    if ist + 1 < ue_source.len() && ist + 1 < full_arc.len() {
        let ds = full_arc[ist + 1] - full_arc[ist];
        if ds.abs() > 1.0e-12 {
            let frac = ((sst - full_arc[ist]) / ds).clamp(0.0, 1.0);
            ue_source[ist] + frac * (ue_source[ist + 1] - ue_source[ist])
        } else {
            ue_source[ist]
        }
    } else {
        ue_source.get(ist).copied().unwrap_or(0.0)
    }
}

fn transition_arc_to_chord_on_surface(
    stagnation: StagnationResult,
    full_arc: &[f64],
    panel_x: &[f64],
    panel_y: &[f64],
    ue_inviscid_full: &[f64],
    is_upper: bool,
    x_transition: Option<f64>,
) -> Option<f64> {
    let xt = x_transition?;
    if full_arc.is_empty() || panel_x.is_empty() || panel_y.is_empty() || ue_inviscid_full.is_empty() {
        return Some(xt);
    }

    let ist_next = (stagnation.ist + 1).min(ue_inviscid_full.len().saturating_sub(1));
    let ue_stag = if ist_next > stagnation.ist && full_arc[ist_next] != full_arc[stagnation.ist] {
        let frac = (stagnation.sst - full_arc[stagnation.ist])
            / (full_arc[ist_next] - full_arc[stagnation.ist]);
        ue_inviscid_full[stagnation.ist]
            + frac * (ue_inviscid_full[ist_next] - ue_inviscid_full[stagnation.ist])
    } else {
        ue_inviscid_full[stagnation.ist]
    };

    let (surface_arc, surface_x, _, _) = extract_surface_xfoil(
        stagnation.ist,
        stagnation.sst,
        ue_stag,
        full_arc,
        panel_x,
        panel_y,
        ue_inviscid_full,
        is_upper,
    );
    if surface_arc.is_empty() || surface_x.is_empty() {
        return Some(xt);
    }
    if xt <= surface_arc[0] {
        return Some(surface_x[0]);
    }

    for i in 1..surface_arc.len() {
        if xt <= surface_arc[i] {
            let ds = surface_arc[i] - surface_arc[i - 1];
            if ds.abs() <= 1.0e-12 {
                return Some(surface_x[i]);
            }
            let frac = ((xt - surface_arc[i - 1]) / ds).clamp(0.0, 1.0);
            return Some(surface_x[i - 1] + frac * (surface_x[i] - surface_x[i - 1]));
        }
    }

    surface_x.last().copied()
}

fn trailing_edge_row_index(rows: &[CanonicalBlRow]) -> usize {
    rows.iter()
        .enumerate()
        .filter(|(_, row)| !row.is_wake)
        .map(|(idx, _)| idx + 1)
        .last()
        .unwrap_or(rows.len())
}

fn build_isys_placeholder(rows: &[CanonicalBlRow]) -> Vec<Option<usize>> {
    let mut next = 0usize;
    rows.iter()
        .map(|row| {
            if row.panel_idx == usize::MAX || row.is_wake {
                None
            } else {
                let current = Some(next);
                next += 1;
                current
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_station(panel_idx: usize, x: f64, u: f64) -> BlStation {
        let mut station = BlStation::new();
        station.panel_idx = panel_idx;
        station.x = x;
        station.x_coord = x;
        station.u = u;
        station.theta = 0.01 + x;
        station.delta_star = 0.02 + x;
        station.refresh_mass_defect();
        station
    }

    #[test]
    fn test_xfoil_like_state_indexes_by_xfoil_coordinates() {
        let upper = vec![sample_station(usize::MAX, 0.0, 0.0), sample_station(7, 0.1, 0.4)];
        let lower = vec![
            sample_station(usize::MAX, 0.0, 0.0),
            sample_station(8, 0.1, 0.5),
            sample_station(9, 0.2, 0.6),
        ];

        let state = XfoilLikeViscousState::from_station_views(&upper, &lower, 12);

        assert_eq!(state.row(1, 1).unwrap().panel_idx, usize::MAX);
        assert_eq!(state.row(2, 1).unwrap().panel_idx, 7);
        assert_eq!(state.row(3, 2).unwrap().panel_idx, 9);
        assert!(state.row(0, 1).is_none());
        assert!(state.row(4, 1).is_none());
        assert!(state.row(1, 3).is_none());
    }

    #[test]
    fn test_xfoil_like_state_round_trips_station_views() {
        let upper = vec![sample_station(usize::MAX, 0.0, 0.0), sample_station(5, 0.1, 0.7)];
        let mut lower = vec![sample_station(usize::MAX, 0.0, 0.0), sample_station(6, 0.1, 0.6)];
        lower[1].is_wake = true;
        lower[1].is_laminar = false;
        lower[1].is_turbulent = true;

        let state = XfoilLikeViscousState::from_station_views(&upper, &lower, 10);
        let upper_view = state.upper_station_view();
        let lower_view = state.lower_station_view();

        assert_eq!(upper_view.len(), upper.len());
        assert_eq!(lower_view.len(), lower.len());
        assert!((upper_view[1].u - upper[1].u).abs() < 1.0e-12);
        assert_eq!(lower_view[1].panel_idx, lower[1].panel_idx);
        assert!(lower_view[1].is_wake);
        assert_eq!(state.wake_panel_indices, vec![6]);
    }

    #[test]
    fn test_xfoil_like_state_canonical_panel_arrays_follow_xfoil_signs() {
        let upper = vec![sample_station(usize::MAX, 0.0, 0.0), sample_station(3, 0.1, 0.8)];
        let lower = vec![sample_station(usize::MAX, 0.0, 0.0), sample_station(4, 0.1, 0.6)];

        let mut state = XfoilLikeViscousState::from_station_views(&upper, &lower, 6);
        state.set_panel_inviscid_arrays(&[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], &[1.0; 6]);
        state.refresh_panel_arrays_from_rows();

        assert!((state.panel_qvis()[3] - 0.8).abs() < 1.0e-12);
        assert!((state.panel_qvis()[4] + 0.6).abs() < 1.0e-12);
        assert!((state.panel_gamma()[3] - 0.8).abs() < 1.0e-12);
        assert!((state.panel_gamma()[4] + 0.6).abs() < 1.0e-12);
        assert_eq!(state.panel_gamma()[5], 0.0);
        assert_eq!(state.panel_gamma_alpha()[2], 1.0);
    }
}
