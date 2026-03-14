use nalgebra::DMatrix;
use rustfoil_bl::state::BlDerivatives;

use crate::config::OperatingMode;

pub type VBlock = [[f64; 2]; 3];
pub type VDelBlock = [[f64; 2]; 3];
pub type VmRow = [f64; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XfoilSurface {
    Upper,
    Lower,
}

impl XfoilSurface {
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

#[derive(Debug, Clone, Default)]
pub struct XfoilBlRow {
    pub x: f64,
    pub x_coord: f64,
    pub y_coord: f64,
    pub panel_idx: usize,
    pub uedg: f64,
    pub uinv: f64,
    pub uinv_a: f64,
    pub theta: f64,
    pub dstr: f64,
    pub mass: f64,
    pub ctau: f64,
    pub ampl: f64,
    pub hk: f64,
    pub h: f64,
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

#[derive(Debug, Clone)]
pub struct XfoilState {
    pub name: String,
    pub alpha_rad: f64,
    pub operating_mode: OperatingMode,
    pub panel_x: Vec<f64>,
    pub panel_y: Vec<f64>,
    pub panel_s: Vec<f64>,
    pub panel_xp: Vec<f64>,
    pub panel_yp: Vec<f64>,
    pub xle: f64,
    pub yle: f64,
    pub xte: f64,
    pub yte: f64,
    pub chord: f64,
    pub sharp: bool,
    pub ante: f64,
    pub qinvu_0: Vec<f64>,
    pub qinvu_90: Vec<f64>,
    pub gamu_0: Vec<f64>,
    pub gamu_90: Vec<f64>,
    pub qinv: Vec<f64>,
    pub qinv_a: Vec<f64>,
    pub qvis: Vec<f64>,
    pub gam: Vec<f64>,
    pub gam_a: Vec<f64>,
    pub wake_x: Vec<f64>,
    pub wake_y: Vec<f64>,
    pub wake_s: Vec<f64>,
    pub wake_nx: Vec<f64>,
    pub wake_ny: Vec<f64>,
    pub wake_apanel: Vec<f64>,
    pub wake_qinvu_0: Vec<f64>,
    pub wake_qinvu_90: Vec<f64>,
    pub wake_qinv: Vec<f64>,
    pub wake_qinv_a: Vec<f64>,
    pub wgap: Vec<f64>,
    pub lwake: bool,
    pub ladij: bool,
    pub lwdij: bool,
    pub awake: f64,
    pub dij: DMatrix<f64>,
    pub ist: usize,
    pub sst: f64,
    pub sst_go: f64,
    pub sst_gp: f64,
    pub nbl_upper: usize,
    pub nbl_lower: usize,
    pub iblte_upper: usize,
    pub iblte_lower: usize,
    pub itran_upper: usize,
    pub itran_lower: usize,
    pub xssitr_upper: f64,
    pub xssitr_lower: f64,
    pub lblini: bool,
    pub ipan_upper: Vec<usize>,
    pub ipan_lower: Vec<usize>,
    pub isys_upper: Vec<Option<usize>>,
    pub isys_lower: Vec<Option<usize>>,
    pub upper_rows: Vec<XfoilBlRow>,
    pub lower_rows: Vec<XfoilBlRow>,
    pub nsys: usize,
    pub va: Vec<VBlock>,
    pub vb: Vec<VBlock>,
    pub vm: Vec<Vec<VmRow>>,
    pub vz: [[f64; 2]; 3],
    pub vdel: Vec<VDelBlock>,
    pub usav_upper: Vec<f64>,
    pub usav_lower: Vec<f64>,
    pub u1_m: Vec<f64>,
    pub u2_m: Vec<f64>,
    pub d1_m: Vec<f64>,
    pub d2_m: Vec<f64>,
    pub ule1_m: Vec<f64>,
    pub ule2_m: Vec<f64>,
    pub ute1_m: Vec<f64>,
    pub ute2_m: Vec<f64>,
    pub u1_a: f64,
    pub u2_a: f64,
    pub d1_a: f64,
    pub d2_a: f64,
    pub due1: f64,
    pub due2: f64,
    pub dds1: f64,
    pub dds2: f64,
    pub u_new: Vec<f64>,
    pub u_ac: Vec<f64>,
    pub q_new: Vec<f64>,
    pub q_ac: Vec<f64>,
    pub dac: f64,
    pub rlx: f64,
    pub cl_new: f64,
    pub cl_a: f64,
    pub cl_ac: f64,
    pub cl_ms: f64,
    pub rmsbl: f64,
    pub rmxbl: f64,
    pub vmxbl: char,
    pub imxbl: usize,
    pub ismxbl: usize,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    pub cdp: f64,
    pub cdf: f64,
    pub converged: bool,
    pub iterations: usize,
    pub residual: f64,
}

impl XfoilState {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: String,
        alpha_rad: f64,
        operating_mode: OperatingMode,
        panel_x: Vec<f64>,
        panel_y: Vec<f64>,
        panel_s: Vec<f64>,
        qinvu_0: Vec<f64>,
        qinvu_90: Vec<f64>,
        gamu_0: Vec<f64>,
        gamu_90: Vec<f64>,
        dij: DMatrix<f64>,
    ) -> Self {
        let n = panel_x.len();
        Self {
            name,
            alpha_rad,
            operating_mode,
            panel_x,
            panel_y,
            panel_s,
            panel_xp: vec![0.0; n],
            panel_yp: vec![0.0; n],
            xle: 0.0,
            yle: 0.0,
            xte: 0.0,
            yte: 0.0,
            chord: 1.0,
            sharp: true,
            ante: 0.0,
            qinvu_0,
            qinvu_90,
            gamu_0,
            gamu_90,
            qinv: vec![0.0; n],
            qinv_a: vec![0.0; n],
            qvis: vec![0.0; n],
            gam: vec![0.0; n],
            gam_a: vec![0.0; n],
            wake_x: Vec::new(),
            wake_y: Vec::new(),
            wake_s: Vec::new(),
            wake_nx: Vec::new(),
            wake_ny: Vec::new(),
            wake_apanel: Vec::new(),
            wake_qinvu_0: Vec::new(),
            wake_qinvu_90: Vec::new(),
            wake_qinv: Vec::new(),
            wake_qinv_a: Vec::new(),
            wgap: Vec::new(),
            lwake: false,
            ladij: false,
            lwdij: false,
            awake: 0.0,
            dij,
            ist: 0,
            sst: 0.0,
            sst_go: 0.0,
            sst_gp: 0.0,
            nbl_upper: 0,
            nbl_lower: 0,
            iblte_upper: 0,
            iblte_lower: 0,
            itran_upper: 0,
            itran_lower: 0,
            xssitr_upper: 0.0,
            xssitr_lower: 0.0,
            lblini: false,
            ipan_upper: Vec::new(),
            ipan_lower: Vec::new(),
            isys_upper: Vec::new(),
            isys_lower: Vec::new(),
            upper_rows: Vec::new(),
            lower_rows: Vec::new(),
            nsys: 0,
            va: Vec::new(),
            vb: Vec::new(),
            vm: Vec::new(),
            vz: [[0.0; 2]; 3],
            vdel: Vec::new(),
            usav_upper: Vec::new(),
            usav_lower: Vec::new(),
            u1_m: Vec::new(),
            u2_m: Vec::new(),
            d1_m: Vec::new(),
            d2_m: Vec::new(),
            ule1_m: Vec::new(),
            ule2_m: Vec::new(),
            ute1_m: Vec::new(),
            ute2_m: Vec::new(),
            u1_a: 0.0,
            u2_a: 0.0,
            d1_a: 0.0,
            d2_a: 0.0,
            due1: 0.0,
            due2: 0.0,
            dds1: 0.0,
            dds2: 0.0,
            u_new: Vec::new(),
            u_ac: Vec::new(),
            q_new: vec![0.0; n],
            q_ac: vec![0.0; n],
            dac: 0.0,
            rlx: 1.0,
            cl_new: 0.0,
            cl_a: 0.0,
            cl_ac: 0.0,
            cl_ms: 0.0,
            rmsbl: 0.0,
            rmxbl: 0.0,
            vmxbl: ' ',
            imxbl: 0,
            ismxbl: 0,
            cl: 0.0,
            cd: 0.0,
            cm: 0.0,
            cdp: 0.0,
            cdf: 0.0,
            converged: false,
            iterations: 0,
            residual: 0.0,
        }
    }

    pub fn n_panel_nodes(&self) -> usize {
        self.panel_x.len()
    }

    pub fn n_total_nodes(&self) -> usize {
        self.n_panel_nodes() + self.wake_x.len()
    }

    pub fn surface_rows(&self, surface: XfoilSurface) -> &[XfoilBlRow] {
        match surface {
            XfoilSurface::Upper => &self.upper_rows,
            XfoilSurface::Lower => &self.lower_rows,
        }
    }

    pub fn surface_rows_mut(&mut self, surface: XfoilSurface) -> &mut Vec<XfoilBlRow> {
        match surface {
            XfoilSurface::Upper => &mut self.upper_rows,
            XfoilSurface::Lower => &mut self.lower_rows,
        }
    }

    pub fn allocate_newton_state(&mut self) {
        self.nsys = self.nbl_upper.saturating_sub(1) + self.nbl_lower.saturating_sub(1);
        self.va = vec![[[0.0; 2]; 3]; self.nsys + 1];
        self.vb = vec![[[0.0; 2]; 3]; self.nsys + 1];
        self.vm = vec![vec![[0.0; 3]; self.nsys + 1]; self.nsys + 1];
        self.vz = [[0.0; 2]; 3];
        self.vdel = vec![[[0.0; 2]; 3]; self.nsys + 1];

        let system_len = self.nsys + 1;
        self.u1_m = vec![0.0; system_len];
        self.u2_m = vec![0.0; system_len];
        self.d1_m = vec![0.0; system_len];
        self.d2_m = vec![0.0; system_len];
        self.ule1_m = vec![0.0; system_len];
        self.ule2_m = vec![0.0; system_len];
        self.ute1_m = vec![0.0; system_len];
        self.ute2_m = vec![0.0; system_len];

        self.usav_upper = vec![0.0; self.upper_rows.len()];
        self.usav_lower = vec![0.0; self.lower_rows.len()];
        self.u_new = vec![0.0; system_len];
        self.u_ac = vec![0.0; system_len];

        let q_len = self.n_total_nodes();
        self.q_new = vec![0.0; q_len];
        self.q_ac = vec![0.0; q_len];
    }
}
