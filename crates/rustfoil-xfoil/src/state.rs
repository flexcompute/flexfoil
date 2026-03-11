use nalgebra::DMatrix;
use rustfoil_bl::state::BlDerivatives;

use crate::config::OperatingMode;

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
    pub sharp: bool,
    pub ante: f64,
    pub qinvu_0: Vec<f64>,
    pub qinvu_90: Vec<f64>,
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
    pub ipan_upper: Vec<usize>,
    pub ipan_lower: Vec<usize>,
    pub isys_upper: Vec<Option<usize>>,
    pub isys_lower: Vec<Option<usize>>,
    pub upper_rows: Vec<XfoilBlRow>,
    pub lower_rows: Vec<XfoilBlRow>,
    pub u_new: Vec<f64>,
    pub u_ac: Vec<f64>,
    pub dac: f64,
    pub rlx: f64,
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
            sharp: true,
            ante: 0.0,
            qinvu_0,
            qinvu_90,
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
            ipan_upper: Vec::new(),
            ipan_lower: Vec::new(),
            isys_upper: Vec::new(),
            isys_lower: Vec::new(),
            upper_rows: Vec::new(),
            lower_rows: Vec::new(),
            u_new: Vec::new(),
            u_ac: Vec::new(),
            dac: 0.0,
            rlx: 1.0,
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
}
