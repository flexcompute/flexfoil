use rustfoil_bl::equations::FlowType;
use rustfoil_bl::state::{BlDerivatives, BlStation};

#[derive(Debug, Clone)]
pub struct CanonicalNewtonRow {
    pub x: f64,
    pub x_coord: f64,
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

impl CanonicalNewtonRow {
    pub fn as_station(&self, flow_type: FlowType) -> BlStation {
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
            is_turbulent,
            is_wake,
            derivs: self.derivs.clone(),
        }
    }

    pub fn overwrite_from_station(&mut self, station: &BlStation) {
        self.uedg = station.u;
        self.theta = station.theta;
        self.dstr = station.delta_star;
        self.ctau = station.ctau;
        self.ampl = station.ampl;
        self.mass = station.mass_defect;
        self.h = station.h;
        self.hk = station.hk;
        self.hs = station.hs;
        self.hc = station.hc;
        self.r_theta = station.r_theta;
        self.cf = station.cf;
        self.cd = station.cd;
        self.us = station.us;
        self.cq = station.cq;
        self.de = station.de;
        self.dw = station.dw;
        self.is_laminar = station.is_laminar;
        self.is_turbulent = station.is_turbulent;
        self.is_wake = station.is_wake;
        self.derivs = station.derivs.clone();
        self.ctau_or_ampl = if station.is_laminar && !station.is_turbulent {
            station.ampl
        } else {
            station.ctau
        };
    }
}

#[derive(Debug, Clone)]
pub struct CanonicalNewtonStateView {
    pub upper_rows: Vec<CanonicalNewtonRow>,
    pub lower_rows: Vec<CanonicalNewtonRow>,
    pub upper_flow_types: Vec<FlowType>,
    pub lower_flow_types: Vec<FlowType>,
    pub upper_ue_current: Vec<f64>,
    pub lower_ue_current: Vec<f64>,
    pub upper_ue_inviscid: Vec<f64>,
    pub lower_ue_inviscid: Vec<f64>,
    pub upper_ue_from_mass: Vec<f64>,
    pub lower_ue_from_mass: Vec<f64>,
    pub upper_ue_operating: Vec<f64>,
    pub lower_ue_operating: Vec<f64>,
    pub sst_go: f64,
    pub sst_gp: f64,
    pub ante: f64,
}

impl CanonicalNewtonStateView {
    pub fn upper_stations(&self) -> Vec<BlStation> {
        self.upper_rows
            .iter()
            .enumerate()
            .map(|(ibl, row)| {
                let flow_type = self
                    .upper_flow_types
                    .get(ibl.saturating_sub(1))
                    .copied()
                    .unwrap_or_else(|| if row.is_wake { FlowType::Wake } else if row.is_laminar { FlowType::Laminar } else { FlowType::Turbulent });
                row.as_station(flow_type)
            })
            .collect()
    }

    pub fn lower_stations(&self) -> Vec<BlStation> {
        self.lower_rows
            .iter()
            .enumerate()
            .map(|(ibl, row)| {
                let flow_type = self
                    .lower_flow_types
                    .get(ibl.saturating_sub(1))
                    .copied()
                    .unwrap_or_else(|| if row.is_wake { FlowType::Wake } else if row.is_laminar { FlowType::Laminar } else { FlowType::Turbulent });
                row.as_station(flow_type)
            })
            .collect()
    }

    pub fn overwrite_from_stations(
        &mut self,
        upper_stations: &[BlStation],
        lower_stations: &[BlStation],
    ) {
        for (row, station) in self.upper_rows.iter_mut().zip(upper_stations.iter()) {
            row.overwrite_from_station(station);
        }
        for (row, station) in self.lower_rows.iter_mut().zip(lower_stations.iter()) {
            row.overwrite_from_station(station);
        }
        self.upper_ue_current = self.upper_rows.iter().map(|row| row.uedg).collect();
        self.lower_ue_current = self.lower_rows.iter().map(|row| row.uedg).collect();
    }
}
