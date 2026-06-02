//! Multi-body inviscid panel method (coupled multi-element).
//!
//! Generalizes the single-body stream-function panel method ([`crate::system`])
//! to N disjoint bodies solved in one coupled system. Each body is its own
//! streamline (one ψ₀ unknown and one Kutta condition per body); the influence
//! matrix sums every body's panel contributions at every control point, which is
//! the element-to-element interaction (slot/circulation coupling).
//!
//! With a single body this reduces *exactly* to [`crate::system::build_and_factorize`]:
//! the diagonal block is the same `psilin(geom, i, …)` per row. Cross-body blocks
//! evaluate one body's panels at a control point on another via `psilin` with an
//! out-of-range node index (the external/wake convention already used for wake
//! points in the DIJ code).
//!
//! Inviscid only (Cd = 0); no wake panels in this first cut.

use crate::geometry::AirfoilGeometry;
use crate::influence::{psilin, psilin_with_dqdm};
use crate::solution::FlowConditions;
use crate::{InviscidError, Result};
use nalgebra::{DMatrix, DVector};

/// Per-body inviscid result.
#[derive(Debug, Clone)]
pub struct BodyInviscidResult {
    /// Vortex strength (= surface tangential velocity) at each node.
    pub gamma: Vec<f64>,
    /// Pressure coefficient at each node.
    pub cp: Vec<f64>,
    /// Lift coefficient (normalized by this body's chord).
    pub cl: f64,
    /// Moment coefficient about this body's quarter-chord.
    pub cm: f64,
}

/// Coupled solution for all bodies at one angle of attack.
#[derive(Debug, Clone)]
pub struct MultiInviscidSolution {
    /// One result per body, in input order.
    pub bodies: Vec<BodyInviscidResult>,
    /// Internal stream-function constant ψ₀ per body.
    pub psi0: Vec<f64>,
}

/// Factorized coupled system, ready for α sweeps (two base solutions combined).
#[derive(Debug, Clone)]
pub struct FactorizedMultiSystem {
    geoms: Vec<AirfoilGeometry>,
    /// Node offset of each body in the global γ vector; `offsets[k]..offsets[k]+geoms[k].n`.
    offsets: Vec<usize>,
    n_total: usize,
    /// Global γ basis for α = 0° and α = 90° (length `n_total`).
    gamu_0: Vec<f64>,
    gamu_90: Vec<f64>,
    /// ψ₀ basis per body for α = 0° and α = 90° (length `K`).
    psi0_0: Vec<f64>,
    psi0_90: Vec<f64>,
}

fn body_offsets(geoms: &[AirfoilGeometry]) -> (Vec<usize>, usize) {
    let mut offsets = Vec::with_capacity(geoms.len());
    let mut acc = 0;
    for g in geoms {
        offsets.push(acc);
        acc += g.n;
    }
    (offsets, acc)
}

/// Build and factorize the coupled influence system for `geoms`.
pub fn build_and_factorize_multi(geoms: &[AirfoilGeometry]) -> Result<FactorizedMultiSystem> {
    if geoms.is_empty() {
        return Err(InviscidError::SingularMatrix);
    }
    let k = geoms.len();
    let (offsets, n_total) = body_offsets(geoms);
    let size = n_total + k;

    let mut a = DMatrix::<f64>::zeros(size, size);
    let mut rhs_0 = DVector::<f64>::zeros(size);
    let mut rhs_90 = DVector::<f64>::zeros(size);

    // Boundary-condition rows: ψ at each node = that body's ψ₀.
    for (a_idx, ga) in geoms.iter().enumerate() {
        let off_a = offsets[a_idx];
        for i in 0..ga.n {
            let row = off_a + i;
            let xi = ga.x[i];
            let yi = ga.y[i];
            for (b_idx, gb) in geoms.iter().enumerate() {
                // Self block handles its own singularity via node index `i`;
                // other blocks evaluate at an external point (index ≥ n).
                let res = if b_idx == a_idx {
                    psilin(gb, i, xi, yi)
                } else {
                    psilin(gb, gb.n + 1, xi, yi)
                };
                let off_b = offsets[b_idx];
                for j in 0..gb.n {
                    a[(row, off_b + j)] = res.dzdg[j];
                }
            }
            a[(row, n_total + a_idx)] = -1.0; // −ψ₀ of this body
            rhs_0[row] = -yi;
            rhs_90[row] = xi;
        }
    }

    // Per body: replace the last node's row with the sharp-TE bisector velocity
    // equation (mirrors single-body GGCALC). Blunt-TE bodies return None and keep
    // their boundary row.
    for (a_idx, ga) in geoms.iter().enumerate() {
        let Some((xbis, ybis, nxbis, nybis)) = ga.sharp_te_bisector_control() else {
            continue;
        };
        let off_a = offsets[a_idx];
        let te_row = off_a + ga.n - 1;
        for c in 0..size {
            a[(te_row, c)] = 0.0;
        }
        for (b_idx, gb) in geoms.iter().enumerate() {
            let res = psilin_with_dqdm(gb, gb.n, xbis, ybis, nxbis, nybis);
            let off_b = offsets[b_idx];
            for j in 0..gb.n {
                a[(te_row, off_b + j)] = res.dqdg[j];
            }
        }
        rhs_0[te_row] = -nybis;
        rhs_90[te_row] = nxbis;
    }

    // One Kutta condition per body: γ(upper TE) + γ(lower TE) = 0.
    for (a_idx, ga) in geoms.iter().enumerate() {
        let row = n_total + a_idx;
        let off_a = offsets[a_idx];
        a[(row, off_a)] = 1.0;
        a[(row, off_a + ga.n - 1)] = 1.0;
    }

    let lu = a.lu();
    let sol_0 = lu.solve(&rhs_0).ok_or(InviscidError::SingularMatrix)?;
    let sol_90 = lu.solve(&rhs_90).ok_or(InviscidError::SingularMatrix)?;

    Ok(FactorizedMultiSystem {
        geoms: geoms.to_vec(),
        offsets,
        n_total,
        gamu_0: sol_0.rows(0, n_total).iter().copied().collect(),
        gamu_90: sol_90.rows(0, n_total).iter().copied().collect(),
        psi0_0: sol_0.rows(n_total, k).iter().copied().collect(),
        psi0_90: sol_90.rows(n_total, k).iter().copied().collect(),
    })
}

impl FactorizedMultiSystem {
    /// Combine the two base solutions for a given angle of attack.
    pub fn solve_alpha(&self, flow: &FlowConditions) -> MultiInviscidSolution {
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();

        let gamma: Vec<f64> = (0..self.n_total)
            .map(|i| cosa * self.gamu_0[i] + sina * self.gamu_90[i])
            .collect();

        let bodies = self
            .geoms
            .iter()
            .enumerate()
            .map(|(k, g)| {
                let off = self.offsets[k];
                let g_body: Vec<f64> = gamma[off..off + g.n].to_vec();
                let cp: Vec<f64> = g_body.iter().map(|&q| 1.0 - (q / flow.v_inf).powi(2)).collect();
                let (cl, cm) = body_forces(g, &cp, flow);
                BodyInviscidResult { gamma: g_body, cp, cl, cm }
            })
            .collect();

        let psi0: Vec<f64> = (0..self.geoms.len())
            .map(|k| cosa * self.psi0_0[k] + sina * self.psi0_90[k])
            .collect();

        MultiInviscidSolution { bodies, psi0 }
    }
}

/// Lift/moment for one body by pressure integration (mirrors `compute_forces`).
fn body_forces(geom: &AirfoilGeometry, cp: &[f64], flow: &FlowConditions) -> (f64, f64) {
    let n = geom.n;
    let cosa = flow.alpha.cos();
    let sina = flow.alpha.sin();
    let x_ref = 0.25 * geom.chord;
    let mut cl = 0.0;
    let mut cm = 0.0;
    for i in 0..n {
        let ip = (i + 1) % n;
        let dx = geom.x[ip] - geom.x[i];
        let dy = geom.y[ip] - geom.y[i];
        let dx_wind = dx * cosa + dy * sina;
        let dy_wind = dy * cosa - dx * sina;
        let cp_avg = 0.5 * (cp[i] + cp[ip]);
        let dg = cp[ip] - cp[i];
        cl += cp_avg * dx_wind;
        let x_mid = 0.5 * (geom.x[i] + geom.x[ip]);
        let y_mid = 0.5 * (geom.y[i] + geom.y[ip]);
        let ax = (x_mid - x_ref) * cosa + y_mid * sina;
        let ay = y_mid * cosa - (x_mid - x_ref) * sina;
        cm -= cp_avg * (ax * dx_wind / geom.chord + ay * dy_wind / geom.chord);
        cm -= dg * dx_wind * dx_wind / (12.0 * geom.chord);
        cm -= dg * dy_wind * dy_wind / (12.0 * geom.chord);
    }
    (cl, cm)
}

/// Convenience: factorize and solve at one angle of attack.
pub fn solve_multi(geoms: &[AirfoilGeometry], flow: &FlowConditions) -> Result<MultiInviscidSolution> {
    Ok(build_and_factorize_multi(geoms)?.solve_alpha(flow))
}
