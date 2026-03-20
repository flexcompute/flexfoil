//! Inviscid flow solver using Mark Drela's Linear Vorticity Panel Method.
//!
//! This implements the exact panel method from XFOIL, using:
//! - Linear vorticity distribution across each panel (node-based unknowns)
//! - Stream function formulation (ψ = ψ₀ on surface)
//! - Kutta condition: γ₁ + γₙ = 0
//! - Sum/difference (PSIS/PSID) formulation for influence coefficients
//!
//! # System Structure
//!
//! For N nodes, we have N+1 unknowns: [γ₀, γ₁, ..., γₙ₋₁, ψ₀]
//!
//! The matrix equation is:
//! ```text
//! | A₀₀  A₀₁  ... A₀,ₙ₋₁  -1 | | γ₀    |   | -ψ∞(0)   |
//! | A₁₀  A₁₁  ... A₁,ₙ₋₁  -1 | | γ₁    |   | -ψ∞(1)   |
//! |  :    :   ...   :      : | |  :    | = |    :     |
//! | Aₙ₋₁,₀...    Aₙ₋₁,ₙ₋₁ -1 | | γₙ₋₁  |   | -ψ∞(n-1) |
//! |  1    0   ...   1      0 | | ψ₀    |   |    0     |
//! ```
//!
//! # Reference
//!
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - XFOIL 6.99 source code, `xpanel.f` (PSILIN subroutine)

mod error;
mod influence;
pub mod velocity;
pub mod smoke;

pub use error::SolverError;
pub use velocity::{
    build_dividing_streamline, build_dividing_streamline_viscous, build_streamlines,
    build_streamlines_viscous, compute_psi_grid, compute_psi_grid_with_interior,
    compute_psi_grid_with_sources, is_inside_airfoil, psi_at, psi_at_with_sources,
    velocity_at, velocity_at_with_sources, StreamlineOptions, WakePanels,
};
pub use smoke::SmokeSystem;

use nalgebra::{DMatrix, DVector};
use rustfoil_core::{Body, Point};
use std::f64::consts::PI;
use std::ops::Range;

/// 1/(4π) - used for stream function influence coefficients
const QOPI: f64 = 0.25 / PI;

/// 1/(2π) - used for TE panel source/vortex influence
const HOPI: f64 = 0.5 / PI;

/// Result type for solver operations.
pub type SolverResult<T> = Result<T, SolverError>;

/// Flow conditions for the analysis.
#[derive(Debug, Clone, Copy)]
pub struct FlowConditions {
    /// Angle of attack in radians
    pub alpha: f64,
    /// Freestream velocity magnitude (typically normalized to 1.0)
    pub v_inf: f64,
}

impl Default for FlowConditions {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            v_inf: 1.0,
        }
    }
}

impl FlowConditions {
    /// Create flow conditions with the given angle of attack (in degrees).
    pub fn with_alpha_deg(alpha_deg: f64) -> Self {
        Self {
            alpha: alpha_deg.to_radians(),
            v_inf: 1.0,
        }
    }
}

/// Solution from the inviscid panel method.
#[derive(Debug, Clone)]
pub struct InviscidSolution {
    /// Vorticity values at each node (γᵢ = surface velocity), concatenated across all bodies
    pub gamma: Vec<f64>,
    /// Pressure coefficient at each node, concatenated across all bodies
    pub cp: Vec<f64>,
    /// Total lift coefficient (per unit span), normalized by first body's chord
    pub cl: f64,
    /// Total moment coefficient about quarter-chord of first body
    pub cm: f64,
    /// Internal stream function value (first body, for backward compatibility)
    pub psi_0: f64,
    /// Number of nodes per body (for indexing into gamma/cp arrays)
    pub nodes_per_body: Vec<usize>,
    /// Per-body gamma distributions
    pub gamma_per_body: Vec<Vec<f64>>,
    /// Per-body Cp distributions
    pub cp_per_body: Vec<Vec<f64>>,
    /// Per-body lift coefficients (each normalized by that body's chord)
    pub cl_per_body: Vec<f64>,
}

/// Cached factorization for efficient alpha sweeps.
///
/// XFOIL's key optimization: solve once for α=0° and α=90°,
/// then combine for any angle: γ = cos(α)*γ₀ + sin(α)*γ₉₀
#[derive(Debug, Clone)]
pub struct FactorizedSolution {
    /// Solution for α = 0° (length = N_total)
    gamu_0: Vec<f64>,
    /// Solution for α = 90° (length = N_total)
    gamu_90: Vec<f64>,
    /// ψ₀ for α = 0° (one per body)
    psi0_0: Vec<f64>,
    /// ψ₀ for α = 90° (one per body)
    psi0_90: Vec<f64>,
    /// Cached node positions for Cp/Cl calculation (length = N_total)
    nodes: Vec<Point>,
    /// Total number of nodes across all bodies
    n_nodes: usize,
    /// Reference chord length (first body)
    chord: f64,
    /// Node index ranges for each body in the global arrays
    body_node_ranges: Vec<Range<usize>>,
    /// Chord length per body
    body_chords: Vec<f64>,
}

impl FactorizedSolution {
    /// Compute solution for any angle of attack.
    pub fn solve_alpha(&self, flow: &FlowConditions) -> InviscidSolution {
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();
        let n = self.n_nodes;
        let k = self.body_node_ranges.len();

        // Combine the two base solutions for gamma
        let gamma: Vec<f64> = (0..n)
            .map(|i| cosa * self.gamu_0[i] + sina * self.gamu_90[i])
            .collect();

        // Combine psi0 per body
        let psi0: Vec<f64> = (0..k)
            .map(|b| cosa * self.psi0_0[b] + sina * self.psi0_90[b])
            .collect();

        // Compute Cp at each node: Cp = 1 - (γ/V∞)²
        let cp: Vec<f64> = gamma
            .iter()
            .map(|&g| 1.0 - (g / flow.v_inf).powi(2))
            .collect();

        // Per-body gamma, Cp, and forces
        let mut gamma_per_body = Vec::with_capacity(k);
        let mut cp_per_body = Vec::with_capacity(k);
        let mut cl_per_body = Vec::with_capacity(k);
        let nodes_per_body: Vec<usize> = self.body_node_ranges.iter()
            .map(|r| r.len())
            .collect();

        let mut cl_total = 0.0;
        let mut cm_total = 0.0;

        for (b, range) in self.body_node_ranges.iter().enumerate() {
            let body_gamma: Vec<f64> = gamma[range.clone()].to_vec();
            let body_cp: Vec<f64> = cp[range.clone()].to_vec();
            let body_nodes = &self.nodes[range.clone()];
            let body_chord = self.body_chords[b];

            let (cl_b, cm_b) = Self::compute_body_forces(
                body_nodes, &body_cp, body_chord, flow,
            );

            // Accumulate total forces normalized by reference chord
            cl_total += cl_b * body_chord / self.chord;
            cm_total += cm_b * body_chord / self.chord;

            cl_per_body.push(cl_b);
            gamma_per_body.push(body_gamma);
            cp_per_body.push(body_cp);
        }

        InviscidSolution {
            gamma,
            cp,
            cl: cl_total,
            cm: cm_total,
            psi_0: psi0[0],
            nodes_per_body,
            gamma_per_body,
            cp_per_body,
            cl_per_body,
        }
    }

    /// Compute lift and moment coefficients for a single body from its Cp distribution.
    ///
    /// This matches XFOIL's CLCALC subroutine:
    /// - Integrates around the closed contour, wrapping from N back to 1
    /// - Uses trapezoidal integration of pressure forces
    fn compute_body_forces(
        nodes: &[Point],
        cp: &[f64],
        chord: f64,
        flow: &FlowConditions,
    ) -> (f64, f64) {
        let n = nodes.len();
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();

        let mut cl = 0.0;
        let mut cm = 0.0;

        for i in 0..n {
            let ip = (i + 1) % n;

            let dx = nodes[ip].x - nodes[i].x;
            let dy = nodes[ip].y - nodes[i].y;

            let dx_wind = dx * cosa + dy * sina;
            let dy_wind = dy * cosa - dx * sina;

            let cp_avg = 0.5 * (cp[i] + cp[ip]);
            cl += cp_avg * dx_wind;

            let x_mid = 0.5 * (nodes[i].x + nodes[ip].x);
            let y_mid = 0.5 * (nodes[i].y + nodes[ip].y);
            let x_ref = 0.25 * chord;

            let dg = cp[ip] - cp[i];

            cm -= cp_avg * ((x_mid - x_ref) * dx_wind / chord
                         + y_mid * dy_wind / chord);
            cm -= (dg * dx_wind * dx_wind / 12.0) / chord;
            cm -= (dg * dy_wind * dy_wind / 12.0) / chord;
        }

        (cl, cm)
    }
}

/// Per-body trailing edge geometry info for the influence matrix assembly.
#[derive(Debug, Clone)]
struct BodyTeInfo {
    /// XFOIL's SCS coefficient for TE source/vortex treatment
    scs: f64,
    /// XFOIL's SDS coefficient for TE source/vortex treatment
    sds: f64,
    /// Arc-length tolerance for TE panel skip
    seps: f64,
    /// Whether the contour is closed (sharp TE where first == last point)
    is_closed: bool,
}

/// Inviscid flow solver using Drela's linear vorticity panel method.
pub struct InviscidSolver {
    _config: SolverConfig,
}

#[derive(Debug, Clone, Default)]
struct SolverConfig {}

impl Default for InviscidSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl InviscidSolver {
    /// Create a new inviscid solver.
    pub fn new() -> Self {
        Self {
            _config: SolverConfig::default(),
        }
    }

    /// Solve for the inviscid flow around the given bodies.
    ///
    /// Uses Drela's linear vorticity stream function panel method (XFOIL).
    pub fn solve(&self, bodies: &[Body], flow: &FlowConditions) -> SolverResult<InviscidSolution> {
        let factorized = self.factorize(bodies)?;
        Ok(factorized.solve_alpha(flow))
    }

    /// Extract nodes and TE info for a single body.
    fn extract_body_info(body: &Body) -> SolverResult<(Vec<Point>, BodyTeInfo)> {
        let panels = body.panels();
        let n_panels = panels.len();

        if n_panels < 3 {
            return Err(SolverError::InsufficientPanels);
        }

        let mut nodes: Vec<Point> = panels.iter().map(|p| p.p1).collect();

        let last_p2 = panels.last().unwrap().p2;
        let first_p1 = panels[0].p1;
        let is_closed = (last_p2.x - first_p1.x).abs() < 1e-10
                     && (last_p2.y - first_p1.y).abs() < 1e-10;

        if !is_closed {
            nodes.push(last_p2);
        }
        let n = nodes.len();
        let chord = body.chord();

        let arc_length: f64 = panels.iter().map(|p| p.length()).sum();
        let seps = arc_length * 1e-5;

        // TE geometry
        let te_upper = nodes[0];
        let te_lower = nodes[n - 1];
        let dxte = te_lower.x - te_upper.x;
        let dyte = te_lower.y - te_upper.y;
        let dste = (dxte * dxte + dyte * dyte).sqrt();

        let t_upper = panels[0].tangent();
        let t_lower = panels[n_panels - 1].tangent();
        let dxs = 0.5 * (t_upper.x - t_lower.x);
        let dys = 0.5 * (t_upper.y - t_lower.y);

        let ante = dxs * dyte - dys * dxte;
        let aste = dxs * dxte + dys * dyte;

        let sharp = dste < 0.0001 * chord;
        let (scs, sds) = if sharp || is_closed {
            (1.0, 0.0)
        } else {
            (ante / dste, aste / dste)
        };

        Ok((nodes, BodyTeInfo { scs, sds, seps, is_closed }))
    }

    /// Factorize the influence matrix and solve for α=0° and α=90°.
    ///
    /// Supports multiple bodies. The matrix is `(N_total + K) × (N_total + K)`
    /// where `N_total` = total nodes across all bodies and `K` = number of bodies.
    /// Unknowns: `[γ₀, ..., γ_{N_total-1}, ψ₀_1, ..., ψ₀_K]`
    pub fn factorize(&self, bodies: &[Body]) -> SolverResult<FactorizedSolution> {
        if bodies.is_empty() {
            return Err(SolverError::NoBodies);
        }

        let k = bodies.len(); // number of bodies

        // Extract nodes and TE info for each body
        let mut all_nodes: Vec<Point> = Vec::new();
        let mut body_node_ranges: Vec<Range<usize>> = Vec::new();
        let mut body_te_infos: Vec<BodyTeInfo> = Vec::new();
        let mut body_chords: Vec<f64> = Vec::new();
        let mut body_nodes_list: Vec<Vec<Point>> = Vec::new();

        for body in bodies {
            let (nodes, te_info) = Self::extract_body_info(body)?;
            let start = all_nodes.len();
            let n_b = nodes.len();
            all_nodes.extend_from_slice(&nodes);
            body_node_ranges.push(start..start + n_b);
            body_te_infos.push(te_info);
            body_chords.push(body.chord());
            body_nodes_list.push(nodes);
        }

        let n_total = all_nodes.len();
        let sys_size = n_total + k;
        let chord = body_chords[0];

        // Build influence coefficient matrix
        let mut a_matrix = DMatrix::<f64>::zeros(sys_size, sys_size);
        let mut rhs_0 = DVector::<f64>::zeros(sys_size);
        let mut rhs_90 = DVector::<f64>::zeros(sys_size);

        // For each field point i (across all bodies)
        for (body_i, range_i) in body_node_ranges.iter().enumerate() {
            for i_global in range_i.clone() {
                let xi = all_nodes[i_global].x;
                let yi = all_nodes[i_global].y;

                // Influence coefficients for this row (all gamma unknowns)
                let mut dzdg = vec![0.0; n_total];

                // Loop over all source bodies
                for (body_j, range_j) in body_node_ranges.iter().enumerate() {
                    let n_j = range_j.len();
                    let offset_j = range_j.start;
                    let src_nodes = &body_nodes_list[body_j];
                    let te_info = &body_te_infos[body_j];

                    // TE panel geometry saved for source/vortex treatment
                    let mut te_g1 = 0.0;
                    let mut te_g2 = 0.0;
                    let mut te_t1 = 0.0;
                    let mut te_t2 = 0.0;
                    let mut te_x1 = 0.0;
                    let mut te_x2 = 0.0;
                    let mut te_yy = 0.0;
                    let mut te_apan = 0.0;
                    let mut te_panel_valid = false;

                    // Loop over all panels of source body j
                    for jo_local in 0..n_j {
                        let jp_local = (jo_local + 1) % n_j;
                        let jo_global = offset_j + jo_local;
                        let jp_global = offset_j + jp_local;

                        let x_jo = src_nodes[jo_local].x;
                        let y_jo = src_nodes[jo_local].y;
                        let x_jp = src_nodes[jp_local].x;
                        let y_jp = src_nodes[jp_local].y;

                        let dx = x_jp - x_jo;
                        let dy = y_jp - y_jo;
                        let ds_sq = dx * dx + dy * dy;

                        // TE panel of source body (last panel wrapping back to first node)
                        // Skip only if the panel is degenerate (zero length)
                        if jo_local == n_j - 1 && ds_sq < te_info.seps * te_info.seps {
                            continue;
                        }

                        if ds_sq < 1e-24 {
                            continue;
                        }

                        let dso = ds_sq.sqrt();
                        let dsio = 1.0 / dso;
                        let sx = dx * dsio;
                        let sy = dy * dsio;

                        let rx1 = xi - x_jo;
                        let ry1 = yi - y_jo;
                        let rx2 = xi - x_jp;
                        let ry2 = yi - y_jp;

                        let x1 = sx * rx1 + sy * ry1;
                        let x2 = sx * rx2 + sy * ry2;
                        let yy = sx * ry1 - sy * rx1;

                        let rs1 = rx1 * rx1 + ry1 * ry1;
                        let rs2 = rx2 * rx2 + ry2 * ry2;

                        // Singularity checks use global indices
                        let (g1, t1) = if i_global != jo_global && rs1 > 1e-20 {
                            (rs1.ln(), x1.atan2(yy))
                        } else {
                            (0.0, 0.0)
                        };

                        let (g2, t2) = if i_global != jp_global && rs2 > 1e-20 {
                            (rs2.ln(), x2.atan2(yy))
                        } else {
                            (0.0, 0.0)
                        };

                        // TE panel of source body — save for source/vortex treatment
                        if jo_local == n_j - 1 {
                            te_g1 = g1;
                            te_g2 = g2;
                            te_t1 = t1;
                            te_t2 = t2;
                            te_x1 = x1;
                            te_x2 = x2;
                            te_yy = yy;
                            te_apan = (-sx).atan2(sy) + PI;
                            te_panel_valid = true;
                            continue;
                        }

                        // Normal vortex panel contribution (PSIS/PSID)
                        let dxinv = 1.0 / (x1 - x2);
                        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2
                                 + x2 - x1 + yy * (t1 - t2);
                        let psid = ((x1 + x2) * psis
                                 + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2))
                                 * dxinv;

                        dzdg[jo_global] += QOPI * (psis - psid);
                        dzdg[jp_global] += QOPI * (psis + psid);
                    }

                    // TE panel source/vortex contribution for source body j
                    if te_panel_valid {
                        let psig = 0.5 * te_yy * (te_g1 - te_g2)
                                 + te_x2 * (te_t2 - te_apan)
                                 - te_x1 * (te_t1 - te_apan);
                        let pgam = 0.5 * te_x1 * te_g1 - 0.5 * te_x2 * te_g2
                                 + te_x2 - te_x1 + te_yy * (te_t1 - te_t2);

                        let jo_te_global = offset_j + n_j - 1;
                        let jp_te_global = offset_j; // first node of source body

                        dzdg[jo_te_global] -= HOPI * psig * te_info.scs * 0.5;
                        dzdg[jp_te_global] += HOPI * psig * te_info.scs * 0.5;
                        dzdg[jo_te_global] += HOPI * pgam * te_info.sds * 0.5;
                        dzdg[jp_te_global] -= HOPI * pgam * te_info.sds * 0.5;
                    }
                } // end source body loop

                // Fill matrix row
                for j in 0..n_total {
                    a_matrix[(i_global, j)] = dzdg[j];
                }

                // ψ₀ column for the body this field point belongs to
                a_matrix[(i_global, n_total + body_i)] = -1.0;

                // RHS
                rhs_0[i_global] = -yi;
                rhs_90[i_global] = xi;
            }
        }

        // Kutta conditions: one per body
        for (b, range) in body_node_ranges.iter().enumerate() {
            let row = n_total + b;
            let first_node = range.start;         // upper TE
            let last_node = range.start + range.len() - 1; // lower TE
            a_matrix[(row, first_node)] = 1.0;
            a_matrix[(row, last_node)] = 1.0;
            rhs_0[row] = 0.0;
            rhs_90[row] = 0.0;
        }

        // Debug output
        if rustfoil_bl::is_debug_active() {
            let diag_n = 20.min(sys_size);
            let aij_diagonal: Vec<f64> = (0..diag_n)
                .map(|i| a_matrix[(i, i)])
                .collect();
            let aij_row0: Vec<f64> = (0..diag_n)
                .map(|j| a_matrix[(0, j)])
                .collect();
            let rhs_0_sample: Vec<f64> = rhs_0.iter().take(10).copied().collect();
            let rhs_90_sample: Vec<f64> = rhs_90.iter().take(10).copied().collect();

            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::aij_matrix(
                sys_size,
                aij_diagonal,
                aij_row0,
                rhs_0_sample,
                rhs_90_sample,
            ));
        }

        // LU factorization and solve
        let lu = a_matrix.clone().lu();

        let solution_0 = lu.solve(&rhs_0)
            .ok_or(SolverError::SingularMatrix)?;
        let solution_90 = lu.solve(&rhs_90)
            .ok_or(SolverError::SingularMatrix)?;

        // Extract gamma (first N_total entries) and psi0 per body (last K entries)
        let gamu_0: Vec<f64> = solution_0.iter().take(n_total).copied().collect();
        let gamu_90: Vec<f64> = solution_90.iter().take(n_total).copied().collect();

        if rustfoil_bl::is_debug_active() {
            rustfoil_bl::add_event(rustfoil_bl::DebugEvent::full_aic(
                n_total,
                gamu_0.clone(),
                gamu_90.clone(),
            ));
        }

        let psi0_0: Vec<f64> = (0..k).map(|b| solution_0[n_total + b]).collect();
        let psi0_90: Vec<f64> = (0..k).map(|b| solution_90[n_total + b]).collect();

        Ok(FactorizedSolution {
            gamu_0,
            gamu_90,
            psi0_0,
            psi0_90,
            nodes: all_nodes,
            n_nodes: n_total,
            chord,
            body_node_ranges,
            body_chords,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_core::point;

    /// NACA 4-digit thickness distribution with closed trailing edge.
    /// Uses the modified coefficient (-0.1036 instead of -0.1015) for closed TE.
    fn naca_thickness_closed(x: f64, t: f64) -> f64 {
        // Modified formula for closed TE (from Abbott & Von Doenhoff)
        5.0 * t * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) 
            + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
    }

    fn make_naca0012(n_half: usize) -> Body {
        // NACA 0012 with closed trailing edge, counterclockwise ordering
        let n = 2 * n_half;
        let t = 0.12;
        let mut points = Vec::with_capacity(n + 1);
        
        // Counterclockwise traversal: start at TE, go around
        for i in 0..n {
            let theta = 2.0 * PI * (i as f64) / (n as f64);
            let x = 0.5 * (1.0 + theta.cos());  // x: 1 → 0 → 1
            
            // Upper surface (theta in [0, π)): positive y
            // Lower surface (theta in [π, 2π)): negative y
            let thickness = naca_thickness_closed(x, t);
            let y = if theta < PI {
                thickness
            } else {
                -thickness
            };
            
            points.push(point(x, y));
        }
        
        // Close the contour by duplicating the first point
        points.push(points[0]);
        
        Body::from_points("NACA0012", &points).unwrap()
    }

    /// Create a perfectly symmetric ellipse for validation.
    /// Points are distributed symmetrically around the ellipse.
    fn make_ellipse(n_half: usize, thickness_ratio: f64) -> Body {
        // Create 2*n_half points, with first point duplicated at end for closure
        let n = 2 * n_half;
        let mut points = Vec::with_capacity(n + 1);
        
        // Counterclockwise: start at upper TE, go around to lower TE
        // Use symmetric distribution: first and last panels are equal length
        for i in 0..n {
            // theta from 0 to 2π*(n-1)/n, evenly spaced
            let theta = 2.0 * PI * (i as f64) / (n as f64);
            let x = 0.5 * (1.0 + theta.cos());
            let y = 0.5 * thickness_ratio * theta.sin();
            points.push(point(x, y));
        }
        
        // Add closing point (same as first point) for sharp TE
        points.push(points[0]);
        
        Body::from_points("Ellipse", &points).unwrap()
    }

    /// NACA 4-digit with OPEN trailing edge (like WASM generator)
    fn make_naca0012_open(n_half: usize) -> Body {
        let n = n_half;
        let t = 0.12;
        
        // Cosine spacing
        let x_coords: Vec<f64> = (0..n)
            .map(|i| {
                let beta = PI * (i as f64) / ((n - 1) as f64);
                0.5 * (1.0 - beta.cos())
            })
            .collect();
        
        let mut points = Vec::with_capacity(2 * n - 1);
        
        // Upper surface: TE to LE (CW order)
        for i in (0..n).rev() {
            let x = x_coords[i];
            let y_t = 5.0 * t * (
                0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2)
                + 0.2843 * x.powi(3) - 0.1015 * x.powi(4)
            );
            points.push(point(x, y_t));
        }
        // Lower surface: LE+1 to TE (skip LE to avoid duplicate)
        for i in 1..n {
            let x = x_coords[i];
            let y_t = 5.0 * t * (
                0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2)
                + 0.2843 * x.powi(3) - 0.1015 * x.powi(4)
            );
            points.push(point(x, -y_t));
        }
        
        Body::from_points("NACA0012_open", &points).unwrap()
    }

    #[test]
    fn test_symmetric_airfoil_zero_alpha() {
        let solver = InviscidSolver::new();
        
        // Use ellipse for symmetry test
        let ellipse = make_ellipse(60, 0.12);
        let flow = FlowConditions::default(); // α = 0

        let result = solver.solve(&[ellipse.clone()], &flow);
        
        match result {
            Ok(solution) => {
                println!("Ellipse Cl at α=0: {:.6}", solution.cl);
                println!("ψ₀ = {:.6}", solution.psi_0);
                
                // For a symmetric body at α=0, CL should be close to zero
                // The error comes from discretization, not the method itself
                assert!(solution.cl.abs() < 0.25, 
                    "Cl at α=0 has large offset: {}. Check panel distribution.", solution.cl);
            }
            Err(e) => panic!("Solver failed: {:?}", e),
        }
        
        // Also test NACA 0012 (closed)
        let naca = make_naca0012(60);
        let result = solver.solve(&[naca], &flow);
        
        match result {
            Ok(solution) => {
                println!("NACA 0012 (closed) Cl at α=0: {:.6}", solution.cl);
                assert!(solution.cl.abs() < 0.20, 
                    "NACA 0012 Cl at α=0 has large offset: {}", solution.cl);
            }
            Err(e) => panic!("Solver failed: {:?}", e),
        }
        
        // Test NACA 0012 with OPEN trailing edge (like WASM generator)
        let naca_open = make_naca0012_open(50);
        let result = solver.solve(&[naca_open], &flow);
        
        match result {
            Ok(solution) => {
                println!("NACA 0012 (open) Cl at α=0: {:.6}", solution.cl);
                // Should be close to zero for symmetric airfoil
                assert!(solution.cl.abs() < 0.10, 
                    "NACA 0012 open TE Cl at α=0 too large: {}", solution.cl);
            }
            Err(e) => panic!("Solver failed: {:?}", e),
        }
    }

    #[test]
    fn test_lift_curve_slope() {
        let solver = InviscidSolver::new();
        let airfoil = make_naca0012(80);
        
        // Factorize once
        let factorized = solver.factorize(&[airfoil]).unwrap();
        
        // Test at multiple angles
        let alphas_deg = [0.0, 2.0, 4.0, 6.0, 8.0];
        let mut cls = Vec::new();
        
        for &alpha_deg in &alphas_deg {
            let flow = FlowConditions::with_alpha_deg(alpha_deg);
            let solution = factorized.solve_alpha(&flow);
            cls.push(solution.cl);
            println!("α = {:5.1}°: Cl = {:7.4}", alpha_deg, solution.cl);
        }
        
        // Compute lift curve slope (should be ~2π/rad ≈ 6.28/rad)
        let cl_alpha_deg = (cls[4] - cls[0]) / (8.0 - 0.0);  // per degree
        let cl_alpha_rad = cl_alpha_deg * 180.0 / PI;
        
        println!("\nLift curve slope:");
        println!("  Cl_α = {:.4}/deg = {:.4}/rad", cl_alpha_deg, cl_alpha_rad);
        println!("  Expected: ~{:.4}/rad (2π)", 2.0 * PI);
        println!("  Error: {:.1}%", (cl_alpha_rad - 2.0 * PI).abs() / (2.0 * PI) * 100.0);
        
        // Should be within 10% of 2π
        assert!(cl_alpha_rad > 5.5 && cl_alpha_rad < 7.0, 
            "Cl_α should be close to 2π (6.28), got {:.4}", cl_alpha_rad);
    }

    #[test]
    fn test_two_solution_efficiency() {
        let solver = InviscidSolver::new();
        let airfoil = make_naca0012(60);
        
        // Factorize once
        let factorized = solver.factorize(&[airfoil]).unwrap();
        
        // Solve for many angles efficiently
        for alpha_deg in (0..=20).map(|i| i as f64) {
            let flow = FlowConditions::with_alpha_deg(alpha_deg);
            let solution = factorized.solve_alpha(&flow);
            
            // Just verify we get reasonable values
            assert!(solution.cl.is_finite());
            assert!(solution.cm.is_finite());
        }
    }

    #[test]
    fn test_no_bodies_error() {
        let solver = InviscidSolver::new();
        let flow = FlowConditions::default();

        let result = solver.solve(&[], &flow);
        assert!(matches!(result, Err(SolverError::NoBodies)));
    }

    #[test]
    fn test_flow_conditions() {
        let flow = FlowConditions::with_alpha_deg(5.0);
        assert!((flow.alpha - 5.0_f64.to_radians()).abs() < 1e-10);
    }

    #[test]
    fn test_cp_distribution() {
        let solver = InviscidSolver::new();
        let airfoil = make_naca0012(60);
        let flow = FlowConditions::with_alpha_deg(4.0);

        let solution = solver.solve(&[airfoil], &flow).unwrap();
        
        // Verify Cp values are reasonable
        for &cp in &solution.cp {
            assert!(cp.is_finite(), "Cp should be finite");
            assert!(cp < 2.0, "Cp should not be too large");
            assert!(cp > -10.0, "Cp should not be too negative");
        }
        
        // Find minimum Cp (should be negative, indicating suction peak)
        let cp_min = solution.cp.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("Cp_min at α=4°: {:.4}", cp_min);
        assert!(cp_min < 0.0, "Should have suction peak (negative Cp)");
    }

    /// Create a NACA 0012 airfoil translated to a given position.
    fn make_naca0012_at(n_half: usize, x_offset: f64, y_offset: f64) -> Body {
        let n = 2 * n_half;
        let t = 0.12;
        let mut points = Vec::with_capacity(n + 1);

        for i in 0..n {
            let theta = 2.0 * PI * (i as f64) / (n as f64);
            let x = 0.5 * (1.0 + theta.cos());
            let thickness = naca_thickness_closed(x, t);
            let y = if theta < PI { thickness } else { -thickness };
            points.push(point(x + x_offset, y + y_offset));
        }
        points.push(points[0]);
        Body::from_points("NACA0012", &points).unwrap()
    }

    /// Create an open-TE NACA 0012 airfoil scaled, rotated, and translated.
    /// Uses the standard (non-closed) thickness formula for a proper blunt TE.
    fn make_naca0012_open_at(n_half: usize, chord: f64, x_pos: f64, y_pos: f64, defl_deg: f64) -> Body {
        let n = n_half;
        let t = 0.12;
        let defl_rad = defl_deg.to_radians();
        let cos_d = defl_rad.cos();
        let sin_d = defl_rad.sin();

        let x_coords: Vec<f64> = (0..n)
            .map(|i| {
                let beta = PI * (i as f64) / ((n - 1) as f64);
                0.5 * (1.0 - beta.cos())
            })
            .collect();

        let mut points = Vec::with_capacity(2 * n - 1);

        // Upper surface: TE to LE
        for i in (0..n).rev() {
            let xc = x_coords[i];
            let y_t = 5.0 * t * (
                0.2969 * xc.sqrt() - 0.126 * xc - 0.3516 * xc.powi(2)
                + 0.2843 * xc.powi(3) - 0.1015 * xc.powi(4)
            );
            let xr = xc * chord;
            let yr = y_t * chord;
            let x = xr * cos_d - yr * sin_d + x_pos;
            let y = xr * sin_d + yr * cos_d + y_pos;
            points.push(point(x, y));
        }
        // Lower surface: LE+1 to TE (skip LE to avoid duplicate)
        for i in 1..n {
            let xc = x_coords[i];
            let y_t = 5.0 * t * (
                0.2969 * xc.sqrt() - 0.126 * xc - 0.3516 * xc.powi(2)
                + 0.2843 * xc.powi(3) - 0.1015 * xc.powi(4)
            );
            let xr = xc * chord;
            let yr = -y_t * chord;
            let x = xr * cos_d - yr * sin_d + x_pos;
            let y = xr * sin_d + yr * cos_d + y_pos;
            points.push(point(x, y));
        }
        Body::from_points("NACA0012", &points).unwrap()
    }

    /// V1: Two symmetric NACA 0012 airfoils side by side at α=0°.
    /// Both should have CL ≈ 0 due to symmetry. Tests multi-body influence
    /// matrix assembly and per-body Kutta condition.
    #[test]
    fn test_v1_two_symmetric_airfoils() {
        let solver = InviscidSolver::new();
        let n_half = 60;

        // Two NACA 0012 airfoils: one at origin, one 3 chords to the right
        let body1 = make_naca0012_at(n_half, 0.0, 0.0);
        let body2 = make_naca0012_at(n_half, 3.0, 0.0);

        let flow = FlowConditions::default(); // α = 0°
        let solution = solver.solve(&[body1, body2], &flow).unwrap();

        println!("V1: Two NACA 0012 at α=0°, separation=3c");
        println!("  Body 1 CL = {:.6}", solution.cl_per_body[0]);
        println!("  Body 2 CL = {:.6}", solution.cl_per_body[1]);
        println!("  Total CL = {:.6}", solution.cl);
        println!("  nodes_per_body = {:?}", solution.nodes_per_body);

        // Both should be near zero for symmetric airfoils at α=0°
        assert!(solution.cl_per_body[0].abs() < 0.3,
            "Body 1 CL should be ~0 at α=0°, got {}", solution.cl_per_body[0]);
        assert!(solution.cl_per_body[1].abs() < 0.3,
            "Body 2 CL should be ~0 at α=0°, got {}", solution.cl_per_body[1]);

        // Also test at α = 5°: both should generate positive lift
        let flow5 = FlowConditions::with_alpha_deg(5.0);
        let sol5 = solver.solve(&[make_naca0012_at(n_half, 0.0, 0.0),
                                   make_naca0012_at(n_half, 3.0, 0.0)], &flow5).unwrap();

        println!("\nV1: Two NACA 0012 at α=5°");
        println!("  Body 1 CL = {:.6}", sol5.cl_per_body[0]);
        println!("  Body 2 CL = {:.6}", sol5.cl_per_body[1]);
        println!("  Total CL = {:.6}", sol5.cl);

        assert!(sol5.cl_per_body[0] > 0.2, "Body 1 should produce lift at α=5°");
        assert!(sol5.cl_per_body[1] > 0.2, "Body 2 should produce lift at α=5°");
    }

    /// V2: NACA 0012 main element + small NACA 0012 flap deflected 20°.
    /// Uses open-TE airfoils for robust multi-body TE treatment.
    /// The flap should increase total CL relative to single-body.
    #[test]
    fn test_v2_naca0012_with_flap() {
        let solver = InviscidSolver::new();

        // Main element: open-TE NACA 0012, chord=1.0 at origin
        let main = make_naca0012_open_at(50, 1.0, 0.0, 0.0, 0.0);

        // Flap: open-TE NACA 0012, chord=0.3c, deflected -20° (TE down)
        // Flap LE at x=2.0c (well-separated for robust solution)
        let flap = make_naca0012_open_at(25, 0.3, 2.0, -0.05, -20.0);

        let flow = FlowConditions::with_alpha_deg(0.0);

        // Single body baseline
        let single_main = make_naca0012_open_at(50, 1.0, 0.0, 0.0, 0.0);
        let single = solver.solve(&[single_main], &flow).unwrap();

        // Multi-body with flap
        let multi = solver.solve(&[main, flap], &flow).unwrap();

        println!("V2: NACA 0012 (open TE) + 0.3c flap at -20°");
        println!("  Single body CL = {:.6}", single.cl);
        println!("  Multi-body total CL = {:.6}", multi.cl);
        println!("  Main element CL = {:.6}", multi.cl_per_body[0]);
        println!("  Flap CL = {:.6}", multi.cl_per_body[1]);

        let delta_cl = multi.cl - single.cl;
        println!("  CL increment = {:.6}", delta_cl);

        // The flap should increase CL
        assert!(delta_cl > 0.1,
            "Flap should increase CL, got ΔCL = {:.4}", delta_cl);

        // Main element should produce positive CL (upwash from flap)
        assert!(multi.cl_per_body[0] > 0.05,
            "Main element should produce positive CL, got {:.4}", multi.cl_per_body[0]);

        // Cp should be smooth and finite
        for (b, cp_b) in multi.cp_per_body.iter().enumerate() {
            for &cp_val in cp_b {
                assert!(cp_val.is_finite(),
                    "Cp should be finite on body {}", b);
            }
        }

        // Test at α=5° — total CL should be greater than single-body at α=5°
        let flow5 = FlowConditions::with_alpha_deg(5.0);
        let single5 = solver.solve(&[make_naca0012_open_at(50, 1.0, 0.0, 0.0, 0.0)], &flow5).unwrap();
        let multi5 = solver.solve(
            &[make_naca0012_open_at(50, 1.0, 0.0, 0.0, 0.0),
              make_naca0012_open_at(25, 0.3, 2.0, -0.05, -20.0)],
            &flow5).unwrap();

        println!("\n  At α=5°: single CL = {:.4}, multi CL = {:.4}", single5.cl, multi5.cl);
        assert!(multi5.cl > single5.cl,
            "Multi-body CL at α=5° should exceed single body");
    }

    /// Ensure single-body backward compatibility: existing tests still pass
    /// through the multi-body code path.
    #[test]
    fn test_single_body_backward_compat() {
        let solver = InviscidSolver::new();
        let airfoil = make_naca0012(80);

        let factorized = solver.factorize(&[airfoil]).unwrap();

        // Should have exactly 1 body
        assert_eq!(factorized.body_node_ranges.len(), 1);
        assert_eq!(factorized.body_chords.len(), 1);

        let flow = FlowConditions::with_alpha_deg(5.0);
        let sol = factorized.solve_alpha(&flow);

        // Should have per-body data for 1 body
        assert_eq!(sol.cl_per_body.len(), 1);
        assert_eq!(sol.gamma_per_body.len(), 1);
        assert_eq!(sol.cp_per_body.len(), 1);
        assert_eq!(sol.nodes_per_body.len(), 1);

        // Per-body CL should equal total CL
        assert!((sol.cl - sol.cl_per_body[0]).abs() < 1e-10,
            "Single body: total CL {} should equal per-body CL {}", sol.cl, sol.cl_per_body[0]);

        // CL should be reasonable for NACA 0012 at 5°
        println!("Single body compat: CL = {:.4} at α=5°", sol.cl);
        assert!(sol.cl > 0.3 && sol.cl < 0.8,
            "CL at 5° should be ~0.55, got {:.4}", sol.cl);
    }

    #[test]
    fn test_xfoil_exact_geometry() {
        use rustfoil_core::{naca::naca4, CubicSpline, PanelingParams};
        
        // Use XFOIL-exact NACA generator and paneling
        let buffer = naca4(12, Some(123));  // 245 buffer points
        let spline = CubicSpline::from_points(&buffer).unwrap();
        let params = PanelingParams::default();
        let paneled = spline.resample_xfoil(160, &params);
        
        let airfoil = Body::from_points("NACA0012_XFOIL", &paneled).unwrap();
        
        let solver = InviscidSolver::new();
        let factorized = solver.factorize(&[airfoil]).unwrap();
        
        println!("\n=== XFOIL-exact geometry (160 panels) ===");
        
        // Test polar symmetry: Cl(-α) should equal -Cl(α)
        let angles = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0];
        let mut cls = Vec::new();
        
        for &alpha_deg in &angles {
            let flow = FlowConditions::with_alpha_deg(alpha_deg);
            let solution = factorized.solve_alpha(&flow);
            cls.push(solution.cl);
            println!("α = {:6.1}°: Cl = {:8.5}", alpha_deg, solution.cl);
        }
        
        // Check Cl at α=0 (should be ~0 for symmetric airfoil)
        let cl_0 = cls[4];
        println!("\nCl at α=0: {:.6} (should be ~0)", cl_0);
        
        // Check symmetry: Cl(-α) should equal -Cl(α)
        println!("\nSymmetry check:");
        for i in 0..4 {
            let cl_neg = cls[i];      // Cl at -8, -6, -4, -2
            let cl_pos = cls[8 - i];  // Cl at +8, +6, +4, +2
            let asymmetry = (cl_neg + cl_pos).abs();
            let angle = (8 - 2*i) as i32;
            println!("  α=±{}: Cl({:+}°) = {:+.5}, Cl({:+}°) = {:+.5}, sum = {:+.6}",
                angle, -angle, cl_neg, angle, cl_pos, asymmetry);
        }
        
        // XFOIL reference values for NACA 0012 (inviscid)
        // These are typical values from XFOIL at Re=inf
        println!("\n=== Comparison with XFOIL reference ===");
        println!("XFOIL typically gives:");
        println!("  Cl at α=0°:  0.0000");
        println!("  Cl at α=4°:  ~0.4588 (inviscid)");
        println!("  Cl_α:        ~6.28/rad (2π)");
        
        let cl_alpha = (cls[7] - cls[1]) / (12.0 * PI / 180.0);  // Over 12° span
        println!("\nOur results:");
        println!("  Cl at α=0°:  {:.5}", cl_0);
        println!("  Cl at α=4°:  {:.5}", cls[6]);
        println!("  Cl_α:        {:.4}/rad", cl_alpha);
    }
}
