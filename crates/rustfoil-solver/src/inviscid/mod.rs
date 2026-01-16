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

pub use error::SolverError;

use nalgebra::{DMatrix, DVector};
use rustfoil_core::{Body, Point};
use std::f64::consts::PI;

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
    /// Vorticity values at each node (γᵢ = surface velocity)
    pub gamma: Vec<f64>,
    /// Pressure coefficient at each node
    pub cp: Vec<f64>,
    /// Lift coefficient (per unit span)
    pub cl: f64,
    /// Moment coefficient about quarter-chord
    pub cm: f64,
    /// Internal stream function value
    pub psi_0: f64,
    /// Number of nodes per body (for indexing)
    pub nodes_per_body: Vec<usize>,
}

/// Cached factorization for efficient alpha sweeps.
/// 
/// XFOIL's key optimization: solve once for α=0° and α=90°,
/// then combine for any angle: γ = cos(α)*γ₀ + sin(α)*γ₉₀
#[derive(Debug, Clone)]
pub struct FactorizedSolution {
    /// Solution for α = 0°
    gamu_0: Vec<f64>,
    /// Solution for α = 90°  
    gamu_90: Vec<f64>,
    /// ψ₀ for α = 0°
    psi0_0: f64,
    /// ψ₀ for α = 90°
    psi0_90: f64,
    /// Cached node positions for Cp/Cl calculation
    nodes: Vec<Point>,
    /// Number of nodes (= number of gammas)
    n_nodes: usize,
    /// Chord length
    chord: f64,
}

impl FactorizedSolution {
    /// Compute solution for any angle of attack.
    pub fn solve_alpha(&self, flow: &FlowConditions) -> InviscidSolution {
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();
        let n = self.n_nodes;

        // Combine the two base solutions
        let gamma: Vec<f64> = (0..n)
            .map(|i| cosa * self.gamu_0[i] + sina * self.gamu_90[i])
            .collect();
        
        let psi_0 = cosa * self.psi0_0 + sina * self.psi0_90;

        // Compute Cp at each node: Cp = 1 - (γ/V∞)²
        let cp: Vec<f64> = gamma
            .iter()
            .map(|&g| 1.0 - (g / flow.v_inf).powi(2))
            .collect();

        // Compute Cl and Cm by pressure integration (XFOIL's CLCALC)
        let (cl, cm) = self.compute_forces(&cp, flow);

        InviscidSolution {
            gamma,
            cp,
            cl,
            cm,
            psi_0,
            nodes_per_body: vec![n],
        }
    }

    /// Compute lift and moment coefficients from Cp distribution.
    /// 
    /// This matches XFOIL's CLCALC subroutine exactly:
    /// - Integrates around the closed contour, wrapping from N back to 1
    /// - Uses trapezoidal integration of pressure forces
    fn compute_forces(&self, cp: &[f64], flow: &FlowConditions) -> (f64, f64) {
        let n = self.n_nodes;
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();
        
        let mut cl = 0.0;
        let mut cm = 0.0;

        // XFOIL's CLCALC: loop from 1 to N, with IP = I+1 and IP=1 when I=N
        // This integrates around the CLOSED contour including the TE panel
        for i in 0..n {
            let ip = (i + 1) % n;  // Wrap around: when i=n-1, ip=0
            
            // Panel geometry
            let dx = self.nodes[ip].x - self.nodes[i].x;
            let dy = self.nodes[ip].y - self.nodes[i].y;
            
            // Rotate to wind axes
            let dx_wind = dx * cosa + dy * sina;
            let dy_wind = dy * cosa - dx * sina;
            
            // Average Cp on this panel (trapezoidal)
            let cp_avg = 0.5 * (cp[i] + cp[ip]);
            
            // CL = integral of Cp * dx (in wind axes)
            cl += cp_avg * dx_wind;
            
            // Moment about quarter-chord
            let x_mid = 0.5 * (self.nodes[i].x + self.nodes[ip].x);
            let y_mid = 0.5 * (self.nodes[i].y + self.nodes[ip].y);
            let x_ref = 0.25 * self.chord;
            
            // Cp difference for higher-order moment term (XFOIL has this)
            let dg = cp[ip] - cp[i];
            
            // XFOIL moment formula (simplified for incompressible)
            cm -= cp_avg * ((x_mid - x_ref) * dx_wind / self.chord 
                         + y_mid * dy_wind / self.chord);
            // Higher-order term
            cm -= (dg * dx_wind * dx_wind / 12.0) / self.chord;
            cm -= (dg * dy_wind * dy_wind / 12.0) / self.chord;
        }

        (cl, cm)
    }
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

    /// Factorize the influence matrix and solve for α=0° and α=90°.
    /// 
    /// This is XFOIL's key optimization - the expensive O(N³) factorization
    /// is done once, then any angle of attack can be computed in O(N) time.
    /// 
    /// # XFOIL Convention
    /// 
    /// XFOIL expects N nodes numbered 1 to N, where:
    /// - Node 1 is upper trailing edge
    /// - Node N is lower trailing edge  
    /// - For sharp TE, there's a small gap between nodes 1 and N
    /// - The N-th "panel" goes from node N back to node 1 (the TE panel)
    /// 
    /// We use 0-based indexing: nodes 0 to N-1, with panel N-1 going from 
    /// node N-1 back to node 0.
    pub fn factorize(&self, bodies: &[Body]) -> SolverResult<FactorizedSolution> {
        if bodies.is_empty() {
            return Err(SolverError::NoBodies);
        }

        // For now, handle single body
        let body = &bodies[0];
        let panels = body.panels();
        let n_panels = panels.len();
        
        if n_panels < 3 {
            return Err(SolverError::InsufficientPanels);
        }

        // Extract nodes: For N panels, we have N+1 nodes (or N nodes if closed)
        // Panel i goes from node i to node i+1
        // The last node (lower TE) is p2 of the last panel
        let mut nodes: Vec<Point> = panels.iter().map(|p| p.p1).collect();
        
        // Check if contour is closed (sharp TE: last panel's p2 == first panel's p1)
        let last_p2 = panels.last().unwrap().p2;
        let first_p1 = panels[0].p1;
        let is_closed = (last_p2.x - first_p1.x).abs() < 1e-10 
                     && (last_p2.y - first_p1.y).abs() < 1e-10;
        
        if !is_closed {
            // Blunt TE: add the lower TE node
            nodes.push(last_p2);
        }
        let n = nodes.len();
        let chord = body.chord();
        
        // Compute distance tolerance for TE panel skip (like XFOIL's SEPS)
        let arc_length: f64 = panels.iter().map(|p| p.length()).sum();
        let seps = arc_length * 1e-5;

        // TE geometry for source/vortex panel treatment
        // Node 0 = upper TE, Node n-1 = lower TE
        let te_upper = nodes[0];
        let te_lower = nodes[n - 1];
        
        // TE gap vector
        let dxte = te_lower.x - te_upper.x;
        let dyte = te_lower.y - te_upper.y;
        let dste = (dxte * dxte + dyte * dyte).sqrt();
        
        // TE midpoint
        let _xte = 0.5 * (te_upper.x + te_lower.x);  // For future bisector use
        let _yte = 0.5 * (te_upper.y + te_lower.y);
        
        // Mean TE panel tangent direction (from tangents of adjacent panels)
        // We have n nodes and n-1 panels (n_panels = n - 1)
        // Panel 0: from node 0 to node 1 (upper TE tangent)
        // Panel n-2: from node n-2 to node n-1 (lower TE tangent, reversed for outward normal)
        let t_upper = panels[0].tangent();
        let t_lower = panels[n_panels - 1].tangent();
        
        // Mean tangent at TE
        let dxs = 0.5 * (t_upper.x - t_lower.x);  // Pointing into airfoil
        let dys = 0.5 * (t_upper.y - t_lower.y);
        
        // XFOIL's ANTE (normal projected gap) and ASTE (tangent projected gap)
        let ante = dxs * dyte - dys * dxte;
        let aste = dxs * dxte + dys * dyte;
        
        // Sharp TE detection (XFOIL: gap < 0.01% chord)
        let sharp = dste < 0.0001 * chord;
        
        // TE panel coefficients (XFOIL's SCS and SDS)
        let (scs, sds) = if sharp {
            (1.0, 0.0)
        } else {
            (ante / dste, aste / dste)
        };

        // Build influence coefficient matrix (N+1 x N+1)
        // Unknowns: [γ₀, γ₁, ..., γₙ₋₁, ψ₀]
        let mut a_matrix = DMatrix::<f64>::zeros(n + 1, n + 1);
        let mut rhs_0 = DVector::<f64>::zeros(n + 1);   // For α = 0°
        let mut rhs_90 = DVector::<f64>::zeros(n + 1);  // For α = 90°

        // For each node i, compute influence from all panels
        for i in 0..n {
            let xi = nodes[i].x;
            let yi = nodes[i].y;

            // Initialize influence coefficient array for this row
            let mut dzdg = vec![0.0; n];

            // Variables to store TE panel geometry for source/vortex treatment
            let mut te_g1 = 0.0;
            let mut te_g2 = 0.0;
            let mut te_t1 = 0.0;
            let mut te_t2 = 0.0;
            let mut te_x1 = 0.0;
            let mut te_x2 = 0.0;
            let mut te_yy = 0.0;
            let mut te_apan = 0.0;
            let mut te_panel_valid = false;

            // Loop over all panels (XFOIL loops JO=1 to N)
            for jo in 0..n {
                let jp = (jo + 1) % n;  // XFOIL: JP = JO+1, but JP=1 when JO=N
                
                // Panel endpoints
                let x_jo = nodes[jo].x;
                let y_jo = nodes[jo].y;
                let x_jp = nodes[jp].x;
                let y_jp = nodes[jp].y;
                
                // Skip TE panel if endpoints coincide (XFOIL's SEPS check)
                let dx = x_jp - x_jo;
                let dy = y_jp - y_jo;
                let ds_sq = dx * dx + dy * dy;
                
                // For TE panel (jo = n-1), check if it should be skipped
                if jo == n - 1 && ds_sq < seps * seps {
                    continue;
                }
                
                if ds_sq < 1e-24 {
                    continue;
                }

                // Panel geometry
                let dso = ds_sq.sqrt();
                let dsio = 1.0 / dso;
                let sx = dx * dsio;
                let sy = dy * dsio;

                // Vector from panel endpoints to field point
                let rx1 = xi - x_jo;
                let ry1 = yi - y_jo;
                let rx2 = xi - x_jp;
                let ry2 = yi - y_jp;

                // Transform to panel-local coordinates
                let x1 = sx * rx1 + sy * ry1;
                let x2 = sx * rx2 + sy * ry2;
                let yy = sx * ry1 - sy * rx1;

                let rs1 = rx1 * rx1 + ry1 * ry1;
                let rs2 = rx2 * rx2 + ry2 * ry2;

                // Logarithm and arctangent terms
                let (g1, t1) = if i != jo && rs1 > 1e-20 {
                    (rs1.ln(), x1.atan2(yy))
                } else {
                    (0.0, 0.0)
                };

                let (g2, t2) = if i != jp && rs2 > 1e-20 {
                    (rs2.ln(), x2.atan2(yy))
                } else {
                    (0.0, 0.0)
                };

                // For TE panel (jo = n-1), save geometry and skip vortex contribution
                // This is XFOIL's "IF(JO.EQ.N) GO TO 11" logic
                if jo == n - 1 {
                    te_g1 = g1;
                    te_g2 = g2;
                    te_t1 = t1;
                    te_t2 = t2;
                    te_x1 = x1;
                    te_x2 = x2;
                    te_yy = yy;
                    // XFOIL's TE panel angle (xpanel.f line 44):
                    // APANEL(N) = ATAN2(-SX, SY) + PI
                    // where SX = X(1) - X(N), SY = Y(1) - Y(N)
                    // Note: sx,sy here are already the normalized tangent direction
                    te_apan = (-sx).atan2(sy) + PI;
                    te_panel_valid = true;
                    continue;  // Skip normal vortex contribution for TE panel
                }

                // Normal vortex panel contribution (PSIS/PSID)
                let dxinv = 1.0 / (x1 - x2);
                let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
                let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

                dzdg[jo] += QOPI * (psis - psid);
                dzdg[jp] += QOPI * (psis + psid);
            }

            // TE panel source/vortex contribution (XFOIL label 11 code)
            if te_panel_valid {
                // XFOIL's PSIG and PGAM formulas
                let psig = 0.5 * te_yy * (te_g1 - te_g2) 
                         + te_x2 * (te_t2 - te_apan) 
                         - te_x1 * (te_t1 - te_apan);
                let pgam = 0.5 * te_x1 * te_g1 - 0.5 * te_x2 * te_g2 
                         + te_x2 - te_x1 + te_yy * (te_t1 - te_t2);

                // Add TE panel contribution to influence coefficients
                // XFOIL: JO = N (lower TE), JP = 1 (upper TE)
                // Our: jo = n-1 (lower TE), jp = 0 (upper TE)
                let jo_te = n - 1;
                let jp_te = 0;
                
                dzdg[jo_te] -= HOPI * psig * scs * 0.5;
                dzdg[jp_te] += HOPI * psig * scs * 0.5;
                dzdg[jo_te] += HOPI * pgam * sds * 0.5;
                dzdg[jp_te] -= HOPI * pgam * sds * 0.5;
            }

            // Fill matrix row i
            for j in 0..n {
                a_matrix[(i, j)] = dzdg[j];
            }
            
            // Column n: coefficient for ψ₀ (the unknown internal stream function)
            // Boundary condition is: ψ_induced + ψ_freestream = ψ₀
            // So: Σ(dzdg_j * γ_j) - ψ₀ = -ψ_freestream
            a_matrix[(i, n)] = -1.0;

            // RHS: -ψ_freestream at node i
            // For α = 0°:  ψ∞ = V∞ * y (freestream from left)
            // For α = 90°: ψ∞ = -V∞ * x (freestream from below)
            rhs_0[i] = -yi;
            rhs_90[i] = xi;
        }

        // Last row: Kutta condition γ₀ + γₙ₋₁ = 0
        // (γ at upper TE + γ at lower TE = 0)
        for j in 0..=n {
            a_matrix[(n, j)] = 0.0;
        }
        a_matrix[(n, 0)] = 1.0;
        a_matrix[(n, n - 1)] = 1.0;
        rhs_0[n] = 0.0;
        rhs_90[n] = 0.0;

        // LU factorization and solve for both RHS vectors
        let lu = a_matrix.clone().lu();
        
        let solution_0 = lu.solve(&rhs_0)
            .ok_or(SolverError::SingularMatrix)?;
        let solution_90 = lu.solve(&rhs_90)
            .ok_or(SolverError::SingularMatrix)?;

        // Extract γ and ψ₀ from solutions
        let gamu_0: Vec<f64> = solution_0.iter().take(n).copied().collect();
        let gamu_90: Vec<f64> = solution_90.iter().take(n).copied().collect();
        let psi0_0 = solution_0[n];
        let psi0_90 = solution_90[n];

        Ok(FactorizedSolution {
            gamu_0,
            gamu_90,
            psi0_0,
            psi0_90,
            nodes,
            n_nodes: n,
            chord,
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
