//! Matrix assembly and solution (XFOIL's GGCALC).
//!
//! This module builds the (N+1)×(N+1) linear system for the inviscid flow
//! and solves for the two base solutions at α=0° and α=90°.
//!
//! # System Structure
//!
//! For N nodes, we have N+1 unknowns: [γ₀, γ₁, ..., γₙ₋₁, ψ₀]
//!
//! ```text
//! | A₀₀  A₀₁  ... A₀,ₙ₋₁  -1 | | γ₀    |   | -ψ∞(0)   |
//! | A₁₀  A₁₁  ... A₁,ₙ₋₁  -1 | | γ₁    |   | -ψ∞(1)   |
//! |  :    :   ...   :      : | |  :    | = |    :     |
//! | Aₙ₋₁,₀...    Aₙ₋₁,ₙ₋₁ -1 | | γₙ₋₁  |   | -ψ∞(n-1) |
//! |  1    0   ...   1      0 | | ψ₀    |   |    0     |
//! ```
//!
//! Where:
//! - A[i][j] = ∂ψᵢ/∂γⱼ (from PSILIN)
//! - Last row: Kutta condition γ₀ + γₙ₋₁ = 0
//!
//! # XFOIL Reference
//!
//! - `xpanel.f`: GGCALC subroutine
//! - `xsolve.f`: LUDCMP, BAKSUB

use crate::geometry::AirfoilGeometry;
use crate::influence::{build_source_influence_matrix, psilin, psilin_with_dqdm};
use crate::solution::{FlowConditions, InviscidSolution};
use crate::{InviscidError, Result};
use nalgebra::{DMatrix, DVector, LU};
use rustfoil_bl::{add_event, is_debug_active, DebugEvent};

/// Factorized system ready for efficient alpha sweeps.
///
/// XFOIL's key optimization: solve the expensive O(N³) factorization once,
/// then any angle of attack can be computed in O(N) time.
#[derive(Debug, Clone)]
pub struct FactorizedSystem {
    /// Solution γ for α = 0° (GAMU(:,1) in XFOIL)
    pub gamu_0: Vec<f64>,
    /// Solution γ for α = 90° (GAMU(:,2) in XFOIL)
    pub gamu_90: Vec<f64>,
    /// Internal stream function ψ₀ for α = 0°
    pub psi0_0: f64,
    /// Internal stream function ψ₀ for α = 90°
    pub psi0_90: f64,
    /// Cached geometry for solution computation
    pub(crate) geom: AirfoilGeometry,
    /// LU factorization of the influence matrix (AIJ)
    pub(crate) lu: LU<f64, nalgebra::Dynamic, nalgebra::Dynamic>,
}

impl FactorizedSystem {
    /// Compute solution for any angle of attack.
    ///
    /// Uses linear combination: γ = cos(α)*γ₀ + sin(α)*γ₉₀
    pub fn solve_alpha(&self, flow: &FlowConditions) -> InviscidSolution {
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();
        let n = self.geom.n;

        // Combine the two base solutions (XFOIL's SPECAL)
        let gamma: Vec<f64> = (0..n)
            .map(|i| cosa * self.gamu_0[i] + sina * self.gamu_90[i])
            .collect();

        // Gamma derivative w.r.t. alpha (XFOIL's QINV_A)
        let gamma_a: Vec<f64> = (0..n)
            .map(|i| -sina * self.gamu_0[i] + cosa * self.gamu_90[i])
            .collect();

        let psi_0 = cosa * self.psi0_0 + sina * self.psi0_90;

        // Compute Cp at each node: Cp = 1 - (γ/V∞)²
        let cp: Vec<f64> = gamma
            .iter()
            .map(|&g| 1.0 - (g / flow.v_inf).powi(2))
            .collect();

        // Compute Cl and Cm by pressure integration (XFOIL's CLCALC)
        let (cl, cm) = self.compute_forces(&cp, flow);

        // Debug output: emit full inviscid solution
        if is_debug_active() {
            // qinv is the same as gamma for inviscid panel method
            add_event(DebugEvent::full_inviscid(
                gamma.clone(),
                gamma.clone(), // qinv = gamma for inviscid
                cp.clone(),
                cl,
                flow.alpha,
            ));
        }

        InviscidSolution {
            gamma,
            gamma_a,
            cp,
            cl,
            cm,
            psi_0,
            n,
        }
    }

    /// Compute the surface inviscid-edge-velocity basis arrays QINVU(:,1:2).
    ///
    /// XFOIL stores these separately from the GAMU circulation basis and uses
    /// them when seeding the wake-edge velocity in QWCALC.
    pub fn surface_qinvu_basis(&self) -> (Vec<f64>, Vec<f64>) {
        // On the airfoil surface, XFOIL's tangential edge velocity matches the
        // node circulation: QINVU(:,1:2) == GAMU(:,1:2). Wake points are handled
        // separately by QWCALC.
        (self.gamu_0.clone(), self.gamu_90.clone())
    }

    /// Compute lift and moment coefficients from Cp distribution.
    ///
    /// This matches XFOIL's CLCALC subroutine exactly:
    /// - Integrates around the closed contour, wrapping from N back to 1
    /// - Uses trapezoidal integration of pressure forces
    fn compute_forces(&self, cp: &[f64], flow: &FlowConditions) -> (f64, f64) {
        let n = self.geom.n;
        let cosa = flow.alpha.cos();
        let sina = flow.alpha.sin();

        let mut cl = 0.0;
        let mut cm = 0.0;

        let x_ref = 0.25 * self.geom.chord;

        // XFOIL's CLCALC: loop from 0 to N-1, with IP = I+1 and IP=0 when I=N-1
        // This integrates around the CLOSED contour including the TE panel
        for i in 0..n {
            let ip = (i + 1) % n;

            // Panel geometry
            let dx = self.geom.x[ip] - self.geom.x[i];
            let dy = self.geom.y[ip] - self.geom.y[i];

            // Rotate to wind axes
            let dx_wind = dx * cosa + dy * sina;
            let dy_wind = dy * cosa - dx * sina;

            // Average Cp on this panel (trapezoidal)
            let cp_avg = 0.5 * (cp[i] + cp[ip]);
            let dg = cp[ip] - cp[i];

            // CL = integral of Cp * dx (in wind axes)
            cl += cp_avg * dx_wind;

            // Moment about quarter-chord
            let x_mid = 0.5 * (self.geom.x[i] + self.geom.x[ip]);
            let y_mid = 0.5 * (self.geom.y[i] + self.geom.y[ip]);

            // XFOIL moment formula
            let ax = (x_mid - x_ref) * cosa + y_mid * sina;
            let ay = y_mid * cosa - (x_mid - x_ref) * sina;

            cm -= cp_avg * (ax * dx_wind / self.geom.chord + ay * dy_wind / self.geom.chord);
            cm -= dg * dx_wind * dx_wind / (12.0 * self.geom.chord);
            cm -= dg * dy_wind * dy_wind / (12.0 * self.geom.chord);
        }

        (cl, cm)
    }

    /// Get a reference to the geometry.
    pub fn geometry(&self) -> &AirfoilGeometry {
        &self.geom
    }

    /// Back-substitute a vector through the influence matrix.
    ///
    /// Solves: AIJ × x = b for x.
    pub fn back_substitute(&self, rhs: &[f64]) -> Result<Vec<f64>> {
        let rhs_vec = DVector::from_column_slice(rhs);
        let solution = self
            .lu
            .solve(&rhs_vec)
            .ok_or(InviscidError::SingularMatrix)?;
        Ok(solution.iter().copied().collect())
    }

    /// Build the mass defect influence matrix DIJ = AIJ⁻¹ × BIJ.
    pub fn build_dij(&self) -> Result<DMatrix<f64>> {
        let n = self.geom.n;
        if n == 0 {
            return Ok(DMatrix::zeros(0, 0));
        }

        let bij = build_source_influence_matrix(&self.geom);
        let mut dij = DMatrix::zeros(n, n);

        for j in 0..n {
            let rhs = bij.column(j).clone_owned();
            let solution = self
                .lu
                .solve(&rhs)
                .ok_or(InviscidError::SingularMatrix)?;

            for i in 0..n {
                dij[(i, j)] = solution[i];
            }
        }
        
        Ok(dij)
    }

    /// Build the mass defect influence matrix including wake panels.
    ///
    /// This extends the airfoil DIJ with wake panels so panel indices N..N+NW-1
    /// are valid for VI coupling (matching XFOIL's QDCALC behavior).
    ///
    /// # Algorithm (XFOIL's QDCALC)
    ///
    /// 1. Airfoil-airfoil block: Use existing LU factorization (DIJ = AIJ^{-1} × BIJ)
    /// 2. Airfoil-wake block: Call PSWLIN for each airfoil node, back-substitute
    /// 3. Wake-airfoil block: Call PSILIN for each wake node  
    /// 4. Wake-wake block: Call PSWLIN for each wake node
    pub fn build_dij_with_wake(&self, wake_x: &[f64], wake_y: &[f64]) -> Result<DMatrix<f64>> {
        let nw = wake_x.len();
        let (wake_apanel, wake_nx, wake_ny) = compute_wake_geometry_from_coordinates(wake_x, wake_y);
        debug_assert_eq!(wake_apanel.len(), nw);
        debug_assert_eq!(wake_nx.len(), nw);
        debug_assert_eq!(wake_ny.len(), nw);
        self.build_dij_with_wake_state(wake_x, wake_y, &wake_nx, &wake_ny, &wake_apanel)
    }

    /// Build the mass defect influence matrix using explicit wake geometry state.
    ///
    /// This matches XFOIL's ownership model where `XYWAKE` defines wake `X/Y/NX/NY/APANEL`
    /// and `QDCALC` consumes that stored state directly instead of reconstructing it from
    /// the wake coordinates.
    pub fn build_dij_with_wake_state(
        &self,
        wake_x: &[f64],
        wake_y: &[f64],
        wake_nx: &[f64],
        wake_ny: &[f64],
        wake_apanel: &[f64],
    ) -> Result<DMatrix<f64>> {
        assert_eq!(
            wake_x.len(),
            wake_y.len(),
            "wake_x and wake_y must have the same length"
        );
        assert_eq!(
            wake_x.len(),
            wake_nx.len(),
            "wake_x and wake_nx must have the same length"
        );
        assert_eq!(
            wake_x.len(),
            wake_ny.len(),
            "wake_x and wake_ny must have the same length"
        );
        assert_eq!(
            wake_x.len(),
            wake_apanel.len(),
            "wake_x and wake_apanel must have the same length"
        );

        let n = self.geom.n;
        let nw = wake_x.len();
        if nw == 0 {
            return self.build_dij();
        }

        let n_total = n + nw;
        let dij_airfoil = self.build_dij()?;
        let mut dij = DMatrix::zeros(n_total, n_total);

        // Copy existing airfoil-airfoil block (most accurate from LU factorization).
        for i in 0..n {
            for j in 0..n {
                dij[(i, j)] = dij_airfoil[(i, j)];
            }
        }

        // Build combined coordinate arrays.
        let mut x_all = Vec::with_capacity(n_total);
        let mut y_all = Vec::with_capacity(n_total);
        x_all.extend_from_slice(&self.geom.x);
        y_all.extend_from_slice(&self.geom.y);
        x_all.extend_from_slice(wake_x);
        y_all.extend_from_slice(wake_y);

        // Use airfoil geometry from SETXY/GGCALC and wake geometry from XYWAKE directly.
        let (apanel, nx, ny) =
            combine_geometry_state(&self.geom, wake_nx, wake_ny, wake_apanel);

        // Step 2: Airfoil-wake block (QDCALC lines 1187-1216).
        // For each wake panel J, compute dPsi/dm at all airfoil nodes, then back-substitute.
        // This follows XFOIL's structure: build BIJ column-by-column for wake panels.
        
        // First, compute DZDM for all airfoil nodes (store as rows temporarily).
        let mut dzdm_matrix = vec![vec![0.0; nw]; n];
        for i in 0..n {
            let xi = self.geom.x[i];
            let yi = self.geom.y[i];
            let nxi = self.geom.nx[i];
            let nyi = self.geom.ny[i];

            let result = pswlin(&x_all, &y_all, &apanel, n, i, xi, yi, nxi, nyi);
            
            // Store dPsi/dm for wake panels in this row.
            for j in 0..nw {
                dzdm_matrix[i][j] = result.dzdm[n + j];
            }
        }
        
        // Now process each wake panel column through back-substitution.
        for j in 0..nw {
            // Build BIJ column for wake panel j: BIJ(I, N+j) = -DZDM_i(N+j).
            let mut bij_column = vec![0.0; n + 1];
            for i in 0..n {
                bij_column[i] = -dzdm_matrix[i][j];
            }
            if self.geom.sharp {
                // XFOIL's sharp-TE bisector row has no direct wake-source influence.
                bij_column[n - 1] = 0.0;
            }
            // Last entry (Kutta condition row) is 0.
            bij_column[n] = 0.0;

            // Back-substitute to get gamma response: AIJ × X = BIJ(:,j).
            let dij_col = self.back_substitute(&bij_column)?;
            
            // Store result in DIJ(:, N+j).
            for i in 0..n {
                dij[(i, n + j)] = dij_col[i];
            }
        }

        // Step 3 & 4: Wake-airfoil and wake-wake blocks (QDCALC lines 1230-1257).
        // For each wake node I, compute influences from all panels.
        // Also store CIJ = DQDG for indirect coupling (Fortran lines 1237-1239).
        let mut cij = vec![vec![0.0; n]; nw]; // CIJ(IW, K) = dQtan/dGam
        for i in 0..nw {
            let i_global = n + i;
            let xi = x_all[i_global];
            let yi = y_all[i_global];
            let nxi = nx[i_global];
            let nyi = ny[i_global];

            // Wake-airfoil: compute source and vortex influence on tangential velocity.
            // XFOIL's PSILIN with SIGLIN=.TRUE. computes DQDM (dQtan/dSigma) and DQDG (dQtan/dGam).
            let psilin_result = psilin_with_dqdm(&self.geom, i_global, xi, yi, nxi, nyi);

            // DIJ(I,J) = DQDM(J) for wake-airfoil block (Fortran lines 1241-1243)
            for j in 0..n {
                dij[(i_global, j)] = psilin_result.dqdm[j];
            }

            // CIJ(IW,J) = DQDG(J) for indirect coupling (Fortran lines 1237-1239)
            for j in 0..n {
                cij[i][j] = psilin_result.dqdg[j];
            }

            // Debug: dump DQDM values for first 2 wake points (matches XFOIL's DBGDQDM_PSILIN)
            if i < 2 && is_debug_active() {
                let dqdm_sample: Vec<f64> = psilin_result.dqdm.iter().take(20).copied().collect();
                add_event(DebugEvent::dqdm_psilin(
                    i_global,      // wake_idx (global index)
                    i + 1,         // iw (1-based to match XFOIL)
                    n,             // n_airfoil
                    xi, yi, nxi, nyi,
                    dqdm_sample,
                ));

                // On first wake point, also dump panel geometry for comparison
                if i == 0 {
                    let panel_x: Vec<f64> = self.geom.x.iter().take(20).copied().collect();
                    let panel_y: Vec<f64> = self.geom.y.iter().take(20).copied().collect();
                    add_event(DebugEvent::panel_geometry(n, panel_x, panel_y));
                }
            }

            // Wake-wake: call PSWLIN for wake panels (tangential velocity influence).
            // DIJ(I,J) = DQDM(J) for wake-wake block (Fortran lines 1251-1256)
            let result = pswlin(&x_all, &y_all, &apanel, n, i_global, xi, yi, nxi, nyi);
            for j in 0..nw {
                dij[(i_global, n + j)] = result.dqdm[n + j];
            }
        }

        // Step 5: CIJ*DIJ indirect coupling (QDCALC lines 1260-1282).
        // Add effect of all sources on airfoil vorticity which affects wake Qtan.
        // DIJ(I,J) += sum_K CIJ(IW,K) * DIJ(K,J)  for K=1..N
        // This is the chain rule: wake velocity depends on vorticity (CIJ),
        // which depends on sources (DIJ airfoil rows).
        for i in 0..nw {
            let i_global = n + i;
            for j in 0..n_total {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += cij[i][k] * dij[(k, j)];
                }
                dij[(i_global, j)] += sum;
            }
        }

        // Step 6: First wake point = trailing edge (QDCALC line 1286).
        // Make sure first wake point has same velocity influence as trailing edge.
        // DIJ(N+1, J) = DIJ(N, J) for all J
        if nw > 0 {
            for j in 0..n_total {
                dij[(n, j)] = dij[(n - 1, j)];
            }
        }

        // Debug output.
        if is_debug_active() {
            let dij_flat: Vec<f64> = dij.iter().copied().collect();
            let diagonal_sample: Vec<f64> = (0..n_total.min(20))
                .map(|i| dij[(i, i)])
                .collect();
            let row1_sample: Vec<f64> = (0..n_total.min(20))
                .map(|j| dij[(0, j)])
                .collect();
            
            add_event(DebugEvent::full_dij(
                n_total,
                dij_flat,
                diagonal_sample,
                row1_sample,
            ));
        }

        Ok(dij)
    }

    /// Build the DIJ matrix using a default XFOIL-style wake.
    pub fn build_dij_with_default_wake(&self) -> Result<DMatrix<f64>> {
        let (wake_x, wake_y) = self.build_wake_coordinates(DEFAULT_WAKLEN);
        self.build_dij_with_wake(&wake_x, &wake_y)
    }

    /// Generate wake coordinates using XFOIL-style spacing.
    pub fn build_wake_coordinates(&self, wake_length: f64) -> (Vec<f64>, Vec<f64>) {
        let n = self.geom.n;
        if n < 2 || wake_length <= 0.0 {
            return (Vec::new(), Vec::new());
        }

        let n_wake = compute_wake_count(n, wake_length);
        if n_wake == 0 {
            return (Vec::new(), Vec::new());
        }

        let ds1 = 0.5 * ((self.geom.s[1] - self.geom.s[0]) + (self.geom.s[n - 1] - self.geom.s[n - 2]));
        let snew = setexp(ds1, wake_length * self.geom.chord, n_wake);

        // XFOIL-style wake normal from TE tangents (XYWAKE).
        let sx = 0.5 * (self.geom.yp[n - 1] - self.geom.yp[0]);
        let sy = 0.5 * (self.geom.xp[0] - self.geom.xp[n - 1]);
        let smod = (sx * sx + sy * sy).sqrt().max(1e-12);
        let nx = sx / smod;
        let ny = sy / smod;

        // Wake tangent direction.
        let tx = -ny;
        let ty = nx;

        let offset = 1.0e-4 * self.geom.chord;
        let mut wake_x = Vec::with_capacity(n_wake);
        let mut wake_y = Vec::with_capacity(n_wake);

        for s in snew {
            let dist = offset + s;
            wake_x.push(self.geom.xte + dist * tx);
            wake_y.push(self.geom.yte + dist * ty);
        }

        (wake_x, wake_y)
    }
}

const DEFAULT_WAKLEN: f64 = 1.0;

pub fn compute_wake_count(n_panels: usize, wake_length: f64) -> usize {
    let wake_blocks = (wake_length.floor() as usize).saturating_mul(10);
    n_panels / 12 + wake_blocks
}

pub fn setexp(ds1: f64, smax: f64, n_points: usize) -> Vec<f64> {
    if n_points == 0 {
        return Vec::new();
    }
    if n_points == 1 {
        return vec![0.0];
    }

    let sigma = smax / ds1.max(1e-12);
    let nex = n_points - 1;
    let rnex = nex as f64;
    let rni = 1.0 / rnex;

    let aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    let bbb = rnex * (rnex - 1.0) / 2.0;
    let ccc = rnex - sigma;

    let mut disc = bbb * bbb - 4.0 * aaa * ccc;
    if disc < 0.0 {
        disc = 0.0;
    }

    let mut ratio = if nex <= 1 {
        1.0
    } else if nex == 2 {
        -ccc / bbb + 1.0
    } else {
        (-bbb + disc.sqrt()) / (2.0 * aaa) + 1.0
    };

    if (ratio - 1.0).abs() > 1e-12 {
        for _ in 0..100 {
            let sigman = (ratio.powi(nex as i32) - 1.0) / (ratio - 1.0);
            let res = sigman.powf(rni) - sigma.powf(rni);
            let dresdr = rni
                * sigman.powf(rni)
                * (rnex * ratio.powi(nex.saturating_sub(1) as i32) - sigman)
                / (ratio.powi(nex as i32) - 1.0);
            let dratio = -res / dresdr;
            ratio += dratio;
            if dratio.abs() < 1.0e-5 {
                break;
            }
        }
    }

    let mut s = vec![0.0; n_points];
    let mut ds = ds1;
    for i in 1..n_points {
        s[i] = s[i - 1] + ds;
        ds *= ratio;
    }
    s
}

fn combine_geometry_state(
    geom: &AirfoilGeometry,
    wake_nx: &[f64],
    wake_ny: &[f64],
    wake_apanel: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n_airfoil = geom.n;
    let n_total = n_airfoil + wake_nx.len();
    let mut apanel = vec![0.0; n_total];
    let mut nx = vec![0.0; n_total];
    let mut ny = vec![0.0; n_total];

    apanel[..n_airfoil].copy_from_slice(&geom.apanel);
    nx[..n_airfoil].copy_from_slice(&geom.nx);
    ny[..n_airfoil].copy_from_slice(&geom.ny);
    apanel[n_airfoil..].copy_from_slice(wake_apanel);
    nx[n_airfoil..].copy_from_slice(wake_nx);
    ny[n_airfoil..].copy_from_slice(wake_ny);

    (apanel, nx, ny)
}

/// Compute wake panel geometry from coordinates only.
///
/// This is used by generic DIJ builders that do not have access to the explicit `XYWAKE`
/// state. The faithful `rustfoil-xfoil` path calls `build_dij_with_wake_state()` instead.
fn compute_wake_geometry_from_coordinates(
    wake_x: &[f64],
    wake_y: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    use std::f64::consts::PI;
    let nw = wake_x.len();
    if nw == 0 {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    let mut apanel = vec![0.0; nw];
    let mut nx = vec![0.0; nw];
    let mut ny = vec![0.0; nw];

    for i in 0..nw {
        let (tx, ty) = if i < nw - 1 {
            (wake_x[i + 1] - wake_x[i], wake_y[i + 1] - wake_y[i])
        } else if nw >= 2 {
            (wake_x[i] - wake_x[i - 1], wake_y[i] - wake_y[i - 1])
        } else {
            (1.0, 0.0)
        };
        let ds = (tx * tx + ty * ty).sqrt().max(1.0e-12);
        nx[i] = ty / ds;
        ny[i] = -tx / ds;
        if i < nw - 1 {
            apanel[i] = ty.atan2(tx);
        } else if nw == 1 {
            apanel[i] = PI;
        } else {
            apanel[i] = apanel[i - 1];
        }
    }

    (apanel, nx, ny)
}

/// Result from PSWLIN computation
struct PswlinResult {
    /// dPsi/dm - streamfunction derivatives
    dzdm: Vec<f64>,
    /// dQtan/dm - tangential velocity derivatives
    dqdm: Vec<f64>,
}

/// XFOIL's PSWLIN: compute wake panel influences on a point.
///
/// Calculates dPsi/dm (DZDM) and dQtan/dm (DQDM) at point (xi, yi)
/// due to wake sources using quadratic vorticity distribution.
///
/// # Arguments
///
/// * `x, y` - Combined airfoil + wake coordinates
/// * `apanel` - Panel angles for all panels
/// * `n_airfoil` - Number of airfoil nodes (wake starts at index n_airfoil)
/// * `i_obs` - Observer point index (for SGN computation)
/// * `xi, yi` - Observer point coordinates
/// * `nxi, nyi` - Observer point normal
///
/// # Returns
///
/// PswlinResult containing both DZDM and DQDM arrays (size = x.len()).
fn pswlin(
    x: &[f64],
    y: &[f64],
    apanel: &[f64],
    n_airfoil: usize,
    i_obs: usize,
    xi: f64,
    yi: f64,
    nxi: f64,
    nyi: f64,
) -> PswlinResult {
    use std::f64::consts::PI;
    let n_total = x.len();
    let nw = n_total - n_airfoil;
    if nw < 2 {
        return PswlinResult {
            dzdm: vec![0.0; n_total],
            dqdm: vec![0.0; n_total],
        };
    }

    let qopi = 0.25 / PI; // 1/(4π)
    let mut dzdm = vec![0.0; n_total];
    let mut dqdm = vec![0.0; n_total];

    // Loop over wake panels (XFOIL PSWLIN lines 829-980).
    // Note: XFOIL loops from N+1 to N+NW-1 (0-indexed: n_airfoil to n_total-2).
    for jo in n_airfoil..n_total - 1 {
        let jp = jo + 1;
        
        // Handle neighboring panel indices with boundary conditions.
        let jm = if jo == n_airfoil { jo } else { jo - 1 };
        let jq = if jo == n_total - 2 { jp } else { jp + 1 };

        // Panel length DSO.
        let dso = ((x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2)).sqrt();
        let dsio = 1.0 / dso.max(1e-12);

        let apan = apanel[jo];

        // Vector from panel nodes to observer point.
        let rx1 = xi - x[jo];
        let ry1 = yi - y[jo];
        let rx2 = xi - x[jp];
        let ry2 = yi - y[jp];

        // Panel coordinate system (SX, SY = tangent).
        let sx = (x[jp] - x[jo]) * dsio;
        let sy = (y[jp] - y[jo]) * dsio;

        // Transform to panel coordinates.
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // Sign convention for branch cuts (XFOIL lines 861-865).
        let sgn = if i_obs >= n_airfoil && i_obs < n_total {
            1.0 // Wake point: no branch cut flipping
        } else {
            yy.signum() // Airfoil point: depends on which side
        };

        // Logarithm and angle terms (XFOIL lines 867-881).
        let (g1, t1) = if i_obs != jo && rs1 > 0.0 {
            let g = rs1.ln();
            let t = (sgn * x1).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;
            (g, t)
        } else {
            (0.0, 0.0)
        };

        let (g2, t2) = if i_obs != jp && rs2 > 0.0 {
            let g = rs2.ln();
            let t = (sgn * x2).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;
            (g, t)
        } else {
            (0.0, 0.0)
        };

        // Transform normal to panel coordinates.
        let x1i = sx * nxi + sy * nyi;
        let x2i = x1i; // Same for both endpoints
        let yyi = sx * nyi - sy * nxi;

        // Midpoint quantities (XFOIL lines 888-891).
        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = rs0.ln();
        let t0 = (sgn * x0).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;

        // ============ First half-panel (1-0) - XFOIL lines 893-934 ============
        {
            let dxinv = 1.0 / (x1 - x0).max(1e-12);
            
            // Panel integrals for streamfunction (PSUM, PDIF).
            let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan)
                + (x0 - x1) * yy) * dxinv;

            // Derivatives for tangential velocity.
            let psx1 = -(t1 - apan);
            let psx0 = t0 - apan;
            let psyy = 0.5 * (g1 - g0);

            let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
            let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
            let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

            // Source strength derivatives (XFOIL lines 907-916).
            let dsm = ((x[jp] - x[jm]).powi(2) + (y[jp] - y[jm]).powi(2)).sqrt();
            let dsim = 1.0 / dsm.max(1e-12);

            // Accumulate dPsi/dm (DZDM) - for streamfunction (XFOIL lines 921-924).
            dzdm[jm] += qopi * (-psum * dsim + pdif * dsim);
            dzdm[jo] += qopi * (-psum * dsio - pdif * dsio);
            dzdm[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

            // dQtan/dm contributions (for tangential velocity, XFOIL lines 927-934).
            let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;

            // Accumulate dQtan/dm (DQDM) - for tangential velocity.
            dqdm[jm] += qopi * (-psni * dsim + pdni * dsim);
            dqdm[jo] += qopi * (-psni * dsio - pdni * dsio);
            dqdm[jp] += qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));
        }

        // ============ Second half-panel (0-2) - XFOIL lines 937-978 ============
        {
            let dxinv = 1.0 / (x0 - x2).max(1e-12);
            
            let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan)
                + (x2 - x0) * yy) * dxinv;

            let psx0 = -(t0 - apan);
            let psx2 = t2 - apan;
            let psyy = 0.5 * (g0 - g2);

            let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
            let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
            let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

            let dsp = ((x[jq] - x[jo]).powi(2) + (y[jq] - y[jo]).powi(2)).sqrt();
            let dsip = 1.0 / dsp.max(1e-12);

            // Accumulate dPsi/dm (DZDM) - for streamfunction (XFOIL lines 965-968).
            dzdm[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            dzdm[jp] += qopi * (psum * dsio - pdif * dsio);
            dzdm[jq] += qopi * (psum * dsip + pdif * dsip);

            // dQtan/dm contributions (XFOIL lines 971-978).
            let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
            let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;

            dqdm[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
            dqdm[jp] += qopi * (psni * dsio - pdni * dsio);
            dqdm[jq] += qopi * (psni * dsip + pdni * dsip);
        }
    }

    PswlinResult { dzdm, dqdm }
}

/// Build and factorize the influence coefficient system.
///
/// This is XFOIL's GGCALC subroutine.
///
/// # Arguments
///
/// * `geom` - Airfoil geometry (from `AirfoilGeometry::from_points`)
///
/// # Returns
///
/// A `FactorizedSystem` ready for alpha sweeps.
pub fn build_and_factorize(geom: &AirfoilGeometry) -> Result<FactorizedSystem> {
    let n = geom.n;

    // Build (N+1)×(N+1) system
    // Unknowns: [γ₀, γ₁, ..., γₙ₋₁, ψ₀]
    let mut a_matrix = DMatrix::<f64>::zeros(n + 1, n + 1);
    let mut rhs_0 = DVector::<f64>::zeros(n + 1);   // For α = 0°
    let mut rhs_90 = DVector::<f64>::zeros(n + 1);  // For α = 90°

    // Build rows 0 to N-1: boundary condition ψ_induced + ψ_freestream = ψ₀
    for i in 0..n {
        let xi = geom.x[i];
        let yi = geom.y[i];

        // Compute influence coefficients from all panels
        let result = psilin(geom, i, xi, yi);

        // Fill matrix row i with dzdg values
        for j in 0..n {
            a_matrix[(i, j)] = result.dzdg[j];
        }

        // Column n: coefficient for ψ₀ (the unknown internal stream function)
        // Boundary condition: Σ(dzdg_j * γ_j) - ψ₀ = -ψ_freestream
        a_matrix[(i, n)] = -1.0;

        // RHS: -ψ_freestream at node i
        // For α = 0°:  ψ∞ = V∞ * y (freestream from left)
        // For α = 90°: ψ∞ = -V∞ * x (freestream from below)
        rhs_0[i] = -yi;
        rhs_90[i] = xi;
    }

    // XFOIL replaces the last airfoil row with a sharp-TE bisector velocity
    // equation before appending the final Kutta row.
    if let Some((xbis, ybis, nxbis, nybis)) = geom.sharp_te_bisector_control() {
        let te_row = n - 1;
        let result = psilin_with_dqdm(geom, n, xbis, ybis, nxbis, nybis);
        for j in 0..n {
            a_matrix[(te_row, j)] = result.dqdg[j];
        }
        a_matrix[(te_row, n)] = 0.0;
        rhs_0[te_row] = -nybis;
        rhs_90[te_row] = nxbis;
    }

    // Row N: Kutta condition γ₀ + γₙ₋₁ = 0
    // (γ at upper TE + γ at lower TE = 0)
    for j in 0..=n {
        a_matrix[(n, j)] = 0.0;
    }
    a_matrix[(n, 0)] = 1.0;
    a_matrix[(n, n - 1)] = 1.0;
    rhs_0[n] = 0.0;
    rhs_90[n] = 0.0;

    // Debug output: dump AIJ matrix elements BEFORE LU factorization
    // This allows direct comparison with XFOIL's influence coefficients
    if is_debug_active() {
        // Extract first 20 diagonal elements
        let aij_diagonal: Vec<f64> = (0..20.min(n + 1))
            .map(|i| a_matrix[(i, i)])
            .collect();

        // Extract first 20 elements of row 0
        let aij_row0: Vec<f64> = (0..20.min(n + 1))
            .map(|j| a_matrix[(0, j)])
            .collect();

        // Extract first 10 elements of RHS vectors
        let rhs_0_sample: Vec<f64> = rhs_0.iter().take(10).copied().collect();
        let rhs_90_sample: Vec<f64> = rhs_90.iter().take(10).copied().collect();

        add_event(DebugEvent::aij_matrix(
            n + 1, // System size is N+1
            aij_diagonal,
            aij_row0,
            rhs_0_sample,
            rhs_90_sample,
        ));
    }

    // LU factorization and solve for both RHS vectors
    let lu = a_matrix.clone().lu();

    let solution_0 = lu
        .solve(&rhs_0)
        .ok_or(InviscidError::SingularMatrix)?;
    let solution_90 = lu
        .solve(&rhs_90)
        .ok_or(InviscidError::SingularMatrix)?;

    // Extract γ and ψ₀ from solutions
    let gamu_0: Vec<f64> = solution_0.iter().take(n).copied().collect();
    let gamu_90: Vec<f64> = solution_90.iter().take(n).copied().collect();
    let psi0_0 = solution_0[n];
    let psi0_90 = solution_90[n];

    // Debug output: emit full AIC base solutions (GGCALC equivalent)
    if is_debug_active() {
        add_event(DebugEvent::full_aic(n, gamu_0.clone(), gamu_90.clone()));
    }

    Ok(FactorizedSystem {
        gamu_0,
        gamu_90,
        psi0_0,
        psi0_90,
        geom: geom.clone(),
        lu,
    })
}

/// Build the system matrix and RHS vectors without solving.
///
/// This is useful for testing matrix values against XFOIL.
pub fn build_system_matrix(geom: &AirfoilGeometry) -> (DMatrix<f64>, DVector<f64>, DVector<f64>) {
    let n = geom.n;

    let mut a_matrix = DMatrix::<f64>::zeros(n + 1, n + 1);
    let mut rhs_0 = DVector::<f64>::zeros(n + 1);
    let mut rhs_90 = DVector::<f64>::zeros(n + 1);

    for i in 0..n {
        let xi = geom.x[i];
        let yi = geom.y[i];

        let result = psilin(geom, i, xi, yi);

        for j in 0..n {
            a_matrix[(i, j)] = result.dzdg[j];
        }

        a_matrix[(i, n)] = -1.0;
        rhs_0[i] = -yi;
        rhs_90[i] = xi;
    }

    for j in 0..=n {
        a_matrix[(n, j)] = 0.0;
    }
    a_matrix[(n, 0)] = 1.0;
    a_matrix[(n, n - 1)] = 1.0;
    rhs_0[n] = 0.0;
    rhs_90[n] = 0.0;

    (a_matrix, rhs_0, rhs_90)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_naca0012(n_panels: usize) -> Vec<(f64, f64)> {
        let n_half = n_panels / 2;
        let t = 0.12;

        let x_coords: Vec<f64> = (0..=n_half)
            .map(|i| {
                let beta = PI * (i as f64) / (n_half as f64);
                0.5 * (1.0 - beta.cos())
            })
            .collect();

        let thickness = |x: f64| -> f64 {
            5.0 * t
                * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
                    - 0.1036 * x.powi(4))
        };

        let mut points = Vec::with_capacity(2 * n_half);

        for i in (0..=n_half).rev() {
            let x = x_coords[i];
            let y = thickness(x);
            points.push((x, y));
        }

        for i in 1..=n_half {
            let x = x_coords[i];
            let y = -thickness(x);
            points.push((x, y));
        }

        points
    }

    #[test]
    fn test_matrix_dimensions() {
        let points = make_naca0012(40);
        let geom = AirfoilGeometry::from_points(&points).unwrap();
        let (a, rhs_0, rhs_90) = build_system_matrix(&geom);

        assert_eq!(a.nrows(), geom.n + 1);
        assert_eq!(a.ncols(), geom.n + 1);
        assert_eq!(rhs_0.len(), geom.n + 1);
        assert_eq!(rhs_90.len(), geom.n + 1);
    }

    #[test]
    fn test_kutta_condition_row() {
        let points = make_naca0012(40);
        let geom = AirfoilGeometry::from_points(&points).unwrap();
        let (a, rhs_0, rhs_90) = build_system_matrix(&geom);

        let n = geom.n;

        // Kutta row should have 1.0 at columns 0 and n-1
        assert!((a[(n, 0)] - 1.0).abs() < 1e-10);
        assert!((a[(n, n - 1)] - 1.0).abs() < 1e-10);

        // All other entries in Kutta row should be 0
        for j in 1..n - 1 {
            assert!(a[(n, j)].abs() < 1e-10);
        }
        assert!(a[(n, n)].abs() < 1e-10);

        // RHS should be 0
        assert!(rhs_0[n].abs() < 1e-10);
        assert!(rhs_90[n].abs() < 1e-10);
    }

    #[test]
    fn test_factorization_succeeds() {
        let points = make_naca0012(40);
        let geom = AirfoilGeometry::from_points(&points).unwrap();

        let result = build_and_factorize(&geom);
        assert!(result.is_ok(), "Factorization should succeed");

        let factorized = result.unwrap();
        assert_eq!(factorized.gamu_0.len(), geom.n);
        assert_eq!(factorized.gamu_90.len(), geom.n);
    }

    #[test]
    fn test_kutta_satisfied() {
        let points = make_naca0012(40);
        let geom = AirfoilGeometry::from_points(&points).unwrap();
        let factorized = build_and_factorize(&geom).unwrap();

        let n = geom.n;

        // Kutta: γ₀ + γₙ₋₁ ≈ 0 for both base solutions
        let kutta_0 = factorized.gamu_0[0] + factorized.gamu_0[n - 1];
        let kutta_90 = factorized.gamu_90[0] + factorized.gamu_90[n - 1];

        assert!(
            kutta_0.abs() < 1e-8,
            "Kutta not satisfied for α=0°: {}",
            kutta_0
        );
        assert!(
            kutta_90.abs() < 1e-8,
            "Kutta not satisfied for α=90°: {}",
            kutta_90
        );
    }

    #[test]
    fn test_symmetric_airfoil_zero_alpha_zero_lift() {
        let points = make_naca0012(60);
        let geom = AirfoilGeometry::from_points(&points).unwrap();
        let factorized = build_and_factorize(&geom).unwrap();

        let flow = FlowConditions::default(); // α = 0
        let solution = factorized.solve_alpha(&flow);

        // Symmetric airfoil at α=0 should have CL ≈ 0
        assert!(
            solution.cl.abs() < 0.1,
            "CL at α=0 should be near 0, got {}",
            solution.cl
        );
    }

    #[test]
    fn test_lift_increases_with_alpha() {
        let points = make_naca0012(60);
        let geom = AirfoilGeometry::from_points(&points).unwrap();
        let factorized = build_and_factorize(&geom).unwrap();

        let flow_0 = FlowConditions::with_alpha_deg(0.0);
        let flow_4 = FlowConditions::with_alpha_deg(4.0);

        let sol_0 = factorized.solve_alpha(&flow_0);
        let sol_4 = factorized.solve_alpha(&flow_4);

        assert!(
            sol_4.cl > sol_0.cl,
            "CL should increase with α: CL(0°)={}, CL(4°)={}",
            sol_0.cl,
            sol_4.cl
        );
    }

}
