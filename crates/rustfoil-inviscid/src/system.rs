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
use crate::influence::{build_source_influence_matrix, psilin};
use crate::solution::{FlowConditions, InviscidSolution};
use crate::{InviscidError, Result};
use nalgebra::{DMatrix, DVector, LU};

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

    // Row N: Kutta condition γ₀ + γₙ₋₁ = 0
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
