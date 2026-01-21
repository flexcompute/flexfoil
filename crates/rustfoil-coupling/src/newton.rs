//! Newton system construction for viscous-inviscid coupling
//!
//! This module implements the Newton system builder for the coupled
//! viscous-inviscid iteration, porting XFOIL's BLSYS subroutine.
//!
//! # System Structure
//!
//! The BL Newton system has the form:
//! ```text
//! |       ||dA1|     |       ||dA2|       |     |
//! |  VS1  ||dT1|  +  |  VS2  ||dT2|   =   |VSREZ|
//! |       ||dD1|     |       ||dD2|       |     |
//!          |dU1|              |dU2|
//!          |dX1|              |dX2|
//!
//!   3x5    5x1         3x5    5x1          3x1
//! ```
//!
//! This forms a block-tridiagonal system where each station has 3 BL equations
//! (amplification/shear-lag, momentum, shape parameter) with coupling to
//! edge velocity through the VS terms.
//!
//! # XFOIL Reference
//! - BLSYS: xblsys.f line 583
//! - SETBL: xbl.f line 21

use nalgebra::{DMatrix, DVector};
use rustfoil_bl::equations::{bldif, FlowType};
use rustfoil_bl::state::BlStation;

/// Block-tridiagonal system for BL Newton iteration
///
/// This structure stores the linearized BL equations in block-tridiagonal form.
/// Each interval between stations i-1 and i contributes:
/// - VA[i]: Diagonal block (derivatives w.r.t. station i variables)
/// - VB[i]: Lower diagonal block (derivatives w.r.t. station i-1 variables)
/// - RHS[i]: Residual vector
/// - VS[i]: Edge velocity coupling (for viscous-inviscid interaction)
///
/// The 3 equations at each station are:
/// 1. Amplification (laminar) or shear-lag (turbulent/wake)
/// 2. Momentum integral
/// 3. Shape parameter (kinetic energy)
///
/// # XFOIL Reference
/// XFOIL xblsys.f BLSYS (line 583)
#[derive(Debug, Clone)]
pub struct BlNewtonSystem {
    /// Number of stations
    pub n: usize,
    /// Diagonal blocks VA[i] (3x3) - derivatives w.r.t. downstream station
    /// Index 0 is unused (boundary condition), indices 1..n are valid
    pub va: Vec<[[f64; 3]; 3]>,
    /// Lower diagonal blocks VB[i] (3x3) - derivatives w.r.t. upstream station
    /// Index 0 is unused, indices 1..n are valid
    pub vb: Vec<[[f64; 3]; 3]>,
    /// Right-hand side residuals (3 per interval)
    /// Index 0 is unused, indices 1..n are valid
    pub rhs: Vec<[f64; 3]>,
    /// Coupling to inviscid (dUe influence)
    /// These are the VS2[:,3] entries from bldif - edge velocity sensitivity
    pub vs: Vec<[f64; 3]>,
}

impl BlNewtonSystem {
    /// Create a new Newton system for n stations
    ///
    /// # Arguments
    /// * `n` - Number of BL stations (including stagnation point)
    ///
    /// # Returns
    /// A new BlNewtonSystem with all blocks initialized to zero
    pub fn new(n: usize) -> Self {
        Self {
            n,
            va: vec![[[0.0; 3]; 3]; n],
            vb: vec![[[0.0; 3]; 3]; n],
            rhs: vec![[0.0; 3]; n],
            vs: vec![[0.0; 3]; n],
        }
    }

    /// Build the Newton system from BL stations
    ///
    /// This method constructs the block-tridiagonal Newton system by calling
    /// `bldif()` for each interval between adjacent stations. The resulting
    /// Jacobian blocks and residuals are extracted and stored.
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval (Laminar, Turbulent, or Wake)
    /// * `msq` - Mach number squared (M²)
    /// * `re` - Reference Reynolds number
    ///
    /// # Note
    /// The `stations` array should have secondary variables already computed
    /// via `blvar()`. This method extracts the 3x3 core blocks from the full
    /// 3x5 Jacobian matrices returned by `bldif()`.
    ///
    /// # Reference
    /// XFOIL xblsys.f BLSYS (line 583)
    pub fn build(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
    ) {
        assert!(
            stations.len() >= 2,
            "Need at least 2 stations to build Newton system"
        );
        assert_eq!(
            stations.len(),
            self.n,
            "Station count must match system size"
        );
        assert_eq!(
            flow_types.len(),
            self.n - 1,
            "Need n-1 flow types for n stations"
        );

        // For each interval (between stations i-1 and i)
        for i in 1..self.n {
            let s1 = &stations[i - 1]; // Upstream station
            let s2 = &stations[i]; // Downstream station
            let flow_type = flow_types[i - 1];

            // Compute residuals and Jacobian for this interval
            let (residuals, jacobian) = bldif(s1, s2, flow_type, msq, re);

            // Store residuals (negated to match Newton iteration convention)
            // Note: bldif already negates residuals, so we use them directly
            self.rhs[i] = [residuals.res_third, residuals.res_mom, residuals.res_shape];

            // Extract 3x3 blocks from the 3x5 Jacobian matrices
            // The first 3 columns correspond to [ampl/ctau, theta, delta_star]
            //
            // VA[i] comes from vs2 (downstream station derivatives)
            // VB[i] comes from vs1 (upstream station derivatives)

            for eq in 0..3 {
                for var in 0..3 {
                    self.va[i][eq][var] = jacobian.vs2[eq][var];
                    self.vb[i][eq][var] = jacobian.vs1[eq][var];
                }
                // Store edge velocity coupling (column 3 = u variable)
                self.vs[i][eq] = jacobian.vs2[eq][3];
            }
        }
    }

    /// Get maximum absolute residual (for convergence checking)
    ///
    /// This returns the infinity norm of the residual vector, which is
    /// commonly used to check Newton iteration convergence.
    ///
    /// # Returns
    /// The maximum absolute value of any residual component
    pub fn max_residual(&self) -> f64 {
        self.rhs
            .iter()
            .skip(1) // Skip index 0 (unused)
            .flat_map(|r| r.iter())
            .map(|&v| v.abs())
            .fold(0.0, f64::max)
    }

    /// Get the L2 norm of residuals (alternative convergence metric)
    ///
    /// # Returns
    /// The Euclidean norm of the residual vector
    pub fn residual_norm(&self) -> f64 {
        self.rhs
            .iter()
            .skip(1)
            .flat_map(|r| r.iter())
            .map(|&v| v * v)
            .sum::<f64>()
            .sqrt()
    }

    /// Check if the system has converged
    ///
    /// # Arguments
    /// * `tolerance` - Maximum allowed residual
    ///
    /// # Returns
    /// True if max_residual() < tolerance
    pub fn is_converged(&self, tolerance: f64) -> bool {
        self.max_residual() < tolerance
    }
}

/// Full coupled Newton system including inviscid interaction
///
/// This structure combines the BL block-tridiagonal system with the
/// inviscid coupling through the DIJ matrix. The DIJ matrix relates
/// changes in mass defect to changes in edge velocity.
///
/// # System Structure
/// The full coupled system solves for:
/// - BL primary variables at each station (ampl/ctau, theta, delta_star)
/// - Edge velocity updates through viscous-inviscid interaction
///
/// # XFOIL Reference
/// XFOIL xbl.f SETBL (line 21)
#[derive(Debug, Clone)]
pub struct CoupledNewtonSystem {
    /// BL system (block tridiagonal)
    pub bl: BlNewtonSystem,
    /// DIJ matrix for mass defect coupling
    /// ΔUe_i = Σ_j DIJ[i,j] × Δ(Ue×δ*)_j
    pub dij: DMatrix<f64>,
    /// Global system matrix (assembled from BL and coupling)
    pub matrix: DMatrix<f64>,
    /// Global RHS vector
    pub rhs: DVector<f64>,
}

impl CoupledNewtonSystem {
    /// Create a new coupled system
    ///
    /// # Arguments
    /// * `n_stations` - Number of BL stations
    /// * `dij` - Pre-computed DIJ influence matrix
    ///
    /// # Returns
    /// A new CoupledNewtonSystem with the BL system and DIJ matrix initialized
    pub fn new(n_stations: usize, dij: DMatrix<f64>) -> Self {
        // Each station has 3 BL variables
        let n_vars = 3 * n_stations;

        Self {
            bl: BlNewtonSystem::new(n_stations),
            dij,
            matrix: DMatrix::zeros(n_vars, n_vars),
            rhs: DVector::zeros(n_vars),
        }
    }

    /// Build the BL portion of the coupled system
    ///
    /// # Arguments
    /// * `stations` - Slice of BlStation with pre-computed secondary variables
    /// * `flow_types` - Flow type for each interval
    /// * `msq` - Mach number squared
    /// * `re` - Reference Reynolds number
    pub fn build_bl_system(
        &mut self,
        stations: &[BlStation],
        flow_types: &[FlowType],
        msq: f64,
        re: f64,
    ) {
        self.bl.build(stations, flow_types, msq, re);
    }

    /// Assemble the full coupled system matrix
    ///
    /// This combines the BL block-tridiagonal system with the inviscid
    /// coupling through the DIJ matrix. The assembled system can then
    /// be solved for the Newton update.
    ///
    /// # Note
    /// The BL system must be built first via `build_bl_system()`.
    ///
    /// # Reference
    /// XFOIL xbl.f SETBL (line 21)
    pub fn assemble(&mut self) {
        let n = self.bl.n;
        let n_vars = 3 * n;

        // Reset matrix and RHS
        self.matrix = DMatrix::zeros(n_vars, n_vars);
        self.rhs = DVector::zeros(n_vars);

        // Assemble block-tridiagonal BL system into global matrix
        for i in 1..n {
            let row_base = 3 * i;
            let col_base_diag = 3 * i; // Diagonal block (station i)
            let col_base_lower = 3 * (i - 1); // Lower diagonal (station i-1)

            // Insert VA[i] (diagonal block)
            for eq in 0..3 {
                for var in 0..3 {
                    self.matrix[(row_base + eq, col_base_diag + var)] = self.bl.va[i][eq][var];
                }
            }

            // Insert VB[i] (lower diagonal block)
            for eq in 0..3 {
                for var in 0..3 {
                    self.matrix[(row_base + eq, col_base_lower + var)] = self.bl.vb[i][eq][var];
                }
            }

            // Copy RHS
            for eq in 0..3 {
                self.rhs[row_base + eq] = self.bl.rhs[i][eq];
            }
        }

        // Add viscous-inviscid coupling through DIJ
        // The coupling adds terms: dR/dUe * dUe/d(mass_defect) * d(mass_defect)/d(primary_vars)
        // This is a simplified version - full implementation would include
        // the complete SETBL coupling logic
        //
        // TODO: Implement full SETBL coupling with DIJ matrix
    }

    /// Get maximum residual from the BL system
    pub fn max_residual(&self) -> f64 {
        self.bl.max_residual()
    }

    /// Check convergence of the coupled system
    pub fn is_converged(&self, tolerance: f64) -> bool {
        self.bl.is_converged(tolerance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfoil_bl::equations::{blvar, FlowType};

    // =========================================================================
    // BlNewtonSystem Creation Tests
    // =========================================================================

    #[test]
    fn test_system_creation() {
        let system = BlNewtonSystem::new(50);

        assert_eq!(system.n, 50);
        assert_eq!(system.va.len(), 50);
        assert_eq!(system.vb.len(), 50);
        assert_eq!(system.rhs.len(), 50);
        assert_eq!(system.vs.len(), 50);

        // All values should be initialized to zero
        for i in 0..50 {
            for eq in 0..3 {
                for var in 0..3 {
                    assert_eq!(system.va[i][eq][var], 0.0);
                    assert_eq!(system.vb[i][eq][var], 0.0);
                }
                assert_eq!(system.rhs[i][eq], 0.0);
                assert_eq!(system.vs[i][eq], 0.0);
            }
        }
    }

    #[test]
    fn test_system_small() {
        // Minimum size system (2 stations)
        let system = BlNewtonSystem::new(2);

        assert_eq!(system.n, 2);
        assert_eq!(system.va.len(), 2);
    }

    // =========================================================================
    // Build Tests - Laminar Flow
    // =========================================================================

    #[test]
    fn test_build_laminar() {
        // Create a simple laminar BL with 5 stations
        let n = 5;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        // Create stations with increasing x and growing BL thickness
        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 + 0.1 * i as f64;
            station.u = 1.0 - 0.01 * i as f64; // Slight deceleration
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);
            station.ampl = 0.5 * i as f64; // Growing amplification

            // Compute secondary variables
            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        // Build the Newton system
        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Check that residuals are finite
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Residual [{},{}] should be finite, got {}",
                    i,
                    eq,
                    system.rhs[i][eq]
                );
            }
        }

        // Check that Jacobian blocks are non-trivial
        let mut has_nonzero_va = false;
        let mut has_nonzero_vb = false;
        for i in 1..n {
            for eq in 0..3 {
                for var in 0..3 {
                    if system.va[i][eq][var].abs() > 1e-15 {
                        has_nonzero_va = true;
                    }
                    if system.vb[i][eq][var].abs() > 1e-15 {
                        has_nonzero_vb = true;
                    }
                }
            }
        }
        assert!(has_nonzero_va, "VA blocks should have non-zero entries");
        assert!(has_nonzero_vb, "VB blocks should have non-zero entries");
    }

    #[test]
    fn test_build_turbulent() {
        // Create a turbulent BL
        let n = 5;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.3 + 0.1 * i as f64;
            station.u = 0.9 - 0.02 * i as f64;
            station.theta = 0.003 * (1.0 + 0.15 * i as f64);
            station.delta_star = 0.005 * (1.0 + 0.15 * i as f64);
            station.ctau = 0.1 + 0.01 * i as f64; // Shear stress coefficient
            station.is_laminar = false;
            station.is_turbulent = true;

            blvar(&mut station, FlowType::Turbulent, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Turbulent; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Check residuals are finite
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Turbulent residual [{},{}] should be finite",
                    i,
                    eq
                );
            }
        }
    }

    #[test]
    fn test_build_wake() {
        // Create a wake region
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 1.0 + 0.1 * i as f64;
            station.u = 0.5 + 0.02 * i as f64; // Recovering velocity
            station.theta = 0.01 * (1.0 + 0.05 * i as f64);
            station.delta_star = 0.02 * (1.0 + 0.05 * i as f64);
            station.ctau = 0.05 - 0.005 * i as f64; // Decaying shear
            station.is_wake = true;
            station.is_laminar = false;

            blvar(&mut station, FlowType::Wake, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Wake; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Verify residuals are finite for wake
        for i in 1..n {
            for eq in 0..3 {
                assert!(
                    system.rhs[i][eq].is_finite(),
                    "Wake residual [{},{}] should be finite",
                    i,
                    eq
                );
            }
        }
    }

    // =========================================================================
    // Jacobian Structure Tests
    // =========================================================================

    #[test]
    fn test_jacobian_block_structure_laminar() {
        // For laminar flow, the amplification equation should have
        // vs1[0][0] = -1 and vs2[0][0] = +1 (from bldif)
        let n = 3;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 * (i + 1) as f64;
            station.u = 1.0;
            station.theta = 0.001;
            station.delta_star = 0.0026;
            station.ampl = i as f64;

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // Check structure of first equation (amplification)
        // VA[i][0][0] should be close to +1 (dres/dampl2)
        // VB[i][0][0] should be close to -1 (dres/dampl1)
        for i in 1..n {
            assert!(
                (system.va[i][0][0] - 1.0).abs() < 0.1,
                "VA[{}][0][0] should be ~1 for amplification, got {}",
                i,
                system.va[i][0][0]
            );
            assert!(
                (system.vb[i][0][0] - (-1.0)).abs() < 0.1,
                "VB[{}][0][0] should be ~-1 for amplification, got {}",
                i,
                system.vb[i][0][0]
            );
        }
    }

    // =========================================================================
    // max_residual Tests
    // =========================================================================

    #[test]
    fn test_max_residual_zero() {
        let system = BlNewtonSystem::new(5);

        // All zeros should give zero max residual
        assert_eq!(system.max_residual(), 0.0);
    }

    #[test]
    fn test_max_residual_finds_max() {
        let mut system = BlNewtonSystem::new(5);

        // Set some residuals
        system.rhs[1] = [0.1, -0.2, 0.05];
        system.rhs[2] = [0.01, 0.5, -0.3];
        system.rhs[3] = [-0.15, 0.02, 0.8]; // Maximum is 0.8

        let max_res = system.max_residual();
        assert!(
            (max_res - 0.8).abs() < 1e-10,
            "max_residual should be 0.8, got {}",
            max_res
        );
    }

    #[test]
    fn test_max_residual_skips_zero_index() {
        let mut system = BlNewtonSystem::new(3);

        // Put large value at index 0 (should be skipped)
        system.rhs[0] = [100.0, 100.0, 100.0];
        system.rhs[1] = [0.1, 0.2, 0.3];
        system.rhs[2] = [0.05, 0.15, 0.25];

        let max_res = system.max_residual();
        assert!(
            (max_res - 0.3).abs() < 1e-10,
            "max_residual should skip index 0, got {}",
            max_res
        );
    }

    #[test]
    fn test_residual_norm() {
        let mut system = BlNewtonSystem::new(3);

        system.rhs[1] = [3.0, 0.0, 0.0];
        system.rhs[2] = [0.0, 4.0, 0.0];

        // L2 norm should be sqrt(9 + 16) = 5
        let norm = system.residual_norm();
        assert!(
            (norm - 5.0).abs() < 1e-10,
            "residual_norm should be 5.0, got {}",
            norm
        );
    }

    #[test]
    fn test_is_converged() {
        let mut system = BlNewtonSystem::new(3);

        system.rhs[1] = [1e-6, 1e-7, 1e-8];
        system.rhs[2] = [1e-7, 1e-6, 1e-7];

        assert!(system.is_converged(1e-5), "Should converge with tolerance 1e-5");
        assert!(!system.is_converged(1e-7), "Should not converge with tolerance 1e-7");
    }

    // =========================================================================
    // Edge Velocity Coupling Tests
    // =========================================================================

    #[test]
    fn test_vs_coupling_nonzero() {
        let n = 4;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);
        let re = 1e6;
        let msq = 0.0;

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = 0.1 * (i + 1) as f64;
            station.u = 1.0 - 0.05 * i as f64; // Decelerating flow
            station.theta = 0.001 * (1.0 + 0.3 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.3 * i as f64);

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let mut system = BlNewtonSystem::new(n);
        let flow_types = vec![FlowType::Laminar; n - 1];
        system.build(&stations, &flow_types, msq, re);

        // VS coupling should have non-zero entries (edge velocity sensitivity)
        let mut has_nonzero_vs = false;
        for i in 1..n {
            for eq in 0..3 {
                if system.vs[i][eq].abs() > 1e-15 {
                    has_nonzero_vs = true;
                }
            }
        }

        assert!(
            has_nonzero_vs,
            "VS coupling terms should have non-zero entries for viscous-inviscid interaction"
        );
    }

    // =========================================================================
    // CoupledNewtonSystem Tests
    // =========================================================================

    #[test]
    fn test_coupled_system_creation() {
        use crate::dij::build_dij_matrix;

        let n = 5;
        let x: Vec<f64> = (0..n).map(|i| 0.1 * (i + 1) as f64).collect();
        let y: Vec<f64> = (0..n).map(|i| 0.01 * i as f64).collect();

        let dij = build_dij_matrix(&x, &y);
        let coupled = CoupledNewtonSystem::new(n, dij);

        assert_eq!(coupled.bl.n, n);
        assert_eq!(coupled.dij.nrows(), n);
        assert_eq!(coupled.dij.ncols(), n);
        assert_eq!(coupled.matrix.nrows(), 3 * n);
        assert_eq!(coupled.matrix.ncols(), 3 * n);
        assert_eq!(coupled.rhs.len(), 3 * n);
    }

    #[test]
    fn test_coupled_system_build_and_assemble() {
        use crate::dij::build_dij_matrix;

        let n = 4;
        let x: Vec<f64> = (0..n).map(|i| 0.1 * (i + 1) as f64).collect();
        let y: Vec<f64> = (0..n).map(|i| 0.01 * i as f64).collect();

        let dij = build_dij_matrix(&x, &y);
        let mut coupled = CoupledNewtonSystem::new(n, dij);

        // Create stations
        let re = 1e6;
        let msq = 0.0;
        let mut stations: Vec<BlStation> = Vec::with_capacity(n);

        for i in 0..n {
            let mut station = BlStation::new();
            station.x = x[i];
            station.u = 1.0;
            station.theta = 0.001 * (1.0 + 0.2 * i as f64);
            station.delta_star = 0.0026 * (1.0 + 0.2 * i as f64);

            blvar(&mut station, FlowType::Laminar, msq, re);
            stations.push(station);
        }

        let flow_types = vec![FlowType::Laminar; n - 1];
        coupled.build_bl_system(&stations, &flow_types, msq, re);
        coupled.assemble();

        // Check that matrix has non-zero entries
        let mut has_nonzero = false;
        for i in 0..coupled.matrix.nrows() {
            for j in 0..coupled.matrix.ncols() {
                if coupled.matrix[(i, j)].abs() > 1e-15 {
                    has_nonzero = true;
                    break;
                }
            }
            if has_nonzero {
                break;
            }
        }

        assert!(has_nonzero, "Assembled matrix should have non-zero entries");
    }

    #[test]
    fn test_coupled_system_max_residual() {
        use crate::dij::build_dij_matrix;

        let n = 3;
        let x = vec![0.1, 0.2, 0.3];
        let y = vec![0.0, 0.01, 0.015];

        let dij = build_dij_matrix(&x, &y);
        let mut coupled = CoupledNewtonSystem::new(n, dij);

        // Manually set BL residuals
        coupled.bl.rhs[1] = [0.1, 0.2, 0.3];
        coupled.bl.rhs[2] = [0.05, 0.4, 0.15];

        let max_res = coupled.max_residual();
        assert!(
            (max_res - 0.4).abs() < 1e-10,
            "Coupled max_residual should be 0.4, got {}",
            max_res
        );
    }
}
