//! Viscous-Inviscid Coupling Module
//!
//! This module implements the viscous-inviscid interaction (VII) solver
//! that couples the inviscid panel method with the boundary layer equations.
//!
//! ## Coupling Methods
//!
//! - **Transpiration** (RECOMMENDED, default): Sequential with transpiration velocity
//!   correction. Achieves <3% accuracy for both Cl and Cd. Use this for all cases.
//! - **FullNewton**: Global Newton-Raphson simultaneous solution. Experimental.
//! - **SemiDirect**: DEPRECATED - has ~180% Cd error due to missing feedback.
//!
//! ## Algorithm (Transpiration)
//! 1. Solve inviscid flow for initial Ue(s)
//! 2. March boundary layer with current Ue → get δ*(s)
//! 3. Compute transpiration velocity: Vn = d(Ue·δ*)/ds
//! 4. Re-solve inviscid with transpiration BC → get new Ue
//! 5. Iterate until δ* converges
//!
//! Reference: Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds
//! Number Airfoils", MIT, 1989.

mod coupling;
/// Global Newton-Raphson VII solver.
pub mod newton;
/// XFOIL-style BL system equations and analytical derivatives.
pub mod blsys;

pub use coupling::{ViscousSolver, ViscousSolution, ViscousConfig, CouplingMethod, compute_transpiration};
pub use newton::{
    NewtonConfig, NewtonGeometry, NewtonState, NewtonStep, Residuals, NewtonResult, NewtonVIISolver,
    // XFOIL-style block structures
    BLBlock, NewtonScaling, VaccelParams, BlockTridiagJacobian, StationState,
};
pub use blsys::{
    BLClosures, BLClosureDerivs, BLSystemResiduals,
    compute_interval_residuals, compute_interval_jacobian,
    LocalNewtonConfig, LocalNewtonResult, solve_station_newton,
    march_bl_surface, build_block_jacobian,
};
