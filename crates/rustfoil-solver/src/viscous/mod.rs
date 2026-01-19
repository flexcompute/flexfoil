//! Viscous-Inviscid Coupling Module
//!
//! This module implements the viscous-inviscid interaction (VII) solver
//! that couples the inviscid panel method with the boundary layer equations.
//!
//! The coupling supports multiple methods:
//! - **SemiDirect**: Sequential iteration without transpiration feedback
//! - **Transpiration**: Sequential with transpiration velocity correction
//! - **FullNewton**: Global Newton-Raphson simultaneous solution
//!
//! ## Algorithm (Newton)
//! 1. Assemble residuals from panel and BL equations
//! 2. Build Jacobian matrix
//! 3. Solve Newton step: J·Δx = -R
//! 4. Apply line search for robustness
//! 5. Iterate until residual converges
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
