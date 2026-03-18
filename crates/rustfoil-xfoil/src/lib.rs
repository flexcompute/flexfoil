//! Standalone XFOIL-faithful side path.
//!
//! This crate owns a separate, array-first solver stack intended to mirror
//! XFOIL's process topology without routing through the modular RustFoil
//! viscous solver.

pub mod assembly;
pub(crate) mod canonical_state;
pub mod config;
pub mod error;
pub mod forces;
pub mod march;
pub mod mdes;
pub mod oper;
pub mod qdes;
pub mod result;
pub mod solve;
pub mod state;
pub mod state_ops;
pub mod update;
pub mod wake_panel;

pub use config::{OperatingMode, XfoilOptions};
pub use error::{Result, XfoilError};
pub use oper::{
    build_state_from_coords, coords_from_body, solve_body_oper_point, solve_coords_oper_point,
    solve_operating_point_from_state, AlphaSpec,
};
pub use qdes::{solve_body_qdes, solve_coords_qdes, InverseDesignSession, QdesOptions};
pub use result::{
    QdesIterationSnapshot, QdesResult, QdesSpec, QdesTarget, QdesTargetKind, SurfaceDistribution,
    XfoilViscousResult,
};
pub use state::{XfoilBlRow, XfoilState, XfoilSurface};
