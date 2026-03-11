//! Standalone XFOIL-faithful side path.
//!
//! This crate owns a separate, array-first solver stack intended to mirror
//! XFOIL's process topology without routing through the modular RustFoil
//! viscous solver.

pub mod assembly;
pub mod config;
pub mod error;
pub mod forces;
pub mod march;
pub mod oper;
pub mod result;
pub mod solve;
pub mod state;
pub mod state_ops;
pub mod update;
pub mod wake_panel;

pub use config::{OperatingMode, XfoilOptions};
pub use error::{Result, XfoilError};
pub use oper::{
    solve_body_oper_point, solve_coords_oper_point, solve_operating_point_from_state, AlphaSpec,
};
pub use result::XfoilViscousResult;
pub use state::{XfoilBlRow, XfoilState, XfoilSurface};
