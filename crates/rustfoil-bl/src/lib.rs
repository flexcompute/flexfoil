//! Boundary layer equations and closures for RustFoil
//!
//! This crate implements XFOIL's integral boundary layer formulation.
//!
//! # Modules
//!
//! - [`closures`] - Closure relations (HKIN, HS, CF, DI, HCT, DAMPL)
//! - [`constants`] - BLPAR constants from XFOIL
//! - [`state`] - BlStation struct for BL state at a single point
//! - [`equations`] - BLVAR/BLDIF for computing secondary variables and residuals
//! - [`debug`] - Debug output for XFOIL comparison

pub mod closures;
pub mod constants;
pub mod debug;
pub mod equations;
pub mod state;

pub use constants::*;
pub use debug::{
    add_event, finalize_debug, init_debug, is_debug_active, BlvarInput, BlvarOutput, DebugEvent,
};
pub use equations::{bldif, bldif_debug, blvar, blvar_debug, BlJacobian, BlResiduals, FlowType};
pub use state::{BlDerivatives, BlStation};
