//! Boundary layer closure relations
//!
//! Each closure function returns both the value and all partial derivatives,
//! which are required for the Newton-Raphson Jacobian.

pub mod cf;
pub mod density;
pub mod dissipation;
pub mod hkin;
pub mod hs;
pub mod transition;

pub use cf::{cf_laminar, cf_turbulent, CfResult};
pub use density::{density_shape, HctResult};
pub use dissipation::{
    dissipation_laminar, dissipation_turbulent, dissipation_wake, DissipationResult,
    TurbDissipationResult,
};
pub use hkin::{hkin, HkinResult};
pub use hs::{hs_laminar, hs_turbulent, HsResult};
pub use transition::{amplification_rate, check_transition, AmplificationResult};
