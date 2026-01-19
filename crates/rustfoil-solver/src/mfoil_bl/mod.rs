//! Mfoil-style boundary layer solver.
//!
//! This module implements a global Newton boundary layer solver that mirrors
//! the approach used in mfoil.py by Krzysztof J. Fidkowski (2023).
//!
//! # State Vector
//!
//! The primary state at each node is `U = [θ, δ*, sa, Ue]`:
//! - `θ` (theta): momentum thickness
//! - `δ*` (delta_star): displacement thickness  
//! - `sa`: amplification factor (laminar) or √(Cτ) shear stress (turbulent)
//! - `Ue`: edge velocity
//!
//! # Residual Equations
//!
//! Three equations per station:
//! 1. **Momentum**: von Kármán momentum integral
//! 2. **Shape**: Kinetic energy shape parameter (Hs) evolution
//! 3. **Lag/Amplification**: Shear lag (turbulent) or amplification (laminar)
//!
//! A fourth equation for V-I coupling:
//!   `Ue = Ue_inviscid + ue_m @ (δ* · Ue)`
//!
//! # Surfaces
//!
//! Three surfaces are tracked:
//! - Lower surface (si=0): From stagnation to lower TE
//! - Upper surface (si=1): From stagnation to upper TE
//! - Wake (si=2): From TE downstream
//!
//! # Reference
//!
//! Based on mfoil.py v2023-06-28 by Krzysztof J. Fidkowski.

mod state;
mod closures;
mod residuals;
mod solver;

pub use state::{MfoilState, MfoilParam, SurfaceIndex, MfoilGeom, MfoilPanel};
pub use closures::{
    get_h, get_hk, get_hs, get_hw, get_ret, get_cf, get_de, get_cteq, get_cttr,
    get_damp, get_upw, upwind,
};
pub use residuals::{residual_station, residual_transition};
pub use solver::{MfoilSolver, MfoilConfig, MfoilSolution, MfoilIsol, MfoilVsol};
