//! Viscous solver module implementing XFOIL's VISCAL main loop.
//!
//! This module provides the complete viscous flow solver that couples the
//! boundary layer equations (from `rustfoil-bl`) with the inviscid solution
//! through the viscous-inviscid interaction scheme (from `rustfoil-coupling`).
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────────┐
//! │                        solve_viscous()                              │
//! ├─────────────────────────────────────────────────────────────────────┤
//! │ 1. Initialize BL from stagnation (setup.rs)                        │
//! │ 2. March BL with initial Ue (march_fixed_ue)                       │
//! │ 3. Newton iteration loop:                                          │
//! │    a. Build Newton system (CoupledNewtonSystem)                    │
//! │    b. Solve block-tridiagonal system (solve_bl_system)             │
//! │    c. Update stations and Ue (update_stations, set_edge_velocities)│
//! │    d. Check convergence                                            │
//! │ 4. Compute forces (forces.rs)                                      │
//! └─────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # XFOIL Reference
//! - VISCAL: xoper.f line 2886

pub mod config;
pub mod circulation;
pub mod forces;
pub mod setup;
pub mod viscal;

// Re-export main types for convenience
pub use config::ViscousSolverConfig;
pub use forces::{compute_forces, AeroForces};
pub use setup::{
    compute_arc_from_stagnation, compute_arc_lengths, extract_surface, extract_surface_xfoil,
    find_stagnation, find_stagnation_by_sign_change, initialize_bl_stations,
    initialize_surface_stations, initialize_surface_stations_with_coords, 
    initialize_surface_stations_with_panel_idx, interpolate_stagnation,
    ViscousSetup,
    // New integration with rustfoil-inviscid
    setup_from_body, setup_from_coords, SetupError, ViscousSetupResult,
};
pub use viscal::{
    solve_viscous, solve_viscous_polar_parallel, solve_viscous_two_surfaces, ViscousResult,
};
