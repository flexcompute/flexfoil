use serde::{Deserialize, Serialize};

/// Physics mode for the CFD solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PhysicsMode {
    /// Inviscid compressible (Euler equations)
    Euler = 0,
    /// Laminar Navier-Stokes
    LaminarNS = 1,
    /// Reynolds-Averaged Navier-Stokes with Spalart-Allmaras
    RansSA = 2,
}

/// Spatial reconstruction scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReconstructionMode {
    /// 2nd-order MUSCL with minmod limiter
    Muscl = 0,
    /// 5th-order WENO5
    Weno5 = 1,
}

/// Time-stepping scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TimeSteppingMode {
    /// Explicit forward Euler (Phase 1)
    ExplicitEuler = 0,
    /// Diagonalized ADI (Phase 2)
    Dadi = 1,
}

/// CFD solver configuration.
///
/// This struct is shared between Rust/WASM and TypeScript (via serde).
/// Fields are laid out to match the GPU uniform buffer `CfdParams`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CfdConfig {
    /// Circumferential grid points (around airfoil)
    pub ni: u32,
    /// Radial grid points (away from airfoil)
    pub nj: u32,
    /// Ratio of specific heats (default 1.4 for air)
    pub gamma: f32,
    /// CFL number for time stepping
    pub cfl: f32,
    /// Freestream Mach number
    pub mach_inf: f32,
    /// Angle of attack in radians
    pub alpha: f32,
    /// Reynolds number (used for NS/RANS)
    pub reynolds: f32,
    /// Prandtl number (default 0.72 for air)
    pub prandtl: f32,
    /// Physics mode
    pub physics: PhysicsMode,
    /// Reconstruction scheme
    pub reconstruction: ReconstructionMode,
    /// Time-stepping scheme
    pub time_stepping: TimeSteppingMode,
    /// Far-field distance in chord lengths
    pub far_field: f32,
    /// First cell wall-normal spacing (for viscous grids)
    pub ds0: f32,
}

impl Default for CfdConfig {
    fn default() -> Self {
        Self {
            ni: 256,
            nj: 128,
            gamma: 1.4,
            cfl: 0.5,
            mach_inf: 0.5,
            alpha: 0.0,
            reynolds: 1e6,
            prandtl: 0.72,
            physics: PhysicsMode::Euler,
            reconstruction: ReconstructionMode::Muscl,
            time_stepping: TimeSteppingMode::ExplicitEuler,
            far_field: 20.0,
            ds0: 1e-4,
        }
    }
}

/// GPU uniform buffer layout for CfdParams.
/// Must match the WGSL struct exactly (std140/std430 layout).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CfdParamsGpu {
    pub ni: u32,
    pub nj: u32,
    pub gamma: f32,
    pub cfl: f32,
    pub mach_inf: f32,
    pub alpha: f32,
    pub reynolds: f32,
    pub prandtl: f32,
    pub dt: f32,
    pub iteration: u32,
    pub physics_mode: u32,
    pub reconstruction: u32,
}

impl CfdParamsGpu {
    pub fn from_config(config: &CfdConfig, dt: f32, iteration: u32) -> Self {
        Self {
            ni: config.ni,
            nj: config.nj,
            gamma: config.gamma,
            cfl: config.cfl,
            mach_inf: config.mach_inf,
            alpha: config.alpha,
            reynolds: config.reynolds,
            prandtl: config.prandtl,
            dt,
            iteration,
            physics_mode: config.physics as u32,
            reconstruction: config.reconstruction as u32,
        }
    }

    /// Return as raw bytes for GPU buffer upload.
    pub fn as_bytes(&self) -> &[u8] {
        unsafe {
            std::slice::from_raw_parts(
                self as *const Self as *const u8,
                std::mem::size_of::<Self>(),
            )
        }
    }
}
