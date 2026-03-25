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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let cfg = CfdConfig::default();
        assert_eq!(cfg.ni, 256);
        assert_eq!(cfg.nj, 128);
        assert!((cfg.gamma - 1.4).abs() < 1e-6);
        assert!((cfg.prandtl - 0.72).abs() < 1e-6);
    }

    #[test]
    fn test_gpu_params_size() {
        // Must be 48 bytes (12 x 4-byte fields) to match WGSL struct
        assert_eq!(std::mem::size_of::<CfdParamsGpu>(), 48);
    }

    #[test]
    fn test_gpu_params_from_config() {
        let cfg = CfdConfig {
            ni: 64,
            nj: 32,
            mach_inf: 0.5,
            alpha: 0.035, // ~2 degrees
            physics: PhysicsMode::Euler,
            reconstruction: ReconstructionMode::Muscl,
            ..CfdConfig::default()
        };
        let params = CfdParamsGpu::from_config(&cfg, 0.001, 42);
        assert_eq!(params.ni, 64);
        assert_eq!(params.nj, 32);
        assert_eq!(params.iteration, 42);
        assert_eq!(params.physics_mode, 0); // Euler
        assert_eq!(params.reconstruction, 0); // MUSCL
        assert!((params.dt - 0.001).abs() < 1e-8);
    }

    #[test]
    fn test_gpu_params_bytes() {
        let cfg = CfdConfig::default();
        let params = CfdParamsGpu::from_config(&cfg, 0.01, 0);
        let bytes = params.as_bytes();
        assert_eq!(bytes.len(), 48);
        // First 4 bytes should be ni=256 as u32 little-endian
        let ni = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        assert_eq!(ni, 256);
    }

    #[test]
    fn test_physics_mode_values() {
        assert_eq!(PhysicsMode::Euler as u32, 0);
        assert_eq!(PhysicsMode::LaminarNS as u32, 1);
        assert_eq!(PhysicsMode::RansSA as u32, 2);
    }

    #[test]
    fn test_serde_roundtrip() {
        let cfg = CfdConfig::default();
        let json = serde_json::to_string(&cfg).unwrap();
        let cfg2: CfdConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(cfg2.ni, cfg.ni);
        assert_eq!(cfg2.physics, cfg.physics);
    }
}
