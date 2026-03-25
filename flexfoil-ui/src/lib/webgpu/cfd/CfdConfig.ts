/**
 * CFD solver configuration types.
 * Must match the Rust CfdConfig struct and WGSL CfdParams struct.
 */

export type PhysicsMode = 'euler' | 'laminar_ns' | 'rans_sa';
export type ReconstructionMode = 'muscl' | 'weno5';
export type TimeSteppingMode = 'explicit_euler' | 'dadi';

export interface CfdConfig {
  ni: number;
  nj: number;
  gamma: number;
  cfl: number;
  machInf: number;
  alphaDeg: number;
  reynolds: number;
  prandtl: number;
  physics: PhysicsMode;
  reconstruction: ReconstructionMode;
  timeStepping: TimeSteppingMode;
  farField: number;
  ds0: number;
}

export const DEFAULT_CFD_CONFIG: CfdConfig = {
  ni: 256,
  nj: 128,
  gamma: 1.4,
  cfl: 0.5,
  machInf: 0.5,
  alphaDeg: 2.0,
  reynolds: 1e6,
  prandtl: 0.72,
  physics: 'euler',
  reconstruction: 'muscl',
  timeStepping: 'explicit_euler',
  farField: 20.0,
  ds0: 1e-4,
};

const PHYSICS_MAP: Record<PhysicsMode, number> = {
  euler: 0,
  laminar_ns: 1,
  rans_sa: 2,
};

const RECON_MAP: Record<ReconstructionMode, number> = {
  muscl: 0,
  weno5: 1,
};

/**
 * Pack CfdConfig into a Float32Array/Uint32Array for GPU uniform upload.
 * Layout matches the WGSL CfdParams struct (48 bytes = 12 x u32/f32).
 */
export function packCfdParams(
  config: CfdConfig,
  dt: number,
  iteration: number,
): ArrayBuffer {
  const buffer = new ArrayBuffer(48);
  const u32 = new Uint32Array(buffer);
  const f32 = new Float32Array(buffer);

  u32[0] = config.ni;
  u32[1] = config.nj;
  f32[2] = config.gamma;
  f32[3] = config.cfl;
  f32[4] = config.machInf;
  f32[5] = (config.alphaDeg * Math.PI) / 180.0;
  f32[6] = config.reynolds;
  f32[7] = config.prandtl;
  f32[8] = dt;
  u32[9] = iteration;
  u32[10] = PHYSICS_MAP[config.physics];
  u32[11] = RECON_MAP[config.reconstruction];

  return buffer;
}

/** Number of conservative variables per cell */
export const NVAR = 5;
