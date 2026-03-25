/**
 * GPU buffer management for the CFD solver.
 * Allocates all storage/uniform buffers needed for the compute pipeline.
 */

import { NVAR } from './CfdConfig';

export interface CfdBufferSet {
  /** Conservative variables Q = [rho, rho*u, rho*v, E, nu_tilde] */
  Q: GPUBuffer;
  /** Primitive variables W = [rho, u, v, p, nu_tilde] */
  W: GPUBuffer;
  /** Mesh x-coordinates */
  meshX: GPUBuffer;
  /** Mesh y-coordinates */
  meshY: GPUBuffer;
  /** Grid metrics [xi_x, xi_y, eta_x, eta_y, J] */
  metrics: GPUBuffer;
  /** Numerical flux in xi-direction */
  fluxXi: GPUBuffer;
  /** Numerical flux in eta-direction */
  fluxEta: GPUBuffer;
  /** Reconstructed left states at xi-faces */
  leftXi: GPUBuffer;
  /** Reconstructed right states at xi-faces */
  rightXi: GPUBuffer;
  /** Reconstructed left states at eta-faces */
  leftEta: GPUBuffer;
  /** Reconstructed right states at eta-faces */
  rightEta: GPUBuffer;
  /** RHS residual */
  residual: GPUBuffer;
  /** Boundary condition types (u32 per cell) */
  bcType: GPUBuffer;
  /** Solver parameters uniform buffer */
  params: GPUBuffer;
  /** Force accumulator [Cl, Cd, Cm] */
  forceAccum: GPUBuffer;
  /** Residual L2 norm */
  residualNorm: GPUBuffer;
  /** Readback buffer (MAP_READ) for async GPU->CPU data transfer */
  readback: GPUBuffer;
}

export function createCfdBuffers(device: GPUDevice, ni: number, nj: number): CfdBufferSet {
  const cellCount = ni * nj;
  const qSize = cellCount * NVAR * 4; // f32 per variable
  const coordSize = cellCount * 4;
  const metricsSize = cellCount * 5 * 4; // 5 metrics per cell

  const STORAGE_RW = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC;
  const STORAGE_R = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST;

  return {
    Q: device.createBuffer({ label: 'CFD Q', size: qSize, usage: STORAGE_RW }),
    W: device.createBuffer({ label: 'CFD W', size: qSize, usage: STORAGE_RW }),
    meshX: device.createBuffer({ label: 'CFD mesh_x', size: coordSize, usage: STORAGE_R }),
    meshY: device.createBuffer({ label: 'CFD mesh_y', size: coordSize, usage: STORAGE_R }),
    metrics: device.createBuffer({ label: 'CFD metrics', size: metricsSize, usage: STORAGE_RW }),
    fluxXi: device.createBuffer({ label: 'CFD flux_xi', size: qSize, usage: STORAGE_RW }),
    fluxEta: device.createBuffer({ label: 'CFD flux_eta', size: qSize, usage: STORAGE_RW }),
    leftXi: device.createBuffer({ label: 'CFD left_xi', size: qSize, usage: STORAGE_RW }),
    rightXi: device.createBuffer({ label: 'CFD right_xi', size: qSize, usage: STORAGE_RW }),
    leftEta: device.createBuffer({ label: 'CFD left_eta', size: qSize, usage: STORAGE_RW }),
    rightEta: device.createBuffer({ label: 'CFD right_eta', size: qSize, usage: STORAGE_RW }),
    residual: device.createBuffer({ label: 'CFD residual', size: qSize, usage: STORAGE_RW }),
    bcType: device.createBuffer({ label: 'CFD bc_type', size: cellCount * 4, usage: STORAGE_R }),
    params: device.createBuffer({
      label: 'CFD params',
      size: 48, // CfdParams struct size
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    }),
    forceAccum: device.createBuffer({
      label: 'CFD force_accum',
      size: 16, // 4 x f32 (Cl, Cd, Cm, _pad) — must be >= 16 for storage
      usage: STORAGE_RW | GPUBufferUsage.COPY_SRC,
    }),
    residualNorm: device.createBuffer({
      label: 'CFD residual_norm',
      size: 16, // 4 x f32
      usage: STORAGE_RW | GPUBufferUsage.COPY_SRC,
    }),
    readback: device.createBuffer({
      label: 'CFD readback',
      size: 32, // 8 x f32 (4 force_accum + 4 residual_norm)
      usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    }),
  };
}

export function destroyCfdBuffers(buffers: CfdBufferSet): void {
  for (const buf of Object.values(buffers)) {
    (buf as GPUBuffer).destroy();
  }
}
