/**
 * WebGPU CFD Solver - Main orchestrator.
 *
 * Manages the GPU compute pipeline for structured-grid CFD.
 * Dispatches per-timestep compute passes and handles async readback.
 */

import type { CfdConfig } from './CfdConfig';
import { packCfdParams } from './CfdConfig';
import type { CfdBufferSet } from './CfdBuffers';
import { createCfdBuffers, destroyCfdBuffers } from './CfdBuffers';
import type { CfdPipelineSet, CfdBindGroupSet } from './CfdPipelines';
import {
  createCfdPipelines,
  createCfdBindGroups,
} from './CfdPipelines';

export interface CfdSolverResult {
  cl: number;
  cd: number;
  cm: number;
  residualL2: number;
  iteration: number;
}

export class CfdSolver {
  private device: GPUDevice;
  private config: CfdConfig;
  private buffers: CfdBufferSet;
  private pipelines: CfdPipelineSet;
  private bindGroups: CfdBindGroupSet;
  private iteration = 0;
  private _destroyed = false;

  private constructor(
    device: GPUDevice,
    config: CfdConfig,
    buffers: CfdBufferSet,
    pipelines: CfdPipelineSet,
    bindGroups: CfdBindGroupSet,
  ) {
    this.device = device;
    this.config = config;
    this.buffers = buffers;
    this.pipelines = pipelines;
    this.bindGroups = bindGroups;
  }

  /**
   * Initialize the CFD solver with mesh data and configuration.
   *
   * @param device - WebGPU device (shared with visualization renderer)
   * @param config - Solver configuration
   * @param meshX - Flat f32 array of mesh x-coordinates (ni*nj)
   * @param meshY - Flat f32 array of mesh y-coordinates (ni*nj)
   * @param initialQ - Flat f32 array of initial conservative variables (ni*nj*5)
   * @param bcTypes - Flat u32 array of boundary condition types (ni*nj)
   */
  static create(
    device: GPUDevice,
    config: CfdConfig,
    meshX: Float32Array,
    meshY: Float32Array,
    initialQ: Float32Array,
    bcTypes: Uint32Array,
  ): CfdSolver {
    const buffers = createCfdBuffers(device, config.ni, config.nj);
    const pipelines = createCfdPipelines(device);
    const bindGroups = createCfdBindGroups(device, pipelines, buffers);

    // Upload mesh data
    device.queue.writeBuffer(buffers.meshX, 0, meshX as Float32Array<ArrayBuffer>);
    device.queue.writeBuffer(buffers.meshY, 0, meshY as Float32Array<ArrayBuffer>);
    device.queue.writeBuffer(buffers.Q, 0, initialQ as Float32Array<ArrayBuffer>);
    device.queue.writeBuffer(buffers.bcType, 0, bcTypes as Uint32Array<ArrayBuffer>);

    // Upload initial params
    const dt = config.cfl * 0.01 / (1.0 + config.machInf);
    const paramsData = packCfdParams(config, dt, 0);
    device.queue.writeBuffer(buffers.params, 0, paramsData);

    const solver = new CfdSolver(device, config, buffers, pipelines, bindGroups);

    // Compute metrics (one-time)
    solver.computeMetrics();

    return solver;
  }

  /** Compute grid metrics (runs once at mesh setup). */
  private computeMetrics(): void {
    const { ni, nj } = this.config;
    const encoder = this.device.createCommandEncoder({ label: 'cfd_metrics' });

    const pass = encoder.beginComputePass({ label: 'metrics' });
    pass.setPipeline(this.pipelines.metrics);
    pass.setBindGroup(0, this.bindGroups.metrics);
    pass.dispatchWorkgroups(Math.ceil(ni / 16), Math.ceil(nj / 16));
    pass.end();

    this.device.queue.submit([encoder.finish()]);
  }

  /**
   * Run N timesteps of the solver.
   * Returns a promise that resolves with forces and residual after readback.
   */
  async step(nSteps: number): Promise<CfdSolverResult> {
    if (this._destroyed) throw new Error('Solver destroyed');

    const { ni, nj } = this.config;
    const dt = this.config.cfl * 0.01 / (1.0 + this.config.machInf);
    const wgXi = Math.ceil(ni / 16);
    const wgEta = Math.ceil(nj / 16);

    for (let s = 0; s < nSteps; s++) {
      this.iteration++;

      // Update params uniform before encoding
      const paramsData = packCfdParams(this.config, dt, this.iteration);
      this.device.queue.writeBuffer(this.buffers.params, 0, paramsData);

      // One encoder per step to avoid writeBuffer/encoder interleaving issues
      const encoder = this.device.createCommandEncoder({ label: `cfd_step_${this.iteration}` });

      // Pass 1: Q -> W (conservative to primitive)
      {
        const pass = encoder.beginComputePass({ label: 'primitives' });
        pass.setPipeline(this.pipelines.primitives);
        pass.setBindGroup(0, this.bindGroups.primitives);
        pass.dispatchWorkgroups(wgXi, wgEta);
        pass.end();
      }

      // Pass 2: Reconstruct xi (one workgroup per j-line)
      {
        const pass = encoder.beginComputePass({ label: 'reconstruct_xi' });
        pass.setPipeline(this.pipelines.reconstructXi);
        pass.setBindGroup(0, this.bindGroups.reconstructXi);
        pass.dispatchWorkgroups(nj);
        pass.end();
      }

      // Pass 3: Roe flux xi
      {
        const pass = encoder.beginComputePass({ label: 'roe_flux_xi' });
        pass.setPipeline(this.pipelines.roeFluxXi);
        pass.setBindGroup(0, this.bindGroups.roeFluxXi);
        pass.dispatchWorkgroups(wgXi, wgEta);
        pass.end();
      }

      // Pass 4: Reconstruct eta (one workgroup per i-column)
      {
        const pass = encoder.beginComputePass({ label: 'reconstruct_eta' });
        pass.setPipeline(this.pipelines.reconstructEta);
        pass.setBindGroup(0, this.bindGroups.reconstructEta);
        pass.dispatchWorkgroups(ni);
        pass.end();
      }

      // Pass 5: Roe flux eta
      {
        const pass = encoder.beginComputePass({ label: 'roe_flux_eta' });
        pass.setPipeline(this.pipelines.roeFluxEta);
        pass.setBindGroup(0, this.bindGroups.roeFluxEta);
        pass.dispatchWorkgroups(wgXi, wgEta);
        pass.end();
      }

      // Pass 6: Compute residual
      {
        const pass = encoder.beginComputePass({ label: 'residual' });
        pass.setPipeline(this.pipelines.residual);
        pass.setBindGroup(0, this.bindGroups.residual);
        pass.dispatchWorkgroups(wgXi, wgEta);
        pass.end();
      }

      // Pass 7: Update solution (explicit forward Euler)
      {
        const pass = encoder.beginComputePass({ label: 'update' });
        pass.setPipeline(this.pipelines.update);
        pass.setBindGroup(0, this.bindGroups.update);
        pass.dispatchWorkgroups(wgXi, wgEta);
        pass.end();
      }

      // Pass 8: Apply boundary conditions
      {
        const pass = encoder.beginComputePass({ label: 'bc' });
        pass.setPipeline(this.pipelines.bc);
        pass.setBindGroup(0, this.bindGroups.bc);
        pass.dispatchWorkgroups(Math.ceil(ni / 256));
        pass.end();
      }

      this.device.queue.submit([encoder.finish()]);
    }

    // Final pass: compute forces and residual norm
    {
      const encoder = this.device.createCommandEncoder({ label: 'cfd_forces' });

      const pass = encoder.beginComputePass({ label: 'forces' });
      pass.setPipeline(this.pipelines.forcesReduce);
      pass.setBindGroup(0, this.bindGroups.forcesReduce);
      pass.dispatchWorkgroups(1); // Single workgroup reduction
      pass.end();

      // Copy forces and residual to readback buffer
      encoder.copyBufferToBuffer(this.buffers.forceAccum, 0, this.buffers.readback, 0, 16);
      encoder.copyBufferToBuffer(this.buffers.residualNorm, 0, this.buffers.readback, 16, 16);

      this.device.queue.submit([encoder.finish()]);
    }

    // Async readback
    await this.buffers.readback.mapAsync(GPUMapMode.READ);
    const data = new Float32Array(this.buffers.readback.getMappedRange().slice(0));
    this.buffers.readback.unmap();

    return {
      cl: data[0],
      cd: data[1],
      cm: data[2],
      residualL2: data[4], // offset 4: first f32 of residualNorm (at readback byte 16)
      iteration: this.iteration,
    };
  }

  /** Get the current iteration count. */
  getIteration(): number {
    return this.iteration;
  }

  /** Get the Q buffer for visualization (e.g., to create velocity texture). */
  getQBuffer(): GPUBuffer {
    return this.buffers.Q;
  }

  /** Get mesh coordinate buffers for mesh visualization. */
  getMeshBuffers(): { x: GPUBuffer; y: GPUBuffer } {
    return { x: this.buffers.meshX, y: this.buffers.meshY };
  }

  /** Update configuration (e.g., CFL change mid-run). */
  updateConfig(config: Partial<CfdConfig>): void {
    Object.assign(this.config, config);
  }

  /** Reset solver to initial conditions. */
  reset(initialQ: Float32Array): void {
    this.iteration = 0;
    this.device.queue.writeBuffer(this.buffers.Q, 0, initialQ as Float32Array<ArrayBuffer>);
  }

  /** Clean up all GPU resources. */
  destroy(): void {
    if (this._destroyed) return;
    this._destroyed = true;
    destroyCfdBuffers(this.buffers);
  }
}
