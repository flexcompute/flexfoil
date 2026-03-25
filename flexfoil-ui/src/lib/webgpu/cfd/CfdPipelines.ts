/**
 * Compute pipeline creation for the CFD solver.
 * Each shader is loaded from WGSL source, with cfd_common.wgsl prepended.
 */

import cfdCommonSrc from '../shaders/cfd/cfd_common.wgsl?raw';
import metricsSrc from '../shaders/cfd/metrics.wgsl?raw';
import primitivesSrc from '../shaders/cfd/primitives.wgsl?raw';
import reconstructMusclSrc from '../shaders/cfd/reconstruct_muscl.wgsl?raw';
import roeFluxSrc from '../shaders/cfd/roe_flux.wgsl?raw';
import residualSrc from '../shaders/cfd/residual.wgsl?raw';
import updateSrc from '../shaders/cfd/update.wgsl?raw';
import bcSrc from '../shaders/cfd/bc.wgsl?raw';
import forcesReduceSrc from '../shaders/cfd/forces_reduce.wgsl?raw';

import type { CfdBufferSet } from './CfdBuffers';

export interface CfdPipelineSet {
  metrics: GPUComputePipeline;
  primitives: GPUComputePipeline;
  reconstructXi: GPUComputePipeline;
  reconstructEta: GPUComputePipeline;
  roeFluxXi: GPUComputePipeline;
  roeFluxEta: GPUComputePipeline;
  residual: GPUComputePipeline;
  update: GPUComputePipeline;
  bc: GPUComputePipeline;
  forcesReduce: GPUComputePipeline;
}

export interface CfdBindGroupSet {
  metrics: GPUBindGroup;
  primitives: GPUBindGroup;
  reconstructXi: GPUBindGroup;
  reconstructEta: GPUBindGroup;
  roeFluxXi: GPUBindGroup;
  roeFluxEta: GPUBindGroup;
  residual: GPUBindGroup;
  update: GPUBindGroup;
  bc: GPUBindGroup;
  forcesReduce: GPUBindGroup;
}

function makeShader(device: GPUDevice, source: string, label: string): GPUShaderModule {
  return device.createShaderModule({
    label,
    code: cfdCommonSrc + '\n' + source,
  });
}

export function createCfdPipelines(device: GPUDevice): CfdPipelineSet {
  const auto = 'auto';

  const metricsModule = makeShader(device, metricsSrc, 'cfd_metrics');
  const primitivesModule = makeShader(device, primitivesSrc, 'cfd_primitives');
  const reconstructModule = makeShader(device, reconstructMusclSrc, 'cfd_reconstruct_muscl');
  const roeFluxModule = makeShader(device, roeFluxSrc, 'cfd_roe_flux');
  const residualModule = makeShader(device, residualSrc, 'cfd_residual');
  const updateModule = makeShader(device, updateSrc, 'cfd_update');
  const bcModule = makeShader(device, bcSrc, 'cfd_bc');
  const forcesModule = makeShader(device, forcesReduceSrc, 'cfd_forces_reduce');

  return {
    metrics: device.createComputePipeline({
      label: 'cfd_metrics_pipeline',
      layout: auto,
      compute: { module: metricsModule, entryPoint: 'compute_metrics' },
    }),
    primitives: device.createComputePipeline({
      label: 'cfd_primitives_pipeline',
      layout: auto,
      compute: { module: primitivesModule, entryPoint: 'conservative_to_primitive' },
    }),
    reconstructXi: device.createComputePipeline({
      label: 'cfd_reconstruct_xi_pipeline',
      layout: auto,
      compute: { module: reconstructModule, entryPoint: 'reconstruct_xi' },
    }),
    reconstructEta: device.createComputePipeline({
      label: 'cfd_reconstruct_eta_pipeline',
      layout: auto,
      compute: { module: reconstructModule, entryPoint: 'reconstruct_eta' },
    }),
    roeFluxXi: device.createComputePipeline({
      label: 'cfd_roe_flux_xi_pipeline',
      layout: auto,
      compute: { module: roeFluxModule, entryPoint: 'roe_flux_xi' },
    }),
    roeFluxEta: device.createComputePipeline({
      label: 'cfd_roe_flux_eta_pipeline',
      layout: auto,
      compute: { module: roeFluxModule, entryPoint: 'roe_flux_eta' },
    }),
    residual: device.createComputePipeline({
      label: 'cfd_residual_pipeline',
      layout: auto,
      compute: { module: residualModule, entryPoint: 'compute_residual' },
    }),
    update: device.createComputePipeline({
      label: 'cfd_update_pipeline',
      layout: auto,
      compute: { module: updateModule, entryPoint: 'update_solution' },
    }),
    bc: device.createComputePipeline({
      label: 'cfd_bc_pipeline',
      layout: auto,
      compute: { module: bcModule, entryPoint: 'apply_boundary_conditions' },
    }),
    forcesReduce: device.createComputePipeline({
      label: 'cfd_forces_reduce_pipeline',
      layout: auto,
      compute: { module: forcesModule, entryPoint: 'reduce_wall_forces' },
    }),
  };
}

export function createCfdBindGroups(
  device: GPUDevice,
  pipelines: CfdPipelineSet,
  buffers: CfdBufferSet,
): CfdBindGroupSet {
  return {
    metrics: device.createBindGroup({
      label: 'cfd_metrics_bg',
      layout: pipelines.metrics.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.meshX } },
        { binding: 2, resource: { buffer: buffers.meshY } },
        { binding: 3, resource: { buffer: buffers.metrics } },
      ],
    }),
    primitives: device.createBindGroup({
      label: 'cfd_primitives_bg',
      layout: pipelines.primitives.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.Q } },
        { binding: 2, resource: { buffer: buffers.W } },
      ],
    }),
    // Separate bind groups for xi/eta because auto-layout only includes
    // bindings actually referenced by each entry point.
    reconstructXi: device.createBindGroup({
      label: 'cfd_reconstruct_xi_bg',
      layout: pipelines.reconstructXi.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.W } },
        { binding: 2, resource: { buffer: buffers.leftXi } },
        { binding: 3, resource: { buffer: buffers.rightXi } },
      ],
    }),
    reconstructEta: device.createBindGroup({
      label: 'cfd_reconstruct_eta_bg',
      layout: pipelines.reconstructEta.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.W } },
        { binding: 4, resource: { buffer: buffers.leftEta } },
        { binding: 5, resource: { buffer: buffers.rightEta } },
      ],
    }),
    roeFluxXi: device.createBindGroup({
      label: 'cfd_roe_flux_xi_bg',
      layout: pipelines.roeFluxXi.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.metrics } },
        { binding: 2, resource: { buffer: buffers.leftXi } },
        { binding: 3, resource: { buffer: buffers.rightXi } },
        { binding: 6, resource: { buffer: buffers.fluxXi } },
      ],
    }),
    roeFluxEta: device.createBindGroup({
      label: 'cfd_roe_flux_eta_bg',
      layout: pipelines.roeFluxEta.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.metrics } },
        { binding: 4, resource: { buffer: buffers.leftEta } },
        { binding: 5, resource: { buffer: buffers.rightEta } },
        { binding: 7, resource: { buffer: buffers.fluxEta } },
      ],
    }),
    residual: device.createBindGroup({
      label: 'cfd_residual_bg',
      layout: pipelines.residual.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.metrics } },
        { binding: 2, resource: { buffer: buffers.fluxXi } },
        { binding: 3, resource: { buffer: buffers.fluxEta } },
        { binding: 4, resource: { buffer: buffers.residual } },
      ],
    }),
    update: device.createBindGroup({
      label: 'cfd_update_bg',
      layout: pipelines.update.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.residual } },
        { binding: 2, resource: { buffer: buffers.Q } },
        { binding: 3, resource: { buffer: buffers.bcType } },
      ],
    }),
    bc: device.createBindGroup({
      label: 'cfd_bc_bg',
      layout: pipelines.bc.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.Q } },
        { binding: 2, resource: { buffer: buffers.metrics } },
        { binding: 3, resource: { buffer: buffers.bcType } },
        { binding: 4, resource: { buffer: buffers.meshX } },
        { binding: 5, resource: { buffer: buffers.meshY } },
      ],
    }),
    forcesReduce: device.createBindGroup({
      label: 'cfd_forces_bg',
      layout: pipelines.forcesReduce.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: buffers.params } },
        { binding: 1, resource: { buffer: buffers.Q } },
        { binding: 2, resource: { buffer: buffers.metrics } },
        { binding: 3, resource: { buffer: buffers.meshX } },
        { binding: 4, resource: { buffer: buffers.meshY } },
        { binding: 5, resource: { buffer: buffers.residual } },
        { binding: 6, resource: { buffer: buffers.forceAccum } },
        { binding: 7, resource: { buffer: buffers.residualNorm } },
      ],
    }),
  };
}
