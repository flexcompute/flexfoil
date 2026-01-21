/**
 * Velocity Texture Module
 * 
 * Computes and caches velocity field as a GPU texture.
 * All GPU visualizations sample from this texture.
 */

import { WebGPUContext } from './context';
import velocityShaderSource from './shaders/velocity.wgsl?raw';

export interface VelocityTextureOptions {
  /** Grid resolution width */
  width?: number;
  /** Grid resolution height */
  height?: number;
}

/*
 * VelocityParams - GPU uniform buffer layout:
 *   xMin, xMax, yMin, yMax (grid bounds)
 *   alpha, vInf (flow conditions)
 *   nx, ny (grid resolution)
 *   nPanels (panel count)
 */

/**
 * GPU-based velocity field texture.
 * 
 * Output format (RGBA32Float):
 *   R = u (x-velocity)
 *   G = v (y-velocity)
 *   B = |V| (speed)
 *   A = ψ (stream function)
 */
export class VelocityTexture {
  private ctx: WebGPUContext;
  private pipeline: GPUComputePipeline;
  private bindGroupLayout: GPUBindGroupLayout;
  
  private paramsBuffer: GPUBuffer;
  private nodesBuffer: GPUBuffer | null = null;
  private gammaBuffer: GPUBuffer | null = null;
  
  private _texture: GPUTexture;
  private _textureView: GPUTextureView;
  private _sampler: GPUSampler;
  
  private _width: number;
  private _height: number;
  private _bounds: [number, number, number, number] = [-1, 2, -1, 1];
  private _psi0: number = 0;
  private _psiMin: number = 0;
  private _psiMax: number = 0;

  get texture(): GPUTexture { return this._texture; }
  get textureView(): GPUTextureView { return this._textureView; }
  get sampler(): GPUSampler { return this._sampler; }
  get width(): number { return this._width; }
  get height(): number { return this._height; }
  get bounds(): [number, number, number, number] { return this._bounds; }
  get psi0(): number { return this._psi0; }
  get psiMin(): number { return this._psiMin; }
  get psiMax(): number { return this._psiMax; }

  constructor(ctx: WebGPUContext, options: VelocityTextureOptions = {}) {
    this.ctx = ctx;
    this._width = options.width ?? 256;
    this._height = options.height ?? 128;

    // Create shader module
    const shaderModule = ctx.createShaderModule(velocityShaderSource, 'velocity-shader');

    // Create bind group layout
    this.bindGroupLayout = ctx.device.createBindGroupLayout({
      label: 'velocity-bind-group-layout',
      entries: [
        {
          binding: 0,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: 'uniform' },
        },
        {
          binding: 1,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: 'read-only-storage' },
        },
        {
          binding: 2,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: 'read-only-storage' },
        },
        {
          binding: 3,
          visibility: GPUShaderStage.COMPUTE,
          storageTexture: {
            access: 'write-only',
            format: 'rgba32float',
          },
        },
      ],
    });

    // Create pipeline
    this.pipeline = ctx.device.createComputePipeline({
      label: 'velocity-pipeline',
      layout: ctx.device.createPipelineLayout({
        bindGroupLayouts: [this.bindGroupLayout],
      }),
      compute: {
        module: shaderModule,
        entryPoint: 'main',
      },
    });

    // Create params buffer (will be updated each frame)
    this.paramsBuffer = ctx.device.createBuffer({
      label: 'velocity-params',
      size: 48, // 12 floats * 4 bytes (padded for alignment)
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    // Create output texture
    this._texture = this.createTexture();
    this._textureView = this._texture.createView();

    // Create sampler - use linear filtering if available, otherwise nearest
    const canFilter = ctx.device.features.has('float32-filterable');
    this._sampler = ctx.device.createSampler({
      label: 'velocity-sampler',
      magFilter: canFilter ? 'linear' : 'nearest',
      minFilter: canFilter ? 'linear' : 'nearest',
      addressModeU: 'clamp-to-edge',
      addressModeV: 'clamp-to-edge',
    });
  }

  private createTexture(): GPUTexture {
    return this.ctx.device.createTexture({
      label: 'velocity-texture',
      size: [this._width, this._height],
      format: 'rgba32float',
      usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING,
    });
  }

  /**
   * Resize the velocity texture.
   */
  resize(width: number, height: number): void {
    if (width === this._width && height === this._height) return;

    this._texture.destroy();
    this._width = width;
    this._height = height;
    this._texture = this.createTexture();
    this._textureView = this._texture.createView();
  }

  /**
   * Update the velocity field with new geometry/flow conditions.
   * 
   * @param panels - Array of panel node positions [{x, y}, ...]
   * @param gamma - Vorticity values at each node
   * @param alpha - Angle of attack in degrees
   * @param bounds - [xMin, xMax, yMin, yMax] domain bounds
   * @param psi0 - Dividing streamline value (for contour rendering)
   */
  update(
    panels: { x: number; y: number }[],
    gamma: number[],
    alpha: number,
    bounds: [number, number, number, number],
    psi0: number = 0
  ): void {
    this._bounds = bounds;
    this._psi0 = psi0;

    const nPanels = panels.length;
    if (nPanels < 3 || gamma.length !== nPanels) {
      return;
    }

    // Update params - must match shader struct layout exactly:
    // struct Params {
    //   x_min, x_max, y_min, y_max: f32  (bytes 0-15)
    //   alpha, v_inf: f32                 (bytes 16-23)
    //   nx, ny, n_panels, _pad: u32       (bytes 24-39)
    // }
    const paramsData = new ArrayBuffer(48);
    const paramsFloat = new Float32Array(paramsData, 0, 6);
    const paramsUint = new Uint32Array(paramsData, 24, 4);
    
    paramsFloat[0] = bounds[0];  // x_min
    paramsFloat[1] = bounds[1];  // x_max
    paramsFloat[2] = bounds[2];  // y_min
    paramsFloat[3] = bounds[3];  // y_max
    paramsFloat[4] = alpha * Math.PI / 180;  // alpha (radians)
    paramsFloat[5] = 1.0;  // v_inf
    
    paramsUint[0] = this._width;   // nx
    paramsUint[1] = this._height;  // ny
    paramsUint[2] = nPanels;       // n_panels
    paramsUint[3] = 0;             // padding
    
    this.ctx.queue.writeBuffer(this.paramsBuffer, 0, paramsData);

    // Update nodes buffer
    const nodesData = new Float32Array(nPanels * 2);
    for (let i = 0; i < nPanels; i++) {
      nodesData[i * 2] = panels[i].x;
      nodesData[i * 2 + 1] = panels[i].y;
    }

    if (!this.nodesBuffer || this.nodesBuffer.size < nodesData.byteLength) {
      this.nodesBuffer?.destroy();
      this.nodesBuffer = this.ctx.device.createBuffer({
        label: 'velocity-nodes',
        size: Math.max(nodesData.byteLength, 256),
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
    }
    this.ctx.queue.writeBuffer(this.nodesBuffer, 0, nodesData);

    // Update gamma buffer
    const gammaData = new Float32Array(gamma);

    if (!this.gammaBuffer || this.gammaBuffer.size < gammaData.byteLength) {
      this.gammaBuffer?.destroy();
      this.gammaBuffer = this.ctx.device.createBuffer({
        label: 'velocity-gamma',
        size: Math.max(gammaData.byteLength, 256),
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
    }
    this.ctx.queue.writeBuffer(this.gammaBuffer, 0, gammaData);

    // Create bind group
    const bindGroup = this.ctx.device.createBindGroup({
      label: 'velocity-bind-group',
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.paramsBuffer } },
        { binding: 1, resource: { buffer: this.nodesBuffer } },
        { binding: 2, resource: { buffer: this.gammaBuffer } },
        { binding: 3, resource: this._textureView },
      ],
    });

    // Dispatch compute
    const encoder = this.ctx.device.createCommandEncoder({
      label: 'velocity-encoder',
    });

    const pass = encoder.beginComputePass({
      label: 'velocity-compute-pass',
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(
      Math.ceil(this._width / 16),
      Math.ceil(this._height / 16)
    );
    pass.end();

    this.ctx.submit([encoder.finish()]);
  }

  /**
   * Create a bind group for consumers (smoke, contours, etc.)
   */
  createConsumerBindGroup(layout: GPUBindGroupLayout, binding: number): GPUBindGroup {
    return this.ctx.device.createBindGroup({
      layout,
      entries: [
        {
          binding,
          resource: this._textureView,
        },
        {
          binding: binding + 1,
          resource: this._sampler,
        },
      ],
    });
  }

  /**
   * Set psi range after computation (typically from CPU analysis).
   */
  setPsiRange(psiMin: number, psiMax: number, psi0: number): void {
    this._psiMin = psiMin;
    this._psiMax = psiMax;
    this._psi0 = psi0;
  }

  /**
   * Clean up GPU resources.
   */
  destroy(): void {
    this._texture.destroy();
    this.paramsBuffer.destroy();
    this.nodesBuffer?.destroy();
    this.gammaBuffer?.destroy();
  }
}
