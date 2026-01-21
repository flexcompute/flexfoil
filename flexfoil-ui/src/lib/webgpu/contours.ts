/**
 * Stream Function Contour Renderer
 * 
 * Renders ψ contours using a full-screen fragment shader.
 * No marching squares needed - direct GPU evaluation.
 */

import { WebGPUContext } from './context';
import { VelocityTexture } from './velocity';
import contoursShaderSource from './shaders/contours.wgsl?raw';

export interface ContourOptions {
  /** Interval between contour lines */
  contourInterval?: number;
}

export class ContourRenderer {
  private ctx: WebGPUContext;
  private velocityTexture: VelocityTexture;
  
  private pipeline: GPURenderPipeline;
  private bindGroupLayout: GPUBindGroupLayout;
  private paramsBuffer: GPUBuffer;
  
  private _contourInterval: number;
  private _isDark: boolean = true;

  constructor(
    ctx: WebGPUContext,
    velocityTexture: VelocityTexture,
    options: ContourOptions = {}
  ) {
    this.ctx = ctx;
    this.velocityTexture = velocityTexture;
    this._contourInterval = options.contourInterval ?? 0.1;

    // Create shader module
    const shaderModule = ctx.createShaderModule(contoursShaderSource, 'contours-shader');

    // Check if float32 filtering is available
    const canFilter = ctx.device.features.has('float32-filterable');
    
    // Create bind group layout - use unfilterable-float if filtering not supported
    this.bindGroupLayout = ctx.device.createBindGroupLayout({
      label: 'contours-bind-group-layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
        { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: canFilter ? 'float' : 'unfilterable-float' } },
        { binding: 2, visibility: GPUShaderStage.FRAGMENT, sampler: { type: canFilter ? 'filtering' : 'non-filtering' } },
      ],
    });

    // Create pipeline
    this.pipeline = ctx.device.createRenderPipeline({
      label: 'contours-pipeline',
      layout: ctx.device.createPipelineLayout({
        bindGroupLayouts: [this.bindGroupLayout],
      }),
      vertex: {
        module: shaderModule,
        entryPoint: 'vs_fullscreen',
      },
      fragment: {
        module: shaderModule,
        entryPoint: 'fs_contour',
        targets: [{
          format: ctx.format,
          blend: {
            color: {
              srcFactor: 'src-alpha',
              dstFactor: 'one-minus-src-alpha',
              operation: 'add',
            },
            alpha: {
              srcFactor: 'one',
              dstFactor: 'one-minus-src-alpha',
              operation: 'add',
            },
          },
        }],
      },
      primitive: {
        topology: 'triangle-list',
      },
    });

    // Create params buffer (72 bytes aligned)
    this.paramsBuffer = ctx.device.createBuffer({
      label: 'contours-params',
      size: 80,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  /**
   * Set theme (affects colors).
   */
  setTheme(isDark: boolean): void {
    this._isDark = isDark;
  }

  /**
   * Set contour interval.
   */
  setContourInterval(interval: number): void {
    this._contourInterval = interval;
  }

  /**
   * Render contours.
   */
  render(
    encoder: GPUCommandEncoder,
    targetView: GPUTextureView,
    viewport: { centerX: number; centerY: number; zoom: number; width: number; height: number },
    alpha: number,
    clearFirst: boolean = true
  ): void {
    const bounds = this.velocityTexture.bounds;
    const aspect = viewport.width / viewport.height;
    const rad = -alpha * Math.PI / 180;
    
    // Compute psi range from velocity texture
    const psi0 = this.velocityTexture.psi0;
    const psiMin = this.velocityTexture.psiMin;
    const psiMax = this.velocityTexture.psiMax;
    
    // Update params buffer
    const params = new Float32Array([
      // tex bounds
      bounds[0], bounds[1], bounds[2], bounds[3],
      // psi range
      psi0, psiMin, psiMax,
      // contour settings
      this._contourInterval,
      // rotation
      Math.cos(rad), Math.sin(rad), 0.25, 0.0, // cos, sin, rot_cx, rot_cy
      // viewport
      viewport.centerX, viewport.centerY,
      viewport.zoom / Math.min(viewport.width, viewport.height) * 2,
      aspect,
      // theme
      this._isDark ? 1.0 : 0.0,
      0, 0, 0, // padding
    ]);
    
    this.ctx.queue.writeBuffer(this.paramsBuffer, 0, params);
    
    // Create bind group
    const bindGroup = this.ctx.device.createBindGroup({
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.paramsBuffer } },
        { binding: 1, resource: this.velocityTexture.textureView },
        { binding: 2, resource: this.velocityTexture.sampler },
      ],
    });
    
    // Render pass
    const pass = encoder.beginRenderPass({
      label: 'contours-render',
      colorAttachments: [{
        view: targetView,
        loadOp: clearFirst ? 'clear' : 'load',
        clearValue: this._isDark 
          ? { r: 0.059, g: 0.059, b: 0.059, a: 1 }  // #0f0f0f
          : { r: 1, g: 1, b: 1, a: 1 },
        storeOp: 'store',
      }],
    });
    
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.draw(3);  // Full-screen triangle
    pass.end();
  }

  /**
   * Clean up resources.
   */
  destroy(): void {
    this.paramsBuffer.destroy();
  }
}
