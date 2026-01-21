/**
 * GPU Streamline Renderer
 * 
 * Renders pre-computed streamlines as line primitives.
 * The streamline computation is done on CPU (via WASM), and this
 * renderer just handles efficient GPU drawing.
 */

import { WebGPUContext } from './context';
import streamlineShaderSource from './shaders/streamlines.wgsl?raw';

export interface LICOptions {
  /** Not used - kept for API compatibility */
  numSteps?: number;
  stepSize?: number;
  kernelLength?: number;
  contrast?: number;
}

export class LICRenderer {
  private ctx: WebGPUContext;
  
  private pipeline: GPURenderPipeline;
  private bindGroupLayout: GPUBindGroupLayout;
  private paramsBuffer: GPUBuffer;
  
  // Streamline vertex data
  private vertexBuffer: GPUBuffer | null = null;
  private vertexCount: number = 0;
  private segmentCounts: number[] = []; // Vertices per streamline for line-strip rendering
  
  private _isDark: boolean = true;
  private _lineAlpha: number = 0.6;

  constructor(
    ctx: WebGPUContext,
    _velocityTexture: unknown, // Not used, kept for API compatibility
    _options: LICOptions = {}
  ) {
    this.ctx = ctx;

    const shaderModule = ctx.createShaderModule(streamlineShaderSource, 'streamline-shader');
    
    this.bindGroupLayout = ctx.device.createBindGroupLayout({
      label: 'streamline-bind-group-layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
      ],
    });

    this.pipeline = ctx.device.createRenderPipeline({
      label: 'streamline-pipeline',
      layout: ctx.device.createPipelineLayout({
        bindGroupLayouts: [this.bindGroupLayout],
      }),
      vertex: {
        module: shaderModule,
        entryPoint: 'vs_streamline',
        buffers: [{
          arrayStride: 8, // 2 floats * 4 bytes
          attributes: [{
            shaderLocation: 0,
            offset: 0,
            format: 'float32x2',
          }],
        }],
      },
      fragment: {
        module: shaderModule,
        entryPoint: 'fs_streamline',
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
        topology: 'line-strip',
      },
    });

    // Params buffer: 32 bytes
    this.paramsBuffer = ctx.device.createBuffer({
      label: 'streamline-params',
      size: 32,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
  }

  /**
   * Update streamline data from pre-computed results.
   * @param streamlines Array of streamlines, each is an array of [x, y] points
   */
  setStreamlines(streamlines: [number, number][][]): void {
    // Count total vertices
    let totalVertices = 0;
    this.segmentCounts = [];
    
    for (const line of streamlines) {
      if (line.length >= 2) {
        totalVertices += line.length;
        this.segmentCounts.push(line.length);
      }
    }
    
    if (totalVertices === 0) {
      this.vertexCount = 0;
      return;
    }
    
    // Create vertex data
    const vertexData = new Float32Array(totalVertices * 2);
    let offset = 0;
    
    for (const line of streamlines) {
      if (line.length >= 2) {
        for (const [x, y] of line) {
          vertexData[offset++] = x;
          vertexData[offset++] = y;
        }
      }
    }
    
    // Create or resize vertex buffer
    if (this.vertexBuffer) {
      this.vertexBuffer.destroy();
    }
    
    this.vertexBuffer = this.ctx.device.createBuffer({
      label: 'streamline-vertices',
      size: vertexData.byteLength,
      usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
    });
    
    this.ctx.queue.writeBuffer(this.vertexBuffer, 0, vertexData);
    this.vertexCount = totalVertices;
  }

  setTheme(isDark: boolean): void {
    this._isDark = isDark;
  }

  setParams(_options: Partial<LICOptions>): void {
    // No-op for compatibility
  }

  render(
    encoder: GPUCommandEncoder,
    targetView: GPUTextureView,
    viewport: { centerX: number; centerY: number; zoom: number; width: number; height: number },
    _alpha: number // Rotation already applied to streamline data
  ): void {
    if (!this.vertexBuffer || this.vertexCount === 0) {
      return;
    }
    
    const aspect = viewport.width / viewport.height;
    
    // Params: view_center_x, view_center_y, zoom, aspect, r, g, b, alpha
    const params = new Float32Array(8);
    params[0] = viewport.centerX;
    params[1] = viewport.centerY;
    params[2] = viewport.zoom / Math.min(viewport.width, viewport.height) * 2;
    params[3] = aspect;
    
    // Line color
    if (this._isDark) {
      params[4] = 0.5;  // R
      params[5] = 0.6;  // G
      params[6] = 0.7;  // B
    } else {
      params[4] = 0.3;  // R
      params[5] = 0.4;  // G
      params[6] = 0.5;  // B
    }
    params[7] = this._lineAlpha;
    
    this.ctx.queue.writeBuffer(this.paramsBuffer, 0, params);
    
    const bindGroup = this.ctx.device.createBindGroup({
      layout: this.bindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.paramsBuffer } },
      ],
    });
    
    const pass = encoder.beginRenderPass({
      label: 'streamline-render',
      colorAttachments: [{
        view: targetView,
        loadOp: 'load',
        storeOp: 'store',
      }],
    });
    
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, bindGroup);
    pass.setVertexBuffer(0, this.vertexBuffer);
    
    // Draw each streamline as a separate line-strip
    let baseVertex = 0;
    for (const count of this.segmentCounts) {
      pass.draw(count, 1, baseVertex, 0);
      baseVertex += count;
    }
    
    pass.end();
  }

  destroy(): void {
    this.paramsBuffer.destroy();
    if (this.vertexBuffer) {
      this.vertexBuffer.destroy();
    }
  }
}
