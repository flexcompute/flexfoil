/**
 * GPU Smoke Particle System
 * 
 * Fully GPU-resident particle system with:
 * - Compute shader for RK2 advection
 * - Instanced rendering for point sprites
 * - No CPU readback per frame
 */

import { WebGPUContext } from './context';
import { VelocityTexture } from './velocity';
import smokeShaderSource from './shaders/smoke.wgsl?raw';

export interface GPUSmokeOptions {
  /** Maximum number of particles */
  maxParticles?: number;
  /** Particles per blob */
  particlesPerBlob?: number;
  /** Number of spawn points */
  spawnCount?: number;
  /** Time between spawns (seconds) */
  spawnInterval?: number;
  /** Maximum particle age (seconds) */
  maxAge?: number;
}

/*
 * SimParams - GPU uniform buffer layout for advection compute shader:
 *   dt, maxAge, spawnInterval, timeSinceSpawn,
 *   texXMin, texXMax, texYMin, texYMax,
 *   psi0, particleCount, spawnCount
 *
 * RenderParams - GPU uniform buffer layout for particle render shader:
 *   viewCenterX, viewCenterY, zoom, aspect,
 *   cosAlpha, sinAlpha, rotCx, rotCy,
 *   colorAbove[3], colorBelow[3], psi0, pointSize, maxAge
 */

export class GPUSmokeSystem {
  private ctx: WebGPUContext;
  private velocityTexture: VelocityTexture;
  
  // Compute pipeline
  private computePipeline: GPUComputePipeline;
  private computeBindGroupLayout: GPUBindGroupLayout;
  
  // Render pipeline
  private renderPipeline: GPURenderPipeline;
  private renderBindGroupLayout: GPUBindGroupLayout;
  
  // Buffers
  private particleBuffer: GPUBuffer;
  private spawnPointsBuffer: GPUBuffer;
  private simParamsBuffer: GPUBuffer;
  private renderParamsBuffer: GPUBuffer;
  
  // State
  private _particleCount: number;
  private _spawnCount: number;
  private _maxAge: number;
  private _spawnInterval: number;
  private _totalTime: number = 0;  // Monotonically increasing global time
  private _psi0: number = 0;
  private _particlesPerBlob: number = 15;  // Particles per spawn point blob
  private _alpha: number = 0;  // Angle of attack (radians) for freestream direction
  
  // Stored spawn points for initialization (CPU-side copy)
  private _spawnPoints: { x: number; y: number }[] = [];
  
  // Colors (can be updated based on theme)
  private colorAbove: [number, number, number] = [1.0, 0.42, 0.42]; // Red-ish
  private colorBelow: [number, number, number] = [0.30, 0.67, 0.97]; // Blue-ish

  get particleCount(): number { return this._particleCount; }

  constructor(
    ctx: WebGPUContext,
    velocityTexture: VelocityTexture,
    options: GPUSmokeOptions = {}
  ) {
    this.ctx = ctx;
    this.velocityTexture = velocityTexture;
    
    this._particleCount = options.maxParticles ?? 2000;
    this._spawnCount = options.spawnCount ?? 30;
    this._spawnInterval = options.spawnInterval ?? 2.0;
    this._maxAge = options.maxAge ?? 6.0;
    
    // Create shader module
    const shaderModule = ctx.createShaderModule(smokeShaderSource, 'smoke-shader');
    
    // ========================================================================
    // Compute Pipeline
    // ========================================================================
    
    this.computeBindGroupLayout = ctx.device.createBindGroupLayout({
      label: 'smoke-compute-bind-group-layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
        { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
        { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
        { binding: 3, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float' } },
      ],
    });
    
    this.computePipeline = ctx.device.createComputePipeline({
      label: 'smoke-compute-pipeline',
      layout: ctx.device.createPipelineLayout({
        bindGroupLayouts: [this.computeBindGroupLayout],
      }),
      compute: {
        module: shaderModule,
        entryPoint: 'advect',
      },
    });
    
    // ========================================================================
    // Render Pipeline
    // ========================================================================
    
    this.renderBindGroupLayout = ctx.device.createBindGroupLayout({
      label: 'smoke-render-bind-group-layout',
      entries: [
        { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
        { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } },
      ],
    });
    
    this.renderPipeline = ctx.device.createRenderPipeline({
      label: 'smoke-render-pipeline',
      layout: ctx.device.createPipelineLayout({
        bindGroupLayouts: [this.renderBindGroupLayout],
      }),
      vertex: {
        module: shaderModule,
        entryPoint: 'vs_particle',
      },
      fragment: {
        module: shaderModule,
        entryPoint: 'fs_particle',
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
    
    // ========================================================================
    // Buffers
    // ========================================================================
    
    // Particle buffer: vec2f pos, f32 age, f32 psi = 16 bytes per particle
    this.particleBuffer = ctx.device.createBuffer({
      label: 'smoke-particles',
      size: this._particleCount * 16,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    
    // Spawn points buffer - initialize with default spawn line upstream
    this.spawnPointsBuffer = ctx.device.createBuffer({
      label: 'smoke-spawn-points',
      size: this._spawnCount * 8, // vec2f per spawn point
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
    
    // Initialize with default spawn points (vertical line upstream)
    const defaultSpawnPoints = new Float32Array(this._spawnCount * 2);
    for (let i = 0; i < this._spawnCount; i++) {
      const t = i / Math.max(1, this._spawnCount - 1);
      defaultSpawnPoints[i * 2] = -0.5;  // x: upstream of LE
      defaultSpawnPoints[i * 2 + 1] = -0.5 + t * 1.0;  // y: -0.5 to 0.5
    }
    ctx.queue.writeBuffer(this.spawnPointsBuffer, 0, defaultSpawnPoints);
    
    // Sim params buffer (64 bytes)
    this.simParamsBuffer = ctx.device.createBuffer({
      label: 'smoke-sim-params',
      size: 64,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    
    // Render params buffer (80 bytes, padded)
    this.renderParamsBuffer = ctx.device.createBuffer({
      label: 'smoke-render-params',
      size: 80,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    
    // Initialize particles to "dead" state (high age)
    this.reset();
  }

  /**
   * Set spawn points.
   */
  setSpawnPoints(points: { x: number; y: number }[]): void {
    this._spawnCount = points.length;
    
    // Resize buffer if needed
    const requiredSize = this._spawnCount * 8;
    if (this.spawnPointsBuffer.size < requiredSize) {
      this.spawnPointsBuffer.destroy();
      this.spawnPointsBuffer = this.ctx.device.createBuffer({
        label: 'smoke-spawn-points',
        size: requiredSize,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
    }
    
    // Store CPU-side copy for use in reset()
    this._spawnPoints = points.slice(0, this._spawnCount);
    
    const data = new Float32Array(this._spawnCount * 2);
    for (let i = 0; i < this._spawnCount; i++) {
      data[i * 2] = points[i].x;
      data[i * 2 + 1] = points[i].y;
    }
    this.ctx.queue.writeBuffer(this.spawnPointsBuffer, 0, data);
  }

  /**
   * Set colors (call when theme changes).
   */
  setColors(isDark: boolean): void {
    if (isDark) {
      this.colorAbove = [1.0, 0.42, 0.42];  // #ff6b6b
      this.colorBelow = [0.30, 0.67, 0.97]; // #4dabf7
    } else {
      this.colorAbove = [0.90, 0.22, 0.27]; // #e63946
      this.colorBelow = [0.10, 0.44, 0.76]; // #1971c2
    }
  }

  /**
   * Set simulation parameters.
   */
  setParams(spawnInterval: number, maxAge: number): void {
    this._spawnInterval = spawnInterval;
    this._maxAge = maxAge;
  }

  /**
   * Set dividing streamline value.
   */
  setPsi0(psi0: number): void {
    this._psi0 = psi0;
  }

  /**
   * Set angle of attack (for freestream direction in compute shader).
   */
  setAlpha(alphaDegrees: number): void {
    this._alpha = alphaDegrees * Math.PI / 180;
  }

  /**
   * Reset all particles and global time.
   * 
   * Particles are organized into WAVES. Each wave contains all spawn points.
   * Waves spawn on a fixed schedule (every spawn_interval seconds).
   * The shader calculates particle age from global time and wave assignment.
   */
  reset(): void {
    const spawnCount = Math.max(1, this._spawnPoints.length);
    const numWaves = Math.max(1, Math.ceil(this._maxAge / this._spawnInterval));
    
    // Calculate particles per blob to fill all waves evenly
    const totalBlobs = spawnCount * numWaves;
    this._particlesPerBlob = Math.max(1, Math.floor(this._particleCount / totalBlobs));
    const particlesPerWave = this._particlesPerBlob * spawnCount;
    
    const data = new Float32Array(this._particleCount * 4);
    
    for (let i = 0; i < this._particleCount; i++) {
      // Which wave does this particle belong to? (fixed assignment)
      const waveId = Math.floor(i / particlesPerWave) % numWaves;
      
      // Which spawn point within the wave?
      const withinWave = i % particlesPerWave;
      const spawnIdx = Math.floor(withinWave / this._particlesPerBlob) % spawnCount;
      const spawnPoint = this._spawnPoints[spawnIdx] ?? { x: -1, y: 0 };
      
      // Random position within blob (uniform disk distribution)
      const angle = Math.random() * Math.PI * 2;
      const r = Math.sqrt(Math.random()) * 0.02;  // Blob radius = 2% of chord
      
      // Initial age: wave 0 is youngest, wave N-1 is oldest
      // This matches what the shader will calculate at t=0
      const initialAge = waveId * this._spawnInterval;
      
      data[i * 4] = spawnPoint.x + Math.cos(angle) * r;
      data[i * 4 + 1] = spawnPoint.y + Math.sin(angle) * r;
      data[i * 4 + 2] = initialAge;
      data[i * 4 + 3] = 0; // psi
    }
    
    this.ctx.queue.writeBuffer(this.particleBuffer, 0, data);
    
    // Reset global time to 0
    this._totalTime = 0;
  }

  /**
   * Update particles (dispatch compute shader).
   */
  update(dt: number): void {
    // Accumulate global time (never resets except on manual reset)
    this._totalTime += dt;
    
    // Calculate number of waves that fit in max_age
    const numWaves = Math.max(1, Math.ceil(this._maxAge / this._spawnInterval));
    
    // Update sim params - 64 byte buffer
    // Layout: dt, maxAge, spawnInterval, totalTime,
    //         texBounds[4], psi0, cosAlpha, sinAlpha,
    //         particleCount, spawnCount, particlesPerBlob, numWaves
    const bounds = this.velocityTexture.bounds;
    const buffer = new ArrayBuffer(64);
    const floatView = new Float32Array(buffer);
    const uintView = new Uint32Array(buffer);
    
    floatView[0] = dt;
    floatView[1] = this._maxAge;
    floatView[2] = this._spawnInterval;
    floatView[3] = this._totalTime;
    floatView[4] = bounds[0];  // tex_x_min
    floatView[5] = bounds[1];  // tex_x_max
    floatView[6] = bounds[2];  // tex_y_min
    floatView[7] = bounds[3];  // tex_y_max
    floatView[8] = this._psi0;
    floatView[9] = Math.cos(this._alpha);   // cos_alpha for freestream
    floatView[10] = Math.sin(this._alpha);  // sin_alpha for freestream
    uintView[11] = this._particleCount;
    uintView[12] = this._spawnCount;
    uintView[13] = this._particlesPerBlob;
    uintView[14] = numWaves;
    uintView[15] = 0;  // padding
    
    this.ctx.queue.writeBuffer(this.simParamsBuffer, 0, buffer);
    
    // Create bind group
    const bindGroup = this.ctx.device.createBindGroup({
      layout: this.computeBindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.simParamsBuffer } },
        { binding: 1, resource: { buffer: this.particleBuffer } },
        { binding: 2, resource: { buffer: this.spawnPointsBuffer } },
        { binding: 3, resource: this.velocityTexture.textureView },
      ],
    });
    
    // Dispatch compute
    const encoder = this.ctx.device.createCommandEncoder({ label: 'smoke-compute' });
    const pass = encoder.beginComputePass({ label: 'smoke-advect' });
    pass.setPipeline(this.computePipeline);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(Math.ceil(this._particleCount / 64));
    pass.end();
    
    this.ctx.submit([encoder.finish()]);
  }

  /**
   * Render particles.
   */
  render(
    encoder: GPUCommandEncoder,
    targetView: GPUTextureView,
    viewport: { centerX: number; centerY: number; zoom: number; width: number; height: number },
    alpha: number
  ): void {
    const aspect = viewport.width / viewport.height;
    const rad = -alpha * Math.PI / 180;
    
    // Update render params
    const renderParams = new Float32Array([
      viewport.centerX,
      viewport.centerY,
      viewport.zoom / Math.min(viewport.width, viewport.height) * 2,
      aspect,
      Math.cos(rad),
      Math.sin(rad),
      0.25, // rot_cx (quarter chord)
      0.0,  // rot_cy
      ...this.colorAbove,
      0, // padding
      ...this.colorBelow,
      0, // padding
      this._psi0,
      0.008, // point size in world units (0.8% of chord)
      this._maxAge,
      0, // padding
    ]);
    
    this.ctx.queue.writeBuffer(this.renderParamsBuffer, 0, renderParams);
    
    // Create bind group
    const bindGroup = this.ctx.device.createBindGroup({
      layout: this.renderBindGroupLayout,
      entries: [
        { binding: 0, resource: { buffer: this.renderParamsBuffer } },
        { binding: 1, resource: { buffer: this.particleBuffer } },
      ],
    });
    
    // Render pass
    const pass = encoder.beginRenderPass({
      label: 'smoke-render',
      colorAttachments: [{
        view: targetView,
        loadOp: 'load',  // Preserve existing content
        storeOp: 'store',
      }],
    });
    
    pass.setPipeline(this.renderPipeline);
    pass.setBindGroup(0, bindGroup);
    // 6 vertices per quad (2 triangles), instanced for each particle
    pass.draw(6, this._particleCount);
    pass.end();
  }

  /**
   * Clean up resources.
   */
  destroy(): void {
    this.particleBuffer.destroy();
    this.spawnPointsBuffer.destroy();
    this.simParamsBuffer.destroy();
    this.renderParamsBuffer.destroy();
  }
}
