/**
 * WebGPU Visualization Renderer
 * 
 * Main API for GPU-accelerated flow visualization.
 * Provides a unified interface for smoke, contours, and LIC streamlines.
 */

export { WebGPUContext, isWebGPUAvailable, getWebGPUInfo } from './context';
export { VelocityTexture } from './velocity';
export { GPUSmokeSystem } from './smoke';
export { ContourRenderer } from './contours';
export { LICRenderer } from './streamlines';

import { WebGPUContext, isWebGPUAvailable } from './context';
import { VelocityTexture } from './velocity';
import { GPUSmokeSystem } from './smoke';
import { ContourRenderer } from './contours';
import { LICRenderer } from './streamlines';

export interface RenderOptions {
  /** Show smoke particles */
  smoke?: boolean;
  /** Show ψ contours */
  contours?: boolean;
  /** Show LIC streamlines */
  streamlines?: boolean;
  /** Angle of attack in degrees */
  alpha?: number;
  /** Theme */
  isDark?: boolean;
}

export interface ViewportParams {
  centerX: number;
  centerY: number;
  zoom: number;
  width: number;
  height: number;
}

export interface PerfMetrics {
  frameTime: number;
  avgFrameTime: number;
  particleCount: number;
  fps: number;
}

/**
 * Main WebGPU renderer for flow visualization.
 */
export class WebGPURenderer {
  private ctx: WebGPUContext;
  private velocity: VelocityTexture;
  private smoke: GPUSmokeSystem;
  private contours: ContourRenderer;
  private lic: LICRenderer;
  
  // Performance tracking
  private frameTimes: number[] = [];
  private lastFrameTime: number = 0;
  private frameCount: number = 0;
  
  // State
  private _initialized: boolean = false;
  private _currentAlpha: number = 0;

  private constructor(ctx: WebGPUContext) {
    this.ctx = ctx;
    
    // Create velocity texture
    this.velocity = new VelocityTexture(ctx, { width: 256, height: 128 });
    
    // Create renderers
    this.smoke = new GPUSmokeSystem(ctx, this.velocity);
    this.contours = new ContourRenderer(ctx, this.velocity);
    this.lic = new LICRenderer(ctx, this.velocity);
    
    this._initialized = true;
  }

  /**
   * Initialize WebGPU renderer.
   * Returns null if WebGPU is not available.
   */
  static async create(canvas: HTMLCanvasElement): Promise<WebGPURenderer | null> {
    if (!isWebGPUAvailable()) {
      console.warn('WebGPU not available');
      return null;
    }

    const ctx = await WebGPUContext.create({ canvas, alphaMode: 'premultiplied' });
    if (!ctx) {
      return null;
    }

    return new WebGPURenderer(ctx);
  }

  get initialized(): boolean {
    return this._initialized;
  }

  /**
   * Update geometry (call when panels change).
   */
  updateGeometry(
    panels: { x: number; y: number }[],
    gamma: number[],
    alpha: number,
    bounds: [number, number, number, number],
    psi0: number = 0,
    psiMin: number = -1,
    psiMax: number = 1
  ): void {
    this._currentAlpha = alpha;
    
    // Update velocity texture
    this.velocity.update(panels, gamma, alpha, bounds, psi0);
    this.velocity.setPsiRange(psiMin, psiMax, psi0);
    
    // Update smoke with new psi0 and alpha (for freestream direction)
    this.smoke.setPsi0(psi0);
    this.smoke.setAlpha(alpha);
  }

  /**
   * Set smoke spawn points.
   */
  setSpawnPoints(points: { x: number; y: number }[]): void {
    this.smoke.setSpawnPoints(points);
  }

  /**
   * Set smoke parameters.
   */
  setSmokeParams(spawnInterval: number, maxAge: number): void {
    this.smoke.setParams(spawnInterval, maxAge);
  }

  /**
   * Reset smoke particles.
   */
  resetSmoke(): void {
    this.smoke.reset();
  }

  /**
   * Set theme for all renderers.
   */
  setTheme(isDark: boolean): void {
    this.smoke.setColors(isDark);
    this.contours.setTheme(isDark);
    this.lic.setTheme(isDark);
  }

  /**
   * Set LIC parameters (no-op for compatibility).
   */
  setLICParams(options: { numSteps?: number; stepSize?: number; kernelLength?: number; contrast?: number }): void {
    this.lic.setParams(options);
  }

  /**
   * Set pre-computed streamline data for GPU rendering.
   * @param streamlines Array of streamlines, each is an array of [x, y] points
   */
  setStreamlines(streamlines: [number, number][][]): void {
    this.lic.setStreamlines(streamlines);
  }

  /**
   * Set contour interval.
   */
  setContourInterval(interval: number): void {
    this.contours.setContourInterval(interval);
  }

  /**
   * Update smoke simulation (call every frame when smoke is visible).
   */
  updateSmoke(dt: number): void {
    this.smoke.update(dt);
  }

  /**
   * Render a frame.
   * The GPU canvas overlays the main canvas, so we use transparent clears
   * to allow the airfoil (drawn on main canvas) to show through.
   */
  render(viewport: ViewportParams, options: RenderOptions = {}): void {
    const startTime = performance.now();
    
    const alpha = options.alpha ?? this._currentAlpha;
    
    // Get render target
    const targetView = this.ctx.getCurrentTextureView();
    
    // Create command encoder
    const encoder = this.ctx.device.createCommandEncoder({ label: 'main-render' });
    
    // Determine what to render and in what order
    const showContours = options.contours ?? false;
    const showSmoke = options.smoke ?? false;
    const showStreamlines = options.streamlines ?? false;
    
    // Always clear to transparent first - this allows the main canvas to show through
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: targetView,
        loadOp: 'clear',
        clearValue: { r: 0, g: 0, b: 0, a: 0 },  // Transparent
        storeOp: 'store',
      }],
    });
    pass.end();
    
    // 1. Contours (drawn first as background - semi-transparent overlay)
    if (showContours) {
      this.contours.render(encoder, targetView, viewport, alpha, false);
    }
    
    // 2. LIC Streamlines (blend on top)
    if (showStreamlines) {
      this.lic.render(encoder, targetView, viewport, alpha);
    }
    
    // 3. Smoke particles (drawn last, on top)
    if (showSmoke) {
      this.smoke.render(encoder, targetView, viewport, alpha);
    }
    
    // Submit
    this.ctx.submit([encoder.finish()]);
    
    // Track performance
    const endTime = performance.now();
    const frameTime = endTime - startTime;
    this.lastFrameTime = frameTime;
    this.frameTimes.push(frameTime);
    if (this.frameTimes.length > 60) {
      this.frameTimes.shift();
    }
    this.frameCount++;
  }

  /**
   * Get performance metrics.
   */
  getPerformanceMetrics(): PerfMetrics {
    const avgFrameTime = this.frameTimes.length > 0
      ? this.frameTimes.reduce((a, b) => a + b, 0) / this.frameTimes.length
      : 0;
    
    return {
      frameTime: this.lastFrameTime,
      avgFrameTime,
      particleCount: this.smoke.particleCount,
      fps: avgFrameTime > 0 ? 1000 / avgFrameTime : 0,
    };
  }

  /**
   * Resize the canvas.
   */
  resize(width: number, height: number): void {
    this.ctx.resize(width, height);
  }

  /**
   * Clean up all resources.
   */
  destroy(): void {
    this.smoke.destroy();
    this.contours.destroy();
    this.lic.destroy();
    this.velocity.destroy();
    this.ctx.destroy();
    this._initialized = false;
  }
}

/**
 * Check if WebGPU is available and get info.
 */
export async function checkWebGPUSupport(): Promise<{
  available: boolean;
  info?: string;
}> {
  if (!isWebGPUAvailable()) {
    return { available: false, info: 'WebGPU not supported in this browser' };
  }

  try {
    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) {
      return { available: false, info: 'No WebGPU adapter found' };
    }

    const info = adapter.info;
    return {
      available: true,
      info: info ? `${info.vendor} ${info.architecture}` : 'WebGPU available',
    };
  } catch (e) {
    return { available: false, info: String(e) };
  }
}
