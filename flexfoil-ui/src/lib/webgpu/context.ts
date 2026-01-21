/**
 * WebGPU Context Management
 * 
 * Handles device/adapter initialization, feature detection, and canvas setup.
 */

export interface WebGPUContextOptions {
  /** Canvas element to render to */
  canvas: HTMLCanvasElement;
  /** Preferred texture format (defaults to navigator.gpu.getPreferredCanvasFormat()) */
  format?: GPUTextureFormat;
  /** Alpha mode for canvas */
  alphaMode?: GPUCanvasAlphaMode;
}

export class WebGPUContext {
  readonly device: GPUDevice;
  readonly queue: GPUQueue;
  readonly canvasContext: GPUCanvasContext;
  readonly format: GPUTextureFormat;
  readonly canvas: HTMLCanvasElement;

  private constructor(
    device: GPUDevice,
    canvasContext: GPUCanvasContext,
    format: GPUTextureFormat,
    canvas: HTMLCanvasElement
  ) {
    this.device = device;
    this.queue = device.queue;
    this.canvasContext = canvasContext;
    this.format = format;
    this.canvas = canvas;
  }

  /**
   * Initialize WebGPU context.
   * Returns null if WebGPU is not available.
   */
  static async create(options: WebGPUContextOptions): Promise<WebGPUContext | null> {
    // Check for WebGPU support
    if (!navigator.gpu) {
      console.warn('WebGPU not supported in this browser');
      return null;
    }

    // Request adapter
    const adapter = await navigator.gpu.requestAdapter({
      powerPreference: 'high-performance',
    });

    if (!adapter) {
      console.warn('No WebGPU adapter found');
      return null;
    }

    // Check for optional features
    const features: GPUFeatureName[] = [];
    if (adapter.features.has('float32-filterable')) {
      features.push('float32-filterable');
    }

    // Request device
    const device = await adapter.requestDevice({
      requiredFeatures: features,
      requiredLimits: {
        maxStorageBufferBindingSize: adapter.limits.maxStorageBufferBindingSize,
        maxComputeWorkgroupsPerDimension: adapter.limits.maxComputeWorkgroupsPerDimension,
      },
    });

    // Handle device lost
    device.lost.then((info) => {
      console.error('WebGPU device lost:', info.message);
      if (info.reason !== 'destroyed') {
        // Could attempt to recreate here
        console.warn('Device lost unexpectedly, may need to reinitialize');
      }
    });

    // Get canvas context
    const canvasContext = options.canvas.getContext('webgpu');
    if (!canvasContext) {
      console.warn('Could not get WebGPU canvas context');
      return null;
    }

    // Configure canvas
    const format = options.format ?? navigator.gpu.getPreferredCanvasFormat();
    canvasContext.configure({
      device,
      format,
      alphaMode: options.alphaMode ?? 'premultiplied',
    });

    return new WebGPUContext(device, canvasContext, format, options.canvas);
  }

  /**
   * Get current render target texture.
   */
  getCurrentTexture(): GPUTexture {
    return this.canvasContext.getCurrentTexture();
  }

  /**
   * Get current render target view.
   */
  getCurrentTextureView(): GPUTextureView {
    return this.getCurrentTexture().createView();
  }

  /**
   * Resize the canvas and reconfigure context.
   */
  resize(width: number, height: number): void {
    this.canvas.width = width;
    this.canvas.height = height;
    // Context is automatically reconfigured on next getCurrentTexture()
  }

  /**
   * Create a buffer with data.
   */
  createBuffer(
    data: Float32Array | Uint32Array | Uint16Array,
    usage: GPUBufferUsageFlags,
    label?: string
  ): GPUBuffer {
    const buffer = this.device.createBuffer({
      label,
      size: data.byteLength,
      usage: usage | GPUBufferUsage.COPY_DST,
      mappedAtCreation: true,
    });

    const mapped = buffer.getMappedRange();
    if (data instanceof Float32Array) {
      new Float32Array(mapped).set(data);
    } else if (data instanceof Uint32Array) {
      new Uint32Array(mapped).set(data);
    } else {
      new Uint16Array(mapped).set(data);
    }
    buffer.unmap();

    return buffer;
  }

  /**
   * Create an empty buffer.
   */
  createEmptyBuffer(size: number, usage: GPUBufferUsageFlags, label?: string): GPUBuffer {
    return this.device.createBuffer({
      label,
      size,
      usage,
    });
  }

  /**
   * Update buffer data.
   */
  writeBuffer(buffer: GPUBuffer, data: Float32Array | Uint32Array, offset = 0): void {
    this.queue.writeBuffer(buffer, offset, data.buffer, data.byteOffset, data.byteLength);
  }

  /**
   * Create a shader module from WGSL source.
   */
  createShaderModule(code: string, label?: string): GPUShaderModule {
    return this.device.createShaderModule({
      label,
      code,
    });
  }

  /**
   * Submit command buffers.
   */
  submit(commandBuffers: GPUCommandBuffer[]): void {
    this.queue.submit(commandBuffers);
  }

  /**
   * Clean up resources.
   */
  destroy(): void {
    this.device.destroy();
  }
}

/**
 * Check if WebGPU is available in this browser.
 */
export function isWebGPUAvailable(): boolean {
  return typeof navigator !== 'undefined' && !!navigator.gpu;
}

/**
 * Get WebGPU adapter info for display.
 */
export async function getWebGPUInfo(): Promise<{
  available: boolean;
  vendor?: string;
  architecture?: string;
  device?: string;
} | null> {
  if (!navigator.gpu) {
    return { available: false };
  }

  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) {
    return { available: false };
  }

  // Note: adapterInfo is available directly on the adapter
  try {
    const info = adapter.info;
    return {
      available: true,
      vendor: info?.vendor,
      architecture: info?.architecture,
      device: info?.device,
    };
  } catch {
    return { available: true };
  }
}
