/**
 * Two-pass LIC flow-field renderer (Three.js), driven by a RANS center-span slice
 * mesh from the rans backend. Ported from the flow360 web app (see licShaders.ts).
 *
 *   pass 1: render the slice mesh (per-vertex velocity `_cfvec`) → flow texture
 *   pass 2: render the mesh again, convolving the flow texture along the flow → LIC,
 *           colored by a per-vertex colormap of Mach.
 *
 * Renders one static frame per RANS point onto a canvas that sits *behind* the 2D
 * airfoil canvas (whose element fills occlude the body interiors).
 */
import * as THREE from 'three';
import { LIC_Dir_VS, LIC_Dir_FS, LIC_VS, LIC_FS } from './licShaders';

export interface FlowMesh {
  points: number[];          // flat [x0,y0,x1,y1,…] (UI coords)
  tris: number[];            // flat [i0,j0,k0,…]
  vel: number[];             // flat [u0,v0,…]
  mach: number[];
  machRange: [number, number];
}

/** Blue→cyan→green→yellow→red ramp for the Mach colormap. */
function colormap(t: number): [number, number, number] {
  t = Math.max(0, Math.min(1, t));
  const stops: Array<[number, [number, number, number]]> = [
    [0.0, [0.23, 0.30, 0.75]], [0.25, [0.0, 0.70, 0.90]], [0.5, [0.0, 0.80, 0.30]],
    [0.75, [0.95, 0.85, 0.10]], [1.0, [0.85, 0.15, 0.10]],
  ];
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [a, ca] = stops[i - 1];
      const [b, cb] = stops[i];
      const f = (t - a) / (b - a);
      return [ca[0] + (cb[0] - ca[0]) * f, ca[1] + (cb[1] - ca[1]) * f, ca[2] + (cb[2] - ca[2]) * f];
    }
  }
  return stops[stops.length - 1][1];
}

export class LicView {
  private renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  private camera = new THREE.OrthographicCamera(-1, 1, 1, -1, -1, 1);
  private flowRT: THREE.WebGLRenderTarget;
  private dirMat: THREE.ShaderMaterial;
  private licMat: THREE.ShaderMaterial;
  private mesh: THREE.Mesh | null = null;
  private w = 1;
  private h = 1;

  constructor(canvas: HTMLCanvasElement) {
    this.renderer = new THREE.WebGLRenderer({ canvas, alpha: true, antialias: true });
    this.renderer.setClearColor(0x000000, 0);
    const dummy = new THREE.DataTexture(new Uint8Array([0, 0, 0, 255]), 1, 1);
    dummy.needsUpdate = true;
    const common = { noiseTexture: { value: dummy }, textureRepeat: { value: 500 } };
    // DoubleSide: VTK triangle winding is arbitrary, so don't let backface culling
    // drop half (or all) of the slice mesh.
    this.dirMat = new THREE.ShaderMaterial({
      vertexShader: LIC_Dir_VS, fragmentShader: LIC_Dir_FS, uniforms: { ...common },
      side: THREE.DoubleSide,
    });
    this.licMat = new THREE.ShaderMaterial({
      vertexShader: LIC_VS, fragmentShader: LIC_FS, vertexColors: true,
      side: THREE.DoubleSide,
      uniforms: {
        ...common, flowTexture: { value: null as THREE.Texture | null },
        sizeK: { value: 1 }, passStep: { value: 0.003 },
      },
    });
    this.flowRT = new THREE.WebGLRenderTarget(2, 2, { minFilter: THREE.LinearFilter, magFilter: THREE.LinearFilter });
  }

  /** World rect (UI coords) covering the full pixel buffer — match the 2D canvas. */
  setCamera(left: number, right: number, bottom: number, top: number): void {
    this.camera.left = left;
    this.camera.right = right;
    this.camera.top = top;
    this.camera.bottom = bottom;
    this.camera.updateProjectionMatrix();
  }

  /** Tune the LIC look (streak length / noise frequency) and re-render. */
  setParams(p: { passStep?: number; textureRepeat?: number }): void {
    if (p.passStep !== undefined) this.licMat.uniforms.passStep.value = p.passStep;
    if (p.textureRepeat !== undefined) {
      this.dirMat.uniforms.textureRepeat.value = p.textureRepeat;
      this.licMat.uniforms.textureRepeat.value = p.textureRepeat;
    }
  }

  resize(wpx: number, hpx: number): void {
    if (wpx === this.w && hpx === this.h) return;
    this.w = wpx;
    this.h = hpx;
    this.renderer.setSize(wpx, hpx, false);
    this.flowRT.setSize(Math.max(2, wpx), Math.max(2, hpx));
    this.licMat.uniforms.sizeK.value = wpx / hpx;
  }

  setFlow(fm: FlowMesh): void {
    if (this.mesh) {
      this.scene.remove(this.mesh);
      this.mesh.geometry.dispose();
    }
    const n = fm.points.length / 2;
    const pos = new Float32Array(n * 3);
    const cfv = new Float32Array(n * 4);
    const col = new Float32Array(n * 3);
    const nor = new Float32Array(n * 3);
    const [m0, m1] = fm.machRange;
    const dm = m1 - m0 || 1;
    for (let i = 0; i < n; i++) {
      pos[i * 3] = fm.points[i * 2];
      pos[i * 3 + 1] = fm.points[i * 2 + 1];
      cfv[i * 4] = fm.vel[i * 2];
      cfv[i * 4 + 1] = fm.vel[i * 2 + 1];
      const [r, g, b] = colormap((fm.mach[i] - m0) / dm);
      col[i * 3] = r; col[i * 3 + 1] = g; col[i * 3 + 2] = b;
      nor[i * 3 + 2] = 1;   // flat slab → +z normal
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geo.setAttribute('_cfvec', new THREE.BufferAttribute(cfv, 4));
    geo.setAttribute('color', new THREE.BufferAttribute(col, 3));
    geo.setAttribute('normal', new THREE.BufferAttribute(nor, 3));
    geo.setIndex(fm.tris);
    this.mesh = new THREE.Mesh(geo, this.dirMat);
    this.scene.add(this.mesh);
  }

  render(): void {
    if (!this.mesh) {
      this.clear();
      return;
    }
    // pass 1: flow texture
    this.mesh.material = this.dirMat;
    this.renderer.setRenderTarget(this.flowRT);
    this.renderer.clear();
    this.renderer.render(this.scene, this.camera);
    // pass 2: LIC to the canvas
    this.licMat.uniforms.flowTexture.value = this.flowRT.texture;
    this.mesh.material = this.licMat;
    this.renderer.setRenderTarget(null);
    this.renderer.clear();
    this.renderer.render(this.scene, this.camera);
  }

  clear(): void {
    this.renderer.setRenderTarget(null);
    this.renderer.clear();
  }

  private colorMat?: THREE.MeshBasicMaterial;

  /** Debug: render the raw mesh flat-colored by Mach (no LIC) — isolates whether the
   *  geometry, camera and per-vertex colors are correct, separately from the shaders. */
  debugColors(): void {
    if (!this.mesh) return;
    if (!this.colorMat) this.colorMat = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
    this.mesh.material = this.colorMat;
    this.renderer.setRenderTarget(null);
    this.renderer.clear();
    this.renderer.render(this.scene, this.camera);
  }
}
