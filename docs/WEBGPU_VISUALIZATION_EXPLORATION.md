# WebGPU Visualization Exploration

## Current Architecture Summary

Your visualization stack uses **Rust/WASM for computation** + **Canvas 2D for rendering**:

| Feature | Computation | Rendering | Bottleneck |
|---------|-------------|-----------|------------|
| **Smoke** | WASM (RK2 per particle) | Overlay Canvas 2D | CPU loop + JS marshal every frame |
| **ψ Contours** | WASM (grid) + JS (marching squares) | Canvas 2D fill | Grid transfer + marching squares |
| **Streamlines** | WASM (RK4 integration) | Canvas 2D stroke | Serial integration |
| **Cp bars** | WASM (panel solve) | Canvas 2D | Minor |

### Current Data Flow (Smoke Example)
```
┌─────────────────────────────────────────────────────────────────────┐
│ Every Frame (~60Hz)                                                 │
├─────────────────────────────────────────────────────────────────────┤
│  WASM                         │   JS                    │  Canvas   │
│  ─────                        │   ──                    │  ──────   │
│  1. For each particle:        │                         │           │
│     - velocity_at() O(n_panels)│                        │           │
│     - RK2 step                │                         │           │
│  2. Copy to Float64Array      │→ 3. Copy to TypedArray  │           │
│                               │   4. For each particle: │→ 5. arc() │
│                               │      rotate, toCanvas() │   fill()  │
└─────────────────────────────────────────────────────────────────────┘
```

**Pain points:**
- `velocity_at()` does O(n_panels) per call — repeated for every particle every frame
- JS←→WASM marshalling overhead for particle arrays
- Canvas 2D `arc()` calls for each particle (not batched)

---

## WebGPU Opportunities (In Priority Order)

### 1. 🚀 Smoke Particles — **Highest Impact**

**Current cost:** `O(particles × panels)` per frame in WASM, plus data transfer

**WebGPU approach:**
```
┌────────────────────────────────────────────────────────────────────┐
│ Setup (once per geometry change)                                   │
├────────────────────────────────────────────────────────────────────┤
│  1. Upload gamma[], panels[] to GPU storage buffers                │
│  2. Create particle position buffer (GPU-resident)                 │
└────────────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────────┐
│ Per Frame                                                          │
├────────────────────────────────────────────────────────────────────┤
│  Compute Pass:                                                     │
│    - workgroup(64,1,1) × ceil(n_particles/64)                      │
│    - Each thread: read position → velocity_at → RK2 → write back   │
│                                                                    │
│  Render Pass:                                                      │
│    - Instanced point sprites                                       │
│    - Vertex shader reads from same position buffer                 │
│    - Fragment shader: soft circle + alpha fade                     │
└────────────────────────────────────────────────────────────────────┘
```

**Benefits:**
- **Zero CPU→GPU copy per frame** — particles live entirely on GPU
- **Parallel velocity evaluation** — 1000 particles computed in parallel
- **Instanced rendering** — single draw call for all particles

**Implementation sketch (compute shader):**
```wgsl
@group(0) @binding(0) var<storage, read_write> positions: array<vec2f>;
@group(0) @binding(1) var<storage, read_write> ages: array<f32>;
@group(0) @binding(2) var<storage, read> panels: array<vec4f>;  // (x1,y1,x2,y2)
@group(0) @binding(3) var<storage, read> gamma: array<f32>;
@group(0) @binding(4) var<uniform> params: Params;  // alpha, v_inf, dt

@compute @workgroup_size(64)
fn advect(@builtin(global_invocation_id) id: vec3u) {
    let idx = id.x;
    if (idx >= arrayLength(&positions)) { return; }
    
    var p = positions[idx];
    var age = ages[idx];
    
    // RK2 midpoint
    let v1 = velocity_at(p);
    let mid = p + 0.5 * params.dt * v1;
    let v2 = velocity_at(mid);
    
    p += params.dt * v2;
    age += params.dt;
    
    // Respawn if too old or inside airfoil
    if (age > params.max_age || is_inside(p)) {
        let spawn_idx = idx % spawn_count;
        p = spawn_positions[spawn_idx];
        age = 0.0;
    }
    
    positions[idx] = p;
    ages[idx] = age;
}
```

---

### 2. 🎨 Stream Function (ψ) Contours — **High Impact**

**Current cost:** WASM grid compute + JS marching squares + Canvas fill

**WebGPU approach — Two options:**

#### Option A: Fragment Shader Direct Visualization (Simpler, Great Results)
Skip marching squares entirely — just sample ψ in the fragment shader:

```wgsl
@fragment
fn fs(@location(0) uv: vec2f) -> @location(0) vec4f {
    // Sample velocity texture or compute inline
    let psi = compute_psi(uv);  // or texture sample
    
    // Color based on distance from ψ₀
    let t = (psi - psi_0) / psi_range;
    let color = mix(blue, red, step(0.0, t)) * abs(t);
    
    // Optional: contour lines via fwidth()
    let line_thickness = 0.02;
    let contour_interval = 0.1;
    let dist_to_line = abs(fract(psi / contour_interval + 0.5) - 0.5);
    let line = 1.0 - smoothstep(0.0, fwidth(psi) * 2.0, dist_to_line);
    
    return vec4f(color.rgb + line * 0.3, 0.8);
}
```

**Benefits:**
- No marching squares needed at all
- Smooth anti-aliased contour lines via `fwidth()`
- Resolution-independent — looks good at any zoom

#### Option B: Compute Shader Marching Squares (If you need exact polylines)
- Parallel marching squares on GPU
- Output to vertex buffer for line rendering

---

### 3. 📈 Streamlines — **Medium Impact, Enables Cool Stuff**

**Current approach works OK**, but WebGPU enables:

#### Option A: Parallel Integration
- Launch 50 workgroups, each integrates one streamline
- Faster initial computation when seeding many lines

#### Option B: LIC (Line Integral Convolution) — **Beautiful Flow Viz**
```wgsl
@fragment
fn lic_fs(@location(0) uv: vec2f) -> @location(0) vec4f {
    var color = 0.0;
    var p = uv;
    
    // Forward integration
    for (var i = 0; i < 20; i++) {
        let v = sample_velocity(p);
        p += normalize(v) * step_size;
        color += noise(p * noise_scale);
    }
    
    // Backward integration
    p = uv;
    for (var i = 0; i < 20; i++) {
        let v = sample_velocity(p);
        p -= normalize(v) * step_size;
        color += noise(p * noise_scale);
    }
    
    return vec4f(vec3f(color / 40.0), 1.0);
}
```

**Result:** Silky smooth flow visualization without discrete lines

---

### 4. 🎯 Shared Optimization: Velocity Texture

The **velocity field** is the common bottleneck. Solution:

```
┌─────────────────────────────────────────────────────────────────────┐
│ On geometry/alpha change:                                          │
│  1. Compute velocity at NxM grid (compute shader, parallel)        │
│  2. Store as RGBA texture: (u, v, |V|, ψ)                          │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│ Per frame:                                                         │
│  - Smoke advection: textureSample(velocity_tex, p) → cheap!        │
│  - Contours: textureSample(velocity_tex, uv).w → ψ value           │
│  - Streamlines: textureSample for RK4 → bilinear interpolation     │
└─────────────────────────────────────────────────────────────────────┘
```

**Implementation note:** The `velocity_at()` function in `velocity.rs` has this kernel:

```rust
// For each panel j:
let (u1, v1) = influence_from_panel_start(j);
let (u2, v2) = influence_from_panel_end(j);
u += gamma[j] * u1 + gamma[j+1] * u2;
v += gamma[j] * v1 + gamma[j+1] * v2;
```

This maps beautifully to a compute shader with shared memory for `gamma[]`.

---

## Implementation Roadmap

### Phase 1: Smoke on GPU (1-2 days)
1. Create WebGPU device/context alongside Canvas
2. Port `velocity_at` to WGSL (straightforward translation)
3. Particle buffers + compute shader
4. Instanced point rendering
5. Keep Canvas fallback for non-WebGPU browsers

### Phase 2: ψ Contours via Fragment Shader (1 day)
1. Full-screen quad with ψ evaluation in fragment
2. Color mapping + fwidth() contour lines
3. Alpha masking for airfoil interior

### Phase 3: Velocity Texture Cache (0.5 days)
1. Compute ψ/velocity grid on geometry change
2. Store as texture
3. Update smoke shader to sample texture instead of recomputing

### Phase 4: LIC Streamlines (Optional, 1 day)
1. Noise texture generation
2. LIC fragment shader
3. Blend with existing viz

---

## Browser Support & Fallback

```typescript
async function initVisualization() {
  if (navigator.gpu) {
    const adapter = await navigator.gpu.requestAdapter();
    if (adapter) {
      // WebGPU path
      return new WebGPUVisualization(adapter);
    }
  }
  // Fallback to current WASM + Canvas 2D
  return new CanvasVisualization();
}
```

**WebGPU support (2024):**
- Chrome 113+ ✅
- Edge 113+ ✅
- Firefox (behind flag, shipping soon)
- Safari 17+ ✅ (macOS Sonoma, iOS 17)

---

## Estimated Performance Gains

| Feature | Current | WebGPU | Speedup |
|---------|---------|--------|---------|
| Smoke (1000 particles) | ~8ms/frame | <1ms/frame | **8-10x** |
| ψ Contours (200×120 grid) | ~12ms compute + ~5ms draw | <2ms total | **8x** |
| Streamlines (50 lines) | ~15ms | ~5ms | **3x** |

**Total frame budget improvement:** From ~25ms (40fps cap) to ~8ms (120fps capable)

---

## Files to Create/Modify

```
flexfoil-ui/
├── src/
│   ├── lib/
│   │   ├── webgpu/
│   │   │   ├── context.ts         # Device/adapter initialization
│   │   │   ├── smoke.ts           # Smoke compute + render pipeline
│   │   │   ├── contours.ts        # ψ contour fragment shader
│   │   │   ├── velocity.ts        # Velocity texture generation
│   │   │   └── shaders/
│   │   │       ├── smoke.wgsl
│   │   │       ├── contours.wgsl
│   │   │       └── velocity.wgsl
│   │   └── wasm.ts                # Existing (kept as fallback)
│   └── components/
│       └── AirfoilCanvas.tsx      # Add WebGPU rendering path
```

---

## Next Steps

1. **Quick Win:** Start with smoke on GPU — most impactful, self-contained
2. **Validate:** Ensure visual parity with current implementation
3. **Iterate:** Add contours, then velocity cache
4. **Polish:** LIC streamlines for "wow factor"

Want me to start implementing Phase 1 (smoke on GPU)?
