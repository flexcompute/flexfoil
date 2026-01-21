/**
 * Smoke Particle Shader
 * 
 * Compute shader: Advects particles using RK2 integration, sampling velocity texture.
 * Render shaders: Instanced point sprites with soft circle and alpha fade.
 */

// ============================================================================
// Shared Structures
// ============================================================================

struct Particle {
    pos: vec2f,
    age: f32,
    psi: f32,  // Stream function value (for coloring)
}

struct SimParams {
    dt: f32,
    max_age: f32,
    spawn_interval: f32,
    total_time: f32,      // Monotonically increasing global time
    // Texture bounds for UV mapping
    tex_x_min: f32,
    tex_x_max: f32,
    tex_y_min: f32,
    tex_y_max: f32,
    // Flow params
    psi_0: f32,          // Dividing streamline value
    cos_alpha: f32,      // Freestream direction (body-fixed)
    sin_alpha: f32,      // Freestream direction (body-fixed)
    particle_count: u32,
    spawn_count: u32,
    particles_per_blob: u32,  // Particles per spawn point blob (for grouping)
    num_waves: u32,      // Number of waves that fit in max_age
}

struct RenderParams {
    // Transform: airfoil coords -> clip space
    view_center_x: f32,
    view_center_y: f32,
    zoom: f32,
    aspect: f32,
    // Rotation for alpha display
    cos_alpha: f32,
    sin_alpha: f32,
    // Quarter chord (rotation center)
    rot_cx: f32,
    rot_cy: f32,
    // Colors
    color_above_r: f32,
    color_above_g: f32,
    color_above_b: f32,
    _pad1: f32,
    color_below_r: f32,
    color_below_g: f32,
    color_below_b: f32,
    _pad2: f32,
    // Dividing streamline
    psi_0: f32,
    point_size: f32,
    _pad3: vec2f,
}

// ============================================================================
// Compute Shader: Particle Advection
// ============================================================================

@group(0) @binding(0) var<uniform> sim_params: SimParams;
@group(0) @binding(1) var<storage, read_write> particles: array<Particle>;
@group(0) @binding(2) var<storage, read> spawn_points: array<vec2f>;
@group(0) @binding(3) var velocity_texture: texture_2d<f32>;

// Get texture dimensions
fn get_tex_dims() -> vec2u {
    return textureDimensions(velocity_texture);
}

// Map world coords to texture pixel coords
fn world_to_texel(pos: vec2f) -> vec2i {
    let dims = vec2f(get_tex_dims());
    let u = (pos.x - sim_params.tex_x_min) / (sim_params.tex_x_max - sim_params.tex_x_min);
    let v = (pos.y - sim_params.tex_y_min) / (sim_params.tex_y_max - sim_params.tex_y_min);
    return vec2i(i32(u * dims.x), i32(v * dims.y));
}

// Check if position is in bounds
fn is_in_bounds(pos: vec2f) -> bool {
    let u = (pos.x - sim_params.tex_x_min) / (sim_params.tex_x_max - sim_params.tex_x_min);
    let v = (pos.y - sim_params.tex_y_min) / (sim_params.tex_y_max - sim_params.tex_y_min);
    return u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0;
}

// Sample velocity from texture using textureLoad (works in compute shader)
fn sample_velocity(pos: vec2f) -> vec2f {
    if (!is_in_bounds(pos)) {
        // Freestream in body-fixed coordinates (comes in at angle alpha)
        return vec2f(sim_params.cos_alpha, sim_params.sin_alpha);
    }
    let texel = world_to_texel(pos);
    let dims = vec2i(get_tex_dims());
    let clamped = clamp(texel, vec2i(0), dims - vec2i(1));
    let sample = textureLoad(velocity_texture, clamped, 0);
    return vec2f(sample.r, sample.g);
}

// Sample stream function from texture
fn sample_psi(pos: vec2f) -> f32 {
    if (!is_in_bounds(pos)) {
        return 0.0;
    }
    let texel = world_to_texel(pos);
    let dims = vec2i(get_tex_dims());
    let clamped = clamp(texel, vec2i(0), dims - vec2i(1));
    let sample = textureLoad(velocity_texture, clamped, 0);
    return sample.a;  // Stream function is in alpha channel
}

// Simple hash for pseudo-random
fn hash(n: u32) -> f32 {
    var x = n;
    x = ((x >> 16u) ^ x) * 0x45d9f3bu;
    x = ((x >> 16u) ^ x) * 0x45d9f3bu;
    x = (x >> 16u) ^ x;
    return f32(x) / f32(0xffffffffu);
}

// Better hash with 2 inputs
fn hash2(a: u32, b: u32) -> f32 {
    var x = (a * 1597334677u) ^ (b * 3812015801u);
    x = ((x >> 16u) ^ x) * 0x45d9f3bu;
    x = ((x >> 16u) ^ x) * 0x45d9f3bu;
    x = (x >> 16u) ^ x;
    return f32(x) / f32(0xffffffffu);
}

@compute @workgroup_size(64)
fn advect(@builtin(global_invocation_id) id: vec3u) {
    let idx = id.x;
    if (idx >= sim_params.particle_count) {
        return;
    }

    var p = particles[idx];
    let dt = sim_params.dt;
    let total_time = sim_params.total_time;
    let spawn_interval = sim_params.spawn_interval;
    let max_age = sim_params.max_age;
    let num_waves = max(1u, sim_params.num_waves);
    
    // Calculate which wave this particle belongs to (fixed assignment)
    let ppb = max(1u, sim_params.particles_per_blob);
    let particles_per_wave = ppb * sim_params.spawn_count;
    let wave_id = (idx / particles_per_wave) % num_waves;
    
    // Calculate particle age from global time (not incremental)
    // Waves spawn on a fixed schedule: wave 0 at t=0, wave 1 at t=spawn_interval, etc.
    // After num_waves spawns, the cycle repeats
    let cycle_duration = f32(num_waves) * spawn_interval;
    
    // Time within current cycle
    let time_in_cycle = total_time % cycle_duration;
    
    // When did this wave spawn in the current cycle?
    let wave_spawn_time = f32(wave_id) * spawn_interval;
    
    // Calculate age: time since this wave spawned
    var age: f32;
    if (time_in_cycle >= wave_spawn_time) {
        // Wave has spawned in current cycle
        age = time_in_cycle - wave_spawn_time;
    } else {
        // Wave hasn't spawned yet in current cycle, use age from end of previous cycle
        age = (cycle_duration - wave_spawn_time) + time_in_cycle;
    }
    
    // Check if particle is "dead" (exceeded max_age) or inside airfoil
    let is_dead = age >= max_age || p.psi > 1e9;
    
    // Check if particle just spawned (age wrapped around)
    // A particle "just spawned" if its new age is less than its stored age (plus some tolerance)
    let just_spawned = age < p.age - dt * 2.0;
    
    // Check if particle hit airfoil (psi > 1e9 indicates interior)
    let hit_airfoil = p.psi > 1e9;
    
    if (just_spawned || (is_dead && age < spawn_interval)) {
        // Respawn at spawn point
        let blob_id = (idx / ppb) % sim_params.spawn_count;
        let spawn_pos = spawn_points[blob_id];
        
        // Use particle index for consistent blob randomization
        let rand1 = hash2(idx, 12345u);
        let rand2 = hash2(idx * 7u, 12345u);
        
        // Random position within blob (uniform disk distribution)
        let angle = rand1 * 6.28318;
        let r = sqrt(rand2) * 0.02;  // Blob radius = 2% of chord
        
        p.pos = spawn_pos + vec2f(cos(angle), sin(angle)) * r;
        p.psi = sample_psi(p.pos);
    } else if (hit_airfoil) {
        // Particle hit airfoil - move off-screen to prevent smearing
        // It will respawn when its wave comes around again
        p.pos = vec2f(-1000.0, -1000.0);
        p.psi = 0.0;  // Reset psi so it's not stuck in "dead" state
    } else if (!is_dead) {
        // RK2 midpoint integration with velocity from texture
        let v1 = sample_velocity(p.pos);
        
        // Use freestream as fallback if velocity is zero (texture not initialized)
        let freestream = vec2f(sim_params.cos_alpha, sim_params.sin_alpha);
        let speed1 = length(v1);
        let vel1 = select(freestream, v1, speed1 > 0.01);
        
        let mid_pos = p.pos + 0.5 * dt * vel1;
        let v2 = sample_velocity(mid_pos);
        let speed2 = length(v2);
        let vel2 = select(freestream, v2, speed2 > 0.01);
        
        p.pos += dt * vel2;

        // Update stream function value
        p.psi = sample_psi(p.pos);
    }
    
    // Store calculated age (for spawn detection and rendering)
    p.age = age;

    particles[idx] = p;
}

// ============================================================================
// Render Shader: Instanced Point Sprites
// ============================================================================

struct VertexOutput {
    @builtin(position) position: vec4f,
    @location(0) alpha: f32,
    @location(1) is_above_dividing: f32,
    @location(2) local_pos: vec2f,
}

@group(0) @binding(0) var<uniform> render_params: RenderParams;
@group(0) @binding(1) var<storage, read> render_particles: array<Particle>;

@vertex
fn vs_particle(
    @builtin(vertex_index) vertex_idx: u32,
    @builtin(instance_index) instance_idx: u32,
) -> VertexOutput {
    let p = render_particles[instance_idx];
    
    // Quad vertices for point sprite
    var quad_pos: array<vec2f, 6> = array<vec2f, 6>(
        vec2f(-1.0, -1.0),
        vec2f( 1.0, -1.0),
        vec2f(-1.0,  1.0),
        vec2f(-1.0,  1.0),
        vec2f( 1.0, -1.0),
        vec2f( 1.0,  1.0),
    );
    
    let local = quad_pos[vertex_idx];
    
    // Rotate particle position from body frame to display frame
    // This matches CPU rotatePoint: x' = cx + dx*cos - dy*sin, y' = cy + dx*sin + dy*cos
    // where cos/sin are computed from -alpha (so cos(-α)=cos(α), sin(-α)=-sin(α))
    let dx = p.pos.x - render_params.rot_cx;
    let dy = p.pos.y - render_params.rot_cy;
    let rotated_x = render_params.rot_cx + dx * render_params.cos_alpha - dy * render_params.sin_alpha;
    let rotated_y = render_params.rot_cy + dx * render_params.sin_alpha + dy * render_params.cos_alpha;
    
    // Transform to clip space
    let world_pos = vec2f(rotated_x, rotated_y);
    let view_offset = (world_pos - vec2f(render_params.view_center_x, render_params.view_center_y)) * render_params.zoom;
    
    // Add quad offset - scale by zoom so particle size is consistent in world units
    // point_size is in world units (e.g., 0.01 = 1% of chord)
    let point_size = render_params.point_size;
    let quad_offset = local * point_size * render_params.zoom / vec2f(render_params.aspect, 1.0);
    
    var out: VertexOutput;
    out.position = vec4f(view_offset.x / render_params.aspect + quad_offset.x, view_offset.y + quad_offset.y, 0.0, 1.0);
    
    // Compute alpha based on age
    let max_age = render_params._pad3.x;
    let t = p.age / max_age;
    
    // Dead particles (age > max_age) should be invisible
    if (t >= 1.0) {
        out.alpha = 0.0;
    } else if (t < 0.05) {
        out.alpha = t / 0.05;
    } else if (t > 0.85) {
        out.alpha = (1.0 - t) / 0.15;
    } else {
        out.alpha = 1.0;
    }
    out.alpha *= 0.6;  // Base opacity
    
    // Determine if above or below dividing streamline
    // Only meaningful for live particles (dead ones are invisible anyway)
    out.is_above_dividing = select(0.0, 1.0, p.psi > render_params.psi_0);
    
    out.local_pos = local;
    
    return out;
}

@fragment
fn fs_particle(in: VertexOutput) -> @location(0) vec4f {
    // Skip dead particles (alpha = 0)
    if (in.alpha <= 0.0) {
        discard;
    }
    
    // Soft circle
    let dist = length(in.local_pos);
    if (dist > 1.0) {
        discard;
    }
    
    let soft_alpha = 1.0 - smoothstep(0.5, 1.0, dist);
    
    // Color based on dividing streamline
    var color: vec3f;
    if (in.is_above_dividing > 0.5) {
        color = vec3f(render_params.color_above_r, render_params.color_above_g, render_params.color_above_b);
    } else {
        color = vec3f(render_params.color_below_r, render_params.color_below_g, render_params.color_below_b);
    }
    
    return vec4f(color, in.alpha * soft_alpha);
}
