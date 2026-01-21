/**
 * Streamline Visualization Shader
 * 
 * Instead of LIC (which produces uniform bands for uniform flow),
 * this shader renders discrete streamlines starting from a grid of seed points.
 * Each streamline is rendered as a thin line that can be followed visually.
 */

struct Params {
    // Velocity texture bounds (world coords -> UV mapping)
    tex_x_min: f32,
    tex_x_max: f32,
    tex_y_min: f32,
    tex_y_max: f32,
    // Streamline settings
    step_size: f32,      // Step size in WORLD units
    num_steps: u32,      // Steps in each direction
    line_spacing: f32,   // Spacing between seed points
    line_width: f32,     // Visual width of lines
    // View transform
    cos_alpha: f32,
    sin_alpha: f32,
    rot_cx: f32,
    rot_cy: f32,
    // Viewport
    view_center_x: f32,
    view_center_y: f32,
    zoom: f32,
    aspect: f32,
    // Theme
    is_dark: f32,
    _pad: vec3f,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var velocity_texture: texture_2d<f32>;
@group(0) @binding(2) var velocity_sampler: sampler;
@group(0) @binding(3) var noise_texture: texture_2d<f32>;
@group(0) @binding(4) var noise_sampler: sampler;

struct VertexOutput {
    @builtin(position) position: vec4f,
    @location(0) world_pos: vec2f,
}

@vertex
fn vs_fullscreen(@builtin(vertex_index) idx: u32) -> VertexOutput {
    var positions = array<vec2f, 3>(
        vec2f(-1.0, -1.0),
        vec2f( 3.0, -1.0),
        vec2f(-1.0,  3.0),
    );
    
    let clip_pos = positions[idx];
    let ndc = clip_pos;
    let view_offset = ndc * vec2f(params.aspect, 1.0) / params.zoom;
    var world = view_offset + vec2f(params.view_center_x, params.view_center_y);
    
    // Inverse rotation
    let dx = world.x - params.rot_cx;
    let dy = world.y - params.rot_cy;
    world.x = params.rot_cx + dx * params.cos_alpha - dy * params.sin_alpha;
    world.y = params.rot_cy + dx * params.sin_alpha + dy * params.cos_alpha;
    
    var out: VertexOutput;
    out.position = vec4f(clip_pos, 0.0, 1.0);
    out.world_pos = world;
    return out;
}

// Map world coords to velocity texture UV [0,1]
fn world_to_uv(world: vec2f) -> vec2f {
    return vec2f(
        (world.x - params.tex_x_min) / (params.tex_x_max - params.tex_x_min),
        (world.y - params.tex_y_min) / (params.tex_y_max - params.tex_y_min)
    );
}

// Sample velocity at world position, return (u, v, speed, inside_flag)
fn get_velocity(world: vec2f) -> vec4f {
    let uv = world_to_uv(world);
    
    // Check if in bounds
    if (uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0) {
        return vec4f(1.0, 0.0, 1.0, 1.0); // Default horizontal flow
    }
    
    let vel = textureSampleLevel(velocity_texture, velocity_sampler, uv, 0.0);
    let u = vel.r;
    let v = vel.g;
    let speed = length(vec2f(u, v));
    
    // Check if inside airfoil (psi > 1e9)
    if (vel.a > 1e9) {
        return vec4f(0.0, 0.0, 0.0, -1.0); // Inside airfoil
    }
    
    if (speed < 0.0001) {
        return vec4f(1.0, 0.0, 0.01, 1.0); // Stagnation, default horizontal
    }
    
    return vec4f(u, v, speed, 1.0);
}

// Check if point p is on a streamline starting from seed point
// Returns distance to nearest streamline segment
fn distance_to_streamline(p: vec2f, seed: vec2f) -> f32 {
    var min_dist = 1000.0;
    let dt = params.step_size;
    let max_steps = params.num_steps;
    
    // Start velocity at seed
    let start_vel = get_velocity(seed);
    if (start_vel.w < 0.0) { return min_dist; } // Seed inside airfoil
    
    // Forward trace
    var pos = seed;
    var prev_pos = seed;
    var vel = vec2f(start_vel.x, start_vel.y);
    let speed = start_vel.z;
    if (speed > 0.001) {
        vel = vel / speed;
    }
    
    for (var i: u32 = 0u; i < max_steps; i++) {
        prev_pos = pos;
        pos = pos + vel * dt;
        
        // Get velocity at new position
        let v = get_velocity(pos);
        if (v.w < 0.0) { break; } // Hit airfoil
        
        // Update velocity direction
        let spd = v.z;
        if (spd > 0.001) {
            vel = vec2f(v.x, v.y) / spd;
        }
        
        // Distance from p to line segment prev_pos -> pos
        let segment = pos - prev_pos;
        let seg_len = length(segment);
        if (seg_len > 0.0001) {
            let t = clamp(dot(p - prev_pos, segment) / (seg_len * seg_len), 0.0, 1.0);
            let proj = prev_pos + t * segment;
            let d = length(p - proj);
            min_dist = min(min_dist, d);
        }
    }
    
    // Backward trace
    pos = seed;
    vel = vec2f(start_vel.x, start_vel.y);
    if (speed > 0.001) {
        vel = vel / speed;
    }
    
    for (var i: u32 = 0u; i < max_steps; i++) {
        prev_pos = pos;
        pos = pos - vel * dt; // Backward
        
        let v = get_velocity(pos);
        if (v.w < 0.0) { break; }
        
        let spd = v.z;
        if (spd > 0.001) {
            vel = vec2f(v.x, v.y) / spd;
        }
        
        let segment = pos - prev_pos;
        let seg_len = length(segment);
        if (seg_len > 0.0001) {
            let t = clamp(dot(p - prev_pos, segment) / (seg_len * seg_len), 0.0, 1.0);
            let proj = prev_pos + t * segment;
            let d = length(p - proj);
            min_dist = min(min_dist, d);
        }
    }
    
    return min_dist;
}

@fragment
fn fs_lic(in: VertexOutput) -> @location(0) vec4f {
    let p = in.world_pos;
    
    // Check if inside airfoil
    let vel = get_velocity(p);
    if (vel.w < 0.0) {
        return vec4f(0.0, 0.0, 0.0, 0.0);
    }
    
    // Grid of seed points
    let spacing = params.line_spacing;
    let half_spacing = spacing * 0.5;
    
    // Find nearest seed points (check 3x3 grid around current position)
    var min_dist = 1000.0;
    let base_seed_x = floor(p.x / spacing) * spacing;
    let base_seed_y = floor(p.y / spacing) * spacing;
    
    for (var dy: i32 = -1; dy <= 1; dy++) {
        for (var dx: i32 = -1; dx <= 1; dx++) {
            let seed = vec2f(
                base_seed_x + f32(dx) * spacing + half_spacing,
                base_seed_y + f32(dy) * spacing + half_spacing
            );
            let d = distance_to_streamline(p, seed);
            min_dist = min(min_dist, d);
        }
    }
    
    // Render as anti-aliased line
    let line_width = params.line_width;
    let edge = line_width * 0.5;
    
    if (min_dist < edge) {
        // Inside line - smooth falloff at edges
        let alpha = 1.0 - smoothstep(edge * 0.5, edge, min_dist);
        
        var color: vec3f;
        if (params.is_dark > 0.5) {
            color = vec3f(0.6, 0.7, 0.8);
        } else {
            color = vec3f(0.2, 0.3, 0.4);
        }
        
        return vec4f(color, alpha * 0.8);
    }
    
    return vec4f(0.0, 0.0, 0.0, 0.0);
}
