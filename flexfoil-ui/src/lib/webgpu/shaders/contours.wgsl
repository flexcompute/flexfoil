/**
 * Stream Function (ψ) Contour Shader
 * 
 * Full-screen quad that samples velocity texture and renders:
 * - Filled color bands (blue for ψ < ψ₀, red for ψ > ψ₀)
 * - Anti-aliased contour lines using fwidth()
 * - Thick dividing streamline (ψ = ψ₀)
 */

struct Params {
    // Transform: world coords -> UV
    tex_x_min: f32,
    tex_x_max: f32,
    tex_y_min: f32,
    tex_y_max: f32,
    // Stream function range
    psi_0: f32,
    psi_min: f32,
    psi_max: f32,
    // Contour settings
    contour_interval: f32,
    // View transform for rotation
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

struct VertexOutput {
    @builtin(position) position: vec4f,
    @location(0) world_pos: vec2f,
}

@vertex
fn vs_fullscreen(@builtin(vertex_index) idx: u32) -> VertexOutput {
    // Full-screen triangle (3 vertices that cover the screen)
    var positions = array<vec2f, 3>(
        vec2f(-1.0, -1.0),
        vec2f( 3.0, -1.0),
        vec2f(-1.0,  3.0),
    );
    
    let clip_pos = positions[idx];
    
    // Convert clip space to world coordinates (accounting for view)
    // Inverse of the standard view transform
    let ndc = clip_pos;
    let view_offset = ndc * vec2f(params.aspect, 1.0) / params.zoom;
    var world = view_offset + vec2f(params.view_center_x, params.view_center_y);
    
    // Inverse rotation to get to airfoil frame
    let dx = world.x - params.rot_cx;
    let dy = world.y - params.rot_cy;
    world.x = params.rot_cx + dx * params.cos_alpha - dy * params.sin_alpha;
    world.y = params.rot_cy + dx * params.sin_alpha + dy * params.cos_alpha;
    
    var out: VertexOutput;
    out.position = vec4f(clip_pos, 0.0, 1.0);
    out.world_pos = world;
    return out;
}

// Map world coords to texture UV
fn world_to_uv(pos: vec2f) -> vec2f {
    return vec2f(
        (pos.x - params.tex_x_min) / (params.tex_x_max - params.tex_x_min),
        (pos.y - params.tex_y_min) / (params.tex_y_max - params.tex_y_min)
    );
}

// Blue gradient for ψ < ψ₀ (flow going under)
fn blue_gradient(t: f32) -> vec3f {
    // t is 0..1 where 1 is furthest from ψ₀
    let intensity = pow(t, 0.5);
    if (params.is_dark > 0.5) {
        // Dark theme: lighter blues
        return vec3f(
            0.86 - intensity * 0.71,  // 220 -> 40
            0.90 - intensity * 0.63,  // 230 -> 70
            1.0 - intensity * 0.22    // 255 -> 200
        );
    } else {
        // Light theme: same
        return vec3f(
            0.86 - intensity * 0.71,
            0.90 - intensity * 0.63,
            1.0 - intensity * 0.22
        );
    }
}

// Red gradient for ψ > ψ₀ (flow going over)
fn red_gradient(t: f32) -> vec3f {
    let intensity = pow(t, 0.5);
    if (params.is_dark > 0.5) {
        return vec3f(
            1.0 - intensity * 0.22,   // 255 -> 200
            0.86 - intensity * 0.63,  // 220 -> 60
            0.82 - intensity * 0.59   // 210 -> 60
        );
    } else {
        return vec3f(
            1.0 - intensity * 0.22,
            0.86 - intensity * 0.63,
            0.82 - intensity * 0.59
        );
    }
}

@fragment
fn fs_contour(in: VertexOutput) -> @location(0) vec4f {
    let uv = world_to_uv(in.world_pos);
    
    // Sample stream function from velocity texture (alpha channel)
    // All texture samples and derivatives MUST happen before any non-uniform control flow
    let clamped_uv = clamp(uv, vec2f(0.0), vec2f(1.0));
    let sample = textureSample(velocity_texture, velocity_sampler, clamped_uv);
    let psi = sample.a;
    
    // Compute derivatives in uniform control flow (before any if/discard)
    let psi_derivative = fwidth(psi);
    
    // Anti-aliased contour lines using fwidth()
    let interval = params.contour_interval;
    let psi_normalized = psi / interval;
    let dist_to_line = abs(fract(psi_normalized + 0.5) - 0.5) * interval;
    let line_width = psi_derivative * 2.0;
    let line = 1.0 - smoothstep(0.0, line_width, dist_to_line);
    
    // Dividing streamline (ψ = ψ₀) - thicker and more prominent
    let psi0_dist = abs(psi - params.psi_0);
    let dividing_width = psi_derivative * 4.0;
    let dividing = 1.0 - smoothstep(0.0, dividing_width, psi0_dist);
    
    // Compute colors in uniform control flow
    let psi_range_below = params.psi_0 - params.psi_min;
    let psi_range_above = params.psi_max - params.psi_0;
    
    // Compute both color options
    let t_below = clamp((params.psi_0 - psi) / max(psi_range_below, 0.001), 0.0, 1.0);
    let t_above = clamp((psi - params.psi_0) / max(psi_range_above, 0.001), 0.0, 1.0);
    let color_below = blue_gradient(t_below);
    let color_above = red_gradient(t_above);
    
    // Select based on psi (using select for uniform flow)
    let color = select(color_above, color_below, psi < params.psi_0);
    
    // Line colors based on theme (use select, not if)
    let line_color = select(0.2, 0.3, params.is_dark > 0.5);
    let dividing_color = select(0.1, 0.95, params.is_dark > 0.5);
    
    // Transparent color for areas where we don't want to draw
    // This lets the main canvas (with airfoil) show through
    let transparent = vec4f(0.0, 0.0, 0.0, 0.0);
    
    // Combine: base color + contour lines + dividing streamline
    let alpha = 0.75;
    var final_color = color;
    final_color = mix(final_color, vec3f(line_color), line * 0.3);
    final_color = mix(final_color, vec3f(dividing_color), dividing * 0.9);
    let result = vec4f(final_color, alpha);
    
    // Check bounds and airfoil interior
    let out_of_bounds = uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0;
    let inside_airfoil = psi > 1e9;
    
    // For out of bounds or inside airfoil, return transparent
    // This lets the main canvas show through
    if (out_of_bounds || inside_airfoil) {
        return transparent;
    }
    
    return result;
}
