/**
 * Simple Streamline Renderer
 * 
 * Renders pre-computed streamlines as line segments.
 * Streamline data is passed via vertex buffer.
 */

struct Params {
    view_center_x: f32,
    view_center_y: f32,
    zoom: f32,
    aspect: f32,
    line_color_r: f32,
    line_color_g: f32,
    line_color_b: f32,
    line_alpha: f32,
}

@group(0) @binding(0) var<uniform> params: Params;

struct VertexInput {
    @location(0) position: vec2f,
}

struct VertexOutput {
    @builtin(position) position: vec4f,
}

@vertex
fn vs_streamline(in: VertexInput) -> VertexOutput {
    // Transform world position to clip space
    let world = in.position;
    
    // Apply view transform
    let view_offset = world - vec2f(params.view_center_x, params.view_center_y);
    let scaled = view_offset * params.zoom;
    let ndc = vec2f(scaled.x / params.aspect, scaled.y);
    
    var out: VertexOutput;
    out.position = vec4f(ndc, 0.0, 1.0);
    return out;
}

@fragment
fn fs_streamline() -> @location(0) vec4f {
    return vec4f(params.line_color_r, params.line_color_g, params.line_color_b, params.line_alpha);
}
