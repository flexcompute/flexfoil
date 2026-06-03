/**
 * Two-pass Line Integral Convolution shaders, ported verbatim from the flow360 web
 * app (flex/.../Visualization/LIC-shader.ts).
 *
 * Pass 1 (LIC_Dir): render the flow mesh with a per-vertex flow vector `_cfvec`;
 *   output a screen-space "flow texture" = (dirX, dirY in [0,1], procedural noise).
 * Pass 2 (LIC): at each fragment, step forward & backward along the flow direction
 *   (read from the flow texture) `passNumber` times, accumulating the noise → the LIC
 *   streak intensity `k`, shaded by the per-vertex color `vColor`.
 *
 * Used by lic.ts with Three.js ShaderMaterials (GLSL ES / WebGL).
 */

export const LIC_Dir_VS = /* glsl */ `
    attribute vec4 _cfvec;
    varying vec2 vCfVec;
    varying vec3 vPos;
    void main() {
        gl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );
        // vCfVec is the flow vector transformed into screen space, mapped [-1,1]->[0,1].
        vCfVec = normalize((modelViewMatrix * vec4( _cfvec.xyz, 0.)).xy) * 0.5 + 0.5;
        vPos = (modelMatrix * vec4(position, 0)).xyz;
    }
`;

export const LIC_Dir_FS = /* glsl */ `
    varying vec2 vCfVec;
    varying vec3 vPos;
    uniform sampler2D noiseTexture;
    uniform float textureRepeat;
    float hash(vec3 p) {
        p  = 50.0*fract( p*0.3183099 + vec3(0.71,0.113,0.419));
        return -1.0+2.0*fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
    }
    float noised(vec3 x ){
        vec3 i = floor(x);
        vec3 w = fract(x);
        vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
        float a = hash(i+vec3(0.0,0.0,0.0));
        float b = hash(i+vec3(1.0,0.0,0.0));
        float c = hash(i+vec3(0.0,1.0,0.0));
        float d = hash(i+vec3(1.0,1.0,0.0));
        float e = hash(i+vec3(0.0,0.0,1.0));
        float f = hash(i+vec3(1.0,0.0,1.0));
        float g = hash(i+vec3(0.0,1.0,1.0));
        float h = hash(i+vec3(1.0,1.0,1.0));
        float k0 =   a;
        float k1 =   b - a;
        float k2 =   c - a;
        float k3 =   e - a;
        float k4 =   a - b - c + d;
        float k5 =   a - c - e + g;
        float k6 =   a - b - e + f;
        float k7 = - a + b + c - d + e - f - g + h;
        return k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z;
    }
    void main() {
        float col = 0.5 + 0.5 * noised(vPos * textureRepeat);
        gl_FragColor = vec4(vCfVec, col, 1.0);
    }
`;

export const LIC_VS = /* glsl */ `
    varying vec3 vCoord;
    varying vec3 vPos;
    varying vec3 vColor;
    varying vec3 vNormal;
    void main() {
        gl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );
        vCoord = gl_Position.xyz / gl_Position.w  * 0.5 + 0.5;
        vPos = position;
        vColor = color;
        vNormal = normalize((modelViewMatrix * vec4(normal, 0.0)).xyz);
    }
`;

export const LIC_FS = /* glsl */ `
    varying vec3 vCoord;
    varying vec3 vPos;
    varying vec3 vColor;
    varying vec3 vNormal;
    uniform sampler2D noiseTexture;
    uniform sampler2D flowTexture;
    uniform float sizeK;
    uniform float textureRepeat;
    uniform float passStep;
    const int passNumber = 40;   // fixed (GLSL ES 1.0 needs a constant loop bound)
    void main() {
        float onestep = passStep;
        vec2 coord1 = vCoord.xy;
        vec2 coord2 = vCoord.xy;
        vec3 currentColor0 = texture2D(flowTexture, coord1 ).xyz;
        vec2 flowStep0 = (currentColor0.xy - 0.5) * 2.0 * onestep;
        flowStep0.y *= sizeK;
        coord1 += flowStep0;
        coord2 -= flowStep0;
        float col = texture2D(flowTexture, coord1 ).z;
        for(int i = 1; i < passNumber; i++) {
            vec3 currentColor1 = texture2D(flowTexture, coord1 ).xyz;
            vec3 currentColor2 = texture2D(flowTexture, coord2 ).xyz;
            vec2 flowStep1 = (currentColor1.xy - 0.5) * 2.0 * onestep;
            vec2 flowStep2 = (currentColor2.xy - 0.5) * 2.0 * onestep;
            flowStep1.y *= sizeK;
            flowStep2.y *= sizeK;
            coord1 += flowStep1;
            coord2 -= flowStep2;
            col += currentColor1.z;
            col += currentColor2.z;
        }
        float k = smoothstep(0.2, 0.7, col / float(passNumber * 2 - 1) );
        float light = 0.25 + 0.75 * dot(vNormal, vec3(0., 0., 1.));
        gl_FragColor = vec4(vColor * k * light, 1.0);
    }
`;
