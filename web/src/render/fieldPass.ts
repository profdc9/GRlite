/* Field render pass: uploads a scalar field (R32F) to a texture and draws
 * it full-screen through a diverging colormap.  Auto-normalizes to the
 * frame's peak |value| so tiny and large fields both show structure. */

import { linkProgram } from './gl';
import { DISPLAY_SCALE } from '../sim/config';

const FIELD_VS = `#version 300 es
in vec2 a_pos;
out vec2 v_uv;
void main() {
    v_uv = a_pos * 0.5 + 0.5;
    gl_Position = vec4(a_pos, 0.0, 1.0);
}`;

const FIELD_FS = `#version 300 es
precision highp float;
in vec2 v_uv;
out vec4 outColor;
uniform sampler2D u_field;
uniform float u_scale;
void main() {
    float v = texture(u_field, v_uv).r * u_scale;
    float pos = clamp(v, 0.0, 1.0);
    float neg = clamp(-v, 0.0, 1.0);
    /* cyan (negative) -> black (zero) -> orange (positive) */
    vec3 c = vec3(pos, pos * 0.5 + neg * 0.6, neg);
    outColor = vec4(c, 1.0);
}`;

export class FieldPass {
    private gl: WebGL2RenderingContext;
    private prog: WebGLProgram;
    private vao: WebGLVertexArrayObject;
    private tex: WebGLTexture;
    private uScale: WebGLUniformLocation | null;
    readonly W: number;
    readonly H: number;

    constructor(gl: WebGL2RenderingContext, w: number, h: number) {
        this.gl = gl; this.W = w; this.H = h;
        this.prog = linkProgram(gl, FIELD_VS, FIELD_FS);

        this.vao = gl.createVertexArray()!;
        gl.bindVertexArray(this.vao);
        const vbo = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
        gl.bufferData(gl.ARRAY_BUFFER,
            new Float32Array([-1,-1, 1,-1, -1,1, -1,1, 1,-1, 1,1]), gl.STATIC_DRAW);
        const aPos = gl.getAttribLocation(this.prog, 'a_pos');
        gl.enableVertexAttribArray(aPos);
        gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

        this.tex = gl.createTexture()!;
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this.tex);
        /* R32F isn't filterable without OES_texture_float_linear -> NEAREST. */
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R32F, w, h);

        gl.useProgram(this.prog);
        gl.uniform1i(gl.getUniformLocation(this.prog, 'u_field'), 0);
        this.uScale = gl.getUniformLocation(this.prog, 'u_scale');
        gl.uniform1f(this.uScale, DISPLAY_SCALE);
    }

    /* Upload + draw.  `scale` overrides auto-normalization when > 0. */
    render(data: Float32Array, scale = 0): void {
        const gl = this.gl;
        gl.bindTexture(gl.TEXTURE_2D, this.tex);
        gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, this.W, this.H, gl.RED, gl.FLOAT, data);

        let s = scale;
        if (s <= 0) {
            let maxAbs = 0;
            for (let k = 0; k < data.length; k++) {
                const a = data[k] < 0 ? -data[k] : data[k];
                if (a > maxAbs) maxAbs = a;
            }
            s = maxAbs > 1e-12 ? 1.0 / maxAbs : DISPLAY_SCALE;
        }
        gl.useProgram(this.prog);
        gl.uniform1f(this.uScale, s);
        gl.bindVertexArray(this.vao);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this.tex);
        gl.drawArrays(gl.TRIANGLES, 0, 6);
    }
}
