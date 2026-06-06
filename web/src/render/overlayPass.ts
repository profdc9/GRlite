/* Overlay render pass: particles as point sprites over the field.
 *
 * Phase 1c will extend this with trails, velocity arrows, and force
 * arrows; for now it draws round particle discs at their grid positions. */

import { linkProgram } from './gl';

const POINT_VS = `#version 300 es
in vec2 a_pos;                 /* sim cells, [0,W) x [0,H) */
uniform vec2 u_gridSize;
uniform float u_pointSize;
void main() {
    vec2 ndc = (a_pos / u_gridSize) * 2.0 - 1.0;
    gl_Position = vec4(ndc, 0.0, 1.0);
    gl_PointSize = u_pointSize;
}`;

const POINT_FS = `#version 300 es
precision highp float;
out vec4 outColor;
uniform vec3 u_color;
void main() {
    vec2 d = gl_PointCoord - vec2(0.5);
    float r2 = dot(d, d);
    if (r2 > 0.25) discard;
    float edge = smoothstep(0.25, 0.16, r2);
    outColor = vec4(u_color, edge);
}`;

export class OverlayPass {
    private gl: WebGL2RenderingContext;
    private prog: WebGLProgram;
    private vao: WebGLVertexArrayObject;
    private vbo: WebGLBuffer;
    private uColor: WebGLUniformLocation | null;

    constructor(gl: WebGL2RenderingContext, w: number, h: number, pointSize: number) {
        this.gl = gl;
        this.prog = linkProgram(gl, POINT_VS, POINT_FS);
        this.vao = gl.createVertexArray()!;
        gl.bindVertexArray(this.vao);
        this.vbo = gl.createBuffer()!;
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo);
        const aPos = gl.getAttribLocation(this.prog, 'a_pos');
        gl.enableVertexAttribArray(aPos);
        gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);
        gl.useProgram(this.prog);
        gl.uniform2f(gl.getUniformLocation(this.prog, 'u_gridSize'), w, h);
        gl.uniform1f(gl.getUniformLocation(this.prog, 'u_pointSize'), pointSize);
        this.uColor = gl.getUniformLocation(this.prog, 'u_color');
        gl.uniform3f(this.uColor, 1.0, 1.0, 1.0);
    }

    /* xy: interleaved [x0,y0, x1,y1, ...] in grid cells; n particles. */
    render(xy: Float32Array, n: number): void {
        if (n <= 0) return;
        const gl = this.gl;
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        gl.useProgram(this.prog);
        gl.bindVertexArray(this.vao);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo);
        gl.bufferData(gl.ARRAY_BUFFER, xy.subarray(0, n * 2), gl.DYNAMIC_DRAW);
        gl.drawArrays(gl.POINTS, 0, n);
        gl.disable(gl.BLEND);
    }
}
