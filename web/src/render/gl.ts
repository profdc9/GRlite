/* Minimal WebGL2 shader helpers. */

export function compileShader(gl: WebGL2RenderingContext, type: number, src: string): WebGLShader {
    const sh = gl.createShader(type);
    if (!sh) throw new Error('createShader failed');
    gl.shaderSource(sh, src);
    gl.compileShader(sh);
    if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
        const log = gl.getShaderInfoLog(sh);
        gl.deleteShader(sh);
        throw new Error(`shader compile error: ${log}`);
    }
    return sh;
}

export function linkProgram(gl: WebGL2RenderingContext, vs: string, fs: string): WebGLProgram {
    const vsh = compileShader(gl, gl.VERTEX_SHADER, vs);
    const fsh = compileShader(gl, gl.FRAGMENT_SHADER, fs);
    const prog = gl.createProgram();
    if (!prog) throw new Error('createProgram failed');
    gl.attachShader(prog, vsh);
    gl.attachShader(prog, fsh);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
        const log = gl.getProgramInfoLog(prog);
        gl.deleteProgram(prog);
        throw new Error(`program link error: ${log}`);
    }
    return prog;
}
