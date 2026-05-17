/* GRlite Stage 1 frontend.
 *
 * Loads the WASM core, instantiates a 2D scalar-wave simulation seeded with the
 * `wave_pulse` scenario (gr_sandbox_v32.tex §12.1 Stage 1), and renders the field
 * as a WebGL2 R32F texture each animation frame. No physics lives here — TS only
 * marshals parameters in and pulls pointers out per memory `grlite-stack`.
 */

import type { GRliteModule } from './grlite';

const GRID_W = 256;
const GRID_H = 256;
const DX = 1.0;
const C_EFF = 1.0;
const CFL = 1.0 / Math.sqrt(2);          // 2D stability limit, §9.2 eq:cfl
const PULSE_SIGMA = 4.0 * DX;
const PULSE_AMP = 1.0;
const STEPS_PER_FRAME = 4;
const DISPLAY_SCALE = 8.0;                // visual gain on |phi| so the wavefront is visible

interface SimAPI {
    sim: number;
    create: (w: number, h: number, dx: number, c: number, cfl: number) => number;
    destroy: (sim: number) => void;
    step: (sim: number) => void;
    stepN: (sim: number, n: number) => void;
    fieldPtr: (sim: number, which: number) => number;
    loadScenario: (sim: number, name: string, paramsPtr: number, n: number) => number;
    stepCount: (sim: number) => number;
    simTime: (sim: number) => number;
}

async function loadWasm(): Promise<GRliteModule> {
    /* Vite serves /grlite/grlite.js from web/public/grlite/, produced by
     * `mingw32-make wasm` in core/. The path is a runtime URL, not a TS or
     * Vite module — both static resolvers need to be told to skip it. */
    // @ts-expect-error  WASM bundle resolves at runtime via Vite's public/ dir
    const mod = (await import(/* @vite-ignore */ '/grlite/grlite.js')) as {
        default: () => Promise<GRliteModule>;
    };
    return mod.default();
}

function bindApi(M: GRliteModule): SimAPI {
    const create = M.cwrap('gr_sim_create', 'number',
        ['number', 'number', 'number', 'number', 'number']) as
        (w: number, h: number, dx: number, c: number, cfl: number) => number;
    const destroy = M.cwrap('gr_sim_destroy', null, ['number']) as (sim: number) => void;
    const step = M.cwrap('gr_sim_step', null, ['number']) as (sim: number) => void;
    const stepN = M.cwrap('gr_sim_step_n', null, ['number', 'number']) as
        (sim: number, n: number) => void;
    const fieldPtr = M.cwrap('gr_sim_field_ptr', 'number', ['number', 'number']) as
        (sim: number, which: number) => number;
    const loadScenario = M.cwrap('gr_sim_load_scenario', 'number',
        ['number', 'string', 'number', 'number']) as
        (sim: number, name: string, paramsPtr: number, n: number) => number;
    const stepCount = M.cwrap('gr_sim_step_count', 'number', ['number']) as
        (sim: number) => number;
    const simTime = M.cwrap('gr_sim_time', 'number', ['number']) as (sim: number) => number;

    const sim = create(GRID_W, GRID_H, DX, C_EFF, CFL);
    if (!sim) throw new Error('gr_sim_create returned NULL');

    // Marshal Gaussian-pulse params via the WASM heap.
    const params = new Float32Array([PULSE_SIGMA, PULSE_AMP]);
    const paramsPtr = M._malloc(params.byteLength);
    M.HEAPF32.set(params, paramsPtr >> 2);
    const rc = loadScenario(sim, 'wave_pulse', paramsPtr, params.length);
    M._free(paramsPtr);
    if (rc !== 0) throw new Error(`load_scenario failed: rc=${rc}`);

    return { sim, create, destroy, step, stepN, fieldPtr, loadScenario, stepCount, simTime };
}

/* WebGL2 setup: one fullscreen quad, one R32F texture, diverging colormap fragment shader. */

function compileShader(gl: WebGL2RenderingContext, type: number, src: string): WebGLShader {
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

function linkProgram(gl: WebGL2RenderingContext, vsSrc: string, fsSrc: string): WebGLProgram {
    const vs = compileShader(gl, gl.VERTEX_SHADER, vsSrc);
    const fs = compileShader(gl, gl.FRAGMENT_SHADER, fsSrc);
    const prog = gl.createProgram();
    if (!prog) throw new Error('createProgram failed');
    gl.attachShader(prog, vs);
    gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
        const log = gl.getProgramInfoLog(prog);
        gl.deleteProgram(prog);
        throw new Error(`program link error: ${log}`);
    }
    return prog;
}

const VS_SRC = `#version 300 es
in vec2 a_pos;
out vec2 v_uv;
void main() {
    v_uv = a_pos * 0.5 + 0.5;
    gl_Position = vec4(a_pos, 0.0, 1.0);
}`;

/* Diverging colormap: blue (negative) -> black (zero) -> red (positive). */
const FS_SRC = `#version 300 es
precision highp float;
in vec2 v_uv;
out vec4 outColor;
uniform sampler2D u_field;
uniform float u_scale;
void main() {
    float v = texture(u_field, v_uv).r * u_scale;
    float pos = clamp(v, 0.0, 1.0);
    float neg = clamp(-v, 0.0, 1.0);
    vec3 c = vec3(pos, 0.0, neg);
    outColor = vec4(c, 1.0);
}`;

async function main(): Promise<void> {
    const statusEl = document.getElementById('status') as HTMLDivElement;
    const canvas = document.getElementById('view') as HTMLCanvasElement;
    const resetBtn = document.getElementById('reset') as HTMLButtonElement;
    const pauseBtn = document.getElementById('pause') as HTMLButtonElement;

    statusEl.textContent = 'loading WASM…';
    const M = await loadWasm();
    const api = bindApi(M);

    const gl = canvas.getContext('webgl2', { antialias: false, premultipliedAlpha: false });
    if (!gl) throw new Error('WebGL2 not available');

    const prog = linkProgram(gl, VS_SRC, FS_SRC);
    gl.useProgram(prog);

    // Fullscreen quad.
    const vao = gl.createVertexArray();
    gl.bindVertexArray(vao);
    const vbo = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
    gl.bufferData(gl.ARRAY_BUFFER,
        new Float32Array([-1, -1,  1, -1, -1, 1,  -1, 1,  1, -1,  1, 1]),
        gl.STATIC_DRAW);
    const aPos = gl.getAttribLocation(prog, 'a_pos');
    gl.enableVertexAttribArray(aPos);
    gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

    // R32F field texture.
    const tex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R32F, GRID_W, GRID_H);

    const uField = gl.getUniformLocation(prog, 'u_field');
    const uScale = gl.getUniformLocation(prog, 'u_scale');
    gl.uniform1i(uField, 0);
    gl.uniform1f(uScale, DISPLAY_SCALE);

    let paused = false;
    pauseBtn.addEventListener('click', () => {
        paused = !paused;
        pauseBtn.textContent = paused ? 'resume' : 'pause';
    });

    function resetPulse(): void {
        const params = new Float32Array([PULSE_SIGMA, PULSE_AMP]);
        const ptr = M._malloc(params.byteLength);
        M.HEAPF32.set(params, ptr >> 2);
        api.loadScenario(api.sim, 'wave_pulse', ptr, params.length);
        M._free(ptr);
    }
    resetBtn.addEventListener('click', resetPulse);

    function frame(): void {
        if (!paused) api.stepN(api.sim, STEPS_PER_FRAME);

        const ptr = api.fieldPtr(api.sim, /* GR_FIELD_PHI_GRAV */ 0);
        // Float32Array view directly over WASM linear memory — no copy.
        const data = new Float32Array(M.HEAPF32.buffer, ptr, GRID_W * GRID_H);
        gl!.texSubImage2D(gl!.TEXTURE_2D, 0, 0, 0, GRID_W, GRID_H, gl!.RED, gl!.FLOAT, data);

        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
        gl!.drawArrays(gl!.TRIANGLES, 0, 6);

        statusEl.textContent =
            `step ${api.stepCount(api.sim).toString().padStart(5)}   ` +
            `t = ${api.simTime(api.sim).toFixed(2)}   ` +
            `grid ${GRID_W}x${GRID_H}   cfl = ${CFL.toFixed(4)}`;
        requestAnimationFrame(frame);
    }

    statusEl.textContent = 'running';
    frame();
}

main().catch((err) => {
    const s = document.getElementById('status');
    if (s) s.textContent = `error: ${err instanceof Error ? err.message : String(err)}`;
    console.error(err);
});
