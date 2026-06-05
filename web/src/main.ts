/* GRlite web frontend.
 *
 * Loads the WASM core, instantiates a 2D linearized-GR + EM PIC sim, and
 * renders the gravitational potential Phi_g as a field texture plus the
 * orbiting particles as point sprites on top.  No physics here -- TS only
 * marshals parameters in and reads pointers out per memory `grlite-stack`.
 *
 * Default scenario: pic_binary (two equal masses in circular orbit via
 * FDTD gravity, gr_sandbox_v32.tex Stage 12).
 */

import type { GRliteModule } from './grlite';

/* ---- Sim configuration -------------------------------------------------- */

const GRID_W = 256;
const GRID_H = 256;
const DX = 1.0;
const C_EFF = 1.0;
const CFL = 1.0 / Math.sqrt(2);
const N_DAMPING = 16;
const STEPS_PER_FRAME = 4;
const DISPLAY_SCALE = 30.0;            /* gain on Phi_g for the colormap */

/* Binary scenario parameters (pic_binary):
 *   params[0] = mass per particle
 *   params[1] = orbital radius
 *   params[2] = v_factor (1.0 = circular for the continuum 2D-log)        */
const BINARY_MASS = 0.01;
const BINARY_R    = 15.0 * DX;
const BINARY_VF   = 1.0;

/* ---- Particle struct layout (gr_particle_t in grlite.h) ----------------
 *  float x, y, px, py, mass, charge, proper_time, [spin, phi_spin, g_factor]
 * The first 7 floats are always there; v40+ adds the three spin floats.
 * We accept either size at runtime. */
const PARTICLE_STRIDE_F32_OLD = 7;     /* before v40 */
const PARTICLE_STRIDE_F32_NEW = 10;    /* v40 and later */

/* ---- TypeScript bindings to the WASM API -------------------------------- */

interface SimAPI {
    sim: number;
    step: (sim: number) => void;
    stepN: (sim: number, n: number) => void;
    fieldPtr: (sim: number, which: number) => number;
    loadScenario: (sim: number, name: string, paramsPtr: number, n: number) => number;
    setDamping: (sim: number, nDamping: number) => void;
    stepCount: (sim: number) => number;
    simTime: (sim: number) => number;
    particleCount: (sim: number) => number;
    getParticle: (sim: number, idx: number) => number;
    clearParticles: (sim: number) => void;
    /* Detected at runtime: byte-size of one particle struct. */
    particleStrideF32: number;
}

async function loadWasm(): Promise<GRliteModule> {
    /* Vite 6 forbids static imports from /public/, so we load grlite.js
     * by injecting a <script type="module"> into the DOM and stashing
     * the factory on the window.  This bypasses Vite's import analyzer
     * entirely and keeps grlite.js + grlite.wasm in /public/ where
     * Emscripten's relative-URL .wasm fetch works correctly. */
    const W = window as unknown as { __GRliteFactory?: () => Promise<GRliteModule> };
    if (!W.__GRliteFactory) {
        await new Promise<void>((resolve, reject) => {
            const script = document.createElement('script');
            script.type = 'module';
            script.textContent =
                `import GRlite from '/grlite/grlite.js';` +
                `window.__GRliteFactory = GRlite;` +
                `window.dispatchEvent(new Event('grlite-ready'));`;
            const onReady = () => { window.removeEventListener('grlite-error', onError); resolve(); };
            const onError = (e: Event) => { window.removeEventListener('grlite-ready', onReady); reject(new Error(String(e))); };
            window.addEventListener('grlite-ready', onReady, { once: true });
            window.addEventListener('grlite-error', onError, { once: true });
            script.onerror = (e) => reject(new Error(`script load failed: ${String(e)}`));
            document.head.appendChild(script);
        });
    }
    const factory = W.__GRliteFactory;
    if (!factory) throw new Error('GRlite factory not set on window');
    return factory();
}

function bindApi(M: GRliteModule): SimAPI {
    type CFn = (...args: number[]) => number;
    type CFnVoid = (...args: number[]) => void;
    const create = M.cwrap('gr_sim_create', 'number',
        ['number','number','number','number','number']) as unknown as
        (w: number, h: number, dx: number, c: number, cfl: number) => number;
    const step = M.cwrap('gr_sim_step', null, ['number']) as unknown as CFnVoid;
    const stepN = M.cwrap('gr_sim_step_n', null, ['number','number']) as unknown as CFnVoid;
    const fieldPtr = M.cwrap('gr_sim_field_ptr', 'number', ['number','number']) as unknown as CFn;
    const loadScenario = M.cwrap('gr_sim_load_scenario', 'number',
        ['number','string','number','number']) as unknown as
        (sim: number, name: string, paramsPtr: number, n: number) => number;
    const setDamping = M.cwrap('gr_sim_set_damping', null,
        ['number','number']) as unknown as CFnVoid;
    const stepCount = M.cwrap('gr_sim_step_count', 'number', ['number']) as unknown as CFn;
    const simTime = M.cwrap('gr_sim_time', 'number', ['number']) as unknown as CFn;
    const particleCount = M.cwrap('gr_sim_particle_count', 'number', ['number']) as unknown as CFn;
    const getParticle = M.cwrap('gr_sim_get_particle', 'number', ['number','number']) as unknown as CFn;
    const clearParticles = M.cwrap('gr_sim_clear_particles', null, ['number']) as unknown as CFnVoid;

    const sim = create(GRID_W, GRID_H, DX, C_EFF, CFL);
    if (!sim) throw new Error('gr_sim_create returned NULL');
    setDamping(sim, N_DAMPING);

    /* Load the pic_binary scenario.  The scenario clears existing particles
     * and sets the two-body initial conditions internally. */
    const params = new Float32Array([BINARY_MASS, BINARY_R, BINARY_VF]);
    const paramsPtr = M._malloc(params.byteLength);
    M.HEAPF32.set(params, paramsPtr >> 2);
    const rc = loadScenario(sim, 'pic_binary', paramsPtr, params.length);
    M._free(paramsPtr);
    if (rc !== 0) throw new Error(`load_scenario 'pic_binary' failed: rc=${rc}`);

    /* Probe particle stride: if any particle exists, distance between
     * successive particle pointers reveals the struct stride. */
    let stride = PARTICLE_STRIDE_F32_NEW;
    const n = particleCount(sim);
    if (n >= 2) {
        const p0 = getParticle(sim, 0);
        const p1 = getParticle(sim, 1);
        if (p1 > p0) stride = (p1 - p0) >> 2;
    }
    if (stride !== PARTICLE_STRIDE_F32_OLD && stride !== PARTICLE_STRIDE_F32_NEW) {
        console.warn(`unexpected particle stride: ${stride} floats`);
    }

    return { sim, step, stepN, fieldPtr, loadScenario, setDamping,
             stepCount, simTime, particleCount, getParticle, clearParticles,
             particleStrideF32: stride };
}

/* ---- WebGL2 plumbing ---------------------------------------------------- */

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

function linkProgram(gl: WebGL2RenderingContext, vs: string, fs: string): WebGLProgram {
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

/* Field shader: fullscreen quad with diverging colormap on R32F texture. */
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
    /* Diverging colormap: cyan (negative) -> black (zero) -> orange (positive). */
    vec3 c = vec3(pos * 1.0, pos * 0.5 + neg * 0.6, neg * 1.0);
    outColor = vec4(c, 1.0);
}`;

/* Particle shader: each particle is a point sprite. */
const POINT_VS = `#version 300 es
in vec2 a_pos;                              /* in sim cells, [0, GRID_W) x [0, GRID_H) */
uniform vec2 u_gridSize;
uniform float u_pointSize;
void main() {
    vec2 ndc = (a_pos / u_gridSize) * 2.0 - 1.0;
    /* Texture is sampled with origin top-left in our texSubImage2D upload;
     * keep particles in the same convention. */
    gl_Position = vec4(ndc, 0.0, 1.0);
    gl_PointSize = u_pointSize;
}`;

const POINT_FS = `#version 300 es
precision highp float;
out vec4 outColor;
uniform vec3 u_color;
void main() {
    /* Round particle: alpha fades to zero at the disc edge for smoother look. */
    vec2 d = gl_PointCoord - vec2(0.5);
    float r2 = dot(d, d);
    if (r2 > 0.25) discard;
    float edge = smoothstep(0.25, 0.16, r2);
    outColor = vec4(u_color, edge);
}`;

/* ---- Main loop ---------------------------------------------------------- */

async function main(): Promise<void> {
    const statusEl = document.getElementById('status') as HTMLDivElement;
    const canvas = document.getElementById('view') as HTMLCanvasElement;
    const resetBtn = document.getElementById('reset') as HTMLButtonElement;
    const pauseBtn = document.getElementById('pause') as HTMLButtonElement;
    const stepBtn  = document.getElementById('stepFrame') as HTMLButtonElement;

    statusEl.textContent = 'loading WASM…';
    const M = await loadWasm();
    const api = bindApi(M);

    const gl = canvas.getContext('webgl2', { antialias: true, premultipliedAlpha: false });
    if (!gl) throw new Error('WebGL2 not available');

    /* Field rendering program */
    const fieldProg = linkProgram(gl, FIELD_VS, FIELD_FS);
    const quadVao = gl.createVertexArray();
    gl.bindVertexArray(quadVao);
    const quadVbo = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quadVbo);
    gl.bufferData(gl.ARRAY_BUFFER,
        new Float32Array([-1,-1, 1,-1, -1,1, -1,1, 1,-1, 1,1]),
        gl.STATIC_DRAW);
    const aPosField = gl.getAttribLocation(fieldProg, 'a_pos');
    gl.enableVertexAttribArray(aPosField);
    gl.vertexAttribPointer(aPosField, 2, gl.FLOAT, false, 0, 0);

    const fieldTex = gl.createTexture();
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, fieldTex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R32F, GRID_W, GRID_H);

    gl.useProgram(fieldProg);
    gl.uniform1i(gl.getUniformLocation(fieldProg, 'u_field'), 0);
    gl.uniform1f(gl.getUniformLocation(fieldProg, 'u_scale'), DISPLAY_SCALE);

    /* Particle rendering program */
    const pointProg = linkProgram(gl, POINT_VS, POINT_FS);
    const pointVao = gl.createVertexArray();
    gl.bindVertexArray(pointVao);
    const pointVbo = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, pointVbo);
    const aPosPoint = gl.getAttribLocation(pointProg, 'a_pos');
    gl.enableVertexAttribArray(aPosPoint);
    gl.vertexAttribPointer(aPosPoint, 2, gl.FLOAT, false, 0, 0);
    gl.useProgram(pointProg);
    gl.uniform2f(gl.getUniformLocation(pointProg, 'u_gridSize'),
                 GRID_W, GRID_H);
    gl.uniform1f(gl.getUniformLocation(pointProg, 'u_pointSize'),
                 Math.min(canvas.width, canvas.height) * 0.018);
    gl.uniform3f(gl.getUniformLocation(pointProg, 'u_color'), 1.0, 1.0, 1.0);

    /* Re-usable float buffer for particle XYs (max 64 particles for now). */
    const PARTICLE_MAX = 64;
    const xyBuf = new Float32Array(PARTICLE_MAX * 2);

    /* ---- Controls ---- */
    let paused = false;
    let singleStep = false;
    pauseBtn.addEventListener('click', () => {
        paused = !paused;
        pauseBtn.textContent = paused ? 'resume' : 'pause';
    });
    stepBtn.addEventListener('click', () => { singleStep = true; });
    resetBtn.addEventListener('click', () => {
        const params = new Float32Array([BINARY_MASS, BINARY_R, BINARY_VF]);
        const ptr = M._malloc(params.byteLength);
        M.HEAPF32.set(params, ptr >> 2);
        api.loadScenario(api.sim, 'pic_binary', ptr, params.length);
        M._free(ptr);
    });

    function frame(): void {
        if (!paused || singleStep) {
            api.stepN(api.sim, singleStep ? 1 : STEPS_PER_FRAME);
            singleStep = false;
        }

        /* Update field texture from WASM heap (zero-copy view). */
        const ptr = api.fieldPtr(api.sim, /* GR_FIELD_PHI_GRAV */ 0);
        const data = new Float32Array(M.HEAPF32.buffer, ptr, GRID_W * GRID_H);
        gl!.bindTexture(gl!.TEXTURE_2D, fieldTex);
        gl!.texSubImage2D(gl!.TEXTURE_2D, 0, 0, 0, GRID_W, GRID_H,
                          gl!.RED, gl!.FLOAT, data);

        /* Render field. */
        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clearColor(0, 0, 0, 1);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
        gl!.useProgram(fieldProg);
        gl!.bindVertexArray(quadVao);
        gl!.activeTexture(gl!.TEXTURE0);
        gl!.bindTexture(gl!.TEXTURE_2D, fieldTex);
        gl!.drawArrays(gl!.TRIANGLES, 0, 6);

        /* Render particles on top. */
        const n = Math.min(api.particleCount(api.sim), PARTICLE_MAX);
        for (let i = 0; i < n; i++) {
            const pptr = api.getParticle(api.sim, i);
            const view = new Float32Array(M.HEAPF32.buffer, pptr, 2);
            xyBuf[2 * i + 0] = view[0];
            xyBuf[2 * i + 1] = view[1];
        }
        if (n > 0) {
            gl!.enable(gl!.BLEND);
            gl!.blendFunc(gl!.SRC_ALPHA, gl!.ONE_MINUS_SRC_ALPHA);
            gl!.useProgram(pointProg);
            gl!.bindVertexArray(pointVao);
            gl!.bindBuffer(gl!.ARRAY_BUFFER, pointVbo);
            gl!.bufferData(gl!.ARRAY_BUFFER,
                xyBuf.subarray(0, n * 2),
                gl!.DYNAMIC_DRAW);
            gl!.drawArrays(gl!.POINTS, 0, n);
            gl!.disable(gl!.BLEND);
        }

        statusEl.textContent =
            `step ${api.stepCount(api.sim).toString().padStart(6)}  ` +
            `t = ${api.simTime(api.sim).toFixed(2).padStart(8)}  ` +
            `particles = ${n}  ` +
            `grid ${GRID_W}x${GRID_H}  cfl=${CFL.toFixed(4)}`;
        requestAnimationFrame(frame);
    }

    statusEl.textContent = 'running pic_binary';
    frame();
}

main().catch((err) => {
    const s = document.getElementById('status');
    if (s) s.textContent = `error: ${err instanceof Error ? err.message : String(err)}`;
    console.error(err);
});
