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

const GRID_W = 384;
const GRID_H = 384;
const DX = 1.0;
const C_EFF = 1.0;
const CFL = 1.0 / Math.sqrt(2);
/* Damping layer chosen to absorb the long-wavelength quasi-static log tail
 * of Phi_g without strong reflection.  Polynomial m=2 with target_reflection
 * = 1e-3 is matched for n_damping = 48 (~12% of the box, ~3x deeper than
 * the typical Stage 02 layer). */
const N_DAMPING = 48;
const STEPS_PER_FRAME = 4;
/* v38 §15.9 "Field initialization": iterate the field leapfrog with the
 * particles' positions and velocities held fixed so the perturbation
 * fields converge to the static Poisson solution at the initial
 * configuration.  The transient needs several domain-crossing times to
 * propagate to the absorber and die.  One crossing at c_eff = 1 with
 * CFL = 1/sqrt(2) takes max(W,H) * sqrt(2) steps; ~10 crossings (≈14x
 * the largest dimension) is comfortably in the doc's 1000-3000-for-N=256
 * recommendation band. */
const WARMUP_CROSSINGS = 10;
const WARMUP_STEPS = Math.ceil(
    WARMUP_CROSSINGS * Math.max(GRID_W, GRID_H) * (1.0 / CFL),
);
/* Colormap gain on Phi_g.  Picked so the log-r far field is visible (not
 * clamped) and the near-zone cyan wells are still bright. */
const DISPLAY_SCALE = 5.0;

/* Per-scenario parameter packages. */
interface ScenarioSpec {
    name: string;
    params: number[];
}
const SCENARIOS: { [k: string]: ScenarioSpec } = {
    pic_binary: {
        name: 'pic_binary',
        /* params[0]=mass, params[1]=orbital radius, params[2]=v_factor */
        params: [0.01, 15.0 * DX, 1.0],
    },
    pic_static: {
        name: 'pic_static',
        /* params[0]=mass.  Use 0.01 so Phi_g is comparable in magnitude
         * to the pic_binary scenario for the same DISPLAY_SCALE. */
        params: [0.01],
    },
};
/* Kept for backwards compat in the diagnostic snapshot code below: the
 * "radius" used to sample at p0/p1 only makes sense for pic_binary. */
const BINARY_MASS = SCENARIOS.pic_binary.params[0];
const BINARY_R    = SCENARIOS.pic_binary.params[1];
const BINARY_VF   = SCENARIOS.pic_binary.params[2];

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
    setParticlesFrozen: (sim: number, frozen: number) => void;
    relaxPhiGPoisson: (sim: number, nIters: number) => void;
    initPotentialsLW: (sim: number) => void;
    setOuterBcNeumann: (sim: number, neumann: number) => void;
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
     * Emscripten's relative-URL .wasm fetch works correctly.
     *
     * Cache-bust the URL on every load so an out-of-date grlite.js
     * cannot get reused across `mingw32-make wasm` rebuilds during
     * dev.  Always force a fresh factory (no window-cache shortcut)
     * for the same reason. */
    const W = window as unknown as { __GRliteFactory?: () => Promise<GRliteModule> };
    W.__GRliteFactory = undefined;
    const cacheBuster = Date.now();
    await new Promise<void>((resolve, reject) => {
        const script = document.createElement('script');
        script.type = 'module';
        script.textContent =
            `import GRlite from '/grlite/grlite.js?v=${cacheBuster}';` +
            `window.__GRliteFactory = GRlite;` +
            `window.dispatchEvent(new Event('grlite-ready'));`;
        const onReady = () => { window.removeEventListener('grlite-error', onError); resolve(); };
        const onError = (e: Event) => { window.removeEventListener('grlite-ready', onReady); reject(new Error(String(e))); };
        window.addEventListener('grlite-ready', onReady, { once: true });
        window.addEventListener('grlite-error', onError, { once: true });
        script.onerror = (e) => reject(new Error(`script load failed: ${String(e)}`));
        document.head.appendChild(script);
    });
    const factory = W.__GRliteFactory;
    if (!factory) throw new Error('GRlite factory not set on window');
    return (factory as () => Promise<GRliteModule>)();
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
    const setParticlesFrozen = M.cwrap('gr_sim_set_particles_frozen', null,
        ['number','number']) as unknown as CFnVoid;
    const relaxPhiGPoisson = M.cwrap('gr_sim_relax_phi_g_poisson', null,
        ['number','number']) as unknown as CFnVoid;
    const initPotentialsLW = M.cwrap('gr_sim_init_potentials_lienard_wiechert',
        null, ['number']) as unknown as CFnVoid;
    const setOuterBcNeumann = M.cwrap('gr_sim_set_outer_bc_neumann', null,
        ['number','number']) as unknown as CFnVoid;
    /* Hard-fail if any expected symbol didn't resolve -- otherwise cwrap
     * silently returns a no-op and we get a confusing "is not a function"
     * deep inside the warmup later. */
    if (typeof initPotentialsLW !== 'function') {
        throw new Error('gr_sim_init_potentials_lienard_wiechert missing from WASM exports — rebuild grlite.wasm');
    }
    if (typeof relaxPhiGPoisson !== 'function') {
        throw new Error('gr_sim_relax_phi_g_poisson missing from WASM exports — rebuild grlite.wasm');
    }
    if (typeof setParticlesFrozen !== 'function') {
        throw new Error('gr_sim_set_particles_frozen missing from WASM exports — rebuild grlite.wasm');
    }
    const stepCount = M.cwrap('gr_sim_step_count', 'number', ['number']) as unknown as CFn;
    const simTime = M.cwrap('gr_sim_time', 'number', ['number']) as unknown as CFn;
    const particleCount = M.cwrap('gr_sim_particle_count', 'number', ['number']) as unknown as CFn;
    const getParticle = M.cwrap('gr_sim_get_particle', 'number', ['number','number']) as unknown as CFn;
    const clearParticles = M.cwrap('gr_sim_clear_particles', null, ['number']) as unknown as CFnVoid;

    const sim = create(GRID_W, GRID_H, DX, C_EFF, CFL);
    console.log('gr_sim_create returned sim=', sim);
    if (!sim) throw new Error('gr_sim_create returned NULL');

    /* Load the scenario FIRST -- pic_binary configures damping internally
     * (the May-17 WASM's set_damping has been observed to throw out-of-
     * bounds when called on a freshly-created sim before any scenario,
     * which suggests the damping allocator depends on scenario state). */
    const params = new Float32Array([BINARY_MASS, BINARY_R, BINARY_VF]);
    const paramsPtr = M._malloc(params.byteLength);
    console.log('alloc params at', paramsPtr);
    M.HEAPF32.set(params, paramsPtr >> 2);
    const rc = loadScenario(sim, 'pic_binary', paramsPtr, params.length);
    M._free(paramsPtr);
    console.log('load_scenario pic_binary returned rc=', rc);
    if (rc !== 0) {
        /* Fall back to wave_pulse if pic_binary isn't registered in this
         * older WASM. */
        const wp = new Float32Array([4.0 * DX, 1.0]);
        const wpPtr = M._malloc(wp.byteLength);
        M.HEAPF32.set(wp, wpPtr >> 2);
        const rc2 = loadScenario(sim, 'wave_pulse', wpPtr, wp.length);
        M._free(wpPtr);
        console.log('fallback load_scenario wave_pulse returned rc=', rc2);
        if (rc2 !== 0) throw new Error(`no scenario could be loaded (binary rc=${rc}, wave rc=${rc2})`);
    }

    /* Apply damping AFTER the scenario.  Wrap in try/catch in case of
     * compatibility issues. */
    try {
        setDamping(sim, N_DAMPING);
        console.log('setDamping ok');
    } catch (e) {
        console.warn('setDamping failed; continuing without damping ring:', e);
    }

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

    /* Enable Neumann outer BC by default for the web demo so that
     * outgoing waves reflect without sign flip and standing-mode
     * antinodes sit at the absorber instead of at the center. */
    setOuterBcNeumann(sim, 1);

    return { sim, step, stepN, fieldPtr, loadScenario, setDamping,
             setParticlesFrozen, relaxPhiGPoisson, initPotentialsLW,
             setOuterBcNeumann, stepCount, simTime, particleCount,
             getParticle, clearParticles, particleStrideF32: stride };
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
    const scenarioSel = document.getElementById('scenario') as HTMLSelectElement;

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
    /* R32F is not natively filterable in WebGL2 without OES_texture_float_linear;
     * use NEAREST sampling which is always supported. */
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R32F, GRID_W, GRID_H);

    gl.useProgram(fieldProg);
    gl.uniform1i(gl.getUniformLocation(fieldProg, 'u_field'), 0);
    const uScale = gl.getUniformLocation(fieldProg, 'u_scale');
    gl.uniform1f(uScale, DISPLAY_SCALE);

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
    /* Start paused so the post-warmup quasi-static field renders without
     * any orbital evolution on top -- press 'resume' to begin dynamics. */
    let paused = true;
    let singleStep = false;
    pauseBtn.textContent = 'resume';
    pauseBtn.addEventListener('click', () => {
        paused = !paused;
        pauseBtn.textContent = paused ? 'resume' : 'pause';
    });
    stepBtn.addEventListener('click', () => { singleStep = true; });
    function loadAndWarmup(specName: string): void {
        const spec = SCENARIOS[specName];
        if (!spec) { console.error(`unknown scenario: ${specName}`); return; }
        const arr = new Float32Array(spec.params);
        const ptr = M._malloc(arr.byteLength);
        M.HEAPF32.set(arr, ptr >> 2);
        const rc = api.loadScenario(api.sim, spec.name, ptr, arr.length);
        M._free(ptr);
        if (rc !== 0) { console.error(`load_scenario ${spec.name} failed rc=${rc}`); return; }
        /* Field initialization: direct-sum 2D Liénard-Wiechert potentials
         * for all six fields, with a damping-ring taper so the IC reaches
         * zero at the outer Dirichlet wall without launching a boundary
         * shock.  Sets prev = curr = next so the wave equation starts from
         * bit-exact zero time derivative.  This is the "true open-BC"
         * initial condition the wave equation would converge toward,
         * computed analytically rather than via iterative relaxation. */
        api.initPotentialsLW(api.sim);
    }

    /* Field-stability probe: run N extra steps with particles still
     * frozen and sample Phi_g(center) along the way.  For a converged
     * static Poisson solution with constant sources this should produce
     * a flat trajectory.  If it oscillates, the wave equation isn't
     * actually at a fixed point. */
    function fieldStabilityProbe(label: string, batches: number, perBatch: number): void {
        api.setParticlesFrozen(api.sim, 1);
        const ptr = api.fieldPtr(api.sim, /* GR_FIELD_PHI_GRAV */ 0);
        const phi = new Float32Array(M.HEAPF32.buffer, ptr, GRID_W * GRID_H);
        const cx = GRID_W >> 1, cy = GRID_H >> 1;
        const series: string[] = [];
        let pmin = Infinity, pmax = -Infinity;
        for (let b = 0; b < batches; b++) {
            const phiCenter = phi[cy * GRID_W + cx];
            let lmin = Infinity, lmax = -Infinity;
            for (let v of phi) { if (v < lmin) lmin = v; if (v > lmax) lmax = v; }
            if (lmin < pmin) pmin = lmin;
            if (lmax > pmax) pmax = lmax;
            series.push(`step+=${(b * perBatch).toString().padStart(5)}  center=${phiCenter.toExponential(3)}  range=[${lmin.toExponential(2)}, ${lmax.toExponential(2)}]`);
            api.stepN(api.sim, perBatch);
        }
        api.setParticlesFrozen(api.sim, 0);
        console.log(`[stability probe: ${label}] frozen, ${batches}x${perBatch}=${batches * perBatch} steps`);
        for (const s of series) console.log('  ' + s);
        console.log(`[stability probe: ${label}] envelope min=${pmin.toExponential(3)} max=${pmax.toExponential(3)} swing=${(pmax - pmin).toExponential(3)}`);
    }
    resetBtn.addEventListener('click', () => {
        loadAndWarmup(scenarioSel.value);
        paused = true;
        pauseBtn.textContent = 'resume';
    });
    scenarioSel.addEventListener('change', () => {
        statusEl.textContent = `loading ${scenarioSel.value}…`;
        requestAnimationFrame(() => {
            loadAndWarmup(scenarioSel.value);
            snapshotParticles(`post-warmup ${scenarioSel.value}`);
            snapshotPhi(`post-warmup ${scenarioSel.value}`);
            fieldStabilityProbe(
                `${scenarioSel.value} post-warmup, particles frozen`, 20, 200);
            paused = true;
            pauseBtn.textContent = 'resume';
            statusEl.textContent = `paused on converged ${scenarioSel.value} (press resume)`;
        });
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

        /* Auto-normalize the colormap to the current frame's |Phi_g| range
         * so a tiny converged Poisson field and a huge over-amplitude
         * transient both render with visible dynamic range. */
        let maxAbs = 0;
        for (let k = 0; k < data.length; k++) {
            const a = data[k] < 0 ? -data[k] : data[k];
            if (a > maxAbs) maxAbs = a;
        }
        gl!.useProgram(fieldProg);
        gl!.uniform1f(uScale, maxAbs > 1e-12 ? 1.0 / maxAbs : DISPLAY_SCALE);

        /* Render field. */
        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clearColor(0, 0, 0, 1);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
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

        let pmin = Infinity, pmax = -Infinity;
        for (let k = 0; k < data.length; k++) {
            const v = data[k];
            if (v < pmin) pmin = v;
            if (v > pmax) pmax = v;
        }
        statusEl.textContent =
            `${scenarioSel.value}  ` +
            `step ${api.stepCount(api.sim).toString().padStart(6)}  ` +
            `t=${api.simTime(api.sim).toFixed(2).padStart(7)}  ` +
            `Phi_g range [${pmin.toExponential(2)}, ${pmax.toExponential(2)}]  ` +
            `parts=${n}  ` +
            `${paused ? '[PAUSED]' : ''}`;
        requestAnimationFrame(frame);
    }

    /* --- Diagnostic helpers ------------------------------------------- */
    function snapshotParticles(label: string): void {
        const n = api.particleCount(api.sim);
        for (let i = 0; i < n; i++) {
            const p = new Float32Array(M.HEAPF32.buffer, api.getParticle(api.sim, i), 4);
            console.log(`[${label}] particle ${i}: x=${p[0].toFixed(4)} y=${p[1].toFixed(4)} px=${p[2].toFixed(4)} py=${p[3].toFixed(4)}`);
        }
    }
    function snapshotPhi(label: string): void {
        const ptr = api.fieldPtr(api.sim, /* GR_FIELD_PHI_GRAV */ 0);
        const phi = new Float32Array(M.HEAPF32.buffer, ptr, GRID_W * GRID_H);
        const at = (i: number, j: number) => phi[j * GRID_W + i];
        const cx = GRID_W >> 1;
        const cy = GRID_H >> 1;
        const r = Math.round(BINARY_R / DX);
        let max = -Infinity, min = Infinity;
        for (let v of phi) { if (v < min) min = v; if (v > max) max = v; }
        console.log(`[${label}] Phi_g min=${min.toExponential(3)} max=${max.toExponential(3)}`);
        console.log(`[${label}]   center=${at(cx, cy).toExponential(3)}`);
        console.log(`[${label}]   at p0 (cx-r,cy)=${at(cx - r, cy).toExponential(3)}`);
        console.log(`[${label}]   at p1 (cx+r,cy)=${at(cx + r, cy).toExponential(3)}`);
        console.log(`[${label}]   at (cx,cy+r)=${at(cx, cy + r).toExponential(3)} (perp axis)`);
        console.log(`[${label}]   at (16,16) corner=${at(16, 16).toExponential(3)}`);
        /* Radial profile along y=cy from i=0 to i=GRID_W-1, every 16 cells.
         * For a converged 2D-log Poisson solution this should be monotone:
         * ~0 at the absorber, growing more negative as we approach a
         * particle, two minima at i = cx +/- r_orb, less-negative bump
         * at the COM between them. */
        const profile: string[] = [];
        for (let i = 0; i < GRID_W; i += 16) {
            profile.push(`(i=${i.toString().padStart(3)}, v=${at(i, cy).toExponential(2)})`);
        }
        console.log(`[${label}]   radial profile y=cy: ${profile.join(' ')}`);
    }

    /* v38 §15.9 convergence iteration: freeze the particles in place
     * (positions AND velocities pinned to scenario IC), iterate the field
     * leapfrog until Phi_g, A_g, etc. converge to their static Poisson
     * solutions for those sources, then unfreeze and begin dynamics.  The
     * displayed field at t=0 of the visible run is already quasi-static --
     * no switch-on wavefront. */
    snapshotParticles('pre-init');
    snapshotPhi('pre-init');
    statusEl.textContent = `direct-sum L-W initialization…`;
    /* Yield to the browser so the status text actually paints, then run. */
    await new Promise<void>((r) => requestAnimationFrame(() => r()));
    api.initPotentialsLW(api.sim);
    snapshotParticles('post-init (should be IDENTICAL to pre)');
    snapshotPhi('post-init (2D L-W direct-sum + damping-ring taper)');

    /* Field-stability probe: confirm or refute whether the converged
     * field is actually a static fixed point of the wave equation. */
    fieldStabilityProbe('post-warmup, particles frozen', 20, 200);

    statusEl.textContent = 'running pic_binary';
    frame();
}

main().catch((err) => {
    const s = document.getElementById('status');
    if (s) s.textContent = `error: ${err instanceof Error ? err.message : String(err)}`;
    console.error(err);
});
