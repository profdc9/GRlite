/* WASM loading + the raw cwrap surface.
 *
 * `loadCore()` loads the Emscripten module and binds every exported
 * gr_sim_* function into a typed `GRliteCore`.  All functions take the
 * sim handle explicitly; the higher-level `World` wrapper (sim/world.ts)
 * owns a handle and presents an ergonomic API.  No physics or rendering
 * here -- this is purely the FFI boundary. */

import type { GRliteModule } from '../grlite';

type CFn = (...args: number[]) => number;
type CFnVoid = (...args: number[]) => void;

export interface GRliteCore {
    M: GRliteModule;
    create: (w: number, h: number, dx: number, c: number, cfl: number) => number;
    step: (sim: number) => void;
    stepN: (sim: number, n: number) => void;
    fieldPtr: (sim: number, which: number) => number;
    loadScenario: (sim: number, name: string, paramsPtr: number, n: number) => number;
    setDamping: (sim: number, nDamping: number) => void;
    setParticlesFrozen: (sim: number, frozen: number) => void;
    relaxPhiGPoisson: (sim: number, nIters: number) => void;
    initPotentialsLW: (sim: number) => void;
    initPotentialsLWTaper: (sim: number, taperInner: number, taperOuter: number) => void;
    initPotentialsSettled: (sim: number, nSettle: number) => void;
    setShapeFunction: (sim: number, shape: number) => void;
    setKernelRadius: (sim: number, radiusCells: number) => void;
    setForceInterp: (sim: number, scheme: number) => void;
    setOuterBcNeumann: (sim: number, neumann: number) => void;
    setVolumeFrictionTaper: (sim: number, uniform: number, taperMax: number, taperDepth: number) => void;
    setZeroMeanScalarPotentials: (sim: number, enabled: number) => void;
    setGEff: (sim: number, gEff: number) => void;
    setKE: (sim: number, kE: number) => void;
    addParticle: (sim: number, x: number, y: number, mass: number, charge: number, vx: number, vy: number) => number;
    clearParticles: (sim: number) => void;
    clearSources: (sim: number) => void;
    clearBackground: (sim: number) => void;
    stepCount: (sim: number) => number;
    simTime: (sim: number) => number;
    simDt: (sim: number) => number;
    particleCount: (sim: number) => number;
    getParticle: (sim: number, idx: number) => number;
    /* Detected at runtime: floats per gr_particle_t. */
    particleStrideF32: number;
}

export async function loadWasmModule(): Promise<GRliteModule> {
    /* Vite forbids static imports from /public/, so load grlite.js by
     * injecting a <script type="module"> and stashing the factory on
     * window.  Cache-bust so a `mingw32-make wasm` rebuild is always
     * picked up during dev. */
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

function bindCore(M: GRliteModule): GRliteCore {
    const num = (name: string, args: string[]): CFn =>
        M.cwrap(name, 'number', args) as unknown as CFn;
    const vfn = (name: string, args: string[]): CFnVoid =>
        M.cwrap(name, null, args) as unknown as CFnVoid;

    const core: GRliteCore = {
        M,
        create: M.cwrap('gr_sim_create', 'number',
            ['number','number','number','number','number']) as unknown as
            (w: number, h: number, dx: number, c: number, cfl: number) => number,
        step: vfn('gr_sim_step', ['number']),
        stepN: vfn('gr_sim_step_n', ['number','number']),
        fieldPtr: num('gr_sim_field_ptr', ['number','number']),
        loadScenario: M.cwrap('gr_sim_load_scenario', 'number',
            ['number','string','number','number']) as unknown as
            (sim: number, name: string, p: number, n: number) => number,
        setDamping: vfn('gr_sim_set_damping', ['number','number']),
        setParticlesFrozen: vfn('gr_sim_set_particles_frozen', ['number','number']),
        relaxPhiGPoisson: vfn('gr_sim_relax_phi_g_poisson', ['number','number']),
        initPotentialsLW: vfn('gr_sim_init_potentials_lienard_wiechert', ['number']),
        initPotentialsLWTaper: vfn('gr_sim_init_potentials_lienard_wiechert_with_taper', ['number','number','number']),
        initPotentialsSettled: vfn('gr_sim_init_potentials_settled', ['number','number']),
        setShapeFunction: vfn('gr_sim_set_shape_function', ['number','number']),
        setKernelRadius: vfn('gr_sim_set_kernel_radius', ['number','number']),
        setForceInterp: vfn('gr_sim_set_force_interp', ['number','number']),
        setOuterBcNeumann: vfn('gr_sim_set_outer_bc_neumann', ['number','number']),
        setVolumeFrictionTaper: vfn('gr_sim_set_volume_friction_taper', ['number','number','number','number']),
        setZeroMeanScalarPotentials: vfn('gr_sim_set_zero_mean_scalar_potentials', ['number','number']),
        setGEff: vfn('gr_sim_set_G_eff', ['number','number']),
        setKE: vfn('gr_sim_set_k_e', ['number','number']),
        addParticle: M.cwrap('gr_sim_add_particle', 'number',
            ['number','number','number','number','number','number','number']) as unknown as
            (sim: number, x: number, y: number, m: number, q: number, vx: number, vy: number) => number,
        clearParticles: vfn('gr_sim_clear_particles', ['number']),
        clearSources: vfn('gr_sim_clear_sources', ['number']),
        clearBackground: vfn('gr_sim_clear_background', ['number']),
        stepCount: num('gr_sim_step_count', ['number']),
        simTime: num('gr_sim_time', ['number']),
        simDt: num('gr_sim_dt', ['number']),
        particleCount: num('gr_sim_particle_count', ['number']),
        getParticle: num('gr_sim_get_particle', ['number','number']),
        particleStrideF32: 10,
    };

    /* Fail loudly if a load-bearing symbol is missing -- otherwise cwrap
     * returns a silent no-op that surfaces as a confusing error later. */
    const required: (keyof GRliteCore)[] = [
        'create','step','stepN','fieldPtr','loadScenario','initPotentialsSettled',
        'setShapeFunction','setKernelRadius','setVolumeFrictionTaper',
        'setOuterBcNeumann','setParticlesFrozen','addParticle','getParticle',
    ];
    for (const name of required) {
        if (typeof core[name] !== 'function') {
            throw new Error(`gr_sim symbol for "${String(name)}" missing from WASM exports — rebuild grlite.wasm`);
        }
    }
    return core;
}

export async function loadCore(): Promise<GRliteCore> {
    const M = await loadWasmModule();
    return bindCore(M);
}
