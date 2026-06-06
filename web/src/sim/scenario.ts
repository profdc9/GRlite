/* GRlite scenario JSON schema — the source of truth for a web simulation.
 *
 * A scenario fully and reproducibly determines a run: grid, all global
 * physics/absorber/init config + switches, background, and particles.  The
 * sim is built from this via the WASM API (sim/build.ts); the C scenarios
 * remain only for native tests. */

export const SCENARIO_PROGRAM = 'grlite';
export const SCENARIO_FORMAT = 'grlite-scenario';
export const SCENARIO_VERSION = 1;

export type ShapeName = 'cic' | 'tsc' | 'bump';
export type ForceInterpName = 'legacy' | 'lewis-birdsall';
export type OuterBCName = 'dirichlet' | 'neumann';
export type InitMethod = 'lw-settled' | 'lw' | 'none';

export interface GridSpec {
    W: number; H: number; dx: number; cEff: number; cfl: number;
}

export interface Switches {
    gravitomagneticForce: boolean;
    gravitomagneticInductive: boolean;
    emLorentz: boolean;
    emInductive: boolean;
    emElectrostatic: boolean;
    emMagnetic: boolean;
    emStressEnergy: boolean;
    emShapiro: boolean;
    fieldEvolution: boolean;
    particleSourceDeposition: boolean;
    esirkepov: boolean;
    periodicBC: boolean;
}

export interface AbsorberSpec {
    outerBC: OuterBCName;
    frictionFloor: number;
    frictionWall: number;
    frictionDepth: number;
    zeroMeanScalar: boolean;
}

export interface InitSpec { method: InitMethod; settleSteps: number; }

export interface GlobalSpec {
    gEff: number;
    kE: number;
    shape: ShapeName;
    kernelRadius: number;
    forceInterp: ForceInterpName;
    absorber: AbsorberSpec;
    init: InitSpec;
    switches: Switches;
    rhoSmooth: number;
    jSmooth: number;
    /* Self-field subtraction (anti-heating) applied to every particle unless
     * the particle overrides via its own `selfField`.  Default on. */
    selfFieldDefault: boolean;
    /* No-radiation mode: forces the inductive (radiation-reaction) pieces
     * off (emInductive + gravitomagneticInductive), to separate genuine
     * radiation-reaction inspiral from numerical heating when debugging. */
    noRadiation: boolean;
}

export interface BackgroundSpec {
    type: 'point-mass' | 'spinning-mass';
    x: number; y: number; GM: number; epsilon: number; Jz?: number;
}

export interface ParticleSpec {
    x: number; y: number; vx: number; vy: number;
    mass: number; charge: number;
    spin?: number; gFactor?: number;
    /* v39 per-particle self-field subtraction: the particle gets its own
     * field set and feels only the OTHER particles' fields (bit-exact
     * self-force cancellation at eps=0).  Removes PIC self-heating at the
     * cost of ~one extra field solve per opted-in particle per step. */
    selfField?: boolean;
    selfFieldEps?: [number, number];
}

export interface ViewSpec { field: number; showTrails: boolean; showVelocity: boolean; }

export interface Scenario {
    program: typeof SCENARIO_PROGRAM;
    format: typeof SCENARIO_FORMAT;
    version: number;
    name: string;
    grid: GridSpec;
    global: GlobalSpec;
    background: BackgroundSpec[];
    particles: ParticleSpec[];
    view: ViewSpec;
}

/* ---- defaults ---- */

export function defaultSwitches(): Switches {
    return {
        gravitomagneticForce: true,
        gravitomagneticInductive: false,
        emLorentz: true,
        emInductive: true,
        emElectrostatic: true,
        emMagnetic: true,
        emStressEnergy: false,
        emShapiro: false,
        fieldEvolution: true,
        particleSourceDeposition: true,
        esirkepov: true,
        periodicBC: false,
    };
}

/* v42 production absorber/deposition defaults (see web/README + addendum). */
export function defaultGlobal(): GlobalSpec {
    return {
        gEff: 1.0,
        kE: 1.0,
        shape: 'bump',
        kernelRadius: 8,
        forceInterp: 'lewis-birdsall',
        absorber: { outerBC: 'dirichlet', frictionFloor: 0.0, frictionWall: 0.02,
                    frictionDepth: 48, zeroMeanScalar: false },
        init: { method: 'lw-settled', settleSteps: 384 },
        switches: defaultSwitches(),
        rhoSmooth: 0,
        jSmooth: 0,
        selfFieldDefault: true,
        noRadiation: false,
    };
}

export function defaultGrid(): GridSpec {
    return { W: 384, H: 384, dx: 1.0, cEff: 1.0, cfl: 1.0 / Math.sqrt(2) };
}

export function emptyScenario(name = 'untitled'): Scenario {
    return {
        program: SCENARIO_PROGRAM, format: SCENARIO_FORMAT, version: SCENARIO_VERSION,
        name,
        grid: defaultGrid(),
        global: defaultGlobal(),
        background: [],
        particles: [],
        view: { field: 0, showTrails: true, showVelocity: true },
    };
}

/* Scenarios are NOT defined in code -- they live as JSON files under
 * public/scenes/ and are the single source of truth (loaded via fetch in
 * main.ts).  emptyScenario() above is only a placeholder / starting point
 * for the Phase-3 editor; it is not a stand-in for a real scenario. */

/* ---- validation / migration ---- */

export function validate(obj: unknown): Scenario {
    const o = obj as Partial<Scenario>;
    if (!o || o.program !== SCENARIO_PROGRAM || o.format !== SCENARIO_FORMAT) {
        throw new Error('not a grlite scenario');
    }
    if (typeof o.version !== 'number' || o.version > SCENARIO_VERSION) {
        throw new Error(`unsupported scenario version ${o.version}`);
    }
    /* Fill any missing sections with defaults (forward-compatible load). */
    const base = emptyScenario(o.name ?? 'loaded');
    return {
        ...base,
        ...o,
        grid: { ...base.grid, ...(o.grid ?? {}) },
        global: {
            ...base.global, ...(o.global ?? {}),
            absorber: { ...base.global.absorber, ...(o.global?.absorber ?? {}) },
            init: { ...base.global.init, ...(o.global?.init ?? {}) },
            switches: { ...base.global.switches, ...(o.global?.switches ?? {}) },
        },
        background: o.background ?? [],
        particles: o.particles ?? [],
        view: { ...base.view, ...(o.view ?? {}) },
    } as Scenario;
}
