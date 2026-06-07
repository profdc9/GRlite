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
/* Field initialization: 'lw-settled' = Lienard-Wiechert + frozen-friction
 * settle to the discrete fixed point; 'lw' = direct-sum L-W only; 'none' = no
 * explicit init (auto: settle if perturbation is active, else leave zero);
 * 'zero' = always leave the field at zero (e.g. a wave-emission/antenna demo
 * that must start from a quiet field even with perturbation on). */
export type InitMethod = 'lw-settled' | 'lw' | 'none' | 'zero';
export type BgModeName = 'sampled' | 'analytic';
/* Which part of a field the visualizations (colormap + vector overlay) show.
 * (Turning the colormap OFF entirely is the 'none' entry of FIELD_VIEWS, not a
 * source -- the source picks WHICH data to draw when something is drawn.) */
export type FieldSourceName = 'perturbation' | 'background' | 'sum';
/* Per-particle trajectory tick markers (§17 Trajectory tracks): circles at
 * equal coordinate-time intervals, triangles at equal proper-time intervals.
 * Showing both makes relativistic time dilation visible as the proper-time
 * ticks spreading out along the track relative to the coordinate-time ones. */
export type TickMode = 'none' | 'coordinate' | 'proper' | 'both';

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
    /* Rebuild the Shapiro c_local each step from the TOTAL Phi_g (background +
     * deposited-mass perturbation), so a moving/deposited mass lenses the EM
     * wave (a 2D log r lens).  Needs emShapiro.  Default off. */
    shapiroDynamic: boolean;
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
    /* How a background field (background[]) is felt by particles: 'sampled'
     * interpolates the field off the grid; 'analytic' evaluates the closed
     * form at the particle's exact position (cleaner for Lense-Thirring /
     * frame-dragging tests where grid sampling smears a sharp source). */
    bgMode: BgModeName;
    /* Self-field subtraction (anti-heating) applied to every particle unless
     * the particle overrides via its own `selfField`.  Default on. */
    selfFieldDefault: boolean;
    /* No-radiation mode: forces the inductive (radiation-reaction) pieces
     * off (emInductive + gravitomagneticInductive), to separate genuine
     * radiation-reaction inspiral from numerical heating when debugging. */
    noRadiation: boolean;
}

/* A single unified "compact body" background, parameterized by the three
 * no-hair quantities.  The classic metrics are special cases:
 *   GM only            -> Schwarzschild
 *   GM + Q             -> Reissner-Nordstrom
 *   GM + Jz            -> Kerr
 *   GM + Q + Jz        -> Kerr-Newman
 * (linearized GEM+EM analog -- see gr_sim_set_background_body).  A scenario
 * holds 0 or 1 of these. */
export interface BackgroundSpec {
    x: number; y: number; epsilon: number;
    GM: number;   // gravitational mass  -> Phi_g
    Q: number;    // electric charge     -> phi_em (Coulomb)
    Jz: number;   // angular momentum    -> A_g (frame dragging)
    /* Gyromagnetic ratio of the spinning charge -> EM magnetic dipole A_em
     * from moment mu = gFactor*(Q/2M)*Jz.  Default 2 (Kerr-Newman); 0 = no
     * EM magnetic field.  Needs GM != 0 (moment undefined for a massless body). */
    gFactor: number;
}

/* The classic-metric name implied by which hairs are nonzero (display only). */
export function backgroundLabel(b: BackgroundSpec | undefined): string {
    if (!b || (b.GM === 0 && b.Q === 0 && b.Jz === 0)) return 'none';
    if (b.Q !== 0 && b.Jz !== 0) return 'Kerr–Newman';
    if (b.Jz !== 0) return 'Kerr';
    if (b.Q !== 0) return 'Reissner–Nordström';
    return 'Schwarzschild';
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
    /* Visualization-only (no physics effect): which trajectory ticks to draw
     * for this particle, and whether to show its proper-time clock.  Edited
     * live (no rebuild) but persisted in the scenario so links reproduce them. */
    ticks?: TickMode;
    clock?: boolean;
    /* Driven oscillating-current source (e.g. an EM dipole antenna): adds a
     * sinusoidal velocity dv = amp*omega*cos(omega t + phase) along `axis` each
     * step (continuity-safe; applied in the pusher).  Displacement amplitude is
     * `amp`.  Use phase = pi on one of a +q/-q pair for an antiphase dipole. */
    drive?: { amp: number; omega: number; phase?: number; axis: [number, number] };
    /* Omit physical forces on this particle (pure source / pinned).  Default
     * true (forces applied; a drive, if any, rides on top). */
    forces?: boolean;
}

export interface ViewSpec {
    field: number; showTrails: boolean; showVelocity: boolean;
    /* Field-visualization source (colormap + vector overlay). */
    source: FieldSourceName;
    /* Vector-arrow overlay: index into VECTOR_FIELDS (0 = off), and the grid
     * spacing in cells between arrows (smaller = denser). */
    vectorField: number;
    vectorSpacing: number;
    /* Trajectory-tick spacing Δ (same value for coordinate- and proper-time
     * ticks, so their differing spacing along the track reveals time dilation).
     * 0 = auto (pick a "nice" Δ giving ~12 ticks across the current trail). */
    tickInterval: number;
    /* Per-particle force-arrow overlay (§17): which force components to draw on
     * every particle.  Arrow lengths share one auto-scale so magnitudes are
     * comparable across components and particles. */
    forceArrows: ForceArrows;
}

export interface ForceArrows {
    grav: boolean;   // gravitational / GEM force (orange)
    em: boolean;     // EM Lorentz force (blue)
    spin: boolean;   // spin-gradient force (green)
    total: boolean;  // total applied force (white)
}

export interface Scenario {
    program: typeof SCENARIO_PROGRAM;
    format: typeof SCENARIO_FORMAT;
    version: number;
    name: string;
    /* Optional human-readable note shown in the inspector (documents a preset). */
    description?: string;
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
        shapiroDynamic: false,
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
        bgMode: 'sampled',
        selfFieldDefault: true,
        noRadiation: false,
    };
}

export function defaultGrid(): GridSpec {
    return { W: 384, H: 384, dx: 1.0, cEff: 1.0, cfl: 1.0 / Math.sqrt(2) };
}

/* Default unified compact body (centered, mass-free).  x/y default to the grid
 * center so a body added without coordinates lands sensibly. */
export function defaultBackground(W: number, H: number): BackgroundSpec {
    return { x: W * 0.5, y: H * 0.5, epsilon: 8, GM: 0, Q: 0, Jz: 0, gFactor: 2 };
}

export function emptyScenario(name = 'untitled'): Scenario {
    return {
        program: SCENARIO_PROGRAM, format: SCENARIO_FORMAT, version: SCENARIO_VERSION,
        name,
        grid: defaultGrid(),
        global: defaultGlobal(),
        background: [],
        particles: [],
        view: { field: 0, showTrails: true, showVelocity: true,
                source: 'perturbation', vectorField: 0, vectorSpacing: 24,
                tickInterval: 0,
                forceArrows: { grav: false, em: false, spin: false, total: false } },
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
        /* Spread-merge each background over the default body, so every JSON key
         * maps straight to the object (and a new BackgroundSpec field round-trips
         * without editing validate) -- same pattern as global/view above. */
        background: (o.background ?? []).map(
            (b) => ({ ...defaultBackground(base.grid.W, base.grid.H), ...b })),
        particles: o.particles ?? [],
        view: {
            ...base.view, ...(o.view ?? {}),
            /* deep-merge so a partial hand-authored forceArrows still gets the
             * other three components defaulted to false (like switches above). */
            forceArrows: { ...base.view.forceArrows, ...(o.view?.forceArrows ?? {}) },
        },
    } as Scenario;
}
