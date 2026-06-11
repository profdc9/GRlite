/* Editable inspectors (Phase 3).
 *
 * The panels edit the canonical Scenario; committing an edit goes through
 * onEdit(), which (in main) snapshots for undo and rebuilds from the JSON
 * (the deterministic edit->rebuild model).  Editable inputs are rebuilt only
 * on scenario/selection change (renderEditors), so typing isn't interrupted;
 * a separate live readout (updateLive) refreshes every frame. */

import type { World } from '../sim/world';
import type { AppState } from '../state';
import { selectedOf } from '../state';
import type { Scenario } from '../sim/scenario';
import { backgroundLabel } from '../sim/scenario';
import { attachTooltip } from './tooltip';

/* Brief descriptions shown as a delayed hover tooltip on each field label
 * (keyed by the exact label string used in the rows below).  Plain-language,
 * one line each -- the labels are terse so this is where the meaning lives. */
const TIPS: Record<string, string> = {
    /* quick toggles */
    'trails': "Draw each particle's past path as a fading trail.",
    'velocity': "Draw each particle's velocity as an arrow (length ∝ speed; brighter when relativistic).",
    'radiation': 'Include radiation-reaction (inductive / retarded) field terms. Off = quasi-static fields only.',
    /* global */
    'G_eff': 'Effective gravitational coupling — overall strength of gravity / GEM. 0 disables gravity.',
    'k_e': 'Effective electrostatic (Coulomb) coupling. 0 disables EM.',
    'force tier': 'newtonian = simple 1/r force; relativistic = full linearized-GR force law (post-Newtonian style).',
    'shape': 'Particle deposition shape on the grid: cic (linear), tsc (quadratic, smoother), bump (compact smooth).',
    'kernel R': 'Radius in cells of the deposit / interpolation kernel. Larger = smoother but slower and more diffuse.',
    'self-field': "Default for new particles: subtract a particle's own field so it feels only the others (kills PIC self-force).",
    'GEM force': "Gravitomagnetic force — the 'magnetic' part of gravity from mass currents (frame-dragging force).",
    'EM Lorentz': 'Apply the electromagnetic Lorentz force (qE + qv×B) to charged particles.',
    'perturbation': 'Evolve and source the dynamic (particle-generated) fields. Off = particles feel only the background body.',
    /* absorber */
    'outer BC': 'Outer-boundary condition for the field solve: dirichlet (value = 0) or neumann (zero normal gradient).',
    'friction floor': 'Baseline derivative-friction damping applied everywhere (gently bleeds off field energy).',
    'friction wall': 'Peak derivative-friction at the domain edge — the absorbing layer that soaks up outgoing waves.',
    'friction depth': 'Thickness in cells of the absorbing friction taper at the boundary.',
    'zero-mean φ': 'Subtract the mean of the scalar potentials each step (removes the unconstrained DC offset in a closed box).',
    /* background body */
    'present': 'Whether a fixed analytic background body (mass / charge / spin) exists at all.',
    'bg mode': 'How the background is applied: sampled onto the grid, or evaluated analytically per particle.',
    'mass GM': 'Background body gravitational mass parameter GM (Schwarzschild-like central mass).',
    'charge Q': 'Background body electric charge (Reissner–Nordström-like).',
    'ang.mom. Jz': 'Background body angular momentum about z — sources frame-dragging (Kerr / Lense–Thirring).',
    'gyromag. g': 'Gyromagnetic ratio of the background body (couples its spin to EM).',
    'softening ε': "Core softening length in cells — smooths the 1/r singularity at the body's center.",
    'x': 'Horizontal position in grid cells (origin at bottom-left).',
    'y': 'Vertical position in grid cells (origin at bottom-left).',
    /* advanced */
    'force interp': 'Force-interpolation scheme: legacy gather, or lewis-birdsall (momentum-conserving, no self-force).',
    'GM inductive': 'Include the inductive (time-derivative) part of the gravitomagnetic field.',
    'EM electrostatic': 'Include the electrostatic (Coulomb) part of the EM field.',
    'EM magnetic': 'Include the magnetostatic (current-sourced B) part of the EM field.',
    'EM inductive': 'Include the inductive (∂A/∂t) part of the EM field — needed for radiation.',
    'EM stress-energy': 'Let EM field energy gravitate (EM stress-energy sources the gravitational field).',
    'EM Shapiro': 'Slow EM propagation in a gravitational potential (Shapiro delay / gravitational lensing of light).',
    'Shapiro dynamic': 'Recompute the Shapiro light speed each step from the total dynamic mass field, so moving mass lenses EM.',
    'Esirkepov': 'Use charge-conserving Esirkepov current deposition (exact continuity for the EM source).',
    'periodic BC': 'Wrap the domain periodically instead of using the absorbing boundary.',
    'ρ smooth': 'Number of smoothing passes on the deposited mass density (reduces grid noise).',
    'J smooth': 'Number of smoothing passes on the deposited current density.',
    'field init': 'How the field is initialized: lw-settled (Liénard–Wiechert + settle), lw, none (auto), or zero.',
    'settle steps': 'Extra frozen-particle relaxation steps used to settle the initial field to its fixed point.',
    /* visualization */
    'tick Δ (0=auto)': 'Spacing between trajectory time-ticks (same Δ for coordinate & proper time; 0 = auto-choose).',
    'force: grav': 'Draw the gravitational / GEM force vector on each particle (orange).',
    'force: EM': 'Draw the EM Lorentz force vector on each particle (blue).',
    'force: spin': 'Draw the spin-gradient force vector on each particle (green).',
    'force: total': 'Draw the total applied force vector on each particle (white).',
    /* particle */
    'mass': 'Particle rest mass (> 0). Sets inertia and gravitational charge.',
    'charge': 'Particle electric charge (may be ±). Sources and feels EM.',
    'spin': 'Intrinsic spin angular momentum about z. Couples to field gradients (spin-gradient force) and precesses (Larmor / Thomas / Lense–Thirring).',
    'gFactor': 'Gyromagnetic ratio g of the particle’s spin. Sets the magnetic moment μ = g·(Q/2M)·spin. Default 2 (Dirac); 0 = no magnetic coupling.',
    'vx': 'Initial x-velocity in units of c (|v| must stay below cEff).',
    'vy': 'Initial y-velocity in units of c (|v| must stay below cEff).',
    'forces': 'Whether physical forces act on this particle. Off = pinned / pure source (a drive still applies).',
    'drive amp': 'Displacement amplitude of a sinusoidal drive (antenna element). 0 = no drive.',
    'drive ω': 'Angular frequency of the sinusoidal drive.',
    'drive phase': 'Phase offset of the drive (use π on one of a +q / −q pair for an antiphase dipole).',
    'drive axis x': 'X-component of the drive oscillation direction.',
    'drive axis y': 'Y-component of the drive oscillation direction.',
    'ticks': 'Trajectory time-tick markers: coordinate-time (circles), proper-time (triangles), both, or none.',
    'clock': 'Show a two-hand clock comparing coordinate vs proper time (visualizes time dilation).',
};

export interface InspectorHandlers {
    getScenario: () => Scenario;
    getState: () => AppState;
    getWorld: () => World;
    onEdit: (mutate: (s: Scenario) => void) => void;
    /* Visualization-only edit: mutate the scenario WITHOUT rebuilding the sim
     * (so the run keeps going) -- used for per-particle ticks/clock + tick Δ. */
    onViewEdit: (mutate: (s: Scenario) => void) => void;
    /* Remove the particle at `index` (rebuilds; undoable via the edit stack). */
    onDeleteParticle: (index: number) => void;
    /* Remove the background body at `index` (rebuilds; undoable). */
    onDeleteBody: (index: number) => void;
    /* Select a background body for editing (e.g. clicking the Global list). */
    onSelectBody: (index: number) => void;
}

function numRow(label: string, value: number, step: number,
                onCommit: (v: number) => void): HTMLElement {
    const row = document.createElement('div'); row.className = 'row';
    const l = document.createElement('label'); l.textContent = label;
    attachTooltip(l, TIPS[label]);
    const inp = document.createElement('input');
    inp.type = 'number'; inp.step = String(step); inp.value = String(value);
    inp.addEventListener('change', () => {
        const v = parseFloat(inp.value);
        if (Number.isFinite(v)) onCommit(v);
    });
    row.append(l, inp);
    return row;
}

function checkRow(label: string, checked: boolean,
                  onCommit: (v: boolean) => void): HTMLElement {
    const row = document.createElement('div'); row.className = 'row';
    const l = document.createElement('label'); l.textContent = label;
    attachTooltip(l, TIPS[label]);
    const inp = document.createElement('input');
    inp.type = 'checkbox'; inp.checked = checked;
    inp.addEventListener('change', () => onCommit(inp.checked));
    row.append(l, inp);
    return row;
}

function header(text: string): HTMLElement {
    const h = document.createElement('div');
    h.textContent = text;
    h.style.cssText = 'margin:8px 0 2px;color:#789;font-size:11px;'
        + 'text-transform:uppercase;letter-spacing:0.04em;';
    return h;
}

function selRow(label: string, value: string, opts: string[],
                onCommit: (v: string) => void): HTMLElement {
    const row = document.createElement('div'); row.className = 'row';
    const l = document.createElement('label'); l.textContent = label;
    attachTooltip(l, TIPS[label]);
    const sel = document.createElement('select');
    for (const o of opts) { const op = document.createElement('option'); op.value = o; op.textContent = o; sel.appendChild(op); }
    sel.value = value;
    sel.addEventListener('change', () => onCommit(sel.value));
    row.append(l, sel);
    return row;
}

export class Inspector {
    private statusEl: HTMLElement;
    private globalEl: HTMLElement;
    private particleEl: HTMLElement;
    private selectionTitle: HTMLElement | null;
    private statusLive: HTMLElement | null = null;
    private particleLive: HTMLElement | null = null;
    private h: InspectorHandlers;

    constructor(h: InspectorHandlers) {
        this.h = h;
        this.statusEl = document.getElementById('statusInfo') as HTMLElement;
        this.globalEl = document.getElementById('globalInfo') as HTMLElement;
        this.particleEl = document.getElementById('particleInfo') as HTMLElement;
        /* The Selection panel's heading -- retitled per what's selected. */
        this.selectionTitle = document.querySelector('#panelParticle h2');
    }

    renderEditors(): void {
        const s = this.h.getScenario();
        const edit = this.h.onEdit;
        const viewEdit = this.h.onViewEdit;

        /* ---- Status (read-only: title, live sim info, then description) ----
         * Live info goes first so the step/time stay visible by default; the
         * (potentially long) description sits below it. */
        this.statusEl.innerHTML = '';
        const title = document.createElement('div');
        title.textContent = s.name;
        title.style.cssText = 'font-weight:600;color:#cde;margin-bottom:4px;';
        this.statusEl.appendChild(title);
        const sl = document.createElement('div'); sl.className = 'live'; sl.textContent = '—';
        sl.style.marginTop = '0';
        this.statusEl.appendChild(sl); this.statusLive = sl;
        if (s.description) {
            const note = document.createElement('div');
            note.className = 'hint';
            note.textContent = s.description;
            this.statusEl.appendChild(note);
        }

        /* ---- Global (editable config) ---- */
        this.globalEl.innerHTML = '';

        /* Quick display/physics toggles up top (the common ones).  trails and
         * velocity are view-only (no rebuild); radiation flips global.noRadiation
         * and so rebuilds.  `radiation` is the user-facing inverse of noRadiation. */
        this.globalEl.appendChild(checkRow('trails', s.view.showTrails, (v) => viewEdit((sc) => { sc.view.showTrails = v; })));
        this.globalEl.appendChild(checkRow('velocity', s.view.showVelocity, (v) => viewEdit((sc) => { sc.view.showVelocity = v; })));
        this.globalEl.appendChild(checkRow('radiation', !s.global.noRadiation, (v) => edit((sc) => { sc.global.noRadiation = !v; })));

        this.globalEl.appendChild(numRow('G_eff', s.global.gEff, 0.1, (v) => edit((sc) => { sc.global.gEff = v; })));
        this.globalEl.appendChild(numRow('k_e', s.global.kE, 0.1, (v) => edit((sc) => { sc.global.kE = v; })));
        this.globalEl.appendChild(selRow('force tier', s.global.forceTier, ['newtonian', 'relativistic'], (v) => edit((sc) => { sc.global.forceTier = v as Scenario['global']['forceTier']; })));
        this.globalEl.appendChild(selRow('shape', s.global.shape, ['cic', 'tsc', 'bump'], (v) => edit((sc) => { sc.global.shape = v as Scenario['global']['shape']; })));
        this.globalEl.appendChild(numRow('kernel R', s.global.kernelRadius, 1, (v) => edit((sc) => { sc.global.kernelRadius = v; })));
        this.globalEl.appendChild(checkRow('self-field', s.global.selfFieldDefault, (v) => edit((sc) => { sc.global.selfFieldDefault = v; })));
        this.globalEl.appendChild(checkRow('GEM force', s.global.switches.gravitomagneticForce, (v) => edit((sc) => { sc.global.switches.gravitomagneticForce = v; })));
        this.globalEl.appendChild(checkRow('EM Lorentz', s.global.switches.emLorentz, (v) => edit((sc) => { sc.global.switches.emLorentz = v; })));

        /* Perturbation fields: one toggle for evolve + deposit.  Off => the
         * particles' own fields are neither sourced nor evolved, so a body
         * responds to ONLY the background (set one below).  Field initialization
         * follows automatically: applyScenario settles the field whenever
         * perturbation is active (see sim/build.ts), so this stays a dumb
         * switch and can't desync from the init method. */
        const pert = s.global.switches.fieldEvolution && s.global.switches.particleSourceDeposition;
        this.globalEl.appendChild(checkRow('perturbation', pert, (v) => edit((sc) => {
            sc.global.switches.fieldEvolution = v;
            sc.global.switches.particleSourceDeposition = v;
        })));

        /* ---- Absorbing boundary (v42 derivative-friction taper) ---- */
        const ab = s.global.absorber;
        this.globalEl.appendChild(header('absorber'));
        this.globalEl.appendChild(selRow('outer BC', ab.outerBC, ['dirichlet', 'neumann'],
            (v) => edit((sc) => { sc.global.absorber.outerBC = v as Scenario['global']['absorber']['outerBC']; })));
        this.globalEl.appendChild(numRow('friction floor', ab.frictionFloor, 0.001,
            (v) => edit((sc) => { sc.global.absorber.frictionFloor = v; })));
        this.globalEl.appendChild(numRow('friction wall', ab.frictionWall, 0.01,
            (v) => edit((sc) => { sc.global.absorber.frictionWall = v; })));
        this.globalEl.appendChild(numRow('friction depth', ab.frictionDepth, 1,
            (v) => edit((sc) => { sc.global.absorber.frictionDepth = v; })));
        this.globalEl.appendChild(checkRow('zero-mean φ', ab.zeroMeanScalar,
            (v) => edit((sc) => { sc.global.absorber.zeroMeanScalar = v; })));

        /* ---- Background bodies (a list; "+ body" is on the toolbar) ----
         * The per-body fields live in the Selection panel; here we show the
         * count, the global bg-eval mode, and a clickable list to select one. */
        const bodySelNow = selectedOf(this.h.getState(), 'body');
        this.globalEl.appendChild(header(`background bodies (${s.background.length})`));
        this.globalEl.appendChild(selRow('bg mode', s.global.bgMode, ['sampled', 'analytic'],
            (v) => edit((sc) => { sc.global.bgMode = v as Scenario['global']['bgMode']; })));
        const bgNote = document.createElement('div');
        bgNote.className = 'hint';
        bgNote.textContent = s.background.length === 0
            ? 'none — use "+ body" to add one.'
            : 'click a body below (or its circle) to edit it:';
        this.globalEl.appendChild(bgNote);
        s.background.forEach((b, i) => {
            const row = document.createElement('div');
            row.className = 'list-row' + (i === bodySelNow ? ' sel' : '');
            row.textContent = `▸ body ${i} — ${backgroundLabel(b)}`;
            row.title = 'select this background body to edit it';
            row.addEventListener('click', () => this.h.onSelectBody(i));
            this.globalEl.appendChild(row);
        });

        /* ---- Advanced switches (less commonly changed) ---- */
        this.globalEl.appendChild(header('advanced'));
        this.globalEl.appendChild(selRow('force interp', s.global.forceInterp, ['legacy', 'lewis-birdsall'],
            (v) => edit((sc) => { sc.global.forceInterp = v as Scenario['global']['forceInterp']; })));
        const swRow = (label: string, key: keyof Scenario['global']['switches']) =>
            this.globalEl.appendChild(checkRow(label, s.global.switches[key],
                (v) => edit((sc) => { sc.global.switches[key] = v; })));
        swRow('GM inductive', 'gravitomagneticInductive');
        swRow('EM electrostatic', 'emElectrostatic');
        swRow('EM magnetic', 'emMagnetic');
        swRow('EM inductive', 'emInductive');
        swRow('EM stress-energy', 'emStressEnergy');
        swRow('EM Shapiro', 'emShapiro');
        swRow('Shapiro dynamic', 'shapiroDynamic');
        swRow('Esirkepov', 'esirkepov');
        swRow('periodic BC', 'periodicBC');
        this.globalEl.appendChild(numRow('ρ smooth', s.global.rhoSmooth, 1, (v) => edit((sc) => { sc.global.rhoSmooth = v; })));
        this.globalEl.appendChild(numRow('J smooth', s.global.jSmooth, 1, (v) => edit((sc) => { sc.global.jSmooth = v; })));
        this.globalEl.appendChild(selRow('field init', s.global.init.method, ['lw-settled', 'lw', 'none', 'zero'],
            (v) => edit((sc) => { sc.global.init.method = v as Scenario['global']['init']['method']; })));
        this.globalEl.appendChild(numRow('settle steps', s.global.init.settleSteps, 1, (v) => edit((sc) => { sc.global.init.settleSteps = v; })));

        /* ---- Visualization (view-only; no rebuild) ---- */
        this.globalEl.appendChild(header('visualization'));
        this.globalEl.appendChild(numRow('tick Δ (0=auto)', s.view.tickInterval, 0.5,
            (v) => viewEdit((sc) => { sc.view.tickInterval = Math.max(0, v); })));

        /* Per-particle force arrows (drawn on every particle). */
        const fa = s.view.forceArrows;
        const faRow = (label: string, key: keyof typeof fa) =>
            this.globalEl.appendChild(checkRow(label, fa[key],
                (v) => viewEdit((sc) => { sc.view.forceArrows[key] = v; })));
        faRow('force: grav', 'grav');
        faRow('force: EM', 'em');
        faRow('force: spin', 'spin');
        faRow('force: total', 'total');

        /* ---- Selection panel: the selected particle OR background body ---- */
        this.particleEl.innerHTML = '';
        this.particleLive = null;
        const st = this.h.getState();
        const pSel = selectedOf(st, 'particle');
        const bSel = selectedOf(st, 'body');
        if (pSel >= 0 && pSel < s.particles.length) {
            if (this.selectionTitle) this.selectionTitle.textContent = `Particle ${pSel}`;
            this.renderParticleEditor(s, pSel, edit, viewEdit);
        } else if (bSel >= 0 && bSel < s.background.length) {
            if (this.selectionTitle) this.selectionTitle.textContent = `Background body ${bSel}`;
            this.renderBodyEditor(s, bSel, edit);
        } else {
            if (this.selectionTitle) this.selectionTitle.textContent = 'Selection';
            this.particleEl.textContent = 'click a particle or a background body';
        }
    }

    /* Editable fields for the selected particle (+ live readout, drive, ticks). */
    private renderParticleEditor(s: Scenario, sel: number,
                                 edit: (m: (s: Scenario) => void) => void,
                                 viewEdit: (m: (s: Scenario) => void) => void): void {
        /* Live readout up top (mirrors the Status panel), so the running values
         * stay visible above the editable fields. */
        const pl = document.createElement('div'); pl.className = 'live'; pl.textContent = '—';
        pl.style.marginTop = '0';
        this.particleEl.appendChild(pl); this.particleLive = pl;

        const mk = (label: string, key: 'x' | 'y' | 'vx' | 'vy' | 'mass' | 'charge', step: number) =>
            this.particleEl.appendChild(numRow(label, s.particles[sel][key], step,
                (v) => edit((sc) => { sc.particles[sel][key] = v; })));
        mk('mass', 'mass', 0.001);
        mk('charge', 'charge', 0.001);
        mk('x', 'x', 1);
        mk('y', 'y', 1);
        mk('vx', 'vx', 0.01);
        mk('vy', 'vy', 0.01);

        /* Intrinsic spin (optional). spin sets the angular momentum that feels
         * the spin-gradient force and precesses; gFactor sets its magnetic
         * moment μ = g·(Q/2M)·spin (default 2, Dirac).  Physics edit => rebuild. */
        this.particleEl.appendChild(numRow('spin', s.particles[sel].spin ?? 0, 0.1,
            (v) => edit((sc) => { sc.particles[sel].spin = v; })));
        this.particleEl.appendChild(numRow('gFactor', s.particles[sel].gFactor ?? 2, 0.1,
            (v) => edit((sc) => { sc.particles[sel].gFactor = v; })));

        /* Source / drive: per-particle sinusoidal oscillation (one antenna
         * element).  `forces` off => pinned / pure source; `amp` 0 => no drive.
         * Displacement amplitude `amp`, angular freq `ω`, `phase` (use π between
         * a +q/-q pair), and drive `axis`.  Physics edits => rebuild. */
        this.particleEl.appendChild(header('source / drive'));
        this.particleEl.appendChild(checkRow('forces', s.particles[sel].forces !== false,
            (v) => edit((sc) => { sc.particles[sel].forces = v; })));
        const drv = s.particles[sel].drive;
        const driveEdit = (mut: (d: NonNullable<Scenario['particles'][number]['drive']>) => void) =>
            edit((sc) => {
                const p = sc.particles[sel];
                if (!p.drive) p.drive = { amp: 0, omega: 0, phase: 0, axis: [0, 1] };
                mut(p.drive);
            });
        this.particleEl.appendChild(numRow('drive amp', drv?.amp ?? 0, 0.1,
            (v) => driveEdit((d) => { d.amp = v; })));
        if ((drv?.amp ?? 0) !== 0 || (drv?.omega ?? 0) !== 0) {
            this.particleEl.appendChild(numRow('drive ω', drv?.omega ?? 0, 0.01,
                (v) => driveEdit((d) => { d.omega = v; })));
            this.particleEl.appendChild(numRow('drive phase', drv?.phase ?? 0, 0.1,
                (v) => driveEdit((d) => { d.phase = v; })));
            this.particleEl.appendChild(numRow('drive axis x', drv?.axis?.[0] ?? 0, 0.1,
                (v) => driveEdit((d) => { d.axis = [v, d.axis[1]]; })));
            this.particleEl.appendChild(numRow('drive axis y', drv?.axis?.[1] ?? 1, 0.1,
                (v) => driveEdit((d) => { d.axis = [d.axis[0], v]; })));
        }

        /* Visualization (view-only): trajectory ticks + proper-time clock. */
        this.particleEl.appendChild(header('time visualization'));
        this.particleEl.appendChild(selRow('ticks', s.particles[sel].ticks ?? 'none',
            ['none', 'coordinate', 'proper', 'both'],
            (v) => viewEdit((sc) => { sc.particles[sel].ticks = v as Scenario['particles'][number]['ticks']; })));
        this.particleEl.appendChild(checkRow('clock', !!s.particles[sel].clock,
            (v) => viewEdit((sc) => { sc.particles[sel].clock = v; })));

        /* Remove this particle (rebuilds; undoable). */
        const del = document.createElement('button');
        del.textContent = 'delete particle';
        del.className = 'del-btn';
        del.title = 'remove this particle (undo restores it)';
        del.addEventListener('click', () => this.h.onDeleteParticle(sel));
        this.particleEl.appendChild(del);
    }

    /* Editable fields for the selected background body (fixed compact source). */
    private renderBodyEditor(s: Scenario, sel: number,
                             edit: (m: (s: Scenario) => void) => void): void {
        const b = s.background[sel];
        const mk = (label: string, key: 'GM' | 'Q' | 'Jz' | 'gFactor' | 'epsilon' | 'x' | 'y', step: number) =>
            this.particleEl.appendChild(numRow(label, b[key], step,
                (v) => edit((sc) => { sc.background[sel][key] = v; })));
        mk('mass GM', 'GM', 0.001);
        mk('charge Q', 'Q', 0.001);
        mk('ang.mom. Jz', 'Jz', 0.1);
        mk('gyromag. g', 'gFactor', 0.1);
        mk('softening ε', 'epsilon', 0.5);
        mk('x', 'x', 1);
        mk('y', 'y', 1);

        const del = document.createElement('button');
        del.textContent = 'delete body';
        del.className = 'del-btn';
        del.title = 'remove this background body (undo restores it)';
        del.addEventListener('click', () => this.h.onDeleteBody(sel));
        this.particleEl.appendChild(del);
    }

    updateLive(): void {
        const w = this.h.getWorld();
        const st = this.h.getState();
        if (this.statusLive) {
            const ps = w.particles();
            let Px = 0, Py = 0, Mx = 0, My = 0, M = 0;
            for (const p of ps) { Px += p.px; Py += p.py; Mx += p.mass * p.x; My += p.mass * p.y; M += p.mass; }
            const com = M > 0 ? `(${(Mx / M).toFixed(2)}, ${(My / M).toFixed(2)})` : '—';
            const sep = ps.length === 2 ? Math.hypot(ps[0].x - ps[1].x, ps[0].y - ps[1].y).toFixed(3) : '—';
            this.statusLive.textContent =
                `step ${w.stepCount()}  t=${w.time().toFixed(2)}  ${st.paused ? 'paused' : 'running'}` +
                `${st.liveModified ? '  [live-modified]' : ''}\n` +
                `P=(${Px.toExponential(2)},${Py.toExponential(2)})  COM=${com}  sep=${sep}`;
        }
        const pSel = selectedOf(st, 'particle');
        if (this.particleLive && pSel >= 0 && pSel < w.particleCount()) {
            const p = w.particle(pSel);
            const g = p.mass > 0 ? Math.sqrt(1 + (p.px * p.px + p.py * p.py) / (p.mass * p.mass)) : 1;
            const v = p.mass > 0 ? Math.hypot(p.px, p.py) / (g * p.mass) : 0;
            this.particleLive.textContent =
                `live: pos=(${p.x.toFixed(2)},${p.y.toFixed(2)})  |v|=${v.toExponential(3)}\n` +
                `p=(${p.px.toExponential(2)},${p.py.toExponential(2)})  τ=${p.properTime.toFixed(2)}`;
        }
    }
}
