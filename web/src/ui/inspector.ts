/* Editable inspectors (Phase 3).
 *
 * The panels edit the canonical Scenario; committing an edit goes through
 * onEdit(), which (in main) snapshots for undo and rebuilds from the JSON
 * (the deterministic edit->rebuild model).  Editable inputs are rebuilt only
 * on scenario/selection change (renderEditors), so typing isn't interrupted;
 * a separate live readout (updateLive) refreshes every frame. */

import type { World } from '../sim/world';
import type { AppState } from '../state';
import type { Scenario } from '../sim/scenario';
import { backgroundLabel, defaultBackground } from '../sim/scenario';

export interface InspectorHandlers {
    getScenario: () => Scenario;
    getState: () => AppState;
    getWorld: () => World;
    onEdit: (mutate: (s: Scenario) => void) => void;
    /* Visualization-only edit: mutate the scenario WITHOUT rebuilding the sim
     * (so the run keeps going) -- used for per-particle ticks/clock + tick Δ. */
    onViewEdit: (mutate: (s: Scenario) => void) => void;
}

function numRow(label: string, value: number, step: number,
                onCommit: (v: number) => void): HTMLElement {
    const row = document.createElement('div'); row.className = 'row';
    const l = document.createElement('label'); l.textContent = label;
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
    const inp = document.createElement('input');
    inp.type = 'checkbox'; inp.checked = checked;
    inp.addEventListener('change', () => onCommit(inp.checked));
    row.append(l, inp);
    return row;
}

/* Add or remove the single (unified) background body.  When adding, seed a
 * Schwarzschild-ish body at the grid center (mass only; user dials in Q/Jz). */
function setBackgroundPresent(sc: Scenario, present: boolean): void {
    if (!present) { sc.background = []; return; }
    if (sc.background[0]) return;
    /* Seed from the canonical default body (mass 0.01 so it's visible). */
    sc.background = [{ ...defaultBackground(sc.grid.W, sc.grid.H), GM: 0.01 }];
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
    private statusLive: HTMLElement | null = null;
    private particleLive: HTMLElement | null = null;
    private h: InspectorHandlers;

    constructor(h: InspectorHandlers) {
        this.h = h;
        this.statusEl = document.getElementById('statusInfo') as HTMLElement;
        this.globalEl = document.getElementById('globalInfo') as HTMLElement;
        this.particleEl = document.getElementById('particleInfo') as HTMLElement;
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
        this.globalEl.appendChild(numRow('G_eff', s.global.gEff, 0.1, (v) => edit((sc) => { sc.global.gEff = v; })));
        this.globalEl.appendChild(numRow('k_e', s.global.kE, 0.1, (v) => edit((sc) => { sc.global.kE = v; })));
        this.globalEl.appendChild(selRow('force tier', s.global.forceTier, ['newtonian', 'relativistic'], (v) => edit((sc) => { sc.global.forceTier = v as Scenario['global']['forceTier']; })));
        this.globalEl.appendChild(selRow('shape', s.global.shape, ['cic', 'tsc', 'bump'], (v) => edit((sc) => { sc.global.shape = v as Scenario['global']['shape']; })));
        this.globalEl.appendChild(numRow('kernel R', s.global.kernelRadius, 1, (v) => edit((sc) => { sc.global.kernelRadius = v; })));
        this.globalEl.appendChild(checkRow('self-field', s.global.selfFieldDefault, (v) => edit((sc) => { sc.global.selfFieldDefault = v; })));
        this.globalEl.appendChild(checkRow('no-radiation', s.global.noRadiation, (v) => edit((sc) => { sc.global.noRadiation = v; })));
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

        /* ---- Background body (unified M, Q, Jz -- see backgroundLabel) ---- */
        const bg = s.background[0];
        this.globalEl.appendChild(header(`background — ${backgroundLabel(bg)}`));
        this.globalEl.appendChild(checkRow('present', !!bg,
            (v) => edit((sc) => setBackgroundPresent(sc, v))));
        if (bg) {
            this.globalEl.appendChild(selRow('bg mode', s.global.bgMode, ['sampled', 'analytic'],
                (v) => edit((sc) => { sc.global.bgMode = v as Scenario['global']['bgMode']; })));
            this.globalEl.appendChild(numRow('mass GM', bg.GM, 0.001, (v) => edit((sc) => { sc.background[0].GM = v; })));
            this.globalEl.appendChild(numRow('charge Q', bg.Q, 0.001, (v) => edit((sc) => { sc.background[0].Q = v; })));
            this.globalEl.appendChild(numRow('ang.mom. Jz', bg.Jz, 0.1, (v) => edit((sc) => { sc.background[0].Jz = v; })));
            this.globalEl.appendChild(numRow('gyromag. g', bg.gFactor, 0.1, (v) => edit((sc) => { sc.background[0].gFactor = v; })));
            this.globalEl.appendChild(numRow('softening ε', bg.epsilon, 0.5, (v) => edit((sc) => { sc.background[0].epsilon = v; })));
            this.globalEl.appendChild(numRow('x', bg.x, 1, (v) => edit((sc) => { sc.background[0].x = v; })));
            this.globalEl.appendChild(numRow('y', bg.y, 1, (v) => edit((sc) => { sc.background[0].y = v; })));
        }

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

        /* ---- Particle ---- */
        this.particleEl.innerHTML = '';
        const sel = this.h.getState().selected;
        if (sel < 0 || sel >= s.particles.length) {
            this.particleEl.textContent = 'click a particle to select';
            this.particleLive = null;
            return;
        }
        const mk = (label: string, key: 'x' | 'y' | 'vx' | 'vy' | 'mass' | 'charge', step: number) =>
            this.particleEl.appendChild(numRow(label, s.particles[sel][key], step,
                (v) => edit((sc) => { sc.particles[sel][key] = v; })));
        mk('mass', 'mass', 0.001);
        mk('charge', 'charge', 0.001);
        mk('x', 'x', 1);
        mk('y', 'y', 1);
        mk('vx', 'vx', 0.01);
        mk('vy', 'vy', 0.01);

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

        const pl = document.createElement('div'); pl.className = 'live'; pl.textContent = '—';
        this.particleEl.appendChild(pl); this.particleLive = pl;
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
        if (this.particleLive && st.selected >= 0 && st.selected < w.particleCount()) {
            const p = w.particle(st.selected);
            const g = p.mass > 0 ? Math.sqrt(1 + (p.px * p.px + p.py * p.py) / (p.mass * p.mass)) : 1;
            const v = p.mass > 0 ? Math.hypot(p.px, p.py) / (g * p.mass) : 0;
            this.particleLive.textContent =
                `live: pos=(${p.x.toFixed(2)},${p.y.toFixed(2)})  |v|=${v.toExponential(3)}\n` +
                `p=(${p.px.toExponential(2)},${p.py.toExponential(2)})  τ=${p.properTime.toFixed(2)}`;
        }
    }
}
