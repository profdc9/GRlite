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

export interface InspectorHandlers {
    getScenario: () => Scenario;
    getState: () => AppState;
    getWorld: () => World;
    onEdit: (mutate: (s: Scenario) => void) => void;
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
    private globalEl: HTMLElement;
    private particleEl: HTMLElement;
    private globalLive: HTMLElement | null = null;
    private particleLive: HTMLElement | null = null;
    private h: InspectorHandlers;

    constructor(h: InspectorHandlers) {
        this.h = h;
        this.globalEl = document.getElementById('globalInfo') as HTMLElement;
        this.particleEl = document.getElementById('particleInfo') as HTMLElement;
    }

    renderEditors(): void {
        const s = this.h.getScenario();
        const edit = this.h.onEdit;

        /* ---- Global ---- */
        this.globalEl.innerHTML = '';
        this.globalEl.appendChild(numRow('G_eff', s.global.gEff, 0.1, (v) => edit((sc) => { sc.global.gEff = v; })));
        this.globalEl.appendChild(numRow('k_e', s.global.kE, 0.1, (v) => edit((sc) => { sc.global.kE = v; })));
        this.globalEl.appendChild(selRow('shape', s.global.shape, ['cic', 'tsc', 'bump'], (v) => edit((sc) => { sc.global.shape = v as Scenario['global']['shape']; })));
        this.globalEl.appendChild(numRow('kernel R', s.global.kernelRadius, 1, (v) => edit((sc) => { sc.global.kernelRadius = v; })));
        this.globalEl.appendChild(checkRow('self-field', s.global.selfFieldDefault, (v) => edit((sc) => { sc.global.selfFieldDefault = v; })));
        this.globalEl.appendChild(checkRow('no-radiation', s.global.noRadiation, (v) => edit((sc) => { sc.global.noRadiation = v; })));
        this.globalEl.appendChild(checkRow('GEM force', s.global.switches.gravitomagneticForce, (v) => edit((sc) => { sc.global.switches.gravitomagneticForce = v; })));
        this.globalEl.appendChild(checkRow('EM Lorentz', s.global.switches.emLorentz, (v) => edit((sc) => { sc.global.switches.emLorentz = v; })));
        const gl = document.createElement('div'); gl.className = 'live'; gl.textContent = '—';
        this.globalEl.appendChild(gl); this.globalLive = gl;

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
        const pl = document.createElement('div'); pl.className = 'live'; pl.textContent = '—';
        this.particleEl.appendChild(pl); this.particleLive = pl;
    }

    updateLive(): void {
        const w = this.h.getWorld();
        const st = this.h.getState();
        if (this.globalLive) {
            const ps = w.particles();
            let Px = 0, Py = 0, Mx = 0, My = 0, M = 0;
            for (const p of ps) { Px += p.px; Py += p.py; Mx += p.mass * p.x; My += p.mass * p.y; M += p.mass; }
            const com = M > 0 ? `(${(Mx / M).toFixed(2)}, ${(My / M).toFixed(2)})` : '—';
            const sep = ps.length === 2 ? Math.hypot(ps[0].x - ps[1].x, ps[0].y - ps[1].y).toFixed(3) : '—';
            this.globalLive.textContent =
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
