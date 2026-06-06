/* Read-only inspector panels (Phase 1c).  Phase 3 makes these editable.
 * Renders the global sim state and the selected particle's values as text. */

import type { World } from '../sim/world';
import type { AppState } from '../state';
import { FIELD_VIEWS } from '../sim/fieldViews';

export class Inspector {
    private globalEl: HTMLElement;
    private particleEl: HTMLElement;

    constructor() {
        this.globalEl = document.getElementById('globalInfo') as HTMLElement;
        this.particleEl = document.getElementById('particleInfo') as HTMLElement;
    }

    update(world: World, state: AppState): void {
        let Px = 0, Py = 0, Mx = 0, My = 0, M = 0;
        const ps = world.particles();
        for (const p of ps) { Px += p.px; Py += p.py; Mx += p.mass * p.x; My += p.mass * p.y; M += p.mass; }
        const com = M > 0 ? `(${(Mx / M).toFixed(2)}, ${(My / M).toFixed(2)})` : '—';
        const sep = ps.length === 2
            ? Math.hypot(ps[0].x - ps[1].x, ps[0].y - ps[1].y).toFixed(3) : '—';

        this.globalEl.textContent = [
            `scenario : ${state.scenario}${state.liveModified ? '  [live-modified]' : ''}`,
            `step/t   : ${world.stepCount()}  /  ${world.time().toFixed(2)}`,
            `grid     : ${world.W}×${world.H}   dt=${world.dt().toFixed(4)}`,
            `field    : ${FIELD_VIEWS[state.viewField]?.label ?? state.viewField}`,
            `particles: ${ps.length}`,
            `P_total  : (${Px.toExponential(2)}, ${Py.toExponential(2)})`,
            `COM      : ${com}    sep: ${sep}`,
            `run      : ${state.paused ? 'paused' : 'playing'}`,
        ].join('\n');

        if (state.selected < 0 || state.selected >= ps.length) {
            this.particleEl.textContent = 'click a particle to select';
            return;
        }
        const p = ps[state.selected];
        const pmag2 = p.px * p.px + p.py * p.py;
        const gamma = p.mass > 0 ? Math.sqrt(1 + pmag2 / (p.mass * p.mass)) : 1;
        const v = p.mass > 0 ? Math.hypot(p.px, p.py) / (gamma * p.mass) : 0;
        this.particleEl.textContent = [
            `index    : ${p.index}`,
            `position : (${p.x.toFixed(3)}, ${p.y.toFixed(3)})`,
            `momentum : (${p.px.toExponential(3)}, ${p.py.toExponential(3)})`,
            `mass     : ${p.mass}`,
            `charge   : ${p.charge}`,
            `|v|      : ${v.toExponential(3)}   γ=${gamma.toFixed(4)}`,
            `spin     : ${p.spin}   g=${p.gFactor}`,
            `proper τ : ${p.properTime.toFixed(3)}`,
        ].join('\n');
    }
}
