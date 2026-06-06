/* Developer diagnostics: text reports computed from the World state.
 *
 * These do the "hard work" in-app and return digestible strings (the same
 * style the debug bridge will serve over the wire).  No rendering. */

import type { World } from '../sim/world';
import { GR_FIELD_PHI_GRAV } from '../sim/config';

export function fieldStats(world: World, which = GR_FIELD_PHI_GRAV): string {
    const f = world.fieldView(which);
    const W = world.W, H = world.H;
    let min = Infinity, max = -Infinity, sum = 0, sumSq = 0;
    let iMin = 0, iMax = 0;
    for (let k = 0; k < f.length; k++) {
        const v = f[k];
        if (v < min) { min = v; iMin = k; }
        if (v > max) { max = v; iMax = k; }
        sum += v; sumSq += v * v;
    }
    const n = f.length;
    const mean = sum / n;
    const rms = Math.sqrt(sumSq / n);
    const xy = (k: number) => `(${k % W},${Math.floor(k / W)})`;
    return [
        `field ${which}: min=${min.toExponential(3)} @${xy(iMin)}`,
        `max=${max.toExponential(3)} @${xy(iMax)}`,
        `mean=${mean.toExponential(3)} rms=${rms.toExponential(3)} (${W}x${H})`,
    ].join('  ');
}

/* 1D slice of a field along a row (y=const) or column (x=const). */
export function fieldProfile(world: World, which: number,
                            axis: 'row' | 'col', at: number, stride = 16): string {
    const f = world.fieldView(which);
    const W = world.W, H = world.H;
    const out: string[] = [];
    if (axis === 'row') {
        for (let i = 0; i < W; i += stride) out.push(`(${i},${f[at * W + i].toExponential(2)})`);
    } else {
        for (let j = 0; j < H; j += stride) out.push(`(${j},${f[j * W + at].toExponential(2)})`);
    }
    return out.join(' ');
}

export function particleReport(world: World): string {
    const ps = world.particles();
    const lines = ps.map((p) => {
        const v = Math.hypot(p.px, p.py) / (p.mass || 1);
        return `  p${p.index}: x=(${p.x.toFixed(3)},${p.y.toFixed(3)}) ` +
               `p=(${p.px.toExponential(3)},${p.py.toExponential(3)}) ` +
               `m=${p.mass} q=${p.charge} |v|=${v.toExponential(3)} tau=${p.properTime.toFixed(3)}`;
    });
    return lines.join('\n');
}

/* Aggregate conserved-ish quantities for a symmetric system. */
export function conservationReport(world: World): string {
    const ps = world.particles();
    let Px = 0, Py = 0, Mx = 0, My = 0, M = 0;
    for (const p of ps) { Px += p.px; Py += p.py; Mx += p.mass * p.x; My += p.mass * p.y; M += p.mass; }
    const com = M > 0 ? `(${(Mx / M).toFixed(3)},${(My / M).toFixed(3)})` : 'n/a';
    let sep = 'n/a';
    if (ps.length === 2) {
        sep = Math.hypot(ps[0].x - ps[1].x, ps[0].y - ps[1].y).toFixed(3);
    }
    return `P_total=(${Px.toExponential(3)},${Py.toExponential(3)})  COM=${com}  sep=${sep}`;
}

export function fullReport(world: World): string {
    return [
        `step=${world.stepCount()} t=${world.time().toFixed(3)} dt=${world.dt().toFixed(4)}`,
        conservationReport(world),
        particleReport(world),
        fieldStats(world, GR_FIELD_PHI_GRAV),
    ].join('\n');
}

/* Frozen-field stability probe (no particle motion): does Phi_g hold? */
export function stabilityProbe(world: World, batches: number, perBatch: number): string {
    world.setParticlesFrozen(true);
    const cx = world.W >> 1, cy = world.H >> 1;
    const series: string[] = [];
    let pmin = Infinity, pmax = -Infinity;
    for (let b = 0; b < batches; b++) {
        const phi = world.fieldView(GR_FIELD_PHI_GRAV);
        const center = phi[cy * world.W + cx];
        let lmin = Infinity, lmax = -Infinity;
        for (const v of phi) { if (v < lmin) lmin = v; if (v > lmax) lmax = v; }
        if (lmin < pmin) pmin = lmin;
        if (lmax > pmax) pmax = lmax;
        series.push(`+${b * perBatch}: center=${center.toExponential(3)} [${lmin.toExponential(2)},${lmax.toExponential(2)}]`);
        world.step(perBatch);
    }
    world.setParticlesFrozen(false);
    return `stability probe ${batches}x${perBatch}: swing=${(pmax - pmin).toExponential(3)}\n  ${series.join('\n  ')}`;
}
