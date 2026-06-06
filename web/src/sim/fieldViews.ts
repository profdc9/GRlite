/* Selectable field views.  Scalar potentials are displayed directly; derived
 * quantities (curl for B_z, gradient magnitude for |g|/|E|) are computed on
 * the CPU from the potential grids -- cheap (one W*H pass for the displayed
 * field) and keeps the colormap auto-normalization (which needs the value
 * range) correct without a multi-texture shader.
 *
 * Each visualization reads a chosen SOURCE of the field: the perturbation
 * (the evolving grid, gr_sim_field_ptr), the background (analytic/sampled
 * grid, gr_sim_background_ptr), their sum, or none.  This lets a scene show
 * e.g. the spinning-mass background A_g that the perturbation grid never
 * carries.
 *
 * Note: derived fields use texel-unit finite differences (the physical 1/dx
 * scale is irrelevant since the colormap auto-normalizes). */

import type { World } from './world';
import type { FieldSourceName } from './scenario';
import {
    GR_FIELD_PHI_GRAV, GR_FIELD_A_GX, GR_FIELD_A_GY,
    GR_FIELD_PHI_EM, GR_FIELD_A_X, GR_FIELD_A_Y,
} from './config';

type ViewMode = 'scalar' | 'curl' | 'gradmag' | 'none';

export interface FieldView {
    label: string;
    mode: ViewMode;
    a: number;       // primary field id
    b?: number;      // secondary (for curl: a=Ax, b=Ay)
}

export const FIELD_VIEWS: FieldView[] = [
    { label: 'Φ_g  (grav potential)', mode: 'scalar', a: GR_FIELD_PHI_GRAV },
    { label: 'B_g,z  (frame drag, curl A_g)', mode: 'curl', a: GR_FIELD_A_GX, b: GR_FIELD_A_GY },
    { label: '|g|  (grav field, |∇Φ_g|)', mode: 'gradmag', a: GR_FIELD_PHI_GRAV },
    { label: 'A_g,x  (grav vector pot.)', mode: 'scalar', a: GR_FIELD_A_GX },
    { label: 'A_g,y  (grav vector pot.)', mode: 'scalar', a: GR_FIELD_A_GY },
    { label: 'φ_em  (electric potential)', mode: 'scalar', a: GR_FIELD_PHI_EM },
    { label: 'B_z  (magnetic, curl A)', mode: 'curl', a: GR_FIELD_A_X, b: GR_FIELD_A_Y },
    { label: '— none (black) —', mode: 'none', a: -1 },
];

/* Whether the selected view draws a colormap at all ('none' => black). */
export function colormapEnabled(viewIndex: number): boolean {
    return (FIELD_VIEWS[viewIndex] ?? FIELD_VIEWS[0]).mode !== 'none';
}

/* 2D vector fields drawable as an arrow overlay.  'components' takes (a,b)
 * directly as (Vx,Vy); 'gradient' forms V = -∇(a) from a scalar.  Index 0 is
 * the off entry. */
export interface VectorField {
    label: string;
    kind: 'components' | 'gradient';
    a: number;       // < 0 means "none"
    b?: number;
}

export const VECTOR_FIELDS: VectorField[] = [
    { label: 'none', kind: 'components', a: -1 },
    { label: 'A_g  (gravitomagnetic)', kind: 'components', a: GR_FIELD_A_GX, b: GR_FIELD_A_GY },
    { label: 'g = −∇Φ_g  (gravity)', kind: 'gradient', a: GR_FIELD_PHI_GRAV },
    { label: 'A  (magnetic vec. pot.)', kind: 'components', a: GR_FIELD_A_X, b: GR_FIELD_A_Y },
    { label: 'E = −∇φ_em  (electric)', kind: 'gradient', a: GR_FIELD_PHI_EM },
];

export interface VectorData { vx: Float32Array; vy: Float32Array; }

/* ---- scratch pools (reused across frames; never escape past the next call) ---- */
let outBuf: Float32Array | null = null;          // computeView curl/gradmag output
function out(n: number): Float32Array {
    if (!outBuf || outBuf.length !== n) outBuf = new Float32Array(n);
    return outBuf;
}
const sumPool: Float32Array[] = [];              // per-slot buffers for 'sum'
function sumBuf(slot: number, n: number): Float32Array {
    if (!sumPool[slot] || sumPool[slot].length !== n) sumPool[slot] = new Float32Array(n);
    return sumPool[slot];
}
const gradPool: Float32Array[] = [];             // gradient-vector components
function gradBuf(slot: number, n: number): Float32Array {
    if (!gradPool[slot] || gradPool[slot].length !== n) gradPool[slot] = new Float32Array(n);
    return gradPool[slot];
}
let zeroBuf: Float32Array | null = null;         // read-only zeros (missing bg)
function zeros(n: number): Float32Array {
    if (!zeroBuf || zeroBuf.length !== n) zeroBuf = new Float32Array(n);
    return zeroBuf;
}

/* Resolve a field id to a concrete array for the chosen source.  `slot`
 * selects a distinct sum-scratch buffer so callers needing two fields (curl,
 * vector components) don't alias.  'none' is treated as perturbation here;
 * callers that mean "draw nothing" check the source before calling. */
function sourced(world: World, id: number, source: FieldSourceName, slot: number): Float32Array {
    const n = world.W * world.H;
    if (source === 'background') return world.backgroundView(id) ?? zeros(n);
    if (source === 'sum') {
        const p = world.fieldView(id);
        const bg = world.backgroundView(id);
        if (!bg) return p;
        const dst = sumBuf(slot, n);
        for (let k = 0; k < n; k++) dst[k] = p[k] + bg[k];
        return dst;
    }
    return world.fieldView(id);
}

/* Returns a Float32Array for the selected colormap view + source.  For scalar
 * views with a single-array source this may be a live heap view (no copy);
 * derived views and 'sum' use reused scratch.  Valid until the next
 * computeView / computeVectorField / sim step. */
export function computeView(world: World, viewIndex: number,
                           source: FieldSourceName = 'perturbation'): Float32Array {
    const view = FIELD_VIEWS[viewIndex] ?? FIELD_VIEWS[0];
    const W = world.W, H = world.H;
    if (view.mode === 'none') return zeros(W * H);
    if (view.mode === 'scalar') return sourced(world, view.a, source, 0);

    const dst = out(W * H);
    const at = (f: Float32Array, i: number, j: number) =>
        f[Math.min(H - 1, Math.max(0, j)) * W + Math.min(W - 1, Math.max(0, i))];

    if (view.mode === 'curl') {
        const Ax = sourced(world, view.a, source, 0);
        const Ay = sourced(world, view.b!, source, 1);
        for (let j = 0; j < H; j++) {
            for (let i = 0; i < W; i++) {
                const dAydx = (at(Ay, i + 1, j) - at(Ay, i - 1, j)) * 0.5;
                const dAxdy = (at(Ax, i, j + 1) - at(Ax, i, j - 1)) * 0.5;
                dst[j * W + i] = dAydx - dAxdy;
            }
        }
    } else { // gradmag
        const P = sourced(world, view.a, source, 0);
        for (let j = 0; j < H; j++) {
            for (let i = 0; i < W; i++) {
                const gx = (at(P, i + 1, j) - at(P, i - 1, j)) * 0.5;
                const gy = (at(P, i, j + 1) - at(P, i, j - 1)) * 0.5;
                dst[j * W + i] = Math.hypot(gx, gy);
            }
        }
    }
    return dst;
}

/* Returns (Vx,Vy) full-grid arrays for the selected vector field + source, or
 * null for the 'none' entry.  The overlay samples these on its own grid and
 * normalizes arrow lengths.  Uses scratch distinct from computeView's. */
export function computeVectorField(world: World, vIndex: number,
                                   source: FieldSourceName): VectorData | null {
    const v = VECTOR_FIELDS[vIndex];
    if (!v || v.a < 0) return null;
    const W = world.W, H = world.H, n = W * H;

    if (v.kind === 'components') {
        return { vx: sourced(world, v.a, source, 2), vy: sourced(world, v.b!, source, 3) };
    }
    /* gradient: V = -∇P */
    const P = sourced(world, v.a, source, 2);
    const at = (i: number, j: number) =>
        P[Math.min(H - 1, Math.max(0, j)) * W + Math.min(W - 1, Math.max(0, i))];
    const vx = gradBuf(0, n), vy = gradBuf(1, n);
    for (let j = 0; j < H; j++) {
        for (let i = 0; i < W; i++) {
            vx[j * W + i] = -(at(i + 1, j) - at(i - 1, j)) * 0.5;
            vy[j * W + i] = -(at(i, j + 1) - at(i, j - 1)) * 0.5;
        }
    }
    return { vx, vy };
}
