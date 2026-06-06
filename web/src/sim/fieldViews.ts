/* Selectable field views.  Scalar potentials are displayed directly; derived
 * quantities (curl for B_z, gradient magnitude for |g|/|E|) are computed on
 * the CPU from the potential grids -- cheap (one W*H pass for the displayed
 * field) and keeps the colormap auto-normalization (which needs the value
 * range) correct without a multi-texture shader.
 *
 * Note: derived fields use texel-unit finite differences (the physical 1/dx
 * scale is irrelevant since the colormap auto-normalizes). */

import type { World } from './world';
import {
    GR_FIELD_PHI_GRAV, GR_FIELD_A_GX, GR_FIELD_A_GY,
    GR_FIELD_PHI_EM, GR_FIELD_A_X, GR_FIELD_A_Y,
} from './config';

type ViewMode = 'scalar' | 'curl' | 'gradmag';

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
];

let scratch: Float32Array | null = null;
function buf(n: number): Float32Array {
    if (!scratch || scratch.length !== n) scratch = new Float32Array(n);
    return scratch;
}

/* Returns a Float32Array for the selected view.  For scalar views this is the
 * live WASM heap view (no copy); for derived views it's a reused scratch
 * buffer.  Either way, valid until the next computeView or sim step. */
export function computeView(world: World, viewIndex: number): Float32Array {
    const view = FIELD_VIEWS[viewIndex] ?? FIELD_VIEWS[0];
    const W = world.W, H = world.H;
    if (view.mode === 'scalar') return world.fieldView(view.a);

    const out = buf(W * H);
    const at = (f: Float32Array, i: number, j: number) =>
        f[Math.min(H - 1, Math.max(0, j)) * W + Math.min(W - 1, Math.max(0, i))];

    if (view.mode === 'curl') {
        const Ax = world.fieldView(view.a);
        const Ay = world.fieldView(view.b!);
        for (let j = 0; j < H; j++) {
            for (let i = 0; i < W; i++) {
                const dAydx = (at(Ay, i + 1, j) - at(Ay, i - 1, j)) * 0.5;
                const dAxdy = (at(Ax, i, j + 1) - at(Ax, i, j - 1)) * 0.5;
                out[j * W + i] = dAydx - dAxdy;
            }
        }
    } else { // gradmag
        const P = world.fieldView(view.a);
        for (let j = 0; j < H; j++) {
            for (let i = 0; i < W; i++) {
                const gx = (at(P, i + 1, j) - at(P, i - 1, j)) * 0.5;
                const gy = (at(P, i, j + 1) - at(P, i, j - 1)) * 0.5;
                out[j * W + i] = Math.hypot(gx, gy);
            }
        }
    }
    return out;
}
