/* World: an ergonomic wrapper around a single sim handle + GRliteCore.
 *
 * Owns the grid dimensions and the WASM heap views.  Presents typed
 * particle reads, field access, config, scenario loading, and stepping.
 * This is the single API that the renderer, the UI, and the debug bridge
 * all call -- nothing else should touch the raw core or HEAP directly. */

import type { GRliteCore } from '../wasm/binding';
import {
    GRID_W, GRID_H, DX, C_EFF, CFL, GR_FIELD_PHI_GRAV,
    PARTICLE_STRIDE_F32_OLD, PARTICLE_STRIDE_F32_NEW,
} from './config';

export interface Particle {
    index: number;
    x: number; y: number;
    px: number; py: number;
    mass: number; charge: number;
    properTime: number;
    spin: number; phiSpin: number; gFactor: number;
    /* v43 last-step force breakdown (0 if the WASM predates the wider stride). */
    fgravX: number; fgravY: number;
    femX: number; femY: number;
    fspinX: number; fspinY: number;
    ftotX: number; ftotY: number;
}

export class World {
    readonly core: GRliteCore;
    readonly sim: number;
    readonly W: number;
    readonly H: number;
    private stride: number;

    constructor(core: GRliteCore, w = GRID_W, h = GRID_H,
                dx = DX, cEff = C_EFF, cfl = CFL) {
        this.core = core;
        this.W = w;
        this.H = h;
        const sim = core.create(w, h, dx, cEff, cfl);
        if (!sim) throw new Error('gr_sim_create returned NULL');
        this.sim = sim;
        this.stride = PARTICLE_STRIDE_F32_NEW;
    }

    /* ---- stepping & time ---- */
    step(n = 1): void { this.core.stepN(this.sim, n); }
    stepCount(): number { return this.core.stepCount(this.sim); }
    time(): number { return this.core.simTime(this.sim); }
    dt(): number { return this.core.simDt(this.sim); }

    setParticlesFrozen(frozen: boolean): void {
        this.core.setParticlesFrozen(this.sim, frozen ? 1 : 0);
    }
    setGEff(g: number): void { this.core.setGEff(this.sim, g); }
    setKE(k: number): void { this.core.setKE(this.sim, k); }

    /* ---- field access ---- */
    fieldView(which = GR_FIELD_PHI_GRAV): Float32Array {
        const ptr = this.core.fieldPtr(this.sim, which);
        return new Float32Array(this.core.M.HEAPF32.buffer, ptr, this.W * this.H);
    }

    /* Background (analytic/sampled) grid for a field id, or null if no
     * background of that kind has been installed (the C side returns NULL). */
    backgroundView(which: number): Float32Array | null {
        const ptr = this.core.backgroundPtr(this.sim, which);
        if (!ptr) return null;
        return new Float32Array(this.core.M.HEAPF32.buffer, ptr, this.W * this.H);
    }

    /* ---- particles ---- */
    /* Re-probe the particle struct stride after a fresh build. */
    refreshStride(): void { this.detectStride(); }

    private detectStride(): void {
        if (this.core.particleCount(this.sim) >= 2) {
            const p0 = this.core.getParticle(this.sim, 0);
            const p1 = this.core.getParticle(this.sim, 1);
            if (p1 > p0) {
                const s = (p1 - p0) >> 2;
                if (s === PARTICLE_STRIDE_F32_OLD || s === PARTICLE_STRIDE_F32_NEW) {
                    this.stride = s;
                    this.core.particleStrideF32 = s;
                }
            }
        }
    }

    particleCount(): number { return this.core.particleCount(this.sim); }

    /* Zero-copy view of one particle's floats. */
    particleView(i: number): Float32Array {
        const ptr = this.core.getParticle(this.sim, i);
        return new Float32Array(this.core.M.HEAPF32.buffer, ptr, this.stride);
    }

    particle(i: number): Particle {
        const v = this.particleView(i);
        const hasF = this.stride >= 18;
        return {
            index: i,
            x: v[0], y: v[1], px: v[2], py: v[3],
            mass: v[4], charge: v[5], properTime: v[6],
            spin: this.stride >= 10 ? v[7] : 0,
            phiSpin: this.stride >= 10 ? v[8] : 0,
            gFactor: this.stride >= 10 ? v[9] : 2,
            fgravX: hasF ? v[10] : 0, fgravY: hasF ? v[11] : 0,
            femX: hasF ? v[12] : 0, femY: hasF ? v[13] : 0,
            fspinX: hasF ? v[14] : 0, fspinY: hasF ? v[15] : 0,
            ftotX: hasF ? v[16] : 0, ftotY: hasF ? v[17] : 0,
        };
    }

    particles(): Particle[] {
        const n = this.particleCount();
        const out: Particle[] = [];
        for (let i = 0; i < n; i++) out.push(this.particle(i));
        return out;
    }

    /* Live-edit a particle by writing the WASM heap view directly.  The
     * gr_particle_t fields are plain floats, so this needs no C export; the
     * pusher picks up the new values on the next step.  Caller is
     * responsible for flagging the run as live-modified (not reproducible). */
    setParticleFields(i: number, f: Partial<{
        x: number; y: number; px: number; py: number;
        mass: number; charge: number;
        spin: number; phiSpin: number; gFactor: number;
    }>): void {
        const v = this.particleView(i);
        if (f.x !== undefined) v[0] = f.x;
        if (f.y !== undefined) v[1] = f.y;
        if (f.px !== undefined) v[2] = f.px;
        if (f.py !== undefined) v[3] = f.py;
        if (f.mass !== undefined) v[4] = f.mass;
        if (f.charge !== undefined) v[5] = f.charge;
        if (this.stride >= 10) {
            if (f.spin !== undefined) v[7] = f.spin;
            if (f.phiSpin !== undefined) v[8] = f.phiSpin;
            if (f.gFactor !== undefined) v[9] = f.gFactor;
        }
    }
}
