/* World: an ergonomic wrapper around a single sim handle + GRliteCore.
 *
 * Owns the grid dimensions and the WASM heap views.  Presents typed
 * particle reads, field access, config, scenario loading, and stepping.
 * This is the single API that the renderer, the UI, and the debug bridge
 * all call -- nothing else should touch the raw core or HEAP directly. */

import type { GRliteCore } from '../wasm/binding';
import {
    GRID_W, GRID_H, DX, C_EFF, CFL, N_DAMPING, N_SETTLE,
    GR_SHAPE_BUMP, GR_FORCE_INTERP_LB, KERNEL_RADIUS_CELLS,
    FRICTION_FLOOR, FRICTION_WALL, GR_FIELD_PHI_GRAV,
    PARTICLE_STRIDE_F32_OLD, PARTICLE_STRIDE_F32_NEW, SCENARIOS,
} from './config';

export interface Particle {
    index: number;
    x: number; y: number;
    px: number; py: number;
    mass: number; charge: number;
    properTime: number;
    spin: number; phiSpin: number; gFactor: number;
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

    /* ---- absorber / deposition config (v42 production defaults) ---- */
    applyAbsorberDefaults(): void {
        const c = this.core, s = this.sim;
        c.setOuterBcNeumann(s, 0);          // Dirichlet (pins DC mode)
        c.setDamping(s, 0);                  // disable legacy multiplicative ring
        c.setVolumeFrictionTaper(s, FRICTION_FLOOR, FRICTION_WALL, N_DAMPING);
        c.setZeroMeanScalarPotentials(s, 0);
        c.setShapeFunction(s, GR_SHAPE_BUMP);
        c.setKernelRadius(s, KERNEL_RADIUS_CELLS);
        c.setForceInterp(s, GR_FORCE_INTERP_LB);
    }

    setParticlesFrozen(frozen: boolean): void {
        this.core.setParticlesFrozen(this.sim, frozen ? 1 : 0);
    }
    setGEff(g: number): void { this.core.setGEff(this.sim, g); }
    setKE(k: number): void { this.core.setKE(this.sim, k); }

    /* ---- scenario loading (Phase 1: still via C scenarios) ---- */
    loadScenario(name: string): number {
        const spec = SCENARIOS[name];
        if (!spec) { console.error(`unknown scenario: ${name}`); return -1; }
        const M = this.core.M;
        const arr = new Float32Array(spec.params);
        const ptr = M._malloc(arr.byteLength);
        M.HEAPF32.set(arr, ptr >> 2);
        const rc = this.core.loadScenario(this.sim, spec.name, ptr, arr.length);
        M._free(ptr);
        if (rc === 0) this.detectStride();
        return rc;
    }

    /* Full (re)build of a scenario: load + production config + settle init.
     * After this the field is at the discrete fixed point and the particles
     * are ready to release. */
    buildScenario(name: string): number {
        const rc = this.loadScenario(name);
        if (rc !== 0) return rc;
        this.applyAbsorberDefaults();
        this.core.initPotentialsSettled(this.sim, N_SETTLE);
        return 0;
    }

    /* ---- field access ---- */
    fieldView(which = GR_FIELD_PHI_GRAV): Float32Array {
        const ptr = this.core.fieldPtr(this.sim, which);
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
        return {
            index: i,
            x: v[0], y: v[1], px: v[2], py: v[3],
            mass: v[4], charge: v[5], properTime: v[6],
            spin: this.stride >= 10 ? v[7] : 0,
            phiSpin: this.stride >= 10 ? v[8] : 0,
            gFactor: this.stride >= 10 ? v[9] : 2,
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
