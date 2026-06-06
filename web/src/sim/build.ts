/* Build a sim from a Scenario via the WASM API.  This is the JSON-driven
 * replacement for gr_sim_load_scenario (which stays for native tests).
 *
 * Assumes the World was created with grid == scenario.grid (the caller
 * recreates the World when the grid changes). */

import type { World } from './world';
import type { Scenario, ShapeName, ForceInterpName } from './scenario';
import { GR_SHAPE_CIC, GR_SHAPE_TSC, GR_SHAPE_BUMP,
         GR_FORCE_INTERP_LEGACY, GR_FORCE_INTERP_LB } from './config';

const shapeId = (s: ShapeName): number =>
    s === 'bump' ? GR_SHAPE_BUMP : s === 'tsc' ? GR_SHAPE_TSC : GR_SHAPE_CIC;
const forceId = (f: ForceInterpName): number =>
    f === 'lewis-birdsall' ? GR_FORCE_INTERP_LB : GR_FORCE_INTERP_LEGACY;
const b = (x: boolean): number => (x ? 1 : 0);

export function applyScenario(world: World, scn: Scenario): void {
    const c = world.core, s = world.sim, g = scn.global;

    c.clearParticles(s);
    c.clearSources(s);
    c.clearBackground(s);

    c.setGEff(s, g.gEff);
    c.setKE(s, g.kE);

    c.setShapeFunction(s, shapeId(g.shape));
    if (g.shape === 'bump') c.setKernelRadius(s, g.kernelRadius);
    c.setForceInterp(s, forceId(g.forceInterp));

    const sw = g.switches;
    /* No-radiation mode forces the inductive (radiation-reaction) pieces off. */
    const gmInductive = sw.gravitomagneticInductive && !g.noRadiation;
    const emInductive = sw.emInductive && !g.noRadiation;
    c.setGravitomagneticForce(s, b(sw.gravitomagneticForce));
    c.setGravitomagneticInductive(s, b(gmInductive));
    c.setEmLorentz(s, b(sw.emLorentz));
    c.setEmInductive(s, b(emInductive));
    c.setEmElectrostatic(s, b(sw.emElectrostatic));
    c.setEmMagnetic(s, b(sw.emMagnetic));
    c.setEmStressEnergy(s, b(sw.emStressEnergy));
    c.setEmShapiro(s, b(sw.emShapiro));
    c.setFieldEvolution(s, b(sw.fieldEvolution));
    c.setParticleSourceDeposition(s, b(sw.particleSourceDeposition));
    c.setEsirkepov(s, b(sw.esirkepov));
    c.setPeriodicBC(s, b(sw.periodicBC));
    c.setRhoSmoothPasses(s, g.rhoSmooth);
    c.setJSmoothPasses(s, g.jSmooth);

    c.setBgMode(s, g.bgMode === 'analytic' ? 1 : 0);
    /* Unified compact body (M, Q, Jz) -- superposes grav + Coulomb + frame
     * dragging in one background.  The classic metrics are special cases. */
    for (const bg of scn.background) {
        c.setBackgroundBody(s, bg.x, bg.y, bg.GM ?? 0, bg.Q ?? 0, bg.Jz ?? 0, bg.gFactor ?? 2, bg.epsilon);
    }

    scn.particles.forEach((p, i) => {
        c.addParticle(s, p.x, p.y, p.mass, p.charge, p.vx, p.vy);
        if (p.spin) c.setParticleSpin(s, i, p.spin, p.gFactor ?? 2.0);
        /* Per-particle selfField overrides the global default. */
        const selfOn = p.selfField !== undefined ? p.selfField : g.selfFieldDefault;
        if (selfOn) {
            c.particleEnableSelfField(s, i);
            if (p.selfFieldEps) c.particleSetSelfFieldEpsilon(s, i, p.selfFieldEps[0], p.selfFieldEps[1]);
        }
    });

    /* Absorber (v42): Dirichlet/Neumann + derivative friction, multiplicative
     * ring off.  Must precede the settle init (which uses the friction). */
    const ab = g.absorber;
    c.setDamping(s, 0);
    c.setOuterBcNeumann(s, b(ab.outerBC === 'neumann'));
    c.setVolumeFrictionTaper(s, ab.frictionFloor, ab.frictionWall, ab.frictionDepth);
    c.setZeroMeanScalarPotentials(s, b(ab.zeroMeanScalar));

    /* Field initialization. */
    if (g.init.method === 'lw-settled') c.initPotentialsSettled(s, g.init.settleSteps);
    else if (g.init.method === 'lw') c.initPotentialsLW(s);
    /* 'none' -> leave fields zero. */

    /* Every (re)build starts a fresh run at t=0.  initPotentialsSettled already
     * zeroes step_count, but 'none'/'lw' inits do not -- without this the sim
     * clock keeps climbing across resets while per-particle proper_time restarts
     * at 0, desyncing the proper-/coordinate-time track ticks. */
    c.resetTime(s);

    world.refreshStride();
}
