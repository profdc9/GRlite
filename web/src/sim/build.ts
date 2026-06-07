/* Build a sim from a Scenario via the WASM API.  This is the JSON-driven
 * replacement for gr_sim_load_scenario (which stays for native tests).
 *
 * Assumes the World was created with grid == scenario.grid (the caller
 * recreates the World when the grid changes). */

import type { World } from './world';
import type { Scenario, ShapeName, ForceInterpName } from './scenario';
import { GR_SHAPE_CIC, GR_SHAPE_TSC, GR_SHAPE_BUMP,
         GR_FORCE_INTERP_LEGACY, GR_FORCE_INTERP_LB,
         GR_FORCE_NEWTONIAN, GR_FORCE_RELATIVISTIC } from './config';

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
    /* Zero the field buffers so every (re)build starts from a quiet field.  The
     * L-W/settled inits overwrite the field anyway, but 'zero'/'none' rely on it
     * being clean -- otherwise reset keeps the previous run's wave. */
    c.clearFields(s);

    c.setGEff(s, g.gEff);
    c.setKE(s, g.kE);

    c.setShapeFunction(s, shapeId(g.shape));
    if (g.shape === 'bump') c.setKernelRadius(s, g.kernelRadius);
    c.setForceInterp(s, forceId(g.forceInterp));
    c.setForceTier(s, g.forceTier === 'relativistic' ? GR_FORCE_RELATIVISTIC : GR_FORCE_NEWTONIAN);

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
    /* Shapiro c_local^2 is built from the background Phi_g, so (re)compute it
     * AFTER the background is installed -- the switch above ran before it. */
    if (sw.emShapiro) c.setEmShapiro(s, 1);
    /* Dynamic Shapiro: per-step c_local from total Phi_g (moving/deposited mass
     * lens).  Set after emShapiro so the recompute has the background. */
    c.setEmShapiroDynamic(s, b(sw.shapiroDynamic));

    scn.particles.forEach((p, i) => {
        c.addParticle(s, p.x, p.y, p.mass, p.charge, p.vx, p.vy);
        if (p.spin) c.setParticleSpin(s, i, p.spin, p.gFactor ?? 2.0);
        /* Driven source / pinned: set the oscillation and the forces-on flag.
         * forces === false also disables self-field below (a source should
         * radiate its full field and feel nothing). */
        const pinned = p.forces === false;
        if (p.drive || pinned) {
            const d = p.drive;
            c.setParticleDrive(s, i, d?.amp ?? 0, d?.omega ?? 0, d?.phase ?? 0,
                               d?.axis?.[0] ?? 1, d?.axis?.[1] ?? 0, pinned ? 0 : 1);
        }
        /* Per-particle selfField overrides the global default; never on a
         * pinned source. */
        const selfOn = !pinned && (p.selfField !== undefined ? p.selfField : g.selfFieldDefault);
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

    /* Field initialization -- derived from the physics config so it can't
     * desync from the switches (this is the single place both scenario load
     * and the reset button initialize through).  When perturbation is active
     * (deposit + evolve) the field MUST start at its settled fixed point, else
     * the run begins from zero with an unphysical transient -- so an explicit
     * 'lw-settled', OR any active-perturbation run, settles the field.  'lw'
     * is the unrelaxed direct sum; 'none' with perturbation off leaves zero. */
    const pertActive = sw.fieldEvolution && sw.particleSourceDeposition;
    const settleSteps = g.init.settleSteps > 0 ? g.init.settleSteps : Math.max(scn.grid.W, scn.grid.H);
    if (g.init.method === 'lw-settled' || (g.init.method === 'none' && pertActive)) {
        c.initPotentialsSettled(s, settleSteps);
    } else if (g.init.method === 'lw') {
        c.initPotentialsLW(s);
    }
    /* else: 'zero' (always), or 'none' + perturbation off -> leave fields zero. */

    /* Every (re)build starts a fresh run at t=0.  initPotentialsSettled already
     * zeroes step_count, but 'none'/'lw' inits do not -- without this the sim
     * clock keeps climbing across resets while per-particle proper_time restarts
     * at 0, desyncing the proper-/coordinate-time track ticks. */
    c.resetTime(s);

    world.refreshStride();
}
