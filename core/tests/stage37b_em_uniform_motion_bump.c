/* Stage 37b -- uniform-motion EM PIC self-force, redone with v38 kernels.
 *
 * Dan, 2026-05-26: "Could you redo the stage where an EM particle moves
 * at a constant velocity and see if heating happens there?"
 *
 * Physics: a single free charge moving at constant v in vacuum CANNOT
 * radiate (Larmor requires acceleration).  Any KE drift in this PIC
 * setup is a self-force artifact.
 *
 * Original Stage 37 (TSC bare, v=0.001): KE rose +1587% over the run
 * -- the "uniform-motion heating" failure the bump kernel and the
 * smoothing work were intended to address.
 *
 * This redo runs the SAME asymmetric-domain setup (128 x 1024 cells,
 * particle launched at low y with +y velocity, no background) with
 * the v38 production kernels:
 *   - TSC bare           (baseline failure)
 *   - TSC + N=16         (smoothing fix)
 *   - BUMP R=6           (wide kernel)
 *   - BUMP R=8           (wider kernel)
 *
 * At three drift velocities to expose v-scaling:
 *   v = 0.1 (Stage 48 Case A regime)
 *   v = 0.01 (Stage 48 Case C regime)
 *   v = 0.001 (original Stage 37's worst case) */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float KE0, KEf, vy0, vyf, y0, yf;
    int   n_steps;
    int   bailed;
} run_t;

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px*p->px + p->py*p->py) / (m*m*c_eff*c_eff));
    return p->py / (gamma * m);
}

static float KE_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px*p->px + p->py*p->py) / (m*m*c_eff*c_eff));
    return (gamma - 1.0f) * m * c_eff * c_eff;
}

static run_t run(gr_shape_function_t shape, float Rk, int smooth, float v_drift) {
    const int W = 128, H = 1024, n_damp = 16;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass = 0.01f, Q = 0.01f;
    const float y_start = (float)(n_damp + 32) * dx;
    const float y_safe = (float)(H - n_damp - 64) * dx;
    const float dt = cfl * dx / c_eff;
    int n_steps = (int)((y_safe - y_start) / (v_drift * dt));
    if (n_steps > 20000) n_steps = 20000;
    if (n_steps < 200)   n_steps = 200;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, smooth);
    gr_sim_set_j_smooth_passes(sim, smooth);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* All EM force pieces ON.  No gravity. */
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);

    const float cx = ((float)(W-1) * 0.5f) * dx;
    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    /* Warmup: let the wake form before measuring (~8 cell traversal times). */
    int n_warmup = (int)(8.0f / (v_drift * dt));
    if (n_warmup > 2000) n_warmup = 2000;
    if (n_warmup < 50)   n_warmup = 50;
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    run_t r = {0};
    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    r.KE0 = KE_of(p0, c_eff);
    r.vy0 = vy_of(p0, c_eff);
    r.y0  = p0->y;

    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        if (p->y > y_safe || p->y < y_start - 8.0f) break;
        if (!isfinite(p->y)) { r.bailed = 1; break; }
        r.n_steps++;
    }
    const gr_particle_t* pf = gr_sim_get_particle(sim, 0);
    r.KEf = KE_of(pf, c_eff);
    r.vyf = vy_of(pf, c_eff);
    r.yf  = pf->y;
    gr_sim_destroy(sim);
    return r;
}

static void show(const char* tag, run_t r) {
    const double dKE = 100.0 * (r.KEf - r.KE0) / r.KE0;
    const double dv  = 100.0 * (r.vyf - r.vy0) / r.vy0;
    printf("  %-20s dKE=%+8.3f%%  dvy=%+8.3f%%   vy: %.5f -> %.5f   y: %.0f -> %.0f   steps=%d%s\n",
           tag, dKE, dv, (double)r.vy0, (double)r.vyf, (double)r.y0, (double)r.yf, r.n_steps,
           r.bailed ? "  [BAILED]" : "");
}

int main(void) {
    printf("=== stage37b_em_uniform_motion_bump ===\n");
    printf("Single +Q in 128 x 1024 asymmetric domain, no background.\n");
    printf("Constant-velocity charge in vacuum MUST NOT heat -- any dKE is artifact.\n\n");

    const float vs[] = {0.1f, 0.01f, 0.001f};
    for (int iv = 0; iv < 3; iv++) {
        printf("--- v_drift = %.3f c ---\n", (double)vs[iv]);
        show("TSC bare",     run(GR_SHAPE_TSC,  0.0f, 0,  vs[iv]));
        show("TSC + N=16",   run(GR_SHAPE_TSC,  0.0f, 16, vs[iv]));
        show("BUMP R=6",     run(GR_SHAPE_BUMP, 6.0f, 0,  vs[iv]));
        show("BUMP R=8",     run(GR_SHAPE_BUMP, 8.0f, 0,  vs[iv]));
        printf("\n");
    }
    return 0;
}
