/* Stage 51 -- v39 per-particle self-field subtraction validation.
 *
 * Sanity test for the new gr_sim_particle_enable_self_field path.
 *
 * Setup mirrors Stage 37e: single +Q charge, EM-only, BUMP R=8, mass=1,
 * equal-distance methodology (8-cell warmup + 2-cell measurement).
 *
 * Three configurations per velocity:
 *   (A) self-field DISABLED (baseline = current PIC behavior, Stage 37e)
 *   (B) self-field ENABLED, eps = 0 (full cancellation predicted)
 *   (C) self-field ENABLED, eps = -1 (no subtraction; should match (A))
 *
 * Prediction: (A) shows the v-dependent spurious force we already know.
 * (B) should drop F to float32 noise at every v (the architecture
 * working as designed).  (C) verifies that eps=-1 gives back the
 * collective-only behavior, so (A) and (C) agree to numerical noise. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px*p->px + p->py*p->py) / (m*m*c_eff*c_eff));
    return p->py / (gamma * m);
}

typedef enum { CFG_DISABLED, CFG_EPS_0, CFG_EPS_NEG1 } cfg_t;

static float measure(cfg_t cfg, float v_drift) {
    const int W = 128, H = 1024, n_damp = 16;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass = 1.0f, Q = 0.01f;
    const float y_start = (float)(n_damp + 32) * dx;
    const float dt = cfl * dx / c_eff;
    const int n_warmup = (int) ceilf(8.0f / (v_drift * dt));
    const int n_measure = (int) ceilf(2.0f / (v_drift * dt));

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 8.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);

    const float cx = ((float)(W-1) * 0.5f) * dx;
    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    if (cfg == CFG_EPS_0 || cfg == CFG_EPS_NEG1) {
        gr_sim_particle_enable_self_field(sim, 0);
        if (cfg == CFG_EPS_NEG1) {
            gr_sim_particle_set_self_field_epsilon(sim, 0, -1.0f, -1.0f);
        }
        /* eps defaults to 0 for CFG_EPS_0. */
    }

    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);
    const float vy_w = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    for (int s = 0; s < n_measure; s++) gr_sim_step(sim);
    const float vy_m = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    const float F = mass * (vy_m - vy_w) / ((float)n_measure * dt);
    gr_sim_destroy(sim);
    return F;
}

int main(void) {
    printf("=== stage51_self_field_validation ===\n");
    printf("Single +Q at constant v.  BUMP R=8.  EM-only.  Equal-distance.\n");
    printf("(A) self-field OFF (baseline).\n");
    printf("(B) self-field ON, eps=0 (predict: F approx 0 across v).\n");
    printf("(C) self-field ON, eps=-1 (predict: F matches (A)).\n\n");

    const float vs[] = {0.005f, 0.010f, 0.020f, 0.040f, 0.050f, 0.075f, 0.100f};
    const int Nv = sizeof(vs)/sizeof(vs[0]);
    printf("  v/c     (A) off         (B) eps=0       (C) eps=-1\n");
    printf("  ------------------------------------------------------\n");
    for (int iv = 0; iv < Nv; iv++) {
        const float FA = measure(CFG_DISABLED, vs[iv]);
        const float FB = measure(CFG_EPS_0,    vs[iv]);
        const float FC = measure(CFG_EPS_NEG1, vs[iv]);
        printf("  %.3f   %+11.4e   %+11.4e   %+11.4e\n",
               (double)vs[iv], (double)FA, (double)FB, (double)FC);
    }
    return 0;
}
