/* Stage 53b -- find the discrete-circular velocity by scanning v_factor.
 *
 * Stage 53 ruled out delayed potential as the rosette source (c-scan
 * showed amplitude is c-independent at ~7.1%).  Hypothesis under test
 * here: the scenario's v_orb = Q sqrt(k_e/m) is the CONTINUUM circular
 * velocity, but the discrete kernel-gather force at d=40 may differ
 * from k_e Q^2/d.  Using the continuum v_orb makes the IC slightly
 * non-circular; Bertrand's theorem then turns the small mismatch into
 * the visible rosette.
 *
 * Method: scan v_factor (multiplier on v_orb) and report the rosette
 * range.  The discrete-circular value is the v_factor that MINIMIZES
 * the rosette amplitude.
 *
 * BUMP R=8, self-field ON (so we're seeing only the IC-mismatch rosette,
 * not any self-force contamination).  EM only (GR mirrored Stage 53). */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run(float v_factor) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass_inertia = 1.0f, Q = 0.01f, r_orb = 20.0f * dx;
    const float k_e = 1.0f;
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float v_orb_continuum = Q * sqrtf(k_e / mass_inertia);
    const float T_ana = 2.0f * (float)M_PI * r_orb / v_orb_continuum;
    const int n_orbits = 2;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    gr_sim_set_k_e(sim, k_e);
    const float params[6] = {mass_inertia, r_orb, v_factor, cx, cy, Q};
    gr_sim_load_scenario(sim, "pic_binary_em", params, 6);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 8.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_particle_enable_self_field(sim, 0);
    gr_sim_particle_enable_self_field(sim, 1);

    const float dt = gr_sim_dt(sim);
    const int n_steps = (int)((float)n_orbits * T_ana / dt);
    float d_min = 40.0f, d_max = 40.0f;
    int unbound = 0;
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        if (!isfinite(p0->x)) { unbound = 1; break; }
        if ((s % 64) == 0) {
            const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
            const float d = sqrtf(dxs*dxs + dys*dys);
            if (d < d_min) d_min = d;
            if (d > d_max) d_max = d;
        }
    }
    if (unbound) {
        printf("  v_factor=%.4f    UNBOUND\n", (double)v_factor);
    } else {
        const float dlow  = 100.0f * (d_min - 40.0f) / 40.0f;
        const float dhigh = 100.0f * (d_max - 40.0f) / 40.0f;
        printf("  v_factor=%.4f   d_min=%6.3f (%+5.2f%%)  d_max=%6.3f (%+5.2f%%)  range=%5.2f%%\n",
               (double)v_factor,
               (double)d_min, (double)dlow,
               (double)d_max, (double)dhigh,
               (double)(dhigh - dlow));
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage53b_v_orb_scan ===\n");
    printf("EM binary, self-field ON, BUMP R=8.  Scan v_factor to find\n");
    printf("the v that minimizes the rosette amplitude.\n\n");

    const float vfs[] = {0.95f, 0.97f, 0.98f, 0.99f, 0.995f, 1.000f,
                         1.005f, 1.01f, 1.02f, 1.03f, 1.05f};
    for (int i = 0; i < (int)(sizeof(vfs)/sizeof(vfs[0])); i++) {
        run(vfs[i]);
    }
    return 0;
}
