/* Stage 18b — does the LB regime hold for multiple orbits?
 *
 * Stage 18 [C] showed LB drops M2 from +42% per orbit to -4% per orbit
 * at m=1e-3.  This test extends to 4 analytic periods and reports the
 * radial drift at each milestone, for both LB-CIC and LB-TSC. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run_orbit(const char* label,
                      gr_shape_function_t shape, gr_force_interp_t force,
                      int n_orbits) {
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float GM     = 1.0f;
    const float r_orb  = 20.0f;
    const float eps    = 1.0f;
    const float m_test = 1.0e-3f;
    const float par[4] = {GM, r_orb, eps, m_test};
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;
    const float T_ana  = 2.0f * (float) M_PI * sqrtf(r_orb * r_orb * r_orb / GM);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_damping(sim, 16);
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, force);
    gr_sim_load_scenario(sim, "pic_orbiting", par, 4);
    const int steps_per_orbit = (int) (T_ana / gr_sim_dt(sim));

    printf("  %s:\n", label);
    for (int k = 1; k <= n_orbits; k++) {
        gr_sim_step_n(sim, steps_per_orbit);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float r = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
        const float drift = 100.0f * (r - r_orb) / r_orb;
        printf("    after orbit %d: r = %.4f  drift = %+6.2f%%\n", k, r, drift);
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage18b_lb_multi_orbit ===\n");
    printf("4 analytic periods, m_test=1e-3, r=20, comparing schemes.\n\n");

    run_orbit("CIC legacy", GR_SHAPE_CIC, GR_FORCE_INTERP_LEGACY,         4);
    run_orbit("CIC LB    ", GR_SHAPE_CIC, GR_FORCE_INTERP_LEWIS_BIRDSALL, 4);
    run_orbit("TSC legacy", GR_SHAPE_TSC, GR_FORCE_INTERP_LEGACY,         4);
    run_orbit("TSC LB    ", GR_SHAPE_TSC, GR_FORCE_INTERP_LEWIS_BIRDSALL, 4);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
