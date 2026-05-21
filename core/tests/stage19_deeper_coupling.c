/* Stage 19 — how far does the new orbital-stability regime extend?
 *
 * After the read-prev timing fix and the TSC+LB pic_orbiting defaults,
 * the orbit at m_test = 1e-3 is essentially stable over 4 analytic
 * periods (stage 18b shows ~-3% per orbit, consistent with physical
 * radiation reaction at v/c ~ 0.22).
 *
 * This sweep pushes m_test up to find where the new regime breaks down.
 * Tests:
 *
 *   m_test in {1e-3, 1e-2, 1e-1, 1.0}.
 *
 * For each m_test, we report the radial drift at orbits 1 and 4 under
 * the production defaults (TSC + LB + 4-pass smoothing).  Also reports
 * CIC legacy at the same m_test for comparison — to see whether the
 * stable regime is genuinely the property of "correct timing" or
 * whether something else degrades at higher coupling.
 *
 * Theoretical scaling: PIC deposit aliasing scales as m_test^2 (the
 * noise injection rate is proportional to m, and the particle sees
 * that noise filtered through its own deposit also proportional to m).
 * So if legacy unbinds in 3 orbits at m=1e-3 (drift ~+140%/3 orbits ~
 * +50%/orbit average), then at m=1e-2 we'd expect ~+5000%/orbit ->
 * instant unbinding.  Anything stable at m=1e-2 is a real regime
 * win. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float orbit1_drift;
    float orbit4_drift;
    int   unbound_by_orbit;  /* 0 if stayed bound through orbit 4 */
    int   nan;
} result_t;

static void run_coupling(float m_test, gr_shape_function_t shape,
                         gr_force_interp_t force, result_t* out) {
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float GM     = 1.0f;
    const float r_orb  = 20.0f;
    const float eps    = 1.0f;
    const float par[4] = {GM, r_orb, eps, m_test};
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;
    const float T_ana  = 2.0f * (float) M_PI * sqrtf(r_orb * r_orb * r_orb / GM);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_damping(sim, 16);
    gr_sim_load_scenario(sim, "pic_orbiting", par, 4);
    /* Override scenario defaults with the variant under test. */
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, force);
    const int steps_per_orbit = (int) (T_ana / gr_sim_dt(sim));

    out->orbit1_drift = 0.0f;
    out->orbit4_drift = 0.0f;
    out->unbound_by_orbit = 0;
    out->nan = 0;
    for (int k = 1; k <= 4; k++) {
        gr_sim_step_n(sim, steps_per_orbit);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float r = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
        if (!isfinite(r)) { out->nan = 1; break; }
        const float drift = 100.0f * (r - r_orb) / r_orb;
        if (k == 1) out->orbit1_drift = drift;
        if (k == 4) out->orbit4_drift = drift;
        if (out->unbound_by_orbit == 0 && fabsf(drift) > 100.0f) out->unbound_by_orbit = k;
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage19_deeper_coupling ===\n");
    printf("Push m_test up under production (TSC+LB+smoothing) vs CIC legacy.\n");
    printf("Reports drift at orbits 1 and 4.  'unbound' = drift > 100%% before orbit 4.\n\n");

    const float ms[] = {1.0e-3f, 1.0e-2f, 1.0e-1f, 1.0f};
    const int   n_m = sizeof(ms) / sizeof(ms[0]);

    printf("%-8s | %-22s | %-22s\n",
           "m_test", "production (TSC+LB)", "CIC legacy");
    printf("%-8s | %-9s %-9s | %-9s %-9s\n",
           "", "orbit 1", "orbit 4", "orbit 1", "orbit 4");
    printf("----------------------------------------------------------\n");
    for (int i = 0; i < n_m; i++) {
        result_t prod, leg;
        run_coupling(ms[i], GR_SHAPE_TSC, GR_FORCE_INTERP_LEWIS_BIRDSALL, &prod);
        run_coupling(ms[i], GR_SHAPE_CIC, GR_FORCE_INTERP_LEGACY,         &leg);

        char prod_o4[32], leg_o4[32];
        if (prod.nan) snprintf(prod_o4, sizeof(prod_o4), "  [NaN]");
        else if (prod.unbound_by_orbit) snprintf(prod_o4, sizeof(prod_o4), "unbound@%d", prod.unbound_by_orbit);
        else snprintf(prod_o4, sizeof(prod_o4), "%+6.2f%%", prod.orbit4_drift);
        if (leg.nan) snprintf(leg_o4, sizeof(leg_o4), "  [NaN]");
        else if (leg.unbound_by_orbit) snprintf(leg_o4, sizeof(leg_o4), "unbound@%d", leg.unbound_by_orbit);
        else snprintf(leg_o4, sizeof(leg_o4), "%+6.2f%%", leg.orbit4_drift);

        printf("%-8.0e | %+7.2f%%  %-9s | %+7.2f%%  %-9s\n",
               (double) ms[i],
               prod.orbit1_drift, prod_o4,
               leg.orbit1_drift, leg_o4);
        fflush(stdout);
    }

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
