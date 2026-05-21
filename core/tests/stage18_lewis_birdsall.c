/* Stage 18 — Lewis-Birdsall energy-conserving force interpolation.
 *
 * Validates the new GR_FORCE_INTERP_LEWIS_BIRDSALL scheme (both LB-CIC and
 * LB-TSC) against the LEGACY (FD-then-interp) scheme on the two existing
 * stage10 phases that probe (1) static self-force and (2) Phase E orbit
 * heating at m=1e-3.
 *
 *   [A] HE self-force cancellation under LB-CIC.  Replicates stage10
 *       Phase C: stationary particle at the box center, 10000 steps,
 *       m=1e-3.  Drift should remain ~0 — LB is in the HE-adjoint family.
 *
 *   [B] HE self-force cancellation under LB-TSC.  Same setup with TSC
 *       shape function.
 *
 *   [C] Phase E orbit comparison.  pic_orbiting, m=1e-3, r=20, 1 T_ana.
 *       Compares M2 (orbit radial drift) across:
 *         legacy CIC, LB-CIC, LB-TSC.
 *       Expectation per Birdsall-Maron 1980: LB reduces moving-particle
 *       heating by 1-2 orders of magnitude vs the FD-then-interp scheme. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...)                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);           \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

/* ----------------------------------------------------------------------- */
/* [A,B] HE self-force cancellation under LB                               */
/* ----------------------------------------------------------------------- */
static int test_lb_self_force(gr_shape_function_t shape, const char* label) {
    printf("\n[%s] HE self-force cancellation under LB\n", label);
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 1.0e-3f;
    const int   N      = 10000;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    gr_sim_set_damping(sim, 16);
    const float params[1] = {mass};
    TEST_ASSERT(gr_sim_load_scenario(sim, "pic_static", params, 1) == 0,
                "pic_static load failed");
    /* Override the scenario's defaults to exercise the variant under test. */
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);

    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float x0 = p->x;
    const float y0 = p->y;
    gr_sim_step_n(sim, N);
    p = gr_sim_get_particle(sim, 0);
    const float drift = sqrtf((p->x - x0) * (p->x - x0)
                            + (p->y - y0) * (p->y - y0)) / dx;
    printf("  %s, %d steps: drift = %.3e cells\n", label, N, drift);
    TEST_ASSERT(drift < 1.0e-3f,
                "LB self-force drift %.3e cells exceeds 1e-3 — HE adjoint broken?", drift);
    gr_sim_destroy(sim);
    return 0;
}

/* ----------------------------------------------------------------------- */
/* [C] Phase E orbit comparison                                            */
/* ----------------------------------------------------------------------- */
typedef struct {
    const char*         label;
    gr_shape_function_t shape;
    gr_force_interp_t   force;
} variant_t;

static float orbit_drift_pct(const variant_t* v, int* nan_out) {
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
    if (!sim) { *nan_out = 1; return 0.0f; }
    gr_sim_set_damping(sim, 16);
    if (gr_sim_load_scenario(sim, "pic_orbiting", par, 4) != 0) {
        gr_sim_destroy(sim);
        *nan_out = 1;
        return 0.0f;
    }
    /* Override the scenario's defaults AFTER loading so the variant
     * under test is what actually runs.  The scenario sets TSC + LB by
     * default (production); this test wants to compare against legacy. */
    gr_sim_set_shape_function(sim, v->shape);
    gr_sim_set_force_interp(sim, v->force);
    const int N_steps = (int) (T_ana / gr_sim_dt(sim));
    gr_sim_step_n(sim, N_steps);
    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float r = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
    *nan_out = !isfinite(r);
    gr_sim_destroy(sim);
    return 100.0f * (r - r_orb) / r_orb;
}

static int test_phase_e_compare(void) {
    printf("\n[C] Phase E orbit, m=1e-3, r=20, 1 T_ana — variant comparison\n");

    const variant_t variants[] = {
        { "CIC legacy ", GR_SHAPE_CIC, GR_FORCE_INTERP_LEGACY         },
        { "CIC LB     ", GR_SHAPE_CIC, GR_FORCE_INTERP_LEWIS_BIRDSALL },
        { "TSC legacy ", GR_SHAPE_TSC, GR_FORCE_INTERP_LEGACY         },
        { "TSC LB     ", GR_SHAPE_TSC, GR_FORCE_INTERP_LEWIS_BIRDSALL },
    };
    const int n_var = sizeof(variants) / sizeof(variants[0]);

    printf("  %-12s   M2 (radial drift %%)\n", "variant");
    printf("  ------------------------------------\n");
    for (int i = 0; i < n_var; i++) {
        int nan = 0;
        const float m2 = orbit_drift_pct(&variants[i], &nan);
        printf("  %-12s   %+8.3f%%%s\n", variants[i].label, m2,
               nan ? "  [NaN]" : "");
        fflush(stdout);
        TEST_ASSERT(!nan, "variant %s produced NaN", variants[i].label);
    }
    return 0;
}

int main(void) {
    printf("=== stage18_lewis_birdsall ===\n");
    if (test_lb_self_force(GR_SHAPE_CIC, "A: LB-CIC self-force") != 0) return 1;
    if (test_lb_self_force(GR_SHAPE_TSC, "B: LB-TSC self-force") != 0) return 1;
    if (test_phase_e_compare()                                   != 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
