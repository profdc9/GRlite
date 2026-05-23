/* Stage 36 -- DIAGNOSTIC: flip the sign on the -q d_t A inductive term.
 *
 * Stage 30 (post-v36 standard-Maxwell flip) shows OUTSPIRAL +25% with the
 * inductive piece +1 (default).  The diagnostic hypothesis: in this
 * convention pair (Maxwell-style field + Maxwell-style force law) the
 * discrete radiation-reaction direction is reversed.  Test by simply
 * flipping the inductive-piece sign via gr_sim_set_em_inductive_sign(-1).
 *
 * Setup: same as Stage 30 (Q=+1 background, q=-1e-3, m=1e-3, r=20, 4 orbits)
 * but at a single CFL (skip the sweep).  Run inductive +1 and inductive
 * -1 side by side. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); \
        return 1; \
    } \
} while (0)

static float run(float inductive_sign, int* n_completed_out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float Q      = +1.0f;
    const float q_test = -1.0e-3f;
    const float m_test = +1.0e-3f;
    const float r_orb  = 20.0f;
    const int   N      = 4;
    const float k_e    = 1.0f;

    const float g_mag  = fabsf(q_test * Q) * k_e * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f) / m_test;
    const float rg     = r_orb * g_mag;
    const float rg2_c2 = rg * rg / (c_eff * c_eff);
    const float u_v2   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ = sqrtf(u_v2);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_circ;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return NAN;
    gr_sim_set_damping(sim, 16);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    gr_sim_set_shape_function(sim, GR_SHAPE_CIC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEGACY);
    gr_sim_set_rho_smooth_passes(sim, 4);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, inductive_sign);   /* +1 or -1 */

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test, 0.0f, v_circ);

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int)(1.2f * (float) N * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    float r_at_end = r_orb;
    for (int s = 0; s < n_max && wraps < N; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) {
            if (n_completed_out) *n_completed_out = wraps;
            gr_sim_destroy(sim);
            return NAN;
        }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            r_at_end = r_now;
        }
        th_prev = th;
    }
    if (n_completed_out) *n_completed_out = wraps;
    gr_sim_destroy(sim);
    if (wraps < N) return NAN;
    return (r_at_end - r_orb) / r_orb;
}

int main(void) {
    printf("=== stage36_em_inductive_signflip ===\n");
    printf("Stage 30 setup at CFL=1/sqrt(2), 4 orbits.\n");
    printf("Compare inductive_sign = +1 (default) vs -1 (flipped).\n\n");

    int n_pos = 0, n_neg = 0;
    const float d_pos = run(+1.0f, &n_pos);
    const float d_neg = run(-1.0f, &n_neg);
    printf("inductive_sign = +1 (default):  drift = %+7.4f%%   (n_completed=%d)\n",
           100.0 * (double) d_pos, n_pos);
    printf("inductive_sign = -1 (flipped):  drift = %+7.4f%%   (n_completed=%d)\n",
           100.0 * (double) d_neg, n_neg);
    printf("\n");
    if (d_pos > 0.0f && d_neg < 0.0f) {
        printf("CONFIRMED: flipping the inductive sign turns OUTSPIRAL into INSPIRAL.\n");
        printf("           Post-v36 standard-Maxwell convention needs the inductive\n");
        printf("           term sign to be FLIPPED for radiation-reaction to give\n");
        printf("           the physically correct direction.\n");
    } else if (d_pos < 0.0f && d_neg > 0.0f) {
        printf("Unexpected: default already inspirals; flip makes it outspiral.\n");
    } else {
        printf("Both same sign -- inductive sign doesn't dominate the drift.\n");
    }
    TEST_ASSERT(n_pos == 4 && n_neg == 4, "didn't complete 4 orbits");
    printf("\nDIAGNOSTIC COMPLETE.\n");
    return 0;
}
