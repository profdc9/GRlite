/* Stage 15 — critical-damping (gamma, gamma) ring sweep.
 *
 * Same axes and design as stage14_sigma_sweep, but the damping kernel
 * is GR_DAMP_TIME_CRITICAL: per cell, gamma = 1 - sigma dt, and the
 * recurrence Phi^{n+1} = 2 gamma Phi^n - gamma^2 Phi^{n-1} + driving
 * places both characteristic-equation roots at gamma (double real
 * root, critical damping), giving per-step decay = gamma exactly.
 *
 * Theory predicts ~2x stronger damping per unit sigma vs the
 * MULTIPLICATIVE form (where the actual decay is sqrt(1 - sigma dt)).
 * Equivalent question: does CRITICAL achieve the same M1 metric at
 * smaller sigma (smaller static-profile distortion)?
 *
 * Interior region held at 256x256 across configurations.  See
 * [[grlite-damping-sweep-result]] memory for the multiplicative
 * baseline numbers. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define INTERIOR 256

static float drift_off_center(const gr_damp_config_t* cfg, int* nan_out) {
    const int   W      = INTERIOR + 2 * cfg->n_damping;
    const int   H      = INTERIOR + 2 * cfg->n_damping;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 1.0e-3f;
    const int   N      = 2000;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float dx_off = 23.5f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { *nan_out = 1; return 0.0f; }
    gr_sim_set_damping_config(sim, cfg);
    gr_sim_set_particle_source_deposition(sim, 1);

    const float x0 = cx + dx_off;
    const float y0 = cy;
    gr_sim_add_particle(sim, x0, y0, mass, 0.0f, 0.0f, 0.0f);
    gr_sim_step_n(sim, N);

    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float dr = sqrtf((p->x - x0) * (p->x - x0) + (p->y - y0) * (p->y - y0));
    *nan_out = !isfinite(dr);
    gr_sim_destroy(sim);
    return dr / dx;
}

static float orbit_drift_pct(const gr_damp_config_t* cfg, int* nan_out) {
    const int   W      = INTERIOR + 2 * cfg->n_damping;
    const int   H      = INTERIOR + 2 * cfg->n_damping;
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
    gr_sim_set_damping_config(sim, cfg);
    if (gr_sim_load_scenario(sim, "pic_orbiting", par, 4) != 0) {
        gr_sim_destroy(sim);
        *nan_out = 1;
        return 0.0f;
    }
    const int N_steps = (int) (T_ana / gr_sim_dt(sim));
    gr_sim_step_n(sim, N_steps);
    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float r = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
    *nan_out = !isfinite(r);
    gr_sim_destroy(sim);
    return 100.0f * (r - r_orb) / r_orb;
}

int main(void) {
    printf("=== stage15_critical_sweep (polynomial m=3, "
           "TIME_CRITICAL, interior=%d) ===\n", INTERIOR);

    const int   Ns[]    = {8, 16, 32, 48};
    const int   n_N     = sizeof(Ns) / sizeof(Ns[0]);
    const float sdts[]  = {0.0f, 0.0001f, 0.001f, 0.01f, 0.05f, 0.1f, 0.2f, 0.4f, 0.7f, 1.0f};
    const int   n_sdt   = sizeof(sdts) / sizeof(sdts[0]);
    const float dt_ref  = 1.0f / sqrtf(2.0f);

    printf("\nM1: HE drift @ +23.5 cells (cells), 2000 steps, m=1e-3.\n");
    printf("M2: pic_orbiting radial drift at 1 T_ana, m=1e-3 (percent).\n");
    printf("sigma*dt at the outer wall as labeled.\n");
    printf("(stage14 baseline for comparison was MULTIPLICATIVE; this is CRITICAL.)\n\n");

    printf("                          M1 (cells, +23.5)                | M2 (pct)\n");
    printf("sigma*dt    \\   N=    %-7d %-7d %-7d %-7d | %-7d %-7d %-7d %-7d\n",
           Ns[0], Ns[1], Ns[2], Ns[3], Ns[0], Ns[1], Ns[2], Ns[3]);
    printf("---------------------------------------------------------------------------\n");

    int total = 0, nan_count = 0;
    for (int is = 0; is < n_sdt; is++) {
        float m1[4], m2[4];
        int   nan_row = 0;
        for (int iN = 0; iN < n_N; iN++) {
            const float sigma_max = sdts[is] / dt_ref;
            const int eff_N = (sdts[is] > 0.0f) ? Ns[iN] : 0;
            const gr_damp_config_t cfg = {
                .n_damping          = eff_N,
                .kind               = GR_DAMP_POLYNOMIAL,
                .poly_order         = 3.0f,
                .exp_beta           = 0.0f,
                .target_reflection  = 0.0f,
                .sigma_max_override = (sdts[is] > 0.0f) ? sigma_max : 0.0f,
                .time_form          = GR_DAMP_TIME_CRITICAL,
            };
            int nan_m1 = 0, nan_m2 = 0;
            m1[iN] = drift_off_center(&cfg, &nan_m1);
            m2[iN] = orbit_drift_pct (&cfg, &nan_m2);
            total++;
            nan_row |= (nan_m1 | nan_m2);
            nan_count += (nan_m1 | nan_m2);
        }
        printf("sigma*dt=%-8.4f| %7.3f %7.3f %7.3f %7.3f | %+6.2f%% %+6.2f%% %+6.2f%% %+6.2f%%%s\n",
               sdts[is],
               m1[0], m1[1], m1[2], m1[3],
               m2[0], m2[1], m2[2], m2[3],
               nan_row ? "  [NaN]" : "");
        fflush(stdout);
    }

    printf("\nconfigurations run: %d, NaN/Inf: %d\n", total, nan_count);
    printf("ALL CHECKS PASSED.\n");
    return 0;
}
