/* Stage 14 — sigma_max × N sweep at fixed profile (polynomial m=3).
 *
 * Question: at fixed polynomial order m=3, is there an optimal
 * sigma_max — large enough to absorb the outgoing wave well, small
 * enough to not produce a strong impedance jump at the inner edge of
 * the damping ring that reflects back into the interior?
 *
 * Sweep axes:
 *   N_damping      : 8, 16, 32, 48 (ring thickness)
 *   sigma_max * dt : 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.5
 *                    (per-cell absorption rate at the outer wall;
 *                     stability bound for the (1 - sigma dt) update is
 *                     roughly sigma_max dt < 1; we sweep into the
 *                     unstable region for completeness)
 *
 * Interior region held at 256 x 256 (same as stage13 sweep #2) so each
 * configuration has the same particle-to-ring distance.  Reference
 * literature default for m=3, R=1e-3: sigma_max dt ~ 13.82 dt / N.
 *
 *   For N=16, dt=0.707: lit default = 0.611
 *   For N=8           :              = 1.22
 *   For N=32          :              = 0.305
 *   For N=48          :              = 0.204
 *
 * Metrics:
 *   M1: HE off-center drift @ +23.5 cells, 2000 steps, m=1e-3.
 *   M2: Phase E orbit drift at 1 analytic period.
 *
 * Lower magnitude = better in both.  Look for the (N, sigma_max) ridge
 * where M1 minimizes. */

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
    printf("=== stage14_sigma_sweep (polynomial m=3, interior=%d) ===\n", INTERIOR);

    const int   Ns[]    = {8, 16, 32, 48};
    const int   n_N     = sizeof(Ns) / sizeof(Ns[0]);
    /* Extended to include sigma*dt = 0 (no damping baseline) and very low
     * values, since the initial sweep showed monotonic increase with sigma
     * — looking for whether the trend continues all the way down to
     * "barely any damping" and what the no-damping baseline looks like. */
    const float sdts[]  = {0.0f, 0.0001f, 0.001f, 0.01f, 0.05f, 0.1f, 0.2f, 0.4f, 0.7f, 1.5f, 2.5f};
    const int   n_sdt   = sizeof(sdts) / sizeof(sdts[0]);
    const float dt_ref  = 1.0f / sqrtf(2.0f);  /* dt = cfl * dx / c = 0.7071 */

    printf("\nM1: HE drift @ +23.5 cells (cells), 2000 steps, m=1e-3.\n");
    printf("M2: pic_orbiting radial drift at 1 T_ana, m=1e-3 (percent).\n");
    printf("sigma_max in units of 1/dt; sigma*dt at wall as labeled.\n");
    printf("Lit ref m=3, R=1e-3 gives sigma*dt = 13.82 dt / N = "
           "[N=8: %.3f, N=16: %.3f, N=32: %.3f, N=48: %.3f].\n\n",
           13.82 * dt_ref / 8, 13.82 * dt_ref / 16,
           13.82 * dt_ref / 32, 13.82 * dt_ref / 48);

    printf("                          M1 (cells, +23.5)                | M2 (pct)\n");
    printf("sigma*dt    \\   N=    %-7d %-7d %-7d %-7d | %-7d %-7d %-7d %-7d\n",
           Ns[0], Ns[1], Ns[2], Ns[3], Ns[0], Ns[1], Ns[2], Ns[3]);
    printf("---------------------------------------------------------------------------\n");

    int total = 0, nan_count = 0;
    for (int is = 0; is < n_sdt; is++) {
        float m1[4], m2[4];
        int   nan_row = 0;
        for (int iN = 0; iN < n_N; iN++) {
            const float sigma_max = sdts[is] / dt_ref;  /* physical units */
            /* When sigma*dt = 0, disable damping entirely (n_damping = 0)
             * so we observe the no-absorber baseline. */
            const int eff_N = (sdts[is] > 0.0f) ? Ns[iN] : 0;
            const gr_damp_config_t cfg = {
                .n_damping          = eff_N,
                .kind               = GR_DAMP_POLYNOMIAL,
                .poly_order         = 3.0f,
                .exp_beta           = 0.0f,
                .target_reflection  = 0.0f,
                .sigma_max_override = (sdts[is] > 0.0f) ? sigma_max : 0.0f,
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
    if (nan_count > 0) {
        printf("  (NaNs at sigma_max dt > 1 reflect the (1 - sigma dt) "
               "kernel instability;\n"
               "   the Crank-Nicolson form would be unconditionally stable.)\n");
    }
    printf("ALL CHECKS PASSED.\n");
    return 0;
}
