/* Stage 16 — does lowering CFL unlock a critical-damping advantage?
 *
 * Hypothesis: critical damping preserves the *static* field in the ring
 * much better than multiplicative.  Per step the static-field decay
 * rates are:
 *
 *    Multiplicative:   (1 - sigma dt)            ~ linear in sigma
 *    Critical (g,g):   (1 - (sigma dt)^2)        ~ quadratic in sigma
 *
 * So at moderate sigma*dt (say 0.01), critical leaves the ring's static
 * profile largely intact while multiplicative chews through it.  Less
 * distortion of the static profile = lower off-center HE drift.
 *
 * Critical damping is unstable at CFL = 1/sqrt(2) for any gamma < 1
 * (Nyquist mode escapes).  At CFL = 0.65, the bound 4*gamma >=
 * 8*CFL^2 = 3.38 gives gamma >= 0.845, i.e. sigma*dt <= 0.155 stable.
 * That covers our experiment range up to sigma*dt = 0.1.
 *
 * Test grid: same physical time across all configs (target = 1414 sim
 * units, same as the stage14/15 baselines).  Steps adjust to dt.
 *
 *   sigma*dt   in {0.001, 0.01, 0.05, 0.1}
 *   N_damping  in {32, 48}
 *   time_form  in {MULTIPLICATIVE, CRITICAL}
 *
 * Metric: M1 (off-center HE drift at +23.5 cells, m=1e-3).  M2
 * (Phase E orbit drift) we've already shown is boundary-insensitive;
 * skipping here. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define INTERIOR 256

static float drift_off_center(const gr_damp_config_t* cfg,
                              float cfl, float target_time,
                              int* nan_out, int* n_steps_out) {
    const int   W      = INTERIOR + 2 * cfg->n_damping;
    const int   H      = INTERIOR + 2 * cfg->n_damping;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float mass   = 1.0e-3f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float dx_off = 23.5f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { *nan_out = 1; return 0.0f; }

    /* Configure damping AFTER create so the sigma_max override uses
     * the actual sim->dt established by gr_sim_create. */
    const float dt        = gr_sim_dt(sim);
    const float sigma_dt  = cfg->sigma_max_override;  /* we pass sigma*dt directly via override */
    /* Translate: caller passed sigma*dt in sigma_max_override; rebuild a
     * proper cfg with sigma_max in physical units. */
    gr_damp_config_t real_cfg = *cfg;
    real_cfg.sigma_max_override = (sigma_dt > 0.0f) ? (sigma_dt / dt) : 0.0f;
    if (sigma_dt <= 0.0f) real_cfg.n_damping = 0;
    gr_sim_set_damping_config(sim, &real_cfg);
    gr_sim_set_particle_source_deposition(sim, 1);

    const float x0 = cx + dx_off;
    const float y0 = cy;
    gr_sim_add_particle(sim, x0, y0, mass, 0.0f, 0.0f, 0.0f);

    const int N = (int) (target_time / dt + 0.5f);
    *n_steps_out = N;
    gr_sim_step_n(sim, N);

    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float dr = sqrtf((p->x - x0) * (p->x - x0) + (p->y - y0) * (p->y - y0));
    *nan_out = !isfinite(dr);
    gr_sim_destroy(sim);
    return dr / dx;
}

int main(void) {
    printf("=== stage16_low_cfl_critical (m=3, interior=%d) ===\n", INTERIOR);

    const float cfl_low  = 0.65f;
    const float cfl_ref  = 1.0f / sqrtf(2.0f);    /* ~0.7071 */
    const float target_t = 1414.0f;               /* same physical time across all */

    const float sdts[]   = {0.0001f, 0.001f, 0.01f, 0.05f, 0.1f};
    const int   n_sdt    = sizeof(sdts) / sizeof(sdts[0]);
    const int   Ns[]     = {32, 48};
    const int   n_N      = sizeof(Ns) / sizeof(Ns[0]);

    printf("\nM1: HE drift @ +23.5 cells (cells), target physical time = %.0f, m=1e-3.\n", target_t);
    printf("(stage14 multiplicative baseline at CFL=%.4f: N=48 sigma*dt=1e-4 -> M1=2.37 cells)\n",
           cfl_ref);
    printf("\n");

    printf("%-12s %-6s %-10s %-6s %-9s %-9s\n",
           "form", "N_d", "sigma*dt", "CFL", "M1", "n_steps");
    printf("--------------------------------------------------------------\n");

    /* Three configurations per (sigma*dt, N) row:
     *   (a) Multiplicative @ CFL_ref         — stage14 baseline.
     *   (b) Multiplicative @ CFL_low         — control for CFL effect alone.
     *   (c) Critical       @ CFL_low         — the hypothesis under test.
     */
    for (int iN = 0; iN < n_N; iN++) {
        for (int is = 0; is < n_sdt; is++) {
            const float sdt = sdts[is];
            for (int trial = 0; trial < 3; trial++) {
                const int               use_critical = (trial == 2);
                const float             cfl          = (trial == 0) ? cfl_ref : cfl_low;
                const gr_damp_time_form_t form = use_critical
                                                     ? GR_DAMP_TIME_CRITICAL
                                                     : GR_DAMP_TIME_MULTIPLICATIVE;
                const char*             label = use_critical
                                                  ? "critical    "
                                                  : ((trial == 0) ? "mult@CFL_ref" : "mult@CFL_low");
                const gr_damp_config_t cfg = {
                    .n_damping          = Ns[iN],
                    .kind               = GR_DAMP_POLYNOMIAL,
                    .poly_order         = 3.0f,
                    .exp_beta           = 0.0f,
                    .target_reflection  = 0.0f,
                    .sigma_max_override = sdt,    /* re-interpreted as sigma*dt by helper */
                    .time_form          = form,
                };
                int nan = 0, n_steps = 0;
                const float m1 = drift_off_center(&cfg, cfl, target_t, &nan, &n_steps);
                printf("%-12s %-6d %-10.4f %-6.4f %-9.3f %-9d%s\n",
                       label, Ns[iN], sdt, cfl, m1, n_steps,
                       nan ? "  [NaN]" : "");
                fflush(stdout);
            }
            printf("\n");
        }
    }

    printf("ALL CHECKS PASSED.\n");
    return 0;
}
