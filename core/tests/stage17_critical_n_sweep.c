/* Stage 17 — push N higher under critical+low-CFL and confirm M2 unchanged.
 *
 * Stage 16 found that critical damping at CFL = 0.65, m=3, N=48 gives
 * M1 = ~2.25 cells essentially independent of sigma*dt across three
 * decades.  Now: how much further does M1 drop as we widen the ring?
 * And does Phase E orbit drift (M2) — which the earlier sweep showed
 * is boundary-insensitive at small N=16/32/48 — remain unchanged at
 * larger N?  We're checking that the new config doesn't introduce a
 * regression on M2.
 *
 * Config: critical + CFL=0.65 + m=3 + sigma*dt = 0.01 (the
 * "sweet-spot" point from stage16 where critical visibly beats
 * multiplicative).  Sweep N in {32, 48, 64, 80}. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define INTERIOR 256

static float drift_off_center(int n_damp, float cfl, float sigma_dt,
                              gr_damp_time_form_t form, float target_time,
                              int* nan_out, int* n_steps_out) {
    const int   W      = INTERIOR + 2 * n_damp;
    const int   H      = INTERIOR + 2 * n_damp;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float mass   = 1.0e-3f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float dx_off = 23.5f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { *nan_out = 1; return 0.0f; }

    const float dt = gr_sim_dt(sim);
    const gr_damp_config_t cfg = {
        .n_damping          = (sigma_dt > 0.0f) ? n_damp : 0,
        .kind               = GR_DAMP_POLYNOMIAL,
        .poly_order         = 3.0f,
        .exp_beta           = 0.0f,
        .target_reflection  = 0.0f,
        .sigma_max_override = (sigma_dt > 0.0f) ? (sigma_dt / dt) : 0.0f,
        .time_form          = form,
    };
    gr_sim_set_damping_config(sim, &cfg);
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

static float orbit_drift_pct(int n_damp, float cfl, float sigma_dt,
                             gr_damp_time_form_t form, int* nan_out) {
    const int   W      = INTERIOR + 2 * n_damp;
    const int   H      = INTERIOR + 2 * n_damp;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
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
    const float dt = gr_sim_dt(sim);
    const gr_damp_config_t cfg = {
        .n_damping          = (sigma_dt > 0.0f) ? n_damp : 0,
        .kind               = GR_DAMP_POLYNOMIAL,
        .poly_order         = 3.0f,
        .exp_beta           = 0.0f,
        .target_reflection  = 0.0f,
        .sigma_max_override = (sigma_dt > 0.0f) ? (sigma_dt / dt) : 0.0f,
        .time_form          = form,
    };
    gr_sim_set_damping_config(sim, &cfg);
    if (gr_sim_load_scenario(sim, "pic_orbiting", par, 4) != 0) {
        gr_sim_destroy(sim);
        *nan_out = 1;
        return 0.0f;
    }
    const int N_steps = (int) (T_ana / dt);
    gr_sim_step_n(sim, N_steps);
    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float r = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
    *nan_out = !isfinite(r);
    gr_sim_destroy(sim);
    return 100.0f * (r - r_orb) / r_orb;
}

int main(void) {
    printf("=== stage17_critical_n_sweep ===\n");
    printf("Config: m=3, sigma*dt=0.01, time_form=CRITICAL, CFL=0.65\n");
    printf("(compare: legacy default = m=2, sigma*dt=0.46, MULTIPLICATIVE, CFL=0.7071, N=16)\n");
    printf("\n");

    const float cfl_low  = 0.65f;
    const float cfl_ref  = 1.0f / sqrtf(2.0f);
    const float target_t = 1414.0f;
    const float sdt      = 0.01f;

    const int Ns[] = {32, 48, 64, 80};
    const int n_N  = sizeof(Ns) / sizeof(Ns[0]);

    printf("%-6s %-22s %-22s\n", "N_d", "M1 (cells, +23.5)", "M2 (orbit drift %)");
    printf("%-6s %-10s %-10s %-10s %-10s\n",
           "", "crit", "mult_ref", "crit", "mult_ref");
    printf("------------------------------------------------------------\n");

    /* Legacy reference at N=16 (the current production setting): */
    int nan; int n_steps;
    const float m1_legacy = drift_off_center(16, cfl_ref, 0.46f,
                                             GR_DAMP_TIME_MULTIPLICATIVE,
                                             target_t, &nan, &n_steps);
    const float m2_legacy = orbit_drift_pct (16, cfl_ref, 0.46f,
                                             GR_DAMP_TIME_MULTIPLICATIVE,
                                             &nan);
    printf("legacy(N=16, m=2, sigma*dt=0.46, mult, CFL=0.7071):\n");
    printf("  M1 = %.3f cells   M2 = %+.2f%%\n\n", m1_legacy, m2_legacy);

    for (int iN = 0; iN < n_N; iN++) {
        int nan_m1c, nan_m1m, nan_m2c, nan_m2m, ns;
        const float m1_crit = drift_off_center(Ns[iN], cfl_low, sdt,
                                               GR_DAMP_TIME_CRITICAL,
                                               target_t, &nan_m1c, &ns);
        const float m1_mult = drift_off_center(Ns[iN], cfl_ref, sdt,
                                               GR_DAMP_TIME_MULTIPLICATIVE,
                                               target_t, &nan_m1m, &ns);
        const float m2_crit = orbit_drift_pct (Ns[iN], cfl_low, sdt,
                                               GR_DAMP_TIME_CRITICAL,
                                               &nan_m2c);
        const float m2_mult = orbit_drift_pct (Ns[iN], cfl_ref, sdt,
                                               GR_DAMP_TIME_MULTIPLICATIVE,
                                               &nan_m2m);
        printf("%-6d %-10.3f %-10.3f %+9.2f%% %+9.2f%%%s\n",
               Ns[iN], m1_crit, m1_mult, m2_crit, m2_mult,
               (nan_m1c | nan_m1m | nan_m2c | nan_m2m) ? "  [NaN]" : "");
        fflush(stdout);
    }

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
