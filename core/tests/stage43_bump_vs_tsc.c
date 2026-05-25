/* Stage 43 -- empirical comparison: TSC vs bump kernel for the
 * Stage 37 spurious-force test.
 *
 * Per gr_sandbox_v38 sec kernel_design (Tier 1): does the C-infinity
 * bump function exp(-1/(1-(d/1.5)^2)) on 3-cell support give a
 * smaller spurious force on a uniformly moving charge than the
 * standard TSC W_3 kernel at the same 3-cell support?
 *
 * Setup: same as Stage 37 / Stage 41 v-scan, mass=1.0 for clean
 * window measurement, rho_smooth=0.  For each v, run twice
 * (shape=TSC, shape=BUMP), report per-step acceleration and ratio.
 *
 * Also runs Stage 39-style reciprocity check at v=0 for both shapes
 * to confirm BUMP preserves the HE adjoint condition (should give
 * |F_1 + F_2| / |F_1| < 1e-6 at every separation, matching TSC). */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static float run_accel_test(float v_drift, gr_shape_function_t shape) {
    const int   W       = 128;
    const int   H       = 1024;
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 16;
    const float Q       = 0.01f;
    const float mass    = 1.0f;
    const float y_start = 16.0f + 32.0f;
    const float dt      = cfl * dx / c_eff;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return 0.0f;
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    /* PHI-ONLY config: inductive and magnetic both off, so only -q grad phi
     * acts.  This is the channel where the deposit/gather kernel matters. */
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 0);
    gr_sim_set_em_magnetic_enabled(sim, 0);

    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    const int n_warmup = 2000;
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    const gr_particle_t* p_post = gr_sim_get_particle(sim, 0);
    const float vy_post = vy_of(p_post, c_eff);

    const int force_window = 200;
    for (int s = 0; s < force_window; s++) gr_sim_step(sim);

    const gr_particle_t* p_after = gr_sim_get_particle(sim, 0);
    const float vy_end = vy_of(p_after, c_eff);
    const float accel = (vy_end - vy_post) / (float) force_window / dt;

    gr_sim_destroy(sim);
    return accel;
}

/* Static reciprocity check at v=0 with HEAVY mass.  Mirrors Stage 39. */
static float run_reciprocity(gr_shape_function_t shape, int sep) {
    const int   W       = 64;
    const int   H       = 64;
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 8;
    const float dt      = cfl * dx / c_eff;
    const float Q       = 0.01f;
    const float mass    = 1.0e6f;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;
    const float cy_base = ((float) (H - 1) * 0.5f) * dx;
    const float y1      = cy_base - 0.5f * (float) sep * dx;
    const float y2      = cy_base + 0.5f * (float) sep * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return -1.0f;
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);

    gr_sim_add_particle(sim, cx, y1, mass, +Q, 0.0f, 0.0f);
    gr_sim_add_particle(sim, cx, y2, mass, +Q, 0.0f, 0.0f);

    for (int s = 0; s < 500; s++) gr_sim_step(sim);
    const gr_particle_t* pp1 = gr_sim_get_particle(sim, 0);
    const gr_particle_t* pp2 = gr_sim_get_particle(sim, 1);
    const float py1_pre = pp1->py;
    const float py2_pre = pp2->py;
    for (int s = 0; s < 100; s++) gr_sim_step(sim);
    pp1 = gr_sim_get_particle(sim, 0);
    pp2 = gr_sim_get_particle(sim, 1);
    const float Fy1 = (pp1->py - py1_pre) / (100.0f * dt);
    const float Fy2 = (pp2->py - py2_pre) / (100.0f * dt);
    const float F1_mag  = fabsf(Fy1);
    const float sum_mag = fabsf(Fy1 + Fy2);
    const float ratio = (F1_mag > 0.0f) ? sum_mag / F1_mag : 0.0f;
    gr_sim_destroy(sim);
    return ratio;
}

int main(void) {
    printf("=== stage43_bump_vs_tsc (v38 Tier-1 empirical test) ===\n");
    printf("Compare bump kernel vs TSC on Stage 37 phi-only spurious force.\n");
    printf("Same 3-cell support; bump is C-infinity (sub-exp Fourier decay),\n");
    printf("TSC is C^1 (1/k^3 polynomial decay).  Mass=1.0, rho_smooth=0,\n");
    printf("inductive and magnetic gates OFF (phi-only isolation).\n\n");

    const float v_vals[] = { 0.001f, 0.01f, 0.05f, 0.1f, 0.2f, 0.3f };
    const int   nv       = (int)(sizeof(v_vals) / sizeof(v_vals[0]));

    printf("--- v-scan: phi-only spurious accel ---\n");
    printf("%-7s %-14s %-14s %-12s\n",
           "v/c", "accel(TSC)", "accel(BUMP)", "BUMP/TSC");
    printf("-----------------------------------------------\n");
    int wins = 0;
    int half = nv / 2;
    for (int i = 0; i < nv; i++) {
        const float aT = run_accel_test(v_vals[i], GR_SHAPE_TSC);
        const float aB = run_accel_test(v_vals[i], GR_SHAPE_BUMP);
        const float ratio = (fabsf(aT) > 0.0f) ? aB / aT : 0.0f;
        const float mag_ratio = (fabsf(aT) > 0.0f) ? fabsf(aB) / fabsf(aT) : 0.0f;
        printf("%-7.4f %+13.3e %+13.3e %+11.3f  |ratio|=%.3f\n",
               (double) v_vals[i], (double) aT, (double) aB,
               (double) ratio, (double) mag_ratio);
        if (mag_ratio < 0.5f) wins++;
    }

    printf("\n--- Stage 39-style static reciprocity (v=0) ---\n");
    printf("%-6s %-14s %-14s\n", "sep", "ratio(TSC)", "ratio(BUMP)");
    printf("------------------------------------\n");
    const int seps[] = { 4, 8, 16 };
    const int ns     = (int)(sizeof(seps) / sizeof(seps[0]));
    float worst_bump = 0.0f;
    for (int i = 0; i < ns; i++) {
        const float rT = run_reciprocity(GR_SHAPE_TSC, seps[i]);
        const float rB = run_reciprocity(GR_SHAPE_BUMP, seps[i]);
        printf("%-6d %12.3e   %12.3e\n", seps[i], (double) rT, (double) rB);
        if (rB > worst_bump) worst_bump = rB;
    }

    printf("\n");
    printf("Acceptance summary (v38 sec kernel_design Tier 1):\n");
    printf("  Threshold A: |BUMP/TSC| < 0.5 at >= %d of %d v values.\n",
           half, nv);
    printf("    Got %d wins.  %s\n", wins,
           (wins >= half) ? "*** KERNEL SHAPE MATTERS -- expand to Tier 2/3 ***"
                          : "Kernel shape does NOT help much at fixed support.");
    printf("  Threshold B: BUMP reciprocity ratio < 1e-6 at every sep.\n");
    printf("    Worst BUMP ratio: %.3e   %s\n", (double) worst_bump,
           (worst_bump < 1.0e-6f) ? "*** PASS ***"
                                  : "*** FAIL -- BUMP breaks HE adjoint ***");

    return 0;
}
