/* Stage 43 -- bump kernel R scan vs TSC + smoothing, with EQUAL
 * DISTANCE TRAVELED methodology.
 *
 * Prior version used a fixed 200-step force window regardless of v,
 * which meant low-v measurements were dominated by a single (nearly
 * fixed) sub-cell offset rather than averaging over the full
 * u-dependence.  Result: spurious v=0.005 "spike."
 *
 * This version measures force over a fixed DISTANCE of 5 cells
 * traveled, with measurement window steps scaled as
 * n_steps = ceil(5 / (v * dt)).  Heavier mass (m=100) keeps the
 * spurious force from materially changing v over the longer windows.
 *
 * Reported quantity: per-unit-time spurious acceleration
 *   accel = (vy_end - vy_start) / total_time
 * where total_time = n_steps * dt.  This is the LOCAL force on a
 * moving deposit averaged over all sub-cell offsets it traverses. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum { CFG_TSC, CFG_BUMP } cfg_kind_t;

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static float run_accel(float v_drift, cfg_kind_t kind, float Rk, int rho_smooth) {
    const int   W       = 128;
    const int   H       = 2048;       /* longer y-axis: 5-cell traversal at high v stays safe */
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 16;
    const float Q       = 0.01f;
    const float mass    = 100.0f;     /* heavy: force barely changes v over the window */
    const float y_start = 16.0f + 32.0f;
    const float dt      = cfl * dx / c_eff;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;

    /* Force-measurement window: ceil(5 / (v * dt)) steps so the
     * particle moves 5 cells regardless of v.  Cap at 200000 to keep
     * runtime bounded; at v=0.001 the cap fires (would otherwise be
     * 7143 steps -- well under the cap).  Min 200. */
    int n_window = (int) ceilf(5.0f / (v_drift * dt));
    if (n_window < 200)   n_window = 200;
    if (n_window > 200000) n_window = 200000;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return 0.0f;
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    if (kind == CFG_BUMP) {
        gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
        gr_sim_set_kernel_radius(sim, Rk);
    } else {
        gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    }
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, rho_smooth);
    gr_sim_set_j_smooth_passes(sim, rho_smooth);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 0);
    gr_sim_set_em_magnetic_enabled(sim, 0);

    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    /* Warmup: 2000 steps (1414 cells of light-travel; wake settles). */
    for (int s = 0; s < 2000; s++) gr_sim_step(sim);
    const gr_particle_t* p_post = gr_sim_get_particle(sim, 0);
    const float vy_post = vy_of(p_post, c_eff);

    /* Measurement window: scaled with 1/v so particle traverses ~5 cells. */
    for (int s = 0; s < n_window; s++) gr_sim_step(sim);
    const gr_particle_t* p_after = gr_sim_get_particle(sim, 0);
    const float vy_end = vy_of(p_after, c_eff);
    const float total_time = (float) n_window * dt;
    const float accel = (vy_end - vy_post) / total_time;

    gr_sim_destroy(sim);
    return accel;
}

int main(void) {
    printf("=== stage43_bump_vs_tsc (R scan + smoothing, equal-distance methodology) ===\n");
    printf("Measurement: 5 cells of trajectory at every v (variable step count).\n");
    printf("Heavy particle (m=100) so spurious force doesn't materially change v.\n");
    printf("Window samples the full sub-cell-offset distribution at every v.\n\n");

    const float v_vals[] = { 0.001f, 0.005f, 0.01f, 0.05f, 0.1f, 0.2f, 0.3f, 0.5f };
    const int   nv       = (int)(sizeof(v_vals) / sizeof(v_vals[0]));

    /* TSC baseline. */
    float aT[8];
    printf("--- TSC baseline (no smoothing) ---\n");
    printf("%-7s %-13s %-10s\n", "v/c", "accel (cells/dt^2)", "n_window");
    printf("------------------------------------------\n");
    for (int i = 0; i < nv; i++) {
        const int nw = (int)(5.0f / (v_vals[i] * (1.0f/sqrtf(2.0f))));
        aT[i] = run_accel(v_vals[i], CFG_TSC, 0.0f, 0);
        printf("%-7.4f %+12.3e   %d\n", (double) v_vals[i], (double) aT[i],
               (nw < 200) ? 200 : (nw > 200000 ? 200000 : nw));
    }
    printf("\n");

    /* Bump R scan. */
    const float Rs[] = { 2.5f, 4.0f, 6.0f, 8.0f };
    const int   nR   = 4;
    printf("--- bump kernel: |accel(BUMP_R) / accel(TSC)| ---\n");
    printf("%-7s", "v/c");
    for (int j = 0; j < nR; j++) printf("  R=%-3.1f       ", (double) Rs[j]);
    printf("\n");
    for (int j = 0; j < 7 + nR * 14; j++) printf("-");
    printf("\n");
    double bump_ratios[nR][8];
    for (int i = 0; i < nv; i++) {
        printf("%-7.4f", (double) v_vals[i]);
        for (int j = 0; j < nR; j++) {
            const float a = run_accel(v_vals[i], CFG_BUMP, Rs[j], 0);
            const double r = (fabs(aT[i]) > 0.0) ? fabs(a) / fabs(aT[i]) : 0.0;
            bump_ratios[j][i] = r;
            printf("  %11.4f  ", r);
        }
        printf("\n");
    }
    printf("\n");

    /* Smoothing scan. */
    const int passes[] = { 4, 16, 64 };
    const int nP       = 3;
    printf("--- TSC + smoothing: |accel(TSC,N passes) / accel(TSC,0)| ---\n");
    printf("%-7s", "v/c");
    for (int j = 0; j < nP; j++) printf("  N=%-4d       ", passes[j]);
    printf("\n");
    for (int j = 0; j < 7 + nP * 14; j++) printf("-");
    printf("\n");
    double smooth_ratios[nP][8];
    for (int i = 0; i < nv; i++) {
        printf("%-7.4f", (double) v_vals[i]);
        for (int j = 0; j < nP; j++) {
            const float a = run_accel(v_vals[i], CFG_TSC, 0.0f, passes[j]);
            const double r = (fabs(aT[i]) > 0.0) ? fabs(a) / fabs(aT[i]) : 0.0;
            smooth_ratios[j][i] = r;
            printf("  %11.4f  ", r);
        }
        printf("\n");
    }
    printf("\n");

    /* Stability metric. */
    printf("--- stability comparison ---\n");
    printf("%-12s %-12s %-12s %-12s %-12s\n",
           "config", "avg ratio", "min ratio", "max ratio", "max/min");
    printf("--------------------------------------------------------------\n");
    for (int j = 0; j < nR; j++) {
        double s = 0.0, mn = 1e30, mx = 0.0;
        for (int i = 0; i < nv; i++) {
            s += bump_ratios[j][i];
            if (bump_ratios[j][i] < mn) mn = bump_ratios[j][i];
            if (bump_ratios[j][i] > mx) mx = bump_ratios[j][i];
        }
        printf("bump R=%-4.1f  %-12.4f %-12.4f %-12.4f %-12.2f\n",
               (double) Rs[j], s / nv, mn, mx, mx / mn);
    }
    for (int j = 0; j < nP; j++) {
        double s = 0.0, mn = 1e30, mx = 0.0;
        for (int i = 0; i < nv; i++) {
            s += smooth_ratios[j][i];
            if (smooth_ratios[j][i] < mn) mn = smooth_ratios[j][i];
            if (smooth_ratios[j][i] > mx) mx = smooth_ratios[j][i];
        }
        printf("TSC N=%-5d  %-12.4f %-12.4f %-12.4f %-12.2f\n",
               passes[j], s / nv, mn, mx, mx / mn);
    }

    return 0;
}
