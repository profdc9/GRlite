/* Stage 45 -- EM PIC Kepler orbit with bump kernel: the downstream
 * verification that ties Stage 37 / 43 / 44 together.
 *
 * Stage 27 reported +25% radius drift over 4 orbits with default
 * settings (CIC + LEGACY, inductive +1).  That was the spurious-
 * heating problem that started the v36/v37/v38 investigation arc.
 *
 * This test runs the same orbital setup four ways:
 *   (1) TSC + rho_smooth=0 (bare; expect heating)
 *   (2) TSC + rho_smooth=16 (smoothing mitigation -- the existing fix)
 *   (3) BUMP + R=6 + rho_smooth=0 (Tier-2 bump kernel mitigation)
 *   (4) BUMP + R=8 + rho_smooth=0 (wider bump)
 *
 * Setup mirrors Stage 27: Q=+1 central charge (analytic background),
 * q_test=-1e-3 orbiting, m_test=+1e-3, r_orb=20.  4 orbits.  Report
 * radius at each pi-wrap. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float r_at_orbit[8];
    int   n_completed;
    int   nan;
    const char* tag;
} orbit_result_t;

static void run_orbit(float Q, float q_test, float m_test, float r_orb,
                       int n_orbits,
                       gr_shape_function_t shape, float Rk, int rho_smooth,
                       orbit_result_t* out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float k_e    = 1.0f;

    const float g_mag  = fabsf(q_test * Q) * k_e * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f) / m_test;
    const float rg     = r_orb * g_mag;
    const float rg2_c2 = rg * rg / (c_eff * c_eff);
    const float u_v2   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ = sqrtf(u_v2);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_circ;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { out->nan = 1; return; }
    gr_sim_set_damping(sim, 16);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* Full force: phi + inductive + magnetic all on (the configuration
     * that gave Stage 27's +25% outspiral). */
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, +1.0f);
    gr_sim_set_em_inductive_disc(sim, GR_INDUCTIVE_CENTERED);
    /* No J time-correction (raw Esirkepov). */
    gr_sim_set_j_time_correction_enabled(sim, 0);
    /* Analytic background point charge. */
    gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    /* Shape function and smoothing per the test condition. */
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, rho_smooth);
    gr_sim_set_j_smooth_passes(sim, rho_smooth);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test,
                        0.0f, v_circ);

    for (int k = 0; k <= n_orbits; k++) out->r_at_orbit[k] = 0.0f;
    {
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        out->r_at_orbit[0] = sqrtf(rx * rx + ry * ry);
    }

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.2f * (float) n_orbits * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    for (int s = 0; s < n_max && wraps < n_orbits; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) { out->nan = 1; gr_sim_destroy(sim); return; }
        if (r_now > 0.9f * (W * 0.5f)) {
            /* Particle escaped the safe interior. */
            out->n_completed = wraps;
            gr_sim_destroy(sim);
            return;
        }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            if (wraps <= n_orbits) {
                out->r_at_orbit[wraps] = r_now;
                out->n_completed = wraps;
            }
        }
        th_prev = th;
    }
    out->nan = 0;
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage45_em_kepler_bump (v38 downstream verification) ===\n");
    printf("Stage 27's +25%% outspiral: does the bump kernel fix it?\n");
    printf("Q=+1, q_test=-1e-3, m_test=+1e-3, r=20, 4 orbits.  Full Lorentz.\n\n");

    const float Q       = +1.0f;
    const float q_test  = -1.0e-3f;
    const float m_test  = +1.0e-3f;
    const float r_orb   = 20.0f;
    const int   N       = 4;

    orbit_result_t r_tsc_bare;   r_tsc_bare.tag    = "TSC bare";
    orbit_result_t r_tsc_smooth; r_tsc_smooth.tag  = "TSC + N=16";
    orbit_result_t r_bump_r6;    r_bump_r6.tag     = "BUMP R=6";
    orbit_result_t r_bump_r8;    r_bump_r8.tag     = "BUMP R=8";

    run_orbit(Q, q_test, m_test, r_orb, N, GR_SHAPE_TSC,  0.0f, 0,  &r_tsc_bare);
    run_orbit(Q, q_test, m_test, r_orb, N, GR_SHAPE_TSC,  0.0f, 16, &r_tsc_smooth);
    run_orbit(Q, q_test, m_test, r_orb, N, GR_SHAPE_BUMP, 6.0f, 0,  &r_bump_r6);
    run_orbit(Q, q_test, m_test, r_orb, N, GR_SHAPE_BUMP, 8.0f, 0,  &r_bump_r8);

    printf("Radius (%% drift from r_orb) at each pi-wrap:\n");
    printf("  %-7s %-15s %-15s %-15s %-15s\n",
           "orbit", "TSC bare", "TSC + N=16", "BUMP R=6", "BUMP R=8");
    printf("  ------------------------------------------------------------------------\n");
    for (int k = 1; k <= N; k++) {
        char c[4][24];
        orbit_result_t* arr[4] = { &r_tsc_bare, &r_tsc_smooth, &r_bump_r6, &r_bump_r8 };
        for (int i = 0; i < 4; i++) {
            if (k > arr[i]->n_completed) snprintf(c[i], 24, "<unbound>");
            else {
                const float pct = 100.0f * (arr[i]->r_at_orbit[k] - r_orb) / r_orb;
                snprintf(c[i], 24, "%+8.3f%%", (double) pct);
            }
        }
        printf("  %-7d %-15s %-15s %-15s %-15s\n", k, c[0], c[1], c[2], c[3]);
    }

    printf("\nFinal drift after %d orbits:\n", N);
    orbit_result_t* arr[4] = { &r_tsc_bare, &r_tsc_smooth, &r_bump_r6, &r_bump_r8 };
    for (int i = 0; i < 4; i++) {
        if (arr[i]->n_completed == N) {
            const float pct = 100.0f * (arr[i]->r_at_orbit[N] - r_orb) / r_orb;
            printf("  %-15s: %+8.3f%%\n", arr[i]->tag, (double) pct);
        } else {
            printf("  %-15s: UNBOUND (%d/%d orbits completed)\n",
                   arr[i]->tag, arr[i]->n_completed, N);
        }
    }

    printf("\nAcceptance:\n");
    printf("  TSC bare: expected unbound or +25%% (Stage 27 baseline).\n");
    printf("  TSC + N=16, BUMP R=6, BUMP R=8: expected bound with small drift\n");
    printf("    (radiation-reaction-like inspiral, or stable within ~1%%/orbit).\n");

    return 0;
}
