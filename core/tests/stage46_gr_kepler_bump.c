/* Stage 46 -- GR Kepler orbit with bump kernel, GRAVITY analog of
 * Stage 45 (EM Kepler).  Same orbital setup, same shape/smoothing
 * comparison, but gravity sector instead of EM.
 *
 * Setup: analytic background point mass (G_eff = 1, M_central = 1,
 * eps=1), test particle m_test = 1e-3 at r=20.  4 orbits.  Compare
 * TSC bare, TSC + smoothing N=16, BUMP R=6, BUMP R=8.
 *
 * Gravitomagnetic inductive defaults OFF (the existing GRlite
 * convention).  Gravitomagnetic v x B_g force defaults ON.  Same
 * shape_function dispatch path as EM (handled by sim.c rho deposit
 * for rho_matter and by the gravity-sector reads in particle.c). */

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

static void run_orbit(float GM, float m_test, float r_orb, int n_orbits,
                       gr_shape_function_t shape, float Rk, int rho_smooth,
                       orbit_result_t* out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;

    /* Newtonian circular velocity for analytic softened point mass. */
    const float g_mag  = GM * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f);
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
    /* Gravity sector on; EM sector irrelevant (no charge on the
     * test particle). */
    gr_sim_set_G_eff(sim, 1.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 1);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);  /* default convention */
    gr_sim_set_background_point_mass(sim, cx, cy, GM, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    /* Shape function and smoothing per the test condition. */
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, rho_smooth);
    gr_sim_set_j_smooth_passes(sim, rho_smooth);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, /*charge=*/0.0f,
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
    printf("=== stage46_gr_kepler_bump (GR analog of Stage 45) ===\n");
    printf("Gravity Kepler with analytic point-mass background.\n");
    printf("GM=+1, m_test=+1e-3, r=20, 4 orbits.  Gravitomagnetic inductive OFF.\n\n");

    const float GM     = +1.0f;
    const float m_test = +1.0e-3f;
    const float r_orb  = 20.0f;
    const int   N      = 4;

    orbit_result_t r_tsc_bare;   r_tsc_bare.tag    = "TSC bare";
    orbit_result_t r_tsc_smooth; r_tsc_smooth.tag  = "TSC + N=16";
    orbit_result_t r_bump_r6;    r_bump_r6.tag     = "BUMP R=6";
    orbit_result_t r_bump_r8;    r_bump_r8.tag     = "BUMP R=8";

    run_orbit(GM, m_test, r_orb, N, GR_SHAPE_TSC,  0.0f, 0,  &r_tsc_bare);
    run_orbit(GM, m_test, r_orb, N, GR_SHAPE_TSC,  0.0f, 16, &r_tsc_smooth);
    run_orbit(GM, m_test, r_orb, N, GR_SHAPE_BUMP, 6.0f, 0,  &r_bump_r6);
    run_orbit(GM, m_test, r_orb, N, GR_SHAPE_BUMP, 8.0f, 0,  &r_bump_r8);

    printf("Radius (%% drift from r_orb) at each pi-wrap:\n");
    printf("  %-7s %-15s %-15s %-15s %-15s\n",
           "orbit", "TSC bare", "TSC + N=16", "BUMP R=6", "BUMP R=8");
    printf("  ------------------------------------------------------------------------\n");
    orbit_result_t* arr[4] = { &r_tsc_bare, &r_tsc_smooth, &r_bump_r6, &r_bump_r8 };
    for (int k = 1; k <= N; k++) {
        char c[4][24];
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
    for (int i = 0; i < 4; i++) {
        if (arr[i]->n_completed == N) {
            const float pct = 100.0f * (arr[i]->r_at_orbit[N] - r_orb) / r_orb;
            printf("  %-15s: %+8.3f%%\n", arr[i]->tag, (double) pct);
        } else {
            printf("  %-15s: UNBOUND (%d/%d orbits)\n",
                   arr[i]->tag, arr[i]->n_completed, N);
        }
    }
    return 0;
}
