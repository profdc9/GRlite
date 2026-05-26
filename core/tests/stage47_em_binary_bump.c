/* Stage 47 -- EM binary PIC (both bodies active sources) with bump
 * kernel.  Parallel to Stage 45 (single test particle around analytic
 * background) but with both particles full PIC participants.  This is
 * the stress test for the self-force problem: at q_test = Q_central
 * (equal-coupling binary), the self-force is comparable to the mutual
 * force, unlike Stage 45's 1e-3 suppression.
 *
 * Setup mirrors Stage 32: two equal-mass (m=0.01), opposite-charge
 * (q = +/-0.01) particles in circular orbit, radius r=8.  Run 4
 * orbital periods.  Compare TSC bare, TSC + smoothing, BUMP R=6,
 * BUMP R=8.
 *
 * Analytic: v_orb = Q * sqrt(k_e / m) = 0.1.
 *           T_orb = 2 * pi * r / v_orb ~ 502 sim units (~711 dt). */

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

static void run_binary_orbit(int n_orbits,
                              gr_shape_function_t shape, float Rk, int rho_smooth,
                              orbit_result_t* out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 0.01f;
    const float Q      = 0.01f;
    const float r_orb  = 20.0f * dx;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float k_e    = 1.0f;
    const float v_orb  = Q * sqrtf(k_e / mass);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_orb;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { out->nan = 1; return; }
    /* Apply pedagogical defaults first; the scenario will set most of
     * its own configuration after, but our explicit overrides below
     * will set the shape/radius/smoothing per the comparison config. */
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);

    const float params[6] = {mass, r_orb, 1.0f, cx, cy, Q};
    if (gr_sim_load_scenario(sim, "pic_binary_em", params, 6) != 0) {
        out->nan = 1; gr_sim_destroy(sim); return;
    }
    /* The scenario set shape=TSC, force_interp=LB.  Override our
     * comparison config AFTER scenario load. */
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, rho_smooth);
    gr_sim_set_j_smooth_passes(sim, rho_smooth);
    /* PURE-EM test: disable gravity (default G_eff=1 would otherwise
     * give 2D-log gravitational attraction comparable to EM at m=Q). */
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

    for (int k = 0; k <= n_orbits; k++) out->r_at_orbit[k] = 0.0f;
    {
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        out->r_at_orbit[0] = sqrtf(rx * rx + ry * ry);
    }

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.4f * (float) n_orbits * T_ana / dt);
    float th_prev = atan2f(0.0f, -r_orb);   /* particle 0 starts at (cx-r, cy) */
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
        /* CW orbit: theta decreases, watch for the +pi -> -pi wrap. */
        if (th_prev < -0.9f * (float) M_PI && th > 0.9f * (float) M_PI) {
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
    printf("=== stage47_em_binary_bump (EM binary PIC, both bodies active) ===\n");
    printf("Equal-mass m=0.01, opposite charge q=+/-0.01.  Circular orbit r=20.\n");
    printf("v_orb = Q sqrt(k_e/m) = 0.1; T_orb ~ 1257 sim units.  4 orbits.\n");
    printf("r=20 chosen so 2r=40 cells > 2*R_max=17 (no kernel overlap).\n");
    printf("\nUnlike Stage 45 (q/Q=1e-3 test particle), here both bodies are\n");
    printf("PIC sources with EQUAL coupling -- self-force is comparable to\n");
    printf("mutual force.  Strongest test of the bump kernel's self-force\n");
    printf("mitigation.\n\n");

    const int N = 4;
    orbit_result_t r_tsc_bare;   r_tsc_bare.tag    = "TSC bare";
    orbit_result_t r_tsc_smooth; r_tsc_smooth.tag  = "TSC + N=16";
    orbit_result_t r_bump_r6;    r_bump_r6.tag     = "BUMP R=6";
    orbit_result_t r_bump_r8;    r_bump_r8.tag     = "BUMP R=8";

    run_binary_orbit(N, GR_SHAPE_TSC,  0.0f, 0,  &r_tsc_bare);
    run_binary_orbit(N, GR_SHAPE_TSC,  0.0f, 16, &r_tsc_smooth);
    run_binary_orbit(N, GR_SHAPE_BUMP, 6.0f, 0,  &r_bump_r6);
    run_binary_orbit(N, GR_SHAPE_BUMP, 8.0f, 0,  &r_bump_r8);

    printf("Radius (%% drift from r_orb=20.0) at each orbit:\n");
    printf("  %-7s %-15s %-15s %-15s %-15s\n",
           "orbit", "TSC bare", "TSC + N=16", "BUMP R=6", "BUMP R=8");
    printf("  ------------------------------------------------------------------------\n");
    orbit_result_t* arr[4] = { &r_tsc_bare, &r_tsc_smooth, &r_bump_r6, &r_bump_r8 };
    const float r_orb = 20.0f;
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
            printf("  %-15s: %+8.3f%% (completed %d/%d orbits)\n",
                   arr[i]->tag, (double) pct, arr[i]->n_completed, N);
        } else {
            printf("  %-15s: UNBOUND (only %d/%d orbits completed)\n",
                   arr[i]->tag, arr[i]->n_completed, N);
        }
    }
    return 0;
}
