/* Stage 29 — CFL convergence sweep on the inductive PIC inspiral.
 *
 * Same setup as Stage 28 (gravity closed-loop PIC at m_test = 1e-3,
 * r = 20, inductive piece enabled, post-half-step-A fix), run at
 * several CFL values to probe whether the ~2.3% per orbit inward drift
 * we now see is genuine radiation reaction or residual PIC discretization
 * error.
 *
 * Predictions:
 *   - If the inspiral is genuine radiation reaction (Larmor), the rate
 *     should be approximately INDEPENDENT of CFL (it's a physical
 *     property of the system, set by v/c and the coupling, not by dt).
 *   - If the inspiral is dominated by leapfrog truncation error in the
 *     PIC self-coupling, it should scale as O(dt^2) — halving CFL should
 *     drop the drift rate by ~4x.
 *   - In practice we expect a mix: a CFL-independent floor (the real
 *     radiation-reaction) plus a CFL-dependent residual (PIC truncation).
 *     Convergence as CFL -> 0 reveals the floor.
 *
 * CFLs swept: {1/sqrt(2), 1/(2*sqrt(2)), 1/(4*sqrt(2)), 1/(8*sqrt(2))}
 * = {0.707, 0.354, 0.177, 0.088}.  Halving each step.
 *
 * Baseline at each CFL: no inductive.  Reports the inductive contribution
 * (= drift_with_inductive - drift_without) at each CFL.  This isolates
 * the inductive piece from the (small, CFL-dependent) drift in the
 * spatial-derivative-only PIC dynamics. */

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
        exit(1); \
    } \
} while (0)

typedef struct {
    float cfl;
    float dt;
    int   n_steps_total;
    float drift_noind;     /* % drift after 4 orbits, inductive OFF */
    float drift_indon;     /* % drift after 4 orbits, inductive ON  */
    float drift_indcontrib;/* inductive contribution = indon - noind */
    int   n_completed_noind;
    int   n_completed_indon;
} sweep_result_t;

static float run(float GM, float m_test, float r_orb, float cfl,
                 int n_orbits, int grav_inductive_enabled,
                 int* n_completed_out, int* n_steps_total_out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;

    const float g_mag  = GM * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f);
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
    gr_sim_set_background_point_mass(sim, cx, cy, GM, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 4);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, grav_inductive_enabled);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/v_circ);

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.2f * (float) n_orbits * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    float r_at_end = r_orb;
    int   s_used = 0;
    for (int s = 0; s < n_max && wraps < n_orbits; s++) {
        gr_sim_step(sim);
        s_used = s + 1;
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) { gr_sim_destroy(sim); return NAN; }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            r_at_end = r_now;
        }
        th_prev = th;
    }
    if (n_completed_out) *n_completed_out = wraps;
    if (n_steps_total_out) *n_steps_total_out = s_used;
    gr_sim_destroy(sim);
    if (wraps < n_orbits) return NAN;
    return (r_at_end - r_orb) / r_orb;
}

int main(void) {
    printf("=== stage29_inductive_cfl_sweep ===\n");
    printf("Does the inductive-piece inspiral rate scale with CFL?\n");
    printf("Setup: gravity, m_test=1e-3, r=20, 4 orbits, half-step A (v36).\n\n");

    const float GM     = 1.0f;
    const float m_test = 1.0e-3f;
    const float r_orb  = 20.0f;
    const int   N      = 4;

    /* Halving CFL each step.  All physically stable (CFL_2D <= 1/sqrt(2)). */
    const float cfls[] = {
        1.0f / sqrtf(2.0f),
        1.0f / (2.0f * sqrtf(2.0f)),
        1.0f / (4.0f * sqrtf(2.0f)),
        1.0f / (8.0f * sqrtf(2.0f)),
    };
    const int n_cfl = (int) (sizeof(cfls) / sizeof(cfls[0]));
    sweep_result_t results[8];

    printf("%-6s %-10s %-10s %-12s %-12s %-12s %-12s\n",
           "i", "CFL", "dt", "steps/4-orb", "no induct", "induct +1", "induct contrib");
    printf("---------------------------------------------------------------------------\n");

    for (int i = 0; i < n_cfl; i++) {
        const float cfl = cfls[i];
        const float dt  = cfl * 1.0f / 1.0f;  /* dx=1, c_eff=1 */
        int n_noind = 0, n_indon = 0;
        int steps_noind = 0, steps_indon = 0;
        const float d_noind = run(GM, m_test, r_orb, cfl, N, 0, &n_noind, &steps_noind);
        const float d_indon = run(GM, m_test, r_orb, cfl, N, 1, &n_indon, &steps_indon);

        results[i].cfl = cfl;
        results[i].dt  = dt;
        results[i].n_steps_total = steps_indon;
        results[i].drift_noind   = d_noind;
        results[i].drift_indon   = d_indon;
        results[i].drift_indcontrib = (isfinite(d_noind) && isfinite(d_indon))
                                          ? (d_indon - d_noind) : NAN;
        results[i].n_completed_noind = n_noind;
        results[i].n_completed_indon = n_indon;

        char no_s[20], in_s[20], dc_s[20];
        if (!isfinite(d_noind)) snprintf(no_s, sizeof(no_s), "<unbound>");
        else snprintf(no_s, sizeof(no_s), "%+7.4f%%", 100.0 * (double) d_noind);
        if (!isfinite(d_indon)) snprintf(in_s, sizeof(in_s), "<unbound>");
        else snprintf(in_s, sizeof(in_s), "%+7.4f%%", 100.0 * (double) d_indon);
        if (!isfinite(results[i].drift_indcontrib)) snprintf(dc_s, sizeof(dc_s), "—");
        else snprintf(dc_s, sizeof(dc_s), "%+7.4f%%", 100.0 * (double) results[i].drift_indcontrib);

        printf("%-6d %-10.5f %-10.5f %-12d %-12s %-12s %-12s\n",
               i, (double) cfl, (double) dt, steps_indon, no_s, in_s, dc_s);
        fflush(stdout);
    }
    printf("\n");

    /* Analysis: report the inductive-contribution scaling between
     * consecutive CFL halvings.  Ratio close to 1 -> CFL-independent
     * (physical).  Ratio close to 4 -> O(dt^2) (numerical truncation).
     * Ratio between 1 and 4 -> mixed. */
    printf("Inductive contribution ratio (consecutive CFL halvings):\n");
    for (int i = 1; i < n_cfl; i++) {
        const float dprev = results[i - 1].drift_indcontrib;
        const float dcur  = results[i].drift_indcontrib;
        if (isfinite(dprev) && isfinite(dcur) && fabsf(dcur) > 1e-8f) {
            const float ratio = dprev / dcur;
            printf("  CFL %.5f -> %.5f:   |contrib_prev / contrib_cur| = %.3f\n",
                   (double) results[i - 1].cfl, (double) results[i].cfl,
                   fabs((double) ratio));
        }
    }
    printf("\nInterpretation:\n");
    printf("  ratio ~ 1.0 : inductive contribution is CFL-independent (physical RR).\n");
    printf("  ratio ~ 4.0 : inductive contribution scales as O(dt^2) (numerical residual).\n");
    printf("  ratio ~ 2.0 : O(dt) -- shouldn't happen for centered leapfrog, would indicate\n");
    printf("                a discretization mismatch.\n");

    /* For reference: analytic Larmor prediction at our parameters.
     *
     * Larmor (3D, gravitoelectric only, applied here as if 2D):
     *   P ~ (2/3) m^2 a^2 / c^3, with a ~ v^2/r at orbit.
     * v_circ ~ 0.22, r=20, a ~ 0.0024.  Orbit energy E ~ (1/2)m v^2.
     * Power ~ (2/3) * 1e-6 * 5.76e-6 = 3.84e-12.
     * Energy E ~ 0.5 * 1e-3 * 0.05 = 2.5e-5.
     * Fractional dE/E per orbit ~ P*T/E = 3.84e-12 * 570 / 2.5e-5 ~ 9e-5,
     * i.e. ~0.01% per orbit.
     * (Order-of-magnitude estimate; 2D-softened gravity gets different prefactors.)
     *
     * We observe a few percent.  Convergence with CFL will say how much of
     * that is physical vs numerical. */

    /* Soft assertion: all four runs should be bound (no <unbound> entries). */
    for (int i = 0; i < n_cfl; i++) {
        TEST_ASSERT(results[i].n_completed_noind == N,
                    "no-inductive run at CFL=%.4f didn't complete %d orbits",
                    (double) cfls[i], N);
        TEST_ASSERT(results[i].n_completed_indon == N,
                    "inductive-on run at CFL=%.4f didn't complete %d orbits",
                    (double) cfls[i], N);
    }

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
