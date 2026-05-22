/* Stage 26 — EM Kepler analog: charged test particle in a softened
 * Coulomb background.
 *
 * Structurally identical to Stage 7 (Kepler in a point-mass background)
 * but with EM force instead of gravity:
 *
 *   Gravity:  Phi_g = -G M / sqrt(r^2 + eps^2),  F = -m grad Phi_g
 *   EM:       phi   = +k_e Q / sqrt(r^2 + eps^2),  F = -q grad phi
 *
 * For opposite-sign Q and q the EM force is attractive; with
 * |q Q| k_e = G M = 1 the orbital dynamics are mathematically
 * identical to Stage 7's gravity case at the same r and eps.  This
 * stage just verifies that the EM force chain (POINT_CHARGE bg ->
 * gr_bg_eval_phi_em -> phi_em_grad_at_total -> em_force_at) is wired
 * correctly end-to-end on a closed-orbit observable, before turning
 * on the perturbation FDTD in Stage 27.
 *
 * Setup:
 *   Q = +1   (central charge), q = -1 (test charge, attractive)
 *   eps = 1, r = 20, c = 1, k_e = 1
 *   field_evolution OFF, particle_source_deposition OFF, ANALYTIC bg
 *   Newtonian tier (no Tier-2 EIH; EM force law is exactly relativistic
 *      already, see v35 §sec:alg_rel line 1001 — no Tier-2 correction
 *      is needed for the EM Lorentz form).
 *
 * Expected: T_orbit = 570.09 (same as Stage 7), orbit stays circular,
 * energy approximately conserved (no perturbation FDTD active). */

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

int main(void) {
    printf("=== stage26_em_kepler ===\n");
    printf("Charged particle in softened-Coulomb background; EM analog of Stage 7.\n");
    printf("Spec: gr_sandbox_v35.tex sec:alg_rel Tier-3 eqbox line 1045.\n\n");

    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float Q      = +1.0f;     /* central charge */
    const float q_test = -1.0f;     /* test charge: attractive */
    const float r_orb  = 20.0f;
    const float eps    = 1.0f;
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;

    /* Relativistic circular velocity: gamma v^2 = g r, where the
     * acceleration magnitude is g = |q Q| k_e r / (r^2 + eps^2)^{3/2}
     * divided by m_test = 1. */
    const float k_e    = 1.0f;
    const float m_test = 1.0f;
    const float g_mag  = fabsf(q_test * Q) * k_e * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f) / m_test;
    const float rg     = r_orb * g_mag;
    const float rg2_c2 = rg * rg / (c_eff * c_eff);
    const float u_v2   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ = sqrtf(u_v2);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_circ;

    printf("Parameters: Q=%+g, q=%+g, r=%g, eps=%g, k_e=%g\n",
           (double) Q, (double) q_test, (double) r_orb,
           (double) eps, (double) k_e);
    printf("Analytic predictions:\n");
    printf("  g_mag   = %.6e\n", (double) g_mag);
    printf("  v_circ  = %.6f  (= sqrt(u) from gamma v^2 = g r)\n", (double) v_circ);
    printf("  gamma   = %.5f\n", (double) (1.0f / sqrtf(1.0f - v_circ * v_circ / (c_eff * c_eff))));
    printf("  T_orbit = %.4f\n\n", (double) T_ana);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test,
                        /*vx=*/0.0f, /*vy=*/v_circ);

    /* Step until two wraps.  Wraps happen at odd multiples of T/2 (the
     * diametrically-opposite side from the starting point), so wrap 1 is
     * at t = T/2 and wrap 2 is at t = 3T/2.  The interval between them
     * is one full orbit, and the particle position at wrap 2 should
     * return to the position at wrap 1 (both at theta = pi). */
    const float dt   = gr_sim_dt(sim);
    const int   n_max = (int) (3.0f * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    float t_first_wrap = 0.0f;
    float t_last_wrap  = 0.0f;
    float x_first_wrap = 0.0f, y_first_wrap = 0.0f;
    float x_last_wrap  = 0.0f, y_last_wrap  = 0.0f;
    float r_min = r_orb, r_max = r_orb;
    for (int s = 0; s < n_max && wraps < 2; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) {
            fprintf(stderr, "FAIL: NaN at step %d\n", s);
            return 1;
        }
        if (r_now < r_min) r_min = r_now;
        if (r_now > r_max) r_max = r_now;
        const float th = atan2f(ry, rx);
        /* Prograde (+v_y at +x): theta increases CCW, wraps +pi -> -pi. */
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            if (wraps == 1) {
                t_first_wrap = gr_sim_time(sim);
                x_first_wrap = p->x;
                y_first_wrap = p->y;
            }
            t_last_wrap = gr_sim_time(sim);
            x_last_wrap = p->x;
            y_last_wrap = p->y;
        }
        th_prev = th;
    }
    TEST_ASSERT(wraps == 2, "didn't see 2 wraps in %d steps (got %d)", n_max, wraps);

    const float T_meas = t_last_wrap - t_first_wrap;
    const float T_err_frac = fabsf(T_meas - T_ana) / T_ana;
    const float radial_range = r_max - r_min;
    const float radial_range_frac = radial_range / r_orb;

    /* Closure: position at wrap 2 should match position at wrap 1 (both
     * at the diametrically-opposite side from the starting point, one
     * full orbital period apart). */
    const float dxw = x_last_wrap - x_first_wrap;
    const float dyw = y_last_wrap - y_first_wrap;
    const float closure = sqrtf(dxw * dxw + dyw * dyw);
    const float closure_frac = closure / r_orb;

    printf("Measured:\n");
    printf("  T_orbit       = %.4f      (err %.4f%%)\n",
           (double) T_meas, 100.0 * (double) T_err_frac);
    printf("  r range       = [%.4f, %.4f]  (range %.4f%% of r)\n",
           (double) r_min, (double) r_max, 100.0 * (double) radial_range_frac);
    printf("  closure err   = %.4e  (%.4f%% of r)\n",
           (double) closure, 100.0 * (double) closure_frac);

    /* Tolerances mirror Stage 7's level (relativistic leapfrog at v/c=0.22,
     * one full orbit). */
    TEST_ASSERT(T_err_frac < 1.0e-3f,
                "T_orbit rel.err %.4f%% exceeds 0.1%%",
                100.0 * (double) T_err_frac);
    TEST_ASSERT(radial_range_frac < 5.0e-4f,
                "radial range %.4f%% exceeds 0.05%% (orbit not circular)",
                100.0 * (double) radial_range_frac);
    TEST_ASSERT(closure_frac < 0.01f,
                "closure error %.4f%% exceeds 1%%",
                100.0 * (double) closure_frac);

    gr_sim_destroy(sim);
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
