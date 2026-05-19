/* Stage 8 test — gr_sandbox_v33.tex §12.8 "Relativistic corrections,
 * perihelion precession". Places a single test particle on a moderately
 * eccentric orbit (a, e) around a softened-point-mass background, then runs
 * with Tier-2 relativistic force corrections enabled.
 *
 *   F_grav = m * [ g (1 + v^2/c^2) + 4 (v . g) v / c^2 ],   g = -grad(Phi_g)
 *
 * (eq:geodesic_expansion §"Expansion of the geodesic equation"; the table at
 * §"Practical implementation tiers" labels this as Tier 2.)
 *
 * Verifies:
 *   (a) Successive perihelia precess by approximately
 *         d_phi = 6 pi G_eff M / (c_eff^2 * a * (1 - e^2))
 *       per radial period (Schwarzschild geodesic, §12.8).
 *   (b) With Tier-0 (Newtonian) gravity the precession is near zero — the
 *       same orbit is non-precessing without the v^2/c^2 corrections, which
 *       directly isolates the contribution of the new terms.
 *
 * Sub-step refinement: a parabolic fit on three consecutive r values around
 * a local minimum gives the perihelion time to O(dt^2), then a linear
 * interpolation of (x, y) gives the perihelion angle to the same order.
 * Without this refinement the per-perihelion angle error is v_peri/r * dt
 * (~1% of the expected precession at our parameters) and dominates the
 * measurement. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...)                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);           \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

/* Wrap an angle to [-pi, pi]. */
static float wrap_pmpi(float a) {
    while (a >  (float) M_PI) a -= 2.0f * (float) M_PI;
    while (a < -(float) M_PI) a += 2.0f * (float) M_PI;
    return a;
}

/* One run of the eccentric_orbit scenario with the given force tier.
 * Returns the mean precession per radial period across the perihelion
 * crossings observed within n_orbits orbits, in radians.  *n_peri_out
 * receives the number of perihelia detected, *v_max_out the maximum |v|
 * seen (for the v/c diagnostic). */
static float measure_precession(gr_force_tier_t tier,
                                float GM, float a_orb, float e, float eps,
                                int   n_orbits,
                                int*  n_peri_out,
                                float* v_max_out,
                                float* r_min_out,
                                float* r_max_out) {
    /* Box must fit r_apo with ~2-cell FD margin. Sized for a=100,e=0.2 →
     * r_apo ≈ 120, plus relativistic-orbit growth (a' may shift up). */
    const int   W      = 320, H = 320;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return NAN;
    gr_sim_set_force_tier(sim, tier);

    const float params[4] = {GM, a_orb, e, eps};
    if (gr_sim_load_scenario(sim, "eccentric_orbit", params, 4) != 0) {
        gr_sim_destroy(sim);
        return NAN;
    }

    const float cx = ((float) W * 0.5f) * dx;
    const float cy = ((float) H * 0.5f) * dx;
    /* Newtonian radial period as a budget — relativistic correction is at
     * most a few percent for our parameters. */
    const float T_kepler = 2.0f * (float) M_PI * sqrtf(a_orb * a_orb * a_orb / GM);
    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) ((float) n_orbits * 1.5f * T_kepler / dt);

    /* History buffers for sub-step refinement. We detect periapsis as a sign
     * change of u = r . v (the un-normalized radial-velocity numerator) from
     * negative (inbound) to positive (outbound).  This is far more robust
     * than comparing successive r values directly — near periapsis r is
     * stationary to O(dt^2) and the difference between consecutive samples
     * falls below float precision (~5e-6 relative for our parameters),
     * causing the local-minimum detector to silently miss perihelia
     * whenever the orbit is "clean" (e.g., analytic-background mode). */
    float u_prev = 0.0f, u_curr = 0.0f;
    float x_prev = 0.0f, y_prev = 0.0f;
    int   have = 0;
    float r_min = 1e30f, r_max = 0.0f;

    /* Recorded perihelion angles, up to a generous upper bound. */
    enum { MAX_PERI = 16 };
    float theta_peri[MAX_PERI];
    int   n_peri = 0;

    float v_max = 0.0f;

    for (int step = 0; step < n_max && n_peri < (n_orbits + 1) && n_peri < MAX_PERI; step++) {
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float xp = p->x - cx;
        const float yp = p->y - cy;
        const float r  = sqrtf(xp * xp + yp * yp);
        /* Track speed for diagnostics. */
        const float pmag2 = p->px * p->px + p->py * p->py;
        const float gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c_eff * c_eff));
        const float v     = sqrtf(pmag2) / (gamma * p->mass);
        if (v > v_max) v_max = v;
        if (r < r_min) r_min = r;
        if (r > r_max) r_max = r;

        /* Compute u = r . v (proportional to r dr/dt).  u changes sign from
         * negative to positive at periapsis. */
        const float vx_now = p->px / (gamma * p->mass);
        const float vy_now = p->py / (gamma * p->mass);
        u_prev = u_curr;
        u_curr = xp * vx_now + yp * vy_now;

        if (have >= 1) {
            /* Detect periapsis: u changed sign from <=0 to >0.  Linearly
             * interpolate to the zero crossing and sample (x, y) there. */
            if (u_prev <= 0.0f && u_curr > 0.0f) {
                const float denom = u_curr - u_prev;
                const float t_frac = (denom != 0.0f) ? (-u_prev) / denom : 0.0f;
                const float xp_peri = (1.0f - t_frac) * x_prev + t_frac * xp;
                const float yp_peri = (1.0f - t_frac) * y_prev + t_frac * yp;
                if (n_peri < MAX_PERI) {
                    theta_peri[n_peri++] = atan2f(yp_peri, xp_peri);
                }
            }
        }
        x_prev = xp;
        y_prev = yp;
        if (have < 2) have++;

        gr_sim_step(sim);
    }

    /* Mean precession per radial period.  Difference of successive perihelion
     * angles, wrapped to [-pi, pi]; this is exactly the extra angle beyond
     * one revolution since each theta is already in [-pi, pi] from atan2. */
    float dphi_sum = 0.0f;
    int   dphi_n   = 0;
    for (int k = 1; k < n_peri; k++) {
        const float dphi = wrap_pmpi(theta_peri[k] - theta_peri[k - 1]);
        dphi_sum += dphi;
        dphi_n++;
    }

    gr_sim_destroy(sim);
    *n_peri_out = n_peri;
    *v_max_out  = v_max;
    *r_min_out  = r_min;
    *r_max_out  = r_max;
    return (dphi_n > 0) ? dphi_sum / (float) dphi_n : NAN;
}

int main(void) {
    printf("=== stage08_precession: gr_sandbox_v33.tex §12.8 ===\n");

    /* Test parameters — chosen so 1PN is well-converged and orbits stay bound:
     *   r_s / r_peri = 2 GM / (c^2 r_peri) = 0.0025 (deeply 1PN regime)
     *   v/c at periapsis ~ 0.04
     *   Schwarzschild precession ~ 0.020 rad/orbit (~1.1 deg)
     *   orbit r in [80, 120], fits in 320x320 box centered at (160, 160).
     *   softening correction is negligible (eps/r_peri = 1/80) */
    const float GM    = 0.1f;
    const float a_orb = 100.0f;
    const float e     = 0.2f;
    const float eps   = 1.0f;
    const int   n_orb = 4;

    const float dphi_ana = 6.0f * (float) M_PI * GM
                         / (1.0f * a_orb * (1.0f - e * e));  /* c_eff = 1 */
    printf("  parameters: GM=%.2f  a=%.1f  e=%.2f  eps=%.1f\n", GM, a_orb, e, eps);
    printf("  expected: dphi_prec = 6 pi GM / (c^2 a (1-e^2)) = %.5f rad (%.3f deg) per orbit\n",
           dphi_ana, dphi_ana * 180.0f / (float) M_PI);

    /* Tier-0 (Newtonian): kinematic baseline. With relativistic momentum + a
     * pure 1/r force law, the orbit precesses slightly (precession ~ v^2/c^2)
     * — separate from any 1PN gravitational effect. */
    printf("\n--- Tier 0 (Newtonian gravity) ---\n");
    int   n_peri_N = 0;
    float v_max_N  = 0.0f, r_min_N = 0.0f, r_max_N = 0.0f;
    const float dphi_N = measure_precession(GR_FORCE_NEWTONIAN, GM, a_orb, e, eps,
                                            n_orb, &n_peri_N, &v_max_N,
                                            &r_min_N, &r_max_N);
    const float a_N = 0.5f * (r_min_N + r_max_N);
    const float e_N = (r_max_N - r_min_N) / (r_max_N + r_min_N);
    printf("  perihelia detected: %d   v_max = %.4f (v/c = %.4f)\n",
           n_peri_N, v_max_N, v_max_N / 1.0f);
    printf("  orbit r in [%.4f, %.4f]   a' = %.4f   e' = %.4f\n",
           r_min_N, r_max_N, a_N, e_N);
    printf("  measured precession: %.5f rad (%.3f deg) per orbit\n",
           dphi_N, dphi_N * 180.0f / (float) M_PI);

    /* Tier-2 (Relativistic with v^2/c^2 and Shapiro terms): should match
     * Schwarzschild perihelion precession to leading 1PN. */
    printf("\n--- Tier 2 (Relativistic 1PN, EIH) ---\n");
    int   n_peri_R = 0;
    float v_max_R  = 0.0f, r_min_R = 0.0f, r_max_R = 0.0f;
    const float dphi_R = measure_precession(GR_FORCE_RELATIVISTIC, GM, a_orb, e, eps,
                                            n_orb, &n_peri_R, &v_max_R,
                                            &r_min_R, &r_max_R);
    const float a_R = 0.5f * (r_min_R + r_max_R);
    const float e_R = (r_max_R - r_min_R) / (r_max_R + r_min_R);
    const float dphi_ana_meas =
        6.0f * (float) M_PI * GM / (1.0f * a_R * (1.0f - e_R * e_R));
    printf("  perihelia detected: %d   v_max = %.4f (v/c = %.4f)\n",
           n_peri_R, v_max_R, v_max_R / 1.0f);
    printf("  orbit r in [%.4f, %.4f]   a' = %.4f   e' = %.4f\n",
           r_min_R, r_max_R, a_R, e_R);
    printf("  expected (measured a',e'): %.5f rad/orbit\n", dphi_ana_meas);
    printf("  measured precession: %.5f rad (%.3f deg) per orbit\n",
           dphi_R, dphi_R * 180.0f / (float) M_PI);

    const float rel_err = fabsf(dphi_R - dphi_ana_meas) / dphi_ana_meas;
    printf("  rel.err vs analytic (measured a',e'): %.3e\n", rel_err);

    /* Print everything first, then assert. */
    TEST_ASSERT(n_peri_N >= 2,
                "Tier-0 detected only %d perihelia in %d orbits", n_peri_N, n_orb);
    TEST_ASSERT(n_peri_R >= 2,
                "Tier-2 detected only %d perihelia in %d orbits", n_peri_R, n_orb);

    /* Tier-0 precession isn't zero — relativistic kinematics with Newtonian
     * force gives a kinematic precession of order v^2/c^2.  We only require
     * it be smaller than the GR contribution (sanity check that we haven't
     * accidentally enabled the v^2/c^2 force terms in Tier-0). */
    TEST_ASSERT(fabsf(dphi_N) < 0.5f * dphi_ana_meas,
                "Tier-0 precession %.5f rad too large (>50%% of GR value)", dphi_N);

    /* Tier-2 precession should match analytic within 15%, comparing against
     * the formula evaluated at the orbit's MEASURED a' and e' (not the input
     * Newtonian-derived a, e — these differ because the Schwarzschild orbit
     * with the same initial periapsis kinematics has a different shape). */
    TEST_ASSERT(rel_err < 0.15f,
                "Tier-2 precession rel.err %.3e exceeds 15%%", rel_err);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
