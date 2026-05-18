/* Stage 7 test — gr_sandbox_v33.tex §12.7 "Single test particle, Boris pusher,
 * Keplerian orbit". Places a single uncharged test particle on a circular
 * orbit around a softened-point-mass background and verifies:
 *
 *   (a) the orbital period matches the softened-force-law prediction
 *       T_soft = 2 pi sqrt(r / |g|),  |g| = GM r / (r^2 + eps^2)^{3/2}
 *       (which reduces to T_Kepler = 2 pi sqrt(r^3 / GM) at r >> eps);
 *   (b) the total energy E = gamma m c^2 + m Phi_g is approximately
 *       conserved over many orbits — no secular drift, consistent with
 *       the Boris-leapfrog symplectic structure (v32 §9.2);
 *   (c) the orbit stays circular — radial drift bounded by tolerance,
 *       no spiral-in or spiral-out from numerical dissipation.
 *
 * Choice of parameters: GM = 1, r = 20, eps = 1 sets r/eps = 20, so the
 * softening correction to Newtonian g is (r^2 / (r^2 + eps^2))^{3/2} ≈
 * 0.9925 — the period differs from Newtonian by ~0.4%, well within our
 * test tolerance. */

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

/* Guard against dt = 0 or NaN when sizing the step budget. */
static float sim_dt_safe(float dt) {
    return (dt > 0.0f && isfinite(dt)) ? dt : 1.0f;
}

int main(void) {
    printf("=== stage07_kepler: gr_sandbox_v33.tex §12.7 ===\n");

    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float GM     = 1.0f;
    const float r_orb  = 20.0f;
    const float eps    = 1.0f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    /* No damping — Stage 7 is just a particle in a fixed background; no
     * waves to absorb. */

    float params[3] = {GM, r_orb, eps};
    TEST_ASSERT(gr_sim_load_scenario(sim, "kepler_orbit", params, 3) == 0,
                "kepler_orbit load failed");

    TEST_ASSERT(gr_sim_particle_count(sim) == 1, "expected exactly 1 particle");

    /* Analytic orbit period for the softened force law, *relativistic*
     * circular-orbit condition gamma * v^2 / r = |g|. At v/c ~ 0.2 the
     * relativistic correction is ~1.3% and matters for sub-percent period
     * agreement. */
    const float g_mag_ana = GM * r_orb / powf(r_orb * r_orb + eps * eps, 1.5f);
    const float rg        = r_orb * g_mag_ana;
    const float rg2_c2    = rg * rg / (c_eff * c_eff);
    const float u_ana     = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ    = sqrtf(u_ana);
    const float T_soft    = 2.0f * (float) M_PI * r_orb / v_circ;
    const float T_kepler  = 2.0f * (float) M_PI * sqrtf(r_orb * r_orb * r_orb / GM);
    printf("  GM=%.2f  r=%.1f  eps=%.1f\n", GM, r_orb, eps);
    printf("  v_circ (relativistic) = %.5f   v/c = %.4f   gamma = %.4f\n",
           v_circ, v_circ / c_eff, 1.0f / sqrtf(1.0f - (v_circ * v_circ) / (c_eff * c_eff)));
    printf("  expected period: T_relativistic = %.3f   (T_Newtonian = %.3f, %.2f%% diff)\n",
           T_soft, T_kepler, 100.0f * (T_kepler - T_soft) / T_soft);

    /* Measure orbit period by tracking when the particle's azimuthal angle
     * (relative to the source center) crosses 0 in the +x direction with
     * positive d_theta/dt. Run for ~5 expected periods, record the times of
     * the first 5 crossings, derive period from the average spacing. */
    const float cx    = ((float) W * 0.5f) * dx;
    const float cy    = ((float) H * 0.5f) * dx;
    const int   n_max = (int) (6.0f * T_soft / sim_dt_safe(gr_sim_dt(sim)));  /* see helper */
    /* Energy snapshot at t=0 for the conservation diagnostic. */
    const float E0 = gr_sim_particle_energy(sim, 0);
    printf("  E(t=0) = %.6f\n", E0);

    float prev_theta = 0.0f;
    int   crossings  = 0;
    float t_cross[6] = {0.0f};
    float r_min = r_orb, r_max = r_orb;
    float E_min = E0,    E_max = E0;
    /* Run until 5 full crossings or step budget exhausted. */
    {
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        prev_theta = atan2f(p->y - cy, p->x - cx);
    }
    int n_steps_done = 0;
    while (crossings < 5 && n_steps_done < n_max) {
        gr_sim_step(sim);
        n_steps_done++;
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float theta = atan2f(p->y - cy, p->x - cx);
        const float r_now = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
        if (r_now < r_min) r_min = r_now;
        if (r_now > r_max) r_max = r_now;
        const float E_now = gr_sim_particle_energy(sim, 0);
        if (E_now < E_min) E_min = E_now;
        if (E_now > E_max) E_max = E_now;
        /* Crossing detection: theta jumps from near +pi to near -pi
         * (counterclockwise orbit completes one revolution). */
        if (prev_theta < -2.5f && theta > 2.5f) {
            /* clockwise crossing — wrong direction; ignore */
        } else if (prev_theta > 2.5f && theta < -2.5f) {
            /* counterclockwise crossing — record */
            if (crossings < 6) {
                t_cross[crossings] = gr_sim_time(sim);
                crossings++;
            }
        }
        prev_theta = theta;
    }

    TEST_ASSERT(crossings >= 5, "only %d period crossings in %d steps (budget %d)",
                crossings, n_steps_done, n_max);

    /* Period from the average inter-crossing spacing. */
    float T_avg = 0.0f;
    for (int i = 1; i < crossings; i++) T_avg += t_cross[i] - t_cross[i - 1];
    T_avg /= (float) (crossings - 1);
    const float T_err = fabsf(T_avg - T_soft) / T_soft;
    printf("  measured period (avg over %d crossings): T = %.3f   rel.err vs T_softened = %.2e\n",
           crossings, T_avg, T_err);

    /* Radial drift over the observation window. */
    const float dr_rel = (r_max - r_min) / r_orb;
    printf("  radial drift: r in [%.4f, %.4f] (rel. range %.3e)\n", r_min, r_max, dr_rel);

    /* Energy drift over the observation window. */
    const float dE_rel = (E_max - E_min) / fabsf(E0);
    printf("  energy   drift: E in [%.6f, %.6f] (rel. range %.3e)\n", E_min, E_max, dE_rel);

    /* Assertions — print everything first, then check. */
    TEST_ASSERT(T_err  < 1.0e-2f, "period rel.err %.3e exceeds 1%%", T_err);
    TEST_ASSERT(dr_rel < 1.0e-2f, "radial drift %.3e exceeds 1%%", dr_rel);
    TEST_ASSERT(dE_rel < 1.0e-3f, "energy drift %.3e exceeds 1e-3", dE_rel);

    gr_sim_destroy(sim);
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
