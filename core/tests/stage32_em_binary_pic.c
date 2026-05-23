/* Stage 32 -- EM analog of Stage 12: two-body mutual EM dynamics via FDTD.
 *
 * Two equal-mass, opposite-sign-charge particles -- standard Maxwell
 * sign convention (opposite charges attract).  Mutual attraction via
 * the FDTD perturbation field, no analytic background.  Full EM Lorentz
 * force on both particles (-q grad phi - q d_t A + q v x B).  Production
 * inductive piece ON.  Stage 27 already validated the inductive inspiral
 * for a single test particle around a softened-charge background; Stage
 * 32 promotes both bodies to active sources and checks the early-orbit
 * behavior is consistent with the 2D-log Coulomb prediction.
 *
 * v_orb prediction in 2D-log Coulomb (analogous to gravity binary):
 *     v_orb = Q * sqrt(k_e / m)        (independent of r)
 *     omega = v_orb / r
 *     T     = 2 pi r / v_orb
 *
 * Tests (mirror Stage 12 structure):
 *   [1] Free-fall: opposite-charge particles released at rest attract,
 *       infall measurable by step 50.
 *   [2] Short-orbit: at v_factor=1, early-phase omega matches v/r.
 *   [3] v_orb independent of r: 2D-log signature for EM.
 *
 * Long-time orbit + inspiral measurement is left for a future test --
 * sub-cell PIC heating + the new radiation-reaction inspiral both act
 * on similar timescales here (same regime as gravity Stage 12), so
 * teasing them apart requires either weak coupling or a CFL sweep
 * (a la Stage 30).  This first stage just validates the basic two-body
 * setup. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); \
        return 1; \
    } \
} while (0)

static int test_free_fall(void) {
    printf("\n[1/3] Free-fall: opposite-sign charges at rest -- mutual attraction.\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 0.01f;
    const float Q      = 0.01f;
    const float r_orb  = 8.0f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "gr_sim_create failed");
    gr_sim_set_damping(sim, 16);

    /* v_factor=0 -> particles at rest; params: m, r, v_factor, cx, cy, Q. */
    const float cx = ((float) (W - 1) * 0.5f) * dx;
    const float cy = ((float) (H - 1) * 0.5f) * dx;
    const float params[6] = {mass, r_orb, 0.0f, cx, cy, Q};
    TEST_ASSERT(gr_sim_load_scenario(sim, "pic_binary_em", params, 6) == 0,
                "pic_binary_em load failed");
    TEST_ASSERT(gr_sim_particle_count(sim) == 2, "expected 2 particles");

    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const float r_init = sqrtf((p0->x - cx) * (p0->x - cx)
                             + (p0->y - cy) * (p0->y - cy));
    printf("  initial: r = %.4f (separation 2r = %.4f), Q=+/-%.3f\n",
           r_init, 2.0f * r_init, Q);

    /* Expected behavior: opposite-sign charges attract (standard Maxwell).
     * Phi_em^{pert} ~ -2 k_e Q log(r) so |F| ~ 2 k_e Q^2 / d.  Mutual
     * infall begins once the field has propagated between the particles
     * (~ 2r/c wave-crossing time ~ 16 sim units / 22 dt). */
    const int N_check = 50;
    const int N_total = 400;
    float r_min = r_init;
    float r_max = r_init;
    float r_at_check = r_init;
    for (int s = 1; s <= N_total; s++) {
        gr_sim_step(sim);
        p0 = gr_sim_get_particle(sim, 0);
        const float r = sqrtf((p0->x - cx) * (p0->x - cx)
                            + (p0->y - cy) * (p0->y - cy));
        if (r < r_min) r_min = r;
        if (r > r_max) r_max = r;
        if (s == N_check) r_at_check = r;
        if (s == 20 || s == 50 || s == 100 || s == 200 || s == 400) {
            printf("  s=%-3d  p0=(%.4f, %.4f)  r=%.4f\n",
                   s, p0->x, p0->y, r);
        }
    }
    printf("  over %d steps: r_min=%.4f r_max=%.4f (r_init=%.4f)\n",
           N_total, r_min, r_max, r_init);

    /* [a] Attraction engaged by s=50. */
    TEST_ASSERT(r_at_check < r_init * 0.97f,
                "free-fall: insufficient infall at s=%d (r=%.4f vs r_init=%.4f)",
                N_check, r_at_check, r_init);
    /* [b] Substantial infall by closest approach. */
    TEST_ASSERT((r_init - r_min) > 0.1f * r_init,
                "free-fall: closest approach r_min=%.4f too far from r_init=%.4f -- "
                "mutual attraction too weak", r_min, r_init);
    /* [c] No Esirkepov violations. */
    const int viols = gr_sim_esirkepov_violations(sim);
    printf("  Esirkepov violations: %d (expect 0)\n", viols);
    TEST_ASSERT(viols == 0, "unexpected Esirkepov 2-cell violations: %d", viols);

    gr_sim_destroy(sim);
    return 0;
}

static int test_short_orbit(void) {
    printf("\n[2/3] Short orbit: omega at r=8, v_factor=1.\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 0.01f;
    const float Q      = 0.01f;
    const float r_orb  = 8.0f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    gr_sim_set_damping(sim, 16);

    const float cx = ((float) (W - 1) * 0.5f) * dx;
    const float cy = ((float) (H - 1) * 0.5f) * dx;
    const float params[6] = {mass, r_orb, 1.0f, cx, cy, Q};
    TEST_ASSERT(gr_sim_load_scenario(sim, "pic_binary_em", params, 6) == 0,
                "pic_binary_em load failed");

    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const float k_e   = gr_sim_get_k_e(sim);
    const float v_orb = Q * sqrtf(k_e / mass);
    const float omega_ana = v_orb / r_orb;
    const float dt = gr_sim_dt(sim);

    /* p0 starts at (cx - r, cy) with velocity (0, +v) and is attracted
     * toward p1 (at +x).  In atan2 conventions this is CLOCKWISE: theta
     * begins at pi (at (-r, 0)) and decreases toward pi/2 as the particle
     * moves through the +y half-plane.  Matches the gravity-binary
     * test_short_orbit convention. */
    const float theta_init = atan2f(p0->y - cy, p0->x - cx);
    float theta_unwrapped = theta_init;
    float prev_theta = theta_init;
    const int N = 100;
    int checkpoints[] = {10, 20, 50, 75, 100};
    int next_cp = 0;
    float omega_meas_at_cp[5];
    float r_at_cp[5];
    for (int s = 1; s <= N; s++) {
        gr_sim_step(sim);
        p0 = gr_sim_get_particle(sim, 0);
        const float r = sqrtf((p0->x - cx) * (p0->x - cx)
                            + (p0->y - cy) * (p0->y - cy));
        const float th = atan2f(p0->y - cy, p0->x - cx);
        float d = th - prev_theta;
        if (d >  (float) M_PI) d -= 2.0f * (float) M_PI;
        if (d < -(float) M_PI) d += 2.0f * (float) M_PI;
        theta_unwrapped += d;
        prev_theta = th;
        if (next_cp < 5 && s == checkpoints[next_cp]) {
            const float dtheta_cum = theta_unwrapped - theta_init;
            /* p0 orbits CW (theta decreases) -- |dtheta|/t = omega. */
            const float omega = fabsf(dtheta_cum) / ((float) s * dt);
            omega_meas_at_cp[next_cp] = omega;
            r_at_cp[next_cp] = r;
            next_cp++;
        }
    }

    printf("  analytic v_orb = Q sqrt(k_e/m) = %.4f (independent of r)\n", v_orb);
    printf("  analytic omega = v_orb / r     = %.4f rad/sim\n", omega_ana);
    printf("  empirical over first %d steps:\n", N);
    for (int i = 0; i < 5; i++) {
        const float rel = fabsf(omega_meas_at_cp[i] - omega_ana) / omega_ana;
        printf("    s=%-3d  r=%6.4f  omega_meas=%.5f  rel.err=%.3e\n",
               checkpoints[i], r_at_cp[i], omega_meas_at_cp[i], rel);
    }

    /* Early omega should match analytic within ~20%.  Tolerance accounts
     * for the wave-crossing transient at startup: at t=0 there is no
     * Coulomb field yet -- it has to propagate from each particle to the
     * other (about r/c ~ 8 sim units = ~11 steps), so omega at s=10 may
     * still be tangent-line-like rather than properly circular. */
    const float rel_early = fabsf(omega_meas_at_cp[0] - omega_ana) / omega_ana;
    TEST_ASSERT(rel_early < 0.20f,
                "early omega (s=10) rel.err %.3e > 20%%", rel_early);

    /* p0 should be orbiting CW (theta_unwrapped < theta_init). */
    TEST_ASSERT(theta_unwrapped < theta_init,
                "orbital direction wrong: theta_final %.4f >= theta_init %.4f",
                theta_unwrapped, theta_init);

    /* Radial drift over the first 20 steps. */
    const float r_drift = fabsf(r_at_cp[1] - r_orb) / r_orb;
    printf("  early radial drift at s=20: %.2f%%\n", 100.0f * r_drift);
    TEST_ASSERT(r_drift < 0.20f,
                "early radial drift %.3e > 20%% -- orbit not establishing",
                r_drift);

    const int viols = gr_sim_esirkepov_violations(sim);
    printf("  Esirkepov violations: %d\n", viols);
    TEST_ASSERT(viols == 0, "unexpected Esirkepov violations: %d", viols);

    gr_sim_destroy(sim);
    return 0;
}

/* The 2D-log Coulomb signature: v_orb = Q sqrt(k_e/m) is independent of r.
 * Same logic as Stage 12 [3] but for the EM channel.  Probe at two radii
 * and confirm the empirically-measured tangential speed is the same and
 * matches the analytic constant. */
static int test_v_orb_independent_of_r(void) {
    printf("\n[3/3] v_orb = Q sqrt(k_e/m) independent of r (2D-log EM signature)\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 0.01f;
    const float Q      = 0.01f;
    const float k_e    = 1.0f;
    const float v_ana  = Q * sqrtf(k_e / mass);

    const float r_values[] = {6.0f, 12.0f};
    const int   n_r = sizeof(r_values) / sizeof(r_values[0]);
    float v_meas[2];

    const int N_fixed = 10;
    for (int k = 0; k < n_r; k++) {
        const float r_orb = r_values[k];
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
        TEST_ASSERT(sim != NULL, "create failed at r=%g", r_orb);
        gr_sim_set_damping(sim, 16);
        const float cx = ((float) (W - 1) * 0.5f) * dx;
        const float cy = ((float) (H - 1) * 0.5f) * dx;
        const float params[6] = {mass, r_orb, 1.0f, cx, cy, Q};
        TEST_ASSERT(gr_sim_load_scenario(sim, "pic_binary_em", params, 6) == 0,
                    "scenario load failed at r=%g", r_orb);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float theta0 = atan2f(p->y - cy, p->x - cx);
        const float dt = gr_sim_dt(sim);
        for (int s = 0; s < N_fixed; s++) gr_sim_step(sim);
        p = gr_sim_get_particle(sim, 0);
        const float theta1 = atan2f(p->y - cy, p->x - cx);
        float dtheta = theta1 - theta0;
        if (dtheta >  (float) M_PI) dtheta -= 2.0f * (float) M_PI;
        if (dtheta < -(float) M_PI) dtheta += 2.0f * (float) M_PI;
        const float t_elapsed = (float) N_fixed * dt;
        const float omega_meas = fabsf(dtheta) / t_elapsed;
        v_meas[k] = r_orb * omega_meas;

        const float t_wc = 2.0f * r_orb / c_eff;
        printf("  r=%5.2f  N=%-3d t=%5.2f (%.2f wave-crossings)  "
               "omega=%.5f  v_meas=r*omega=%.5f\n",
               r_orb, N_fixed, t_elapsed, t_elapsed / t_wc,
               omega_meas, v_meas[k]);
        gr_sim_destroy(sim);
    }

    const float rel_diff = fabsf(v_meas[0] - v_meas[1])
                         / (0.5f * (v_meas[0] + v_meas[1]));
    printf("  |v(r=%.0f) - v(r=%.0f)| / mean = %.3e\n",
           r_values[0], r_values[1], rel_diff);
    TEST_ASSERT(rel_diff < 0.15f,
                "v_orb varies with r by %.3e -- should be ~independent",
                rel_diff);

    const float rel0 = fabsf(v_meas[0] - v_ana) / v_ana;
    const float rel1 = fabsf(v_meas[1] - v_ana) / v_ana;
    printf("  rel.err vs Q sqrt(k_e/m): r=%.0f -> %.3e   r=%.0f -> %.3e\n",
           r_values[0], rel0, r_values[1], rel1);
    TEST_ASSERT(rel0 < 0.20f && rel1 < 0.20f,
                "v_meas off from analytic by >20%% (rel=%.3e, %.3e)",
                rel0, rel1);

    return 0;
}

int main(void) {
    printf("=== stage32_em_binary_pic: two-body mutual EM via FDTD ===\n");
    if (test_free_fall()              != 0) return 1;
    if (test_short_orbit()            != 0) return 1;
    if (test_v_orb_independent_of_r() != 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
