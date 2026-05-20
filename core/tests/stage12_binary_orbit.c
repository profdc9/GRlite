/* Stage 12 — two-body mutual FDTD gravity (free-fall + bound-orbit checks).
 *
 * gr_sandbox_vNN.tex §sec:yee_migration_plan S7.  This is the milestone
 * scenario the Yee + Esirkepov refactor was built to enable: two equal-mass
 * particles experience mutual gravity through the FDTD perturbation field,
 * with no analytic background.  v34 cell-centered's PIC grid-heating
 * unbound any mutual orbit within a fraction of a wave-crossing time.
 * After S6's Esirkepov continuity and post-S5's HE-adjoint corner chain,
 * the mutual attraction is empirically present and stable on physically
 * meaningful timescales.
 *
 * Two checks:
 *   [1] Free-fall (v_factor = 0): two equal masses released at rest at
 *       separation 2r should accelerate toward each other under their
 *       mutual gravity.  We measure that r (each particle's distance from
 *       COM) decreases monotonically over the early-time evolution.
 *   [2] Circular orbit (v_factor = 1, short horizon): the particles
 *       complete a measurable arc of an orbit before the FDTD-vs-Newton
 *       retardation transient distorts the trajectory.  This is a
 *       robustness check, not a precise dynamics test — full bound orbits
 *       require either pre-equilibrating the field or much weaker
 *       coupling than the time horizon affords. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

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

static int test_free_fall(void) {
    printf("\n[1/2] Free-fall: two equal masses released at rest\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 0.01f;
    const float r_orb  = 8.0f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "gr_sim_create failed");
    gr_sim_set_damping(sim, 16);

    /* v_factor = 0 -> particles at rest. */
    const float params[3] = {mass, r_orb, 0.0f};
    TEST_ASSERT(gr_sim_load_scenario(sim, "pic_binary", params, 3) == 0,
                "pic_binary load failed");
    TEST_ASSERT(gr_sim_particle_count(sim) == 2, "expected 2 particles");

    const float cx = (float) (W / 2) * dx;
    const float cy = (float) (H / 2) * dx;
    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const float r_init = sqrtf((p0->x - cx) * (p0->x - cx)
                             + (p0->y - cy) * (p0->y - cy));
    printf("  initial: r = %.4f (separation 2r = %.4f)\n", r_init, 2.0f * r_init);

    /* Free-fall produces a head-on infall (along the initial separation
     * axis), the particles pass close at some minimum r, then continue
     * outward on the other side.  We sample at a point well past the
     * field-buildup transient (s=50, t~35) but before the closest-approach
     * "slingshot" — the half-period of the linear infall in the 2D log
     * potential for our parameters is ~50-70 steps. */
    const int N_check  = 50;   /* must show infall by here */
    const int N_min    = 100;  /* track to/near closest approach */
    float r_min = r_init;
    float r_at_check = r_init;
    for (int s = 1; s <= N_min; s++) {
        gr_sim_step(sim);
        p0 = gr_sim_get_particle(sim, 0);
        const float r = sqrtf((p0->x - cx) * (p0->x - cx)
                            + (p0->y - cy) * (p0->y - cy));
        if (r < r_min) r_min = r;
        if (s == N_check) r_at_check = r;
        if (s == 20 || s == 50 || s == 100) {
            printf("  s=%-3d  p0=(%.4f, %.4f)  r=%.4f  (infall %.4f)\n",
                   s, p0->x, p0->y, r, r_init - r);
        }
    }
    printf("  min r over %d steps = %.4f (infall %.4f, %.2f%% of r_init)\n",
           N_min, r_min, r_init - r_min, 100.0f * (r_init - r_min) / r_init);

    /* [a] Mutual gravity engages — particles measurably approach within
     * ~50 steps (~35 sim units, ~2 wave-crossing times). */
    TEST_ASSERT(r_at_check < r_init * 0.97f,
                "free-fall: insufficient infall at s=%d (r=%.4f vs r_init=%.4f)",
                N_check, r_at_check, r_init);

    /* [b] Substantial infall by closest approach (>10% of initial r). */
    TEST_ASSERT((r_init - r_min) > 0.1f * r_init,
                "free-fall: closest approach (r_min=%.4f) didn't get within "
                "10%% of r_init=%.4f — mutual gravity too weak", r_min, r_init);

    /* [c] Esirkepov violations: should be 0 (motion well under 1 cell/step). */
    const int viols = gr_sim_esirkepov_violations(sim);
    printf("  Esirkepov violations: %d (expect 0)\n", viols);
    TEST_ASSERT(viols == 0, "unexpected Esirkepov 2-cell violations: %d", viols);

    gr_sim_destroy(sim);
    return 0;
}

static int test_short_orbit(void) {
    printf("\n[2/2] Short-horizon orbit: arc of an attempted circular orbit\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 0.01f;
    const float r_orb  = 8.0f * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    gr_sim_set_damping(sim, 16);

    /* Default v_factor = 1.0 (analytic circular-orbit speed). */
    const float params[2] = {mass, r_orb};
    TEST_ASSERT(gr_sim_load_scenario(sim, "pic_binary", params, 2) == 0,
                "pic_binary load failed");

    const float cx = (float) (W / 2) * dx;
    const float cy = (float) (H / 2) * dx;
    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const float G_eff = gr_sim_get_G_eff(sim);
    const float v_orb = sqrtf(G_eff * mass);
    const float T_ana = 2.0f * (float) M_PI * r_orb / v_orb;
    const float dt = gr_sim_dt(sim);

    /* Run for a few wave-crossing times so the field equilibrates and
     * the orbital arc is measurable. */
    const int N = 30;  /* ~0.06 of an analytic period */
    const float theta_init = atan2f(p0->y - cy, p0->x - cx);
    for (int s = 0; s < N; s++) {
        gr_sim_step(sim);
    }
    p0 = gr_sim_get_particle(sim, 0);
    const float r = sqrtf((p0->x - cx) * (p0->x - cx)
                        + (p0->y - cy) * (p0->y - cy));
    const float theta = atan2f(p0->y - cy, p0->x - cx);
    float dtheta = theta - theta_init;
    if (dtheta >  (float) M_PI) dtheta -= 2.0f * (float) M_PI;
    if (dtheta < -(float) M_PI) dtheta += 2.0f * (float) M_PI;
    /* Expected angular advance over N*dt of an ideal orbit:
     * dtheta_ana = 2 pi * (N * dt / T_ana). */
    const float dtheta_ana = 2.0f * (float) M_PI * (float) N * dt / T_ana;
    /* p0 starts at theta = pi (left side of COM) moving in +y direction,
     * so the unwrapped angle DECREASES toward pi/2.  Absolute value
     * comparison. */
    printf("  after %d steps (t=%.2f, %.1f%% of T_ana=%.2f):\n",
           N, N * dt, 100.0f * N * dt / T_ana, T_ana);
    printf("  p0=(%.4f, %.4f)  r=%.4f (init %.4f)\n",
           p0->x, p0->y, r, r_orb);
    printf("  arc dtheta = %.4f rad (analytic %.4f)\n",
           dtheta, -dtheta_ana);

    /* The orbital arc should be in the expected direction (decreasing
     * theta) and of approximately the right magnitude (within factor 2). */
    TEST_ASSERT(dtheta < 0.0f,
                "orbital direction wrong: dtheta = %.4f (expected < 0)", dtheta);
    const float arc_rel = fabsf(fabsf(dtheta) - dtheta_ana) / dtheta_ana;
    TEST_ASSERT(arc_rel < 1.0f,
                "orbital arc off by >100%%: measured %.4f vs analytic %.4f",
                fabsf(dtheta), dtheta_ana);

    /* Radial drift should be modest at this short horizon. */
    const float r_drift = fabsf(r - r_orb) / r_orb;
    printf("  radial drift = %.3e (%.2f%%)\n", r_drift, 100.0f * r_drift);
    TEST_ASSERT(r_drift < 0.5f,
                "radial drift %.3e too large for short horizon", r_drift);

    const int viols = gr_sim_esirkepov_violations(sim);
    printf("  Esirkepov violations: %d\n", viols);
    TEST_ASSERT(viols == 0, "unexpected Esirkepov violations: %d", viols);

    gr_sim_destroy(sim);
    return 0;
}

int main(void) {
    printf("=== stage12_binary_orbit: two-body mutual gravity via FDTD ===\n");
    if (test_free_fall()   != 0) return 1;
    if (test_short_orbit() != 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
