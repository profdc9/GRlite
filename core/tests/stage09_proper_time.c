/* Stage 9 test — gr_sandbox_v34.tex §12.9 "Proper time accumulation".
 *
 * Two-phase test:
 *
 *   Phase A — non-spinning point mass: circular orbits at several radii
 *   verify the test particle's accumulated proper time per orbit matches
 *
 *       d_tau/dt = sqrt( 1 + 2 Phi/c^2 - (1 - 2 Phi/c^2) v^2/c^2 )
 *
 *   evaluated at the (Phi, v) of the relativistic circular orbit.  This
 *   exercises both the gravitational ($2 Phi/c^2$) and kinematic ($v^2/c^2$)
 *   time-dilation contributions in a regime where they're comparable.  The
 *   two contributions stand in a 2:1 ratio for a Newtonian circular orbit
 *   with $v^2 = GM/r$, so changing $r$ slides the balance: gravity dominates
 *   close in, kinematic farther out (since $v$ falls off slower than $\Phi$).
 *
 *   Phase B — spinning point mass: same circular orbit, once prograde and
 *   once retrograde, verifies the gravitomagnetic clock effect (Eq. 75 with
 *   v34 sign + factor corrections):
 *
 *       Delta_tau_{pro - retro} = -8 pi G J_z / (c^4 r)   per orbit
 *
 *   derived from the $-8 (v . A_g)/c^2$ cross-term in d_tau/dt and the
 *   $|A_g| = G J_z / (2 c^2 r^2)$ dipole field at the orbit.
 *
 * Both phases use the analytic-background path (§sec:bg_mode) to eliminate
 * the CIC+FD tangential force artifact. */

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

/* Phase A — non-spinning. */
static int test_static_circular_at_radius(float GM, float r, float eps,
                                          float tol_rel) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");

    const float params[3] = {GM, r, eps};
    TEST_ASSERT(gr_sim_load_scenario(sim, "kepler_orbit", params, 3) == 0,
                "kepler_orbit load failed (r=%g)", r);
    TEST_ASSERT(gr_sim_particle_count(sim) == 1, "particle count != 1");

    /* Relativistic v_circ (matches kepler_orbit.c). */
    const float r2_e2_15 = powf(r * r + eps * eps, 1.5f);
    const float g_mag    = GM * r / r2_e2_15;
    const float rg       = r * g_mag;
    const float rg2_c2   = rg * rg / (c_eff * c_eff);
    const float u_v2     = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ   = sqrtf(u_v2);
    const float Phi      = -GM / sqrtf(r * r + eps * eps);
    const float two_phi_c2 = 2.0f * Phi / (c_eff * c_eff);
    const float v2_c2      = v_circ * v_circ / (c_eff * c_eff);
    const float radicand   = 1.0f + two_phi_c2 - (1.0f - two_phi_c2) * v2_c2;
    const float dtau_dt_ana = sqrtf(radicand);
    const float T_orbit     = 2.0f * (float) M_PI * r / v_circ;

    /* Run for one orbital period.  Track when angle theta crosses through 0
     * with theta increasing (the particle started at theta = 0, going +y) —
     * that's one full revolution. */
    const float cx = ((float) W * 0.5f) * dx;
    const float cy = ((float) H * 0.5f) * dx;
    const int   n_max = (int) (1.2f * T_orbit / gr_sim_dt(sim));
    float prev_theta = 0.0f;
    int   n_done = 0;
    float t_completion = 0.0f;
    int   completed = 0;
    for (int s = 0; s < n_max; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float th = atan2f(p->y - cy, p->x - cx);
        /* Detect counterclockwise crossing from +x: prev theta > 2.5 and
         * current theta < -2.5, OR a linear crossing from -small to +small. */
        if (s >= 10 && prev_theta < -2.5f && th > 2.5f) {
            /* came back around */
        } else if (s >= 10 && prev_theta > 2.5f && th < -2.5f) {
            t_completion = gr_sim_time(sim);
            completed = 1;
            n_done = s + 1;
            break;
        }
        prev_theta = th;
    }
    TEST_ASSERT(completed,
                "particle didn't complete one orbit in %d steps (r=%g)",
                n_max, r);

    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    const float tau_measured = p->proper_time;
    /* The simulation's "one orbit" may end slightly off from analytic
     * T_orbit due to integration error; we report tau and the elapsed
     * coordinate time at completion, then compare d_tau/dt ratio. */
    const float dtau_dt_measured = tau_measured / t_completion;
    const float rel_err = fabsf(dtau_dt_measured - dtau_dt_ana) / dtau_dt_ana;
    /* Break d_tau/dt down: a Newtonian circular orbit at r has v^2 = GM/r,
     * so the SR contribution is -v^2/(2c^2) and the gravitational is +Phi/c^2
     * ≈ -GM/(rc^2), in a 1:2 ratio (gravity wins by 2x).  Report both at
     * leading order for context. */
    const float dt_sr_lin   = -0.5f * v2_c2;
    const float dt_grav_lin = Phi / (c_eff * c_eff);

    printf("  r=%5.1f  v/c=%.4f  Phi/c^2=%+.5f  v^2/2c^2=%.5f  ratio g:SR = %.2f\n",
           r, v_circ / c_eff, two_phi_c2 * 0.5f, 0.5f * v2_c2,
           fabsf(dt_grav_lin / dt_sr_lin));
    printf("           T_orb=%.3f  steps=%d\n", T_orbit, n_done);
    printf("           predicted dtau/dt = %.7f  (full sqrt form)\n", dtau_dt_ana);
    printf("           measured  dtau/dt = %.7f   rel.err = %.3e\n",
           dtau_dt_measured, rel_err);

    TEST_ASSERT(rel_err < tol_rel,
                "r=%g: dtau/dt rel.err %.3e exceeds tol %.3e",
                r, rel_err, tol_rel);
    gr_sim_destroy(sim);
    return 0;
}

/* Phase B — spinning, prograde vs retrograde clock effect. */
static int test_spinning_clock_effect(float GM, float r, float eps,
                                      float Jz, float tol_rel) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);

    /* Two sims: prograde and retrograde, otherwise identical.  Run each for
     * the same number of timesteps (= one orbital period at the
     * non-spinning v_circ) and compare accumulated proper time. */
    float taus[2] = {0.0f, 0.0f};
    int   completed_steps = -1;
    for (int dir_idx = 0; dir_idx < 2; dir_idx++) {
        const float sign = (dir_idx == 0) ? +1.0f : -1.0f;
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
        TEST_ASSERT(sim != NULL, "create failed");
        const float params[5] = {GM, r, eps, Jz, sign};
        TEST_ASSERT(gr_sim_load_scenario(sim, "spinning_orbit", params, 5) == 0,
                    "spinning_orbit load failed");
        /* Disable the Tier-1 gravitomagnetic Lorentz force (+4 m v x B_g)
         * for this CLOCK-only isolation test.  When enabled, the v x B_g
         * piece perturbs prograde and retrograde orbits oppositely (Lense-
         * Thirring frame-dragging on the orbit shape), which adds an O(J_z)
         * contribution to Delta_tau on top of the clock-effect formula
         * being verified here.  Stage 20+ tests the v x B_g force directly. */
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);

        /* Orbital period (non-spinning approximation; the spinning shift is
         * a 1PN correction to the orbit). */
        const float r2_e2_15 = powf(r * r + eps * eps, 1.5f);
        const float g_mag    = GM * r / r2_e2_15;
        const float rg       = r * g_mag;
        const float rg2_c2   = rg * rg / (c_eff * c_eff);
        const float u_v2     = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
        const float v_circ   = sqrtf(u_v2);
        const float T_orbit  = 2.0f * (float) M_PI * r / v_circ;
        const int   n_steps  = (int) (T_orbit / gr_sim_dt(sim) + 0.5f);
        if (completed_steps < 0) completed_steps = n_steps;
        gr_sim_step_n(sim, completed_steps);

        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        taus[dir_idx] = p->proper_time;
        gr_sim_destroy(sim);
    }
    const float tau_pro   = taus[0];
    const float tau_retro = taus[1];
    const float dtau_measured = tau_pro - tau_retro;

    /* Two analytic predictions for the same observable:
     *
     *   (1) Linearized form (textbook clock-effect):
     *         Delta_tau = -8 pi G J_z / (c^4 r)   per orbit
     *       derived from -4 (v.A_g)/c^2 in d_tau/dt and integrated around
     *       one orbit of period T = 2 pi r / v.  Accurate when the
     *       A_g.v term is much smaller than (Phi/c^2) and (v^2/c^2).
     *
     *   (2) Exact-sqrt form (what the integrator actually uses):
     *         d_tau/dt = sqrt(1 + 2 Phi/c^2 - (1-2 Phi/c^2) v^2/c^2
     *                            -+ 8 (v.A_g)/c^2)
     *       evaluated at the prograde and retrograde signs of (v.A_g);
     *       Delta_tau = T_coord * (d_tau/dt_+ - d_tau/dt_-).
     *
     * At our test parameters (Jz=2, r=20, near-extremal Kerr in toy units)
     * the two predictions differ at the ~10% level — small but visible.
     * We report both and tolerance against the exact form, which is the
     * intrinsic accuracy of the integrator. */
    const float G_eff   = 1.0f;
    const float c2      = c_eff * c_eff;
    const float c4      = c2 * c2;
    const float dtau_lin = -8.0f * (float) M_PI * G_eff * Jz / (c4 * r);

    /* Reconstruct the exact prediction.  T_orbit, v_circ, and the (Phi,
     * |A_g|) at the orbit's radius are deterministic from the input. */
    const float r2_e2_15 = powf(r * r + eps * eps, 1.5f);
    const float g_mag    = G_eff * GM * r / r2_e2_15;
    const float rg       = r * g_mag;
    const float rg2_c2   = rg * rg / c2;
    const float u_v2     = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ   = sqrtf(u_v2);
    const float Phi      = -GM / sqrtf(r * r + eps * eps);
    const float Ag_mag   = G_eff * Jz * r / (2.0f * c2 * r2_e2_15);
    const float v_dot_Ag = v_circ * Ag_mag;
    const float two_phi_c2 = 2.0f * Phi / c2;
    const float v2_c2      = v_circ * v_circ / c2;
    const float radicand_pro = 1.0f + two_phi_c2 - (1.0f - two_phi_c2) * v2_c2
                             - 8.0f * v_dot_Ag / c2;
    const float radicand_ret = 1.0f + two_phi_c2 - (1.0f - two_phi_c2) * v2_c2
                             + 8.0f * v_dot_Ag / c2;
    /* Reconstruct total coordinate time from CFL (must match the sim's dt). */
    const float cfl_local  = 1.0f / sqrtf(2.0f);
    const float dt_local   = cfl_local * 1.0f / c_eff;  /* dx = 1, c_eff above */
    const float T_coord    = (float) completed_steps * dt_local;
    const float dtau_exact = T_coord * (sqrtf(radicand_pro) - sqrtf(radicand_ret));

    const float rel_err_lin   = fabsf(dtau_measured - dtau_lin)   / fabsf(dtau_lin);
    const float rel_err_exact = fabsf(dtau_measured - dtau_exact) / fabsf(dtau_exact);

    printf("  r=%5.1f  Jz=%.3f  |A_g|=%.5f  v.A_g=%.5e\n",
           r, Jz, Ag_mag, v_dot_Ag);
    printf("    tau_prograde   = %.6f\n", tau_pro);
    printf("    tau_retrograde = %.6f\n", tau_retro);
    printf("    measured delta = %+.6f\n", dtau_measured);
    printf("    linearized ana = %+.6f  ( -8 pi G Jz / (c^4 r) )   rel.err = %.3e\n",
           dtau_lin, rel_err_lin);
    printf("    exact-sqrt ana = %+.6f  (full radicand difference)  rel.err = %.3e\n",
           dtau_exact, rel_err_exact);

    TEST_ASSERT(rel_err_exact < tol_rel,
                "Gravitomagnetic clock-effect (exact form) rel.err %.3e exceeds tol %.3e",
                rel_err_exact, tol_rel);
    return 0;
}

int main(void) {
    printf("=== stage09_proper_time: gr_sandbox_v34.tex §12.9 ===\n");

    /* Phase A: non-spinning at multiple radii, tolerance 1e-3 (FP +
     * leapfrog integration). */
    printf("\n[A] Static point mass, circular orbits at several radii\n");
    const float radii[] = {10.0f, 20.0f, 40.0f, 80.0f};
    const int   n_r = sizeof(radii) / sizeof(radii[0]);
    for (int i = 0; i < n_r; i++) {
        if (test_static_circular_at_radius(/*GM=*/1.0f, radii[i],
                                           /*eps=*/1.0f,
                                           /*tol_rel=*/1.0e-3f) != 0) return 1;
    }

    /* Phase B: spinning point mass, prograde vs retrograde.
     * Moderate spin Jz=2 at r=20 (near-extremal Kerr in our toy units but
     * well inside the GEM linearized regime: |A_g|/(c^2) ~ GJ/(2 c^4 r^2)
     * ~ 0.0025, comfortably small).  Expected clock effect magnitude is
     * ~2.5 sim time units per orbit, ~0.5% of the orbital period. */
    printf("\n[B] Spinning point mass, gravitomagnetic clock effect\n");
    if (test_spinning_clock_effect(/*GM=*/1.0f, /*r=*/20.0f,
                                   /*eps=*/1.0f, /*Jz=*/2.0f,
                                   /*tol_rel=*/5.0e-3f) != 0) return 1;

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
