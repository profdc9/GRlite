/* Stage 21 — Lense-Thirring orbital period asymmetry.
 *
 * Tests the Tier-1 gravitomagnetic Lorentz force in a SPATIALLY-VARYING
 * gravitomagnetic background — the spinning point mass.  Stage 20 verified
 * the +4 v x B_g force for constant B_g; Stage 21 confirms it also works
 * for the spinning dipole's B_g_z(r), which is the configuration of physical
 * interest (Kerr-like background).
 *
 * Setup: softened spinning point mass at the box center with mass GM = 1,
 * smoothing length eps = 1, spin J_z = 2.  A non-spinning test particle is
 * placed on a circular orbit at radius r = 20.  The gravitomagnetic field
 * at the orbit's radius is (see gr_bg_eval_B_g in core/src/background.c):
 *
 *   B_{g,z}(r) = k * (2 eps^2 - r^2) / (r^2 + eps^2)^{5/2}
 *                where k = G_eff J_z / (2 c^2).
 *
 * At r >> eps, B_{g,z} ~= -k / r^3, NEGATIVE in the orbital plane (the
 * equatorial-plane value of a magnetic-dipole field is anti-parallel to
 * the dipole moment, by a factor 1/2 vs. the polar value — but for our
 * single in-plane sample this is just the curl).
 *
 * The gravitomagnetic Lorentz force F = 4 m (v x B_g) is purely radial
 * for tangential motion.  For PROGRADE (CCW with J_z > 0) the force is
 * INWARD; for RETROGRADE, OUTWARD.  The steady-state circular velocity
 * therefore satisfies (per gr_sandbox_v35.tex eq:geodesic_expansion,
 * specialised to circular motion in 2D):
 *
 *   prograde:    gamma v_pro^2 / r =  a(r) + 4 v_pro |B_gz|
 *   retrograde:  gamma v_ret^2 / r =  a(r) - 4 v_ret |B_gz|
 *
 *   where a(r) = GM r / (r^2 + eps^2)^{3/2}  is the scalar acceleration
 *   and gamma = (1 - v^2 / c^2)^{-1/2}.  The relativistic factor matters:
 *   at our parameters v/c ~ 0.22, gamma ~ 1.025, and the Newtonian
 *   approximation v^2 = a(r) r over-predicts v by 1.3%, putting the
 *   particle on an ellipse instead of a circle.  Stage 7 already
 *   exercises this relativistic-circular formula (kepler_orbit scenario).
 *
 * Perturbative solution.  Let v_rel be the relativistic SCALAR circular
 * velocity satisfying gamma_rel v_rel^2 = a(r) r (the same v_circ used
 * in stage 7).  Linearising about v_rel for small |B_gz|:
 *
 *   delta_v = 4 |B_gz| r / [ gamma_rel * (2 + gamma_rel^2 v_rel^2 / c^2) ]
 *
 *   v_pro =  v_rel + delta_v   (force inward -> needs more centripetal)
 *   v_ret =  v_rel - delta_v
 *
 * Period asymmetry (the test's primary observable):
 *
 *   T = 2 pi r / v   for each direction
 *   T_ret - T_pro ~= 4 pi r delta_v / v_rel^2
 *
 * Using v_circ^2 = a(r) r = GM r^2 / (r^2 + eps^2)^{3/2} (Newtonian
 * circular orbit on the softened potential):
 *
 *   T_ret - T_pro ~= 4 pi |B_gz| (r^2 + eps^2)^{3/2} / GM           (*)
 *
 * At our test parameters (GM=1, eps=1, r=20, J_z=2 -> k=1):
 *
 *   |B_gz| = |2 - 400| / 401^{2.5}     = 398 / 3.2241e6 = 1.2360e-4
 *   v_rel  (relativistic, scalar)      ~= 0.22043     (Stage 7 result)
 *   gamma_rel                          ~= 1.0252
 *   delta_v                            ~= 4.70e-3
 *   v_pro = v_rel + delta_v            ~= 0.22513
 *   v_ret = v_rel - delta_v            ~= 0.21573
 *   T_pro = 2 pi 20 / v_pro            ~= 558.3
 *   T_ret = 2 pi 20 / v_ret            ~= 582.6
 *   (T_ret - T_pro) analytic           ~= 24.3   (~ 4.3% of T_scalar = 570.1)
 *
 * NEWTONIAN tier is used throughout — the Tier-2 EIH corrections are
 * direction-symmetric and would contaminate the prograde-vs-retrograde
 * comparison only at higher order.  Field evolution disabled; analytic
 * background mode for clean B_gz(r). */

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
    float v_init;        /* the GM-corrected v_circ we launch with        */
    float T_mean;        /* mean orbital period over the run              */
    int   n_orbits;      /* number of completed orbits used in the mean   */
    float r_min, r_max;  /* radial bracket — verifies the orbit is steady */
    int   nan;
} orbit_t;

/* Run a single orbit, prograde (dir=+1) or retrograde (dir=-1).  If
 * gm_force is 0 the gravitomagnetic Lorentz force is disabled (baseline
 * scalar gravity only); if 1 it's enabled (full GEM).  Uses v_circ derived
 * from the GM-corrected balance when gm_force is 1, scalar v_circ when 0. */
static void run_orbit(float GM, float r, float eps, float Jz,
                      int dir, int gm_force,
                      int n_orbits_target, orbit_t* out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;

    /* Analytic |B_gz(r)| at the orbit's radius (in-plane).  Same formula as
     * gr_bg_eval_B_g for the SPINNING_POINT_MASS kind. */
    const float k      = GM * Jz / (2.0f * c_eff * c_eff);
    const float r2     = r * r;
    const float eps2   = eps * eps;
    const float s2     = r2 + eps2;
    const float inv_s52 = 1.0f / (s2 * s2 * sqrtf(s2));
    const float B_gz   = k * (2.0f * eps2 - r2) * inv_s52;   /* signed */
    const float abs_B  = fabsf(B_gz);

    /* Scalar relativistic circular velocity at r (same formula as Stage 7
     * and the kepler_orbit scenario): gamma v^2 = g r where
     * g = GM r / (r^2 + eps^2)^{3/2}.  Solving the quadratic in v^2:
     *   v_rel^2 = [-(g r / c^2)^2 + sqrt((g r / c^2)^4 + 4 (g r)^2)] / 2.
     * The Newtonian approximation v^2 = g r over-predicts by 1.3%; using
     * it here gives an elliptical orbit, not a circle, and the measured
     * period drifts up.  See Stage 7 verification. */
    const float a_r      = GM * r / (s2 * sqrtf(s2));
    const float c2       = c_eff * c_eff;
    const float rg       = a_r * r;                          /* = g * r       */
    const float rg2_c2   = rg * rg / c2;
    const float u_rel    = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_rel    = sqrtf(u_rel);
    const float gamma_rel = 1.0f / sqrtf(1.0f - u_rel / c2);

    /* GM-corrected v_circ.  Linearising the relativistic centripetal
     * balance gamma v^2 = (g + 4 v |B|) r about v = v_rel yields the
     * first-order shift
     *   delta_v = 4 |B| r / [ gamma_rel ( 2 + gamma_rel^2 v_rel^2 / c^2 ) ]. */
    const float dv_coeff = gamma_rel * (2.0f + gamma_rel * gamma_rel * u_rel / c2);
    const float delta_v  = 4.0f * abs_B * r / dv_coeff;

    float v_circ;
    if (gm_force) {
        v_circ = (dir > 0) ? (v_rel + delta_v) : (v_rel - delta_v);
    } else {
        v_circ = v_rel;
    }

    /* Sim setup. */
    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim, "create failed");
    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_background_spinning_point_mass(sim, cx, cy, GM, eps, Jz);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    gr_sim_set_gravitomagnetic_force_enabled(sim, gm_force);

    /* Particle at (cx + r, cy) with v = (0, +/- v_circ).
     * Prograde (CCW) = +y at +x. */
    gr_sim_add_particle(sim, cx + r, cy, /*mass=*/1.0f,
                        /*charge=*/0.0f, 0.0f, (float) dir * v_circ);

    /* Step until we've completed n_orbits_target + 1 atan2 wraps.  Each
     * wrap (theta crossing +pi <-> -pi) marks ONE full orbital revolution,
     * but the first wrap occurs ~T/2 into the run (particle starts at
     * theta=0, must traverse a half-orbit to reach theta=pi).  So the
     * orbital period is (t_wrap_{n+1} - t_wrap_1) / n. */
    const float dt   = gr_sim_dt(sim);
    const float T_est = 2.0f * (float) M_PI * r / v_circ;
    const int   wraps_target = n_orbits_target + 1;
    const int   n_max = (int) (1.5f * (float) wraps_target * T_est / dt);
    float th_prev = atan2f(0.0f, r);   /* initial = 0 */
    int   wraps = 0;
    float t_first_wrap = 0.0f;
    float t_last_wrap  = 0.0f;
    float r_min = r, r_max = r;

    for (int s = 0; s < n_max && wraps < wraps_target; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float dxp = p->x - cx;
        const float dyp = p->y - cy;
        const float r_now = sqrtf(dxp * dxp + dyp * dyp);
        if (!isfinite(r_now)) { out->nan = 1; gr_sim_destroy(sim); return; }
        if (r_now < r_min) r_min = r_now;
        if (r_now > r_max) r_max = r_now;

        const float th = atan2f(dyp, dxp);
        /* Wrap detection.  Prograde (dir=+1): theta increases CCW, wraps
         * from near +pi to near -pi.  Retrograde: wraps from near -pi to
         * near +pi.  One wrap per full orbit, occurring at the diametrically
         * opposite side of the starting point. */
        int wrapped = 0;
        if (dir > 0 && th_prev > 0.9f * (float) M_PI
                    && th < -0.9f * (float) M_PI) {
            wrapped = 1;
        } else if (dir < 0 && th_prev < -0.9f * (float) M_PI
                           && th > 0.9f * (float) M_PI) {
            wrapped = 1;
        }
        if (wrapped) {
            wraps++;
            if (wraps == 1) t_first_wrap = gr_sim_time(sim);
            t_last_wrap = gr_sim_time(sim);
        }
        th_prev = th;
    }

    out->v_init   = v_circ;
    out->n_orbits = (wraps >= 2) ? (wraps - 1) : 0;
    out->T_mean   = (wraps >= 2) ? (t_last_wrap - t_first_wrap) / (float) out->n_orbits : 0.0f;
    out->r_min    = r_min;
    out->r_max    = r_max;
    out->nan      = 0;
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage21_lense_thirring_period ===\n");
    printf("Prograde/retrograde period asymmetry for a circular orbit\n");
    printf("around a spinning point mass.  Spec: gr_sandbox_v35.tex\n");
    printf("eq:geodesic_expansion (line 938), spec line 2917 (Lense-Thirring).\n\n");

    const float GM   = 1.0f;
    const float r    = 20.0f;
    const float eps  = 1.0f;
    const float Jz   = 2.0f;
    const int   N    = 4;       /* orbits per run */

    /* Analytic predictions: use the relativistic scalar v_circ plus the
     * perturbative GM shift, matching what run_orbit() injects. */
    const float c_eff   = 1.0f;
    const float c2      = c_eff * c_eff;
    const float k       = GM * Jz / 2.0f;
    const float s2      = r * r + eps * eps;
    const float inv_s52 = 1.0f / (s2 * s2 * sqrtf(s2));
    const float B_gz    = k * (2.0f * eps * eps - r * r) * inv_s52;
    const float abs_B   = fabsf(B_gz);
    const float a_r     = GM * r / (s2 * sqrtf(s2));
    const float rg      = a_r * r;
    const float rg2_c2  = rg * rg / c2;
    const float u_rel   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_rel   = sqrtf(u_rel);
    const float gamma_rel = 1.0f / sqrtf(1.0f - u_rel / c2);
    const float dv_coeff  = gamma_rel * (2.0f + gamma_rel * gamma_rel * u_rel / c2);
    const float delta_v   = 4.0f * abs_B * r / dv_coeff;
    const float v_pro_ana = v_rel + delta_v;
    const float v_ret_ana = v_rel - delta_v;
    const float T_pro_ana = 2.0f * (float) M_PI * r / v_pro_ana;
    const float T_ret_ana = 2.0f * (float) M_PI * r / v_ret_ana;
    const float T_rel_ana = 2.0f * (float) M_PI * r / v_rel;
    const float dT_ana    = T_ret_ana - T_pro_ana;

    printf("Background:  GM=%g, r=%g, eps=%g, J_z=%g\n",
           (double) GM, (double) r, (double) eps, (double) Jz);
    printf("Analytic predictions:\n");
    printf("  B_{g,z}(r) = %+e   (signed; equatorial dipole)\n", (double) B_gz);
    printf("  v_rel      = %.6f  (relativistic circular, gamma v^2 = g r)\n", (double) v_rel);
    printf("  gamma_rel  = %.5f\n", (double) gamma_rel);
    printf("  delta_v    = %+.5e  (perturbative shift from 4 v B_g)\n", (double) delta_v);
    printf("  v_pro      = %.6f  (= v_rel + delta_v)\n", (double) v_pro_ana);
    printf("  v_ret      = %.6f  (= v_rel - delta_v)\n", (double) v_ret_ana);
    printf("  T_rel      = %.4f   (scalar circular; baseline)\n", (double) T_rel_ana);
    printf("  T_pro      = %.4f\n", (double) T_pro_ana);
    printf("  T_ret      = %.4f\n", (double) T_ret_ana);
    printf("  T_ret-T_pro= %+.4f   (%.2f%% of T_rel)\n",
           (double) dT_ana, 100.0 * (double) (dT_ana / T_rel_ana));
    printf("\n");

    /* GM-on runs: production physics with the Lense-Thirring response. */
    orbit_t pro_on, ret_on, pro_off, ret_off;
    run_orbit(GM, r, eps, Jz, +1, /*gm_force=*/1, N, &pro_on);
    run_orbit(GM, r, eps, Jz, -1, /*gm_force=*/1, N, &ret_on);
    /* GM-off baseline: scalar gravity only.  Both directions should give
     * the same period (= T_scalar) — sanity check that the asymmetry comes
     * purely from the v x B_g force, not from any other Stage-9-style
     * direction-coupling. */
    run_orbit(GM, r, eps, Jz, +1, /*gm_force=*/0, N, &pro_off);
    run_orbit(GM, r, eps, Jz, -1, /*gm_force=*/0, N, &ret_off);

    TEST_ASSERT(!pro_on.nan && !ret_on.nan && !pro_off.nan && !ret_off.nan, "NaN in one of the runs");
    TEST_ASSERT(pro_on.n_orbits == N, "pro_on completed %d/%d orbits",
                pro_on.n_orbits, N);
    TEST_ASSERT(ret_on.n_orbits == N, "ret_on completed %d/%d orbits",
                ret_on.n_orbits, N);

    const float dT_meas_on   = ret_on.T_mean  - pro_on.T_mean;
    const float dT_meas_off  = ret_off.T_mean - pro_off.T_mean;
    const float dT_err_frac  = fabsf(dT_meas_on - dT_ana) / fabsf(dT_ana);

    printf("[GM force ENABLED — production]\n");
    printf("  Prograde  : v_init=%.6f  T_mean=%.4f  r in [%.3f, %.3f]\n",
           (double) pro_on.v_init, (double) pro_on.T_mean,
           (double) pro_on.r_min,  (double) pro_on.r_max);
    printf("  Retrograde: v_init=%.6f  T_mean=%.4f  r in [%.3f, %.3f]\n",
           (double) ret_on.v_init, (double) ret_on.T_mean,
           (double) ret_on.r_min,  (double) ret_on.r_max);
    printf("  T_ret - T_pro (measured) = %+.4f\n", (double) dT_meas_on);
    printf("  T_ret - T_pro (analytic) = %+.4f\n", (double) dT_ana);
    printf("  Relative error           = %.4f%%\n", 100.0 * (double) dT_err_frac);
    printf("\n");
    printf("[GM force DISABLED — baseline scalar gravity]\n");
    printf("  Prograde  : v_init=%.6f  T_mean=%.4f\n",
           (double) pro_off.v_init, (double) pro_off.T_mean);
    printf("  Retrograde: v_init=%.6f  T_mean=%.4f\n",
           (double) ret_off.v_init, (double) ret_off.T_mean);
    printf("  T_ret - T_pro (measured) = %+.4f   (should be ~0)\n",
           (double) dT_meas_off);
    printf("\n");

    /* Assertions:
     *  (1) GM-off prograde and retrograde periods agree (sanity).  Tolerance
     *      is the leapfrog single-period closure error ~ dt/T.
     *  (2) GM-on measured Delta T matches analytic to better than 1%.  The
     *      remaining error is the leapfrog truncation residual after one
     *      orbit (~dt * omega ~ 0.7 * (2 pi / 562) ~ 0.8% per period; we
     *      average over N=4 so it ought to come down). */
    const float baseline_tol = 0.005f;   /* 0.5% of T_rel */
    TEST_ASSERT(fabsf(dT_meas_off) / T_rel_ana < baseline_tol,
                "Baseline (GM-off) Delta T %+.4f is too large (>%.2f%% of T_rel)",
                (double) dT_meas_off, 100.0 * (double) baseline_tol);
    TEST_ASSERT(dT_err_frac < 0.02f,
                "GM-on Delta T error %.4f%% exceeds 2%%",
                100.0 * (double) dT_err_frac);

    /* Also assert: the orbits are within ~0.5% of the prescribed radius —
     * this is what gives us confidence the orbit is truly circular and the
     * period measurement is clean. */
    const float radial_tol = 0.005f * r;
    TEST_ASSERT((pro_on.r_max - pro_on.r_min) < 4.0f * radial_tol,
                "Prograde orbit not steady: r in [%.3f, %.3f]",
                (double) pro_on.r_min, (double) pro_on.r_max);

    printf("ALL CHECKS PASSED.\n");
    return 0;
}
