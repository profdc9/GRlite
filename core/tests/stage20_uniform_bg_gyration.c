/* Stage 20 — gravitomagnetic Lorentz force, unit isolation test.
 *
 * Sets up a spatially-uniform gravitomagnetic field B_g_z = B_0 (via the
 * symmetric-gauge background, see gr_sim_set_background_uniform_gravitomagnetic
 * in core/src/background.c) with Phi_g = 0 (no scalar gravity).  The force on
 * a test particle reduces to the pure Tier-1 gravitomagnetic Lorentz piece
 * (gr_sandbox_v35.tex eq:geodesic_expansion, line 938; algorithmic Tier-3
 * eqbox, line 1040):
 *
 *     F = m * 4 (v x B_g)        with    (v x B_g_z z)_x = +v_y B_g_z,
 *                                       (v x B_g_z z)_y = -v_x B_g_z.
 *
 * The motion is exactly cyclotron-like.  For initial velocity v_0 = (v_0, 0)
 * with B_g_z > 0, the particle gyrates CLOCKWISE (viewed from +z) at
 *
 *     gyrofrequency      omega_gm  = 4 |B_g_z| / gamma
 *     period             T         = 2 pi / omega_gm   =   pi / (2 B_0 / gamma)
 *     Larmor radius      r_L       = v_0 / omega_gm    =   v_0 gamma / (4 B_0)
 *     orbit center       (x0, y0 - r_L)
 *
 * In the non-relativistic limit (gamma ~= 1) these become T = pi/(2 B_0) and
 * r_L = v_0/(4 B_0).
 *
 * Checks performed at v_0 = 0.1 c, B_0 = 0.01 (so v/c = 0.1 -> gamma ~= 1.005,
 * Newtonian limit good to about 0.5 percent):
 *
 *   (1) After one full analytic period the particle returns within 1% of r_L
 *       of its starting position (closure test).
 *   (2) |v|^2 is conserved to better than 1% over one period (v x B does no
 *       work, modulo Boris-leapfrog corrector error).
 *   (3) Measured Larmor radius (max(|r - r_center|) over one period) matches
 *       analytic prediction to within 2%.
 *   (4) Reversing the sign of B_0 reverses the gyration direction (the
 *       particle moves to +y at the start of the orbit instead of -y).
 *
 * This is a fixed-background test: field evolution is disabled, no particle
 * sources are deposited, no perturbation A_g is present.  ANALYTIC bg mode
 * is used so B_g_z reads from the closed form (not a Yee curl). */

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
    float B0;
    float v0;
    float r_L;
    float T;
    float gamma;
    float closure_err_frac;    /* |r_end - r_start| / r_L                      */
    float v2_drift_frac;       /* |v^2_end - v^2_start| / v^2_start            */
    float r_L_measured;        /* max(|r - r_center|) over one orbit           */
    float r_L_err_frac;        /* |r_L_measured - r_L| / r_L                   */
    int   sign_B0;             /* +1 or -1                                     */
    int   gyration_sign;       /* +1 = clockwise (B0>0 expectation), -1 = ccw  */
} run_t;

static void run_cyclotron(int sign_B0, run_t* out) {
    /* Grid sized so the orbit fits well inside the interior.  r_L ~ 25
     * cells, centered at (W/2, H/2 - r_L), so 128x128 keeps the entire
     * orbit ~10 cells away from each edge.  Out-of-bounds CIC interp
     * returns 0 (matches the analytic Phi_g = 0 background), so the test
     * also works at smaller grid sizes — the bump is for visual clarity. */
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    /* B_0 chosen so dt * omega_gm << 1, otherwise the closure error is
     * dominated by rounding T_ana to integer steps (the leapfrog itself
     * has tiny truncation error here).  At B_0 = 1e-3, omega ~= 4e-3,
     * dt ~= 0.7, so dt * omega ~= 2.8e-3 -> well below 1%. */
    const float B0     = (float) sign_B0 * 1.0e-3f;
    const float v0     = 0.1f;
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;

    /* Analytic predictions (relativistic momentum -> use the Boris pusher's
     * lagged-momentum convention; v0 = p0 / (gamma m), so for a particle
     * launched with speed v0 the conserved 3-velocity magnitude is v0 and
     * gamma = 1/sqrt(1 - v0^2/c^2). */
    const float gamma = 1.0f / sqrtf(1.0f - (v0 * v0) / (c_eff * c_eff));
    const float omega = 4.0f * fabsf(B0) / gamma;
    const float T_ana = 2.0f * (float) M_PI / omega;
    const float r_L   = v0 / omega;     /* = v0 gamma / (4 |B0|) */

    /* Sim setup: no field evolution, no sources, no damping, analytic bg. */
    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "gr_sim_create failed");

    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_background_uniform_gravitomagnetic(sim, cx, cy, B0);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

    /* Particle at the gauge origin with v = (+v0, 0).  The relativistic
     * momentum is p = gamma m v.  gr_sim_add_particle takes vx, vy and
     * stores p = gamma m v internally, so we pass v0 directly. */
    const int idx = gr_sim_add_particle(sim, cx, cy, /*mass=*/1.0f,
                                        /*charge=*/0.0f, v0, 0.0f);
    TEST_ASSERT(idx == 0, "gr_sim_add_particle expected idx 0, got %d", idx);

    const float dt = gr_sim_dt(sim);
    const int   steps_per_period = (int) ceilf(T_ana / dt);

    /* For B0 > 0, v starts as (+v0, 0), so F_y = -4 m v_x B0 < 0,
     * pushing the particle DOWN -> clockwise gyration in (x, y) plane
     * (viewed from +z).  The orbit center sits at (cx, cy - r_L) for B0 > 0
     * and (cx, cy + r_L) for B0 < 0. */
    const float yc_orbit = (B0 > 0.0f) ? (cy - r_L) : (cy + r_L);

    /* Step and probe.  Track max radius from orbit center to estimate r_L. */
    const float v0_2 = v0 * v0;
    float r_max = 0.0f;
    float min_x = cx + 10.0f, max_x = cx - 10.0f;
    float min_y = cy + 10.0f, max_y = cy - 10.0f;
    /* Early gyration direction probe — sample after one short step. */
    gr_sim_step(sim);
    const gr_particle_t* p1 = gr_sim_get_particle(sim, 0);
    out->gyration_sign = (p1->py < 0.0f) ? +1 : -1;   /* +1 = clockwise */
    /* (gyration_sign agrees with the sign of B0 by the convention above:
     *  B0 > 0 -> initial F_y < 0 -> p_y goes negative -> "clockwise". */

    for (int s = 1; s < steps_per_period; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - yc_orbit;
        const float r  = sqrtf(rx * rx + ry * ry);
        if (r > r_max) r_max = r;
        if (p->x < min_x) min_x = p->x;
        if (p->x > max_x) max_x = p->x;
        if (p->y < min_y) min_y = p->y;
        if (p->y > max_y) max_y = p->y;
    }

    const gr_particle_t* p_end = gr_sim_get_particle(sim, 0);
    const float dx_end = p_end->x - cx;
    const float dy_end = p_end->y - cy;
    const float closure_err = sqrtf(dx_end * dx_end + dy_end * dy_end);

    /* Final |v|^2.  Reconstruct v from p via p = gamma m v. */
    const float pm2 = p_end->px * p_end->px + p_end->py * p_end->py;
    const float gamma_end = sqrtf(1.0f + pm2 / (p_end->mass * p_end->mass * c_eff * c_eff));
    const float vx_end = p_end->px / (gamma_end * p_end->mass);
    const float vy_end = p_end->py / (gamma_end * p_end->mass);
    const float v2_end = vx_end * vx_end + vy_end * vy_end;

    out->B0               = B0;
    out->v0               = v0;
    out->r_L              = r_L;
    out->T                = T_ana;
    out->gamma            = gamma;
    out->closure_err_frac = closure_err / r_L;
    out->v2_drift_frac    = fabsf(v2_end - v0_2) / v0_2;
    out->r_L_measured     = r_max;
    out->r_L_err_frac     = fabsf(r_max - r_L) / r_L;
    out->sign_B0          = sign_B0;

    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage20_uniform_bg_gyration ===\n");
    printf("Gravitomagnetic cyclotron, F = 4 m (v x B_g), B_g_z = constant.\n");
    printf("Spec: gr_sandbox_v35.tex eq:geodesic_expansion (line 938).\n\n");

    /* Run with B_0 > 0 (clockwise expected). */
    run_t r_pos, r_neg;
    run_cyclotron(+1, &r_pos);
    run_cyclotron(-1, &r_neg);

    printf("Analytic predictions (B_0 = +%g, v_0 = %g c):\n", (double) r_pos.B0, (double) r_pos.v0);
    printf("  gamma         = %.5f\n", (double) r_pos.gamma);
    printf("  omega_gm      = %.5f  (= 4 |B_0| / gamma)\n",
           (double) (4.0f * fabsf(r_pos.B0) / r_pos.gamma));
    printf("  period T      = %.5f\n", (double) r_pos.T);
    printf("  Larmor r_L    = %.5f  (= v_0 / omega_gm)\n", (double) r_pos.r_L);
    printf("\n");

    printf("B_0 = +%g:\n", (double) r_pos.B0);
    printf("  closure error (|r_end - r_start| / r_L)    = %.4f%%\n",
           100.0 * (double) r_pos.closure_err_frac);
    printf("  v^2 drift     (|v^2_end - v^2_0| / v^2_0)   = %.4f%%\n",
           100.0 * (double) r_pos.v2_drift_frac);
    printf("  r_L measured  (max |r - r_center|)          = %.5f  (err %.4f%%)\n",
           (double) r_pos.r_L_measured,
           100.0 * (double) r_pos.r_L_err_frac);
    printf("  gyration sign (+1=CW, -1=CCW)               = %+d\n", r_pos.gyration_sign);

    printf("\nB_0 = -%g (sign-reversal test):\n", fabs((double) r_neg.B0));
    printf("  gyration sign                               = %+d\n", r_neg.gyration_sign);

    /* Assertions. */
    TEST_ASSERT(r_pos.closure_err_frac < 0.01f,
                "B_0 > 0: closure error %.4f%% exceeds 1%%",
                100.0 * (double) r_pos.closure_err_frac);
    TEST_ASSERT(r_pos.v2_drift_frac < 0.01f,
                "B_0 > 0: |v|^2 drift %.4f%% exceeds 1%%",
                100.0 * (double) r_pos.v2_drift_frac);
    TEST_ASSERT(r_pos.r_L_err_frac < 0.02f,
                "B_0 > 0: Larmor radius error %.4f%% exceeds 2%%",
                100.0 * (double) r_pos.r_L_err_frac);
    TEST_ASSERT(r_pos.gyration_sign == +1,
                "B_0 > 0: expected clockwise (+1) gyration, got %+d",
                r_pos.gyration_sign);
    TEST_ASSERT(r_neg.gyration_sign == -1,
                "B_0 < 0: expected counter-clockwise (-1) gyration, got %+d",
                r_neg.gyration_sign);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
