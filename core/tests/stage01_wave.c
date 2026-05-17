/* Stage 1 test — gr_sandbox_v32.tex §12.1 "Stage 1: Scalar wave equation, free
 * propagation". Verifies:
 *   (a) wavefront arrival time at radius r matches r / c_eff to within a few dt
 *   (b) amplitude falls off as 1/sqrt(r) in 2D (§10.3 "The 2D Green's function")
 *   (c) the scheme is stable at cfl = 1/sqrt(2), unstable just above it
 *   (d) field remains finite (no NaN/Inf) at end of measurement window
 *
 * Note: the dispersion-relation check from §12.1 bullet 4 (against eq:dispersion_physical
 * in §9.4) is deferred — that requires a Fourier analysis of a propagating sinusoid
 * and will be added as a separate test once the harness supports it. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TEST_ASSERT(cond, fmt, ...)                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);           \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

static int test_propagation_and_falloff(void) {
    /* Domain & scenario chosen so the wavefront has time to reach the outermost
     * sample point (r=80) and the peak is captured before the reflection from
     * the far wall returns. With c=1, dx=1, cfl=1/sqrt(2), dt~=0.7071, sample at
     * r=80 peaks at step ~113; reflection from the +x wall returns at step ~248.
     * We measure for 180 steps — safely inside the no-reflection window. */
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);  /* 2D stability limit, §9.2 eq:cfl */
    const float sigma  = 4.0f * dx;
    const float amp    = 1.0f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "gr_sim_create returned NULL");

    float params[2] = {sigma, amp};
    int   rc        = gr_sim_load_scenario(sim, "wave_pulse", params, 2);
    TEST_ASSERT(rc == 0, "load_scenario(wave_pulse) returned %d", rc);

    const int   cx       = W / 2;
    const int   cy       = H / 2;
    const int   radii[3] = {20, 40, 80};
    const int   n_smp    = 3;
    const float dt       = gr_sim_dt(sim);

    float peak_amp[3] = {0};
    float peak_t[3]   = {0};

    const int max_steps = 180;
    for (int n = 0; n < max_steps; n++) {
        gr_sim_step(sim);
        const float  t   = gr_sim_time(sim);
        const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
        for (int s = 0; s < n_smp; s++) {
            const float v = fabsf(phi[cy * W + (cx + radii[s])]);
            if (v > peak_amp[s]) {
                peak_amp[s] = v;
                peak_t[s]   = t;
            }
        }
    }

    /* (a) Wavefront arrival times.
     * The peak of a Gaussian-launched 2D wavefront at radius r lags the leading
     * edge by ~sigma/c, so r/c_eff matches the peak time to within a few timesteps
     * (sigma=4, c=1 => ~5.66 dt lag). We allow 8 dt tolerance. */
    const float tol_steps = 8.0f;
    for (int s = 0; s < n_smp; s++) {
        const float r          = (float) radii[s] * dx;
        const float t_expected = r / c_eff;
        const float err        = fabsf(peak_t[s] - t_expected);
        printf("  r=%-3d : peak |phi|=%.5f at t=%6.3f, expected ~%6.3f (err = %4.2f dt)\n",
               radii[s], peak_amp[s], peak_t[s], t_expected, err / dt);
        TEST_ASSERT(err < tol_steps * dt,
                    "wavefront arrival at r=%d off by %.2f dt (tolerance %.1f dt)",
                    radii[s], err / dt, tol_steps);
    }

    /* (b) 1/sqrt(r) amplitude falloff in 2D — eq:gem_metric / §10.3.
     * peak(r) * sqrt(r) should be approximately constant for r >> sigma. */
    float A[3];
    for (int s = 0; s < n_smp; s++) A[s] = peak_amp[s] * sqrtf((float) radii[s] * dx);
    printf("  amp*sqrt(r): %.4f, %.4f, %.4f\n", A[0], A[1], A[2]);

    /* Inner radii (r=20 = 5 sigma) carry transient finite-r corrections; 25%
     * tolerance on the ratio across the three samples is the published-result
     * regime. */
    const float ratio01 = A[1] / A[0];
    const float ratio12 = A[2] / A[1];
    TEST_ASSERT(ratio01 > 0.75f && ratio01 < 1.25f,
                "1/sqrt(r) falloff violated r=20 -> r=40: ratio = %.3f", ratio01);
    TEST_ASSERT(ratio12 > 0.75f && ratio12 < 1.25f,
                "1/sqrt(r) falloff violated r=40 -> r=80: ratio = %.3f", ratio12);

    /* (d) No non-finite values anywhere in the field. */
    const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
    for (int k = 0; k < W * H; k++) {
        TEST_ASSERT(isfinite(phi[k]), "non-finite phi[%d] after %d steps", k, gr_sim_step_count(sim));
    }

    gr_sim_destroy(sim);
    return 0;
}

static int test_cfl_stability(void) {
    /* (c) Stability boundary check — gr_sandbox_v32.tex §9.2 eq:cfl.
     * In 2D, the leapfrog 5-point scheme requires cfl <= 1/sqrt(2). At
     * cfl = 1/sqrt(2) * 1.05 (~5% over the limit) the highest-frequency mode
     * grows exponentially and the field blows up within tens of steps. */
    const int   W     = 64, H = 64;
    const float dx    = 1.0f;
    const float c_eff = 1.0f;
    const float cfl_stable   = 1.0f / sqrtf(2.0f);
    const float cfl_unstable = cfl_stable * 1.05f;

    /* Stable run: should remain finite. */
    {
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl_stable);
        TEST_ASSERT(sim != NULL, "create at stable CFL returned NULL");
        float p[2] = {2.0f, 1.0f};
        TEST_ASSERT(gr_sim_load_scenario(sim, "wave_pulse", p, 2) == 0, "scenario load failed");
        gr_sim_step_n(sim, 200);
        const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
        float max_abs = 0.0f;
        for (int k = 0; k < W * H; k++) {
            const float v = fabsf(phi[k]);
            if (!isfinite(v)) { max_abs = INFINITY; break; }
            if (v > max_abs) max_abs = v;
        }
        printf("  CFL=%.4f (stable):   max|phi| after 200 steps = %.4f\n", cfl_stable, max_abs);
        TEST_ASSERT(isfinite(max_abs) && max_abs < 10.0f * 1.0f /* amplitude */,
                    "stable run grew unexpectedly: max|phi| = %g", max_abs);
        gr_sim_destroy(sim);
    }

    /* Unstable run: should blow up. */
    {
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl_unstable);
        TEST_ASSERT(sim != NULL, "create at unstable CFL returned NULL");
        float p[2] = {2.0f, 1.0f};
        TEST_ASSERT(gr_sim_load_scenario(sim, "wave_pulse", p, 2) == 0, "scenario load failed");
        int blew_up_at = -1;
        for (int n = 0; n < 1000; n++) {
            gr_sim_step(sim);
            const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
            float max_abs = 0.0f;
            for (int k = 0; k < W * H; k++) {
                const float v = fabsf(phi[k]);
                if (!isfinite(v) || v > max_abs) max_abs = v;
            }
            if (!isfinite(max_abs) || max_abs > 1e6f) { blew_up_at = n + 1; break; }
        }
        printf("  CFL=%.4f (unstable): field blew up at step %d (max=1000 budget)\n",
               cfl_unstable, blew_up_at);
        TEST_ASSERT(blew_up_at > 0, "unstable CFL did NOT blow up — instability mechanism missing");
        gr_sim_destroy(sim);
    }
    return 0;
}

int main(void) {
    printf("=== stage01_wave: gr_sandbox_v32.tex §12.1 ===\n");

    printf("\n[1/2] propagation & 1/sqrt(r) falloff\n");
    if (test_propagation_and_falloff() != 0) return 1;

    printf("\n[2/2] CFL stability boundary\n");
    if (test_cfl_stability() != 0) return 1;

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
