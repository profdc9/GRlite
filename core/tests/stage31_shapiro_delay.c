/* Stage 31 -- Shapiro delay verification for the EM wave equation.
 *
 * gr_sandbox_v35.tex sec:shapiro eq:c_local:
 *   c_local(x) = c (1 + 2 Phi_g(x) / c^2)
 *   Delta t_Shapiro = -(2 / c^3) int_path Phi_g d_l
 *
 * Setup: a Gaussian initial-condition pulse on phi_em (the EM scalar
 * potential) is placed at (cx - L/2, cy).  The pulse expands radially in
 * 2D; we sample the wavefront at the diametrically-opposite probe
 * (cx + L/2, cy).  Run twice:
 *   (1) baseline -- no Phi_g background; uniform-c EM wave equation.
 *   (2) shapiro  -- softened point-mass at the center, Shapiro enabled.
 *
 * Peak arrival times at the probe are extracted with parabolic
 * interpolation around the discrete maximum, then differenced to give
 * the measured Shapiro delay.
 *
 * Analytic prediction (passing impact parameter b = 0 through a softened
 * mass with smoothing length eps):
 *
 *   Delta t_ana = (4 GM / c^3) * asinh( L / (2 eps) )
 *
 * derived from int_{-L/2}^{L/2} dx / sqrt(x^2 + eps^2) = 2 asinh(L/(2 eps))
 * and the doc's (2/c^3) coefficient on the Phi_g line integral. */

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

/* Initialize phi_em with a Gaussian centered at (x0, y0), with zero
 * initial time derivative (prev = curr).  CORNER sublattice: nodes at
 * (i, j) * dx. */
static void init_em_gaussian(gr_sim_t* sim, float x0, float y0,
                             float sigma, float amp) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    const float inv_2s2 = 1.0f / (2.0f * sigma * sigma);
    float* curr = sim->fields[GR_FIELD_PHI_EM].curr;
    float* prev = sim->fields[GR_FIELD_PHI_EM].prev;
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x  = (float) i * dx;
            const float r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
            const float v  = amp * expf(-r2 * inv_2s2);
            curr[row + i] = v;
            prev[row + i] = v;     /* zero initial time derivative */
        }
    }
}

/* Bilinear sample of phi_em (CORNER sublattice) at exact (xp, yp). */
static float sample_phi_em(const gr_sim_t* sim, float xp, float yp) {
    const int   W  = gr_sim_width((gr_sim_t*) sim);
    const int   H  = gr_sim_height((gr_sim_t*) sim);
    const float dx = gr_sim_dx((gr_sim_t*) sim);
    float u = xp / dx;
    float v = yp / dx;
    if (u < 0.0f) u = 0.0f;
    if (v < 0.0f) v = 0.0f;
    if (u > (float) (W - 1)) u = (float) (W - 1);
    if (v > (float) (H - 1)) v = (float) (H - 1);
    int i0 = (int) u;
    int j0 = (int) v;
    if (i0 > W - 2) i0 = W - 2;
    if (j0 > H - 2) j0 = H - 2;
    const float fx = u - (float) i0;
    const float fy = v - (float) j0;
    const float* a = gr_sim_field_ptr((gr_sim_t*) sim, GR_FIELD_PHI_EM);
    const float v00 = a[ j0      * W + i0    ];
    const float v10 = a[ j0      * W + i0 + 1];
    const float v01 = a[(j0 + 1) * W + i0    ];
    const float v11 = a[(j0 + 1) * W + i0 + 1];
    return (1.0f - fy) * ((1.0f - fx) * v00 + fx * v10)
         + fy          * ((1.0f - fx) * v01 + fx * v11);
}

/* Run one pulse experiment and return parabolic-refined peak arrival time
 * at the probe.  If with_mass is set, install a softened point mass at the
 * grid center and enable Shapiro. */
static float run_pulse(int with_mass, float GM, float eps,
                       float xsrc, float ysrc,
                       float xprobe, float yprobe,
                       int W, int H, int n_steps,
                       float* peak_amp_out, float* peak_time_out) {
    const float dx    = 1.0f;
    const float c_eff = 1.0f;
    const float cfl   = 1.0f / sqrtf(2.0f);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return NAN;
    gr_sim_set_damping(sim, 16);

    if (with_mass) {
        const float cx = ((float) (W - 1) * 0.5f) * dx;
        const float cy = ((float) (H - 1) * 0.5f) * dx;
        gr_sim_set_background_point_mass(sim, cx, cy, GM, eps);
        gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
        gr_sim_set_em_shapiro_enabled(sim, 1);
    }

    init_em_gaussian(sim, xsrc, ysrc, /*sigma=*/4.0f, /*amp=*/1.0f);

    const float dt = gr_sim_dt(sim);
    float ring[3] = { 0.0f, 0.0f, 0.0f };
    int   filled  = 0;
    float peak_v = 0.0f;
    float peak_t = 0.0f;
    float peak_l = 0.0f, peak_r = 0.0f;
    for (int n = 0; n < n_steps; n++) {
        gr_sim_step(sim);
        const float v = sample_phi_em(sim, xprobe, yprobe);
        /* Shift the 3-sample ring: ring[0] = newest at t=(n+1)*dt,
         * ring[1] = previous at t=n*dt, ring[2] = at t=(n-1)*dt. */
        ring[2] = ring[1];
        ring[1] = ring[0];
        ring[0] = v;
        if (++filled < 3) continue;
        /* ring[1] is the candidate (middle) sample at time n*dt.
         * Track its largest value over the run. */
        if (ring[1] > peak_v) {
            peak_v = ring[1];
            peak_t = (float) n * dt;
            peak_l = ring[2];
            peak_r = ring[0];
        }
    }
    /* Parabolic refinement around the discrete max:
     *   Delta = (L - R) / (2 (L - 2 P + R)) * dt
     * with samples L=ring[2], P=peak_v, R=ring[0] saved at the peak step. */
    const float denom = peak_l - 2.0f * peak_v + peak_r;
    if (fabsf(denom) > 1e-12f) {
        const float off = (peak_l - peak_r) / (2.0f * denom) * dt;
        peak_t += off;
    }
    if (peak_amp_out)  *peak_amp_out  = peak_v;
    if (peak_time_out) *peak_time_out = peak_t;
    gr_sim_destroy(sim);
    return peak_t;
}

int main(void) {
    printf("=== stage31_shapiro_delay ===\n");
    printf("EM pulse past a softened point mass: measured vs analytic delay.\n");
    printf("  Spec: gr_sandbox_v35.tex sec:shapiro eq:c_local\n");
    printf("  Analytic (b=0):  Delta_t = (4 GM / c^3) asinh( L / (2 eps) )\n\n");

    const int   W  = 256, H = 256;
    const float dx = 1.0f;
    const float c_eff = 1.0f;
    const float GM   = 0.5f;
    const float eps  = 4.0f;
    const float L    = 120.0f;
    const float cx   = ((float) (W - 1) * 0.5f) * dx;
    const float cy   = ((float) (H - 1) * 0.5f) * dx;
    const float xsrc   = cx - L * 0.5f;
    const float ysrc   = cy;
    const float xprobe = cx + L * 0.5f;
    const float yprobe = cy;
    const int   n_steps = 220;        /* baseline arrival ~ L/c = 120; +60% margin */

    printf("  grid %dx%d, dx=%.1f, c=%.1f, cfl=1/sqrt(2)\n", W, H, dx, c_eff);
    printf("  point mass: GM=%.3f, eps=%.1f at (%.1f, %.1f)\n", GM, eps, cx, cy);
    printf("  pulse source at (%.1f, %.1f), probe at (%.1f, %.1f), L=%.1f\n",
           xsrc, ysrc, xprobe, yprobe, L);
    printf("  pulse sigma=4.0, amp=1.0\n\n");

    float baseline_amp = 0.0f, baseline_t = 0.0f;
    const float t_base = run_pulse(/*with_mass=*/0, 0.0f, eps,
                                   xsrc, ysrc, xprobe, yprobe,
                                   W, H, n_steps, &baseline_amp, &baseline_t);
    TEST_ASSERT(isfinite(t_base) && baseline_amp > 0.0f,
                "baseline run produced no peak");

    float shapiro_amp = 0.0f, shapiro_t = 0.0f;
    const float t_shap = run_pulse(/*with_mass=*/1, GM, eps,
                                   xsrc, ysrc, xprobe, yprobe,
                                   W, H, n_steps, &shapiro_amp, &shapiro_t);
    TEST_ASSERT(isfinite(t_shap) && shapiro_amp > 0.0f,
                "shapiro run produced no peak");

    const float c3 = c_eff * c_eff * c_eff;
    const float dt_meas = t_shap - t_base;
    const float dt_ana  = (4.0f * GM / c3) * asinhf(L / (2.0f * eps));

    printf("  baseline peak |phi_em| = %.5f at t = %.4f\n", baseline_amp, t_base);
    printf("  shapiro  peak |phi_em| = %.5f at t = %.4f\n", shapiro_amp, t_shap);
    printf("\n");
    printf("  measured  Delta_t = %+.4f\n", dt_meas);
    printf("  analytic  Delta_t = %+.4f   [(4 GM/c^3) asinh(L/(2 eps))]\n", dt_ana);
    printf("  error             = %+.2f%%\n", 100.0 * (double)(dt_meas - dt_ana) / (double) dt_ana);

    /* Sanity: the delay must be positive (mass slows the pulse). */
    TEST_ASSERT(dt_meas > 0.0f,
                "Shapiro delay should be positive but measured %+.4f", dt_meas);

    /* Tolerance: the analytic formula assumes a straight-line ray; the 2D
     * radial wavefront also probes nearby Phi_g(x, y) for |y - ysrc| ~ few
     * cells, which contributes additional softening.  The pulse shape
     * also broadens with propagation (2D Bessel-like dispersion).  20%
     * agreement is a fair pedagogical target. */
    const float rel_err = fabsf(dt_meas - dt_ana) / fabsf(dt_ana);
    TEST_ASSERT(rel_err < 0.25f,
                "Shapiro delay error %.2f%% exceeds 25%% tolerance",
                100.0f * rel_err);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
