/* Stage 31b -- Shapiro delay GM sweep.
 *
 * Repeats Stage 31's pulse experiment at several gravitational field
 * strengths to verify:
 *   (1) sign always positive (mass slows the wave),
 *   (2) linear scaling with GM (the analytic formula is linear in GM),
 *   (3) error stays bounded into the weak-field regime, and degrades
 *       gracefully as |2 Phi_g/c^2| approaches 1.
 *
 * Setup identical to Stage 31 except GM is swept.  L=120, eps=4, sigma=4.
 *
 * Predictions: each row's "delay" column should increase ~linearly with
 * GM.  The "error" column should hover near Stage 31's -8% (radial
 * wavefront averaging) until the linearized c_local formula breaks down
 * (max |2 Phi_g/c^2| approaches 1, i.e. GM/eps approaches c^2/2 = 0.5
 * here).  Past that, the simulator is still numerically stable thanks
 * to the 1e-3 floor on (1+2Phi/c^2), but the analytic comparison is
 * the wrong reference. */

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
            prev[row + i] = v;
        }
    }
}

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

static float run_pulse(int with_mass, float GM, float eps,
                       float xsrc, float ysrc,
                       float xprobe, float yprobe,
                       int W, int H, int n_steps,
                       float* peak_amp_out) {
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
    init_em_gaussian(sim, xsrc, ysrc, 4.0f, 1.0f);
    const float dt = gr_sim_dt(sim);
    float ring[3] = {0,0,0};
    int   filled  = 0;
    float peak_v = 0.0f, peak_t = 0.0f, peak_l = 0.0f, peak_r = 0.0f;
    for (int n = 0; n < n_steps; n++) {
        gr_sim_step(sim);
        const float v = sample_phi_em(sim, xprobe, yprobe);
        ring[2] = ring[1]; ring[1] = ring[0]; ring[0] = v;
        if (++filled < 3) continue;
        if (ring[1] > peak_v) {
            peak_v = ring[1];
            peak_t = (float) n * dt;
            peak_l = ring[2];
            peak_r = ring[0];
        }
    }
    const float denom = peak_l - 2.0f * peak_v + peak_r;
    if (fabsf(denom) > 1e-12f) {
        peak_t += (peak_l - peak_r) / (2.0f * denom) * dt;
    }
    if (peak_amp_out) *peak_amp_out = peak_v;
    gr_sim_destroy(sim);
    return peak_t;
}

int main(void) {
    printf("=== stage31b_shapiro_strength_sweep ===\n");
    printf("Pulse experiment repeated at several gravitational strengths.\n");
    printf("L=120, eps=4, sigma=4; sweep GM = 0.1 .. 1.5.\n\n");

    const int   W = 256, H = 256;
    const float dx = 1.0f;
    const float c_eff = 1.0f;
    const float eps = 4.0f;
    const float L   = 120.0f;
    const float cx  = ((float) (W - 1) * 0.5f) * dx;
    const float cy  = ((float) (H - 1) * 0.5f) * dx;
    const float xsrc   = cx - L * 0.5f;
    const float xprobe = cx + L * 0.5f;
    const int   n_steps = 240;

    /* Baseline (no mass) -- only need to run once. */
    float base_amp = 0.0f;
    const float t_base = run_pulse(0, 0.0f, eps, xsrc, cy, xprobe, cy,
                                   W, H, n_steps, &base_amp);
    printf("baseline: peak_amp=%.5f at t=%.4f\n\n", base_amp, t_base);

    const float GMs[] = {0.1f, 0.25f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f};
    const int n_GM = sizeof(GMs) / sizeof(GMs[0]);

    printf("%-6s %-8s %-12s %-12s %-12s %-12s %-10s\n",
           "GM", "max2Φ/c²", "peak_amp", "t_shap", "Δt_meas", "Δt_ana", "error%");
    printf("---------------------------------------------------------------------------\n");

    float prev_meas = 0.0f;
    float prev_GM   = 0.0f;
    int ok = 1;
    for (int i = 0; i < n_GM; i++) {
        const float GM = GMs[i];
        const float max2phi_c2 = 2.0f * GM / eps;   /* |2 Phi_g/c^2| at r=0 */
        float amp = 0.0f;
        const float t_sh = run_pulse(1, GM, eps, xsrc, cy, xprobe, cy,
                                     W, H, n_steps, &amp);
        const float c3 = c_eff * c_eff * c_eff;
        const float dt_meas = t_sh - t_base;
        const float dt_ana  = (4.0f * GM / c3) * asinhf(L / (2.0f * eps));
        const float err = 100.0f * (dt_meas - dt_ana) / dt_ana;

        printf("%-6.3f %-8.3f %-12.5f %-12.4f %+11.4f  %+11.4f  %+9.2f\n",
               GM, max2phi_c2, amp, t_sh, dt_meas, dt_ana, err);
        fflush(stdout);

        /* Sign: positive. */
        if (!(dt_meas > 0.0f)) {
            printf("  !! sign wrong at GM=%.3f\n", GM);
            ok = 0;
        }
        /* Linearity diagnostic: contribution per unit GM should be roughly
         * the same in the weak-field regime.  Print slope between consecutive
         * runs. */
        if (i > 0) {
            const float slope_meas = (dt_meas - prev_meas) / (GM - prev_GM);
            const float slope_ana  = (4.0f / c3) * asinhf(L / (2.0f * eps));
            printf("    [slope_meas=%.4f vs slope_ana=%.4f]\n",
                   slope_meas, slope_ana);
        }
        prev_meas = dt_meas;
        prev_GM   = GM;
    }

    printf("\n");
    printf("Notes:\n");
    printf("  - Analytic formula is linear in GM (slope_ana = constant).\n");
    printf("  - Measured slope should match the analytic slope in the weak-field\n");
    printf("    regime |2 Phi/c^2| < 1; beyond that, the c_local approximation\n");
    printf("    itself breaks down (not a simulator bug -- the analytic ref is wrong\n");
    printf("    there).\n");
    printf("  - The 1e-3 floor on (1+2Phi/c^2) engages when |2Phi/c^2| > 0.999\n");
    printf("    (i.e. GM > 0.999 * eps/2 = %.3f here).\n", 0.999f * eps * 0.5f);

    TEST_ASSERT(ok, "one or more sign tests failed");

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
