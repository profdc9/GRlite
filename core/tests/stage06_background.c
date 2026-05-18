/* Stage 6 test — gr_sandbox_v33.tex §12.6 "Stage 6: Sampled background field
 * arrays". Verifies the softened-point-mass generator (eq:bg_softened_point_mass)
 * produces the expected sampled values and gradients, that the perturbation
 * FDTD is unaffected by an installed background, and that the damping layer
 * does not touch background arrays. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(cond, fmt, ...)                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__);           \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

static float phi_at_cell(const float* phi, int W, int i, int j) {
    return phi[j * W + i];
}

static int test_softened_point_mass_sampling(void) {
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float GM     = 5.0f;
    const float eps    = 4.0f * dx;
    const float x0     = ((float) W * 0.5f) * dx;
    const float y0     = ((float) H * 0.5f) * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "gr_sim_create returned NULL");

    /* Background pointer is NULL before any set_background_* call. */
    TEST_ASSERT(gr_sim_background_ptr(sim, GR_FIELD_PHI_GRAV) == NULL,
                "background ptr should be NULL before set");

    gr_sim_set_background_point_mass(sim, x0, y0, GM, eps);
    const float* phi_bg = gr_sim_background_ptr(sim, GR_FIELD_PHI_GRAV);
    TEST_ASSERT(phi_bg != NULL, "background ptr is NULL after set");

    /* (a) Sampled values match the analytic formula to float precision.
     * Sample at center, 4eps away on +x, and at a few off-axis points. */
    const struct { int i, j; } pts[] = {
        {W / 2, H / 2}, {W / 2 + 16, H / 2}, {W / 2 + 32, H / 2},
        {W / 2 + 16, H / 2 + 16}, {W / 2 + 24, H / 2 - 24},
    };
    const int n_pts = sizeof(pts) / sizeof(pts[0]);
    for (int k = 0; k < n_pts; k++) {
        const int i = pts[k].i, j = pts[k].j;
        const float x = ((float) i + 0.5f) * dx;
        const float y = ((float) j + 0.5f) * dx;
        const float r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        const float expected = -GM / sqrtf(r2 + eps * eps);
        const float observed = phi_at_cell(phi_bg, W, i, j);
        const float err = fabsf(observed - expected);
        const float rel = err / fabsf(expected);
        TEST_ASSERT(rel < 1.0e-5f,
                    "cell (%d,%d): observed=%.6g expected=%.6g rel.err=%.2e",
                    i, j, observed, expected, rel);
    }
    printf("  sampling matches analytic at %d sample points (rel.err < 1e-5)\n", n_pts);

    /* (b) Finite-difference field component g_x = -dphi/dx at r >> eps matches
     * the analytic signed gravitational field
     *     g_x = -G*M * (x - x0) / ( (x-x0)^2 + (y-y0)^2 + eps^2 )^{3/2}
     * which at our +x-axis sample points in the -x direction (toward the
     * source). At r/eps = 6 this should be within ~0.1% of the Newton limit
     * -G*M / r^2 and the centered finite-difference should reproduce it to
     * within the O((dx)^2) truncation error of the FD stencil. */
    const int   r_cells = 24;  /* r/eps = 6 */
    const int   icx = W / 2 + r_cells, jcy = H / 2;
    const float xs = ((float) icx + 0.5f) * dx;
    const float ys = ((float) jcy + 0.5f) * dx;
    const float dxr = xs - x0;
    const float dyr = ys - y0;
    const float r_actual2 = dxr * dxr + dyr * dyr;
    const float gx_observed = -(phi_at_cell(phi_bg, W, icx + 1, jcy)
                              - phi_at_cell(phi_bg, W, icx - 1, jcy)) / (2.0f * dx);
    const float gx_analytic = -GM * dxr / powf(r_actual2 + eps * eps, 1.5f);
    const float gx_newton   = -GM * dxr / powf(r_actual2, 1.5f);
    printf("  g_x at sample (i=%d,j=%d, r=%.2f):  observed=%.6f  softened=%.6f  Newton=%.6f\n",
           icx, jcy, sqrtf(r_actual2), gx_observed, gx_analytic, gx_newton);
    const float fd_rel = fabsf(gx_observed - gx_analytic) / fabsf(gx_analytic);
    TEST_ASSERT(fd_rel < 5.0e-3f, "FD-vs-softened-analytic rel.err %.3e exceeds 5e-3", fd_rel);
    const float newton_rel = fabsf(gx_analytic - gx_newton) / fabsf(gx_newton);
    TEST_ASSERT(newton_rel < 5.0e-2f,
                "softened differs from Newton by %.3e at r/eps=%g — smoothing too coarse",
                newton_rel, sqrtf(r_actual2) / eps);

    /* (c) Symmetry: phi_bg(W/2 + dr, H/2) == phi_bg(W/2 - dr, H/2) to fp prec. */
    for (int dr = 1; dr <= 30; dr += 5) {
        const float left  = phi_at_cell(phi_bg, W, W / 2 - dr, H / 2);
        const float right = phi_at_cell(phi_bg, W, W / 2 + dr - 1, H / 2);
        TEST_ASSERT(fabsf(left - right) < 1.0e-6f,
                    "asymmetric at dr=%d: left=%.6g right=%.6g", dr, left, right);
    }
    printf("  symmetry across x-axis verified for dr in [1,30]\n");

    gr_sim_destroy(sim);
    return 0;
}

static int test_perturbation_unaffected_by_background(void) {
    /* Background must not perturb the leapfrog. Run wave_pulse with and
     * without a background installed; the perturbation arrays must agree
     * bit-for-bit after the same number of steps. */
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const int   n_steps = 100;
    float       p[2]    = {4.0f, 1.0f};

    gr_sim_t* a = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_t* b = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(a && b, "create failed");

    /* Install a non-trivial background in `a` only. */
    gr_sim_set_background_point_mass(a, (W * 0.5f) * dx, (H * 0.5f) * dx, 5.0f, 4.0f);

    TEST_ASSERT(gr_sim_load_scenario(a, "wave_pulse", p, 2) == 0, "load a failed");
    TEST_ASSERT(gr_sim_load_scenario(b, "wave_pulse", p, 2) == 0, "load b failed");

    gr_sim_step_n(a, n_steps);
    gr_sim_step_n(b, n_steps);

    const float* phi_a = gr_sim_field_ptr(a, GR_FIELD_PHI_GRAV);
    const float* phi_b = gr_sim_field_ptr(b, GR_FIELD_PHI_GRAV);
    int diffs = 0;
    float max_abs_diff = 0.0f;
    for (int k = 0; k < W * H; k++) {
        const float d = fabsf(phi_a[k] - phi_b[k]);
        if (d > max_abs_diff) max_abs_diff = d;
        if (d != 0.0f) diffs++;
    }
    printf("  with-bg vs no-bg perturbation: %d differing cells, max|diff| = %.2e\n",
           diffs, max_abs_diff);
    TEST_ASSERT(diffs == 0, "perturbation differs in %d cells — background is leaking", diffs);

    gr_sim_destroy(a);
    gr_sim_destroy(b);
    return 0;
}

static int test_damping_does_not_touch_background(void) {
    /* Install a background, snapshot it, enable damping, step many times,
     * verify the background array is byte-identical at the end. */
    const int   W = 128, H = 128;
    gr_sim_t* sim = gr_sim_create(W, H, 1.0f, 1.0f, 1.0f / sqrtf(2.0f));
    TEST_ASSERT(sim != NULL, "create failed");

    gr_sim_set_background_point_mass(sim, 64.0f, 64.0f, 5.0f, 4.0f);
    const float* phi_bg = gr_sim_background_ptr(sim, GR_FIELD_PHI_GRAV);

    /* Snapshot. */
    float* snap = (float*) malloc((size_t) W * H * sizeof(float));
    TEST_ASSERT(snap != NULL, "snapshot alloc failed");
    memcpy(snap, phi_bg, (size_t) W * H * sizeof(float));

    /* Enable damping and step. */
    gr_sim_set_damping(sim, 16);
    float p[2] = {4.0f, 1.0f};
    TEST_ASSERT(gr_sim_load_scenario(sim, "wave_pulse", p, 2) == 0, "scenario failed");
    gr_sim_step_n(sim, 200);

    /* Background array must be unchanged. */
    const int cmp = memcmp(snap, phi_bg, (size_t) W * H * sizeof(float));
    TEST_ASSERT(cmp == 0, "background array was modified during stepping with damping");
    printf("  background array unchanged after 200 steps with damping enabled\n");

    free(snap);
    gr_sim_destroy(sim);
    return 0;
}

int main(void) {
    printf("=== stage06_background: gr_sandbox_v33.tex §12.6 ===\n");

    printf("\n[1/3] softened point-mass sampling and gradient\n");
    if (test_softened_point_mass_sampling() != 0) return 1;

    printf("\n[2/3] perturbation FDTD unaffected by installed background\n");
    if (test_perturbation_unaffected_by_background() != 0) return 1;

    printf("\n[3/3] damping does not touch background array\n");
    if (test_damping_does_not_touch_background() != 0) return 1;

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
