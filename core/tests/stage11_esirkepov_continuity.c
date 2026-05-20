/* Stage 11 — gr_sandbox_vNN.tex §sec:yee_migration_plan answer B2.
 *
 * Verifies that Esirkepov current deposition (deposit.c) gives exact
 * discrete continuity
 *   (rho^n - rho^{n-1}) / dt + (J_x^{n-1/2}[i, j] - J_x^{n-1/2}[i-1, j]) / dx
 *                            + (J_y^{n-1/2}[i, j] - J_y^{n-1/2}[i, j-1]) / dx = 0
 * cell-by-cell, while direct CIC deposition (rho deposited at x^n, J at x^n
 * with the instantaneous velocity) violates continuity by exactly v . grad rho.
 *
 * Two trajectories:
 *   1. Straight line: rigid motion across several CFL-bounded steps.
 *   2. Circular: angular motion; J_x and J_y both change non-trivially.
 *
 * For each, we directly call gr_cic_deposit_corner (rho) and
 * gr_esirkepov_deposit_jxy (J) — no auto-deposit, no leapfrog, no particle
 * pusher.  Pure unit test of the discrete continuity property. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* Compute max |drho/dt + div J| over the interior of a W x H corner grid.
 * Corners at index (i, j); J_x at storage (i, j) is the x-edge between
 * corners (i, j) and (i+1, j); J_y at (i, j) is the y-edge between corners
 * (i, j) and (i, j+1).  Hence the divergence stencil reads J_x[i, j] and
 * J_x[i-1, j] (the edges to the right and left of corner (i, j)), and
 * similarly for J_y. */
static float max_continuity_residual(
        const float* rho_old, const float* rho_new,
        const float* Jx,      const float* Jy,
        int W, int H, float dx, float dt) {
    float maxres = 0.0f;
    for (int j = 1; j < H - 1; j++) {
        for (int i = 1; i < W - 1; i++) {
            const float drho_dt = (rho_new[j * W + i] - rho_old[j * W + i]) / dt;
            const float divJx   = (Jx[j * W + i] - Jx[j * W + i - 1]) / dx;
            const float divJy   = (Jy[j * W + i] - Jy[(j - 1) * W + i]) / dx;
            const float r = drho_dt + divJx + divJy;
            const float ar = (r < 0.0f) ? -r : r;
            if (ar > maxres) maxres = ar;
        }
    }
    return maxres;
}

/* Straight-line trajectory: a particle moves rigidly from x0 to x1 = x0 + v*dt
 * in a single step. */
static int test_straight_line(void) {
    printf("\n[1/3] Straight-line trajectory — single step continuity\n");
    const int   W  = 32, H = 32;
    const float dx = 1.0f;
    const float dt = 0.5f;
    const float m  = 1.0f;

    float* rho_old = (float*) calloc((size_t) W * H, sizeof(float));
    float* rho_new = (float*) calloc((size_t) W * H, sizeof(float));
    float* Jx      = (float*) calloc((size_t) W * H, sizeof(float));
    float* Jy      = (float*) calloc((size_t) W * H, sizeof(float));
    TEST_ASSERT(rho_old && rho_new && Jx && Jy, "alloc failed");

    /* Several positions / velocities to exercise the various sub-cell cases. */
    const struct { float x0, y0, vx, vy; const char* label; } cases[] = {
        {16.3f, 16.7f, 0.4f,  0.2f, "in-cell motion"},
        {16.0f, 16.0f, 0.4f,  0.2f, "from integer corner"},
        {16.5f, 16.5f, 0.4f,  0.2f, "from cell midpoint"},
        {16.8f, 16.8f, 0.5f,  0.5f, "crossing both x and y boundaries"},
        {16.3f, 16.7f,-0.4f, -0.2f, "negative-velocity in-cell"},
    };
    const int n_cases = sizeof(cases) / sizeof(cases[0]);

    float worst = 0.0f;
    for (int k = 0; k < n_cases; k++) {
        memset(rho_old, 0, (size_t) W * H * sizeof(float));
        memset(rho_new, 0, (size_t) W * H * sizeof(float));
        memset(Jx,      0, (size_t) W * H * sizeof(float));
        memset(Jy,      0, (size_t) W * H * sizeof(float));
        const float x0 = cases[k].x0, y0 = cases[k].y0;
        const float x1 = x0 + cases[k].vx * dt;
        const float y1 = y0 + cases[k].vy * dt;
        gr_cic_deposit_corner(rho_old, W, H, dx, x0, y0, m);
        gr_cic_deposit_corner(rho_new, W, H, dx, x1, y1, m);
        const int ok = gr_esirkepov_deposit_jxy(Jx, Jy, W, H, dx, dt,
                                                x0, y0, x1, y1, m);
        TEST_ASSERT(ok == 1, "case '%s': Esirkepov returned 0", cases[k].label);
        const float r = max_continuity_residual(rho_old, rho_new, Jx, Jy,
                                                W, H, dx, dt);
        printf("    %-50s max|residual| = %.3e\n", cases[k].label, r);
        if (r > worst) worst = r;
    }
    printf("  worst single-step residual = %.3e\n", worst);
    TEST_ASSERT(worst < 1e-5f,
                "max single-step continuity residual %.3e exceeds 1e-5", worst);

    free(rho_old); free(rho_new); free(Jx); free(Jy);
    return 0;
}

/* Circular trajectory: particle marches around a small circle, hitting all
 * eight octants.  Verify per-step continuity holds throughout. */
static int test_circular(void) {
    printf("\n[2/3] Circular trajectory — per-step continuity\n");
    const int   W  = 64, H = 64;
    const float dx = 1.0f;
    const float dt = 0.5f;
    const float m  = 1.0f;
    const float cx = (float) W * 0.5f * dx;
    const float cy = (float) H * 0.5f * dx;
    const float r0 = 10.0f * dx;
    const int   N_steps = 64;  /* full loop split into 64 small chord steps */

    float* rho_old = (float*) calloc((size_t) W * H, sizeof(float));
    float* rho_new = (float*) calloc((size_t) W * H, sizeof(float));
    float* Jx      = (float*) calloc((size_t) W * H, sizeof(float));
    float* Jy      = (float*) calloc((size_t) W * H, sizeof(float));
    TEST_ASSERT(rho_old && rho_new && Jx && Jy, "alloc failed");

    float max_step_residual = 0.0f;
    for (int s = 0; s < N_steps; s++) {
        const float th0 = 2.0f * (float) M_PI * (float)  s      / (float) N_steps;
        const float th1 = 2.0f * (float) M_PI * (float) (s + 1) / (float) N_steps;
        const float x0 = cx + r0 * cosf(th0);
        const float y0 = cy + r0 * sinf(th0);
        const float x1 = cx + r0 * cosf(th1);
        const float y1 = cy + r0 * sinf(th1);
        memset(rho_old, 0, (size_t) W * H * sizeof(float));
        memset(rho_new, 0, (size_t) W * H * sizeof(float));
        memset(Jx,      0, (size_t) W * H * sizeof(float));
        memset(Jy,      0, (size_t) W * H * sizeof(float));
        gr_cic_deposit_corner(rho_old, W, H, dx, x0, y0, m);
        gr_cic_deposit_corner(rho_new, W, H, dx, x1, y1, m);
        const int ok = gr_esirkepov_deposit_jxy(Jx, Jy, W, H, dx, dt,
                                                x0, y0, x1, y1, m);
        TEST_ASSERT(ok == 1, "step %d: Esirkepov returned 0", s);
        const float r = max_continuity_residual(rho_old, rho_new, Jx, Jy,
                                                W, H, dx, dt);
        if (r > max_step_residual) max_step_residual = r;
    }
    printf("  worst per-step residual over %d circle steps = %.3e\n",
           N_steps, max_step_residual);
    TEST_ASSERT(max_step_residual < 1e-5f,
                "max circular continuity residual %.3e exceeds 1e-5",
                max_step_residual);

    free(rho_old); free(rho_new); free(Jx); free(Jy);
    return 0;
}

/* Direct CIC (no Esirkepov) violates continuity by ~|v . grad rho|.  This
 * regression check ensures that disabling Esirkepov measurably affects the
 * discrete continuity property — so a future "did we accidentally bypass
 * Esirkepov?" question can be answered empirically. */
static int test_cic_violates(void) {
    printf("\n[3/3] Direct CIC violates discrete continuity (control)\n");
    const int   W  = 32, H = 32;
    const float dx = 1.0f;
    const float dt = 0.5f;
    const float m  = 1.0f;
    const float x0 = 16.3f, y0 = 16.7f;
    const float vx = 0.4f,  vy = 0.2f;
    const float x1 = x0 + vx * dt;
    const float y1 = y0 + vy * dt;

    float* rho_old = (float*) calloc((size_t) W * H, sizeof(float));
    float* rho_new = (float*) calloc((size_t) W * H, sizeof(float));
    float* Jx      = (float*) calloc((size_t) W * H, sizeof(float));
    float* Jy      = (float*) calloc((size_t) W * H, sizeof(float));
    TEST_ASSERT(rho_old && rho_new && Jx && Jy, "alloc failed");

    gr_cic_deposit_corner(rho_old, W, H, dx, x0, y0, m);
    gr_cic_deposit_corner(rho_new, W, H, dx, x1, y1, m);

    /* Direct CIC deposit: rho * v at the current particle position. */
    gr_cic_deposit_xedge(Jx, W, H, dx, x1, y1, m * vx);
    gr_cic_deposit_yedge(Jy, W, H, dx, x1, y1, m * vy);

    const float cic_residual = max_continuity_residual(rho_old, rho_new,
                                                       Jx, Jy, W, H, dx, dt);
    printf("  direct-CIC max|residual| = %.3e (expected > 1e-3, well above Esirkepov)\n",
           cic_residual);
    TEST_ASSERT(cic_residual > 1e-3f,
                "direct-CIC residual %.3e unexpectedly small — control test broken",
                cic_residual);

    /* Same case via Esirkepov for comparison. */
    memset(Jx, 0, (size_t) W * H * sizeof(float));
    memset(Jy, 0, (size_t) W * H * sizeof(float));
    gr_esirkepov_deposit_jxy(Jx, Jy, W, H, dx, dt, x0, y0, x1, y1, m);
    const float esirkepov_residual = max_continuity_residual(rho_old, rho_new,
                                                             Jx, Jy, W, H, dx, dt);
    printf("  Esirkepov max|residual|  = %.3e (should be float-precision small)\n",
           esirkepov_residual);
    printf("  ratio Esirkepov/CIC      = %.3e (lower is better, smaller violation)\n",
           esirkepov_residual / cic_residual);
    TEST_ASSERT(esirkepov_residual < cic_residual * 1e-3f,
                "Esirkepov should be at least 1000x better than direct CIC "
                "(got %.3e vs %.3e)",
                esirkepov_residual, cic_residual);

    free(rho_old); free(rho_new); free(Jx); free(Jy);
    return 0;
}

int main(void) {
    printf("=== stage11_esirkepov_continuity: discrete charge conservation ===\n");
    if (test_straight_line() != 0) return 1;
    if (test_circular()      != 0) return 1;
    if (test_cic_violates()  != 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
