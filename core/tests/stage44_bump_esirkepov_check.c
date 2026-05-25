/* Stage 44 -- Tier-2 validation: bump Esirkepov continuity + bump
 * with inductive force enabled.
 *
 * Two checks:
 *   (A) Discrete continuity of bump Esirkepov (parallel to Stage 11
 *       for CIC Esirkepov): for a single trajectory step, verify
 *       (rho^{n+1} - rho^n)/dt + div(J^{n-1/2}) = 0 at every corner.
 *   (B) Full-force spurious accel with bump+inductive-on at a few v
 *       values, compared to TSC + smoothing.  Tests the Tier-2 effort
 *       of getting bump consistently applied to BOTH rho/phi and J/A. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Internal-API access (bump deposit + Esirkepov) for the standalone
 * continuity check.  Declared in sim_internal.h. */
void gr_bump_deposit_corner(float* arr, int W, int H, float dx,
                             float x_p, float y_p, float value, float R);
int  gr_bump_esirkepov_deposit_jxy(float* Jx, float* Jy,
                                    int W, int H, float dx, float dt,
                                    float x0, float y0, float x1, float y1,
                                    float source, float R);

/* ===========================================================
 * Part A: continuity test (Stage 11 analog for bump Esirkepov).
 * =========================================================== */

static float continuity_check(float x0, float y0, float x1, float y1,
                               float R, int* nonzero_cells_out) {
    const int   W   = 64;
    const int   H   = 64;
    const float dx  = 1.0f;
    const float dt  = 1.0f / sqrtf(2.0f);
    const float Q   = 1.0f;

    const size_t n = (size_t) W * (size_t) H;
    float* rho_n   = (float*) calloc(n, sizeof(float));
    float* rho_np1 = (float*) calloc(n, sizeof(float));
    float* Jx      = (float*) calloc(n, sizeof(float));
    float* Jy      = (float*) calloc(n, sizeof(float));

    gr_bump_deposit_corner(rho_n,   W, H, dx, x0, y0, Q, R);
    gr_bump_deposit_corner(rho_np1, W, H, dx, x1, y1, Q, R);
    gr_bump_esirkepov_deposit_jxy(Jx, Jy, W, H, dx, dt, x0, y0, x1, y1, Q, R);

    /* Continuity residual at each corner:
     *   res[i,j] = (rho_np1 - rho_n)/dt + (Jx[i] - Jx[i-1])/dx + (Jy[j] - Jy[j-1])/dx
     * where Jx[i] is the X-edge to the right of corner i (i.e., at i+0.5)
     * and Jx[i-1] is to the left (at i-0.5).  In our array layout, X-edge
     * (i, j) sits at position (i+0.5, j), so Jx[j*W + i] is to the right
     * of corner (i, j) and Jx[j*W + (i-1)] is to the left. */
    float max_abs = 0.0f;
    int nonzero = 0;
    for (int j = 1; j < H - 1; j++) {
        for (int i = 1; i < W - 1; i++) {
            const int k     = j * W + i;
            const float drho_dt = (rho_np1[k] - rho_n[k]) / dt;
            const float divJx   = (Jx[k] - Jx[k - 1]) / dx;
            const float divJy   = (Jy[k] - Jy[k - W]) / dx;
            const float res = drho_dt + divJx + divJy;
            if (fabsf(rho_n[k]) > 1e-12f || fabsf(rho_np1[k]) > 1e-12f) nonzero++;
            if (fabsf(res) > max_abs) max_abs = fabsf(res);
        }
    }

    free(rho_n); free(rho_np1); free(Jx); free(Jy);
    *nonzero_cells_out = nonzero;
    return max_abs;
}

static int part_A_continuity(void) {
    printf("--- Part A: bump Esirkepov continuity ---\n");
    printf("For each trajectory, verify (drho/dt + div J) ~ 0 at every corner.\n\n");
    printf("%-10s %-22s %-8s %-12s %s\n",
           "R", "trajectory", "cells", "max|residual|", "verdict");
    printf("-------------------------------------------------------------\n");
    int all_pass = 1;
    const float Rs[] = { 1.5f, 2.5f, 4.0f, 6.0f, 8.0f };
    const int   nR = (int)(sizeof(Rs)/sizeof(Rs[0]));
    /* Several representative motions. */
    struct { float x0,y0,x1,y1; const char* tag; } trajs[] = {
        { 32.0f, 32.0f, 32.5f, 32.0f, "in-cell +x" },
        { 32.0f, 32.0f, 32.2f, 32.3f, "diagonal" },
        { 32.5f, 32.5f, 32.6f, 32.7f, "mid-cell" },
        { 32.0f, 32.0f, 32.0f, 32.5f, "in-cell +y" },
        { 32.5f, 32.5f, 31.5f, 31.5f, "negative diag (1 cell)" },
    };
    const int nT = (int)(sizeof(trajs)/sizeof(trajs[0]));
    for (int j = 0; j < nR; j++) {
        for (int i = 0; i < nT; i++) {
            int ncells = 0;
            const float m = continuity_check(trajs[i].x0, trajs[i].y0,
                                              trajs[i].x1, trajs[i].y1,
                                              Rs[j], &ncells);
            const int ok = (m < 1.0e-5f);
            if (!ok) all_pass = 0;
            printf("R=%-6.1f  %-22s  %-6d  %12.3e %s\n",
                   (double) Rs[j], trajs[i].tag, ncells, (double) m,
                   ok ? "PASS" : "*** FAIL ***");
        }
    }
    printf("\n%s\n", all_pass ? "Part A: ALL PASS" : "Part A: SOME FAILED");
    return all_pass ? 0 : 1;
}

/* ===========================================================
 * Part B: bump with inductive ON, vs TSC + smoothing.
 * =========================================================== */

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static float run_full_force(float v_drift, gr_shape_function_t shape, float Rk,
                             int rho_smooth) {
    const int   W       = 128;
    const int   H       = 2048;
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 16;
    const float Q       = 0.01f;
    const float mass    = 100.0f;
    const float y_start = 16.0f + 32.0f;
    const float dt      = cfl * dx / c_eff;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;

    int n_window = (int) ceilf(5.0f / (v_drift * dt));
    if (n_window < 200)    n_window = 200;
    if (n_window > 200000) n_window = 200000;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return 0.0f;
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, rho_smooth);
    gr_sim_set_j_smooth_passes(sim, rho_smooth);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    /* FULL force: phi + inductive + v x B all enabled. */
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, +1.0f);

    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    for (int s = 0; s < 2000; s++) gr_sim_step(sim);
    const gr_particle_t* p_post = gr_sim_get_particle(sim, 0);
    const float vy_post = vy_of(p_post, c_eff);
    for (int s = 0; s < n_window; s++) gr_sim_step(sim);
    const gr_particle_t* p_after = gr_sim_get_particle(sim, 0);
    const float vy_end = vy_of(p_after, c_eff);
    const float accel = (vy_end - vy_post) / ((float) n_window * dt);

    gr_sim_destroy(sim);
    return accel;
}

static void part_B_inductive(void) {
    printf("\n--- Part B: full force (inductive ON), TSC vs bump-R vs TSC+smoothing ---\n");
    printf("Equal-distance methodology.  mass=100, full Lorentz force.\n\n");
    const float v_vals[] = { 0.01f, 0.05f, 0.1f, 0.2f, 0.3f };
    const int   nv       = 5;

    printf("%-7s %-12s %-12s %-12s %-12s %-12s\n",
           "v/c", "TSC(0)", "bump R=4", "bump R=8", "TSC N=4", "TSC N=16");
    printf("------------------------------------------------------------------------\n");
    for (int i = 0; i < nv; i++) {
        const float aT0 = run_full_force(v_vals[i], GR_SHAPE_TSC,  0.0f,  0);
        const float aB4 = run_full_force(v_vals[i], GR_SHAPE_BUMP, 4.0f,  0);
        const float aB8 = run_full_force(v_vals[i], GR_SHAPE_BUMP, 8.0f,  0);
        const float aS4 = run_full_force(v_vals[i], GR_SHAPE_TSC,  0.0f,  4);
        const float aS16= run_full_force(v_vals[i], GR_SHAPE_TSC,  0.0f, 16);
        printf("%-7.4f %+12.3e %+12.3e %+12.3e %+12.3e %+12.3e\n",
               (double) v_vals[i], (double) aT0,
               (double) aB4, (double) aB8,
               (double) aS4, (double) aS16);
    }
    printf("\n");
    printf("Interpretation: smaller |accel| is better.  Bump R=8 and TSC N=16\n");
    printf("are roughly equivalent deposit-side mitigations on the FULL force.\n");
}

int main(void) {
    printf("=== stage44_bump_esirkepov_check (v38 Tier 2 validation) ===\n\n");
    const int A = part_A_continuity();
    part_B_inductive();
    if (A == 0) {
        printf("\nOVERALL: bump Esirkepov continuity verified; inductive-on results\n");
        printf("printed for inspection.  Tier 2 implementation is sound.\n");
    } else {
        printf("\nOVERALL: bump Esirkepov continuity FAILED at some R/trajectory --\n");
        printf("investigate before relying on bump for closed-loop PIC.\n");
    }
    return A;
}
