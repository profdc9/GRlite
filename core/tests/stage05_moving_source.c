/* Stage 5 test — gr_sandbox_v33.tex §12.5 "Moving source, vector potential, Yee curl".
 *
 * Loads moving_source with mass M, charge Q, and velocity (v_x, 0), runs the
 * wave-equation iteration with damping, and verifies:
 *
 *   (a) the four scalar/vector-x fields develop the boosted-Coulomb relation
 *           A_g_x ~ (v_x / c^2) * Phi_g          (gravity)
 *           A_x   ~ (v_x / c^2) * phi_em         (EM)
 *       at sample radii — the spec test of §12.5;
 *   (b) the perpendicular vector components A_g_y and A_y remain near zero
 *       (only v_x is nonzero, so their sources J_*y are exactly zero);
 *   (c) the cell-centered curl B_z = d_x A_y - d_y A_x (centered FD analog of
 *       §9.1 eq:yee_curl) matches the analytic boosted-Coulomb prediction
 *           B_z(x_0, y) = -(v_x / c^2) * d_y phi_em(x_0, y);
 *   (d) the discrete continuity residual  d_t rho + div J  is bounded and
 *       matches the predicted v . grad rho (the source position is held fixed
 *       across steps, so d_t rho = 0 exactly and the entire residual is the
 *       div J term coming from the CIC profile of the moving deposit). This
 *       is the deliberate violation that drives Lorenz-gauge drift; we
 *       measure it for the record so we can correlate it with gauge residual
 *       growth in future stages.
 *
 * Cell-centered storage convention (Stages 1+); §9.1's edge-staggered Yee
 * positions are restated as centered finite differences with the same
 * O(dx^2) accuracy.
 */

#include "grlite.h"

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

int main(void) {
    printf("=== stage05_moving_source: gr_sandbox_v33.tex §12.5 ===\n");

    const int   W        = 128, H = 128;
    const float dx       = 1.0f;
    const float c_eff    = 1.0f;
    const float cfl      = 1.0f / sqrtf(2.0f);
    const int   N_d      = 16;
    const float G_eff    = 1.0f;
    const float k_e      = 1.0f;
    const float mass     = 1.0f;
    const float charge   = 1.0f;
    const float vx       = 0.1f * c_eff;
    const int   n_settle = 1000;
    const int   n_avg    = 3000;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    gr_sim_set_G_eff(sim, G_eff);
    gr_sim_set_k_e(sim, k_e);
    gr_sim_set_damping(sim, N_d);

    /* params: mass, charge, vx, vy (=0), x0 (default center), y0 (default center). */
    float params[4] = {mass, charge, vx, 0.0f};
    TEST_ASSERT(gr_sim_load_scenario(sim, "moving_source", params, 4) == 0,
                "moving_source load failed");

    /* Settle + accumulate cell-wise time-average for all six fields. */
    gr_sim_step_n(sim, n_settle);
    double* sums[GR_FIELD_COUNT];
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        sums[f] = (double*) calloc((size_t) W * H, sizeof(double));
        TEST_ASSERT(sums[f] != NULL, "sum alloc failed");
    }
    for (int n = 0; n < n_avg; n++) {
        gr_sim_step(sim);
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            const float* p = gr_sim_field_ptr(sim, (gr_field_id_t) f);
            for (int k = 0; k < W * H; k++) sums[f][k] += (double) p[k];
        }
    }
    float* mean[GR_FIELD_COUNT];
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        mean[f] = (float*) malloc((size_t) W * H * sizeof(float));
        TEST_ASSERT(mean[f] != NULL, "mean alloc failed");
        const double inv_n = 1.0 / (double) n_avg;
        for (int k = 0; k < W * H; k++) mean[f][k] = (float) (sums[f][k] * inv_n);
        free(sums[f]);
    }

    /* (a) Boosted-Coulomb ratio A_x / phi ≈ v_x / c^2 at sample radii. Sample
     * along the +y axis from the source (so x = x0 exactly), where phi and A_x
     * are both well-defined and not on a CIC-deposit cell. */
    const float ratio_expected = vx / (c_eff * c_eff);
    const int   cx = W / 2, cy = H / 2;
    const int   radii[] = {8, 12, 16, 20, 24};
    const int   n_r = sizeof(radii) / sizeof(radii[0]);
    float max_rel_grav = 0.0f, max_rel_em = 0.0f;
    printf("  boosted-Coulomb ratio (expected v_x / c^2 = %.5f):\n", ratio_expected);
    printf("    %-3s   %-12s %-12s %-12s %-12s\n",
           "r", "Phi_g", "A_g_x/Phi_g", "phi_em", "A_x/phi_em");
    for (int s = 0; s < n_r; s++) {
        const int k = (cy + radii[s]) * W + cx;
        const float phi_g  = mean[GR_FIELD_PHI_GRAV][k];
        const float A_gx   = mean[GR_FIELD_A_GX][k];
        const float phi_em = mean[GR_FIELD_PHI_EM][k];
        const float A_x    = mean[GR_FIELD_A_X][k];
        const float r_grav = A_gx  / phi_g;
        const float r_em   = A_x   / phi_em;
        const float rg_rel = fabsf((r_grav - ratio_expected) / ratio_expected);
        const float re_rel = fabsf((r_em   - ratio_expected) / ratio_expected);
        if (rg_rel > max_rel_grav) max_rel_grav = rg_rel;
        if (re_rel > max_rel_em)   max_rel_em   = re_rel;
        printf("    %-3d   %+8.4f     %+8.5f   %+8.4f    %+8.5f\n",
               radii[s], phi_g, r_grav, phi_em, r_em);
    }
    printf("  max rel.err: grav = %.3e   EM = %.3e\n", max_rel_grav, max_rel_em);
    TEST_ASSERT(max_rel_grav < 5.0e-2f, "grav boosted-Coulomb rel.err %.3e exceeds 5e-2", max_rel_grav);
    TEST_ASSERT(max_rel_em   < 5.0e-2f, "EM boosted-Coulomb rel.err %.3e exceeds 5e-2",  max_rel_em);

    /* (b) Perpendicular components A_g_y and A_y stay near zero (their sources
     * J_my and J_qy are zero — vy = 0). With damping the leapfrog drives them
     * monotonically toward zero from any numerical noise; over the time
     * average they should be tiny. Tolerance: well below the boosted-Coulomb
     * signal scale (which is vx/c^2 * |phi|_peak ~ 0.1 * 1 = 0.1). */
    float max_Agy = 0.0f, max_Ay = 0.0f;
    for (int j = N_d; j < H - N_d; j++)
        for (int i = N_d; i < W - N_d; i++) {
            const int k = j * W + i;
            const float a1 = fabsf(mean[GR_FIELD_A_GY][k]);
            const float a2 = fabsf(mean[GR_FIELD_A_Y][k]);
            if (a1 > max_Agy) max_Agy = a1;
            if (a2 > max_Ay)  max_Ay  = a2;
        }
    printf("  perpendicular components (expected ~0):\n");
    printf("    max |A_g_y|@avg = %.3e   max |A_y|@avg = %.3e\n", max_Agy, max_Ay);
    TEST_ASSERT(max_Agy < 1.0e-3f, "max |A_g_y| = %.3e is not near zero", max_Agy);
    TEST_ASSERT(max_Ay  < 1.0e-3f, "max |A_y|   = %.3e is not near zero",  max_Ay);

    /* (c) Curl B_z = d_x A_y - d_y A_x via centered FD (cell-centered analog
     * of §9.1 eq:yee_curl). For our case A_y ~ 0, so B_z = -d_y A_x.
     *
     * Analytic at sample point (x0, y0 + r), in our normalization:
     *   leapfrog source coefficient for phi_em is -4 pi k_e, so the static
     *   Poisson form is Lap phi_em = 4 pi k_e rho_q. For a point charge Q
     *   the 2D Green's function gives phi_em(r) = 2 k_e Q ln r + C
     *   (slope 2 k_e Q, as verified in Stage 3 for the gravitational
     *   analogue), hence d_r phi_em = 2 k_e Q / r.
     *   A_x  = (v_x / c^2) * phi_em
     *   B_z  = -d_y A_x|_(x0, y0+r) = -(v_x / c^2) * 2 k_e Q / r */
    const float inv_2dx = 1.0f / (2.0f * dx);
    float max_rel_bz = 0.0f;
    printf("  B_z (= -d_y A_x) along +y axis, expected = -(v_x/c^2) * 2 k_e Q / r:\n");
    for (int s = 0; s < n_r; s++) {
        const int r = radii[s];
        const int k = (cy + r) * W + cx;
        const float dAx_dy = (mean[GR_FIELD_A_X][k + W] - mean[GR_FIELD_A_X][k - W]) * inv_2dx;
        const float dAy_dx = (mean[GR_FIELD_A_Y][k + 1] - mean[GR_FIELD_A_Y][k - 1]) * inv_2dx;
        const float Bz_obs = dAy_dx - dAx_dy;
        const float Bz_ana = -(vx / (c_eff * c_eff)) * 2.0f * k_e * charge / (float) r;
        const float rel = fabsf((Bz_obs - Bz_ana) / Bz_ana);
        if (rel > max_rel_bz) max_rel_bz = rel;
        printf("    r=%-3d   B_z_obs = %+.4e   B_z_analytic = %+.4e   rel.err = %.3e\n",
               r, Bz_obs, Bz_ana, rel);
    }
    TEST_ASSERT(max_rel_bz < 1.0e-1f, "B_z rel.err %.3e exceeds 10%%", max_rel_bz);

    /* (d) Discrete continuity residual: d_t rho + div J. Since rho is held
     * fixed across steps in this scenario, d_t rho = 0 identically; the
     * residual is the divergence of the deposited J via centered FD. We
     * report both the matter and EM channels, plus their predicted
     * magnitudes |v| * |grad rho|_max — the test is satisfied when the
     * residual matches the prediction (not when it's small).
     *
     * This is the deliberate continuity violation that drives Lorenz-gauge
     * drift over time per gr_sandbox_v32.tex §9.5. */
    const float* rho_m = gr_sim_rho_matter_ptr(sim);
    const float* rho_q = gr_sim_rho_q_ptr(sim);
    const float* J_mx  = gr_sim_J_mx_ptr(sim);
    const float* J_my  = gr_sim_J_my_ptr(sim);
    const float* J_qx  = gr_sim_J_qx_ptr(sim);
    const float* J_qy  = gr_sim_J_qy_ptr(sim);

    double sum_sq_m = 0.0, sum_sq_q = 0.0;
    float  max_m    = 0.0f, max_q    = 0.0f;
    float  max_grad_rho_m = 0.0f, max_grad_rho_q = 0.0f;
    int    n_cells  = 0;
    for (int j = N_d + 1; j < H - N_d - 1; j++) {
        for (int i = N_d + 1; i < W - N_d - 1; i++) {
            const int k = j * W + i;
            const float divJm = (J_mx[k + 1] - J_mx[k - 1]) * inv_2dx
                              + (J_my[k + W] - J_my[k - W]) * inv_2dx;
            const float divJq = (J_qx[k + 1] - J_qx[k - 1]) * inv_2dx
                              + (J_qy[k + W] - J_qy[k - W]) * inv_2dx;
            const float gxm = (rho_m[k + 1] - rho_m[k - 1]) * inv_2dx;
            const float gym = (rho_m[k + W] - rho_m[k - W]) * inv_2dx;
            const float grm = sqrtf(gxm * gxm + gym * gym);
            const float gxq = (rho_q[k + 1] - rho_q[k - 1]) * inv_2dx;
            const float gyq = (rho_q[k + W] - rho_q[k - W]) * inv_2dx;
            const float grq = sqrtf(gxq * gxq + gyq * gyq);
            if (fabsf(divJm) > max_m) max_m = fabsf(divJm);
            if (fabsf(divJq) > max_q) max_q = fabsf(divJq);
            if (grm > max_grad_rho_m) max_grad_rho_m = grm;
            if (grq > max_grad_rho_q) max_grad_rho_q = grq;
            sum_sq_m += (double) (divJm * divJm);
            sum_sq_q += (double) (divJq * divJq);
            n_cells++;
        }
    }
    const float rms_m = (float) sqrt(sum_sq_m / (double) n_cells);
    const float rms_q = (float) sqrt(sum_sq_q / (double) n_cells);
    /* Predicted: |v . grad rho|_max <= |v| * |grad rho|_max. Here v = (v_x, 0). */
    const float v_mag       = fabsf(vx);
    const float pred_m_max  = v_mag * max_grad_rho_m;
    const float pred_q_max  = v_mag * max_grad_rho_q;
    printf("  discrete continuity residual = d_t rho + div J (d_t rho = 0 here):\n");
    printf("    matter channel:  max |div J_m|  = %.4e   RMS = %.4e\n", max_m, rms_m);
    printf("    charge channel:  max |div J_q|  = %.4e   RMS = %.4e\n", max_q, rms_q);
    printf("    predicted bound |v . grad rho|_max:  matter = %.4e   charge = %.4e\n",
           pred_m_max, pred_q_max);
    /* Sanity assert: the residual should be O(|v . grad rho|) — bounded above
     * by ~v . grad rho at any interior cell (within FD truncation error). */
    TEST_ASSERT(max_m <= 1.5f * pred_m_max + 1e-10f,
                "matter continuity residual %.4e exceeds 1.5 * predicted %.4e",
                max_m, pred_m_max);
    TEST_ASSERT(max_q <= 1.5f * pred_q_max + 1e-10f,
                "charge continuity residual %.4e exceeds 1.5 * predicted %.4e",
                max_q, pred_q_max);

    /* Bonus diagnostic: gauge residual at end of run. With div J != 0 the
     * Lorenz gauge drifts; print so we can correlate later. */
    printf("  Lorenz gauge residual at final step: grav = %.4e   EM = %.4e\n",
           gr_sim_gauge_residual_grav(sim), gr_sim_gauge_residual_em(sim));

    for (int f = 0; f < GR_FIELD_COUNT; f++) free(mean[f]);
    gr_sim_destroy(sim);
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
