/* Stage 34 -- end-to-end sign verification of the EM chain.
 *
 * Verifies each piece of the EM machinery independently, with the
 * minimum setup needed for each test.  This is the post-v36 audit
 * the user requested after the sign-convention flip: confirm at every
 * step that the simulator implements standard Maxwell.
 *
 *   34a  Static Poisson sign:
 *        Lap phi - (1/c^2) d^2 phi/dt^2 = -rho_q / epsilon_0
 *        Static limit:  Lap phi = -4 pi k_e rho_q.
 *        => For +Q point deposit, phi_em peaks POSITIVELY at the source
 *           and decreases with distance.  Discrete Lap(phi) at source
 *           has sign OPPOSITE to rho.
 *
 *   34b  J = rho * v deposition:
 *        Particle of charge Q moving at (vx, vy) must deposit total
 *        J_qx = Q*vx, J_qy = Q*vy (signs and magnitudes both match).
 *
 *   34c  A wave-equation sign:
 *        Lap A - (1/c^2) d^2 A/dt^2 = -J / (epsilon_0 c^2)
 *        Static limit:  Lap A = -4 pi k_e J / c^2.
 *        => Manually set J_qy at a cell and watch A_y develop.  Peak
 *           POSITIVE at source; Lap(A_y) sign OPPOSITE to J_y.
 *
 *   34d  Lorentz force direction:
 *        F_E = -q grad phi   (force on +q points DOWN the phi gradient).
 *        F_B = q v x B       (right-hand rule: +q with v=+x, B=+z gives
 *                             F in -y direction).
 *        Verified via momentum kick from analytic uniform-E and
 *        uniform-B backgrounds. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); \
        return 1; \
    } \
} while (0)

/* ============================================================
 * 34a: phi_em Poisson sign
 * ============================================================ */
static int test_phi_poisson_sign(void) {
    printf("\n[34a] Poisson sign: Lap(phi_em) = -4 pi k_e rho_q.\n");
    const int   W      = 64, H = 64;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create");
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 0);   /* manual sources */
    gr_sim_set_damping(sim, 8);

    /* Manually inject a +Q point charge at the center cell. */
    const int   ci = W / 2;
    const int   cj = H / 2;
    const float Q  = 1.0f;
    float* rho_q = (float*) gr_sim_array_ptr(sim, GR_ARR_RHO_Q);
    memset(rho_q, 0, (size_t) W * H * sizeof(float));
    rho_q[cj * W + ci] = Q;

    /* Run the leapfrog ~80 steps to let phi_em approach quasi-static.
     * Source array stays nonzero across steps because particle_source_
     * deposition is OFF (no gr_sim_clear_sources call). */
    for (int s = 0; s < 80; s++) gr_sim_step(sim);

    const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_EM);
    const float phi_at_source = phi[cj * W + ci];
    const float phi_at_d8     = phi[cj * W + (ci + 8)];
    const float phi_at_d16    = phi[cj * W + (ci + 16)];
    printf("  phi_em(at +Q source) = %+.5e\n", phi_at_source);
    printf("  phi_em(at +8 cells)  = %+.5e\n", phi_at_d8);
    printf("  phi_em(at +16 cells) = %+.5e\n", phi_at_d16);

    /* For +Q: phi has POSITIVE local peak at source, monotonically
     * decreasing in the damped box. */
    TEST_ASSERT(phi_at_source > 0.0f,
                "phi_em at +Q source should be POSITIVE, got %+.3e", phi_at_source);
    TEST_ASSERT(phi_at_d8 > 0.0f && phi_at_d8 < phi_at_source,
                "phi_em should decrease away from +Q: source=%+.3e, +8=%+.3e",
                phi_at_source, phi_at_d8);
    TEST_ASSERT(phi_at_d16 > 0.0f && phi_at_d16 < phi_at_d8,
                "phi_em should be monotonic: +8=%+.3e, +16=%+.3e",
                phi_at_d8, phi_at_d16);

    /* Discrete Laplacian at the source cell.  In the damped box with
     * a delta source, the dynamic field hasn't perfectly converged in
     * 80 steps, but the SIGN of Lap(phi) at the source cell should
     * be NEGATIVE (Lap < 0 at the deposit) -- this is the standard
     * Maxwell Poisson sign.
     *
     * Discrete Laplacian uses 5-point stencil:
     *   Lap(phi)[k] = (phi[k-1] + phi[k+1] + phi[k-W] + phi[k+W] - 4 phi[k]) / dx^2 */
    const int k = cj * W + ci;
    const float lap_phi = (phi[k - 1] + phi[k + 1] + phi[k - W] + phi[k + W]
                          - 4.0f * phi[k]) / (dx * dx);
    printf("  Lap(phi_em) at source cell = %+.5e (expect < 0 for +Q)\n", lap_phi);
    TEST_ASSERT(lap_phi < 0.0f,
                "Lap(phi_em) at +Q source should be NEGATIVE (standard Poisson "
                "Lap phi = -4 pi k_e rho), got %+.3e", lap_phi);

    gr_sim_destroy(sim);
    return 0;
}

/* ============================================================
 * 34b: J = rho * v deposition
 * ============================================================ */
static int test_J_equals_rho_v(void) {
    printf("\n[34b] Deposition: J_q = rho_q * v.\n");
    const int   W      = 64, H = 64;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create");
    gr_sim_set_field_evolution(sim, 0);   /* no wave evolution */
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_shape_function(sim, GR_SHAPE_CIC);

    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;
    const float Q      = 0.5f;
    const float vx     = +0.1f;
    const float vy     = -0.2f;
    /* Heavy mass so the particle barely moves during the single step. */
    gr_sim_add_particle(sim, cx, cy, 1.0e6f, Q, vx, vy);

    /* One step: clears sources, then deposits.  rho_q and J_qx/J_qy
     * arrays now hold the deposit for this particle. */
    gr_sim_step(sim);

    /* Sum total deposited charge and current. */
    const float* rho_q = gr_sim_rho_q_ptr(sim);
    const float* J_qx  = gr_sim_J_qx_ptr(sim);
    const float* J_qy  = gr_sim_J_qy_ptr(sim);
    double Q_sum = 0.0, Jx_sum = 0.0, Jy_sum = 0.0;
    for (int k = 0; k < W * H; k++) {
        Q_sum  += (double) rho_q[k];
        Jx_sum += (double) J_qx[k];
        Jy_sum += (double) J_qy[k];
    }
    printf("  Sum(rho_q)  = %+.5e (expect %+.3e for Q=%.2f)\n",
           Q_sum, Q, (double) Q);
    printf("  Sum(J_qx)   = %+.5e (expect %+.3e for Q*vx=%.3e)\n",
           Jx_sum, (double) (Q * vx), (double) (Q * vx));
    printf("  Sum(J_qy)   = %+.5e (expect %+.3e for Q*vy=%.3e)\n",
           Jy_sum, (double) (Q * vy), (double) (Q * vy));

    /* Total charge conservation: deposit should sum to Q. */
    TEST_ASSERT(fabs(Q_sum - (double) Q) < 1e-5,
                "Sum(rho_q) = %.4e does not match deposited Q = %.4e",
                Q_sum, (double) Q);
    /* J_qx and J_qy should match Q*v_x and Q*v_y in sign and magnitude. */
    TEST_ASSERT(Jx_sum * (double) (Q * vx) > 0.0,
                "Sign mismatch: Sum(J_qx)=%.4e but Q*vx=%.4e",
                Jx_sum, (double) (Q * vx));
    TEST_ASSERT(Jy_sum * (double) (Q * vy) > 0.0,
                "Sign mismatch: Sum(J_qy)=%.4e but Q*vy=%.4e",
                Jy_sum, (double) (Q * vy));
    /* Tolerance on Esirkepov-deposited J magnitude.  The Esirkepov
     * decomposition guarantees the discrete continuity equation, not
     * total-J = Q*v exactly, but for a single-cell-crossing trajectory
     * the totals match closely. */
    TEST_ASSERT(fabs(Jx_sum - (double) (Q * vx)) < 1e-3 * fabs((double) (Q * vx)),
                "Sum(J_qx)=%.4e off from Q*vx=%.4e", Jx_sum, (double) (Q * vx));
    TEST_ASSERT(fabs(Jy_sum - (double) (Q * vy)) < 1e-3 * fabs((double) (Q * vy)),
                "Sum(J_qy)=%.4e off from Q*vy=%.4e", Jy_sum, (double) (Q * vy));

    gr_sim_destroy(sim);
    return 0;
}

/* ============================================================
 * 34c: A_y wave-equation sign
 * ============================================================ */
static int test_A_wave_eq_sign(void) {
    printf("\n[34c] A_y wave eq: Lap(A_y) = -4 pi k_e J_qy / c^2.\n");
    const int   W      = 64, H = 64;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create");
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 0);   /* manual sources */
    gr_sim_set_damping(sim, 8);

    /* Manually inject a positive J_y at a single Y_EDGE cell. */
    const int   ci = W / 2;
    const int   cj = H / 2;
    const float Jy = 1.0f;
    float* J_qy = (float*) gr_sim_array_ptr(sim, GR_ARR_J_QY);
    memset(J_qy, 0, (size_t) W * H * sizeof(float));
    J_qy[cj * W + ci] = Jy;

    /* Settle the A_y field.  Source persists since deposition is OFF. */
    for (int s = 0; s < 80; s++) gr_sim_step(sim);

    const float* A_y = gr_sim_field_ptr(sim, GR_FIELD_A_Y);
    const float A_at_source = A_y[cj * W + ci];
    const float A_at_d8     = A_y[cj * W + (ci + 8)];
    const float A_at_d16    = A_y[cj * W + (ci + 16)];
    printf("  A_y(at +J source) = %+.5e\n", A_at_source);
    printf("  A_y(at +8 cells)  = %+.5e\n", A_at_d8);
    printf("  A_y(at +16 cells) = %+.5e\n", A_at_d16);

    /* For positive J_y: A_y peaks POSITIVELY at the source (in the
     * damped box) -- same structural relationship as phi to rho. */
    TEST_ASSERT(A_at_source > 0.0f,
                "A_y at +J source should be POSITIVE, got %+.3e", A_at_source);
    TEST_ASSERT(A_at_d8 > 0.0f && A_at_d8 < A_at_source,
                "A_y should decrease away from +J source: %+.3e -> %+.3e",
                A_at_source, A_at_d8);

    /* Discrete Laplacian of A_y at the source cell.  Should be NEGATIVE
     * for positive J_y (Lap A = -4 pi k_e J/c^2 standard Maxwell). */
    const int k = cj * W + ci;
    const float lap_A = (A_y[k - 1] + A_y[k + 1] + A_y[k - W] + A_y[k + W]
                        - 4.0f * A_y[k]) / (dx * dx);
    printf("  Lap(A_y) at source = %+.5e (expect < 0 for +J_y)\n", lap_A);
    TEST_ASSERT(lap_A < 0.0f,
                "Lap(A_y) at +J source should be NEGATIVE (standard Maxwell), got %+.3e",
                lap_A);

    gr_sim_destroy(sim);
    return 0;
}

/* ============================================================
 * 34d: Lorentz force direction
 * ============================================================ */
static int test_lorentz_force_directions(void) {
    printf("\n[34d] Lorentz force: F_E = -q grad phi, F_B = q v x B.\n");

    /* Sub-test (i): F_E direction with uniform E background.
     * Uniform E in +x direction:  phi_bg = -E_0*(x-x_0),  -grad phi = +E_0 x_hat.
     * Positive charge at rest experiences F = +q*E_0 in +x.  Velocity
     * gain after one step is in +x direction. */
    {
        const int   W = 64, H = 64;
        const float dx = 1.0f, c = 1.0f, cfl = 1.0f / sqrtf(2.0f);
        gr_sim_t* sim = gr_sim_create(W, H, dx, c, cfl);
        TEST_ASSERT(sim != NULL, "create");
        gr_sim_set_field_evolution(sim, 0);
        gr_sim_set_particle_source_deposition(sim, 0);
        gr_sim_set_damping(sim, 0);
        gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

        const float cx = ((float) W * 0.5f) * dx;
        const float cy = ((float) H * 0.5f) * dx;
        const float Ex = +1.0e-3f;
        gr_sim_set_background_uniform_electric(sim, cx, cy, Ex, 0.0f);
        gr_sim_add_particle(sim, cx, cy, /*mass=*/1.0f, /*charge=*/+1.0f,
                            0.0f, 0.0f);

        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        printf("  [i] +q in E=+x: p_x = %+.5e (expect > 0), p_y = %+.5e (expect ~0)\n",
               (double) p->px, (double) p->py);
        TEST_ASSERT(p->px > 0.0f,
                    "F_E direction wrong: +q in +E_x should give +p_x, got %+.3e", p->px);
        TEST_ASSERT(fabsf(p->py) < 1e-6f,
                    "+q in pure +E_x should not gain p_y, got %+.3e", p->py);
        gr_sim_destroy(sim);
    }

    /* Sub-test (ii): F_B direction with uniform B_z background.
     * Standard right-hand rule:  +q with v = (v_0, 0, 0) and B = (0, 0, +B_0)
     * gives v x B = (0, -v_0 B_0, 0), so F = q*(v x B) has F_y < 0.
     * After one step, p_y < 0 for the test particle. */
    {
        const int   W = 64, H = 64;
        const float dx = 1.0f, c = 1.0f, cfl = 1.0f / sqrtf(2.0f);
        gr_sim_t* sim = gr_sim_create(W, H, dx, c, cfl);
        TEST_ASSERT(sim != NULL, "create");
        gr_sim_set_field_evolution(sim, 0);
        gr_sim_set_particle_source_deposition(sim, 0);
        gr_sim_set_damping(sim, 0);
        gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

        const float cx = ((float) W * 0.5f) * dx;
        const float cy = ((float) H * 0.5f) * dx;
        const float B0 = +1.0e-3f;
        gr_sim_set_background_uniform_magnetic(sim, cx, cy, B0);
        gr_sim_add_particle(sim, cx, cy, /*mass=*/1.0f, /*charge=*/+1.0f,
                            /*vx=*/+0.1f, /*vy=*/0.0f);

        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        printf("  [ii] +q v=+x in B=+z: p_y = %+.5e (expect < 0)\n",
               (double) p->py);
        TEST_ASSERT(p->py < 0.0f,
                    "F_B = q(v x B) direction wrong: +q with v=+x in B=+z should "
                    "give -p_y deflection, got p_y = %+.3e", p->py);
        gr_sim_destroy(sim);
    }

    /* Sub-test (iii): opposite-charge attractive Coulomb in analytic
     * background.  +Q central, -q test charge at +x: should accelerate
     * INWARD (-x direction). */
    {
        const int   W = 64, H = 64;
        const float dx = 1.0f, c = 1.0f, cfl = 1.0f / sqrtf(2.0f);
        gr_sim_t* sim = gr_sim_create(W, H, dx, c, cfl);
        TEST_ASSERT(sim != NULL, "create");
        gr_sim_set_field_evolution(sim, 0);
        gr_sim_set_particle_source_deposition(sim, 0);
        gr_sim_set_damping(sim, 0);
        gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

        const float cx = ((float) W * 0.5f) * dx;
        const float cy = ((float) H * 0.5f) * dx;
        const float Q  = +1.0f;
        const float eps = 1.0f;
        gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
        gr_sim_add_particle(sim, cx + 10.0f, cy, /*mass=*/1.0f,
                            /*charge=*/-1.0f, 0.0f, 0.0f);

        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        printf("  [iii] -q at +x of +Q: p_x = %+.5e (expect < 0, attractive)\n",
               (double) p->px);
        TEST_ASSERT(p->px < 0.0f,
                    "Opposite charges should attract: -q at +x of +Q should get "
                    "p_x < 0, got %+.3e", p->px);
        gr_sim_destroy(sim);
    }

    /* Sub-test (iv): like-charge repulsive Coulomb in analytic background.
     * +Q central, +q test at +x: should accelerate OUTWARD (+x). */
    {
        const int   W = 64, H = 64;
        const float dx = 1.0f, c = 1.0f, cfl = 1.0f / sqrtf(2.0f);
        gr_sim_t* sim = gr_sim_create(W, H, dx, c, cfl);
        TEST_ASSERT(sim != NULL, "create");
        gr_sim_set_field_evolution(sim, 0);
        gr_sim_set_particle_source_deposition(sim, 0);
        gr_sim_set_damping(sim, 0);
        gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

        const float cx = ((float) W * 0.5f) * dx;
        const float cy = ((float) H * 0.5f) * dx;
        const float Q  = +1.0f;
        const float eps = 1.0f;
        gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
        gr_sim_add_particle(sim, cx + 10.0f, cy, /*mass=*/1.0f,
                            /*charge=*/+1.0f, 0.0f, 0.0f);

        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        printf("  [iv] +q at +x of +Q: p_x = %+.5e (expect > 0, repulsive)\n",
               (double) p->px);
        TEST_ASSERT(p->px > 0.0f,
                    "Like charges should repel: +q at +x of +Q should get "
                    "p_x > 0, got %+.3e", p->px);
        gr_sim_destroy(sim);
    }

    return 0;
}

int main(void) {
    printf("=== stage34_em_chain_signs ===\n");
    printf("End-to-end sign verification: Poisson, J=rho*v, A wave eq, Lorentz force.\n");
    if (test_phi_poisson_sign()        != 0) return 1;
    if (test_J_equals_rho_v()          != 0) return 1;
    if (test_A_wave_eq_sign()          != 0) return 1;
    if (test_lorentz_force_directions()!= 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
