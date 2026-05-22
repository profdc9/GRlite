/* Stage 25 — Inductive EM force (-q d_t A piece), unit isolation.
 *
 * Wires the -q d_t A inductive piece of the EM Lorentz force into the
 * pusher.  Combined with Stages 23 (q v x B) and 24 (-q grad phi), the
 * full EM Lorentz force F = q (-grad phi - d_t A + v x B) is now
 * implemented.  Spec: gr_sandbox_v35.tex sec:alg_rel Tier-3 eqbox
 * line 1045.
 *
 * The inductive piece is the trickiest to test in isolation because the
 * "natural" way to drive d_t A is to run the wave equation with a
 * time-varying source, which couples the inductive force to a complex
 * spatial profile.  Stage 25 instead drives d_t A directly by
 * MANUALLY pre-populating the fields[A_X].curr and fields[A_X].next
 * buffers to mimic the post-rotation state of the leapfrog
 *     curr = A^{n+1},  next = A^{n-1}
 * with field_evolution OFF (so the buffers don't rotate and the
 * manually-set values persist).  Specifically:
 *
 *   fields[A_X].curr[k] = +A_0 dt   (i.e., A_x at t = +dt)
 *   fields[A_X].next[k] = -A_0 dt   (i.e., A_x at t = -dt)
 *   centered difference  = (curr - next) / (2 dt) = A_0    (uniform).
 *
 * No scalar EM background, no magnetic field, so the only force on the
 * particle is the inductive piece:
 *
 *     F_x = -q d_t A_x = -q A_0,     F_y = 0.
 *
 * A charged particle at rest then accelerates uniformly along -x at
 * acceleration q A_0 / (gamma m), exactly the kinematics of Stage 24's
 * electrostatic test but driven by the OPPOSITE force-law piece.  At
 * Boris half-step convention, p_x^{N-1/2} = -q A_0 (N - 1/2) dt.
 *
 * This is a direct unit verification of the dt_A_em_at_total kernel and
 * its wiring into em_force_at, independent of any wave-equation
 * dynamics. */

#include "grlite.h"
#include "sim_internal.h"   /* direct fields[].curr / .next access for setup */

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
    float px_end;
    float py_end;
    float px_ana;
    float px_err_frac;
    int   nan;
} run_t;

/* Run with a uniform d_t A_x = A_0 imposed by direct buffer manipulation. */
static void run_inductive(int sign_A, float charge, int n_steps, run_t* out) {
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float A0     = (float) sign_A * 1.0e-3f;     /* d_t A_x target */
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");

    gr_sim_set_field_evolution(sim, 0);    /* no leapfrog, no rotation */
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);   /* no bg fields installed */

    const float dt = gr_sim_dt(sim);
    /* Manually populate fields[A_X].curr and .next so that the centered
     * difference (curr - next) / (2 dt) = A_0 at every interior cell.
     * Set curr = +A_0 dt, next = -A_0 dt uniformly over the grid. */
    float* Ax_curr = sim->fields[GR_FIELD_A_X].curr;
    float* Ax_next = sim->fields[GR_FIELD_A_X].next;
    const float curr_val = +A0 * dt;
    const float next_val = -A0 * dt;
    for (int k = 0; k < W * H; k++) {
        Ax_curr[k] = curr_val;
        Ax_next[k] = next_val;
    }
    /* A_y curr/next stay at zero -> no d_t A_y. */

    gr_sim_add_particle(sim, cx, cy, /*mass=*/1.0f, charge,
                        /*vx=*/0.0f, /*vy=*/0.0f);
    gr_sim_step_n(sim, n_steps);

    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    if (!isfinite(p->px) || !isfinite(p->py)) {
        out->nan = 1;
        gr_sim_destroy(sim);
        return;
    }

    /* Analytic prediction: a uniform d_t A_x = A_0 produces inductive
     * E_x = -A_0, hence F_x = -q A_0 (NOTE: this is opposite-signed from
     * Stage 24's q E_0 with the same E_0 sign, because here E is induced
     * by d_t A, not by a static phi gradient).  Boris half-step
     * convention -> p_x^{N - 1/2} = -q A_0 (N - 1/2) dt. */
    const float qA       = charge * A0;
    const float t_half   = ((float) n_steps - 0.5f) * dt;
    const float px_ana   = -qA * t_half;

    out->px_end      = p->px;
    out->py_end      = p->py;
    out->px_ana      = px_ana;
    out->px_err_frac = (fabsf(px_ana) > 0.0f)
                          ? fabsf(p->px - px_ana) / fabsf(px_ana) : 0.0f;
    out->nan         = 0;
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage25_inductive_e ===\n");
    printf("EM inductive force F_ind = -q d_t A.  Uniform d_t A_x = 1e-3.\n");
    printf("Spec: gr_sandbox_v35.tex sec:alg_rel Tier-3 eqbox line 1045.\n\n");

    const int N = 200;
    run_t r_pos, r_neg, r_neutral;
    run_inductive(+1, +1.0f, N, &r_pos);
    run_inductive(-1, +1.0f, N, &r_neg);
    run_inductive(+1,  0.0f, N, &r_neutral);

    printf("[d_t A_x = +1e-3, q = +1]\n");
    printf("  analytic p_x = %+.6e\n", (double) r_pos.px_ana);
    printf("  measured p_x = %+.6e   (err %.4e)\n",
           (double) r_pos.px_end, (double) r_pos.px_err_frac);
    printf("  measured p_y = %+.6e   (expect 0)\n", (double) r_pos.py_end);

    printf("\n[d_t A_x = -1e-3, q = +1]\n");
    printf("  measured p_x = %+.6e   (expect %+.6e)\n",
           (double) r_neg.px_end, (double) r_neg.px_ana);

    printf("\n[neutral particle q=0]\n");
    printf("  measured p_x = %+.6e   (expect 0)\n", (double) r_neutral.px_end);

    /* Assertions. */
    TEST_ASSERT(!r_pos.nan && !r_neg.nan && !r_neutral.nan, "NaN in run");
    TEST_ASSERT(r_pos.px_err_frac < 1.0e-5f,
                "d_t A > 0: p_x rel.err %.4e exceeds 1e-5",
                (double) r_pos.px_err_frac);
    TEST_ASSERT(fabsf(r_pos.py_end) < 1.0e-10f,
                "d_t A > 0: p_y is nonzero (%.4e)",
                (double) r_pos.py_end);
    TEST_ASSERT(r_neg.px_end > 0.0f,
                "d_t A < 0, q > 0: p_x expected positive (force in +x), got %.4e",
                (double) r_neg.px_end);
    TEST_ASSERT(fabsf(r_neg.px_end + r_pos.px_end) < 1.0e-8f,
                "Sign reversal: |p_x(dtA>0) + p_x(dtA<0)| = %.4e should be ~0",
                (double) fabsf(r_neg.px_end + r_pos.px_end));
    TEST_ASSERT(fabsf(r_neutral.px_end) < 1.0e-10f,
                "neutral q=0: p_x %.4e should be 0",
                (double) r_neutral.px_end);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
