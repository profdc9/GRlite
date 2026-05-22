/* Stage 24 — Electrostatic Lorentz force (q E piece), unit isolation.
 *
 * Wires the -q grad phi_em piece of the EM Lorentz force into the pusher
 * (Stage 23 covered the q v x B piece).  A spatially uniform electric
 * field E_0 along +x is installed via the new
 * gr_sim_set_background_uniform_electric generator:
 *
 *     phi^{bg}(x, y) = -E_0 (x - x_0)
 *     -grad phi^{bg} = (E_0, 0)              (uniform, by construction)
 *
 * The force on a charged particle reduces to F = q E_0 hat x.  Starting
 * from rest, the relativistic momentum integrates exactly to
 *     p_x(t) = q E_0 t,  p_y(t) = 0,
 * with the Boris half-step convention storing p^{n - 1/2}; at sim step N
 * the recorded momentum is p^{N - 1/2} = q E_0 (N - 1/2) dt.  Position
 * follows from
 *     x(t) - x(0) = (m c^2 / q E_0) [sqrt(1 + (q E_0 t / m c)^2) - 1].
 *
 * Checks at q = 1, m = 1, E_0 = 1e-3, c = 1, N = 200 steps:
 *   (1) p_x matches q E_0 (N - 1/2) dt to part-in-1e5.
 *   (2) p_y = 0 exactly (no transverse force).
 *   (3) x(t) matches the relativistic prediction.
 *   (4) sign reversal: E_0 -> -E_0 reverses the acceleration direction.
 *   (5) neutral particle (q = 0) feels no force (baseline isolation).
 *   (6) ANALYTIC and SAMPLED bg modes give identical results (linear phi
 *       -> centered FD is exact).
 *
 * No magnetic field is installed here, so the v x B piece is zero and
 * the only EM Lorentz contribution under test is -q grad phi. */

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
    float px_end;
    float py_end;
    float x_drift;       /* x(t) - x(0) at end */
    float px_ana;
    float x_ana;
    float px_err_frac;
    float x_err_frac;
    int   nan;
} run_t;

static void run_uniform_e(int sign_E, float charge, gr_bg_mode_t bg_mode,
                          int n_steps, run_t* out) {
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float E0     = (float) sign_E * 1.0e-3f;
    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");

    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_background_uniform_electric(sim, cx, cy, E0, 0.0f);
    gr_sim_set_bg_mode(sim, bg_mode);

    gr_sim_add_particle(sim, cx, cy, /*mass=*/1.0f, charge,
                        /*vx=*/0.0f, /*vy=*/0.0f);
    const float dt = gr_sim_dt(sim);
    gr_sim_step_n(sim, n_steps);

    const gr_particle_t* p = gr_sim_get_particle(sim, 0);
    if (!isfinite(p->x) || !isfinite(p->y)
        || !isfinite(p->px) || !isfinite(p->py)) {
        out->nan = 1;
        gr_sim_destroy(sim);
        return;
    }

    /* Analytic predictions (Boris half-step convention: at step N, the
     * stored momentum is p^{N - 1/2}).  Position follows from the
     * relativistic kinematic integral. */
    const float qE     = charge * E0;
    const float t_half = ((float) n_steps - 0.5f) * dt;
    const float px_ana = qE * t_half;
    const float t_pos  = (float) n_steps * dt;
    /* x(t) - x(0) = (m c^2 / qE) [sqrt(1 + (qE t / m c)^2) - 1].  For
     * neutral particle (qE = 0) this is 0 by limit. */
    float x_ana = 0.0f;
    if (fabsf(qE) > 0.0f) {
        const float u   = qE * t_pos / (1.0f * c_eff);   /* qE t / (m c), m=1 */
        x_ana = (1.0f * c_eff * c_eff / qE) * (sqrtf(1.0f + u * u) - 1.0f);
    }

    out->px_end      = p->px;
    out->py_end      = p->py;
    out->x_drift     = p->x - cx;
    out->px_ana      = px_ana;
    out->x_ana       = x_ana;
    out->px_err_frac = (fabsf(px_ana) > 0.0f)
                          ? fabsf(p->px - px_ana) / fabsf(px_ana) : 0.0f;
    out->x_err_frac  = (fabsf(x_ana) > 0.0f)
                          ? fabsf(out->x_drift - x_ana) / fabsf(x_ana) : 0.0f;
    out->nan         = 0;
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage24_uniform_e_field ===\n");
    printf("EM electrostatic force F = -q grad phi.  Constant E_0 = 1e-3 hat x.\n");
    printf("Spec: gr_sandbox_v35.tex sec:alg_rel Tier-3 eqbox line 1045.\n\n");

    const int N = 200;
    run_t r_ana_pos, r_ana_neg, r_samp_pos, r_neutral;
    run_uniform_e(+1, +1.0f, GR_BG_MODE_ANALYTIC, N, &r_ana_pos);
    run_uniform_e(-1, +1.0f, GR_BG_MODE_ANALYTIC, N, &r_ana_neg);
    run_uniform_e(+1, +1.0f, GR_BG_MODE_SAMPLED,  N, &r_samp_pos);
    run_uniform_e(+1,  0.0f, GR_BG_MODE_ANALYTIC, N, &r_neutral);

    printf("[ANALYTIC bg, E_0 = +1e-3, q = +1]\n");
    printf("  analytic p_x(t = (N-1/2) dt) = %+.6e\n", (double) r_ana_pos.px_ana);
    printf("  measured p_x                = %+.6e   (err %.4e)\n",
           (double) r_ana_pos.px_end, (double) r_ana_pos.px_err_frac);
    printf("  measured p_y                = %+.6e   (expect 0)\n",
           (double) r_ana_pos.py_end);
    printf("  analytic x(N dt) - x(0)     = %+.6e\n", (double) r_ana_pos.x_ana);
    printf("  measured drift              = %+.6e   (err %.4e)\n",
           (double) r_ana_pos.x_drift, (double) r_ana_pos.x_err_frac);

    printf("\n[ANALYTIC bg, E_0 = -1e-3, q = +1]\n");
    printf("  measured p_x                = %+.6e   (expect %+.6e)\n",
           (double) r_ana_neg.px_end, (double) r_ana_neg.px_ana);

    printf("\n[SAMPLED bg, E_0 = +1e-3, q = +1]\n");
    printf("  measured p_x                = %+.6e   (err %.4e)\n",
           (double) r_samp_pos.px_end, (double) r_samp_pos.px_err_frac);

    printf("\n[neutral particle q=0, E_0 = +1e-3]\n");
    printf("  measured p_x                = %+.6e   (expect 0)\n",
           (double) r_neutral.px_end);
    printf("  measured drift              = %+.6e   (expect 0)\n",
           (double) r_neutral.x_drift);

    /* Assertions. */
    TEST_ASSERT(!r_ana_pos.nan && !r_ana_neg.nan && !r_samp_pos.nan && !r_neutral.nan,
                "NaN in one of the runs");
    TEST_ASSERT(r_ana_pos.px_err_frac < 1.0e-5f,
                "ANALYTIC q=+1 E>0: p_x rel.err %.4e exceeds 1e-5",
                (double) r_ana_pos.px_err_frac);
    TEST_ASSERT(fabsf(r_ana_pos.py_end) < 1.0e-10f,
                "ANALYTIC q=+1 E>0: p_y is nonzero (%.4e)",
                (double) r_ana_pos.py_end);
    /* Position is integrated over the v drift; leapfrog truncation enters
     * here at O(dt^2 * a) per step, so we allow a looser bound. */
    TEST_ASSERT(r_ana_pos.x_err_frac < 1.0e-3f,
                "ANALYTIC q=+1 E>0: position rel.err %.4e exceeds 1e-3",
                (double) r_ana_pos.x_err_frac);
    TEST_ASSERT(r_ana_neg.px_end < 0.0f,
                "ANALYTIC q=+1 E<0: p_x expected negative, got %.4e",
                (double) r_ana_neg.px_end);
    TEST_ASSERT(fabsf(r_ana_neg.px_end + r_ana_pos.px_end) < 1.0e-8f,
                "Sign reversal: |p_x_pos + p_x_neg| = %.4e should be ~0",
                (double) fabsf(r_ana_neg.px_end + r_ana_pos.px_end));
    TEST_ASSERT(r_samp_pos.px_err_frac < 1.0e-5f,
                "SAMPLED q=+1 E>0: p_x rel.err %.4e exceeds 1e-5",
                (double) r_samp_pos.px_err_frac);
    TEST_ASSERT(fabsf(r_samp_pos.px_end - r_ana_pos.px_end) < 1.0e-8f,
                "ANALYTIC vs SAMPLED disagreement on linear phi: %.4e vs %.4e",
                (double) r_samp_pos.px_end, (double) r_ana_pos.px_end);
    TEST_ASSERT(fabsf(r_neutral.px_end) < 1.0e-10f,
                "neutral q=0: p_x %.4e should be 0",
                (double) r_neutral.px_end);
    TEST_ASSERT(fabsf(r_neutral.x_drift) < 1.0e-10f,
                "neutral q=0: x drift %.4e should be 0",
                (double) r_neutral.x_drift);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
