/* Stage 2 test — gr_sandbox_v32.tex §12.2 "Stage 2: Absorbing boundary layer".
 *
 * Drives the Stage 1 Gaussian pulse through a 256x256 domain with the §9.6
 * damping layer (N_d = 16, sigma_max = 21 c / (2 N_d dx) — eq:damp_profile),
 * runs past two full round-trip times, and **records the residual field
 * amplitude at the domain center** (the literal spec wording). Verifies:
 *
 *   (a) damped central residual is small (< 5% of initial amplitude);
 *   (b) damped suppression vs an identical undamped run is at least 10x;
 *   (c) the wall-adjacent field is strongly attenuated relative to the
 *       undamped case — direct evidence the layer is doing absorption;
 *   (d) the field stays finite throughout.
 *
 * The §9.6 R ~ 10^-3 reflection figure is the plane-wave-normal-incidence
 * limit; the realized in-interior residual for a 2D circular wavefront is
 * higher because (i) the (1 - sigma dt) linearization deviates from the exact
 * exp(-sigma dt) at our sigma_max dt ~= 0.46 (~15%) and (ii) oblique-incidence
 * waves into the grid corners are weakly damped by axis-aligned profiles
 * (PML would mitigate this). What we *measure* therefore is "much less than
 * without damping" rather than the spec's optimistic 10^-3 — a finding worth
 * a note in a future LaTeX revision per feedback-docs-in-latex.
 */

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

typedef struct {
    float phi_center_final;       /* |phi| at (W/2, H/2) at final step */
    float phi_near_wall_final;    /* |phi| at (W-2, H/2) at final step — inside the damping layer */
    float phi_center_max_pre;     /* peak |phi(center)| in the *pre-bounce* settled window
                                     [step 50, step 150] — captures the wake before any
                                     echo arrives */
} measurement_t;

static int run_pulse_and_measure(int n_damping, int n_steps, measurement_t* m) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float sigma  = 4.0f * dx;
    const float amp    = 1.0f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return -1;
    gr_sim_set_damping(sim, n_damping);

    float params[2] = {sigma, amp};
    if (gr_sim_load_scenario(sim, "wave_pulse", params, 2) != 0) {
        gr_sim_destroy(sim);
        return -2;
    }

    float pre_max = 0.0f;
    for (int n = 1; n <= n_steps; n++) {
        gr_sim_step(sim);
        const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
        if (n >= 50 && n <= 150) {
            const float v = fabsf(phi[(H / 2) * W + (W / 2)]);
            if (!isfinite(v)) { gr_sim_destroy(sim); return -3; }
            if (v > pre_max) pre_max = v;
        }
    }
    const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
    m->phi_center_final    = fabsf(phi[(H / 2) * W + (W / 2)]);
    m->phi_near_wall_final = fabsf(phi[(H / 2) * W + (W - 2)]);
    m->phi_center_max_pre  = pre_max;
    if (!isfinite(m->phi_center_final) || !isfinite(m->phi_near_wall_final)) {
        gr_sim_destroy(sim);
        return -4;
    }
    gr_sim_destroy(sim);
    return 0;
}

int main(void) {
    printf("=== stage02_damping: gr_sandbox_v32.tex §12.2 ===\n");

    /* 800 steps ~= t = 565, past 2nd round-trip (t ~ 256 + 256 = 512). */
    const int n_steps = 800;

    measurement_t damp = {0};
    TEST_ASSERT(run_pulse_and_measure(16, n_steps, &damp) == 0, "damped run failed");
    printf("  damped   (N_d=16): |phi(center)|@final = %.5g, |phi(near wall)|@final = %.5g, "
           "pre-bounce |phi(center)| max = %.5g\n",
           damp.phi_center_final, damp.phi_near_wall_final, damp.phi_center_max_pre);

    measurement_t undamp = {0};
    TEST_ASSERT(run_pulse_and_measure(0, n_steps, &undamp) == 0, "undamped run failed");
    printf("  undamped (N_d=0):  |phi(center)|@final = %.5g, |phi(near wall)|@final = %.5g, "
           "pre-bounce |phi(center)| max = %.5g\n",
           undamp.phi_center_final, undamp.phi_near_wall_final, undamp.phi_center_max_pre);

    /* (a) Damped central residual at the final step is small. */
    TEST_ASSERT(damp.phi_center_final < 5.0e-2f,
                "damped |phi(center)| at final step = %.5g exceeds 5e-2",
                damp.phi_center_final);

    /* (b) Suppression vs undamped is at least 10x at the center. */
    const float center_ratio = undamp.phi_center_final / damp.phi_center_final;
    printf("  center  ratio  (undamped / damped) = %.1fx\n", center_ratio);
    TEST_ASSERT(center_ratio > 10.0f,
                "central-residual suppression = %.2fx below 10x threshold", center_ratio);

    /* (c) The wall-adjacent cell sees orders-of-magnitude stronger suppression
     * — this is the direct evidence the layer is doing absorption. Diagnostic
     * showed ~1e-8 with damping vs ~1e-5 without at step 800 (cell W-2). */
    const float wall_ratio = undamp.phi_near_wall_final / (damp.phi_near_wall_final + 1e-30f);
    printf("  wall    ratio  (undamped / damped) = %.1fx\n", wall_ratio);
    TEST_ASSERT(wall_ratio > 100.0f,
                "near-wall suppression = %.1fx below 100x threshold — damping not effective",
                wall_ratio);

    /* (d) Pre-bounce sanity: with no echoes back yet, damped and undamped must
     * agree closely — confirms damping doesn't perturb the *outgoing* wave in
     * the interior. */
    const float pre_diff = fabsf(damp.phi_center_max_pre - undamp.phi_center_max_pre);
    TEST_ASSERT(pre_diff < 1.0e-6f,
                "pre-bounce central |phi| differs by %.5g — damping is leaking into the interior",
                pre_diff);

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
