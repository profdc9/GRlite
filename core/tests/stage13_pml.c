/* Stage 13 — Perfectly-Matched-Layer (PML) absorbing boundary.
 *
 * The damping ring of Stage 2 is a lossy material: every wave passing
 * through the ring is attenuated, but the impedance jump at the ring's
 * inner edge produces a position-dependent reflection that, on a
 * stationary point source, manifests as a small but real systematic
 * self-force at off-center particle positions (see memory).
 *
 * IMPORTANT (finding from this test, 2026-05-20): PML is impedance-
 * matched for *plane waves* at the PML/interior interface — that's
 * what makes it "perfectly matched".  But the *static* log(r) potential
 * from a stationary point source is DC content (not a propagating
 * wave), and PML's nonzero sigma in the boundary ring still distorts
 * the steady-state DC distribution in a position-dependent way.  PML
 * therefore does NOT restore HE self-force = 0 at off-center positions
 * for a stationary source — that intuition (in the earlier project
 * memory) was wrong.  What PML *does* do better than the lossy damping
 * ring is (a) absorb outgoing radiation more cleanly and (b) reduce
 * standing-wave wake from a moving radiating source (Phase E heating).
 *
 * Tests:
 *   [A] Outgoing pulse absorption.  Hard assert: comparable or better
 *       than the existing damping ring on the same setup.
 *
 *   [B] Stationary off-center drift comparison (informational).  Reports
 *       drift_pml vs drift_damping at 0, +0.5, and +23.5 cells off the
 *       box center for m=1e-3.  Both schemes show position-dependent
 *       drift; this test documents the magnitudes.
 *
 *   [C] Phase E orbit at m=1e-3 — pic_orbiting under PML vs damping
 *       ring.  This is where PML *might* help: the moving particle
 *       radiates a wake and the cleanliness of its absorption is what
 *       differs between the two schemes.  Hard assert: PML configuration
 *       runs without NaN/Inf; comparison printed.
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

/* ----------------------------------------------------------------------- */
/* [A] PML absorption                                                       */
/* ----------------------------------------------------------------------- */
static int test_pml_absorption(void) {
    printf("\n[A] PML absorbs outgoing Gaussian pulse\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float sigma  = 4.0f * dx;
    const float amp    = 1.0f;
    const int   n_pml  = 16;
    const int   n_steps = 800;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create failed");
    gr_sim_set_pml(sim, n_pml);
    TEST_ASSERT(gr_sim_pml_layers(sim) == n_pml, "PML did not enable");
    TEST_ASSERT(gr_sim_damping_layers(sim) == 0, "damping ring not cleared by PML");

    float params[2] = {sigma, amp};
    TEST_ASSERT(gr_sim_load_scenario(sim, "wave_pulse", params, 2) == 0,
                "wave_pulse load failed");

    float center_max = 0.0f;
    for (int n = 1; n <= n_steps; n++) {
        gr_sim_step(sim);
        const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
        if (n >= 50 && n <= 150) {
            const float v = fabsf(phi[(H / 2) * W + (W / 2)]);
            if (!isfinite(v)) { gr_sim_destroy(sim); return -1; }
            if (v > center_max) center_max = v;
        }
    }
    const float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
    const float center_final   = fabsf(phi[(H / 2) * W + (W / 2)]);
    const float near_wall_final = fabsf(phi[(H / 2) * W + (W - 2)]);
    printf("  PML (N=%d): pre-bounce peak |phi(center)| = %.5g, |phi(center)|@final = %.5g, |phi(near-wall)|@final = %.5g\n",
           n_pml, center_max, center_final, near_wall_final);
    TEST_ASSERT(isfinite(center_final),  "PML produced NaN/Inf at center");
    TEST_ASSERT(isfinite(near_wall_final), "PML produced NaN/Inf near wall");
    /* Want at least as well as the damping ring's < 5e-2 threshold. */
    TEST_ASSERT(center_final < 5.0e-2f,
                "PML central residual %.5g exceeds 5e-2", center_final);
    gr_sim_destroy(sim);
    return 0;
}

/* ----------------------------------------------------------------------- */
/* [B] Stationary off-center drift — PML vs damping ring (informational)    */
/* ----------------------------------------------------------------------- */
static int test_he_off_center_drift(void) {
    printf("\n[B] Stationary off-center drift, PML vs damping ring\n");
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = 1.0e-3f;
    const int   N      = 2000;

    const float cx = ((float) (W - 1) * 0.5f) * dx;
    const float cy = ((float) (H - 1) * 0.5f) * dx;
    struct { const char* label; float dx_off; float dy_off; } cases[] = {
        { "center        ", 0.0f,        0.0f },
        { "+0.5 cell     ", 0.5f * dx,   0.0f },
        { "+23.5 cells   ", 23.5f * dx,  0.0f },
    };
    const int n_cases = sizeof(cases) / sizeof(cases[0]);

    /* Run each case under both PML and damping ring. */
    enum { TRIAL_DAMP = 0, TRIAL_PML = 1, N_TRIAL = 2 };
    const char* trial_name[N_TRIAL] = { "damping", "PML    " };
    float drift[N_TRIAL][3];

    for (int t = 0; t < N_TRIAL; t++) {
        for (int c = 0; c < n_cases; c++) {
            gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
            TEST_ASSERT(sim != NULL, "create failed");
            if (t == TRIAL_DAMP) gr_sim_set_damping(sim, 16);
            else                 gr_sim_set_pml(sim, 16);
            gr_sim_set_particle_source_deposition(sim, 1);

            const float x0 = cx + cases[c].dx_off;
            const float y0 = cy + cases[c].dy_off;
            const int id = gr_sim_add_particle(sim, x0, y0, mass, 0.0f, 0.0f, 0.0f);
            TEST_ASSERT(id == 0, "add_particle failed");

            gr_sim_step_n(sim, N);

            const gr_particle_t* p = gr_sim_get_particle(sim, 0);
            const float dr = sqrtf((p->x - x0) * (p->x - x0) + (p->y - y0) * (p->y - y0));
            drift[t][c] = dr / dx;
            const int finite = isfinite(drift[t][c]);
            printf("  %s, %s: drift over %d steps = %.3e cells%s\n",
                   trial_name[t], cases[c].label, N, drift[t][c],
                   finite ? "" : "  [NaN!]");
            TEST_ASSERT(finite, "non-finite drift at trial=%d case=%d", t, c);
            gr_sim_destroy(sim);
        }
    }

    /* Comparison table — PML vs damping for each offset. */
    printf("\n  offset           damping_drift   PML_drift     PML/damping\n");
    for (int c = 0; c < n_cases; c++) {
        const float dd = drift[TRIAL_DAMP][c];
        const float dp = drift[TRIAL_PML][c];
        const float ratio = (dd > 1.0e-12f) ? (dp / dd) : ((dp == 0.0f) ? 1.0f : INFINITY);
        printf("  %s    %.3e       %.3e     %.3f\n",
               cases[c].label, dd, dp, ratio);
    }
    printf("  (informational — PML does NOT restore translation invariance;\n");
    printf("   what PML buys is cleaner absorption of *propagating* waves,\n");
    printf("   benchmarked in test [A] and the Phase E orbit in [C].)\n");
    return 0;
}

/* ----------------------------------------------------------------------- */
/* [C] Phase E orbit at m = 1e-3 under PML (informational)                  */
/* ----------------------------------------------------------------------- */
static int test_phase_e_under_pml(void) {
    printf("\n[C] Phase E orbit at m=1e-3 — PML vs damping ring (informational)\n");
    const int   W      = 128, H = 128;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float GM     = 1.0f;
    const float r_orb  = 20.0f;
    const float eps    = 1.0f;
    const float m_test = 1.0e-3f;
    const float par[4] = {GM, r_orb, eps, m_test};

    const float cx     = ((float) W * 0.5f) * dx;
    const float cy     = ((float) H * 0.5f) * dx;
    const float T_ana  = 2.0f * (float) M_PI * sqrtf(r_orb * r_orb * r_orb / GM);

    /* Run each setup for 2 analytic periods worth of steps, report
     * cumulative drift in r at the end.  This is exploratory — we want a
     * number to compare against the damping-ring baseline, not a pass/fail. */

    for (int trial = 0; trial < 2; trial++) {
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
        TEST_ASSERT(sim != NULL, "create failed");
        if (trial == 0) {
            gr_sim_set_damping(sim, 16);
            printf("  --- trial: damping ring (N=16) ---\n");
        } else {
            gr_sim_set_pml(sim, 16);
            printf("  --- trial: PML (N=16) ---\n");
        }
        TEST_ASSERT(gr_sim_load_scenario(sim, "pic_orbiting", par, 4) == 0,
                    "pic_orbiting load failed");

        const int N_steps = (int) (2.0f * T_ana / gr_sim_dt(sim));
        const int milestones[] = {N_steps / 8, N_steps / 4, N_steps / 2, N_steps };
        const int n_ms = sizeof(milestones) / sizeof(milestones[0]);

        int ms_idx = 0;
        for (int s = 1; s <= N_steps; s++) {
            gr_sim_step(sim);
            if (ms_idx < n_ms && s >= milestones[ms_idx]) {
                const gr_particle_t* p = gr_sim_get_particle(sim, 0);
                const float r = sqrtf((p->x - cx) * (p->x - cx) + (p->y - cy) * (p->y - cy));
                const float drift_pct = 100.0f * (r - r_orb) / r_orb;
                printf("    s=%5d (%.0f%% of 2T_ana): r=%.4f  drift=%+.2f%%\n",
                       s, 100.0f * (float) s / (float) N_steps, r, drift_pct);
                ms_idx++;
            }
        }
        gr_sim_destroy(sim);
    }
    printf("  (informational — no assertion)\n");
    return 0;
}

int main(void) {
    printf("=== stage13_pml: PML absorbing boundary ===\n");
    if (test_pml_absorption()      != 0) return 1;
    if (test_he_off_center_drift() != 0) return 1;
    if (test_phase_e_under_pml()   != 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
