/* Stage 32b -- EM sign convention verification.
 *
 * After the v36 EM source-sign flip (sc_em = +4 pi k_e), this stage
 * verifies the two unambiguous sign tests on PIC EM:
 *
 *   [A] Like charges REPEL (Coulomb sign of PIC phi_em).
 *   [B] Opposite charges ATTRACT (mirror of [A] -- catches sign convention
 *       errors that flip both [A] and [B] together).
 *   [C] phi_em profile from a single +Q source decreases monotonically
 *       with distance (verifies the SHAPE of the box-bounded 2D Green's
 *       function -- phi positive everywhere in the damped box, peaked
 *       at source).
 *
 * Why no v x B PIC-Ampere dynamics test here:  the magnetic force is a
 * (v/c)^2 correction to Coulomb (~9% at v/c=0.3), much smaller than
 * the v-dependent PIC moving-source self-heating + Cherenkov-like
 * artifacts at non-trivial v.  Cleanly isolating the small Ampere
 * signal from the large PIC artifacts requires further investigation
 * (a future stage).  The v x B FORCE LAW itself is already verified
 * with the analytic uniform-B background in Stage 23 (cyclotron). */

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

static float read_corner(const float* arr, int W, int i, int j) {
    return arr[j * W + i];
}

/* Run two static charges (q1, q2) at separation d_sep; return widening. */
static float static_pair_widening(float q1, float q2, float mass,
                                  float d_sep, int n_steps) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return NAN;
    gr_sim_set_damping(sim, 16);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    const float cx = ((float) (W - 1) * 0.5f) * dx;
    const float cy = ((float) (H - 1) * 0.5f) * dx;
    gr_sim_add_particle(sim, cx - d_sep * 0.5f, cy, mass, q1, 0.0f, 0.0f);
    gr_sim_add_particle(sim, cx + d_sep * 0.5f, cy, mass, q2, 0.0f, 0.0f);
    for (int s = 0; s < n_steps; s++) gr_sim_step(sim);
    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
    const float r_final = fabsf(p1->x - p0->x);
    gr_sim_destroy(sim);
    return r_final - d_sep;
}

static int test_static_coulomb_signs(void) {
    printf("\n[A,B] Static Coulomb sign tests:\n");
    const float Q = 0.01f, mass = 0.01f, d_sep = 16.0f;
    const int   n_steps = 80;

    const float w_like     = static_pair_widening(+Q, +Q, mass, d_sep, n_steps);
    const float w_opposite = static_pair_widening(+Q, -Q, mass, d_sep, n_steps);
    printf("  like     +Q,+Q: widening = %+.4f (expect > 0 -- repulsion)\n", w_like);
    printf("  opposite +Q,-Q: widening = %+.4f (expect < 0 -- attraction)\n", w_opposite);

    TEST_ASSERT(w_like > 0.0f,
                "Like charges did not repel: %+.3f", w_like);
    TEST_ASSERT(w_opposite < 0.0f,
                "Opposite charges did not attract: %+.3f", w_opposite);
    return 0;
}

static int test_phi_em_profile(void) {
    printf("\n[C] PIC phi_em profile around a single static +Q.\n");
    printf("    Expect monotonically DECREASING with distance from source.\n");
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float Q      = 0.01f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    TEST_ASSERT(sim != NULL, "create");
    gr_sim_set_damping(sim, 16);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    const int   ci = W / 2;
    const int   cj = H / 2;
    const float cx = (float) ci * dx;
    const float cy = (float) cj * dx;
    gr_sim_add_particle(sim, cx, cy, 1.0e6f, +Q, 0.0f, 0.0f);   /* heavy, static */
    for (int s = 0; s < 80; s++) gr_sim_step(sim);

    const float* phi_arr = sim->fields[GR_FIELD_PHI_EM].curr;
    const int dxs[] = {0, 4, 8, 12, 16, 24, 32};
    const int n_dx = (int)(sizeof(dxs)/sizeof(dxs[0]));
    float phi_at[7];
    printf("    %-6s %-12s\n", "+dx", "phi_em");
    for (int k = 0; k < n_dx; k++) {
        phi_at[k] = read_corner(phi_arr, W, ci + dxs[k], cj);
        printf("    %-6d %+11.5e\n", dxs[k], phi_at[k]);
    }

    TEST_ASSERT(phi_at[0] > 0.0f,
                "phi_em at source should be POSITIVE for +Q, got %+e",
                phi_at[0]);
    /* Strict monotonic decrease. */
    for (int k = 1; k < n_dx; k++) {
        TEST_ASSERT(phi_at[k] < phi_at[k - 1],
                    "phi_em not monotonically decreasing at +dx=%d (%.4e >= %.4e)",
                    dxs[k], phi_at[k], phi_at[k - 1]);
    }
    printf("    => phi_em peaks at source and decreases with distance. CORRECT.\n");
    gr_sim_destroy(sim);
    return 0;
}

int main(void) {
    printf("=== stage32b_em_sign_verification ===\n");
    printf("Verifies the v36 c_em sign flip puts PIC EM in standard Maxwell.\n");
    if (test_static_coulomb_signs() != 0) return 1;
    if (test_phi_em_profile()       != 0) return 1;
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
