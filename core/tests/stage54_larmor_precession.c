/* Stage 54 -- Larmor precession of a charged particle's spin in uniform B.
 *
 * v40 spin phase 1 acceptance gate.  A charged spinning particle at rest
 * in uniform background B_z precesses its phi_spin at the Larmor rate:
 *
 *   d phi_spin/dt = (g_s q / 2m) B_z + B_g,z + Omega_Thomas,z
 *
 * Setup: particle at rest (no v x B force, no Omega_Thomas), no gravity,
 * no perturbation field evolution.  Only the analytic uniform B
 * background acts.  Expected: phi_spin(t) = (g_s q / 2m) B_z * t.
 *
 * Test compares measured phi_spin after N steps to the analytic value.
 * Gate: |Delta phi_spin / phi_analytic - 1| < 1e-3 (limited by the
 * leapfrog integration of a constant-rate ODE -- should be at float32
 * floor). */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== stage54_larmor_precession ===\n");
    printf("Charged spinning particle at rest in uniform B_z; expect\n");
    printf("phi_spin(t) = (g_s q / 2m) B_z t.\n\n");

    const int   W = 64, H = 64;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float mass     = 1.0f;
    const float charge   = 1.0f;
    const float g_factor = 2.0f;
    const float spin     = 0.5f;          /* magnitude; doesn't affect rate */
    const float B0       = 0.01f;         /* uniform B_z */

    /* Test a few step counts to confirm linear scaling. */
    const int Ns[] = {100, 1000, 10000};
    const int nN = (int)(sizeof(Ns)/sizeof(Ns[0]));

    for (int it = 0; it < nN; it++) {
        const int N = Ns[it];
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
        if (!sim) { fprintf(stderr, "sim_create failed\n"); return 1; }
        gr_sim_set_damping(sim, 8);
        gr_sim_set_G_eff(sim, 0.0f);
        gr_sim_set_k_e(sim, 1.0f);
        /* No perturbation field evolution -- uniform B is static background. */
        gr_sim_set_field_evolution(sim, 0);
        gr_sim_set_particle_source_deposition(sim, 0);
        /* All EM force pieces on (so v x B would matter if v != 0). */
        gr_sim_set_em_lorentz_force_enabled(sim, 1);
        gr_sim_set_em_electrostatic_enabled(sim, 1);
        gr_sim_set_em_magnetic_enabled(sim, 1);
        gr_sim_set_em_inductive_enabled(sim, 1);
        /* Gravity off (no B_g_z contribution). */
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

        gr_sim_set_background_uniform_magnetic(sim, cx, cy, B0);

        /* Particle at rest at the center. */
        const int idx = gr_sim_add_particle(sim, cx, cy, mass, charge, 0.0f, 0.0f);
        if (idx < 0) { fprintf(stderr, "add_particle failed\n"); gr_sim_destroy(sim); return 1; }
        gr_sim_set_particle_spin(sim, idx, spin, g_factor);

        const float dt = gr_sim_dt(sim);
        for (int s = 0; s < N; s++) gr_sim_step(sim);

        const gr_particle_t* p = gr_sim_get_particle(sim, idx);
        const float t = (float)N * dt;
        const float omega_analytic = g_factor * charge / (2.0f * mass) * B0;
        const float phi_analytic   = omega_analytic * t;
        const float rel_err        = (p->phi_spin - phi_analytic) / phi_analytic;
        printf("  N=%-6d  t=%9.4f  phi_spin = %+10.6f  expected = %+10.6f  rel_err = %+.3e\n",
               N, (double)t,
               (double)p->phi_spin,
               (double)phi_analytic,
               (double)rel_err);
        gr_sim_destroy(sim);

        if (fabsf(rel_err) > 1e-3f) {
            fprintf(stderr, "  FAIL: relative error exceeds 1e-3\n");
            return 1;
        }
    }
    printf("\nPASS: Larmor precession rate matches analytic to <1e-3.\n");
    return 0;
}
