/* Stage 60 -- EM effective GEM sources (gr_sandbox_v38.tex sec:alg_2d step 3).
 *
 * Verifies the v41 phase 1 implementation: the EM perturbation field's
 * energy density and Poynting flux are added to rho_matter / J_m before
 * each GEM wave step, so the EM field gravitates.
 *
 * Test design: deposit a STATIC charge distribution that generates a
 * known phi_em near-field.  With EM stress-energy disabled, Phi_g
 * stays at zero (no matter, no current).  With EM stress-energy
 * enabled, the |E|^2 term in rho_EM = eps_0 / (2 c^2) |E|^2 + ... acts
 * as a positive gravitating source -- Phi_g should acquire a small
 * NEGATIVE value (since the gravity source coupling is
 *   sc_grav = -4 pi G_eff
 * giving Lap Phi_g = +4 pi G_eff rho_eff in static limit, so for
 * rho_eff > 0 Phi_g has a maximum where the source is -- but in our
 * convention "attractive" gravity has Phi_g > 0 at the source location
 * for positive masses, so EM stress-energy as positive source produces
 * Phi_g > 0 at the EM-energy peak).
 *
 * Gate: with stress-energy ON, Phi_g at the EM-energy peak is nonzero
 * and SAME sign as the matter-induced Phi_g would be for a positive
 * mass deposit at the same location. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float Phi_g_at_center(gr_sim_t* sim, int W, int H) {
    const float* phi = sim->fields[GR_FIELD_PHI_GRAV].curr;
    const int j = H / 2, i = W / 2;
    return phi[j * W + i];
}

int main(void) {
    printf("=== stage60_em_stress_energy ===\n");
    printf("Verify that EM field energy sources Phi_g (linearized-GR coupling).\n\n");

    const int   W = 64, H = 64;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;

    /* Same setup twice -- once with stress-energy OFF, once ON. */
    float Phi_g_off = 0.0f, Phi_g_on  = 0.0f;
    for (int run = 0; run < 2; run++) {
        gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
        gr_sim_set_damping(sim, 16);
        gr_sim_set_G_eff(sim, 1.0f);          /* large G_eff to make Phi_g detectable */
        gr_sim_set_k_e(sim, 1.0f);
        gr_sim_set_field_evolution(sim, 1);
        gr_sim_set_particle_source_deposition(sim, 0);
        /* No EM force on the particle -- we just want to deposit charge
         * and let the EM field develop. */
        gr_sim_set_em_lorentz_force_enabled(sim, 1);
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_em_stress_energy_enabled(sim, run == 1 ? 1 : 0);

        /* Deposit a static positive charge at the center to build up an
         * E field there.  Use the API for a direct point-charge deposit. */
        gr_sim_clear_sources(sim);
        gr_sim_deposit_point_charge(sim, cx, cy, 1.0f);

        /* Run several steps to let phi_em propagate (~10 cells of light
         * travel = ~14 sim units = 20 steps). */
        for (int s = 0; s < 100; s++) gr_sim_step(sim);

        const float phi_g = Phi_g_at_center(sim, W, H);
        if (run == 0) Phi_g_off = phi_g;
        else          Phi_g_on  = phi_g;
        gr_sim_destroy(sim);
    }

    printf("  Phi_g at center, EM stress-energy OFF: %+.6e\n", (double) Phi_g_off);
    printf("  Phi_g at center, EM stress-energy ON : %+.6e\n", (double) Phi_g_on);
    printf("  Difference (EM-gravity coupling)     : %+.6e\n",
           (double)(Phi_g_on - Phi_g_off));

    /* Sanity checks:
     *   (1) Phi_g_off should be ~zero (no mass deposited, only charge).
     *   (2) Phi_g_on should be DIFFERENT from Phi_g_off by a finite amount
     *       (the EM-gravity coupling is doing something).
     *   (3) Sign: in 2D linearized gravity, a positive mass at origin
     *       gives Phi_g = 2 G_eff m ln(r), which goes to -infinity at the
     *       source (ln r -> -inf as r -> 0).  Force is -m grad Phi_g; for
     *       attractive force inward Phi_g must increase outward, i.e.,
     *       Phi_g is LARGE NEGATIVE at the source.  Adding positive
     *       rho_EM makes Phi_g MORE negative at the source.  So we
     *       expect Phi_g_on < Phi_g_off. */
    const int pass_off_small  = (fabsf(Phi_g_off) < 1e-10f);
    const int pass_on_nonzero = (fabsf(Phi_g_on)  > 10.0f * fabsf(Phi_g_off) + 1e-12f);
    const int pass_sign       = (Phi_g_on < Phi_g_off);
    printf("\n  Phi_g_off near zero (no matter)   : %s\n", pass_off_small ? "yes" : "NO");
    printf("  Phi_g_on differs from off         : %s\n", pass_on_nonzero ? "yes" : "NO");
    printf("  Sign: EM energy decreases Phi_g  : %s\n", pass_sign ? "yes" : "NO");

    const int all_pass = pass_off_small && pass_on_nonzero && pass_sign;
    printf("\n%s: EM stress-energy sources Phi_g as expected.\n",
           all_pass ? "PASS" : "FAIL");
    return all_pass ? 0 : 1;
}
