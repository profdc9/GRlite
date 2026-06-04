/* Stage 61 -- EM field energy gravitates a separate test mass.
 *
 * Dynamic acceptance test for v41 phase 1: a static positive charge
 * creates a phi_em near-field; the E = -grad phi field has finite
 * energy density which (with em_stress_energy ON) contributes to
 * rho_eff and sources Phi_g via the gravitational wave equation.
 * A neutral test mass placed at a distance from the charge should
 * feel a gravitational attraction toward the EM-energy peak.
 *
 * Setup:
 *   - Static +Q charge at the box center (no mass; pure EM source).
 *   - Neutral test mass at (cx + r, cy), no charge, initially at rest.
 *   - G_eff large so the tiny EM-energy gravity is measurable.
 *   - k_e small so the EM near-field is a manageable magnitude.
 *
 * Predictions:
 *   (A) With em_stress_energy=OFF: nothing sources Phi_g (the test mass's
 *       own rho is at corner sublattice and produces a self-force we
 *       can subtract via self-field; or we simply note its dvx is
 *       dominated by self-force and roughly zero net inward).  Use
 *       self-field ON with eps=0 to cancel.
 *   (B) With em_stress_energy=ON: the EM-energy distribution near the
 *       +Q charge gravitates -- the test mass feels small inward
 *       attraction toward the charge.
 *
 * Gate: dvx_on < dvx_off  (with self-field on so own contribution is
 * subtracted) -- i.e., turning EM stress-energy ON adds a finite
 * inward (negative) acceleration. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float run(int em_stress_energy_on) {
    const int   W = 128, H = 128;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float mass    = 1.0f;
    const float r_off   = 15.0f;   /* test mass at +r from the static charge */

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_damping(sim, 16);
    gr_sim_set_G_eff(sim, 1.0f);
    gr_sim_set_k_e(sim, 1.0f);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* No EM force on the test particle; it has no charge anyway.  Keep the
     * gravitomagnetic gate ON so B_g (zero in this test) doesn't break
     * anything; standard force tier. */
    gr_sim_set_em_lorentz_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_em_stress_energy_enabled(sim, em_stress_energy_on);

    /* Static charge: add as a particle (heavy enough to stay at rest) with
     * essentially no mass impact, AT THE CENTER.  Charge but ~no mass: its
     * matter contribution to gravity is negligible compared to its EM-field
     * energy contribution. */
    const float src_mass = 1e-8f;          /* near-zero matter mass */
    const int idx_src = gr_sim_add_particle(sim, cx, cy,
                                             src_mass, /*charge=*/1.0f,
                                             0.0f, 0.0f);
    /* Test particle: neutral mass at distance r from the static charge. */
    const int idx_test = gr_sim_add_particle(sim,
                                              cx + r_off, cy,
                                              mass, /*charge=*/0.0f,
                                              0.0f, 0.0f);
    /* Self-field on the test mass so its own gravity self-force is
     * cancelled, isolating the cross contribution from the EM field
     * energy near the central charge. */
    gr_sim_particle_enable_self_field(sim, idx_test);
    gr_sim_particle_enable_self_field(sim, idx_src);

    const int n_warmup = 50;
    const int n_meas   = 200;

    for (int s = 0; s < n_warmup + n_meas; s++) gr_sim_step(sim);

    const gr_particle_t* p = gr_sim_get_particle(sim, idx_test);
    return p->px;
}

int main(void) {
    printf("=== stage61_em_energy_attraction ===\n");
    printf("Neutral test mass placed near a static +Q charge.\n");
    printf("With em_stress_energy OFF: pure self-field test mass should feel no force.\n");
    printf("With em_stress_energy ON:  EM-field energy gravitates; expect inward (-x) drift.\n\n");

    const float px_off = run(/*em_stress_energy_on=*/0);
    const float px_on  = run(/*em_stress_energy_on=*/1);
    printf("  test mass px after warmup+meas:\n");
    printf("    em_stress_energy OFF: %+11.4e\n", (double) px_off);
    printf("    em_stress_energy ON : %+11.4e\n", (double) px_on);
    printf("    difference (ON - OFF) = %+11.4e   (expect NEGATIVE -- inward force)\n",
           (double)(px_on - px_off));

    const int pass_attractive  = ((px_on - px_off) < 0.0f);
    const int pass_finite      = (fabsf(px_on - px_off) > 1e-10f);
    printf("\n  Sign: ON gives more negative px (inward force) : %s\n",
           pass_attractive ? "yes" : "NO");
    printf("  Magnitude finite                                : %s\n",
           pass_finite ? "yes" : "NO");

    const int pass = pass_attractive && pass_finite;
    printf("\n%s: EM field energy gravitates a separate test mass.\n",
           pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
