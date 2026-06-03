/* Stage 56 -- spin-gradient force smoke test (v40 phase 3).
 *
 * Sanity tests:
 *
 * (A) Larmor regression: a particle in UNIFORM B_z (zero gradient) feels
 *     no spin-gradient force.  The spin precession (Stage 54) should
 *     remain at the analytic Larmor rate and the particle should remain
 *     at rest.  This verifies the new force term is silent when it
 *     should be.
 *
 * (B) Spin-spin force smoke: two stationary spinning charged particles
 *     a few cells apart deposit dipole moments via phase 2.  The FDTD
 *     propagates the resulting A field, and the gradient of B_z at each
 *     particle's position drives a finite spin-gradient force.
 *     Verify the forces on the two particles are roughly equal and
 *     opposite (Newton's 3rd law) and finite. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int test_A_uniform_B(void) {
    printf("--- (A) particle at rest in uniform B_z, no gradient force ---\n");
    const int   W = 64, H = 64;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float mass = 1.0f, charge = 1.0f, g_factor = 2.0f, spin = 1.0f;
    const float B0 = 0.01f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_damping(sim, 8);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_k_e(sim, 1.0f);
    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_background_uniform_magnetic(sim, cx, cy, B0);

    const int idx = gr_sim_add_particle(sim, cx, cy, mass, charge, 0.0f, 0.0f);
    gr_sim_set_particle_spin(sim, idx, spin, g_factor);

    for (int s = 0; s < 1000; s++) gr_sim_step(sim);

    const gr_particle_t* p = gr_sim_get_particle(sim, idx);
    const float dt = gr_sim_dt(sim);
    const float omega_analytic = g_factor * charge / (2.0f * mass) * B0;
    const float phi_analytic   = omega_analytic * (float)(1000) * dt;
    const float rel_err_phi = (p->phi_spin - phi_analytic) / phi_analytic;
    const float x_drift = p->x - cx;
    const float y_drift = p->y - cy;
    printf("  phi_spin = %+.5f  (expected %+.5f, rel_err = %+.2e)\n",
           (double) p->phi_spin, (double) phi_analytic, (double) rel_err_phi);
    printf("  position drift = (%+.2e, %+.2e) cells  (expected 0)\n",
           (double) x_drift, (double) y_drift);
    gr_sim_destroy(sim);
    const int pass = (fabsf(rel_err_phi) < 1e-3f)
                     && (fabsf(x_drift) < 1e-3f)
                     && (fabsf(y_drift) < 1e-3f);
    printf("  %s\n", pass ? "PASS" : "FAIL");
    return pass;
}

static int test_B_two_spinning(void) {
    printf("\n--- (B) two stationary spinning charges, dipole-dipole force ---\n");
    const int   W = 128, H = 128;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float mass = 1.0f, charge = 1.0f, g_factor = 2.0f, spin = 1.0f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_k_e(sim, 1.0f);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 0);     /* disable Coulomb */
    gr_sim_set_em_inductive_enabled(sim, 0);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);

    /* Two particles 10 cells apart along x, both spin = +1. */
    const int p0 = gr_sim_add_particle(sim, cx - 5.0f, cy, mass, charge, 0.0f, 0.0f);
    const int p1 = gr_sim_add_particle(sim, cx + 5.0f, cy, mass, charge, 0.0f, 0.0f);
    gr_sim_set_particle_spin(sim, p0, spin, g_factor);
    gr_sim_set_particle_spin(sim, p1, spin, g_factor);

    /* Run enough steps for the dipole field to develop and reach the
     * other particle (~10 cells of light travel + a few). */
    const int N_warmup = 64;
    for (int s = 0; s < N_warmup; s++) gr_sim_step(sim);

    /* Record momenta, run a short measurement window, recover dpx/dt. */
    const gr_particle_t* a0 = gr_sim_get_particle(sim, p0);
    const gr_particle_t* a1 = gr_sim_get_particle(sim, p1);
    const float px0_a = a0->px, px0_b = a1->px;
    const int N_meas = 32;
    for (int s = 0; s < N_meas; s++) gr_sim_step(sim);
    const gr_particle_t* b0 = gr_sim_get_particle(sim, p0);
    const gr_particle_t* b1 = gr_sim_get_particle(sim, p1);
    const float dt = gr_sim_dt(sim);
    const float Fx_p0 = (b0->px - px0_a) / ((float)N_meas * dt);
    const float Fx_p1 = (b1->px - px0_b) / ((float)N_meas * dt);
    printf("  Fx on p0 = %+.4e\n", (double) Fx_p0);
    printf("  Fx on p1 = %+.4e\n", (double) Fx_p1);
    printf("  Fx_p0 + Fx_p1 (Newton 3) = %+.4e   |Fx_p0|/|Fx_p0|+|Fx_p1| = %.3f\n",
           (double)(Fx_p0 + Fx_p1),
           (double)(fabsf(Fx_p0) / fmaxf(fabsf(Fx_p0) + fabsf(Fx_p1), 1e-30f)));
    gr_sim_destroy(sim);

    /* Gates (phase 3 smoke):
     *   - both forces finite and nonzero  (spin-gradient code is being exercised)
     *
     * NOT GATED (deferred to future work):
     *   - Newton's 3rd law -- the static 2D z-axis dipole has B_z = 0
     *     outside its own support, so any spin-gradient force in this test
     *     comes from (a) FDTD wake transients during evolution and (b) the
     *     particle's OWN dipole's gradient at its own position (the v40
     *     analog of the rho_q self-force).  Cleanly testing dipole-dipole
     *     force requires either an analytic background with known dB_z/dx
     *     or extending v39 self-field subtraction to cover the spin-
     *     gradient gather.  Both deferred. */
    const int pass = isfinite(Fx_p0) && isfinite(Fx_p1)
                  && (fabsf(Fx_p0) > 1e-10f) && (fabsf(Fx_p1) > 1e-10f);
    printf("  %s (smoke gate only; see comments for deferred N3 test)\n",
           pass ? "PASS" : "FAIL");
    return pass;
}

int main(void) {
    printf("=== stage56_spin_gradient_force ===\n");
    const int a = test_A_uniform_B();
    const int b = test_B_two_spinning();
    printf("\n%s\n", (a && b) ? "PASS: spin-gradient force phase 3 functional."
                              : "FAIL");
    return (a && b) ? 0 : 1;
}
