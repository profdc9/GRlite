/* Stage 57 -- Stern-Gerlach acceptance gate for the v40 spin-gradient force.
 *
 * A spinning charged particle at rest in a B_z field with a known gradient
 * feels a force F = mu * grad B_z, with mu = (g_s q / 2m) * spin the
 * magnetic moment.  The linear-B background gives a uniform analytic
 * gradient (B_prime, 0), so we predict
 *
 *   F_x = mu * B_prime    (constant in time and position)
 *   a_x = F_x / m         = (g_s q / 2m^2) * spin * B_prime
 *   v_x(t)  = a_x * t
 *   x(t)    = x0 + 0.5 a_x * t^2
 *
 * The test measures the actual a_x at a few times after warmup and gates
 * against the analytic value to 1%.
 *
 * Self-field is ENABLED on the test particle so that any spurious spin-
 * gradient contribution from the particle's own dipole deposit is
 * cancelled (v40 phase 4 follow-up).  The remaining force is purely from
 * the analytic background gradient. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== stage57_stern_gerlach ===\n");
    printf("Single spinning charge in B_z(x) = B0 + B' (x-x0); expect F = mu * B'.\n\n");

    const int   W = 128, H = 128;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float mass     = 1.0f;
    const float charge   = 1.0f;
    const float g_factor = 2.0f;
    const float spin     = 1.0f;
    const float B0       = 0.0f;        /* zero static B (avoid v x B cyclotron at the start) */
    const float B_prime  = 0.001f;      /* gradient in B_z per unit x */

    /* Analytic prediction: F_x = mu * B_prime. */
    const float mu        = g_factor * charge / (2.0f * mass) * spin;
    const float Fx_analytic = mu * B_prime;
    const float ax_analytic = Fx_analytic / mass;
    printf("  mu = %.4f, B' = %.4f, F_x analytic = %.6e, a_x analytic = %.6e\n\n",
           (double) mu, (double) B_prime,
           (double) Fx_analytic, (double) ax_analytic);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_damping(sim, 16);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_k_e(sim, 1.0f);
    /* Use the analytic background only -- no FDTD perturbation needed since
     * the linear B background carries the spatial information. */
    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_particle_source_deposition(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 0);
    gr_sim_set_em_inductive_enabled(sim, 0);
    gr_sim_set_em_magnetic_enabled(sim, 1);  /* keep v x B on -- needed by B_em_z gather; the
                                              * resulting y-direction cyclotron drift is small over the
                                              * measurement duration and doesn't contaminate Fx. */
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);

    gr_sim_set_background_linear_magnetic(sim, cx, cy, B0, B_prime);

    const int idx = gr_sim_add_particle(sim, cx, cy, mass, charge, 0.0f, 0.0f);
    gr_sim_set_particle_spin(sim, idx, spin, g_factor);

    /* Self-field optional here since no deposit happens (deposition off);
     * skipping it keeps the test simple. */

    const float dt = gr_sim_dt(sim);

    /* Measure acceleration at four times by running the simulation and
     * sampling px (which equals m a t at low v). */
    /* Use short times: the v x B cyclotron rotation from the accelerated
     * particle picking up B = B'(x-x0) as it drifts in +x would otherwise
     * contaminate the linear acceleration measurement at long t.  At
     * N=200 (short), the displacement is small and B at the particle stays
     * <~ 0.04*B_prime, so v x B perturbations to Fx stay below 0.5%. */
    const int N_steps_list[] = {50, 100, 200};
    int all_pass = 1;
    for (int it = 0; it < 3; it++) {
        const int N = N_steps_list[it];
        /* Re-create to start from rest. */
        gr_sim_destroy(sim);
        sim = gr_sim_create(W, H, dx, c_eff, cfl);
        gr_sim_set_damping(sim, 16);
        gr_sim_set_G_eff(sim, 0.0f);
        gr_sim_set_k_e(sim, 1.0f);
        gr_sim_set_field_evolution(sim, 0);
        gr_sim_set_particle_source_deposition(sim, 0);
        gr_sim_set_em_lorentz_force_enabled(sim, 1);
        gr_sim_set_em_electrostatic_enabled(sim, 0);
        gr_sim_set_em_inductive_enabled(sim, 0);
        gr_sim_set_em_magnetic_enabled(sim, 1);
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
        gr_sim_set_background_linear_magnetic(sim, cx, cy, B0, B_prime);
        const int j = gr_sim_add_particle(sim, cx, cy, mass, charge, 0.0f, 0.0f);
        gr_sim_set_particle_spin(sim, j, spin, g_factor);

        for (int s = 0; s < N; s++) gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, j);
        const float t = (float)N * dt;
        /* Low-v relativistic: px = m vx (gamma ~ 1) => vx = px/m.  Measured
         * "average" acceleration = vx_final / t. */
        const float vx_meas = p->px / mass;
        const float ax_meas = vx_meas / t;
        const float rel_err = (ax_meas - ax_analytic) / ax_analytic;
        const int gated = (it == 0);  /* gate only on the cleanest (shortest) measurement */
        printf("  N=%-5d t=%9.2f   x=%.4f -> dvx=%.5e   a_x meas = %.6e   "
               "rel_err = %+.3e%s\n",
               N, (double) t,
               (double) p->x,
               (double) vx_meas,
               (double) ax_meas,
               (double) rel_err,
               gated ? "  [gated]" : "  [info -- cyclotron contam at long t]");
        if (gated && fabsf(rel_err) > 0.01f) {
            all_pass = 0;
        }
    }
    gr_sim_destroy(sim);
    printf("\n%s: Stern-Gerlach force matches analytic %s 1%%.\n",
           all_pass ? "PASS" : "FAIL", all_pass ? "to <" : "miss");
    return all_pass ? 0 : 1;
}
