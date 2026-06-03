/* Stage 58 -- Thomas precession of a charged spinning particle in
 * uniform B (cyclotron motion).
 *
 * A particle in circular motion accumulates a kinematic spin precession
 * from acceleration of its rest frame -- the Thomas precession.  Our
 * simulator (gr_sandbox_v38.tex sec:alg_2d step 11) integrates
 *
 *   d phi_spin / dt = (g_s q / 2m) B_z + Omega_Thomas,z
 *   Omega_Thomas,z  = -gamma^2/(gamma+1) (v_x a_y - v_y a_x) / c^2
 *
 * For a charged particle in uniform B_z, the Lorentz force produces
 * circular motion with cyclotron angular frequency
 *
 *   omega_c = q B / (gamma m)
 *
 * and (v x a)_z = v^2 omega_c (constant magnitude over the circular
 * orbit).  The combined precession rate is
 *
 *   d phi_spin/dt = (g_s/2 - 1) omega_c + omega_c / gamma
 *                 = omega_c [(g_s/2 - 1) + 1/gamma]
 *
 * which for g_s = 2 (Dirac) reduces to omega_c / gamma -- the BMT result.
 * Equivalently, after substituting v^2/c^2 = 1 - 1/gamma^2,
 *
 *   d phi_spin/dt = omega_c (1 + gamma - gamma^2) / gamma     [g_s = 2]
 *
 * This stage measures the actual phi_spin accumulation in a relativistic
 * cyclotron orbit and compares to that closed form. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int run_at(float v_init) {
    const int   W = 64, H = 64;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float mass     = 1.0f;
    const float charge   = 1.0f;
    const float g_factor = 2.0f;
    const float spin     = 1.0f;
    /* Use NEGATIVE B_z so q v x B with v = v_init y_hat gives -q v_init B_z x_hat
     * = +q v_init B0 x_hat, an inward (centripetal) force toward the cyclotron
     * center.  Choose B0 small enough that the orbital radius (r = gamma m v
     * / |q B|) fits inside the active grid region. */
    const float B0  = -0.05f;
    /* Run for ONE cyclotron period and compare phi_spin accumulation to
     * analytic.  Cyclotron rate: omega_c = |q B| / (gamma m).  Period
     * T = 2 pi / omega_c. */
    const float gamma_init = 1.0f / sqrtf(fmaxf(1.0f - v_init*v_init/(c_eff*c_eff), 1e-12f));
    const float omega_c    = fabsf(charge * B0) / (gamma_init * mass);
    const float T_orb      = 2.0f * (float)M_PI / omega_c;

    /* Analytic BMT precession rate for g=2 in 2D:
     *   d phi_spin/dt = omega_c (1 + gamma - gamma^2) / gamma
     * Sign: B0 < 0 so q v x B at v = +y_hat gives + x_hat force; the
     * cyclotron rotates CCW.  The Larmor term contributes
     *   (g q / 2m) B_z = q B_z / m = (charge * B0) / m
     * which is negative since B0 < 0.  Sign of phi_spin accumulation is
     * thus negative.  */
    const float phi_dot_analytic = (charge * B0) / mass
                                  - (gamma_init * gamma_init) / (gamma_init + 1.0f)
                                    * v_init * v_init * omega_c
                                    / (c_eff * c_eff);

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_damping(sim, 8);
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
    gr_sim_set_background_uniform_magnetic(sim, cx, cy, B0);

    const int idx = gr_sim_add_particle(sim, cx, cy, mass, charge, 0.0f, v_init);
    gr_sim_set_particle_spin(sim, idx, spin, g_factor);

    const float dt = gr_sim_dt(sim);
    const int N = (int) (T_orb / dt);
    for (int s = 0; s < N; s++) gr_sim_step(sim);

    const gr_particle_t* p = gr_sim_get_particle(sim, idx);
    const float t = (float) N * dt;
    const float phi_meas    = p->phi_spin;
    const float phi_analytic = phi_dot_analytic * t;
    const float rel_err     = (fabsf(phi_analytic) > 0.0f)
        ? (phi_meas - phi_analytic) / phi_analytic
        : (phi_meas - phi_analytic);

    printf("  v/c=%.3f  gamma=%.4f   omega_c=%.4f  T_orb=%.2f  N=%d\n",
           (double) v_init, (double) gamma_init,
           (double) omega_c, (double) T_orb, N);
    printf("    Larmor only (g=2)              : %+.4e\n",
           (double)((charge * B0) / mass) * (double) t);
    printf("    Thomas contribution            : %+.4e\n",
           (double)(- (gamma_init * gamma_init) / (gamma_init + 1.0f)
                    * v_init * v_init * omega_c / (c_eff * c_eff)) * (double) t);
    printf("    analytic total phi_spin(t)     : %+.4e\n", (double) phi_analytic);
    printf("    measured phi_spin(t)           : %+.4e   rel_err = %+.3e\n",
           (double) phi_meas, (double) rel_err);
    gr_sim_destroy(sim);
    return fabsf(rel_err) < 0.01f;
}

int main(void) {
    printf("=== stage58_thomas_precession ===\n");
    printf("Charged spinning particle in uniform B; expect phi_spin = Larmor + Thomas.\n\n");

    int all_pass = 1;
    /* Three velocities from mildly to substantially relativistic.  Thomas
     * scales as v^2/c^2 so it grows rapidly with v. */
    const float vs[] = {0.1f, 0.3f, 0.5f};
    for (int i = 0; i < 3; i++) {
        printf("\n");
        if (!run_at(vs[i])) all_pass = 0;
    }
    printf("\n%s: phi_spin matches analytic Larmor+Thomas %s 1%%.\n",
           all_pass ? "PASS" : "FAIL", all_pass ? "to <" : "miss");
    return all_pass ? 0 : 1;
}
