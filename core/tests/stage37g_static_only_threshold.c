/* Stage 37g -- does the v-dependent threshold survive without
 * the relativistic (velocity-dependent) force terms?
 *
 * Dan, 2026-06-01: "What if we turn off the relativistic force terms
 * to see if that changes the v dependence?"
 *
 * EM Lorentz: F = q*(-grad phi - d_t A + v x B)
 *   Static (Coulomb-like):     -q grad phi
 *   Relativistic / vel-dep:    -q d_t A,  +q v x B
 *
 * GR (Newtonian tier + GEM): F = -m grad Phi_g + (-m d_t A_g) + 4 m v x B_g
 *   Static:                    -m grad Phi_g
 *   Velocity-dependent:        -m d_t A_g,  +4 m v x B_g
 *
 * If the heating/drag sign flips in F(v) come from -q grad phi alone,
 * disabling the velocity-dependent terms should KEEP the same threshold
 * pattern.  If they come from d_t A / v x B, the threshold should
 * disappear or shift.
 *
 * Single-particle, equal-distance methodology, mass=1, BUMP R=8.
 * G_eff*m^2 = k_e*Q^2 = 1e-4 (same coupling magnitude as 37e/f). */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px*p->px + p->py*p->py) / (m*m*c_eff*c_eff));
    return p->py / (gamma * m);
}

static float measure(int is_em, int veldep_on, float v_drift) {
    const int W = 128, H = 1024, n_damp = 16;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass = 1.0f, Q = 0.01f;
    const float y_start = (float)(n_damp + 32) * dx;
    const float dt = cfl * dx / c_eff;
    const int n_warmup = (int) ceilf(8.0f / (v_drift * dt));
    const int n_measure = (int) ceilf(2.0f / (v_drift * dt));

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 8.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);

    if (is_em) {
        gr_sim_set_G_eff(sim, 0.0f);
        gr_sim_set_k_e(sim, 1.0f);
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
        gr_sim_set_em_lorentz_force_enabled(sim, 1);
        gr_sim_set_em_electrostatic_enabled(sim, 1);
        gr_sim_set_em_magnetic_enabled    (sim, veldep_on);
        gr_sim_set_em_inductive_enabled   (sim, veldep_on);
        gr_sim_add_particle(sim, ((float)(W-1)*0.5f)*dx, y_start, mass, Q, 0.0f, v_drift);
    } else {
        gr_sim_set_G_eff(sim, 1.0e-4f);
        gr_sim_set_k_e(sim, 0.0f);
        gr_sim_set_em_lorentz_force_enabled(sim, 0);
        gr_sim_set_em_inductive_enabled(sim, 0);
        gr_sim_set_em_electrostatic_enabled(sim, 0);
        gr_sim_set_em_magnetic_enabled(sim, 0);
        gr_sim_set_gravitomagnetic_force_enabled    (sim, veldep_on);
        gr_sim_set_gravitomagnetic_inductive_enabled(sim, veldep_on);
        gr_sim_add_particle(sim, ((float)(W-1)*0.5f)*dx, y_start, mass, /*charge=*/0.0f, 0.0f, v_drift);
    }

    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);
    const float vy_w = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    for (int s = 0; s < n_measure; s++) gr_sim_step(sim);
    const float vy_m = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    const float F = mass * (vy_m - vy_w) / ((float)n_measure * dt);
    gr_sim_destroy(sim);
    return F;
}

int main(void) {
    printf("=== stage37g_static_only_threshold ===\n");
    printf("Static-only force (grad phi or grad Phi_g) vs full force.\n");
    printf("BUMP R=8, mass=1, equal-distance 8-cell warmup + 2-cell measurement.\n\n");

    const float vs[] = {0.005f, 0.010f, 0.015f, 0.020f, 0.025f, 0.030f,
                         0.035f, 0.040f, 0.045f, 0.050f, 0.060f, 0.075f, 0.100f};
    const int Nv = sizeof(vs)/sizeof(vs[0]);

    printf("EM single +Q charge:\n");
    printf("  v/c     full          static-only   diff (vel-dep)\n");
    printf("  -----------------------------------------------------\n");
    for (int iv = 0; iv < Nv; iv++) {
        const float F_full   = measure(/*is_em=*/1, /*veldep=*/1, vs[iv]);
        const float F_static = measure(/*is_em=*/1, /*veldep=*/0, vs[iv]);
        printf("  %.3f   %+11.4e   %+11.4e   %+11.4e\n",
               (double)vs[iv], (double)F_full, (double)F_static, (double)(F_full - F_static));
    }
    printf("\n");

    printf("GR single mass:\n");
    printf("  v/c     full          static-only   diff (vel-dep)\n");
    printf("  -----------------------------------------------------\n");
    for (int iv = 0; iv < Nv; iv++) {
        const float F_full   = measure(/*is_em=*/0, /*veldep=*/1, vs[iv]);
        const float F_static = measure(/*is_em=*/0, /*veldep=*/0, vs[iv]);
        printf("  %.3f   %+11.4e   %+11.4e   %+11.4e\n",
               (double)vs[iv], (double)F_full, (double)F_static, (double)(F_full - F_static));
    }
    return 0;
}
