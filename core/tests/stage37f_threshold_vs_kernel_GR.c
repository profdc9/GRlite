/* Stage 37f -- GR analog of Stage 37e.
 *
 * Dan, 2026-06-01: "Try the same on GR to see if we have just one
 * sign of self-force or it inverts at different v/c for different
 * kernels."
 *
 * Setup mirrors Stage 37e but with EM disabled, gravity enabled:
 *   mass = 1, charge = 0
 *   G_eff = 1e-4 (so G_eff * m^2 = k_e * Q^2 = 1e-4 -- same coupling
 *       magnitude as the EM test for direct comparison)
 *   gravitomagnetic force ON (-4 m v x B_g, GEM spin-2 enhancement)
 *   gravitomagnetic inductive ON (-m d_t A_g)
 *   em_lorentz_force_enabled = 0 (no EM)
 *
 * Note: gravity has spin-2 +4x coefficient on v x B_g vs EM's +1x on
 * v x B.  So the velocity-dependent piece is 4x larger in GR.  If the
 * spurious-force sign flips are dominated by the v x B contribution
 * this could shift the threshold; if they come from the static
 * grad-phi/dt-A discretization the threshold should match EM.
 *
 * Equal-distance methodology (8-cell warmup + 2-cell measurement),
 * mass = 1 so v changes little during measurement.  CFL = 1/sqrt2. */

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

static float measure_F(gr_shape_function_t shape, float Rk, int smooth, float v_drift) {
    const int W = 128, H = 1024, n_damp = 16;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass = 1.0f;
    const float y_start = (float)(n_damp + 32) * dx;
    const float dt = cfl * dx / c_eff;
    const int n_warmup = (int) ceilf(8.0f / (v_drift * dt));
    const int n_measure = (int) ceilf(2.0f / (v_drift * dt));

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, smooth);
    gr_sim_set_j_smooth_passes(sim, smooth);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* GR setup: gravity ON, EM OFF.  G_eff*m^2 = 1e-4 = same coupling
     * magnitude as EM (k_e*Q^2) in stage 37e. */
    gr_sim_set_G_eff(sim, 1.0e-4f);
    gr_sim_set_k_e(sim, 0.0f);
    gr_sim_set_em_lorentz_force_enabled(sim, 0);
    gr_sim_set_em_inductive_enabled(sim, 0);
    gr_sim_set_em_electrostatic_enabled(sim, 0);
    gr_sim_set_em_magnetic_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 1);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 1);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);

    const float cx = ((float)(W-1) * 0.5f) * dx;
    gr_sim_add_particle(sim, cx, y_start, mass, /*charge=*/0.0f, 0.0f, v_drift);
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);
    const float vy_w = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    for (int s = 0; s < n_measure; s++) gr_sim_step(sim);
    const float vy_m = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    const float F = mass * (vy_m - vy_w) / ((float)n_measure * dt);
    gr_sim_destroy(sim);
    return F;
}

typedef struct {
    const char* label;
    gr_shape_function_t shape;
    float Rk;
    int   smooth;
} kernel_t;

int main(void) {
    printf("=== stage37f_threshold_vs_kernel_GR ===\n");
    printf("GR analog: single mass (no charge) at constant v in self-gravity field.\n");
    printf("G_eff=1e-4, m=1, all gravity force pieces ON, EM OFF.\n");
    printf("Equal-distance 8-cell warmup + 2-cell measurement.\n\n");

    const float vs[] = {0.005f, 0.010f, 0.015f, 0.020f, 0.025f, 0.030f,
                         0.035f, 0.040f, 0.045f, 0.050f, 0.060f, 0.075f, 0.100f};
    const int Nv = sizeof(vs)/sizeof(vs[0]);

    const kernel_t kernels[] = {
        {"TSC bare",   GR_SHAPE_TSC,  0.0f, 0},
        {"TSC + N=4",  GR_SHAPE_TSC,  0.0f, 4},
        {"TSC + N=16", GR_SHAPE_TSC,  0.0f, 16},
        {"TSC + N=64", GR_SHAPE_TSC,  0.0f, 64},
        {"BUMP R=2.5", GR_SHAPE_BUMP, 2.5f, 0},
        {"BUMP R=4",   GR_SHAPE_BUMP, 4.0f, 0},
        {"BUMP R=6",   GR_SHAPE_BUMP, 6.0f, 0},
        {"BUMP R=8",   GR_SHAPE_BUMP, 8.0f, 0},
    };
    const int Nk = sizeof(kernels)/sizeof(kernels[0]);

    printf("  v/c   ");
    for (int ik = 0; ik < Nk; ik++) printf(" %-11s", kernels[ik].label);
    printf("\n  ");
    for (int ik = 0; ik < Nk; ik++) printf("------------");
    printf("\n");

    for (int iv = 0; iv < Nv; iv++) {
        printf("  %.3f", (double)vs[iv]);
        for (int ik = 0; ik < Nk; ik++) {
            const float F = measure_F(kernels[ik].shape, kernels[ik].Rk, kernels[ik].smooth, vs[iv]);
            printf("  %+10.3e", (double)F);
        }
        printf("\n");
    }
    return 0;
}
