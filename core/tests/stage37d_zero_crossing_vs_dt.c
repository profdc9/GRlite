/* Stage 37d -- does the F(v) zero-crossing depend on dt (CFL)?
 *
 * Dan, 2026-05-31: "Would the velocity for which heating to drag
 * change if the timestep was different?"
 *
 * The original Stage 37 v-scan with CFL = 1/sqrt(2) showed
 *   F(v=0.001) = +1.18e-6 (heating)
 *   F(v=0.010) = +1.68e-6 (heating)
 *   F(v=0.050) = -6.6e-7  (drag)
 * with zero crossing somewhere in (0.01, 0.05).
 *
 * Memory note from prior convergence work says the spurious-force
 * MAGNITUDE is dt-independent but dx-dependent.  Here we focus on the
 * specific question of whether the ZERO CROSSING velocity v_zero shifts
 * with dt.  If v_zero is purely a spatial-grid property, halving dt
 * (smaller CFL) should leave v_zero unchanged.  If v_zero shifts with
 * dt, that points to a temporal-discretization piece.
 *
 * Method: measure F(v) at several v straddling the expected zero
 * crossing, at three CFL values (1/sqrt(2), 1/2, 1/4).  Equal-distance
 * 2-cell measurement window after 8-cell warmup, mass = 1.0 so v
 * barely changes during measurement.  BUMP R=8 production kernel
 * plus TSC bare for reference. */

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

static float measure_F(gr_shape_function_t shape, float Rk, int smooth,
                        float v_drift, float cfl) {
    const int W = 128, H = 1024, n_damp = 16;
    const float dx = 1.0f, c_eff = 1.0f;
    const float mass = 1.0f, Q = 0.01f;
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
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);

    const float cx = ((float)(W-1) * 0.5f) * dx;
    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);
    const float vy_w = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    for (int s = 0; s < n_measure; s++) gr_sim_step(sim);
    const float vy_m = vy_of(gr_sim_get_particle(sim, 0), c_eff);
    const float F = mass * (vy_m - vy_w) / ((float)n_measure * dt);
    gr_sim_destroy(sim);
    return F;
}

static void scan(const char* label, gr_shape_function_t shape, float Rk, int smooth) {
    const float vs[] = {0.005f, 0.010f, 0.020f, 0.030f, 0.040f, 0.050f, 0.075f, 0.100f};
    const int Nv = sizeof(vs)/sizeof(vs[0]);
    const float cfls[]  = {1.0f/sqrtf(2.0f), 0.5f, 0.25f};
    const char* lab[]   = {"CFL=1/sqrt2", "CFL=1/2    ", "CFL=1/4    "};
    const int Nc = 3;
    printf("=== %s ===\n", label);
    printf("  v/c         CFL=1/sqrt2  CFL=1/2      CFL=1/4    \n");
    printf("  -----------------------------------------------\n");
    for (int iv = 0; iv < Nv; iv++) {
        printf("  v=%.3f", (double)vs[iv]);
        for (int ic = 0; ic < Nc; ic++) {
            const float F = measure_F(shape, Rk, smooth, vs[iv], cfls[ic]);
            (void)lab;
            printf("   %+11.4e", (double)F);
        }
        printf("\n");
    }
    printf("\n");
}

int main(void) {
    printf("=== stage37d_zero_crossing_vs_dt ===\n");
    printf("Does F_spurious(v) zero-crossing depend on dt (CFL)?\n");
    printf("Equal-distance measurement (8-cell warmup + 2-cell window) at mass=1.0.\n\n");
    scan("TSC bare", GR_SHAPE_TSC,  0.0f, 0);
    scan("BUMP R=8", GR_SHAPE_BUMP, 8.0f, 0);
    return 0;
}
