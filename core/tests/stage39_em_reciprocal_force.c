/* Stage 39 -- BENCHMARK B-3: Reciprocal force action-reaction (gr_sandbox_v36 sec gempic_benchmark).
 *
 * Two like charges (+Q, +Q), both at rest at $(c_x, c_y)$ and
 * $(c_x, c_y + r)$.  Static configuration.  Measure the force on each
 * particle by reading the momentum change per step.  Verify
 *      |F_1 + F_2| / |F_1| < 1e-6
 * at all separations r >= 4 cells.
 *
 * Physics: continuous Newton's third law.  Each particle sees the
 * other's static Coulomb field, plus its own self-force (which is
 * zero by Hockney-Eastwood at v = 0 with matched deposit/interp).
 *
 * Why this probe (per v36 sec gempic_benchmark): GUARDRAIL test.  The
 * current scheme already passes static reciprocity at v=0; the point
 * is to confirm GEMPIC's modified force law does not break it
 * accidentally.  A "fix" that solves Stage 37 but breaks Stage 39 is
 * a regression.
 *
 * Method: use a HEAVY particle mass (1e6) so the particles do not
 * appreciably move during the test.  Charge is q = 0.01, same as the
 * other benchmark stages.  Run warmup so phi reaches its static
 * Coulomb pattern, then read each particle's dp / dt over a few steps.
 *
 * Acceptance (per v36):
 *   |F_1 + F_2| / |F_1| < 1e-6 at all separations r >= 4 cells. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float Fx_1, Fy_1;
    float Fx_2, Fy_2;
    float Fx_sum, Fy_sum;
    float F1_mag, sum_mag;
    float reciprocity_ratio;  /* |F_1 + F_2| / |F_1| */
    int   sep_cells;
    int   bailed_out;
} reciprocity_result_t;

static reciprocity_result_t run_reciprocity(int sep_cells, int n_warmup, int n_force_avg) {
    const int   W       = 64;
    const int   H       = 64;
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 8;
    const float dt      = cfl * dx / c_eff;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;
    const float cy_base = ((float) (H - 1) * 0.5f) * dx;
    /* Particles separated along y, symmetric about cy_base. */
    const float y1 = cy_base - 0.5f * (float) sep_cells * dx;
    const float y2 = cy_base + 0.5f * (float) sep_cells * dx;
    const float Q  = 0.01f;
    /* Very heavy mass so the particle does not move during the test. */
    const float mass = 1.0e6f;

    reciprocity_result_t r = {0};
    r.sep_cells = sep_cells;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { r.bailed_out = 1; return r; }
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);

    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, +1.0f);

    gr_sim_add_particle(sim, cx, y1, mass, +Q, 0.0f, 0.0f);
    gr_sim_add_particle(sim, cx, y2, mass, +Q, 0.0f, 0.0f);

    /* Warmup to let static phi develop. */
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    /* Sample F = dp/dt over n_force_avg steps after warmup.  Average for noise. */
    const gr_particle_t* pp1 = gr_sim_get_particle(sim, 0);
    const gr_particle_t* pp2 = gr_sim_get_particle(sim, 1);
    const float px1_pre = pp1->px;
    const float py1_pre = pp1->py;
    const float px2_pre = pp2->px;
    const float py2_pre = pp2->py;
    for (int s = 0; s < n_force_avg; s++) gr_sim_step(sim);
    pp1 = gr_sim_get_particle(sim, 0);
    pp2 = gr_sim_get_particle(sim, 1);
    const float Fx1 = (pp1->px - px1_pre) / ((float) n_force_avg * dt);
    const float Fy1 = (pp1->py - py1_pre) / ((float) n_force_avg * dt);
    const float Fx2 = (pp2->px - px2_pre) / ((float) n_force_avg * dt);
    const float Fy2 = (pp2->py - py2_pre) / ((float) n_force_avg * dt);

    r.Fx_1 = Fx1; r.Fy_1 = Fy1;
    r.Fx_2 = Fx2; r.Fy_2 = Fy2;
    r.Fx_sum = Fx1 + Fx2;
    r.Fy_sum = Fy1 + Fy2;
    r.F1_mag  = sqrtf(Fx1*Fx1 + Fy1*Fy1);
    r.sum_mag = sqrtf(r.Fx_sum*r.Fx_sum + r.Fy_sum*r.Fy_sum);
    r.reciprocity_ratio = (r.F1_mag > 0.0f) ? r.sum_mag / r.F1_mag : 0.0f;

    gr_sim_destroy(sim);
    return r;
}

int main(void) {
    printf("=== stage39_em_reciprocal_force (BENCHMARK B-3) ===\n");
    printf("Two like charges at rest, separated by r cells along y.  Heavy mass\n");
    printf("(1e6) so they do not move.  Read F = dp/dt on each particle after\n");
    printf("static phi has developed.  Newton's third law: F_1 + F_2 = 0.\n\n");
    printf("Setup: 64x64 grid, BARE deposit (rho_smooth=0, j_smooth=0),\n");
    printf("mass=1e6, Q=0.01, TSC+LB.  Warmup 500 steps, force avg over 100 steps.\n\n");

    const int n_warmup   = 500;
    const int n_force_avg = 100;
    const int sep_vals[] = { 4, 6, 8, 10, 16, 24 };
    const int nsep       = (int)(sizeof(sep_vals) / sizeof(sep_vals[0]));

    printf("%-6s %-13s %-13s %-13s %-13s %-13s %-13s\n",
           "sep", "F1_y", "F2_y", "F1+F2_y", "|F1|", "|F1+F2|", "ratio");
    printf("--------------------------------------------------------------------------------------\n");
    int worst_violation_sep = sep_vals[0];
    double worst_violation_ratio = 0.0;
    for (int i = 0; i < nsep; i++) {
        reciprocity_result_t r = run_reciprocity(sep_vals[i], n_warmup, n_force_avg);
        printf("%-6d %+12.4e %+12.4e %+12.4e %12.4e %12.4e %12.4e%s\n",
               sep_vals[i],
               (double) r.Fy_1, (double) r.Fy_2,
               (double) r.Fy_sum,
               (double) r.F1_mag, (double) r.sum_mag,
               (double) r.reciprocity_ratio,
               r.bailed_out ? " [BAIL]" : "");
        if ((double) r.reciprocity_ratio > worst_violation_ratio) {
            worst_violation_ratio = (double) r.reciprocity_ratio;
            worst_violation_sep = sep_vals[i];
        }
    }

    printf("\nInterpretation:\n");
    printf("  F1_y, F2_y: y-component of force on each particle.  Like charges\n");
    printf("    repel; with particle 1 below and particle 2 above on the y-axis,\n");
    printf("    F1_y should be NEGATIVE and F2_y POSITIVE (opposite signs).\n");
    printf("  F1+F2: vector sum.  Continuum third law => 0.  Should be NUMERICAL\n");
    printf("    NOISE, not physical signal.\n");
    printf("  ratio = |F1+F2| / |F1|: dimensionless reciprocity violation.\n\n");
    printf("Acceptance (per v36 sec gempic_benchmark B-3):\n");
    printf("  ratio < 1e-6 at every sep >= 4.\n");
    printf("  Worst observed ratio: %.4e at sep=%d.\n",
           worst_violation_ratio, worst_violation_sep);
    if (worst_violation_ratio < 1.0e-6) {
        printf("  *** PASSES *** the current scheme satisfies static reciprocity.\n");
        printf("  Guardrail set: GEMPIC must not regress this.\n");
    } else {
        printf("  *** FAILS *** at the 1e-6 threshold.  This is unexpected --\n");
        printf("  Hockney-Eastwood adjoint condition should give exact zero at\n");
        printf("  v=0.  Investigate before treating as a guardrail.\n");
    }

    printf("\nBASELINE COMPLETE.\n");
    return 0;
}
