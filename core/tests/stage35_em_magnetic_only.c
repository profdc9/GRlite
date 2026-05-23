/* Stage 35 -- DIAGNOSTIC: pure-magnetic force between parallel like-currents.
 *
 * Two like-sign +Q particles separated in x and drifting in +y (so
 * each has a +y current).  The full EM Lorentz force has three pieces:
 *   F_E    = -q grad phi          (electrostatic)
 *   F_ind  = -q d_t A             (inductive)
 *   F_B    = q v x B = q v x (grad x A)  (magnetic)
 *
 * Here we GATE OFF F_E and F_ind, leaving ONLY F_B active in the
 * particle pusher.  The fields phi_em and A_em are still updated
 * normally by the wave equation -- the deposit, leapfrog, and rotation
 * all proceed as usual.  Only the FORCE EVALUATION on the particle
 * uses just q v x B.
 *
 * Question to answer: do the two parallel like-currents move TOWARD
 * each other (Ampere attraction, standard Maxwell) or AWAY from each
 * other (anti-Ampere)?
 *
 * This isolates the sign of the magnetic force from the much-larger
 * electrostatic repulsion and from the inductive piece's
 * radiation-reaction contribution. */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); \
        return 1; \
    } \
} while (0)

static float run_pair(float v_drift, int n_steps, int* viols_out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float Q      = 0.01f;
    const float mass   = 0.01f;
    const float d_sep  = 16.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return NAN;
    gr_sim_set_damping(sim, 16);
    gr_sim_set_field_evolution(sim, 1);                /* fields update normally */
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);

    /* Force law: ONLY v x B active.  Disable -grad phi and -d_t A. */
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 0);
    gr_sim_set_em_inductive_enabled(sim, 0);
    /* Disable GRAVITY entirely so the gravitational attraction between
     * the two massive particles doesn't dominate the small v x B signal. */
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

    gr_sim_add_particle(sim, cx - d_sep * 0.5f, cy, mass, +Q, 0.0f, v_drift);
    gr_sim_add_particle(sim, cx + d_sep * 0.5f, cy, mass, +Q, 0.0f, v_drift);

    for (int s = 0; s < n_steps; s++) gr_sim_step(sim);

    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
    const float r_final = fabsf(p1->x - p0->x);
    if (viols_out) *viols_out = gr_sim_esirkepov_violations(sim);
    gr_sim_destroy(sim);
    return r_final - d_sep;   /* >0: separated; <0: attracted */
}

int main(void) {
    printf("=== stage35_em_magnetic_only ===\n");
    printf("Two like +Q charges drifting in +y, ONLY v x B force active.\n");
    printf("(Fields phi_em and A_em update normally via wave equation.)\n\n");

    const float Q      = 0.01f;
    const float mass   = 0.01f;
    const float d_sep  = 16.0f;
    const float k_e    = 1.0f;
    (void) Q; (void) mass; (void) d_sep; (void) k_e;

    /* Sweep v_drift.  v=0 should give zero magnetic force (no current).
     * v>0 should give a magnetic force whose sign tells us if Ampere
     * (attractive) or anti-Ampere. */
    const float vs[]      = { 0.00f, 0.10f, 0.20f, 0.30f, 0.40f };
    const int   n_v       = (int)(sizeof(vs)/sizeof(vs[0]));
    const int   n_steps   = 200;

    printf("%-8s %-14s %-7s\n", "v_drift", "widening", "viols");
    printf("---------------------------------------\n");
    for (int i = 0; i < n_v; i++) {
        int viols = 0;
        const float w = run_pair(vs[i], n_steps, &viols);
        const char* dir;
        if (fabsf(w) < 1e-4f) dir = "(zero)";
        else if (w > 0.0f)    dir = "REPULSIVE";
        else                  dir = "ATTRACTIVE";
        printf("%-8.3f %+11.5e   %-7d  %s\n", (double) vs[i], (double) w, viols, dir);
    }

    printf("\n");
    printf("Interpretation:\n");
    printf("  v=0 baseline: widening should be ~0 (no current -> no magnetic force).\n");
    printf("  v>0: sign of widening tells us the v x B direction:\n");
    printf("    - NEGATIVE widening = attractive = standard Ampere (parallel\n");
    printf("      like-currents attract).\n");
    printf("    - POSITIVE widening = repulsive = anti-Ampere (suspicious).\n");
    printf("  Magnitude should scale as ~(v/c)^2.\n");

    /* Stage 35 is purely diagnostic: no assertions on direction. */
    printf("\nDIAGNOSTIC COMPLETE (no pass/fail on direction).\n");
    return 0;
}
