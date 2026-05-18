/* Scenario "static_source" — CIC-deposited point mass and/or point charge at
 * the domain center for Stage 3 "Static CIC source, Poisson convergence"
 * (§12.3) and Stage 4 "All six wave equations, gauge monitoring" (§12.4).
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: mass M    (default 1.0)
 *   params[1]: charge Q  (default 0.0; deposits nothing if zero)
 *   params[2]: x0        (default width  * dx / 2)
 *   params[3]: y0        (default height * dx / 2)
 *
 * Both deposits are at zero velocity, so J_matter = J_q = 0 and the four
 * vector potentials (A_g_{x,y}, A_{x,y}) have no source. They start at zero
 * (from gr_sim_create's calloc) and stay exactly zero through every leapfrog
 * step — verified by the Stage 4 test.
 */

#include "grlite.h"
#include "sim_internal.h"

#include <string.h>

static int build_static_source(gr_sim_t* sim, const float* params, int n_params) {
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx  = sim->dx;

    const float mass   = (n_params >= 1) ? params[0] : 1.0f;
    const float charge = (n_params >= 2) ? params[1] : 0.0f;
    const float x0     = (n_params >= 3) ? params[2] : ((float) W * 0.5f) * dx;
    const float y0     = (n_params >= 4) ? params[3] : ((float) H * 0.5f) * dx;

    /* All six fields' prev/curr to zero (§9.7 zero-derivative initialization). */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }

    gr_sim_clear_sources(sim);
    if (mass   != 0.0f) gr_sim_deposit_point_mass(sim, x0, y0, mass);
    if (charge != 0.0f) gr_sim_deposit_point_charge(sim, x0, y0, charge);

    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_STATIC_SOURCE = {
    .name  = "static_source",
    .build = build_static_source,
};

void gr_scenario_register_static_source(void) {
    gr_scenario_register(&SCENARIO_STATIC_SOURCE);
}
