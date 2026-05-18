/* Scenario "moving_source" — point particle with mass, charge, and (held-fixed)
 * velocity deposited at the domain center. Used by the Stage 5 test (§12.5):
 * verifies that nonzero J sources drive the vector potentials toward the
 * boosted-Coulomb relation A ~ v phi / c^2.
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: mass M      (default 1.0)
 *   params[1]: charge Q    (default 0.0)
 *   params[2]: vx          (default 0.1 * c_eff)
 *   params[3]: vy          (default 0.0)
 *   params[4]: x0          (default center)
 *   params[5]: y0          (default center)
 *
 * The velocity here is a SCENARIO-TIME constant — the particle does not
 * physically move during stepping. This is the §9.7 "convergence iteration"
 * usage: position and velocity held fixed so the wave equations settle to a
 * static field that *would* exist around the moving particle.
 *
 * Holding rho static while J = rho * v is nonzero deliberately violates the
 * continuum continuity equation d_t rho + div J = 0 by exactly v . grad rho.
 * Stage 5's test measures this violation directly so we can track its effect
 * on gauge drift over time.
 */

#include "grlite.h"
#include "sim_internal.h"

#include <string.h>

static int build_moving_source(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float mass   = (n_params >= 1) ? params[0] : 1.0f;
    const float charge = (n_params >= 2) ? params[1] : 0.0f;
    const float vx     = (n_params >= 3) ? params[2] : 0.1f * sim->c_eff;
    const float vy     = (n_params >= 4) ? params[3] : 0.0f;
    const float x0     = (n_params >= 5) ? params[4] : ((float) W * 0.5f) * dx;
    const float y0     = (n_params >= 6) ? params[5] : ((float) H * 0.5f) * dx;

    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }

    gr_sim_clear_sources(sim);
    gr_sim_deposit_point_particle(sim, x0, y0, mass, charge, vx, vy);

    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_MOVING_SOURCE = {
    .name  = "moving_source",
    .build = build_moving_source,
};

void gr_scenario_register_moving_source(void) {
    gr_scenario_register(&SCENARIO_MOVING_SOURCE);
}
