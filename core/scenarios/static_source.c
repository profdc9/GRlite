/* Scenario "static_source" — single CIC-deposited point source for Stage 3
 * "Static CIC source, Poisson convergence" (gr_sandbox_v33.tex §12.3).
 *
 * Parameters (defaults applied when <= 0):
 *   params[0]: mass M           (default: 1.0)
 *   params[1]: x0               (default: width  * dx / 2)
 *   params[2]: y0               (default: height * dx / 2)
 *
 * Initial field is zero; the wave-equation source -4*pi*G_eff*rho at the
 * deposited cells drives Phi outward as the simulation steps. With the
 * absorbing layer engaged (Section §9.6) the transient radiates out and the
 * field converges to the discrete-Poisson static solution Lap Phi = 4 pi
 * G_eff rho_matter.
 */

#include "grlite.h"
#include "sim_internal.h"

static int build_static_source(gr_sim_t* sim, const float* params, int n_params) {
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx  = sim->dx;

    const float mass = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    const float x0   = (n_params >= 2) ? params[1] : ((float) W * 0.5f) * dx;
    const float y0   = (n_params >= 3) ? params[2] : ((float) H * 0.5f) * dx;

    /* Phi^{-1} = Phi^0 = 0 (gr_sim_create's calloc has already done this; reset
     * explicitly in case the scenario is loaded into a re-used sim). */
    const size_t n = (size_t) W * (size_t) H;
    for (size_t k = 0; k < n; k++) {
        sim->phi_prev[k] = 0.0f;
        sim->phi_curr[k] = 0.0f;
    }

    /* Clear any prior sources, then CIC-deposit the single point mass. */
    gr_sim_clear_sources(sim);
    gr_sim_deposit_point_mass(sim, x0, y0, mass);

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
