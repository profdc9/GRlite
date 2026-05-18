/* Scenario "kepler_orbit" — single test particle on a circular orbit around a
 * softened-point-mass background (gr_sandbox_v33.tex §12.7 "Single test
 * particle, Boris pusher, Keplerian orbit").
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: GM (central mass * G_eff)  (default 1.0)
 *   params[1]: orbital radius r           (default 20.0)
 *   params[2]: smoothing length eps       (default 1.0)
 *   params[3]: x0  (default width *dx / 2)
 *   params[4]: y0  (default height*dx / 2)
 *
 * Initial conditions: particle at (x0 + r, y0) with velocity in the +y
 * direction, magnitude chosen from the SOFTENED force law
 *   |g| = GM * r / (r^2 + eps^2)^{3/2}
 *   v   = sqrt(r * |g|)  (centripetal balance)
 * so the orbit is circular for the softened potential, not the Newtonian
 * 1/r^2 force law. At r/eps >> 1 these converge.
 *
 * Damping is NOT enabled here — the particle is in vacuum, and no waves
 * radiate (Stage 7 has no perturbation sources). */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_kepler_orbit(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float GM    = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    const float r     = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 20.0f;
    const float eps   = (n_params >= 3 && params[2] > 0.0f) ? params[2] : 1.0f;
    const float x0    = (n_params >= 4) ? params[3] : ((float) W * 0.5f) * dx;
    const float y0    = (n_params >= 5) ? params[4] : ((float) H * 0.5f) * dx;

    /* Zero out all perturbation fields and sources. */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    gr_sim_clear_sources(sim);

    /* Install the softened point mass as background, then clear any particles
     * from a prior scenario and add the single orbiter. */
    gr_sim_clear_background(sim);
    gr_sim_set_background_point_mass(sim, x0, y0, GM, eps);
    gr_sim_clear_particles(sim);

    /* Relativistic circular-orbit velocity for the softened force law.
     *
     * The Boris pusher is fully relativistic, so the centripetal condition is
     *     gamma * v^2 / r = g(r),    gamma = 1 / sqrt(1 - v^2/c^2)
     * not the Newtonian v^2 / r = g(r). At v/c ~ 0.2 the relativistic v_circ
     * is ~1.3% below the Newtonian value — initializing with the Newtonian
     * formula puts the particle on an eccentric orbit with apoapsis above r.
     *
     * Closed form (solve for u = v^2):
     *     u = ( sqrt( (rg)^4 / c^4 + 4 (rg)^2 ) - (rg)^2 / c^2 ) / 2
     * which reduces to u = rg in the non-relativistic limit c -> infinity. */
    const float r2_e2_15 = powf(r * r + eps * eps, 1.5f);
    const float g_mag    = GM * r / r2_e2_15;
    const float rg       = r * g_mag;
    const float rg2      = rg * rg;
    const float c2       = sim->c_eff * sim->c_eff;
    const float rg2_c2   = rg2 / c2;
    const float u        = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg2) - rg2_c2) * 0.5f;
    const float v_circ   = sqrtf(u);
    gr_sim_add_particle(sim, x0 + r, y0, /*mass=*/1.0f, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/v_circ);
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_KEPLER_ORBIT = {
    .name  = "kepler_orbit",
    .build = build_kepler_orbit,
};

void gr_scenario_register_kepler_orbit(void) {
    gr_scenario_register(&SCENARIO_KEPLER_ORBIT);
}
