/* Scenario "eccentric_orbit" — single test particle on an eccentric orbit
 * around a softened-point-mass background, for the Stage 8 perihelion
 * precession test (gr_sandbox_v33.tex §12.8).
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: GM (central mass * G_eff)  (default 1.0)
 *   params[1]: semi-major axis a          (default 30.0)
 *   params[2]: eccentricity e             (default 0.3)
 *   params[3]: smoothing length eps       (default 1.0)
 *   params[4]: x0  (default width *dx / 2)
 *   params[5]: y0  (default height*dx / 2)
 *
 * Initial conditions: particle at periapsis (x0 + a(1-e), y0) with velocity
 * in the +y direction, magnitude from the Newtonian vis-viva relation
 *     v_peri^2 = GM * (1+e) / (a * (1-e)).
 * This defines a/e in the Kepler-limit sense: at v/c -> 0 the orbit is the
 * standard ellipse.  With Tier-2 relativistic corrections enabled, the orbit
 * precesses at rate dphi = 6 pi GM / (c^2 a (1-e^2)) per radial period
 * (Schwarzschild result, eq:dphi_prec). The softening correction is
 * negligible whenever eps << r_peri, which is the practical regime.
 *
 * Damping is NOT enabled — Stage 8 is just a particle in a fixed background;
 * no radiation, no waves. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_eccentric_orbit(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float GM    = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    const float a     = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 30.0f;
    const float e     = (n_params >= 3 && params[2] >= 0.0f && params[2] < 1.0f)
                            ? params[2] : 0.3f;
    const float eps   = (n_params >= 4 && params[3] > 0.0f) ? params[3] : 1.0f;
    const float x0    = (n_params >= 5) ? params[4] : ((float) W * 0.5f) * dx;
    const float y0    = (n_params >= 6) ? params[5] : ((float) H * 0.5f) * dx;

    /* Zero perturbation fields and sources. */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    gr_sim_clear_sources(sim);

    gr_sim_clear_background(sim);
    gr_sim_set_background_point_mass(sim, x0, y0, GM, eps);
    gr_sim_clear_particles(sim);

    /* No perturbation dynamics in Stage 8 — fields stay at zero. Skip the
     * wave-equation leapfrog to save compute. */
    gr_sim_set_field_evolution(sim, 0);
    /* Analytic-background evaluation: the precession test is sensitive to
     * tangential force error from the CIC+FD-of-sampled-bg chain, which
     * would otherwise contaminate the measurement. The doc explicitly
     * permits this fallback for Stage 8 (§12.8 caveat at line 2240). */
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

    /* Newtonian periapsis kinematics — for moderately eccentric, moderately
     * relativistic orbits the actual (a, e) deviate from these inputs by
     * O(v^2/c^2), and the precession formula treats (a, e) in the
     * weak-field-limit sense as the orbit-averaged geometry. */
    const float r_peri = a * (1.0f - e);
    /* vis-viva at periapsis: v^2 = GM (2/r - 1/a) = GM (1+e)/(a(1-e)) */
    const float v_peri = sqrtf(GM * (1.0f + e) / (a * (1.0f - e)));

    gr_sim_add_particle(sim, x0 + r_peri, y0,
                        /*mass=*/1.0f, /*charge=*/0.0f,
                        /*vx=*/0.0f,   /*vy=*/v_peri);
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_ECCENTRIC_ORBIT = {
    .name  = "eccentric_orbit",
    .build = build_eccentric_orbit,
};

void gr_scenario_register_eccentric_orbit(void) {
    gr_scenario_register(&SCENARIO_ECCENTRIC_ORBIT);
}
