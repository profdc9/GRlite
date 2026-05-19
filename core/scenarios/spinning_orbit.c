/* Scenario "spinning_orbit" — single test particle on a circular orbit
 * around a softened SPINNING point mass.  Used by Stage 9's gravitomagnetic
 * clock-effect test (gr_sandbox_v34.tex §12.9).
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: GM (central mass * G_eff)             (default 1.0)
 *   params[1]: orbital radius r                      (default 20.0)
 *   params[2]: smoothing length eps                  (default 1.0)
 *   params[3]: J_z (central spin along +z)           (default 0.0)
 *   params[4]: direction (+1 prograde, -1 retrograde)(default +1.0)
 *   params[5]: x0                                    (default width *dx/2)
 *   params[6]: y0                                    (default height*dx/2)
 *
 * Initial conditions: same relativistic v_circ derivation as kepler_orbit
 * (Newtonian-softened gradient + relativistic momentum); the gravitomagnetic
 * force shift in v_circ is a 1PN correction to the orbit shape that doesn't
 * affect the leading-order clock-effect measurement.  The particle is placed
 * at (x0 + r, y0) with velocity (0, +/- v_circ) according to the direction
 * parameter.
 *
 * Damping is NOT enabled; the perturbation FDTD is skipped via
 * gr_sim_set_field_evolution(sim, 0); the background is evaluated
 * analytically (GR_BG_MODE_ANALYTIC) so that the clock-effect signal isn't
 * contaminated by the CIC+FD tangential force error. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_spinning_orbit(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float GM    = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    const float r     = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 20.0f;
    const float eps   = (n_params >= 3 && params[2] > 0.0f) ? params[2] : 1.0f;
    const float Jz    = (n_params >= 4) ? params[3] : 0.0f;
    float       sign  = (n_params >= 5) ? params[4] : 1.0f;
    /* Snap direction to ±1. */
    sign = (sign >= 0.0f) ? 1.0f : -1.0f;
    const float x0    = (n_params >= 6) ? params[5] : ((float) W * 0.5f) * dx;
    const float y0    = (n_params >= 7) ? params[6] : ((float) H * 0.5f) * dx;

    /* Zero perturbation fields and sources. */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    gr_sim_clear_sources(sim);

    gr_sim_clear_background(sim);
    gr_sim_set_background_spinning_point_mass(sim, x0, y0, GM, eps, Jz);
    gr_sim_clear_particles(sim);

    /* Stage 9 setup: static-background-only, evaluate analytic. */
    gr_sim_set_field_evolution(sim, 0);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);

    /* Relativistic circular-orbit velocity from the softened (non-spinning)
     * force law, identical to kepler_orbit.c.  The spin-induced 1PN shift in
     * v_circ is small at our parameters and doesn't affect the leading
     * gravitomagnetic-clock measurement. */
    const float r2_e2_15 = powf(r * r + eps * eps, 1.5f);
    const float g_mag    = GM * r / r2_e2_15;
    const float rg       = r * g_mag;
    const float rg2      = rg * rg;
    const float c2       = sim->c_eff * sim->c_eff;
    const float rg2_c2   = rg2 / c2;
    const float u        = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg2) - rg2_c2) * 0.5f;
    const float v_circ   = sqrtf(u);

    /* Particle at (x0 + r, y0) with tangential v = (0, +/- v_circ).
     * Prograde (sign = +1) co-rotates with J_z > 0; retrograde counter-rotates. */
    gr_sim_add_particle(sim, x0 + r, y0,
                        /*mass=*/1.0f, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/sign * v_circ);
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_SPINNING_ORBIT = {
    .name  = "spinning_orbit",
    .build = build_spinning_orbit,
};

void gr_scenario_register_spinning_orbit(void) {
    gr_scenario_register(&SCENARIO_SPINNING_ORBIT);
}
