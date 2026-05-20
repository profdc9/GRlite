/* Scenario "pic_static" — single stationary particle in empty space, with
 * the perturbation FDTD running and auto-deposition enabled.  Stage 10
 * static-field-shape and self-force-cancellation tests use this setup.
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: mass m  (default 1.0)
 *   params[1]: x0      (default width *dx/2)
 *   params[2]: y0      (default height*dx/2)
 *
 * No background is installed (vacuum); the perturbation field is initially
 * zero and is sourced solely by the particle's deposited mass. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_pic_static(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float mass = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    /* Default to the TRUE box center ((W-1)/2, (H-1)/2) so the absorbing
     * damping layer is exactly symmetric around the particle.  Using
     * W*0.5*dx as the center (which puts the particle at the integer
     * corner W/2 for even W) puts it 1 cell closer to the right/top
     * damping than to the left/bottom, breaking the symmetry the HE
     * argument relies on.  At W=128, 0.82-cell static drift over 2000
     * steps becomes EXACTLY 0 when this is fixed.  Caller can still
     * override by passing params[1], params[2]. */
    const float x0   = (n_params >= 2) ? params[1] : ((float) (W - 1) * 0.5f) * dx;
    const float y0   = (n_params >= 3) ? params[2] : ((float) (H - 1) * 0.5f) * dx;

    /* Zero perturbation fields. */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    gr_sim_clear_sources(sim);
    gr_sim_clear_background(sim);
    gr_sim_clear_particles(sim);

    /* Vacuum: no analytic background, so any nonzero force on the particle
     * must come from its own perturbation field (the case we're testing
     * against zero in Test C). */
    gr_sim_set_bg_mode(sim, GR_BG_MODE_SAMPLED);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);

    /* Stationary particle. */
    gr_sim_add_particle(sim, x0, y0,
                        mass, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/0.0f);
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_PIC_STATIC = {
    .name  = "pic_static",
    .build = build_pic_static,
};

void gr_scenario_register_pic_static(void) {
    gr_scenario_register(&SCENARIO_PIC_STATIC);
}
