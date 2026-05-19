/* Scenario "pic_constant_v" — single particle moving at constant velocity
 * through empty space, with the perturbation FDTD running and auto-
 * deposition enabled.  Stage 10's boosted-Coulomb field-shape test uses
 * this.
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: mass m  (default 1.0)
 *   params[1]: vx      (default 0.05)
 *   params[2]: vy      (default 0.0)
 *   params[3]: x0      (default width *dx/2)
 *   params[4]: y0      (default height*dx/2)
 *
 * No background is installed; the particle is initialized at (x0,y0) with
 * velocity (vx,vy).  Use the damping layer if you want to suppress
 * back-reflection from the box walls while measuring fields. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_pic_constant_v(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float mass = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    const float vx   = (n_params >= 2) ? params[1] : 0.05f;
    const float vy   = (n_params >= 3) ? params[2] : 0.0f;
    const float x0   = (n_params >= 4) ? params[3] : ((float) W * 0.5f) * dx;
    const float y0   = (n_params >= 5) ? params[4] : ((float) H * 0.5f) * dx;

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

    gr_sim_set_bg_mode(sim, GR_BG_MODE_SAMPLED);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);

    gr_sim_add_particle(sim, x0, y0,
                        mass, /*charge=*/0.0f,
                        vx, vy);
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_PIC_CONSTANT_V = {
    .name  = "pic_constant_v",
    .build = build_pic_constant_v,
};

void gr_scenario_register_pic_constant_v(void) {
    gr_scenario_register(&SCENARIO_PIC_CONSTANT_V);
}
