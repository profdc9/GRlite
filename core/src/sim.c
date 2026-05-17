/* Simulation lifecycle and scenario dispatch.
 * Spec reference: gr_sandbox_v32.tex §9 "Numerical Implementation". */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

gr_sim_t* gr_sim_create(int width, int height, float dx, float c_eff, float cfl) {
    if (width <= 2 || height <= 2 || !(dx > 0.0f) || !(c_eff > 0.0f) || !(cfl > 0.0f)) {
        return NULL;
    }
    gr_sim_t* sim = (gr_sim_t*) calloc(1, sizeof(*sim));
    if (!sim) return NULL;

    sim->width  = width;
    sim->height = height;
    sim->dx     = dx;
    sim->c_eff  = c_eff;
    sim->cfl    = cfl;
    /* dt from CFL — gr_sandbox_v32.tex §9.2: dt = cfl * dx / c_eff, with stability requiring
     * cfl <= 1/sqrt(d) in d spatial dimensions. We do not enforce here so callers can
     * deliberately exceed the limit for the Stage 1 instability test (§12.1). */
    sim->dt = cfl * dx / c_eff;

    const size_t n = (size_t) width * (size_t) height;
    sim->phi_prev = (float*) calloc(n, sizeof(float));
    sim->phi_curr = (float*) calloc(n, sizeof(float));
    sim->phi_next = (float*) calloc(n, sizeof(float));
    if (!sim->phi_prev || !sim->phi_curr || !sim->phi_next) {
        gr_sim_destroy(sim);
        return NULL;
    }
    return sim;
}

void gr_sim_destroy(gr_sim_t* sim) {
    if (!sim) return;
    free(sim->phi_prev);
    free(sim->phi_curr);
    free(sim->phi_next);
    free(sim);
}

void gr_sim_step(gr_sim_t* sim) {
    if (!sim) return;
    gr_field_leapfrog_step(sim);
    /* Three-pointer rotation: the old prev buffer becomes the next scratch.
     * gr_sandbox_v32.tex §9.2 — only two time levels need be retained between steps. */
    float* tmp     = sim->phi_prev;
    sim->phi_prev  = sim->phi_curr;
    sim->phi_curr  = sim->phi_next;
    sim->phi_next  = tmp;
    sim->step_count++;
}

void gr_sim_step_n(gr_sim_t* sim, int n) {
    if (!sim || n <= 0) return;
    for (int i = 0; i < n; i++) gr_sim_step(sim);
}

int   gr_sim_step_count(const gr_sim_t* sim) { return sim ? sim->step_count : 0; }
float gr_sim_time(const gr_sim_t* sim)       { return sim ? sim->dt * (float) sim->step_count : 0.0f; }
float gr_sim_dt(const gr_sim_t* sim)         { return sim ? sim->dt : 0.0f; }
float gr_sim_dx(const gr_sim_t* sim)         { return sim ? sim->dx : 0.0f; }
int   gr_sim_width(const gr_sim_t* sim)      { return sim ? sim->width : 0; }
int   gr_sim_height(const gr_sim_t* sim)     { return sim ? sim->height : 0; }

float* gr_sim_field_ptr(gr_sim_t* sim, gr_field_id_t which) {
    if (!sim) return NULL;
    if (which == GR_FIELD_PHI_GRAV) return sim->phi_curr;
    return NULL;
}

int gr_sim_load_scenario(gr_sim_t* sim, const char* name, const float* params, int n_params) {
    if (!sim || !name) return -1;
    gr_scenarios_init();
    const gr_scenario_t* s = gr_scenario_find(name);
    if (!s) return -2;
    return s->build(sim, params, n_params);
}
