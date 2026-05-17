/* Internal simulation struct — visible only to .c files inside core/. */
#ifndef GRLITE_SIM_INTERNAL_H
#define GRLITE_SIM_INTERNAL_H

#include "grlite.h"

/* Three time levels of the scalar potential field, rotated after each leapfrog step.
 * Layout: row-major, index k = j * width + i, for i in [0, width), j in [0, height).
 *
 * The leapfrog (gr_sandbox_v32.tex §9.2 eq:leapfrog_field) needs Phi at the previous
 * and current timesteps to produce the next:
 *
 *     Phi^{n+1} = 2 Phi^n - Phi^{n-1} + (c dt)^2 (Lap Phi^n + source^n)
 *
 * After computing next, we rotate so the buffer that was "prev" becomes the
 * scratch for the following step's "next" (zero-copy three-pointer rotation). */
struct gr_sim {
    int   width, height;
    float dx;
    float c_eff;
    float dt;
    float cfl;
    int   step_count;

    float* phi_prev;   /* Phi^{n-1} */
    float* phi_curr;   /* Phi^n     */
    float* phi_next;   /* Phi^{n+1} (scratch) */

    /* Damping layer (Stage 2, gr_sandbox_v32.tex §9.6).
     *   n_damping       : layer thickness on each edge (0 = no damping)
     *   damping_d       : precomputed d_{i,j} = max(d_x(i), d_y(j)) array (NULL if disabled)
     * Applied per step as Phi^{n+1}_{i,j} <- Phi^{n+1}_{i,j} * (1 - d_{i,j}) inside the
     * leapfrog kernel — see eq:damp_profile and the implementation block below it. */
    int    n_damping;
    float* damping_d;
};

/* Defined in field.c — fills sim->phi_next from sim->phi_curr and sim->phi_prev. */
void gr_field_leapfrog_step(struct gr_sim* sim);

/* Defined in scenarios/registry.c. */
const gr_scenario_t* gr_scenario_find(const char* name);

#endif /* GRLITE_SIM_INTERNAL_H */
