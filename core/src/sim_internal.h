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
    float G_eff;        /* effective gravitational constant (Stage 3+); default 1.0 */
    int   step_count;

    float* phi_prev;   /* Phi^{n-1} */
    float* phi_curr;   /* Phi^n     */
    float* phi_next;   /* Phi^{n+1} (scratch) */

    /* Stage 3 static source — same staggering as phi_g (cell-centered).
     * Always allocated by gr_sim_create (zero-filled); deposited into by
     * scenarios via the CIC kernel gr_cic_deposit_scalar (eq:cic_deposit). */
    float* rho_matter;

    /* Damping layer (Stage 2, gr_sandbox_v32.tex §9.6).
     *   n_damping       : layer thickness on each edge (0 = no damping)
     *   damping_d       : precomputed d_{i,j} = max(d_x(i), d_y(j)) array (NULL if disabled)
     * Applied per step as Phi^{n+1}_{i,j} <- Phi^{n+1}_{i,j} * (1 - d_{i,j}) inside the
     * leapfrog kernel — see eq:damp_profile and the implementation block below it. */
    int    n_damping;
    float* damping_d;

    /* Background field arrays (Stage 6, gr_sandbox_v33.tex §12.6).
     * Lazily allocated by gr_sim_set_background_*; NULL until first use.
     * Read-only after scenario setup. Never touched by the leapfrog or the
     * damping layer — see field.c (leapfrog only modifies the perturbation
     * arrays phi_prev/curr/next). Force evaluation (future Stage 7+) will read
     * Phi_total = Phi_bg + Phi_pert via a single CIC interpolation. */
    float* phi_g_bg;
    float* Agx_bg;
    float* Agy_bg;
    float* phi_bg;
    float* Ax_bg;
    float* Ay_bg;
};

/* Defined in field.c — fills sim->phi_next from sim->phi_curr and sim->phi_prev. */
void gr_field_leapfrog_step(struct gr_sim* sim);

/* Defined in deposit.c — CIC deposition of a scalar value at sub-cell position.
 * Used by scenarios to set up rho_matter (Stage 3) and, later, particle
 * source deposition (Stage 9+). */
void gr_cic_deposit_scalar(float* rho, int W, int H, float dx,
                           float x_p, float y_p, float value);

/* Defined in scenarios/registry.c. */
const gr_scenario_t* gr_scenario_find(const char* name);

#endif /* GRLITE_SIM_INTERNAL_H */
