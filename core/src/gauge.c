/* Lorenz gauge residual monitoring.
 * Spec reference: gr_sandbox_v32.tex §9.5 eq:lorenz_gem, eq:lorenz_em,
 * and the "Monitoring" subsection. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>

/* Eq. (eq:lorenz_gem) / eq:lorenz_em — discrete Lorenz gauge residual.
 *
 *   G = (1/c^2) d_t Phi_scalar + div A_vector
 *
 * In our cell-centered convention with rotation order applied AFTER each
 * step (so .prev holds Phi^n, .curr holds Phi^{n+1}), the one-sided backward
 * difference (curr - prev)/dt approximates d_t Phi at the half-step n+1/2.
 * The divergence of A uses a centered finite difference on .curr arrays
 * (which are at time n+1 after rotation). RMS is taken over the interior,
 * excluding the damping layer where the modified equation perturbs the
 * gauge condition. */
static float gauge_residual_rms(const gr_sim_t* sim,
                                gr_field_id_t scalar_id,
                                gr_field_id_t Ax_id,
                                gr_field_id_t Ay_id) {
    if (!sim) return 0.0f;
    const int   W      = sim->width;
    const int   H      = sim->height;
    const float dx     = sim->dx;
    const float dt     = sim->dt;
    const float c2     = sim->c_eff * sim->c_eff;
    const int   skip   = sim->n_damping > 0 ? sim->n_damping + 1 : 1;

    const float* phi_prev = sim->fields[scalar_id].prev;
    const float* phi_curr = sim->fields[scalar_id].curr;
    const float* Ax_curr  = sim->fields[Ax_id].curr;
    const float* Ay_curr  = sim->fields[Ay_id].curr;

    const float inv_c2dt = 1.0f / (c2 * dt);
    const float inv_2dx  = 1.0f / (2.0f * dx);

    double sum_sq = 0.0;
    long   count  = 0;
    for (int j = skip; j < H - skip; j++) {
        const int row = j * W;
        for (int i = skip; i < W - skip; i++) {
            const int   k     = row + i;
            const float dphi  = (phi_curr[k] - phi_prev[k]) * inv_c2dt;
            const float divA  = (Ax_curr[k + 1] - Ax_curr[k - 1]) * inv_2dx
                              + (Ay_curr[k + W] - Ay_curr[k - W]) * inv_2dx;
            const float g     = dphi + divA;
            sum_sq += (double) (g * g);
            count++;
        }
    }
    if (count == 0) return 0.0f;
    return (float) sqrt(sum_sq / (double) count);
}

float gr_sim_gauge_residual_grav(const gr_sim_t* sim) {
    return gauge_residual_rms(sim, GR_FIELD_PHI_GRAV, GR_FIELD_A_GX, GR_FIELD_A_GY);
}

float gr_sim_gauge_residual_em(const gr_sim_t* sim) {
    return gauge_residual_rms(sim, GR_FIELD_PHI_EM, GR_FIELD_A_X, GR_FIELD_A_Y);
}
