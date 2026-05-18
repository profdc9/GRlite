/* Background field array management and generators.
 * Spec reference: gr_sandbox_v33.tex §12.6 (sec:stage_bg)
 *                 "Stage 6: Sampled background field arrays". */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>

/* Lookup the storage slot for a given field id. Returns a pointer-to-pointer
 * so callers can both read and (re)assign the array. NULL for unknown ids. */
static float** background_slot(gr_sim_t* sim, gr_field_id_t which) {
    if (!sim) return NULL;
    switch (which) {
    case GR_FIELD_PHI_GRAV: return &sim->phi_g_bg;
    /* Future field IDs (vector potentials, EM scalar) plug in here as those
     * generators are added in later stages — Stage 12 onward per v33 §12. */
    default: return NULL;
    }
}

float* gr_sim_background_ptr(gr_sim_t* sim, gr_field_id_t which) {
    float** slot = background_slot(sim, which);
    return slot ? *slot : NULL;
}

void gr_sim_clear_background(gr_sim_t* sim) {
    if (!sim) return;
    free(sim->phi_g_bg); sim->phi_g_bg = NULL;
    free(sim->Agx_bg);   sim->Agx_bg   = NULL;
    free(sim->Agy_bg);   sim->Agy_bg   = NULL;
    free(sim->phi_bg);   sim->phi_bg   = NULL;
    free(sim->Ax_bg);    sim->Ax_bg    = NULL;
    free(sim->Ay_bg);    sim->Ay_bg    = NULL;
}

/* Lazily allocate and zero a background array if not already present. */
static float* ensure_bg_alloc(gr_sim_t* sim, gr_field_id_t which) {
    float** slot = background_slot(sim, which);
    if (!slot) return NULL;
    if (*slot) return *slot;
    const size_t n = (size_t) sim->width * (size_t) sim->height;
    *slot = (float*) calloc(n, sizeof(float));
    return *slot;
}

/* Eq. (eq:bg_softened_point_mass) — gr_sandbox_v33.tex §12.6
 * "Stage 6: Sampled background field arrays".
 *
 *   Phi_g^{bg}(x) = -G_eff * M / sqrt(|x - x0|^2 + epsilon^2)
 *
 * Cell-centered sampling: cell (i,j) covers [i, i+1) x [j, j+1) in cell units,
 * sampled at its center (i + 0.5, j + 0.5) * dx, matching the wave_pulse
 * scenario's convention (see scenarios/wave_pulse.c). The smoothing length
 * epsilon avoids the 1/r singularity; ~few cells is typical. */
void gr_sim_set_background_point_mass(gr_sim_t* sim,
                                      float x0, float y0,
                                      float GM, float epsilon) {
    if (!sim) return;
    float* phi_bg = ensure_bg_alloc(sim, GR_FIELD_PHI_GRAV);
    if (!phi_bg) return;

    const int   W    = sim->width;
    const int   H    = sim->height;
    const float dx   = sim->dx;
    const float eps2 = epsilon * epsilon;

    for (int j = 0; j < H; j++) {
        const float y  = ((float) j + 0.5f) * dx;
        const float dy = y - y0;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = ((float) i + 0.5f) * dx;
            const float dxi = x - x0;
            const float r2  = dxi * dxi + dy * dy;
            phi_bg[row + i] = -GM / sqrtf(r2 + eps2);
        }
    }
}
