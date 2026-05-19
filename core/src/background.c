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
    /* Clear analytic generator parameters too — they're paired with the
     * sampled array.  Mode is left alone (user-set). */
    sim->bg_kind   = GR_BG_KIND_NONE;
    sim->bg_x0     = 0.0f;
    sim->bg_y0     = 0.0f;
    sim->bg_GM     = 0.0f;
    sim->bg_eps    = 0.0f;
    sim->bg_charge = 0.0f;
    sim->bg_spin   = 0.0f;
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
    /* Store generator parameters for the analytic-mode evaluator. */
    sim->bg_kind = GR_BG_KIND_POINT_MASS;
    sim->bg_x0   = x0;
    sim->bg_y0   = y0;
    sim->bg_GM   = GM;
    sim->bg_eps  = epsilon;
}

/* Analytic-mode evaluation of the installed background generator at the
 * particle's exact (x, y).  See gr_sandbox §12.6 / §12.8 for the rationale:
 * the sampled CIC+FD path introduces an O((dx/r)^2) tangential force error
 * that breaks angular-momentum conservation in test-particle orbits.  The
 * analytic path replaces that with the closed-form expression, eliminating
 * the grid-induced precession in Stage 7/8 tests entirely. */
int gr_bg_eval_analytic(const struct gr_sim* sim, float x, float y,
                        float* phi_out, float* gx_out, float* gy_out) {
    if (!sim || sim->bg_kind == GR_BG_KIND_NONE) return 0;

    switch (sim->bg_kind) {
    case GR_BG_KIND_POINT_MASS: {
        /* Softened Newtonian potential (eq:bg_softened_point_mass): */
        /*   Phi(x,y) = -G*M / sqrt(r^2 + eps^2)                    */
        /*   grad     =  G*M * (r - r0) / (r^2 + eps^2)^{3/2}       */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        const float r2  = dxi * dxi + dyi * dyi;
        const float s2  = r2 + sim->bg_eps * sim->bg_eps;
        const float inv_s  = 1.0f / sqrtf(s2);
        const float inv_s3 = inv_s / s2;            /* 1 / (s2)^{3/2} */
        *phi_out = -sim->bg_GM * inv_s;
        *gx_out  =  sim->bg_GM * dxi * inv_s3;
        *gy_out  =  sim->bg_GM * dyi * inv_s3;
        return 1;
    }
    default:
        return 0;
    }
}
