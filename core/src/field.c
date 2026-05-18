/* Leapfrog FDTD update for all six potentials.
 * Spec reference: gr_sandbox_v32.tex §9.2 + §9.7 (source coupling). */

#include "grlite.h"
#include "sim_internal.h"

/* Eq. (eq:leapfrog_field) with the §9.7 source term:
 *   Phi^{n+1} = 2 Phi^n - Phi^{n-1} + (c dt)^2 (Lap Phi^n + source_coeff * source^n)
 *
 * Applied identically to all six potentials each step (gr_sandbox_v32.tex §9.6:
 * "the same damping array d_{i,j} applies to all six because they all satisfy
 * the same wave equation with the same propagation speed"). The per-field
 * variation lives entirely in source_coeff (e.g. -4 pi G_eff for Phi_g vs
 * -4 pi G_eff / c^2 for A_g) and the source array bound to f->source. */
static void leapfrog_field_damped(gr_field_state_t* f, const float* damp,
                                  int W, int H, float c2dt2, float inv_dx2) {
    const float* prev = f->prev;
    const float* curr = f->curr;
    float*       next = f->next;
    const float* src  = f->source;
    const float  sc   = f->source_coeff;
    for (int j = 1; j < H - 1; j++) {
        const int row = j * W;
        for (int i = 1; i < W - 1; i++) {
            const int   k   = row + i;
            const float lap = (curr[k - 1] + curr[k + 1] + curr[k - W] + curr[k + W]
                              - 4.0f * curr[k]) * inv_dx2;
            next[k] = (2.0f * curr[k] - prev[k]
                      + c2dt2 * (lap + sc * src[k])) * (1.0f - damp[k]);
        }
    }
    /* Zero-Dirichlet boundary; same justification as Stage 1/2 — the damping
     * layer absorbs incident energy before it reaches the wall (§9.6). */
    for (int i = 0; i < W; i++) {
        next[i] = 0.0f;
        next[(H - 1) * W + i] = 0.0f;
    }
    for (int j = 0; j < H; j++) {
        next[j * W] = 0.0f;
        next[j * W + (W - 1)] = 0.0f;
    }
}

static void leapfrog_field_undamped(gr_field_state_t* f,
                                    int W, int H, float c2dt2, float inv_dx2) {
    const float* prev = f->prev;
    const float* curr = f->curr;
    float*       next = f->next;
    const float* src  = f->source;
    const float  sc   = f->source_coeff;
    for (int j = 1; j < H - 1; j++) {
        const int row = j * W;
        for (int i = 1; i < W - 1; i++) {
            const int   k   = row + i;
            const float lap = (curr[k - 1] + curr[k + 1] + curr[k - W] + curr[k + W]
                              - 4.0f * curr[k]) * inv_dx2;
            next[k] = 2.0f * curr[k] - prev[k] + c2dt2 * (lap + sc * src[k]);
        }
    }
    for (int i = 0; i < W; i++) {
        next[i] = 0.0f;
        next[(H - 1) * W + i] = 0.0f;
    }
    for (int j = 0; j < H; j++) {
        next[j * W] = 0.0f;
        next[j * W + (W - 1)] = 0.0f;
    }
}

void gr_field_leapfrog_step_all(struct gr_sim* sim) {
    const int   W       = sim->width;
    const int   H       = sim->height;
    const float c2dt2   = sim->c_eff * sim->c_eff * sim->dt * sim->dt;
    const float inv_dx2 = 1.0f / (sim->dx * sim->dx);
    const float* damp   = sim->damping_d;  /* may be NULL */
    if (damp) {
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            leapfrog_field_damped(&sim->fields[f], damp, W, H, c2dt2, inv_dx2);
        }
    } else {
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            leapfrog_field_undamped(&sim->fields[f], W, H, c2dt2, inv_dx2);
        }
    }
}
