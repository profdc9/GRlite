/* Leapfrog FDTD update for all six potentials.
 * Spec reference: gr_sandbox_vNN.tex §9.2 + §9.7 (source coupling) +
 * §sec:yee_pivot (v35 layout).
 *
 * Per §9 each potential lives on its own Yee sublattice:
 *   Phi_g, phi_em : GR_LATTICE_CORNER     - nodes at (i,  j  ) * dx
 *   A_{g,x}, A_x  : GR_LATTICE_X_EDGE     - nodes at (i+0.5, j) * dx
 *   A_{g,y}, A_y  : GR_LATTICE_Y_EDGE     - nodes at (i, j+0.5) * dx
 *
 * The 5-point discrete Laplacian at each field's own index (i,j) sums
 * neighbors at (i±1, j) and (i, j±1) — all on the SAME sublattice.  So
 * the leapfrog kernel is identical in form regardless of which sublattice
 * the field lives on; only the physical interpretation of (i,j) changes.
 *
 * Boundary cells in storage:
 *   - For CORNER fields, the cells at i=0, i=W-1, j=0, j=H-1 ARE the
 *     physical box corners; zero-Dirichlet there is correct.
 *   - For X_EDGE fields, the cell at i=W-1 is a "ghost" (no physical
 *     position, outside box); zero is correct.  The cell at i=0 is at
 *     position (0.5)*dx — actually a half-cell inside the box — but the
 *     leapfrog treats it as zero-Dirichlet too.  This is a half-cell
 *     discretization error confined within the damping layer (which
 *     starts at i ~ N_d cells and damps the field there anyway), and is
 *     below other numerical noise.  Same for Y_EDGE at j=0.
 *   - Proper per-sublattice boundary handling (evolving the i=0 / j=0
 *     edge-field cells via a virtual zero ghost at i=-1 / j=-1) is
 *     deferred to S5 if force-evaluation accuracy requires it.
 *
 * Applied identically to all six potentials each step (§9.6: "the same
 * damping array d_{i,j} applies to all six because they all satisfy the
 * same wave equation with the same propagation speed").  The per-field
 * variation lives entirely in source_coeff (e.g. -4 pi G_eff for Phi_g vs
 * -4 pi G_eff / c^2 for A_g) and the source array bound to f->source. */

#include "grlite.h"
#include "sim_internal.h"

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
            /* Sum x-pair and y-pair separately, then combine.  This preserves
             * the lattice's reflection-and-transpose (D4) symmetry under float
             * arithmetic: cells related by the symmetry produce bit-identical
             * lap values because the pairwise sums avoid the x-before-y bias
             * of a left-to-right (((curr[k-1]+curr[k+1])+curr[k-W])+curr[k+W])
             * evaluation.  Critical for HE self-force stability over many
             * steps (v35 §sec:yee_pivot, post-S6 stability fix). */
            const float sum_x = curr[k - 1] + curr[k + 1];
            const float sum_y = curr[k - W] + curr[k + W];
            const float lap   = ((sum_x + sum_y) - 4.0f * curr[k]) * inv_dx2;
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
            /* Sum x-pair and y-pair separately, then combine.  This preserves
             * the lattice's reflection-and-transpose (D4) symmetry under float
             * arithmetic: cells related by the symmetry produce bit-identical
             * lap values because the pairwise sums avoid the x-before-y bias
             * of a left-to-right (((curr[k-1]+curr[k+1])+curr[k-W])+curr[k+W])
             * evaluation.  Critical for HE self-force stability over many
             * steps (v35 §sec:yee_pivot, post-S6 stability fix). */
            const float sum_x = curr[k - 1] + curr[k + 1];
            const float sum_y = curr[k - W] + curr[k + W];
            const float lap   = ((sum_x + sum_y) - 4.0f * curr[k]) * inv_dx2;
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
            /* Sublattice consulted via gr_array_lattice((gr_array_id_t) f) —
             * not used by the kernel currently (loop bounds and stencil are
             * sublattice-invariant), but stays available for downstream
             * stages that may want sublattice-aware boundary handling. */
            (void) gr_array_lattice((gr_array_id_t) f);
            leapfrog_field_damped(&sim->fields[f], damp, W, H, c2dt2, inv_dx2);
        }
    } else {
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            (void) gr_array_lattice((gr_array_id_t) f);
            leapfrog_field_undamped(&sim->fields[f], W, H, c2dt2, inv_dx2);
        }
    }
}
