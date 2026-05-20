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

/* Generic per-cell leapfrog body, parameterized on the four neighbor reads.
 * Used to share the leapfrog logic between zero-Dirichlet (interior loop +
 * boundary-zeroing) and periodic (full-grid loop with wrap-around indexing)
 * variants. */
static inline void leapfrog_cell(float* next, const float* prev, const float* curr,
                                 const float* src, float sc,
                                 int k, float xm, float xp, float ym, float yp,
                                 float damp_k, float c2dt2, float inv_dx2,
                                 int damped) {
    const float sum_x = xm + xp;
    const float sum_y = ym + yp;
    const float lap   = ((sum_x + sum_y) - 4.0f * curr[k]) * inv_dx2;
    const float upd   = 2.0f * curr[k] - prev[k] + c2dt2 * (lap + sc * src[k]);
    next[k] = damped ? upd * (1.0f - damp_k) : upd;
}

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
            const int k = row + i;
            leapfrog_cell(next, prev, curr, src, sc, k,
                          curr[k - 1], curr[k + 1], curr[k - W], curr[k + W],
                          damp[k], c2dt2, inv_dx2, /*damped=*/1);
        }
    }
    /* Zero-Dirichlet boundary; the damping layer absorbs incident energy
     * before it reaches the wall (§9.6). */
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
            const int k = row + i;
            leapfrog_cell(next, prev, curr, src, sc, k,
                          curr[k - 1], curr[k + 1], curr[k - W], curr[k + W],
                          0.0f, c2dt2, inv_dx2, /*damped=*/0);
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

/* Periodic-BC leapfrog: wrap-around indexing at all four edges.  Updates
 * every cell (no boundary-zeroing), restoring translation invariance of
 * the discrete Laplacian. */
static void leapfrog_field_periodic(gr_field_state_t* f, const float* damp,
                                    int W, int H, float c2dt2, float inv_dx2) {
    const float* prev = f->prev;
    const float* curr = f->curr;
    float*       next = f->next;
    const float* src  = f->source;
    const float  sc   = f->source_coeff;
    const int    damped = (damp != NULL);
    for (int j = 0; j < H; j++) {
        const int jm = (j == 0) ? (H - 1) : (j - 1);
        const int jp = (j == H - 1) ? 0 : (j + 1);
        const int row  = j  * W;
        const int rowm = jm * W;
        const int rowp = jp * W;
        for (int i = 0; i < W; i++) {
            const int im = (i == 0) ? (W - 1) : (i - 1);
            const int ip = (i == W - 1) ? 0 : (i + 1);
            const int k = row + i;
            leapfrog_cell(next, prev, curr, src, sc, k,
                          curr[row + im], curr[row + ip],
                          curr[rowm + i], curr[rowp + i],
                          damped ? damp[k] : 0.0f, c2dt2, inv_dx2, damped);
        }
    }
}

void gr_field_leapfrog_step_all(struct gr_sim* sim) {
    const int   W       = sim->width;
    const int   H       = sim->height;
    const float c2dt2   = sim->c_eff * sim->c_eff * sim->dt * sim->dt;
    const float inv_dx2 = 1.0f / (sim->dx * sim->dx);
    const float* damp   = sim->damping_d;  /* may be NULL */
    if (sim->periodic_bc) {
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            (void) gr_array_lattice((gr_array_id_t) f);
            leapfrog_field_periodic(&sim->fields[f], damp, W, H, c2dt2, inv_dx2);
        }
        return;
    }
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
