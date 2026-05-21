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

/* ----------------------------------------------------------------------------
 * PML (Hu split-field) leapfrog
 *
 * Each scalar potential Phi is split into Phi_x + Phi_y, evolved separately:
 *
 *   (1 + sigma_x dt/2) Phi_x^{n+1} = 2 Phi_x^n - (1 - sigma_x dt/2) Phi_x^{n-1}
 *                                  + c^2 dt^2 (d_x^2 Phi^n + S^n / 2)
 *   (1 + sigma_y dt/2) Phi_y^{n+1} = 2 Phi_y^n - (1 - sigma_y dt/2) Phi_y^{n-1}
 *                                  + c^2 dt^2 (d_y^2 Phi^n + S^n / 2)
 *
 * with Phi = Phi_x + Phi_y, sigma_x = sigma_x(i), sigma_y = sigma_y(j) ramped
 * cubically to zero at the PML/interior interface.  In the interior (both
 * sigma=0) summing the two recovers the standard wave-equation leapfrog
 * bit-for-bit (subject to floating-point associativity).
 *
 * Interior reads d_x^2 / d_y^2 from the assembled sum field f->curr, which
 * gr_field_pml_rotate_and_assemble keeps in sync after each step. */
static void leapfrog_field_pml(gr_field_state_t* f,
                               const float* sigma_dt_x, const float* sigma_dt_y,
                               int W, int H, float c2dt2, float inv_dx2) {
    const float* curr_sum = f->curr;      /* Phi^n = phi_x^n + phi_y^n */
    const float* px_prev  = f->phi_x_prev;
    const float* px_curr  = f->phi_x_curr;
    float*       px_next  = f->phi_x_next;
    const float* py_prev  = f->phi_y_prev;
    const float* py_curr  = f->phi_y_curr;
    float*       py_next  = f->phi_y_next;
    const float* src      = f->source;
    const float  sc       = f->source_coeff;

    for (int j = 1; j < H - 1; j++) {
        const int   row  = j * W;
        const float sy   = sigma_dt_y[j];
        const float ay   = 1.0f / (1.0f + sy);
        const float by   = (1.0f - sy);
        for (int i = 1; i < W - 1; i++) {
            const int   k     = row + i;
            const float sx    = sigma_dt_x[i];
            const float ax    = 1.0f / (1.0f + sx);
            const float bx    = (1.0f - sx);
            /* Direction-resolved second differences of the assembled field. */
            const float lap_x = (curr_sum[k - 1] + curr_sum[k + 1] - 2.0f * curr_sum[k]) * inv_dx2;
            const float lap_y = (curr_sum[k - W] + curr_sum[k + W] - 2.0f * curr_sum[k]) * inv_dx2;
            const float half_src = 0.5f * sc * src[k];
            px_next[k] = ax * (2.0f * px_curr[k] - bx * px_prev[k]
                               + c2dt2 * (lap_x + half_src));
            py_next[k] = ay * (2.0f * py_curr[k] - by * py_prev[k]
                               + c2dt2 * (lap_y + half_src));
        }
    }
    /* Zero-Dirichlet at the literal outer wall.  In a properly tuned PML the
     * field at the wall is already attenuated to ~R; the wall BC is a
     * formality, not a physical absorber. */
    for (int i = 0; i < W; i++) {
        px_next[i]               = 0.0f;
        py_next[i]               = 0.0f;
        px_next[(H - 1) * W + i] = 0.0f;
        py_next[(H - 1) * W + i] = 0.0f;
    }
    for (int j = 0; j < H; j++) {
        px_next[j * W]             = 0.0f;
        py_next[j * W]             = 0.0f;
        px_next[j * W + (W - 1)]   = 0.0f;
        py_next[j * W + (W - 1)]   = 0.0f;
    }
}

void gr_field_pml_rotate_and_assemble(struct gr_sim* sim) {
    const size_t n = (size_t) sim->width * (size_t) sim->height;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        gr_field_state_t* fs = &sim->fields[f];
        /* Rotate split-field time levels. */
        float* tx        = fs->phi_x_prev;
        fs->phi_x_prev   = fs->phi_x_curr;
        fs->phi_x_curr   = fs->phi_x_next;
        fs->phi_x_next   = tx;
        float* ty        = fs->phi_y_prev;
        fs->phi_y_prev   = fs->phi_y_curr;
        fs->phi_y_curr   = fs->phi_y_next;
        fs->phi_y_next   = ty;
        /* Refresh assembled Phi^n = phi_x^n + phi_y^n into the existing
         * curr slot so all downstream readers (force, gauge, web renderer)
         * keep working unchanged.  Also rotate prev -> next for callers
         * that still expect a time-derivative via (curr - prev) / dt
         * (notably the gauge residual diagnostic). */
        float* tmp_sum   = fs->prev;
        fs->prev         = fs->curr;
        fs->curr         = tmp_sum;
        for (size_t k = 0; k < n; k++) {
            fs->curr[k] = fs->phi_x_curr[k] + fs->phi_y_curr[k];
        }
    }
}

void gr_field_leapfrog_step_all(struct gr_sim* sim) {
    const int   W       = sim->width;
    const int   H       = sim->height;
    const float c2dt2   = sim->c_eff * sim->c_eff * sim->dt * sim->dt;
    const float inv_dx2 = 1.0f / (sim->dx * sim->dx);
    const float* damp   = sim->damping_d;  /* may be NULL */
    if (sim->n_pml > 0) {
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            (void) gr_array_lattice((gr_array_id_t) f);
            leapfrog_field_pml(&sim->fields[f],
                               sim->pml_sigma_dt_x, sim->pml_sigma_dt_y,
                               W, H, c2dt2, inv_dx2);
        }
        return;
    }
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
