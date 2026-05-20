/* Source deposition kernels.
 * Spec reference: gr_sandbox_vNN.tex §9.5 "Source deposition: approximating
 * the delta function" + §sec:yee_pivot (v35 sublattice layout).
 *
 * Per §9 each source array lives on its own Yee sublattice:
 *   rho_matter, rho_q       : GR_LATTICE_CORNER     - nodes at (i,  j  ) * dx
 *   J_{m,x}, J_{q,x}        : GR_LATTICE_X_EDGE     - nodes at (i+0.5, j) * dx
 *   J_{m,y}, J_{q,y}        : GR_LATTICE_Y_EDGE     - nodes at (i, j+0.5) * dx
 *
 * CIC ($W_2$) weights are the same bilinear form on every sublattice; only
 * the sub-cell fraction is measured from the sublattice's own node positions.
 * Given a particle at $(x_p, y_p)$:
 *
 *   alpha = x_p/dx - offset_x - ic   in [0, 1)
 *   beta  = y_p/dx - offset_y - jc   in [0, 1)
 *
 * with (offset_x, offset_y) = (0,0) for CORNER, (0.5, 0) for X_EDGE,
 * (0, 0.5) for Y_EDGE.  The four-cell deposit:
 *
 *   arr[ic    , jc    ] += value * (1 - alpha) * (1 - beta) * inv_area
 *   arr[ic + 1, jc    ] += value *      alpha  * (1 - beta) * inv_area
 *   arr[ic    , jc + 1] += value * (1 - alpha) *      beta  * inv_area
 *   arr[ic + 1, jc + 1] += value *      alpha  *      beta  * inv_area
 *
 * Total deposited sum_k arr[k] * dx^2 == value (the four bilinear weights
 * sum to 1 by construction).
 *
 * Out-of-range positions (ic, jc outside the valid stencil range for the
 * sublattice) are silently dropped.  Valid ranges in storage:
 *   CORNER: ic in [0, W-2], jc in [0, H-2]
 *   X_EDGE: ic in [0, W-3], jc in [0, H-2]   (i = W-1 is ghost, never deposit)
 *   Y_EDGE: ic in [0, W-2], jc in [0, H-3]   (j = H-1 is ghost, never deposit) */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>

/* CORNER sublattice: sample at (i, j) * dx.  Used for rho_matter, rho_q. */
void gr_cic_deposit_corner(float* arr, int W, int H, float dx,
                           float x_p, float y_p, float value) {
    if (!arr) return;
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 1) return;
    const float alpha    = xn - (float) ic;
    const float beta     = yn - (float) jc;
    const float inv_area = 1.0f / (dx * dx);

    arr[jc       * W + ic    ] += value * (1.0f - alpha) * (1.0f - beta) * inv_area;
    arr[jc       * W + ic + 1] += value *         alpha  * (1.0f - beta) * inv_area;
    arr[(jc + 1) * W + ic    ] += value * (1.0f - alpha) *         beta  * inv_area;
    arr[(jc + 1) * W + ic + 1] += value *         alpha  *         beta  * inv_area;
}

/* X_EDGE sublattice: sample at (i+0.5, j) * dx.  Used for J_mx, J_qx, and
 * the implicit deposition target for d/dx of corner-sublattice fields. */
void gr_cic_deposit_xedge(float* arr, int W, int H, float dx,
                          float x_p, float y_p, float value) {
    if (!arr) return;
    const float xn = x_p / dx - 0.5f;     /* subtract sublattice offset */
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    /* Last valid x-edge cell is i = W-2; cell W-1 is the ghost.  The +1
     * neighbor in x must stay within [0, W-2], so ic in [0, W-3]. */
    if (ic < 0 || ic >= W - 2 || jc < 0 || jc >= H - 1) return;
    const float alpha    = xn - (float) ic;
    const float beta     = yn - (float) jc;
    const float inv_area = 1.0f / (dx * dx);

    arr[jc       * W + ic    ] += value * (1.0f - alpha) * (1.0f - beta) * inv_area;
    arr[jc       * W + ic + 1] += value *         alpha  * (1.0f - beta) * inv_area;
    arr[(jc + 1) * W + ic    ] += value * (1.0f - alpha) *         beta  * inv_area;
    arr[(jc + 1) * W + ic + 1] += value *         alpha  *         beta  * inv_area;
}

/* Y_EDGE sublattice: sample at (i, j+0.5) * dx.  Used for J_my, J_qy. */
void gr_cic_deposit_yedge(float* arr, int W, int H, float dx,
                          float x_p, float y_p, float value) {
    if (!arr) return;
    const float xn = x_p / dx;
    const float yn = y_p / dx - 0.5f;     /* subtract sublattice offset */
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    /* Last valid y-edge cell is j = H-2; cell H-1 is the ghost. */
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 2) return;
    const float alpha    = xn - (float) ic;
    const float beta     = yn - (float) jc;
    const float inv_area = 1.0f / (dx * dx);

    arr[jc       * W + ic    ] += value * (1.0f - alpha) * (1.0f - beta) * inv_area;
    arr[jc       * W + ic + 1] += value *         alpha  * (1.0f - beta) * inv_area;
    arr[(jc + 1) * W + ic    ] += value * (1.0f - alpha) *         beta  * inv_area;
    arr[(jc + 1) * W + ic + 1] += value *         alpha  *         beta  * inv_area;
}
