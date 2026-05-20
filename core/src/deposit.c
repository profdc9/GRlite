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

/* 1D CIC (W_2) kernel: w(u) = max(1 - |u|, 0). */
static inline float cic_w1(float u) {
    const float au = (u < 0.0f) ? -u : u;
    return (au < 1.0f) ? (1.0f - au) : 0.0f;
}

/* Esirkepov 2D current deposit (W_2 / CIC shape), 2-cell case.
 * Esirkepov, J. Comp. Phys. Comm. 135 (2001) 144-153, eq. 30.
 *
 * Splits the 2D CIC shape change DS_2D = S^n_x*S^n_y - S^{n-1}_x*S^{n-1}_y
 * into orthogonal x-flux and y-flux parts:
 *   DS_x[i,j] = (S^n_x - S^{n-1}_x) * (S^{n-1}_y + 0.5 * (S^n_y - S^{n-1}_y))
 *   DS_y[i,j] = (S^n_y - S^{n-1}_y) * (S^{n-1}_x + 0.5 * (S^n_x - S^{n-1}_x))
 * which satisfy DS_x + DS_y = DS_2D exactly.  The currents are then
 *   J_x(x-edge i, j) = -(source/(dt dx)) * sum_{k <= i} DS_x[k, j]
 *   J_y(i, y-edge j) = -(source/(dt dx)) * sum_{l <= j} DS_y[i, l]
 * yielding discrete continuity (rho^n - rho^{n-1})/dt + div(J^{n-1/2}) = 0
 * at every corner cell. */
int gr_esirkepov_deposit_jxy(float* Jx, float* Jy,
                             int W, int H, float dx, float dt,
                             float x0, float y0, float x1, float y1,
                             float source) {
    if (!Jx || !Jy || dx <= 0.0f || dt <= 0.0f) return 0;
    if (source == 0.0f) return 1;
    const float xn0 = x0 / dx;
    const float yn0 = y0 / dx;
    const float xn1 = x1 / dx;
    const float yn1 = y1 / dx;

    /* 2-cell case: motion strictly less than 1 cell in each direction. */
    {
        const float ax = xn1 - xn0;
        const float ay = yn1 - yn0;
        if ((ax >  1.0f) || (ax < -1.0f) ||
            (ay >  1.0f) || (ay < -1.0f)) return 0;
    }

    /* Anchor at leftmost/bottommost cell so the 4-corner patch covers both
     * endpoints' CIC support (each is 2 wide; union <= 3; patch of 4 is safe). */
    const int im_lo = (int) floorf(fminf(xn0, xn1));
    const int jm_lo = (int) floorf(fminf(yn0, yn1));

    float S0x[4], S1x[4], S0y[4], S1y[4];
    for (int k = 0; k < 4; k++) {
        const float ix = (float) (im_lo + k);
        const float jy = (float) (jm_lo + k);
        S0x[k] = cic_w1(xn0 - ix);
        S1x[k] = cic_w1(xn1 - ix);
        S0y[k] = cic_w1(yn0 - jy);
        S1y[k] = cic_w1(yn1 - jy);
    }

    const float prefactor = -source / (dt * dx);

    /* J_x on X_EDGE: cumulative sum along i for each j-corner row. */
    for (int kj = 0; kj < 4; kj++) {
        const int jc = jm_lo + kj;
        if (jc < 0 || jc >= H) continue;
        const float Wy = S0y[kj] + 0.5f * (S1y[kj] - S0y[kj]);
        float cumsum = 0.0f;
        for (int ki = 0; ki < 4; ki++) {
            const float DSx = (S1x[ki] - S0x[ki]) * Wy;
            cumsum += DSx;
            const int ie = im_lo + ki;
            if (ie >= 0 && ie < W - 1) {
                Jx[jc * W + ie] += prefactor * cumsum;
            }
        }
    }
    /* J_y on Y_EDGE: cumulative sum along j for each i-corner column. */
    for (int ki = 0; ki < 4; ki++) {
        const int ic = im_lo + ki;
        if (ic < 0 || ic >= W) continue;
        const float Wx = S0x[ki] + 0.5f * (S1x[ki] - S0x[ki]);
        float cumsum = 0.0f;
        for (int kj = 0; kj < 4; kj++) {
            const float DSy = (S1y[kj] - S0y[kj]) * Wx;
            cumsum += DSy;
            const int je = jm_lo + kj;
            if (je >= 0 && je < H - 1) {
                Jy[je * W + ic] += prefactor * cumsum;
            }
        }
    }
    return 1;
}
