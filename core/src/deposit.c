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

/* TSC (W_3 / quadratic B-spline) deposit + interp on the CORNER sublattice.
 *
 * 1D shape, anchored at nearest-integer corner inew with u = x_p/dx - inew
 * in [-0.5, 0.5):
 *   w_left   = 0.5 * (0.5 - u)^2
 *   w_center = 0.75 - u^2
 *   w_right  = 0.5 * (0.5 + u)^2
 * Weights sum to 1 by construction.  2D shape is the outer product, giving
 * a 3x3 = 9-cell footprint smoother than the CIC 2x2 = 4-cell footprint.
 *
 * Same kernel used for both deposit and interp preserves the Hockney-
 * Eastwood adjoint condition end-to-end at any sub-cell particle position
 * (combined with centered FD on corners for gradient evaluation). */
static inline void tsc_weights_1d(float u, float w[3]) {
    const float a = 0.5f - u;  /* in (0, 1] */
    const float b = 0.5f + u;  /* in [0, 1) */
    w[0] = 0.5f * a * a;
    w[1] = 0.75f - u * u;
    w[2] = 0.5f * b * b;
}

void gr_tsc_deposit_corner(float* arr, int W, int H, float dx,
                           float x_p, float y_p, float value) {
    if (!arr) return;
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);  /* nearest corner */
    const int   jc = (int) floorf(yn + 0.5f);
    /* Need indices ic-1, ic, ic+1 (resp. jc) in range. */
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) return;
    const float u = xn - (float) ic;  /* in [-0.5, 0.5) */
    const float v = yn - (float) jc;
    float wx[3], wy[3];
    tsc_weights_1d(u, wx);
    tsc_weights_1d(v, wy);
    const float inv_area = 1.0f / (dx * dx);
    for (int dj = -1; dj <= 1; dj++) {
        const int j = jc + dj;
        const int row = j * W;
        for (int di = -1; di <= 1; di++) {
            arr[row + ic + di] += value * wx[di + 1] * wy[dj + 1] * inv_area;
        }
    }
}

float gr_tsc_interp_corner(const float* arr, int W, int H, float dx,
                           float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], wy[3];
    tsc_weights_1d(u, wx);
    tsc_weights_1d(v, wy);
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += wx[di + 1] * wy[dj + 1] * arr[row + ic + di];
        }
    }
    return result;
}

/* TSC interp on X_EDGE sublattice -- nodes at (i+0.5, j) * dx.
 * Used to gather edge-staggered fields (A_x, J_x) at a particle position
 * with the W_3 kernel matched to gr_tsc_deposit_*.  Hockney-Eastwood
 * adjoint pair to the X_EDGE deposit (which Esirkepov produces). */
float gr_tsc_interp_xedge(const float* arr, int W, int H, float dx,
                          float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx - 0.5f;     /* node index space for X_EDGE */
    const float yn = y_p / dx;            /* corner-like in y */
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    /* Need ic-1, ic, ic+1 and jc-1, jc, jc+1.  X_EDGE has W-1 valid columns
     * (column W-1 is the ghost), so ic in [1, W-3].  jc in [1, H-2]. */
    if (ic < 1 || ic > W - 3 || jc < 1 || jc > H - 2) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], wy[3];
    tsc_weights_1d(u, wx);
    tsc_weights_1d(v, wy);
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += wx[di + 1] * wy[dj + 1] * arr[row + ic + di];
        }
    }
    return result;
}

/* TSC interp on Y_EDGE sublattice -- nodes at (i, j+0.5) * dx. */
float gr_tsc_interp_yedge(const float* arr, int W, int H, float dx,
                          float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx;            /* corner-like in x */
    const float yn = y_p / dx - 0.5f;     /* node index space for Y_EDGE */
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 3) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], wy[3];
    tsc_weights_1d(u, wx);
    tsc_weights_1d(v, wy);
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += wx[di + 1] * wy[dj + 1] * arr[row + ic + di];
        }
    }
    return result;
}

/* Derivative of the 1D W_3 kernel (analytic gradient of the quadratic
 * B-spline).  For u in [-0.5, 0.5) the 1D weights are
 *   w0(u) = 0.5 * (0.5 - u)^2,   w1(u) = 0.75 - u^2,   w2(u) = 0.5 * (0.5 + u)^2
 * which differentiate to
 *   dw0/du = -(0.5 - u) = u - 0.5
 *   dw1/du = -2 u
 *   dw2/du =  (0.5 + u) = u + 0.5
 * Used by the Lewis-Birdsall-style gather on edge sublattices (e.g. the
 * direct curl A read for the magnetic force at the particle). */
static inline void tsc_dw_1d(float u, float dw[3]) {
    dw[0] = u - 0.5f;
    dw[1] = -2.0f * u;
    dw[2] = u + 0.5f;
}

/* ===================================================================
 * Bump (C-infinity compactly-supported) kernel on cell-distance
 * support [-1.5, 1.5].  Sub-exponential Fourier decay; 3-cell
 * footprint matched to TSC.  See gr_sandbox_v38 sec kernel_design.
 * ===================================================================
 *
 * Raw form: b_raw(d) = exp(-1/(1 - (d/1.5)^2)) for |d| < 1.5, else 0.
 * Discrete normalization: w_i = b_raw(i - u) / sum_j b_raw(j - u)
 * so cell weights sum to 1 at every sub-cell offset u (charge
 * conservation).  Derivative includes the renormalization correction
 * so the HE-adjoint paring with the deposit is preserved. */

static inline float bump_raw(float d) {
    const float q = d / 1.5f;
    const float one_minus_q2 = 1.0f - q * q;
    if (one_minus_q2 <= 0.0f) return 0.0f;
    return expf(-1.0f / one_minus_q2);
}

/* d/dd of raw bump function. */
static inline float bump_raw_dd(float d) {
    const float q = d / 1.5f;
    const float one_minus_q2 = 1.0f - q * q;
    if (one_minus_q2 <= 0.0f) return 0.0f;
    const float b = expf(-1.0f / one_minus_q2);
    return b * (-2.0f * q) / (1.5f * one_minus_q2 * one_minus_q2);
}

/* Normalized discrete bump weights at cell offsets {-1, 0, +1} from
 * the nearest cell.  u in [-0.5, 0.5).  Sum is 1 by construction. */
static inline void bump_weights_1d(float u, float w[3]) {
    const float b0 = bump_raw(-1.0f - u);
    const float b1 = bump_raw( 0.0f - u);
    const float b2 = bump_raw( 1.0f - u);
    const float S = b0 + b1 + b2;
    if (S > 0.0f) {
        const float invS = 1.0f / S;
        w[0] = b0 * invS;
        w[1] = b1 * invS;
        w[2] = b2 * invS;
    } else {
        w[0] = w[1] = w[2] = 0.0f;
    }
}

/* d/du of the normalized discrete bump weights, including the
 * renormalization correction.  Used by the LB-style gradient gather
 * so the discrete adjoint of the deposit is exact. */
static inline void bump_dw_1d(float u, float dw[3]) {
    const float b0 = bump_raw(-1.0f - u);
    const float b1 = bump_raw( 0.0f - u);
    const float b2 = bump_raw( 1.0f - u);
    /* d/du of b_raw(i - u) = -b_raw'(i - u). */
    const float db0 = -bump_raw_dd(-1.0f - u);
    const float db1 = -bump_raw_dd( 0.0f - u);
    const float db2 = -bump_raw_dd( 1.0f - u);
    const float S = b0 + b1 + b2;
    const float dSdu = db0 + db1 + db2;
    if (S > 0.0f) {
        const float invS = 1.0f / S;
        const float W0 = b0 * invS, W1 = b1 * invS, W2 = b2 * invS;
        dw[0] = (db0 - W0 * dSdu) * invS;
        dw[1] = (db1 - W1 * dSdu) * invS;
        dw[2] = (db2 - W2 * dSdu) * invS;
    } else {
        dw[0] = dw[1] = dw[2] = 0.0f;
    }
}

/* Bump deposit on the CORNER sublattice.  Mirrors gr_tsc_deposit_corner
 * exactly except for the weight function. */
void gr_bump_deposit_corner(float* arr, int W, int H, float dx,
                             float x_p, float y_p, float value) {
    if (!arr) return;
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) return;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], wy[3];
    bump_weights_1d(u, wx);
    bump_weights_1d(v, wy);
    const float inv_area = 1.0f / (dx * dx);
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            arr[row + ic + di] += value * wx[di + 1] * wy[dj + 1] * inv_area;
        }
    }
}

/* Bump interp on the CORNER sublattice (used for diagnostics). */
float gr_bump_interp_corner(const float* arr, int W, int H, float dx,
                             float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], wy[3];
    bump_weights_1d(u, wx);
    bump_weights_1d(v, wy);
    float r = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            r += wx[di + 1] * wy[dj + 1] * arr[row + ic + di];
        }
    }
    return r;
}

/* Bump-LB gather of d/dx and d/dy of a CORNER field at the particle
 * position.  Adjoint of gr_bump_deposit_corner -- HE static self-force
 * should be zero at v=0 when both deposit and gather use bump. */
void gr_bump_lb_grad_corner(const float* arr, int W, int H, float dx,
                             float x_p, float y_p,
                             float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!arr) return;
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) return;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], wy[3], dwx[3], dwy[3];
    bump_weights_1d(u, wx);
    bump_weights_1d(v, wy);
    bump_dw_1d(u, dwx);
    bump_dw_1d(v, dwy);
    const float inv_dx = 1.0f / dx;
    float gx = 0.0f, gy = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            const float val = arr[row + ic + di];
            /* LB convention (matches TSC path in particle.c
             * phi_em_grad_at_lb): bump_dw_1d returns d/du of w_i(u),
             * and d phi / d x_p = (1/dx) * sum dw_i/du * phi_i. */
            gx += dwx[di + 1] * wy[dj + 1] * val;
            gy += wx[di + 1] * dwy[dj + 1] * val;
        }
    }
    *gx_out = inv_dx * gx;
    *gy_out = inv_dx * gy;
}

/* TSC-LB direct gather of d/dx of an X_EDGE field at the particle:
 *   dA_x/dx (x_p) = sum_{i,j} dW_x(x_p - x_{i+1/2,j}) * W_y(y_p - y_j) * A_{i+1/2,j} / dx
 * Used by the LB-style v x B path so the magnetic force inherits the same
 * adjoint pairing as the LB phi gradient. */
float gr_tsc_lb_dx_xedge(const float* arr, int W, int H, float dx,
                         float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx - 0.5f;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 3 || jc < 1 || jc > H - 2) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float dwx[3], wy[3];
    tsc_dw_1d(u, dwx);
    tsc_weights_1d(v, wy);
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += dwx[di + 1] * wy[dj + 1] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
}

/* TSC-LB direct gather of d/dy of an X_EDGE field at the particle.
 * Needed by GEMPIC for the mixed-component gradient d A_x / d y_p
 * (see gr_sandbox_v37 sec gempic_derivation). */
float gr_tsc_lb_dy_xedge(const float* arr, int W, int H, float dx,
                         float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx - 0.5f;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 3 || jc < 1 || jc > H - 2) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], dwy[3];
    tsc_weights_1d(u, wx);
    tsc_dw_1d(v, dwy);
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += wx[di + 1] * dwy[dj + 1] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
}

/* TSC-LB direct gather of d/dx of a Y_EDGE field at the particle.
 * Needed by GEMPIC for the mixed-component gradient d A_y / d x_p
 * (see gr_sandbox_v37 sec gempic_derivation). */
float gr_tsc_lb_dx_yedge(const float* arr, int W, int H, float dx,
                         float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx;
    const float yn = y_p / dx - 0.5f;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 3) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float dwx[3], wy[3];
    tsc_dw_1d(u, dwx);
    tsc_weights_1d(v, wy);
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += dwx[di + 1] * wy[dj + 1] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
}

/* TSC-LB direct gather of d/dy of a Y_EDGE field at the particle. */
float gr_tsc_lb_dy_yedge(const float* arr, int W, int H, float dx,
                         float x_p, float y_p) {
    if (!arr) return 0.0f;
    const float xn = x_p / dx;
    const float yn = y_p / dx - 0.5f;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 3) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[3], dwy[3];
    tsc_weights_1d(u, wx);
    tsc_dw_1d(v, dwy);
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -1; dj <= 1; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -1; di <= 1; di++) {
            result += wx[di + 1] * dwy[dj + 1] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
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
