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
 * Bump (C-infinity compactly-supported) kernel.  PARAMETERIZED by
 * radius R in cell units; cells touched per axis = 2*ceil(R+0.5)+1
 * (capped at BUMP_MAX_HW).  Default R=1.5 matches TSC's 3-cell support.
 * Larger R = wider macroparticle, alternative to rho_smooth_passes.
 * See gr_sandbox_v38 sec kernel_design.
 * ===================================================================
 *
 * Raw form (R-parameterized):
 *   b_raw(d, R) = exp(-1/(1 - (d/R)^2)) for |d| < R, else 0
 * Discrete weights: w_i = b_raw(i - u, R) / sum_j b_raw(j - u, R)
 * (per-axis renormalization so weights sum to 1 -> charge conservation).
 * Derivative includes the renormalization correction so the HE-adjoint
 * pairing with the deposit is preserved. */

#define BUMP_MAX_HW 8                       /* max half-width in cells */
#define BUMP_MAX_W  (2 * BUMP_MAX_HW + 1)   /* max cells/axis = 17 */

static inline float bump_raw_R(float d, float R) {
    const float q = d / R;
    const float one_minus_q2 = 1.0f - q * q;
    if (one_minus_q2 <= 0.0f) return 0.0f;
    return expf(-1.0f / one_minus_q2);
}

/* d/dd of raw bump (used for LB derivative). */
static inline float bump_raw_dd_R(float d, float R) {
    const float q = d / R;
    const float one_minus_q2 = 1.0f - q * q;
    if (one_minus_q2 <= 0.0f) return 0.0f;
    const float b = expf(-1.0f / one_minus_q2);
    return b * (-2.0f * q) / (R * one_minus_q2 * one_minus_q2);
}

/* Half-width in cells given R.  Cells in window: 2*hw+1. */
static inline int bump_half_width(float R) {
    int hw = (int) ceilf(R + 0.5f);
    if (hw < 1) hw = 1;
    if (hw > BUMP_MAX_HW) hw = BUMP_MAX_HW;
    return hw;
}

/* Fill weight array w[0..2*hw] for cells at offsets [-hw, hw] from
 * the nearest cell.  Per-axis renormalized to sum to 1. */
static inline void bump_weights_1d_R(float u, float R, int hw, float w[BUMP_MAX_W]) {
    float S = 0.0f;
    for (int i = -hw; i <= hw; i++) {
        const float b = bump_raw_R((float) i - u, R);
        w[i + hw] = b;
        S += b;
    }
    if (S > 0.0f) {
        const float invS = 1.0f / S;
        for (int i = -hw; i <= hw; i++) w[i + hw] *= invS;
    } else {
        for (int i = -hw; i <= hw; i++) w[i + hw] = 0.0f;
    }
}

/* Fill derivative array dw[0..2*hw] = d/du of normalized weights,
 * including the renormalization correction. */
static inline void bump_dw_1d_R(float u, float R, int hw,
                                  float w[BUMP_MAX_W], float dw[BUMP_MAX_W]) {
    float b_raw_arr[BUMP_MAX_W], db_raw_arr[BUMP_MAX_W];
    float S = 0.0f, dSdu = 0.0f;
    for (int i = -hw; i <= hw; i++) {
        const float b  =  bump_raw_R   ((float) i - u, R);
        const float db = -bump_raw_dd_R((float) i - u, R);  /* d/du of b((i-u)) */
        b_raw_arr [i + hw] = b;
        db_raw_arr[i + hw] = db;
        S    += b;
        dSdu += db;
    }
    if (S > 0.0f) {
        const float invS = 1.0f / S;
        for (int i = -hw; i <= hw; i++) {
            const float W = b_raw_arr[i + hw] * invS;
            w[i + hw]  = W;
            dw[i + hw] = (db_raw_arr[i + hw] - W * dSdu) * invS;
        }
    } else {
        for (int i = -hw; i <= hw; i++) { w[i + hw] = 0.0f; dw[i + hw] = 0.0f; }
    }
}

/* Bump deposit on the CORNER sublattice. */
void gr_bump_deposit_corner(float* arr, int W, int H, float dx,
                             float x_p, float y_p, float value, float R) {
    if (!arr) return;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 1 - hw || jc < hw || jc > H - 1 - hw) return;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W];
    bump_weights_1d_R(u, R, hw, wx);
    bump_weights_1d_R(v, R, hw, wy);
    const float inv_area = 1.0f / (dx * dx);
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            arr[row + ic + di] += value * wx[di + hw] * wy[dj + hw] * inv_area;
        }
    }
}

/* Bump interp on the CORNER sublattice (used for diagnostics). */
float gr_bump_interp_corner(const float* arr, int W, int H, float dx,
                             float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 1 - hw || jc < hw || jc > H - 1 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W];
    bump_weights_1d_R(u, R, hw, wx);
    bump_weights_1d_R(v, R, hw, wy);
    float r = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            r += wx[di + hw] * wy[dj + hw] * arr[row + ic + di];
        }
    }
    return r;
}

/* Bump interp on X_EDGE sublattice -- nodes at (i+0.5, j) * dx.
 * Used by Boris EM force evaluation: A_x reads for d_t A and curl A. */
float gr_bump_interp_xedge(const float* arr, int W, int H, float dx,
                            float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx - 0.5f;       /* X_EDGE: half-cell offset in x */
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    /* X_EDGE valid range: i in [hw, W-2-hw]; jc in [hw, H-1-hw]. */
    if (ic < hw || ic > W - 2 - hw || jc < hw || jc > H - 1 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W];
    bump_weights_1d_R(u, R, hw, wx);
    bump_weights_1d_R(v, R, hw, wy);
    float r = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            r += wx[di + hw] * wy[dj + hw] * arr[row + ic + di];
        }
    }
    return r;
}

/* Bump interp on Y_EDGE sublattice -- nodes at (i, j+0.5) * dx. */
float gr_bump_interp_yedge(const float* arr, int W, int H, float dx,
                            float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx;
    const float yn = y_p / dx - 0.5f;       /* Y_EDGE: half-cell offset in y */
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 1 - hw || jc < hw || jc > H - 2 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W];
    bump_weights_1d_R(u, R, hw, wx);
    bump_weights_1d_R(v, R, hw, wy);
    float r = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            r += wx[di + hw] * wy[dj + hw] * arr[row + ic + di];
        }
    }
    return r;
}

/* Bump Esirkepov current deposit on X_EDGE / Y_EDGE.  Variable-window
 * version of gr_esirkepov_deposit_jxy.  Uses the renormalized bump
 * weights so that the discrete continuity equation
 *   (rho^{n+1} - rho^n)/dt + div(J^{n-1/2}) = 0
 * holds exactly at every corner, provided rho is deposited with the
 * same bump weights (which it is when GR_SHAPE_BUMP is selected).
 *
 * Algebra parallel to gr_esirkepov_deposit_jxy:
 *   DS_x[i, j] = (S^1_x - S^0_x)[i] * (S^0_y + 0.5 * (S^1_y - S^0_y))[j]
 *   DS_y[i, j] = (S^1_y - S^0_y)[j] * (S^0_x + 0.5 * (S^1_x - S^0_x))[i]
 *   J_x at X_edge i, j = -(source/(dt*dx)) * sum_{k<=i} DS_x[k, j]
 *   J_y at Y_edge i, j = -(source/(dt*dx)) * sum_{l<=j} DS_y[i, l]
 *
 * The cell window covers the union of bump support at x^n and x^{n+1}.
 * For sub-CFL motion (|v*dt| < dx) and bump half-width hw, the union
 * is at most hw + 1 cells in each direction beyond the central cell. */
int gr_bump_esirkepov_deposit_jxy(float* Jx, float* Jy,
                                   int W, int H, float dx, float dt,
                                   float x0, float y0, float x1, float y1,
                                   float source, float R) {
    if (!Jx || !Jy || dx <= 0.0f || dt <= 0.0f) return 0;
    if (source == 0.0f) return 1;
    const int hw = bump_half_width(R);

    const float xn0 = x0 / dx;
    const float yn0 = y0 / dx;
    const float xn1 = x1 / dx;
    const float yn1 = y1 / dx;

    /* Reject motion larger than one cell per axis (CFL ought to keep it
     * sub-cell anyway; this is a safety check matching the CIC version). */
    const float ax = xn1 - xn0;
    const float ay = yn1 - yn0;
    if (ax > 1.0f || ax < -1.0f || ay > 1.0f || ay < -1.0f) return 0;

    /* Window: cover both endpoints' bump support.  Anchor at leftmost
     * cell touched by either endpoint.  Endpoint 0 touches cells
     * [floor(xn0) - hw, floor(xn0) + hw]; endpoint 1 similarly.  Union
     * extends from min(floor(xn0), floor(xn1)) - hw to max(...) + hw.
     * Total width <= 2*hw + 2 cells (since motion is sub-cell). */
    const int ic0 = (int) floorf(xn0 + 0.5f);
    const int ic1 = (int) floorf(xn1 + 0.5f);
    const int jc0 = (int) floorf(yn0 + 0.5f);
    const int jc1 = (int) floorf(yn1 + 0.5f);
    const int ic_min = (ic0 < ic1) ? ic0 : ic1;
    const int jc_min = (jc0 < jc1) ? jc0 : jc1;
    const int ic_max = (ic0 > ic1) ? ic0 : ic1;
    const int jc_max = (jc0 > jc1) ? jc0 : jc1;
    const int im_lo = ic_min - hw;
    const int jm_lo = jc_min - hw;
    const int im_hi = ic_max + hw;
    const int jm_hi = jc_max + hw;
    const int N_x = im_hi - im_lo + 1;
    const int N_y = jm_hi - jm_lo + 1;
    /* Width upper bound: 2*hw + 2 (one extra column for the +1 cell shift). */
    if (N_x > 2 * BUMP_MAX_HW + 2 || N_y > 2 * BUMP_MAX_HW + 2) return 0;

    /* Compute renormalized weights at endpoints 0 and 1, for the
     * uniform cell window [im_lo..im_hi] x [jm_lo..jm_hi]. */
    float S0x[2 * BUMP_MAX_HW + 2], S1x[2 * BUMP_MAX_HW + 2];
    float S0y[2 * BUMP_MAX_HW + 2], S1y[2 * BUMP_MAX_HW + 2];
    /* For endpoint 0: nearest cell is ic0, sub-cell u = xn0 - ic0.  But
     * here we want weights at the SAME window cells [im_lo..im_hi] for
     * both endpoints, so we compute raw bump_raw_R values at d = i - xn
     * for i in window, then renormalize so the sum over the window is 1. */
    float S0x_raw[2 * BUMP_MAX_HW + 2], S1x_raw[2 * BUMP_MAX_HW + 2];
    float S0y_raw[2 * BUMP_MAX_HW + 2], S1y_raw[2 * BUMP_MAX_HW + 2];
    float Sx0_sum = 0.0f, Sx1_sum = 0.0f, Sy0_sum = 0.0f, Sy1_sum = 0.0f;
    for (int k = 0; k < N_x; k++) {
        const float ix = (float) (im_lo + k);
        S0x_raw[k] = bump_raw_R(ix - xn0, R);
        S1x_raw[k] = bump_raw_R(ix - xn1, R);
        Sx0_sum += S0x_raw[k];
        Sx1_sum += S1x_raw[k];
    }
    for (int k = 0; k < N_y; k++) {
        const float jy = (float) (jm_lo + k);
        S0y_raw[k] = bump_raw_R(jy - yn0, R);
        S1y_raw[k] = bump_raw_R(jy - yn1, R);
        Sy0_sum += S0y_raw[k];
        Sy1_sum += S1y_raw[k];
    }
    const float invSx0 = (Sx0_sum > 0.0f) ? 1.0f / Sx0_sum : 0.0f;
    const float invSx1 = (Sx1_sum > 0.0f) ? 1.0f / Sx1_sum : 0.0f;
    const float invSy0 = (Sy0_sum > 0.0f) ? 1.0f / Sy0_sum : 0.0f;
    const float invSy1 = (Sy1_sum > 0.0f) ? 1.0f / Sy1_sum : 0.0f;
    for (int k = 0; k < N_x; k++) {
        S0x[k] = S0x_raw[k] * invSx0;
        S1x[k] = S1x_raw[k] * invSx1;
    }
    for (int k = 0; k < N_y; k++) {
        S0y[k] = S0y_raw[k] * invSy0;
        S1y[k] = S1y_raw[k] * invSy1;
    }

    const float prefactor = -source / (dt * dx);

    /* J_x on X_EDGE: cumulative sum along i for each j-corner row. */
    for (int kj = 0; kj < N_y; kj++) {
        const int jc = jm_lo + kj;
        if (jc < 0 || jc >= H) continue;
        const float Wy = S0y[kj] + 0.5f * (S1y[kj] - S0y[kj]);
        float cumsum = 0.0f;
        for (int ki = 0; ki < N_x; ki++) {
            const float DSx = (S1x[ki] - S0x[ki]) * Wy;
            cumsum += DSx;
            const int ie = im_lo + ki;
            if (ie >= 0 && ie < W - 1) {
                Jx[jc * W + ie] += prefactor * cumsum;
            }
        }
    }
    /* J_y on Y_EDGE: cumulative sum along j for each i-corner column. */
    for (int ki = 0; ki < N_x; ki++) {
        const int ic = im_lo + ki;
        if (ic < 0 || ic >= W) continue;
        const float Wx = S0x[ki] + 0.5f * (S1x[ki] - S0x[ki]);
        float cumsum = 0.0f;
        for (int kj = 0; kj < N_y; kj++) {
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

/* Bump-LB gather of d/dx and d/dy of a CORNER field at the particle. */
void gr_bump_lb_grad_corner(const float* arr, int W, int H, float dx,
                             float x_p, float y_p, float R,
                             float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!arr) return;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 1 - hw || jc < hw || jc > H - 1 - hw) return;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W], dwx[BUMP_MAX_W], dwy[BUMP_MAX_W];
    bump_dw_1d_R(u, R, hw, wx, dwx);
    bump_dw_1d_R(v, R, hw, wy, dwy);
    const float inv_dx = 1.0f / dx;
    float gx = 0.0f, gy = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            const float val = arr[row + ic + di];
            gx += dwx[di + hw] * wy[dj + hw] * val;
            gy += wx[di + hw] * dwy[dj + hw] * val;
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

/* Bump-LB mixed-axis gathers on edge sublattices.  Four functions for the
 * full Tier-2 coverage: d/dx and d/dy of X_EDGE and Y_EDGE fields at the
 * particle.  Used for the bump-shaped curl A and (in case of future
 * use) GEMPIC-style force evaluation.  All use the renormalized bump
 * weights with the corresponding derivative kernel. */

float gr_bump_lb_dx_xedge(const float* arr, int W, int H, float dx,
                           float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx - 0.5f;       /* X_EDGE */
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 2 - hw || jc < hw || jc > H - 1 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W], dwx[BUMP_MAX_W], dwy_unused[BUMP_MAX_W];
    bump_dw_1d_R(u, R, hw, wx, dwx);
    bump_weights_1d_R(v, R, hw, wy);
    (void) dwy_unused;
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            result += dwx[di + hw] * wy[dj + hw] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
}

float gr_bump_lb_dy_xedge(const float* arr, int W, int H, float dx,
                           float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx - 0.5f;
    const float yn = y_p / dx;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 2 - hw || jc < hw || jc > H - 1 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W], dwx_unused[BUMP_MAX_W], dwy[BUMP_MAX_W];
    bump_weights_1d_R(u, R, hw, wx);
    bump_dw_1d_R(v, R, hw, wy, dwy);
    (void) dwx_unused;
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            result += wx[di + hw] * dwy[dj + hw] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
}

float gr_bump_lb_dx_yedge(const float* arr, int W, int H, float dx,
                           float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx;
    const float yn = y_p / dx - 0.5f;       /* Y_EDGE */
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 1 - hw || jc < hw || jc > H - 2 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W], dwx[BUMP_MAX_W], dwy_unused[BUMP_MAX_W];
    bump_dw_1d_R(u, R, hw, wx, dwx);
    bump_weights_1d_R(v, R, hw, wy);
    (void) dwy_unused;
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            result += dwx[di + hw] * wy[dj + hw] * arr[row + ic + di];
        }
    }
    return result * inv_dx;
}

float gr_bump_lb_dy_yedge(const float* arr, int W, int H, float dx,
                           float x_p, float y_p, float R) {
    if (!arr) return 0.0f;
    const int   hw = bump_half_width(R);
    const float xn = x_p / dx;
    const float yn = y_p / dx - 0.5f;
    const int   ic = (int) floorf(xn + 0.5f);
    const int   jc = (int) floorf(yn + 0.5f);
    if (ic < hw || ic > W - 1 - hw || jc < hw || jc > H - 2 - hw) return 0.0f;
    const float u = xn - (float) ic;
    const float v = yn - (float) jc;
    float wx[BUMP_MAX_W], wy[BUMP_MAX_W], dwx_unused[BUMP_MAX_W], dwy[BUMP_MAX_W];
    bump_weights_1d_R(u, R, hw, wx);
    bump_dw_1d_R(v, R, hw, wy, dwy);
    (void) dwx_unused;
    const float inv_dx = 1.0f / dx;
    float result = 0.0f;
    for (int dj = -hw; dj <= hw; dj++) {
        const int row = (jc + dj) * W;
        for (int di = -hw; di <= hw; di++) {
            result += wx[di + hw] * dwy[dj + hw] * arr[row + ic + di];
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
