/* Source deposition kernels.
 * Spec reference: gr_sandbox_v32.tex §9.5 "Source deposition: approximating
 * the delta function". */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>

/* Eq. (eq:cic_deposit) — §9.5 "Source deposition" (CIC / W_2 tent kernel).
 *
 * Distribute `value` (e.g. a mass M, or charge q) into the four cells
 * surrounding sub-cell position (x_p, y_p), with bilinear weights. Cell
 * (i, j) has its center at ((i + 0.5) * dx, (j + 0.5) * dx) in our
 * cell-centered convention (matches wave_pulse.c and background.c). Letting
 *
 *   alpha = x_p/dx - 0.5 - ic        in [0, 1)
 *   beta  = y_p/dx - 0.5 - jc        in [0, 1)
 *
 * the deposits are
 *   rho[ic    , jc    ] += value * (1 - alpha) * (1 - beta) / dx^2
 *   rho[ic + 1, jc    ] += value *      alpha  * (1 - beta) / dx^2
 *   rho[ic    , jc + 1] += value * (1 - alpha) *      beta  / dx^2
 *   rho[ic + 1, jc + 1] += value *      alpha  *      beta  / dx^2
 *
 * The total integral sum_k rho[k] * dx^2 equals `value` exactly (the four
 * weights sum to 1 by construction). Positions outside the depositable
 * interior (ic in [0, W-2], jc in [0, H-2]) are silently dropped — this
 * matches the Stage 1 design choice of confining sources to the interior
 * (§9.6 absorbing layer doesn't damp deposited rho). */
void gr_cic_deposit_scalar(float* rho, int W, int H, float dx,
                           float x_p, float y_p, float value) {
    if (!rho) return;
    const float xn = x_p / dx - 0.5f;
    const float yn = y_p / dx - 0.5f;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 1) return;
    const float alpha    = xn - (float) ic;
    const float beta     = yn - (float) jc;
    const float inv_area = 1.0f / (dx * dx);

    rho[jc       * W + ic    ] += value * (1.0f - alpha) * (1.0f - beta) * inv_area;
    rho[jc       * W + ic + 1] += value *         alpha  * (1.0f - beta) * inv_area;
    rho[(jc + 1) * W + ic    ] += value * (1.0f - alpha) *         beta  * inv_area;
    rho[(jc + 1) * W + ic + 1] += value *         alpha  *         beta  * inv_area;
}
