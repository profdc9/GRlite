/* Scenario "wave_pulse": narrow Gaussian initial condition for the scalar potential.
 * Used by the Stage 1 test (gr_sandbox_v32.tex §12.1 "Stage 1: Scalar wave equation,
 * free propagation") and by the web frontend's first end-to-end visualization.
 *
 * Parameters (all in simulation units; defaults applied when value <= 0):
 *   params[0]: sigma     — Gaussian width             (default: 4 * dx)
 *   params[1]: amplitude — peak value at t = 0        (default: 1.0)
 *   params[2]: x0        — center x coordinate        (default: width  * dx / 2)
 *   params[3]: y0        — center y coordinate        (default: height * dx / 2)
 *
 * Initial condition: Phi(x, t=0) = amplitude * exp(-|x - x0|^2 / (2 sigma^2)) with
 * zero initial time derivative (gr_sandbox_v32.tex §9.7 "Field initialization":
 * for a static IC, set Phi^{-1} = Phi^0 so the leapfrog free-propagation
 * extrapolation 2 Phi^0 - Phi^{-1} is consistent). The Cauchy problem in 2D
 * then produces an outgoing wavefront whose amplitude falls as 1/sqrt(r) for
 * r >> sigma (§10.3 "The 2D Green's function").
 */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>

static int build_wave_pulse(gr_sim_t* sim, const float* params, int n_params) {
    const int W   = sim->width;
    const int H   = sim->height;
    const float dx = sim->dx;

    const float sigma     = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 4.0f * dx;
    const float amplitude = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 1.0f;
    const float x0        = (n_params >= 3) ? params[2] : ((float) W * 0.5f) * dx;
    const float y0        = (n_params >= 4) ? params[3] : ((float) H * 0.5f) * dx;

    const float inv_2sigma2 = 1.0f / (2.0f * sigma * sigma);

    float* phi_curr = sim->fields[GR_FIELD_PHI_GRAV].curr;
    float* phi_prev = sim->fields[GR_FIELD_PHI_GRAV].prev;
    for (int j = 0; j < H; j++) {
        const float y = ((float) j + 0.5f) * dx;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float x  = ((float) i + 0.5f) * dx;
            const float r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
            const float v  = amplitude * expf(-r2 * inv_2sigma2);
            phi_curr[row + i] = v;
            phi_prev[row + i] = v;  /* zero initial time derivative — §9.7 */
        }
    }
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_WAVE_PULSE = {
    .name  = "wave_pulse",
    .build = build_wave_pulse,
};

void gr_scenario_register_wave_pulse(void) {
    gr_scenario_register(&SCENARIO_WAVE_PULSE);
}
