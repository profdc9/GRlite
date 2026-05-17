/* Leapfrog FDTD update for the scalar potential.
 * Spec reference: gr_sandbox_v32.tex §9.2 "Time stepping: the leapfrog scheme". */

#include "grlite.h"
#include "sim_internal.h"

/* Eq. (eq:leapfrog_field) — §9.2 "Time stepping: the leapfrog scheme".
 *
 * Three-point leapfrog (Stoermer-Verlet) update of the homogeneous d'Alembert
 * equation, Box Phi = 0:
 *
 *     Phi^{n+1} = 2 Phi^n - Phi^{n-1} + (c dt)^2 ( Lap Phi^n + source^n )
 *
 * Stage 1 has no source term — source = 0. The Laplacian uses the standard
 * 5-point isotropic stencil (gr_sandbox_v32.tex §9.4, the dispersion analysis
 * around eq:dispersion_exact assumes this same 5-point form):
 *
 *     (Lap Phi)_{i,j} = ( Phi_{i-1,j} + Phi_{i+1,j} + Phi_{i,j-1} + Phi_{i,j+1}
 *                          - 4 Phi_{i,j} ) / dx^2
 *
 * Boundary treatment for Stage 1: hard reflecting (zero-Dirichlet) edges. The
 * design doc §9.6 prescribes an absorbing damping layer for production runs,
 * but Stage 1 deliberately omits it (§12.1) — tests must be sized so the
 * outgoing pulse does not reach the wall during the measurement window.
 */
void gr_field_leapfrog_step(struct gr_sim* sim) {
    const int W = sim->width;
    const int H = sim->height;

    const float c2dt2   = sim->c_eff * sim->c_eff * sim->dt * sim->dt;
    const float inv_dx2 = 1.0f / (sim->dx * sim->dx);

    const float* prev = sim->phi_prev;
    const float* curr = sim->phi_curr;
    float* next       = sim->phi_next;

    for (int j = 1; j < H - 1; j++) {
        const int row = j * W;
        for (int i = 1; i < W - 1; i++) {
            const int k    = row + i;
            const float lap = (curr[k - 1] + curr[k + 1] + curr[k - W] + curr[k + W]
                              - 4.0f * curr[k]) * inv_dx2;
            next[k] = 2.0f * curr[k] - prev[k] + c2dt2 * lap;
        }
    }

    /* Zero-Dirichlet boundary (Stage 1: hard wall — pulse reflection beyond
     * the measurement window is tolerated). */
    for (int i = 0; i < W; i++) {
        next[i] = 0.0f;
        next[(H - 1) * W + i] = 0.0f;
    }
    for (int j = 0; j < H; j++) {
        next[j * W] = 0.0f;
        next[j * W + (W - 1)] = 0.0f;
    }
}
