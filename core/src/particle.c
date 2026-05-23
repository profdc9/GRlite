/* Particle integrator — relativistic Boris-leapfrog (kick-drift) in 2D.
 * Spec reference: gr_sandbox_v32.tex §9.2 (eq:leapfrog_field neighborhood),
 * §9.5 (CIC adjoint condition for force interpolation), and v33 §12.7. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>

/* Sublattice-aware CIC interpolators (v35 §sec:yee_pivot, answer C1):
 * three inline variants matching the deposit kernels gr_cic_deposit_*.
 * Pairing deposit+interp on the same sublattice with the same bilinear
 * weights satisfies the §9.5 adjoint condition for Hockney-Eastwood
 * self-force cancellation.  Each variant computes its sub-cell offset
 * against its sublattice's own node positions:
 *   CORNER:  (i,    j   ) * dx  - no offset
 *   X_EDGE:  (i+0.5, j  ) * dx  - x offset 0.5
 *   Y_EDGE:  (i,    j+0.5) * dx - y offset 0.5
 */
static inline float cic_interp_corner(const float* f, int W, int H, float dx,
                                      float x, float y) {
    if (!f) return 0.0f;
    const float xn = x / dx;
    const float yn = y / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 1) return 0.0f;
    const float alpha = xn - (float) ic;
    const float beta  = yn - (float) jc;
    return (1.0f - alpha) * (1.0f - beta) * f[jc       * W + ic]
         +         alpha  * (1.0f - beta) * f[jc       * W + ic + 1]
         + (1.0f - alpha) *         beta  * f[(jc + 1) * W + ic]
         +         alpha  *         beta  * f[(jc + 1) * W + ic + 1];
}

static inline float cic_interp_xedge(const float* f, int W, int H, float dx,
                                     float x, float y) {
    if (!f) return 0.0f;
    const float xn = x / dx - 0.5f;
    const float yn = y / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 2 || jc < 0 || jc >= H - 1) return 0.0f;
    const float alpha = xn - (float) ic;
    const float beta  = yn - (float) jc;
    return (1.0f - alpha) * (1.0f - beta) * f[jc       * W + ic]
         +         alpha  * (1.0f - beta) * f[jc       * W + ic + 1]
         + (1.0f - alpha) *         beta  * f[(jc + 1) * W + ic]
         +         alpha  *         beta  * f[(jc + 1) * W + ic + 1];
}

static inline float cic_interp_yedge(const float* f, int W, int H, float dx,
                                     float x, float y) {
    if (!f) return 0.0f;
    const float xn = x / dx;
    const float yn = y / dx - 0.5f;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 2) return 0.0f;
    const float alpha = xn - (float) ic;
    const float beta  = yn - (float) jc;
    return (1.0f - alpha) * (1.0f - beta) * f[jc       * W + ic]
         +         alpha  * (1.0f - beta) * f[jc       * W + ic + 1]
         + (1.0f - alpha) *         beta  * f[(jc + 1) * W + ic]
         +         alpha  *         beta  * f[(jc + 1) * W + ic + 1];
}

float gr_phi_g_total_at(const struct gr_sim* sim, float x, float y) {
    if (!sim) return 0.0f;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    /* Phi_g lives on the CORNER sublattice (v35 §9).  Use TSC interp when
     * selected to match the TSC deposit kernel (HE adjoint pairing). */
    const int use_tsc = (sim->shape_function == GR_SHAPE_TSC);
    const float pert = use_tsc
        ? gr_tsc_interp_corner(sim->fields[GR_FIELD_PHI_GRAV].curr, W, H, dx, x, y)
        :     cic_interp_corner(sim->fields[GR_FIELD_PHI_GRAV].curr, W, H, dx, x, y);
    float bg = 0.0f;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC && sim->bg_kind != GR_BG_KIND_NONE) {
        float phi_a, gx_a, gy_a;
        if (gr_bg_eval_analytic(sim, x, y, &phi_a, &gx_a, &gy_a)) {
            bg = phi_a;
        }
    } else {
        bg = use_tsc
            ? gr_tsc_interp_corner(sim->phi_g_bg, W, H, dx, x, y)
            :     cic_interp_corner(sim->phi_g_bg, W, H, dx, x, y);
    }
    return pert + bg;
}

/* Lewis-Birdsall variant: gradient of Phi at the particle computed as
 *     dPhi/dx_p = sum_g (dW/dx)(x_p - x_g) * Phi_g
 * using the analytic gradient of the deposit kernel itself.  This is the
 * variationally-consistent force from the discrete Lagrangian — the
 * discrete work-energy theorem holds exactly for the particle moving
 * under this force in a self-consistently deposited field.
 *
 * For CIC (W_2):
 *   particle at (x_p, y_p), sub-cell offset (fx, fy) in [0, 1).
 *   4 corners involved (ic, jc), (ic+1, jc), (ic, jc+1), (ic+1, jc+1).
 *   x-weights:  -(1-fy)/dx, +(1-fy)/dx, -fy/dx, +fy/dx
 *   y-weights:  -(1-fx)/dx, -fx/dx, +(1-fx)/dx, +fx/dx
 *
 * For TSC (W_3, quadratic B-spline) anchored at the NEAREST corner ic =
 *   floor(xn + 0.5), with u = xn - ic in [-1/2, 1/2]:
 *   d/ds W_3 at the three x-cells (di = -1, 0, +1):
 *       di = -1:   (u - 1/2)/dx
 *       di =  0:  -2u/dx
 *       di = +1:   (u + 1/2)/dx
 *   In 2D the x-weight at cell (ic+di, jc+dj) is
 *       w_x(di, dj) = (dW_3/ds at di) * W_3(v at dj) / dx
 *   where W_3(v at dj) = {0.5*(0.5-v)^2, 0.75-v^2, 0.5*(0.5+v)^2}.  Symmetric
 *   for w_y.
 *
 * Both LB-CIC and LB-TSC satisfy the HE adjoint condition (F_self = 0 at
 * any sub-cell position for stationary sources) AND the discrete energy
 * conservation theorem (no per-step work imbalance under motion). */
static void grav_grad_at_lb(const struct gr_sim* sim, float x, float y,
                            float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!sim) return;
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx  = sim->dx;

    /* Background contribution. */
    float gx_bg_ana = 0.0f, gy_bg_ana = 0.0f;
    const float* bg = NULL;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC && sim->bg_kind != GR_BG_KIND_NONE) {
        float phi_a;
        gr_bg_eval_analytic(sim, x, y, &phi_a, &gx_bg_ana, &gy_bg_ana);
    } else if (sim->bg_mode == GR_BG_MODE_SAMPLED) {
        bg = sim->phi_g_bg;  /* may be NULL */
    }
    /* Read the field at time n (matched to the deposit at the start of
     * the current step), not Phi^{n+1}.  After the field leapfrog +
     * rotation in gr_sim_step, Phi^n now lives in prev; curr holds
     * Phi^{n+1}.  For LB's energy-conservation guarantee the force must
     * be evaluated from the SAME field-time as the deposit.  When the
     * field evolution is disabled (Stage 7/8 fixed-background tests),
     * prev is untouched (zero-initialized), so fall back to curr. */
    const float* pert = sim->field_evolution_enabled
                            ? sim->fields[GR_FIELD_PHI_GRAV].prev
                            : sim->fields[GR_FIELD_PHI_GRAV].curr;

    if (sim->shape_function == GR_SHAPE_TSC) {
        /* TSC anchor: nearest corner.  u, v in [-1/2, 1/2]. */
        const float xn = x / dx;
        const float yn = y / dx;
        const int   ic = (int) floorf(xn + 0.5f);
        const int   jc = (int) floorf(yn + 0.5f);
        if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) {
            *gx_out = gx_bg_ana;
            *gy_out = gy_bg_ana;
            return;
        }
        const float u  = xn - (float) ic;
        const float v  = yn - (float) jc;
        const float a  = 0.5f - v, b = 0.5f + v;
        const float c  = 0.5f - u, d  = 0.5f + u;
        /* Value weights (TSC W_3) along x and y. */
        const float wx[3] = { 0.5f * c * c,  0.75f - u * u, 0.5f * d * d };
        const float wy[3] = { 0.5f * a * a,  0.75f - v * v, 0.5f * b * b };
        /* Derivative weights dW_3/ds along x and y. */
        const float dwx[3] = { u - 0.5f, -2.0f * u, u + 0.5f };
        const float dwy[3] = { v - 0.5f, -2.0f * v, v + 0.5f };
        const float inv_dx = 1.0f / dx;
        float gx_sum = 0.0f, gy_sum = 0.0f;
        for (int dj = -1; dj <= 1; dj++) {
            const int row = (jc + dj) * W;
            for (int di = -1; di <= 1; di++) {
                const int k = row + ic + di;
                float phi = pert[k];
                if (bg) phi += bg[k];
                gx_sum += (dwx[di + 1] * wy[dj + 1]) * phi;
                gy_sum += (wx[di + 1]  * dwy[dj + 1]) * phi;
            }
        }
        *gx_out = inv_dx * gx_sum + gx_bg_ana;
        *gy_out = inv_dx * gy_sum + gy_bg_ana;
        return;
    }

    /* CIC: anchor at lower-left corner, fx, fy in [0, 1). */
    const float xn = x / dx;
    const float yn = y / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 1) {
        *gx_out = gx_bg_ana;
        *gy_out = gy_bg_ana;
        return;
    }
    const float fx = xn - (float) ic;
    const float fy = yn - (float) jc;
    const float inv_dx = 1.0f / dx;
    const int k00 =  jc      * W + ic;
    const int k10 =  jc      * W + ic + 1;
    const int k01 = (jc + 1) * W + ic;
    const int k11 = (jc + 1) * W + ic + 1;
    float p00 = pert[k00], p10 = pert[k10], p01 = pert[k01], p11 = pert[k11];
    if (bg) { p00 += bg[k00]; p10 += bg[k10]; p01 += bg[k01]; p11 += bg[k11]; }
    /* F_x = -m dPhi/dx_p; here we return dPhi/dx, dPhi/dy. */
    const float gx_lb = inv_dx * (-(1.0f - fy) * p00 + (1.0f - fy) * p10
                                  -        fy  * p01 +        fy  * p11);
    const float gy_lb = inv_dx * (-(1.0f - fx) * p00 -        fx  * p10
                                  +(1.0f - fx) * p01 +        fx  * p11);
    *gx_out = gx_lb + gx_bg_ana;
    *gy_out = gy_lb + gy_bg_ana;
}

/* Gradient of Phi_g_total at (x, y) on the CORNER sublattice.
 *
 * Centered finite difference on corners + corner-CIC interpolation to the
 * particle.  This is the Hockney-Eastwood-preserving chain: deposit kernel
 * (corner-CIC) and force-interp kernel (corner-CIC) are identical, the
 * discrete corner Laplacian Green's function is symmetric, and the centered
 * FD matrix is antisymmetric in (i, j).  Self-force on a stationary
 * particle vanishes identically for ALL sub-cell positions, not just at
 * symmetric points.
 *
 *   (d Phi / d x)_corner[c]  = (Phi[c+1] - Phi[c-1]) / (2 dx)
 *   (d Phi / d y)_corner[c]  = (Phi[c+W] - Phi[c-W]) / (2 dx)
 *
 *   particle -> rho deposit (corner-CIC) -> Phi (corner Laplacian)
 *           -> grad (centered FD on corners)
 *           -> force interp (corner-CIC back to particle).
 *
 * Trade-off vs the natural Yee corner->edge forward FD: we give up the
 * sublattice-elegance of having grad components live on edges, but in
 * exchange the adjoint condition is restored end-to-end.  The Esirkepov
 * current (J on edges) is independent and unaffected — currents still
 * sit on edges to give exact continuity at corners; only the *force*
 * evaluation goes corner->corner here. */
static void grav_grad_at(const struct gr_sim* sim, float x, float y,
                         float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!sim) return;
    /* Energy-conserving force pairing — dispatched here so the rest of
     * the pusher (Tier-2 corrections, proper-time, leapfrog init) is
     * unchanged.  Both schemes return dPhi/dx, dPhi/dy at the particle;
     * only the kernel differs. */
    if (sim->force_interp == GR_FORCE_INTERP_LEWIS_BIRDSALL) {
        grav_grad_at_lb(sim, x, y, gx_out, gy_out);
        return;
    }
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    const float inv_2dx = 1.0f / (2.0f * dx);

    /* Background contribution. */
    float gx_bg = 0.0f, gy_bg = 0.0f;
    const float* bg_corner = NULL;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC && sim->bg_kind != GR_BG_KIND_NONE) {
        float phi_a;
        gr_bg_eval_analytic(sim, x, y, &phi_a, &gx_bg, &gy_bg);
    } else if (sim->bg_mode == GR_BG_MODE_SAMPLED) {
        bg_corner = sim->phi_g_bg;  /* may be NULL — handled inside the kernel */
    }

    /* Read Phi^n (the field at the time matched to the deposit done at
     * the start of this gr_sim_step), NOT Phi^{n+1}.  After the field
     * leapfrog + rotation in gr_sim_step, Phi^n now lives in prev; curr
     * holds Phi^{n+1}.  The legacy FD-then-interp path previously read
     * curr, giving a half-step time mismatch that masked the LB benefit
     * and contributed to moving-particle heating.  Phi_g is always on
     * the corner sublattice (v35 §9).
     * When field_evolution_enabled is 0 (Stage 7/8 fixed-background runs)
     * the field is static and prev is zero-initialized, so fall back to
     * curr to preserve those tests' behavior. */
    const float* pert = sim->field_evolution_enabled
                            ? sim->fields[GR_FIELD_PHI_GRAV].prev
                            : sim->fields[GR_FIELD_PHI_GRAV].curr;

    if (sim->shape_function == GR_SHAPE_TSC) {
        /* TSC corner interp: 3x3 stencil anchored at nearest corner. */
        const float xn = x / dx;
        const float yn = y / dx;
        const int   ic = (int) floorf(xn + 0.5f);
        const int   jc = (int) floorf(yn + 0.5f);
        /* Need ic-1..ic+1 corners plus their ±1 neighbors for centered FD,
         * so ic in [2, W-3], jc in [2, H-3]. */
        if (ic >= 2 && ic < W - 2 && jc >= 2 && jc < H - 2) {
            const float u = xn - (float) ic;
            const float v = yn - (float) jc;
            const float a = 0.5f - u, b = 0.5f + u;
            const float c = 0.5f - v, d = 0.5f + v;
            const float wx[3] = {0.5f * a * a, 0.75f - u * u, 0.5f * b * b};
            const float wy[3] = {0.5f * c * c, 0.75f - v * v, 0.5f * d * d};
            float gx_sum = 0.0f, gy_sum = 0.0f;
            for (int dj = -1; dj <= 1; dj++) {
                const int row = (jc + dj) * W;
                for (int di = -1; di <= 1; di++) {
                    const int k = row + ic + di;
                    float vgx = (pert[k + 1] - pert[k - 1]) * inv_2dx;
                    float vgy = (pert[k + W] - pert[k - W]) * inv_2dx;
                    if (bg_corner) {
                        vgx += (bg_corner[k + 1] - bg_corner[k - 1]) * inv_2dx;
                        vgy += (bg_corner[k + W] - bg_corner[k - W]) * inv_2dx;
                    }
                    const float w = wx[di + 1] * wy[dj + 1];
                    gx_sum += w * vgx;
                    gy_sum += w * vgy;
                }
            }
            *gx_out = gx_sum + gx_bg;
            *gy_out = gy_sum + gy_bg;
        } else {
            *gx_out = gx_bg;
            *gy_out = gy_bg;
        }
        return;
    }

    /* Corner-CIC interp at (x, y): xn = x/dx, no sublattice offset. */
    const float xn = x / dx;
    const float yn = y / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);

    /* Need 4 corners (ic..ic+1) x (jc..jc+1), each reading a centered-FD
     * neighbor stencil one cell in each direction.  Interior margin:
     * ic in [1, W-3], jc in [1, H-3]. */
    if (ic >= 1 && ic < W - 2 && jc >= 1 && jc < H - 2) {
        const float alpha = xn - (float) ic;
        const float beta  = yn - (float) jc;
        float gx[4], gy[4];
        for (int dq = 0; dq < 2; dq++) {
            const int jj  = jc + dq;
            const int row = jj * W;
            for (int dp = 0; dp < 2; dp++) {
                const int ii = ic + dp;
                float vgx = (pert[row + ii + 1] - pert[row + ii - 1]) * inv_2dx;
                float vgy = (pert[row + W + ii] - pert[row - W + ii]) * inv_2dx;
                if (bg_corner) {
                    vgx += (bg_corner[row + ii + 1] - bg_corner[row + ii - 1]) * inv_2dx;
                    vgy += (bg_corner[row + W + ii] - bg_corner[row - W + ii]) * inv_2dx;
                }
                gx[dq * 2 + dp] = vgx;
                gy[dq * 2 + dp] = vgy;
            }
        }
        const float w00 = (1.0f - alpha) * (1.0f - beta);
        const float w10 =         alpha  * (1.0f - beta);
        const float w01 = (1.0f - alpha) *         beta;
        const float w11 =         alpha  *         beta;
        *gx_out = w00 * gx[0] + w10 * gx[1] + w01 * gx[2] + w11 * gx[3] + gx_bg;
        *gy_out = w00 * gy[0] + w10 * gy[1] + w01 * gy[2] + w11 * gy[3] + gy_bg;
    } else {
        *gx_out = gx_bg;
        *gy_out = gy_bg;
    }
}

/* Gravitational force F at a particle with mass m, velocity (vx, vy), in the
 * total Phi_g (perturbation + background, value phi at the particle), its
 * gradient (grad_x, grad_y) = grad(Phi_g_total), and the local gravitomagnetic
 * field B_g_z = (curl A_g)_z (added Stage 20 — Tier-1 gravitomagnetic Lorentz).
 * The local gravity vector is g = -grad(Phi_g).
 *
 * Tiers (gr_sandbox_v35.tex §"Practical implementation tiers" line 985):
 *   NEWTONIAN:    F = m * g                              + 4 m v x B_g
 *   RELATIVISTIC: F = m * [(1 + v^2/c^2 + 4 phi/c^2) g
 *                          - 4 (v . g) v / c^2 ]         + 4 m v x B_g
 *
 * Both tiers include the Tier-1 gravitomagnetic Lorentz piece
 *   F_gm = +4 m (v x B_g),   (v x B_g_z z)_x = v_y B_g_z, (...)_y = -v_x B_g_z
 * per gr_sandbox_v35.tex eq:geodesic_expansion (line 938) and the Tier-3
 * algorithmic eqbox (line 1040, in potential form
 *   4 grad(v . A_g) - 4 (v . grad) A_g  =  4 v x (curl A_g)
 * via the vector identity).  The factor of 4 is the spin-2 enhancement
 * relative to the EM Lorentz force (line 458, 989, 2125 in v35); v34
 * corrected v33's sign on the velocity-coupling term (line 948).  When
 * A_g = 0 (no spinning / uniform-B_g background and no perturbation A_g)
 * this term vanishes identically, recovering the older Tier-0/Tier-2 code
 * paths bit-exactly.
 *
 * The relativistic scalar expression is the Einstein-Infeld-Hoffmann 1PN
 * equation of motion for a test particle in a static field, harmonic gauge,
 * with the v33 doc's isotropic-form metric (g_{ij} = (1-2 phi/c^2) delta_ij).
 * See Ali-Haïmoud, GR Fall 2019 lecture 25, eq. 37, with psi = xi = 0 and
 * dt phi = 0:
 *
 *   dv/dt = -(1 + v^2/c^2 + 4 phi/c^2) grad(phi) + 4 (v . grad(phi)) v / c^2
 *
 * Substituting g = -grad(phi) and (v . g) = -(v . grad(phi)) gives the form
 * above.  Differs from v33 eq:geodesic_expansion (line 668) in two ways:
 *   (a) the 4*phi*g/c^2 ("Shapiro") term is added — comes from the g_{00}
 *       expansion at O(v^4) in Gamma^i_{00} (lecture eq. 30).
 *   (b) the sign on the velocity-coupling term is NEGATIVE, not positive as
 *       v33 has.  Empirical confirmation: with the v33 sign the test orbit
 *       precesses retrograde; with the EIH sign it precesses prograde. */
static inline void grav_force_at(const struct gr_sim* sim,
                                 float mass, float vx, float vy,
                                 float phi,
                                 float grad_x, float grad_y,
                                 float Bg_z,
                                 float dAg_x, float dAg_y,
                                 float* Fx, float* Fy) {
    if (sim->force_tier == GR_FORCE_RELATIVISTIC) {
        /* EIH 1PN coordinate acceleration (Ali-Haïmoud, GR Fall 2019 lecture
         * 25, eq. 37, static test particle with psi = xi = 0, dt phi = 0):
         *
         *   dv/dt = g (1 + v^2/c^2 + 4 phi/c^2)  -  4 (v . g) v / c^2
         *
         * We apply this as F/m in the Boris-leapfrog dp/dt update.  The
         * exact conversion adds a gamma factor and a (v.a) term, both small
         * (~v^2/c^2) at the velocities of interest; numerical tests at deep
         * 1PN regime (v/c ~ 0.04) recover the Schwarzschild precession to
         * ~12 percent without it, well within 2PN truncation error. */
        const float inv_c2 = 1.0f / (sim->c_eff * sim->c_eff);
        const float v2     = vx * vx + vy * vy;
        /* v . g = -(v . grad Phi_g) */
        const float v_dot_g = -(vx * grad_x + vy * grad_y);
        const float k1 = -(1.0f + (v2 + 4.0f * phi) * inv_c2);   /* (1 + v^2/c^2 + 4 phi/c^2) g */
        const float k2 = -4.0f * v_dot_g * inv_c2;               /* -4 (v.g) v / c^2 */
        *Fx = mass * (k1 * grad_x + k2 * vx);
        *Fy = mass * (k1 * grad_y + k2 * vy);
    } else {
        *Fx = -mass * grad_x;
        *Fy = -mass * grad_y;
    }
    /* Tier-1 gravitomagnetic Lorentz force, additive in all scalar tiers
     * (gr_sandbox_v35.tex eq:geodesic_expansion line 938; spin-2 factor of 4
     * at line 458/989/2125).  In 2D with B_g_z along +z:
     *   (v x B_g_z z)_x = +v_y B_g_z
     *   (v x B_g_z z)_y = -v_x B_g_z
     * Contributes only when an A_g background or perturbation is active;
     * when Bg_z = 0 the scalar tiers above are recovered bit-exactly. */
    *Fx += 4.0f * mass * vy * Bg_z;
    *Fy -= 4.0f * mass * vx * Bg_z;
    /* GM inductive piece -m d_t A_g per Eq.1040 (Tier-3 eqbox).
     * Coefficient is 1 (not 4) on the time-derivative piece, per the
     * doc convention -- the factor of 4 lives on the v x B_g spatial
     * cross-product, not on the d_t A_g term itself.  Default OFF; gated
     * by gravitomagnetic_inductive_enabled (see Stage 28). */
    *Fx -= mass * dAg_x;
    *Fy -= mass * dAg_y;
}

/* EM Lorentz force on a charged particle.  Per gr_sandbox_v35.tex
 * eq:eih_full (line 921) and the Tier-3 eqbox (line 1045):
 *
 *   F_em = q ( -grad phi - d_t A + v x B )                  (B = curl A)
 *
 * All three pieces are now implemented (Stages 23/24/25):
 *   -grad phi  — electrostatic (Stage 24).  Passed as (phi_grad_x,
 *                phi_grad_y); F_E = -q (phi_grad_x, phi_grad_y).
 *   -d_t A     — inductive E (Stage 25).  Passed as (dAx, dAy);
 *                F_ind = -q (dAx, dAy).
 *   v x B      — magnetic (Stage 23).  Passed as B_em_z (scalar in 2D);
 *                F_B = q (+v_y B_z, -v_x B_z).
 *
 * The full coefficient on v x B is +1 (no spin-2 enhancement, unlike
 * the GEM coefficient +4 on v x B_g). */
static inline void em_force_at(float charge, float vx, float vy,
                               float phi_grad_x, float phi_grad_y,
                               float dAx, float dAy,
                               float B_em_z,
                               float* Fx, float* Fy) {
    *Fx = charge * (-phi_grad_x - dAx + vy * B_em_z);
    *Fy = charge * (-phi_grad_y - dAy - vx * B_em_z);
}

/* Yee curl: B_g_z = d/dx A_{g,y} - d/dy A_{g,x} on the cell-center
 * sublattice (i+0.5, j+0.5) * dx from edge-staggered A_g arrays
 * (gr_sandbox_v35.tex §9.1, Yee staggering).  A_{g,x} lives on X_EDGE
 * (nodes at (i+0.5, j) * dx); A_{g,y} on Y_EDGE (nodes at (i, j+0.5) * dx).
 * Forward differences:
 *   d/dx A_{g,y}     = (A_gy[j, i+1] - A_gy[j, i]) / dx
 *   d/dy A_{g,x}     = (A_gx[j+1, i] - A_gx[j, i]) / dx
 * yielding B_g_z at cell center (i+0.5, j+0.5).  Bilinear interpolation of
 * the four cell-center values around the particle gives B_g_z at (x, y).
 * Returns 0 if either array is NULL or the stencil straddles the box edge. */
static float curl_Agz_sampled_at(const float* Agx, const float* Agy,
                                 int W, int H, float dx,
                                 float x, float y) {
    if (!Agx || !Agy) return 0.0f;
    const float inv_dx = 1.0f / dx;
    /* Cell-center coordinates of the particle. */
    const float u = x * inv_dx - 0.5f;
    const float v = y * inv_dx - 0.5f;
    const int   ic = (int) floorf(u);
    const int   jc = (int) floorf(v);
    /* Need cell centers (ic, jc), (ic+1, jc), (ic, jc+1), (ic+1, jc+1).
     * Curl at (ic+1, jc+1) reads A_gy at column ic+2 and A_gx at row jc+2,
     * so interior bounds are ic in [0, W-3], jc in [0, H-3]. */
    if (ic < 0 || ic > W - 3 || jc < 0 || jc > H - 3) return 0.0f;
    const float fx = u - (float) ic;
    const float fy = v - (float) jc;
    /* Inline the curl at each of the four cell centers (ic+di, jc+dj). */
    float bz[4];
    for (int dj = 0; dj < 2; dj++) {
        const int j = jc + dj;
        for (int di = 0; di < 2; di++) {
            const int i = ic + di;
            const float dAgy_dx = (Agy[j * W + (i + 1)] - Agy[j * W + i]) * inv_dx;
            const float dAgx_dy = (Agx[(j + 1) * W + i] - Agx[j * W + i]) * inv_dx;
            bz[dj * 2 + di] = dAgy_dx - dAgx_dy;
        }
    }
    const float w00 = (1.0f - fx) * (1.0f - fy);
    const float w10 =         fx  * (1.0f - fy);
    const float w01 = (1.0f - fx) *         fy;
    const float w11 =         fx  *         fy;
    return w00 * bz[0] + w10 * bz[1] + w01 * bz[2] + w11 * bz[3];
}

/* Total gravitomagnetic field B_g_z at the particle position, combining
 * background (analytic OR sampled) and perturbation (always sampled).
 * Returns 0 if the gravitomagnetic force is gated off.
 *
 * Background:
 *   - bg_mode = ANALYTIC: closed form via gr_bg_eval_B_g (handles
 *     UNIFORM_GRAVITOMAGNETIC, SPINNING_POINT_MASS; zero otherwise).
 *   - bg_mode = SAMPLED:  Yee curl of Agx_bg, Agy_bg arrays (zero when
 *     the kind has no A_g background, since those arrays remain NULL).
 *
 * Perturbation: Yee curl of fields[A_GX] and fields[A_GY] arrays.  Reads
 * .prev when field_evolution_enabled (same read-prev convention as Phi_g
 * for time-consistency with the deposit done at the start of this step);
 * falls back to .curr when field evolution is off (static fields).  When
 * the perturbation A_g potentials are unsourced and zero this contributes
 * nothing, preserving Stage 7-21 behavior bit-exactly. */
static float B_g_z_at_total(const struct gr_sim* sim, float x, float y) {
    if (!sim || !sim->gravitomagnetic_force_enabled) return 0.0f;
    float bg = 0.0f;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC) {
        gr_bg_eval_B_g(sim, x, y, &bg);
    } else {
        bg = curl_Agz_sampled_at(sim->Agx_bg, sim->Agy_bg,
                                 sim->width, sim->height, sim->dx, x, y);
    }
    /* Half-step A_g convention: perturbation A_g lives at half-integer times,
     * so .prev = A_g^{n-1/2} and .curr = A_g^{n+1/2}.  To evaluate B_g at
     * the integer time t^n where the particle lives, average the curls of
     * the two surrounding half-step slices:
     *    B_g^n ≈ ( curl A_g^{n-1/2} + curl A_g^{n+1/2} ) / 2.
     * When field_evolution is off, prev and curr are both at their static
     * values (no rotation, both buffers untouched), so the average reduces
     * to a single-slice read. */
    float pert = 0.0f;
    if (sim->field_evolution_enabled) {
        const float c1 = curl_Agz_sampled_at(sim->fields[GR_FIELD_A_GX].prev,
                                             sim->fields[GR_FIELD_A_GY].prev,
                                             sim->width, sim->height, sim->dx, x, y);
        const float c2 = curl_Agz_sampled_at(sim->fields[GR_FIELD_A_GX].curr,
                                             sim->fields[GR_FIELD_A_GY].curr,
                                             sim->width, sim->height, sim->dx, x, y);
        pert = 0.5f * (c1 + c2);
    } else {
        pert = curl_Agz_sampled_at(sim->fields[GR_FIELD_A_GX].curr,
                                   sim->fields[GR_FIELD_A_GY].curr,
                                   sim->width, sim->height, sim->dx, x, y);
    }
    return bg + pert;
}

/* Lewis-Birdsall variant of the EM scalar gradient.  Structurally identical
 * to grav_grad_at_lb but reads phi_em (and phi_bg for sampled-mode) instead
 * of Phi_g.  See grav_grad_at_lb's doc-block above for the derivation
 * details; this is the EM mirror.  Gated by em_lorentz_force_enabled. */
static void phi_em_grad_at_lb(const struct gr_sim* sim, float x, float y,
                              float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!sim || !sim->em_lorentz_force_enabled) return;
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx  = sim->dx;

    float gx_bg_ana = 0.0f, gy_bg_ana = 0.0f;
    const float* bg = NULL;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC && sim->bg_kind != GR_BG_KIND_NONE) {
        float phi_a;
        gr_bg_eval_phi_em(sim, x, y, &phi_a, &gx_bg_ana, &gy_bg_ana);
    } else if (sim->bg_mode == GR_BG_MODE_SAMPLED) {
        bg = sim->phi_bg;
    }
    const float* pert = sim->field_evolution_enabled
                            ? sim->fields[GR_FIELD_PHI_EM].prev
                            : sim->fields[GR_FIELD_PHI_EM].curr;

    if (sim->shape_function == GR_SHAPE_TSC) {
        const float xn = x / dx;
        const float yn = y / dx;
        const int   ic = (int) floorf(xn + 0.5f);
        const int   jc = (int) floorf(yn + 0.5f);
        if (ic < 1 || ic > W - 2 || jc < 1 || jc > H - 2) {
            *gx_out = gx_bg_ana;
            *gy_out = gy_bg_ana;
            return;
        }
        const float u  = xn - (float) ic;
        const float v  = yn - (float) jc;
        const float a  = 0.5f - v, b = 0.5f + v;
        const float c  = 0.5f - u, d  = 0.5f + u;
        const float wx[3]  = { 0.5f * c * c,  0.75f - u * u, 0.5f * d * d };
        const float wy[3]  = { 0.5f * a * a,  0.75f - v * v, 0.5f * b * b };
        const float dwx[3] = { u - 0.5f, -2.0f * u, u + 0.5f };
        const float dwy[3] = { v - 0.5f, -2.0f * v, v + 0.5f };
        const float inv_dx = 1.0f / dx;
        float gx_sum = 0.0f, gy_sum = 0.0f;
        for (int dj = -1; dj <= 1; dj++) {
            const int row = (jc + dj) * W;
            for (int di = -1; di <= 1; di++) {
                const int k = row + ic + di;
                float phi = pert[k];
                if (bg) phi += bg[k];
                gx_sum += (dwx[di + 1] * wy[dj + 1]) * phi;
                gy_sum += (wx[di + 1]  * dwy[dj + 1]) * phi;
            }
        }
        *gx_out = inv_dx * gx_sum + gx_bg_ana;
        *gy_out = inv_dx * gy_sum + gy_bg_ana;
        return;
    }

    /* CIC-LB */
    const float xn = x / dx;
    const float yn = y / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 1) {
        *gx_out = gx_bg_ana;
        *gy_out = gy_bg_ana;
        return;
    }
    const float fx = xn - (float) ic;
    const float fy = yn - (float) jc;
    const float inv_dx = 1.0f / dx;
    const int k00 =  jc      * W + ic;
    const int k10 =  jc      * W + ic + 1;
    const int k01 = (jc + 1) * W + ic;
    const int k11 = (jc + 1) * W + ic + 1;
    float p00 = pert[k00], p10 = pert[k10], p01 = pert[k01], p11 = pert[k11];
    if (bg) { p00 += bg[k00]; p10 += bg[k10]; p01 += bg[k01]; p11 += bg[k11]; }
    const float gx_lb = inv_dx * (-(1.0f - fy) * p00 + (1.0f - fy) * p10
                                  -        fy  * p01 +        fy  * p11);
    const float gy_lb = inv_dx * (-(1.0f - fx) * p00 -        fx  * p10
                                  +(1.0f - fx) * p01 +        fx  * p11);
    *gx_out = gx_lb + gx_bg_ana;
    *gy_out = gy_lb + gy_bg_ana;
}

/* EM scalar gradient at the particle: grad phi_em^{total} = grad phi^{bg}
 * + grad phi^{pert}.  Returns (gx, gy) such that the electrostatic force
 * piece on charge q is F = -q * (gx, gy).  Mirrors grav_grad_at structurally.
 *
 * Dispatch:
 *   force_interp == LEWIS_BIRDSALL -> phi_em_grad_at_lb (CIC-LB or TSC-LB).
 *   force_interp == LEGACY:
 *     shape_function == TSC -> TSC corner interp + centered FD.
 *     shape_function == CIC -> CIC corner interp + centered FD (the basic path).
 *
 * Gated by em_lorentz_force_enabled. */
static void phi_em_grad_at_total(const struct gr_sim* sim, float x, float y,
                                 float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!sim || !sim->em_lorentz_force_enabled) return;

    if (sim->force_interp == GR_FORCE_INTERP_LEWIS_BIRDSALL) {
        phi_em_grad_at_lb(sim, x, y, gx_out, gy_out);
        return;
    }

    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    const float inv_2dx = 1.0f / (2.0f * dx);

    float gx_bg = 0.0f, gy_bg = 0.0f;
    const float* bg_corner = NULL;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC && sim->bg_kind != GR_BG_KIND_NONE) {
        float phi_a;
        gr_bg_eval_phi_em(sim, x, y, &phi_a, &gx_bg, &gy_bg);
    } else if (sim->bg_mode == GR_BG_MODE_SAMPLED) {
        bg_corner = sim->phi_bg;
    }
    const float* pert = sim->field_evolution_enabled
                            ? sim->fields[GR_FIELD_PHI_EM].prev
                            : sim->fields[GR_FIELD_PHI_EM].curr;

    if (sim->shape_function == GR_SHAPE_TSC) {
        /* TSC corner interp: 3x3 stencil + centered FD at each neighbor. */
        const float xn = x / dx;
        const float yn = y / dx;
        const int   ic = (int) floorf(xn + 0.5f);
        const int   jc = (int) floorf(yn + 0.5f);
        if (ic >= 2 && ic < W - 2 && jc >= 2 && jc < H - 2) {
            const float u = xn - (float) ic;
            const float v = yn - (float) jc;
            const float a = 0.5f - u, b = 0.5f + u;
            const float c = 0.5f - v, d = 0.5f + v;
            const float wx[3] = {0.5f * a * a, 0.75f - u * u, 0.5f * b * b};
            const float wy[3] = {0.5f * c * c, 0.75f - v * v, 0.5f * d * d};
            float gx_sum = 0.0f, gy_sum = 0.0f;
            for (int dj = -1; dj <= 1; dj++) {
                const int row = (jc + dj) * W;
                for (int di = -1; di <= 1; di++) {
                    const int k = row + ic + di;
                    float vgx = (pert[k + 1] - pert[k - 1]) * inv_2dx;
                    float vgy = (pert[k + W] - pert[k - W]) * inv_2dx;
                    if (bg_corner) {
                        vgx += (bg_corner[k + 1] - bg_corner[k - 1]) * inv_2dx;
                        vgy += (bg_corner[k + W] - bg_corner[k - W]) * inv_2dx;
                    }
                    const float w = wx[di + 1] * wy[dj + 1];
                    gx_sum += w * vgx;
                    gy_sum += w * vgy;
                }
            }
            *gx_out = gx_sum + gx_bg;
            *gy_out = gy_sum + gy_bg;
        } else {
            *gx_out = gx_bg;
            *gy_out = gy_bg;
        }
        return;
    }

    /* CIC + centered-FD (original LEGACY path). */
    const float xn = x / dx;
    const float yn = y / dx;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 1 || ic >= W - 2 || jc < 1 || jc >= H - 2) {
        *gx_out = gx_bg;
        *gy_out = gy_bg;
        return;
    }
    const float alpha = xn - (float) ic;
    const float beta  = yn - (float) jc;
    float gx_pert = 0.0f, gy_pert = 0.0f;
    for (int dq = 0; dq < 2; dq++) {
        const int jj  = jc + dq;
        const int row = jj * W;
        for (int dp = 0; dp < 2; dp++) {
            const int ii = ic + dp;
            float vgx = (pert[row + ii + 1] - pert[row + ii - 1]) * inv_2dx;
            float vgy = (pert[row + W + ii] - pert[row - W + ii]) * inv_2dx;
            if (bg_corner) {
                vgx += (bg_corner[row + ii + 1] - bg_corner[row + ii - 1]) * inv_2dx;
                vgy += (bg_corner[row + W + ii] - bg_corner[row - W + ii]) * inv_2dx;
            }
            const float wx = (dp == 0) ? (1.0f - alpha) : alpha;
            const float wy = (dq == 0) ? (1.0f - beta)  : beta;
            const float w = wx * wy;
            gx_pert += w * vgx;
            gy_pert += w * vgy;
        }
    }
    *gx_out = gx_pert + gx_bg;
    *gy_out = gy_pert + gy_bg;
}

/* Time derivative of the GEM vector potential A_g at the particle, for
 * the GM inductive piece -m d_t A_g of the Tier-3 gravity force.
 *
 * Half-step A_g convention (v36): the A_g buffers store half-step time
 * slices, so .prev = A^{n-1/2}, .curr = A^{n+1/2} after the leapfrog +
 * rotation.  Then the 1-step centered difference
 *     d_t A_g^n = (curr - prev) / dt = (A^{n+1/2} - A^{n-1/2}) / dt
 * is centered cleanly at t^n with a tight 1-step stencil.  The source
 * J^{n-1/2} fed into the wave equation now matches A's half-step time
 * naturally, dissolving the source-time staggering issue.  Both pieces
 * resolve in one architectural change. */
static void dt_A_g_at_total(const struct gr_sim* sim, float x, float y,
                            float* dAx_out, float* dAy_out) {
    *dAx_out = 0.0f;
    *dAy_out = 0.0f;
    if (!sim || !sim->gravitomagnetic_inductive_enabled) return;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    const float inv_dt = 1.0f / sim->dt;

    const float Ax_curr = cic_interp_xedge(sim->fields[GR_FIELD_A_GX].curr,
                                           W, H, dx, x, y);
    const float Ax_prev = cic_interp_xedge(sim->fields[GR_FIELD_A_GX].prev,
                                           W, H, dx, x, y);
    const float Ay_curr = cic_interp_yedge(sim->fields[GR_FIELD_A_GY].curr,
                                           W, H, dx, x, y);
    const float Ay_prev = cic_interp_yedge(sim->fields[GR_FIELD_A_GY].prev,
                                           W, H, dx, x, y);
    const float s = sim->grav_inductive_sign;
    *dAx_out = s * (Ax_curr - Ax_prev) * inv_dt;
    *dAy_out = s * (Ay_curr - Ay_prev) * inv_dt;
}

/* Time derivative of the EM vector potential A_em at the particle, for
 * the inductive piece -q d_t A of the EM Lorentz force.
 *
 * The leapfrog rotation in gr_sim_step cycles the three time-buffers per
 * field.  After each step's rotation:
 *     fields[*].prev = A^n        (matched to particle's t = n dt)
 *     fields[*].curr = A^{n+1}
 *     fields[*].next = A^{n-1}    (recycled to hold next leapfrog's output)
 * So the CENTERED difference at t^n is
 *     d_t A^n = (A^{n+1} - A^{n-1}) / (2 dt) = (curr - next) / (2 dt),
 * a 2nd-order accurate temporal derivative.  Spatial interpolation uses
 * the sublattice-aware kernels (X_EDGE for A_x, Y_EDGE for A_y) -- the
 * SAME interpolators that already serve the magnetic curl path.
 *
 * Background A:  for the static UNIFORM_MAGNETIC background and any
 * other current analytic generator, A^{bg} is time-independent, so its
 * contribution to d_t A vanishes by construction (sampled-bg arrays are
 * never written by the leapfrog).  This routine therefore reads only the
 * PERTURBATION arrays.  If a future time-varying analytic background is
 * added, this is the place to add an extra term.
 *
 * Gated by em_lorentz_force_enabled.  When field_evolution_enabled is 0
 * and the buffers haven't been manipulated by hand, all three buffers are
 * zero so this evaluator returns zero -- bit-exact baseline.  Tests can
 * manually pre-populate curr and next to drive a known d_t A. */
static void dt_A_em_at_total(const struct gr_sim* sim, float x, float y,
                             float* dAx_out, float* dAy_out) {
    *dAx_out = 0.0f;
    *dAy_out = 0.0f;
    if (!sim || !sim->em_lorentz_force_enabled) return;
    /* Fine-grained gate: skip the centered-time-difference if the
     * inductive piece is disabled (diagnostic isolation for Stage 27). */
    if (!sim->em_inductive_enabled) return;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float s = sim->em_inductive_sign;
    /* Half-step A_em convention (v36): .prev = A^{n-1/2}, .curr = A^{n+1/2}.
     * d_t A^n = (curr - prev) / dt = 1-step centered diff at t^n.  The
     * BACKWARD discretization option is now structurally equivalent and
     * retained for API compat (legacy stages that didn't drive A.next).
     * CENTERED (default) reads the natural 1-step diff. */
    const float inv_dt = 1.0f / sim->dt;
    const float Ax_curr = cic_interp_xedge(sim->fields[GR_FIELD_A_X].curr,
                                           W, H, dx, x, y);
    const float Ax_prev = cic_interp_xedge(sim->fields[GR_FIELD_A_X].prev,
                                           W, H, dx, x, y);
    const float Ay_curr = cic_interp_yedge(sim->fields[GR_FIELD_A_Y].curr,
                                           W, H, dx, x, y);
    const float Ay_prev = cic_interp_yedge(sim->fields[GR_FIELD_A_Y].prev,
                                           W, H, dx, x, y);
    *dAx_out = s * (Ax_curr - Ax_prev) * inv_dt;
    *dAy_out = s * (Ay_curr - Ay_prev) * inv_dt;
}

/* Total EM magnetic field B_z = (curl A)_z at the particle position.
 * Structurally identical to B_g_z_at_total but points at the EM arrays
 * (Ax_bg / Ay_bg for the background, fields[A_X] / fields[A_Y] for the
 * perturbation).  The curl evaluator is shared (curl_Agz_sampled_at is
 * sublattice-aware via the X_EDGE / Y_EDGE conventions and applies to any
 * pair of edge-staggered components).  Gated by em_lorentz_force_enabled. */
static float B_em_z_at_total(const struct gr_sim* sim, float x, float y) {
    if (!sim || !sim->em_lorentz_force_enabled) return 0.0f;
    float bg = 0.0f;
    if (sim->bg_mode == GR_BG_MODE_ANALYTIC) {
        gr_bg_eval_B_em(sim, x, y, &bg);
    } else {
        bg = curl_Agz_sampled_at(sim->Ax_bg, sim->Ay_bg,
                                 sim->width, sim->height, sim->dx, x, y);
    }
    /* Half-step A convention: same averaging as B_g_z_at_total. */
    float pert = 0.0f;
    if (sim->field_evolution_enabled) {
        const float c1 = curl_Agz_sampled_at(sim->fields[GR_FIELD_A_X].prev,
                                             sim->fields[GR_FIELD_A_Y].prev,
                                             sim->width, sim->height, sim->dx, x, y);
        const float c2 = curl_Agz_sampled_at(sim->fields[GR_FIELD_A_X].curr,
                                             sim->fields[GR_FIELD_A_Y].curr,
                                             sim->width, sim->height, sim->dx, x, y);
        pert = 0.5f * (c1 + c2);
    } else {
        pert = curl_Agz_sampled_at(sim->fields[GR_FIELD_A_X].curr,
                                   sim->fields[GR_FIELD_A_Y].curr,
                                   sim->width, sim->height, sim->dx, x, y);
    }
    return bg + pert;
}

/* Boris-leapfrog kick-drift for one timestep.
 *
 * gr_sandbox_v32.tex §9.2:
 *   p_{n+1/2} = p_{n-1/2} + F^n * dt
 *   v_{n+1/2} = p_{n+1/2} / sqrt(m^2 + |p|^2/c^2)
 *   x_{n+1}   = x_n + v_{n+1/2} * dt
 *
 * For Tier-0 (Newtonian) the force depends only on position, so plain
 * kick-drift is fine. For Tier-2 the force depends on v as well, and we
 * approximate v^n by the lagged half-step velocity from p^{n-1/2}. This is
 * 2nd-order in dt and preserves the symmetric kick-drift-kick structure when
 * the corrector below is used. */
void gr_particle_push_all(struct gr_sim* sim) {
    if (!sim || sim->n_particles <= 0) return;
    const float dt = sim->dt;
    const float c2 = sim->c_eff * sim->c_eff;
    for (int i = 0; i < sim->n_particles; i++) {
        gr_particle_t* p = &sim->particles[i];
        float grad_x, grad_y;
        grav_grad_at(sim, p->x, p->y, &grad_x, &grad_y);
        const float phi = gr_phi_g_total_at(sim, p->x, p->y);

        /* Total gravitomagnetic field B_g_z at the particle (Tier-1 source
         * for the v x B_g term in grav_force_at).  Combines background
         * (analytic OR sampled, depending on bg_mode) and perturbation
         * (always sampled, via Yee curl).  Gated by the
         * gravitomagnetic_force_enabled flag so that stage09 (clock-only
         * isolation) can disable the v x B_g piece. */
        const float Bg_z = B_g_z_at_total(sim, p->x, p->y);
        /* EM magnetic field B_em_z at the particle for the q v x B Lorentz
         * piece (Stage 23+).  Same combined bg+perturbation evaluator as
         * B_g_z but pointing at the EM arrays.  Gated by em_lorentz_force_
         * enabled.  When the particle is neutral or no EM A is present
         * this contributes nothing. */
        const float B_em_z = B_em_z_at_total(sim, p->x, p->y);
        /* EM scalar gradient at the particle for the -q grad phi
         * electrostatic piece (Stage 24+).  Combines bg + perturbation;
         * gated by em_lorentz_force_enabled. */
        float E_phi_grad_x, E_phi_grad_y;
        phi_em_grad_at_total(sim, p->x, p->y, &E_phi_grad_x, &E_phi_grad_y);
        /* EM vector-potential time derivative at the particle for the
         * -q d_t A inductive piece (Stage 25+).  Centered difference at
         * t^n using fields[A_X/A_Y].curr (=A^{n+1}) and .next (=A^{n-1})
         * after rotation.  Returns 0 when buffers are quiescent. */
        float E_dAx, E_dAy;
        dt_A_em_at_total(sim, p->x, p->y, &E_dAx, &E_dAy);
        /* GM vector-potential time derivative for the analogous gravity
         * inductive piece -m d_t A_g (Stage 28+).  Default OFF; gated by
         * gravitomagnetic_inductive_enabled. */
        float G_dAx, G_dAy;
        dt_A_g_at_total(sim, p->x, p->y, &G_dAx, &G_dAy);

        /* v from lagged p^{n-1/2} — used to evaluate the velocity-dependent
         * Tier-2 force terms.  Negligible cost for Newtonian since grav_force_at
         * ignores v in that mode. */
        float pmag2 = p->px * p->px + p->py * p->py;
        float gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
        const float vx_pre = p->px / (gamma * p->mass);
        const float vy_pre = p->py / (gamma * p->mass);

        float Fx, Fy;
        grav_force_at(sim, p->mass, vx_pre, vy_pre, phi, grad_x, grad_y, Bg_z, G_dAx, G_dAy, &Fx, &Fy);
        /* Additive EM Lorentz contribution: -q grad phi - q d_t A + q v x B
         * (Stages 23/24/25; full static-and-dynamic EM Lorentz force). */
        {
            float Fx_em, Fy_em;
            em_force_at(p->charge, vx_pre, vy_pre,
                        E_phi_grad_x, E_phi_grad_y,
                        E_dAx, E_dAy,
                        B_em_z, &Fx_em, &Fy_em);
            Fx += Fx_em;
            Fy += Fy_em;
        }

        /* Corrector iteration for velocity-dependent terms (Tier-2 scalar
         * coupling + Tier-1 gravitomagnetic v x B_g + EM v x B_em):
         * midpoint velocity gives 2nd-order time accuracy.  Engage it
         * whenever the force depends on v — RELATIVISTIC tier always,
         * plus any tier when Bg_z != 0 or B_em_z != 0 (and the particle
         * carries charge for the EM piece).  The -q grad phi piece is
         * v-independent so it doesn't drive the corrector. */
        const int em_velocity_dep = (B_em_z != 0.0f) && (p->charge != 0.0f);
        if (sim->force_tier == GR_FORCE_RELATIVISTIC
            || Bg_z != 0.0f || em_velocity_dep) {
            const float px_pred = p->px + Fx * dt;
            const float py_pred = p->py + Fy * dt;
            const float pmag2_pred = px_pred * px_pred + py_pred * py_pred;
            const float gamma_pred = sqrtf(1.0f + pmag2_pred / (p->mass * p->mass * c2));
            const float vx_post = px_pred / (gamma_pred * p->mass);
            const float vy_post = py_pred / (gamma_pred * p->mass);
            const float vx_mid = 0.5f * (vx_pre + vx_post);
            const float vy_mid = 0.5f * (vy_pre + vy_post);
            grav_force_at(sim, p->mass, vx_mid, vy_mid, phi, grad_x, grad_y, Bg_z, G_dAx, G_dAy, &Fx, &Fy);
            float Fx_em, Fy_em;
            em_force_at(p->charge, vx_mid, vy_mid,
                        E_phi_grad_x, E_phi_grad_y,
                        E_dAx, E_dAy,
                        B_em_z, &Fx_em, &Fy_em);
            Fx += Fx_em;
            Fy += Fy_em;
        }

        p->px += Fx * dt;
        p->py += Fy * dt;
        pmag2 = p->px * p->px + p->py * p->py;
        gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
        const float vx = p->px / (gamma * p->mass);
        const float vy = p->py / (gamma * p->mass);

        /* Stage 9: accumulate proper time over [t_n, t_{n+1}] using v^{n+1/2}
         * (just computed) and the (Phi, A_g) at the particle's current
         * position x_n.  This is 2nd-order accurate since v^{n+1/2} sits at
         * the midpoint of the drift interval.  Formula (v34 Eq. 75 with
         * v34 sign + factor corrections):
         *   d_tau = dt * sqrt(1 + 2 Phi/c^2 - (1 - 2 Phi/c^2) v^2/c^2
         *                       - 8 (v . A_g)/c^2)
         * A_g is zero unless a SPINNING_POINT_MASS background is installed. */
        {
            const float inv_c2 = 1.0f / c2;
            float Agx = 0.0f, Agy = 0.0f;
            gr_bg_eval_A_g(sim, p->x, p->y, &Agx, &Agy);
            const float v2      = vx * vx + vy * vy;
            const float two_phi = 2.0f * phi * inv_c2;
            const float v_dot_A = vx * Agx + vy * Agy;
            const float radicand = 1.0f + two_phi
                                 - (1.0f - two_phi) * v2 * inv_c2
                                 - 8.0f * v_dot_A * inv_c2;
            const float dtau = (radicand > 0.0f) ? dt * sqrtf(radicand) : 0.0f;
            p->proper_time += dtau;
        }

        p->x += vx * dt;
        p->y += vy * dt;
    }
}

/* ----------------------------------------------------------------------------
 * Particle API
 * --------------------------------------------------------------------------*/

int gr_sim_add_particle(gr_sim_t* sim, float x, float y,
                        float mass, float charge,
                        float vx, float vy) {
    if (!sim || mass <= 0.0f) return -1;
    if (sim->n_particles >= sim->particles_capacity) {
        int new_cap = sim->particles_capacity > 0 ? sim->particles_capacity * 2 : 16;
        gr_particle_t* new_arr = (gr_particle_t*) realloc(sim->particles,
                                                          (size_t) new_cap * sizeof(*new_arr));
        if (!new_arr) return -1;
        sim->particles          = new_arr;
        sim->particles_capacity = new_cap;
    }
    /* Initialize half-step-back momentum from a single half-kick of the
     * force at t=0 (gr_sandbox_v32.tex §9.7 "leapfrog initialization"). The
     * leapfrog needs p^{-1/2} to advance to p^{+1/2} on the first step. */
    const float c2    = sim->c_eff * sim->c_eff;
    /* Compute force at the initial position, using the current tier setting
     * so the half-step-back kick matches what the pusher will apply on step 0. */
    float grad_x, grad_y;
    grav_grad_at(sim, x, y, &grad_x, &grad_y);
    const float phi_init = gr_phi_g_total_at(sim, x, y);
    /* Tier-1 gravitomagnetic B_g_z + EM (B_z, grad phi, d_t A) at the
     * initial position.  Uses the same total evaluators as
     * gr_particle_push_all so the half-step-back kick stays consistent
     * with what the pusher will apply on step 0. */
    const float Bg_z_init   = B_g_z_at_total(sim, x, y);
    const float B_em_z_init = B_em_z_at_total(sim, x, y);
    float Egx_init, Egy_init;
    phi_em_grad_at_total(sim, x, y, &Egx_init, &Egy_init);
    float dAx_init, dAy_init;
    dt_A_em_at_total(sim, x, y, &dAx_init, &dAy_init);
    float G_dAx_init, G_dAy_init;
    dt_A_g_at_total(sim, x, y, &G_dAx_init, &G_dAy_init);
    float Fx, Fy;
    grav_force_at(sim, mass, vx, vy, phi_init, grad_x, grad_y, Bg_z_init,
                  G_dAx_init, G_dAy_init, &Fx, &Fy);
    {
        float Fx_em, Fy_em;
        em_force_at(charge, vx, vy, Egx_init, Egy_init,
                    dAx_init, dAy_init, B_em_z_init,
                    &Fx_em, &Fy_em);
        Fx += Fx_em;
        Fy += Fy_em;
    }
    /* p^0 = gamma_0 m v (relativistic). */
    const float gamma0 = 1.0f / sqrtf(fmaxf(1.0f - (vx * vx + vy * vy) / c2, 1e-12f));
    const float px0    = gamma0 * mass * vx;
    const float py0    = gamma0 * mass * vy;
    /* p^{-1/2} = p^0 - F^0 * dt/2 */
    const float half_dt = 0.5f * sim->dt;
    const int   idx = sim->n_particles++;
    sim->particles[idx].x           = x;
    sim->particles[idx].y           = y;
    sim->particles[idx].px          = px0 - Fx * half_dt;
    sim->particles[idx].py          = py0 - Fy * half_dt;
    sim->particles[idx].mass        = mass;
    sim->particles[idx].charge      = charge;
    sim->particles[idx].proper_time = 0.0f;
    return idx;
}

int gr_sim_particle_count(const gr_sim_t* sim) {
    return sim ? sim->n_particles : 0;
}

const gr_particle_t* gr_sim_get_particle(const gr_sim_t* sim, int idx) {
    if (!sim || idx < 0 || idx >= sim->n_particles) return NULL;
    return &sim->particles[idx];
}

void gr_sim_clear_particles(gr_sim_t* sim) {
    if (!sim) return;
    sim->n_particles = 0;
    /* keep capacity allocated for reuse */
}

float gr_sim_particle_energy(const gr_sim_t* sim, int idx) {
    if (!sim || idx < 0 || idx >= sim->n_particles) return 0.0f;
    const gr_particle_t* p = &sim->particles[idx];
    const float c2     = sim->c_eff * sim->c_eff;
    const float pmag2  = p->px * p->px + p->py * p->py;
    const float gamma  = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
    const float Etot   = gamma * p->mass * c2 + p->mass * gr_phi_g_total_at(sim, p->x, p->y);
    return Etot;
}
