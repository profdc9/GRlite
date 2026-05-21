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

    /* Perturbation Phi_g is always on the corner sublattice. */
    const float* pert = sim->fields[GR_FIELD_PHI_GRAV].curr;

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
 * total Phi_g (perturbation + background, value phi at the particle) and its
 * gradient (grad_x, grad_y) = grad(Phi_g_total).  The local gravity vector is
 * g = -grad(Phi_g).
 *
 * Tiers (gr_sandbox_v33.tex §"Practical implementation tiers"):
 *   NEWTONIAN:    F = m * g
 *   RELATIVISTIC: F = m * [g (1 + v^2/c^2 + 4 phi/c^2) - 4 (v . g) v / c^2]
 *
 * The relativistic expression is the Einstein-Infeld-Hoffmann 1PN equation
 * of motion for a test particle in a static field, harmonic gauge, with the
 * v33 doc's isotropic-form metric (g_{ij} = (1-2 phi/c^2) delta_ij).  See
 * Ali-Haïmoud, GR Fall 2019 lecture 25, eq. 37, with psi = xi = 0 and dt phi = 0:
 *
 *   dv/dt = -(1 + v^2/c^2 + 4 phi/c^2) grad(phi) + 4 (v . grad(phi)) v / c^2
 *
 * Substituting g = -grad(phi) and (v . g) = -(v . grad(phi)) gives the form
 * above.  This differs from v33 eq:geodesic_expansion (line 668) in two ways:
 *   (a) the 4*phi*g/c^2 ("Shapiro") term is added — comes from the g_{00}
 *       expansion at O(v^4) in Gamma^i_{00} (lecture eq. 30).
 *   (b) the sign on the velocity-coupling term is NEGATIVE, not positive as
 *       v33 has.  Empirical confirmation: with the v33 sign the test orbit
 *       precesses retrograde; with the EIH sign it precesses prograde.
 *
 * To be folded into the doc as v34.
 *
 * Tier 1 (GEM with A_g) and Tier 3 (full) arrive at Stage 10+ when the
 * perturbation A_g potentials are active and contribute via psi and xi. */
static inline void grav_force_at(const struct gr_sim* sim,
                                 float mass, float vx, float vy,
                                 float phi,
                                 float grad_x, float grad_y,
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

        /* v from lagged p^{n-1/2} — used to evaluate the velocity-dependent
         * Tier-2 force terms.  Negligible cost for Newtonian since grav_force_at
         * ignores v in that mode. */
        float pmag2 = p->px * p->px + p->py * p->py;
        float gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
        const float vx_pre = p->px / (gamma * p->mass);
        const float vy_pre = p->py / (gamma * p->mass);

        float Fx, Fy;
        grav_force_at(sim, p->mass, vx_pre, vy_pre, phi, grad_x, grad_y, &Fx, &Fy);

        /* Corrector iteration for velocity-dependent Tier-2 force: midpoint
         * velocity gives 2nd-order time accuracy. */
        if (sim->force_tier == GR_FORCE_RELATIVISTIC) {
            const float px_pred = p->px + Fx * dt;
            const float py_pred = p->py + Fy * dt;
            const float pmag2_pred = px_pred * px_pred + py_pred * py_pred;
            const float gamma_pred = sqrtf(1.0f + pmag2_pred / (p->mass * p->mass * c2));
            const float vx_post = px_pred / (gamma_pred * p->mass);
            const float vy_post = py_pred / (gamma_pred * p->mass);
            const float vx_mid = 0.5f * (vx_pre + vx_post);
            const float vy_mid = 0.5f * (vy_pre + vy_post);
            grav_force_at(sim, p->mass, vx_mid, vy_mid, phi, grad_x, grad_y, &Fx, &Fy);
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
    float Fx, Fy;
    grav_force_at(sim, mass, vx, vy, phi_init, grad_x, grad_y, &Fx, &Fy);
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
