/* Particle integrator — relativistic Boris-leapfrog (kick-drift) in 2D.
 * Spec reference: gr_sandbox_v32.tex §9.2 (eq:leapfrog_field neighborhood),
 * §9.5 (CIC adjoint condition for force interpolation), and v33 §12.7. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>

/* CIC interpolation of a cell-centered scalar field at (x, y).
 * Same sub-cell math as gr_cic_deposit_scalar (eq:cic_deposit) — bilinear
 * weights satisfy the §9.5 adjoint condition automatically. */
static float cic_interpolate(const float* f, int W, int H, float dx, float x, float y) {
    if (!f) return 0.0f;
    const float xn = x / dx - 0.5f;
    const float yn = y / dx - 0.5f;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    if (ic < 0 || ic >= W - 1 || jc < 0 || jc >= H - 1) return 0.0f;
    const float alpha = xn - (float) ic;
    const float beta  = yn - (float) jc;
    const float w00 = (1.0f - alpha) * (1.0f - beta);
    const float w10 =         alpha  * (1.0f - beta);
    const float w01 = (1.0f - alpha) *         beta;
    const float w11 =         alpha  *         beta;
    return w00 * f[jc * W + ic]
         + w10 * f[jc * W + ic + 1]
         + w01 * f[(jc + 1) * W + ic]
         + w11 * f[(jc + 1) * W + ic + 1];
}

float gr_phi_g_total_at(const struct gr_sim* sim, float x, float y) {
    if (!sim) return 0.0f;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    const float pert = cic_interpolate(sim->fields[GR_FIELD_PHI_GRAV].curr, W, H, dx, x, y);
    const float bg   = cic_interpolate(sim->phi_g_bg,                       W, H, dx, x, y);
    return pert + bg;
}

/* CIC-interpolated centered-FD gradient of Phi_g_total at (x, y). The CIC
 * weights match those used for deposit (W_2 / bilinear), satisfying the
 * §9.5 adjoint pairing so that a stationary self-deposit produces no
 * self-force on the particle (relevant once Stage 10 adds source coupling;
 * for Stage 7 in a fixed background it's just the natural force law). */
static void grav_grad_at(const struct gr_sim* sim, float x, float y,
                         float* gx_out, float* gy_out) {
    *gx_out = 0.0f;
    *gy_out = 0.0f;
    if (!sim) return;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    const float xn = x / dx - 0.5f;
    const float yn = y / dx - 0.5f;
    const int   ic = (int) floorf(xn);
    const int   jc = (int) floorf(yn);
    /* Need ic-1 and ic+2 in bounds for the centered FD at all four corners. */
    if (ic < 1 || ic >= W - 2 || jc < 1 || jc >= H - 2) return;
    const float alpha = xn - (float) ic;
    const float beta  = yn - (float) jc;
    const float inv_2dx = 1.0f / (2.0f * dx);

    const float* pert = sim->fields[GR_FIELD_PHI_GRAV].curr;
    const float* bg   = sim->phi_g_bg;

    /* Evaluate d_x Phi_g and d_y Phi_g at the four surrounding cells via
     * centered FD on (pert + bg), then bilinear-combine to the particle. */
    float gx[4], gy[4];
    const int idx[4] = {jc * W + ic,
                        jc * W + ic + 1,
                        (jc + 1) * W + ic,
                        (jc + 1) * W + ic + 1};
    for (int q = 0; q < 4; q++) {
        const int k = idx[q];
        float dphi_dx = (pert[k + 1] - pert[k - 1]) * inv_2dx;
        float dphi_dy = (pert[k + W] - pert[k - W]) * inv_2dx;
        if (bg) {
            dphi_dx += (bg[k + 1] - bg[k - 1]) * inv_2dx;
            dphi_dy += (bg[k + W] - bg[k - W]) * inv_2dx;
        }
        gx[q] = dphi_dx;
        gy[q] = dphi_dy;
    }
    const float w00 = (1.0f - alpha) * (1.0f - beta);
    const float w10 =         alpha  * (1.0f - beta);
    const float w01 = (1.0f - alpha) *         beta;
    const float w11 =         alpha  *         beta;
    *gx_out = w00 * gx[0] + w10 * gx[1] + w01 * gx[2] + w11 * gx[3];
    *gy_out = w00 * gy[0] + w10 * gy[1] + w01 * gy[2] + w11 * gy[3];
}

/* Boris-leapfrog kick-drift for one timestep.
 *
 * gr_sandbox_v32.tex §9.2:
 *   p_{n+1/2} = p_{n-1/2} + F^n * dt
 *   v_{n+1/2} = p_{n+1/2} / sqrt(m^2 + |p|^2/c^2)
 *   x_{n+1}   = x_n + v_{n+1/2} * dt
 *
 * For Stage 7 the force is purely the gravitational gradient,
 *   F = -m * grad(Phi_g_total) (eq:force_gem_pot with v=0 and A_g=0).
 * Later stages add Lorentz/gravitomagnetic terms via the Boris rotation. */
void gr_particle_push_all(struct gr_sim* sim) {
    if (!sim || sim->n_particles <= 0) return;
    const float dt    = sim->dt;
    const float c2    = sim->c_eff * sim->c_eff;
    for (int i = 0; i < sim->n_particles; i++) {
        gr_particle_t* p = &sim->particles[i];
        float gx, gy;
        grav_grad_at(sim, p->x, p->y, &gx, &gy);
        const float Fx = -p->mass * gx;
        const float Fy = -p->mass * gy;
        p->px += Fx * dt;
        p->py += Fy * dt;
        const float pmag2 = p->px * p->px + p->py * p->py;
        const float gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
        const float vx    = p->px / (gamma * p->mass);
        const float vy    = p->py / (gamma * p->mass);
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
    /* Compute force at the initial position. */
    float gx, gy;
    grav_grad_at(sim, x, y, &gx, &gy);
    const float Fx = -mass * gx;
    const float Fy = -mass * gy;
    /* p^0 = m * v (linear approximation; the relativistic correction at low
     * v is negligible and corrected by the integrator over many steps). */
    const float gamma0 = 1.0f / sqrtf(fmaxf(1.0f - (vx * vx + vy * vy) / c2, 1e-12f));
    const float px0    = gamma0 * mass * vx;
    const float py0    = gamma0 * mass * vy;
    /* p^{-1/2} = p^0 - F^0 * dt/2 */
    const float half_dt = 0.5f * sim->dt;
    const int   idx = sim->n_particles++;
    sim->particles[idx].x      = x;
    sim->particles[idx].y      = y;
    sim->particles[idx].px     = px0 - Fx * half_dt;
    sim->particles[idx].py     = py0 - Fy * half_dt;
    sim->particles[idx].mass   = mass;
    sim->particles[idx].charge = charge;
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
