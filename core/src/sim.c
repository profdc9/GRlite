/* Simulation lifecycle, source management, scenario dispatch.
 * Spec reference: gr_sandbox_v32.tex §9 + gr_sandbox_v33.tex §12.4. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void gr_sim_recompute_source_coeffs(struct gr_sim* sim) {
    if (!sim) return;
    const float inv_c2  = 1.0f / (sim->c_eff * sim->c_eff);
    /* Gravity: leapfrog Lap + sc*src form yields static limit
     *   Lap Phi_g = -sc * rho_m = +4 pi G * rho_m   (Newton Poisson)
     * Positive masses attract.  sc_grav = -4 pi G.
     *
     * EM (standard Maxwell, post-v36 sign fix):
     *   Lap phi_em - (1/c^2) d^2phi/dt^2 = -rho_q / epsilon_0
     *                                   = -4 pi k_e rho_q
     * which in leapfrog form gives static Lap phi_em = -4 pi k_e rho_q.
     * Yielding standard-EM behavior:
     *   like charges REPEL, opposite charges ATTRACT,
     *   parallel like-currents ATTRACT (Ampere).
     *
     * sc_em = +4 pi k_e   (OPPOSITE sign from gravity).
     *
     * The originally-derived EM code shared the gravity-style negative
     * sign (matching how the EM scaffolding was lifted from GEM); the
     * result was an "anti-Maxwell" PIC channel in which same-sign
     * charges attracted.  v36 corrects this; Stage 34 verifies the
     * full EM chain (Poisson sign, J=rho*v, A wave-eq sign, Lorentz
     * force direction). */
    const float c_grav  = -4.0f * (float) M_PI * sim->G_eff;
    const float c_em    = +4.0f * (float) M_PI * sim->k_e;
    sim->fields[GR_FIELD_PHI_GRAV].source_coeff = c_grav;
    sim->fields[GR_FIELD_A_GX    ].source_coeff = c_grav * inv_c2;
    sim->fields[GR_FIELD_A_GY    ].source_coeff = c_grav * inv_c2;
    sim->fields[GR_FIELD_PHI_EM  ].source_coeff = c_em;
    sim->fields[GR_FIELD_A_X     ].source_coeff = c_em   * inv_c2;
    sim->fields[GR_FIELD_A_Y     ].source_coeff = c_em   * inv_c2;
    /* v39: propagate to per-particle self-field sets. */
    if (sim->self_field_sets) {
        for (int i = 0; i < sim->particles_capacity; i++) {
            if (sim->self_field_sets[i])
                gr_self_field_set_bind_source_coeffs(sim, sim->self_field_sets[i]);
        }
    }
}

gr_sim_t* gr_sim_create(int width, int height, float dx, float c_eff, float cfl) {
    if (width <= 2 || height <= 2 || !(dx > 0.0f) || !(c_eff > 0.0f) || !(cfl > 0.0f)) {
        return NULL;
    }
    gr_sim_t* sim = (gr_sim_t*) calloc(1, sizeof(*sim));
    if (!sim) return NULL;

    sim->width  = width;
    sim->height = height;
    sim->dx     = dx;
    sim->c_eff  = c_eff;
    sim->cfl    = cfl;
    sim->G_eff                   = 1.0f;
    sim->k_e                     = 1.0f;
    sim->field_evolution_enabled = 1;
    /* Esirkepov current deposition on by default (v35 answer B1). */
    sim->esirkepov_enabled       = 1;
    /* Tier-1 gravitomagnetic Lorentz force on by default (Stage 20+). */
    sim->gravitomagnetic_force_enabled = 1;
    /* GM inductive piece OFF by default for backward compatibility
     * (Stage 28+).  Existing gravity tests pre-date this piece. */
    sim->gravitomagnetic_inductive_enabled = 0;
    /* EM Lorentz force on by default (Stage 23+); all sub-pieces on. */
    sim->em_lorentz_force_enabled      = 1;
    sim->em_inductive_enabled          = 1;
    sim->em_electrostatic_enabled      = 1;
    sim->em_magnetic_enabled           = 1;
    /* v40 parabolic Lorenz-gauge cleaning (gr_sandbox_v38.tex sec:gauge).
     * Default OFF: the doc claims cleaning is "purely cosmetic" but
     * empirically (2026-06-03 regression) it shifts orbital periods
     * by 6.8%% (stage10), causes CFL-sweep failures (stage29), and
     * flips field signs in PIC profiles (stage32b).  The discrete
     * scheme inherits gauge invariance only approximately, so
     * modifying phi (the cleaning operation) leaks into the gradient
     * force on charged particles.  When enabled the cleaning IS
     * effective at bounding the Lorenz-residual drift; turn on via
     * gr_sim_set_parabolic_gauge_cleaning_enabled(sim, 1) if that
     * matters more than dynamic fidelity. */
    sim->parabolic_gauge_cleaning_enabled = 0;
    sim->parabolic_gauge_R                = cfl * cfl * 0.25f;
    /* v41 EM stress-energy contribution to GEM sources: default OFF. */
    sim->em_stress_energy_enabled         = 0;
    sim->j_smooth_passes               = 0;
    sim->kernel_radius                 = 1.5f;
    /* J time-correction off by default (raw Esirkepov J^{n-1/2}). */
    sim->esirkepov_violations    = 0;
    sim->rho_smooth_passes       = 0;
    sim->shape_function          = GR_SHAPE_CIC;
    sim->periodic_bc             = 0;
    /* dt from CFL — gr_sandbox_v32.tex §9.2 eq:cfl. Not enforced to allow
     * the Stage 1 instability test (§12.1) to deliberately exceed the limit. */
    sim->dt = cfl * dx / c_eff;

    const size_t n = (size_t) width * (size_t) height;

    /* Allocate six source arrays (always present, zero-filled). */
    sim->rho_matter = (float*) calloc(n, sizeof(float));
    sim->J_mx       = (float*) calloc(n, sizeof(float));
    sim->J_my       = (float*) calloc(n, sizeof(float));
    sim->rho_q      = (float*) calloc(n, sizeof(float));
    sim->J_qx       = (float*) calloc(n, sizeof(float));
    sim->J_qy       = (float*) calloc(n, sizeof(float));
    if (!sim->rho_matter || !sim->J_mx || !sim->J_my
        || !sim->rho_q   || !sim->J_qx || !sim->J_qy) {
        gr_sim_destroy(sim);
        return NULL;
    }

    /* Allocate three time levels for each of six fields (18 arrays). */
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        sim->fields[f].prev = (float*) calloc(n, sizeof(float));
        sim->fields[f].curr = (float*) calloc(n, sizeof(float));
        sim->fields[f].next = (float*) calloc(n, sizeof(float));
        if (!sim->fields[f].prev || !sim->fields[f].curr || !sim->fields[f].next) {
            gr_sim_destroy(sim);
            return NULL;
        }
    }

    /* Bind each field to its source array. */
    sim->fields[GR_FIELD_PHI_GRAV].source = sim->rho_matter;
    sim->fields[GR_FIELD_A_GX    ].source = sim->J_mx;
    sim->fields[GR_FIELD_A_GY    ].source = sim->J_my;
    sim->fields[GR_FIELD_PHI_EM  ].source = sim->rho_q;
    sim->fields[GR_FIELD_A_X     ].source = sim->J_qx;
    sim->fields[GR_FIELD_A_Y     ].source = sim->J_qy;

    gr_sim_recompute_source_coeffs(sim);
    return sim;
}

/* v39 self-field set lifecycle ------------------------------------------- */

gr_self_field_set_t* gr_self_field_set_alloc(int W, int H) {
    if (W <= 0 || H <= 0) return NULL;
    const size_t n = (size_t) W * (size_t) H;
    gr_self_field_set_t* s = (gr_self_field_set_t*) calloc(1, sizeof(*s));
    if (!s) return NULL;
    /* 6 source arrays. */
    s->rho_matter = (float*) calloc(n, sizeof(float));
    s->J_mx       = (float*) calloc(n, sizeof(float));
    s->J_my       = (float*) calloc(n, sizeof(float));
    s->rho_q      = (float*) calloc(n, sizeof(float));
    s->J_qx       = (float*) calloc(n, sizeof(float));
    s->J_qy       = (float*) calloc(n, sizeof(float));
    /* 6 fields x 3 time slices. */
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        s->fields[f].prev = (float*) calloc(n, sizeof(float));
        s->fields[f].curr = (float*) calloc(n, sizeof(float));
        s->fields[f].next = (float*) calloc(n, sizeof(float));
    }
    /* Verify and bind sources. */
    int ok = (s->rho_matter && s->J_mx && s->J_my
              && s->rho_q && s->J_qx && s->J_qy);
    for (int f = 0; f < GR_FIELD_COUNT && ok; f++) {
        if (!s->fields[f].prev || !s->fields[f].curr || !s->fields[f].next) ok = 0;
    }
    if (!ok) {
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            free(s->fields[f].prev);
            free(s->fields[f].curr);
            free(s->fields[f].next);
        }
        free(s->rho_matter); free(s->J_mx); free(s->J_my);
        free(s->rho_q);      free(s->J_qx); free(s->J_qy);
        free(s);
        return NULL;
    }
    s->fields[GR_FIELD_PHI_GRAV].source = s->rho_matter;
    s->fields[GR_FIELD_A_GX    ].source = s->J_mx;
    s->fields[GR_FIELD_A_GY    ].source = s->J_my;
    s->fields[GR_FIELD_PHI_EM  ].source = s->rho_q;
    s->fields[GR_FIELD_A_X     ].source = s->J_qx;
    s->fields[GR_FIELD_A_Y     ].source = s->J_qy;
    return s;
}

void gr_self_field_set_free(gr_self_field_set_t* s) {
    if (!s) return;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        free(s->fields[f].prev);
        free(s->fields[f].curr);
        free(s->fields[f].next);
    }
    free(s->rho_matter); free(s->J_mx); free(s->J_my);
    free(s->rho_q);      free(s->J_qx); free(s->J_qy);
    free(s);
}

void gr_self_field_set_bind_source_coeffs(gr_sim_t* sim, gr_self_field_set_t* s) {
    if (!sim || !s) return;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        s->fields[f].source_coeff = sim->fields[f].source_coeff;
    }
}

/* v39 self-field public API ----------------------------------------------- */

static int ensure_self_field_tracking(gr_sim_t* sim) {
    if (sim->particles_capacity <= 0) return -1;
    const int cap = sim->particles_capacity;
    if (!sim->self_field_sets) {
        sim->self_field_sets = (gr_self_field_set_t**) calloc((size_t) cap, sizeof(*sim->self_field_sets));
        if (!sim->self_field_sets) return -1;
    }
    if (!sim->self_field_eps_x) {
        sim->self_field_eps_x = (float*) calloc((size_t) cap, sizeof(float));
        if (!sim->self_field_eps_x) return -1;
    }
    if (!sim->self_field_eps_y) {
        sim->self_field_eps_y = (float*) calloc((size_t) cap, sizeof(float));
        if (!sim->self_field_eps_y) return -1;
    }
    return 0;
}

int gr_sim_particle_enable_self_field(gr_sim_t* sim, int idx) {
    if (!sim) return -1;
    if (idx < 0 || idx >= sim->n_particles) return -1;
    if (ensure_self_field_tracking(sim) != 0) return -1;
    if (sim->self_field_sets[idx]) return 0;  /* already enabled */
    gr_self_field_set_t* s = gr_self_field_set_alloc(sim->width, sim->height);
    if (!s) return -1;
    gr_self_field_set_bind_source_coeffs(sim, s);
    sim->self_field_sets[idx] = s;
    return 0;
}

void gr_sim_particle_disable_self_field(gr_sim_t* sim, int idx) {
    if (!sim || idx < 0 || idx >= sim->particles_capacity) return;
    if (!sim->self_field_sets) return;
    gr_self_field_set_free(sim->self_field_sets[idx]);
    sim->self_field_sets[idx] = NULL;
    if (sim->self_field_eps_x) sim->self_field_eps_x[idx] = 0.0f;
    if (sim->self_field_eps_y) sim->self_field_eps_y[idx] = 0.0f;
}

int gr_sim_particle_has_self_field(const gr_sim_t* sim, int idx) {
    if (!sim || idx < 0 || idx >= sim->particles_capacity) return 0;
    if (!sim->self_field_sets) return 0;
    return sim->self_field_sets[idx] != NULL;
}

void gr_sim_particle_set_self_field_epsilon(gr_sim_t* sim, int idx,
                                              float eps_x, float eps_y) {
    if (!sim || idx < 0 || idx >= sim->particles_capacity) return;
    if (ensure_self_field_tracking(sim) != 0) return;
    sim->self_field_eps_x[idx] = eps_x;
    sim->self_field_eps_y[idx] = eps_y;
}

void gr_sim_particle_get_self_field_epsilon(const gr_sim_t* sim, int idx,
                                              float* eps_x_out, float* eps_y_out) {
    if (eps_x_out) *eps_x_out = 0.0f;
    if (eps_y_out) *eps_y_out = 0.0f;
    if (!sim || idx < 0 || idx >= sim->particles_capacity) return;
    if (sim->self_field_eps_x && eps_x_out) *eps_x_out = sim->self_field_eps_x[idx];
    if (sim->self_field_eps_y && eps_y_out) *eps_y_out = sim->self_field_eps_y[idx];
}

void gr_sim_destroy(gr_sim_t* sim) {
    if (!sim) return;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        free(sim->fields[f].prev);
        free(sim->fields[f].curr);
        free(sim->fields[f].next);
    }
    free(sim->rho_matter);
    free(sim->J_mx);
    free(sim->J_my);
    free(sim->rho_q);
    free(sim->J_qx);
    free(sim->J_qy);
    free(sim->J_mx_prev);
    free(sim->J_my_prev);
    free(sim->J_qx_prev);
    free(sim->J_qy_prev);
    free(sim->damping_d);
    free(sim->friction_d);
    free(sim->phi_g_bg);
    free(sim->Agx_bg);
    free(sim->Agy_bg);
    free(sim->phi_bg);
    free(sim->Ax_bg);
    free(sim->Ay_bg);
    free(sim->c_local2_corner);
    free(sim->c_local2_xedge);
    free(sim->c_local2_yedge);
    /* v39: free all per-particle self-field sets. */
    if (sim->self_field_sets) {
        for (int i = 0; i < sim->particles_capacity; i++) {
            gr_self_field_set_free(sim->self_field_sets[i]);
        }
        free(sim->self_field_sets);
    }
    free(sim->self_field_eps_x);
    free(sim->self_field_eps_y);
    free(sim->particles);
    free(sim);
}

void gr_sim_step(gr_sim_t* sim) {
    if (!sim) return;

    /* Stage 10: deposit every particle's sources onto the grid before the
     * field leapfrog reads them.  Velocity for the current density is taken
     * from p^{n-1/2}, matching the force-evaluation convention. */
    if (sim->particle_source_deposition && sim->n_particles > 0) {
        gr_sim_clear_sources(sim);
        const float c2 = sim->c_eff * sim->c_eff;
        const float dt = sim->dt;
        const int   W  = sim->width;
        const int   H  = sim->height;
        const float dx = sim->dx;
        for (int i = 0; i < sim->n_particles; i++) {
            const gr_particle_t* p = &sim->particles[i];
            const float pmag2 = p->px * p->px + p->py * p->py;
            const float gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
            const float vx    = p->px / (gamma * p->mass);
            const float vy    = p->py / (gamma * p->mass);

            /* v39: if this particle has a self-field, also clear its
             * sources and deposit there.  Same kernel, same Esirkepov path. */
            gr_self_field_set_t* self = (sim->self_field_sets) ? sim->self_field_sets[i] : NULL;
            if (self) {
                const size_t n_cells = (size_t) W * (size_t) H;
                memset(self->rho_matter, 0, n_cells * sizeof(float));
                memset(self->rho_q,      0, n_cells * sizeof(float));
                memset(self->J_mx,       0, n_cells * sizeof(float));
                memset(self->J_my,       0, n_cells * sizeof(float));
                memset(self->J_qx,       0, n_cells * sizeof(float));
                memset(self->J_qy,       0, n_cells * sizeof(float));
            }

            /* Deposit rho^n at the current particle position.  Use TSC
             * (3x3, smoother) if selected; otherwise CIC (2x2).  BUMP
             * is v38 Tier-1: same 3x3 footprint as TSC but C-infinity
             * kernel (sub-exponential Fourier decay). */
            if (sim->shape_function == GR_SHAPE_BUMP) {
                const float Rk = sim->kernel_radius;
                if (p->mass   != 0.0f) gr_bump_deposit_corner(sim->rho_matter, W, H, dx, p->x, p->y, p->mass,   Rk);
                if (p->charge != 0.0f) gr_bump_deposit_corner(sim->rho_q,      W, H, dx, p->x, p->y, p->charge, Rk);
                if (self) {
                    if (p->mass   != 0.0f) gr_bump_deposit_corner(self->rho_matter, W, H, dx, p->x, p->y, p->mass,   Rk);
                    if (p->charge != 0.0f) gr_bump_deposit_corner(self->rho_q,      W, H, dx, p->x, p->y, p->charge, Rk);
                }
            } else if (sim->shape_function == GR_SHAPE_TSC) {
                if (p->mass   != 0.0f) gr_tsc_deposit_corner(sim->rho_matter, W, H, dx, p->x, p->y, p->mass);
                if (p->charge != 0.0f) gr_tsc_deposit_corner(sim->rho_q,      W, H, dx, p->x, p->y, p->charge);
                if (self) {
                    if (p->mass   != 0.0f) gr_tsc_deposit_corner(self->rho_matter, W, H, dx, p->x, p->y, p->mass);
                    if (p->charge != 0.0f) gr_tsc_deposit_corner(self->rho_q,      W, H, dx, p->x, p->y, p->charge);
                }
            } else {
                if (p->mass   != 0.0f) gr_cic_deposit_corner(sim->rho_matter, W, H, dx, p->x, p->y, p->mass);
                if (p->charge != 0.0f) gr_cic_deposit_corner(sim->rho_q,      W, H, dx, p->x, p->y, p->charge);
                if (self) {
                    if (p->mass   != 0.0f) gr_cic_deposit_corner(self->rho_matter, W, H, dx, p->x, p->y, p->mass);
                    if (p->charge != 0.0f) gr_cic_deposit_corner(self->rho_q,      W, H, dx, p->x, p->y, p->charge);
                }
            }

            /* Deposit J^{n-1/2} for the trajectory x^{n-1} -> x^n.
             * x^{n-1} = x^n - v^{n-1/2} * dt is exact under leapfrog drift. */
            const float x0     = p->x - vx * dt;
            const float y0     = p->y - vy * dt;
            const float x1     = p->x;
            const float y1     = p->y;
            int violated = 0;
            if (sim->esirkepov_enabled) {
                const int use_bump = (sim->shape_function == GR_SHAPE_BUMP);
                const float Rk_j  = sim->kernel_radius;
                if (p->mass != 0.0f) {
                    const int ok = use_bump
                        ? gr_bump_esirkepov_deposit_jxy(sim->J_mx, sim->J_my, W, H, dx, dt,
                                                        x0, y0, x1, y1, p->mass, Rk_j)
                        : gr_esirkepov_deposit_jxy     (sim->J_mx, sim->J_my, W, H, dx, dt,
                                                        x0, y0, x1, y1, p->mass);
                    if (!ok) violated = 1;
                    if (self) {
                        if (use_bump) gr_bump_esirkepov_deposit_jxy(self->J_mx, self->J_my, W, H, dx, dt,
                                                                    x0, y0, x1, y1, p->mass, Rk_j);
                        else          gr_esirkepov_deposit_jxy     (self->J_mx, self->J_my, W, H, dx, dt,
                                                                    x0, y0, x1, y1, p->mass);
                    }
                }
                if (p->charge != 0.0f) {
                    const int ok = use_bump
                        ? gr_bump_esirkepov_deposit_jxy(sim->J_qx, sim->J_qy, W, H, dx, dt,
                                                        x0, y0, x1, y1, p->charge, Rk_j)
                        : gr_esirkepov_deposit_jxy     (sim->J_qx, sim->J_qy, W, H, dx, dt,
                                                        x0, y0, x1, y1, p->charge);
                    if (!ok) violated = 1;
                    if (self) {
                        if (use_bump) gr_bump_esirkepov_deposit_jxy(self->J_qx, self->J_qy, W, H, dx, dt,
                                                                    x0, y0, x1, y1, p->charge, Rk_j);
                        else          gr_esirkepov_deposit_jxy     (self->J_qx, self->J_qy, W, H, dx, dt,
                                                                    x0, y0, x1, y1, p->charge);
                    }
                }
            }
            if (!sim->esirkepov_enabled || violated) {
                /* Direct CIC fallback (also used when Esirkepov is opted out
                 * for regression testing).  Deposit J at the trajectory
                 * midpoint -- which for nonzero shift is (x^n + jdx/?). */
                const float xc = 0.5f * (x0 + x1);
                const float yc = 0.5f * (y0 + y1);
                if (p->mass   != 0.0f && vx != 0.0f) gr_cic_deposit_xedge(sim->J_mx, W, H, dx, xc, yc, p->mass   * vx);
                if (p->mass   != 0.0f && vy != 0.0f) gr_cic_deposit_yedge(sim->J_my, W, H, dx, xc, yc, p->mass   * vy);
                if (p->charge != 0.0f && vx != 0.0f) gr_cic_deposit_xedge(sim->J_qx, W, H, dx, xc, yc, p->charge * vx);
                if (p->charge != 0.0f && vy != 0.0f) gr_cic_deposit_yedge(sim->J_qy, W, H, dx, xc, yc, p->charge * vy);
                if (violated) sim->esirkepov_violations++;
            }

            /* v40 spin dipole curl deposit: a spinning particle contributes
             * a curl-of-dipole-density on top of the regular J.  Magnetic
             * moment mu = (g q / 2m) S for EM; gravitomagnetic moment
             * sigma = 2 S.  Same kernel as the rho/J deposit.  Inherently
             * divergence-free so rho continuity is preserved. */
            if (p->spin != 0.0f) {
                const float mu_em   = (p->mass > 0.0f)
                    ? (p->g_factor * p->charge / (2.0f * p->mass)) * p->spin
                    : 0.0f;
                const float sigma_g = 2.0f * p->spin;
                if (sim->shape_function == GR_SHAPE_BUMP) {
                    const float Rk = sim->kernel_radius;
                    if (mu_em   != 0.0f) gr_bump_curl_dipole_deposit_jxy(sim->J_qx, sim->J_qy, W, H, dx, p->x, p->y, mu_em,   Rk);
                    if (sigma_g != 0.0f) gr_bump_curl_dipole_deposit_jxy(sim->J_mx, sim->J_my, W, H, dx, p->x, p->y, sigma_g, Rk);
                    if (self) {
                        if (mu_em   != 0.0f) gr_bump_curl_dipole_deposit_jxy(self->J_qx, self->J_qy, W, H, dx, p->x, p->y, mu_em,   Rk);
                        if (sigma_g != 0.0f) gr_bump_curl_dipole_deposit_jxy(self->J_mx, self->J_my, W, H, dx, p->x, p->y, sigma_g, Rk);
                    }
                } else if (sim->shape_function == GR_SHAPE_TSC) {
                    if (mu_em   != 0.0f) gr_tsc_curl_dipole_deposit_jxy(sim->J_qx, sim->J_qy, W, H, dx, p->x, p->y, mu_em);
                    if (sigma_g != 0.0f) gr_tsc_curl_dipole_deposit_jxy(sim->J_mx, sim->J_my, W, H, dx, p->x, p->y, sigma_g);
                    if (self) {
                        if (mu_em   != 0.0f) gr_tsc_curl_dipole_deposit_jxy(self->J_qx, self->J_qy, W, H, dx, p->x, p->y, mu_em);
                        if (sigma_g != 0.0f) gr_tsc_curl_dipole_deposit_jxy(self->J_mx, self->J_my, W, H, dx, p->x, p->y, sigma_g);
                    }
                } else {
                    if (mu_em   != 0.0f) gr_cic_curl_dipole_deposit_jxy(sim->J_qx, sim->J_qy, W, H, dx, p->x, p->y, mu_em);
                    if (sigma_g != 0.0f) gr_cic_curl_dipole_deposit_jxy(sim->J_mx, sim->J_my, W, H, dx, p->x, p->y, sigma_g);
                    if (self) {
                        if (mu_em   != 0.0f) gr_cic_curl_dipole_deposit_jxy(self->J_qx, self->J_qy, W, H, dx, p->x, p->y, mu_em);
                        if (sigma_g != 0.0f) gr_cic_curl_dipole_deposit_jxy(self->J_mx, self->J_my, W, H, dx, p->x, p->y, sigma_g);
                    }
                }
            }
        }
        /* Binomial smoothing on rho_matter (and rho_q) if enabled.
         * Applied N_smooth times, where each pass is the 3x3 [[1,2,1],
         * [2,4,2],[1,2,1]]/16 stencil.  Reduces high-spatial-frequency
         * aliasing in the deposited rho that the wave-equation leapfrog
         * would otherwise turn into a moving-particle self-force wake.
         * This is the canonical PIC noise-reduction technique. */
        if (sim->rho_smooth_passes > 0) {
            const size_t n = (size_t) W * (size_t) H;
            float* scratch = (float*) malloc(n * sizeof(float));
            if (scratch) {
                for (int pass = 0; pass < sim->rho_smooth_passes; pass++) {
                    /* rho_matter */
                    memcpy(scratch, sim->rho_matter, n * sizeof(float));
                    for (int j = 1; j < H - 1; j++) {
                        for (int i = 1; i < W - 1; i++) {
                            const int k = j * W + i;
                            sim->rho_matter[k] = (1.0f / 16.0f) * (
                                  1.0f * scratch[k - W - 1] + 2.0f * scratch[k - W] + 1.0f * scratch[k - W + 1]
                                + 2.0f * scratch[k - 1]     + 4.0f * scratch[k]     + 2.0f * scratch[k + 1]
                                + 1.0f * scratch[k + W - 1] + 2.0f * scratch[k + W] + 1.0f * scratch[k + W + 1]);
                        }
                    }
                    /* rho_q */
                    memcpy(scratch, sim->rho_q, n * sizeof(float));
                    for (int j = 1; j < H - 1; j++) {
                        for (int i = 1; i < W - 1; i++) {
                            const int k = j * W + i;
                            sim->rho_q[k] = (1.0f / 16.0f) * (
                                  1.0f * scratch[k - W - 1] + 2.0f * scratch[k - W] + 1.0f * scratch[k - W + 1]
                                + 2.0f * scratch[k - 1]     + 4.0f * scratch[k]     + 2.0f * scratch[k + 1]
                                + 1.0f * scratch[k + W - 1] + 2.0f * scratch[k + W] + 1.0f * scratch[k + W + 1]);
                        }
                    }
                }
                free(scratch);
            }
        }
        /* Optional binomial smoothing on the deposited J arrays (Stage 37
         * diagnostic).  Same 3x3 [[1,2,1],[2,4,2],[1,2,1]]/16 stencil
         * applied to the X_EDGE / Y_EDGE arrays in turn.  This damps
         * high-k content of the J deposit so the A wake propagated by
         * the wave equation is also damped -- a test of whether the
         * `-q d_t A` artifact on uniform motion is dispersion-driven. */
        if (sim->j_smooth_passes > 0) {
            const size_t n = (size_t) W * (size_t) H;
            float* scratch = (float*) malloc(n * sizeof(float));
            if (scratch) {
                float* targets[4] = { sim->J_mx, sim->J_my, sim->J_qx, sim->J_qy };
                for (int pass = 0; pass < sim->j_smooth_passes; pass++) {
                    for (int t = 0; t < 4; t++) {
                        memcpy(scratch, targets[t], n * sizeof(float));
                        for (int j = 1; j < H - 1; j++) {
                            for (int i = 1; i < W - 1; i++) {
                                const int k = j * W + i;
                                targets[t][k] = (1.0f / 16.0f) * (
                                      1.0f * scratch[k - W - 1] + 2.0f * scratch[k - W] + 1.0f * scratch[k - W + 1]
                                    + 2.0f * scratch[k - 1]     + 4.0f * scratch[k]     + 2.0f * scratch[k + 1]
                                    + 1.0f * scratch[k + W - 1] + 2.0f * scratch[k + W] + 1.0f * scratch[k + W + 1]);
                            }
                        }
                    }
                }
                free(scratch);
            }
        }

        /* v39: apply the SAME smoothing passes to each opted-in particle's
         * self-field sources so the self-gather mirrors what the collective
         * gather receives from this particle. */
        if (sim->self_field_sets && (sim->rho_smooth_passes > 0 || sim->j_smooth_passes > 0)) {
            const size_t n = (size_t) W * (size_t) H;
            float* scratch = (float*) malloc(n * sizeof(float));
            if (scratch) {
                for (int pi = 0; pi < sim->n_particles; pi++) {
                    gr_self_field_set_t* sf = sim->self_field_sets[pi];
                    if (!sf) continue;
                    if (sim->rho_smooth_passes > 0) {
                        float* rho_targets[2] = { sf->rho_matter, sf->rho_q };
                        for (int pass = 0; pass < sim->rho_smooth_passes; pass++) {
                            for (int t = 0; t < 2; t++) {
                                memcpy(scratch, rho_targets[t], n * sizeof(float));
                                for (int j = 1; j < H - 1; j++) {
                                    for (int i = 1; i < W - 1; i++) {
                                        const int k = j * W + i;
                                        rho_targets[t][k] = (1.0f / 16.0f) * (
                                              1.0f * scratch[k - W - 1] + 2.0f * scratch[k - W] + 1.0f * scratch[k - W + 1]
                                            + 2.0f * scratch[k - 1]     + 4.0f * scratch[k]     + 2.0f * scratch[k + 1]
                                            + 1.0f * scratch[k + W - 1] + 2.0f * scratch[k + W] + 1.0f * scratch[k + W + 1]);
                                    }
                                }
                            }
                        }
                    }
                    if (sim->j_smooth_passes > 0) {
                        float* j_targets[4] = { sf->J_mx, sf->J_my, sf->J_qx, sf->J_qy };
                        for (int pass = 0; pass < sim->j_smooth_passes; pass++) {
                            for (int t = 0; t < 4; t++) {
                                memcpy(scratch, j_targets[t], n * sizeof(float));
                                for (int j = 1; j < H - 1; j++) {
                                    for (int i = 1; i < W - 1; i++) {
                                        const int k = j * W + i;
                                        j_targets[t][k] = (1.0f / 16.0f) * (
                                              1.0f * scratch[k - W - 1] + 2.0f * scratch[k - W] + 1.0f * scratch[k - W + 1]
                                            + 2.0f * scratch[k - 1]     + 4.0f * scratch[k]     + 2.0f * scratch[k + 1]
                                            + 1.0f * scratch[k + W - 1] + 2.0f * scratch[k + W] + 1.0f * scratch[k + W + 1]);
                                    }
                                }
                            }
                        }
                    }
                }
                free(scratch);
            }
        }
    }

    /* v41 EM stress-energy contribution to gravity sources.  Added on
     * top of the particle-derived rho_matter / J_m before the GEM
     * wave-equation leapfrog so that the EM field energy and Poynting
     * flux gravitate via the standard linearized-GR coupling. */
    gr_sim_apply_em_stress_energy_sources(sim);

    /* Dynamic Shapiro: refresh c_local^2 from the CURRENT total Phi_g (bg +
     * perturbation) so a moving/deposited mass lenses the EM wave this step. */
    if (sim->em_shapiro_enabled && sim->em_shapiro_dynamic) {
        gr_em_shapiro_recompute_c_local2(sim);
    }

    if (sim->field_evolution_enabled) {
        gr_field_leapfrog_step_all(sim);

        /* v42 derivative-friction: CENTERED-IMPLICIT discretization
         * applied AFTER the standard leapfrog computes next, BEFORE
         * the rotation.
         *
         * Continuous PDE:  d^2 Phi/dt^2 + 2 gamma d Phi/dt = c^2 Lap Phi + ...
         * Centered discrete (gives ~implicit pull on Phi_next):
         *   (1+beta) Phi_next = 2 Phi_curr - (1-beta) Phi_prev + dt^2(...)
         * with beta = gamma * dt.  Equivalent to the post-leapfrog
         * correction:
         *   next[k] := (next_standard[k] + beta[k] * prev[k]) / (1 + beta[k])
         *
         * Eigenvalues of the recurrence at the 2D Nyquist mode (s^2=4)
         * are -1 (marginal) and -(1-beta)/(1+beta) (stable), so unlike
         * a post-step bleed on prev (which destabilizes the Nyquist
         * mode at CFL = 1/sqrt(2)), this form is stable for any beta. */
        const size_t ncells_fr = (size_t) sim->width * (size_t) sim->height;
        const float* beta = sim->friction_d;   /* may be NULL */

        /* Friction on the collective's freshly-computed `next`, before rotate. */
        if (beta) {
            for (int f = 0; f < GR_FIELD_COUNT; f++) {
                float* prev = sim->fields[f].prev;
                float* next = sim->fields[f].next;
                for (size_t k = 0; k < ncells_fr; k++)
                    next[k] = (next[k] + beta[k] * prev[k]) / (1.0f + beta[k]);
            }
        }

        /* Three-pointer rotation per field (collective). */
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            float* tmp           = sim->fields[f].prev;
            sim->fields[f].prev  = sim->fields[f].curr;
            sim->fields[f].curr  = sim->fields[f].next;
            sim->fields[f].next  = tmp;
        }

        /* v39: step each per-particle self-field set, applying the SAME
         * friction to its freshly-computed `next` before rotating.
         *
         * FIX: previously the self-set friction ran in the block above --
         * before gr_field_leapfrog_step_self had computed `next` -- so it
         * damped a stale buffer that step_self then overwrote.  The self
         * sets therefore got NO friction while the collective did, so the
         * two diverged and a lone particle kept a residual self-force.
         * Applying friction here (post-step_self) keeps collective == self. */
        if (sim->self_field_sets) {
            for (int pi = 0; pi < sim->n_particles; pi++) {
                gr_self_field_set_t* sf = sim->self_field_sets[pi];
                if (!sf) continue;
                gr_field_leapfrog_step_self(sim, sf);
                if (beta) {
                    for (int f = 0; f < GR_FIELD_COUNT; f++) {
                        float* prev = sf->fields[f].prev;
                        float* next = sf->fields[f].next;
                        for (size_t k = 0; k < ncells_fr; k++)
                            next[k] = (next[k] + beta[k] * prev[k]) / (1.0f + beta[k]);
                    }
                }
                for (int f = 0; f < GR_FIELD_COUNT; f++) {
                    float* tmp          = sf->fields[f].prev;
                    sf->fields[f].prev  = sf->fields[f].curr;
                    sf->fields[f].curr  = sf->fields[f].next;
                    sf->fields[f].next  = tmp;
                }
            }
        }

        /* v42 zero-mean gauge fix.  Under Neumann outer BC the DC mode
         * of the scalar potentials is unconstrained, and the wave equation
         * accelerates it whenever the spatial source mean is nonzero
         * (Phi_DC_dotdot = sc * rho_avg).  Subtract the mean from prev
         * and curr (same value, so the encoded time derivative
         * curr - prev is preserved) to keep the DC mode at zero. */
        if (sim->zero_mean_scalar_potentials) {
            const size_t ncells = (size_t) sim->width * (size_t) sim->height;
            const int scalar_ids[2] = { GR_FIELD_PHI_GRAV, GR_FIELD_PHI_EM };
            for (int idx = 0; idx < 2; idx++) {
                const int f = scalar_ids[idx];
                float* curr = sim->fields[f].curr;
                float* prev = sim->fields[f].prev;
                double sum = 0.0;
                for (size_t k = 0; k < ncells; k++) sum += (double) curr[k];
                const float mean = (float)(sum / (double) ncells);
                if (mean != 0.0f) {
                    for (size_t k = 0; k < ncells; k++) {
                        curr[k] -= mean;
                        prev[k] -= mean;
                    }
                }
            }
        }

        /* v40 parabolic Lorenz gauge cleaning: applied to phi and Phi_g
         * after each leapfrog + buffer rotation.  Bounds the gauge drift
         * from inconsistent source deposition (eg the spin dipole's
         * approximately-but-not-exactly div-free J).  Default ON. */
        gr_field_parabolic_gauge_clean(sim);
    }
    /* Stage 7+: push particles each step (Boris-leapfrog kick-drift).
     * Skipped while particles_frozen is set so the v38 §15.9 convergence
     * iteration can drive the fields to the static Poisson solution at
     * the initial particle configuration. */
    if (sim->n_particles > 0 && !sim->particles_frozen) gr_particle_push_all(sim);
    sim->step_count++;
}

/* v38 §15.9: freeze particle motion so the field leapfrog can iterate to
 * the static (Poisson) solution at the initial configuration before
 * dynamics begin.  Sources still deposit (positions and velocities are
 * read from the frozen particle state), the field still evolves, only the
 * Boris kick-drift on each particle is suppressed. */
void gr_sim_set_particles_frozen(gr_sim_t* sim, int frozen) {
    if (!sim) return;
    sim->particles_frozen = frozen ? 1 : 0;
}
int gr_sim_get_particles_frozen(const gr_sim_t* sim) {
    return sim ? sim->particles_frozen : 0;
}

void gr_sim_set_outer_bc_neumann(gr_sim_t* sim, int neumann) {
    if (!sim) return;
    sim->outer_bc_neumann = neumann ? 1 : 0;
}
int gr_sim_get_outer_bc_neumann(const gr_sim_t* sim) {
    return sim ? sim->outer_bc_neumann : 0;
}

void gr_sim_set_zero_mean_scalar_potentials(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->zero_mean_scalar_potentials = enabled ? 1 : 0;
}
int gr_sim_get_zero_mean_scalar_potentials(const gr_sim_t* sim) {
    return sim ? sim->zero_mean_scalar_potentials : 0;
}

/* v42 derivative-friction taper builder.  Per-cell value stored in
 * friction_d[k] is beta = gamma * dt, used in the centered-implicit
 * post-leapfrog correction:
 *   next[k] := (next_std[k] + beta[k] * prev[k]) / (1 + beta[k])
 *
 *   uniform_beta:   friction floor everywhere (catches the lowest
 *                   standing mode, which has nonzero d/dt content
 *                   throughout the box).  Typical: 1e-3 .. 5e-3.
 *
 *   taper_max_beta: friction at the outer wall.  Typical: 1e-2 ..
 *                   5e-2 so radiative wave packets attenuate by
 *                   several e-folds across the taper width.
 *
 *   taper_depth:    width of the taper from the wall inward, in cells.
 *                   Beyond taper_depth into the interior, beta ==
 *                   uniform.
 *
 * Profile: polynomial m=2 in the taper region.  Pass (0, 0, 0) to
 * deallocate. */
void gr_sim_set_volume_friction_taper(gr_sim_t* sim,
                                       float uniform_beta,
                                       float taper_max_beta,
                                       int   taper_depth) {
    if (!sim) return;
    free(sim->friction_d);
    sim->friction_d = NULL;
    if (uniform_beta <= 0.0f && (taper_max_beta <= 0.0f || taper_depth <= 0)) {
        return;
    }
    const int W = sim->width;
    const int H = sim->height;
    const size_t ncells = (size_t) W * (size_t) H;
    sim->friction_d = (float*) calloc(ncells, sizeof(float));
    if (!sim->friction_d) return;

    const float floor = (uniform_beta > 0.0f) ? uniform_beta : 0.0f;
    const float wall  = (taper_max_beta > 0.0f) ? taper_max_beta : 0.0f;
    const int   Nd    = (taper_depth > 0) ? taper_depth : 0;

    for (int j = 0; j < H; j++) {
        int dy = j; if (H - 1 - j < dy) dy = H - 1 - j;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            int dxi = i; if (W - 1 - i < dxi) dxi = W - 1 - i;
            const int depth = (dxi < dy) ? dxi : dy;
            float g;
            if (Nd > 0 && depth < Nd && wall > floor) {
                const float u = (float)(Nd - depth) / (float) Nd;
                g = floor + (wall - floor) * u * u;
            } else {
                g = floor;
            }
            if (g < 0.0f)   g = 0.0f;
            if (g > 0.999f) g = 0.999f;
            sim->friction_d[row + i] = g;
        }
    }
}
int gr_sim_volume_friction_enabled(const gr_sim_t* sim) {
    return (sim && sim->friction_d) ? 1 : 0;
}

/* Direct Poisson solver for Phi_g via SOR (Successive Over-Relaxation).
 *
 * The wave-equation convergence iteration from v38 §15.9 does not
 * actually damp the lowest standing mode of the absorber-bounded box --
 * that mode has nodes at the wall and an antinode at the center, so the
 * absorbing ring sees essentially no energy from it.  Running for many
 * mode periods leaves residual oscillation of order |Phi_static| (Dan,
 * 2026-06-04: stability probe shows swing ~50% of Phi_static even after
 * 5400-step warmup).
 *
 * This routine bypasses the wave equation entirely and solves the
 * static Poisson equation directly:
 *
 *     Lap Phi_g = -source_coeff * rho_matter        (per leapfrog form)
 *
 * 5-point Laplacian, zero-Dirichlet on the outer boundary, SOR with
 * omega = 2 / (1 + sin(pi / max(W,H))) which converges in ~max(W,H)
 * iterations.  After convergence prev := curr so the wave-equation
 * leapfrog continues from a true static state (zero time-derivative,
 * no excitation of standing modes).
 *
 * rho_matter must already be deposited (call gr_sim_step with frozen
 * particles or gr_sim_deposit_point_mass first). */
void gr_sim_relax_phi_g_poisson(gr_sim_t* sim, int n_iters) {
    if (!sim || n_iters <= 0) return;
    gr_field_state_t* f = &sim->fields[GR_FIELD_PHI_GRAV];
    const float* src = f->source;
    const float  sc  = f->source_coeff;
    if (!src) return;
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx2 = sim->dx * sim->dx;
    float* phi = f->curr;
    const size_t ncells = (size_t) W * (size_t) H;

    /* Zero-Dirichlet on the outer edge (matches leapfrog_field_critical). */
    for (int i = 0; i < W; i++) { phi[i] = 0.0f; phi[(H - 1) * W + i] = 0.0f; }
    for (int j = 0; j < H; j++) { phi[j * W] = 0.0f; phi[j * W + (W - 1)] = 0.0f; }

    const int N_for_omega = (W > H) ? W : H;
    const float omega = 2.0f / (1.0f + sinf((float) M_PI / (float) N_for_omega));

    /* Red-black Gauss-Seidel ordering would be friendlier to vectorization;
     * for clarity use plain row-major Gauss-Seidel with in-place updates. */
    for (int it = 0; it < n_iters; it++) {
        for (int j = 1; j < H - 1; j++) {
            const int row = j * W;
            for (int i = 1; i < W - 1; i++) {
                const int k = row + i;
                /* 5-pt Poisson: Lap Phi = -sc*src
                 *   (phi[k-1]+phi[k+1]+phi[k-W]+phi[k+W] - 4 phi[k]) / dx^2 = -sc * src[k]
                 * => phi[k] = (sum_neighbors + sc * src[k] * dx^2) / 4
                 * SOR update applied. */
                const float gauss = 0.25f * (phi[k - 1] + phi[k + 1]
                                            + phi[k - W] + phi[k + W]
                                            + sc * src[k] * dx2);
                phi[k] += omega * (gauss - phi[k]);
            }
        }
    }

    /* Lock as static: zero the wave-equation time derivative. */
    memcpy(f->prev, f->curr, ncells * sizeof(float));
    memcpy(f->next, f->curr, ncells * sizeof(float));
}

/* 2D Liénard-Wiechert initialization of all six perturbation potentials.
 *
 * For each particle p, every cell gets a contribution
 *
 *   Phi_g(x)  += -(sc_grav / 2pi)         * m_p     * ln |x - x_p|
 *   A_g_x(x)  += -(sc_grav / 2pi c^2)     * m_p v_x * ln |x - x_p|
 *   A_g_y(x)  += -(sc_grav / 2pi c^2)     * m_p v_y * ln |x - x_p|
 *   phi(x)    += -(sc_em   / 2pi)         * q_p     * ln |x - x_p|
 *   A_x(x)    += -(sc_em   / 2pi c^2)     * q_p v_x * ln |x - x_p|
 *   A_y(x)    += -(sc_em   / 2pi c^2)     * q_p v_y * ln |x - x_p|
 *
 * which is the 2D fundamental-solution Poisson result for each field's
 * wave-equation static limit ( Lap Phi = -sc * src ).  In the static
 * gauge with v = const, the vector potential reduces to v/c^2 times the
 * scalar potential -- which is the 2D Liénard-Wiechert form.
 *
 * Each particle's contribution is shifted by the mean of ln r over the
 * boundary ring (a pure gauge constant -- no gradient, hence no force,
 * is changed) so the summed field is ~0 at the Dirichlet wall and is
 * consistent with the discrete Dirichlet fixed point, leaving the wave
 * equation almost no boundary shock to launch.  This replaces the old
 * multiplicative taper (which distorted gradients in the ring); the
 * taper_inner/taper_outer args are retained for diagnostics but the
 * default path passes 0 (taper off).
 *
 * .prev = .curr = .next is set for every field so the wave equation
 * starts with bit-exact zero time derivative.
 *
 * Velocity is taken as v = p / (gamma m), matching the deposit code's
 * convention.  Particles with m == 0 contribute nothing to Phi_g/A_g;
 * particles with q == 0 contribute nothing to phi/A. */
/* Apply the boundary taper (if any) and lock prev=curr=next on a set of
 * potential fields.  Used for BOTH the collective sim->fields and each
 * per-particle self-field set, so a single particle's self set carries the
 * SAME L-W as its collective contribution -> exact self-force cancellation
 * (the difference field collective-self starts at 0 and, being source-free,
 * stays 0). */
static void lw_taper_and_lock(gr_field_state_t* fields, int W, int H,
                              int taper_inner, int taper_outer) {
    const size_t ncells = (size_t) W * (size_t) H;
    if (taper_inner > 0 && taper_inner > taper_outer) {
        const float tw = (float)(taper_inner - taper_outer);
        for (int j = 0; j < H; j++) {
            int dy = j; if (H - 1 - j < dy) dy = H - 1 - j;
            const int row = j * W;
            for (int i = 0; i < W; i++) {
                int dxi = i; if (W - 1 - i < dxi) dxi = W - 1 - i;
                const int depth = (dxi < dy) ? dxi : dy;
                if (depth >= taper_inner) continue;
                float taper;
                if (depth <= taper_outer) taper = 0.0f;
                else { const float u = (float)(taper_inner - depth) / tw; taper = 1.0f - u * u; }
                const int k = row + i;
                for (int f = 0; f < GR_FIELD_COUNT; f++)
                    if (fields[f].curr) fields[f].curr[k] *= taper;
            }
        }
    }
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        if (fields[f].curr && fields[f].prev) memcpy(fields[f].prev, fields[f].curr, ncells * sizeof(float));
        if (fields[f].curr && fields[f].next) memcpy(fields[f].next, fields[f].curr, ncells * sizeof(float));
    }
}

void gr_sim_init_potentials_lienard_wiechert_with_taper(gr_sim_t* sim,
                                                         int taper_inner,
                                                         int taper_outer) {
    if (!sim) return;
    if (taper_inner < 0) taper_inner = 0;
    if (taper_outer < 0) taper_outer = 0;
    /* taper_inner == 0 (or taper_inner <= taper_outer) disables the
     * taper -- L-W extends all the way to the wall.  Required when
     * using Neumann outer BC + derivative friction (the taper-induced
     * Phi(wall)=0 is incompatible with Neumann's Phi(wall)=Phi(wall-1)
     * and the discontinuity launches a wave inward each step). */
    const int   W       = sim->width;
    const int   H       = sim->height;
    const float dx      = sim->dx;
    const float c_eff   = sim->c_eff;
    const float inv_c2  = 1.0f / (c_eff * c_eff);
    const float inv_2pi = 1.0f / (2.0f * (float) M_PI);
    const float r_min   = 0.5f * dx;
    const size_t ncells = (size_t) W * (size_t) H;

    const float sc_grav = sim->fields[GR_FIELD_PHI_GRAV].source_coeff;
    const float sc_em   = sim->fields[GR_FIELD_PHI_EM  ].source_coeff;

    for (int f = 0; f < 6; f++) {
        if (sim->fields[f].curr) memset(sim->fields[f].curr, 0, ncells * sizeof(float));
    }

    float* phi_g = sim->fields[GR_FIELD_PHI_GRAV].curr;
    float* a_gx  = sim->fields[GR_FIELD_A_GX].curr;
    float* a_gy  = sim->fields[GR_FIELD_A_GY].curr;
    float* phi_e = sim->fields[GR_FIELD_PHI_EM].curr;
    float* a_x   = sim->fields[GR_FIELD_A_X].curr;
    float* a_y   = sim->fields[GR_FIELD_A_Y].curr;

    /* Direct sum over particles. */
    for (int pi = 0; pi < sim->n_particles; pi++) {
        const gr_particle_t* p = &sim->particles[pi];
        const float xp = p->x, yp = p->y;
        const float pmag2 = p->px * p->px + p->py * p->py;
        const float gamma = (p->mass > 0.0f)
            ? sqrtf(1.0f + pmag2 / (p->mass * p->mass * c_eff * c_eff))
            : 1.0f;
        const float vx = (p->mass > 0.0f) ? (p->px / (gamma * p->mass)) : 0.0f;
        const float vy = (p->mass > 0.0f) ? (p->py / (gamma * p->mass)) : 0.0f;

        const float kg   = -sc_grav * inv_2pi          * p->mass;
        const float kgA  = -sc_grav * inv_2pi * inv_c2 * p->mass;
        const float kqe  = -sc_em   * inv_2pi          * p->charge;
        const float kqA  = -sc_em   * inv_2pi * inv_c2 * p->charge;
        const int do_m = (p->mass   != 0.0f);
        const int do_q = (p->charge != 0.0f);

        /* If this particle has a self-field set, seed it with this particle's
         * OWN L-W contribution (zero it first, then accumulate just this
         * particle below).  Makes collective == self for the lone particle. */
        gr_self_field_set_t* self = (sim->self_field_sets) ? sim->self_field_sets[pi] : NULL;
        float *sg = NULL, *sgx = NULL, *sgy = NULL, *se = NULL, *sx = NULL, *sy = NULL;
        if (self) {
            sg  = self->fields[GR_FIELD_PHI_GRAV].curr;
            sgx = self->fields[GR_FIELD_A_GX].curr;
            sgy = self->fields[GR_FIELD_A_GY].curr;
            se  = self->fields[GR_FIELD_PHI_EM].curr;
            sx  = self->fields[GR_FIELD_A_X].curr;
            sy  = self->fields[GR_FIELD_A_Y].curr;
            memset(sg,  0, ncells * sizeof(float));
            memset(sgx, 0, ncells * sizeof(float));
            memset(sgy, 0, ncells * sizeof(float));
            memset(se,  0, ncells * sizeof(float));
            memset(sx,  0, ncells * sizeof(float));
            memset(sy,  0, ncells * sizeof(float));
        }

        /* Per-particle boundary-mean subtraction.  The 2D log GROWS toward
         * the wall: ln|x-x_p| ~ ln(R_wall) >> the near-field structure, so
         * the bare L-W is large and positive at the Dirichlet wall (forced
         * to 0), which launches the dominant startup transient.  Subtract
         * the mean of ln r over the boundary ring from this particle's
         * contribution, so the summed field is ~0 (mean) at the wall and
         * consistent with the Dirichlet fixed point.  This is a pure
         * constant (gauge) shift -- it changes no gradient, hence no force,
         * anywhere -- applied uniformly to all six components (the scalar
         * shift maps to A via the same v/c^2 factor).  Replaces the old
         * multiplicative taper, which distorted gradients in the ring. */
        double sum_lr = 0.0;
        long   n_b    = 0;
        for (int i = 0; i < W; i++) {
            const float dxx = (float) i * dx - xp;
            /* top row j=0 */
            {
                const float dyt = -yp;
                float r = sqrtf(dxx * dxx + dyt * dyt);
                if (r < r_min) r = r_min;
                sum_lr += (double) logf(r); n_b++;
            }
            /* bottom row j=H-1 */
            {
                const float dyb = (float)(H - 1) * dx - yp;
                float r = sqrtf(dxx * dxx + dyb * dyb);
                if (r < r_min) r = r_min;
                sum_lr += (double) logf(r); n_b++;
            }
        }
        for (int j = 1; j < H - 1; j++) {
            const float dyj = (float) j * dx - yp;
            /* left col i=0 */
            {
                const float dxl = -xp;
                float r = sqrtf(dxl * dxl + dyj * dyj);
                if (r < r_min) r = r_min;
                sum_lr += (double) logf(r); n_b++;
            }
            /* right col i=W-1 */
            {
                const float dxr = (float)(W - 1) * dx - xp;
                float r = sqrtf(dxr * dxr + dyj * dyj);
                if (r < r_min) r = r_min;
                sum_lr += (double) logf(r); n_b++;
            }
        }
        const float bmean = (n_b > 0) ? (float)(sum_lr / (double) n_b) : 0.0f;

        for (int j = 0; j < H; j++) {
            const float yc = (float) j * dx;
            const float dy = yc - yp;
            const float dy2 = dy * dy;
            const int row = j * W;
            for (int i = 0; i < W; i++) {
                const float xc  = (float) i * dx;
                const float dxx = xc - xp;
                float r2 = dxx * dxx + dy2;
                float r  = sqrtf(r2);
                if (r < r_min) r = r_min;
                const float lr = logf(r) - bmean;   /* boundary-mean-subtracted */
                const int k = row + i;
                if (do_m) {
                    const float vg = kg * lr, vgx = kgA * vx * lr, vgy = kgA * vy * lr;
                    phi_g[k] += vg;  a_gx[k] += vgx;  a_gy[k] += vgy;
                    if (self) { sg[k] += vg;  sgx[k] += vgx;  sgy[k] += vgy; }
                }
                if (do_q) {
                    const float ve = kqe * lr, vax = kqA * vx * lr, vay = kqA * vy * lr;
                    phi_e[k] += ve;  a_x[k] += vax;  a_y[k] += vay;
                    if (self) { se[k] += ve;  sx[k] += vax;  sy[k] += vay; }
                }
            }
        }
    }

    /* Boundary taper: smoothly drop the field to zero so the IC matches
     * the wave equation's wall BC and the L-W log-tail doesn't have to
     * be supported by the leapfrog all the way to the absorber.
     *
     *   taper(d) = 1 - u^2 where d = depth-from-nearest-wall in cells,
     *   u = (taper_inner - d) / (taper_inner - taper_outer) clipped to
     *   [0, 1].
     *
     *   d >= taper_inner            -> full L-W
     *   taper_outer < d < taper_inner -> smooth ramp
     *   d <= taper_outer            -> zero field
     *
     * Default (call to gr_sim_init_potentials_lienard_wiechert with no
     * taper args): taper_inner = n_damping, taper_outer = 0.  That
     * tapers across the damping ring, matching the wall to zero.
     *
     * Pulling the taper inward (larger taper_outer) gives a buffer of
     * zero field between the taper and the absorber -- useful as a
     * diagnostic for whether the boundary handling is responsible for
     * residual oscillation. */
    /* Taper (if any) + lock prev=curr=next, for the collective field AND
     * every per-particle self-field set, so collective and self carry the
     * identical L-W and the self-force cancels exactly. */
    lw_taper_and_lock(sim->fields, W, H, taper_inner, taper_outer);
    if (sim->self_field_sets) {
        for (int pi = 0; pi < sim->n_particles; pi++) {
            gr_self_field_set_t* self = sim->self_field_sets[pi];
            if (self) lw_taper_and_lock(self->fields, W, H, taper_inner, taper_outer);
        }
    }
}

/* Default-args entry: no taper.  Per-particle boundary-mean subtraction
 * (always applied inside _with_taper) makes the L-W field ~0 at the wall,
 * which is what consistency with the Dirichlet fixed point needs; the old
 * multiplicative taper is no longer used on the default path. */
void gr_sim_init_potentials_lienard_wiechert(gr_sim_t* sim) {
    if (!sim) return;
    gr_sim_init_potentials_lienard_wiechert_with_taper(sim, 0, 0);
}

/* v42: boundary-mean L-W init + a frozen-particle friction settle to the
 * exact discrete fixed point.
 *
 * The boundary-mean L-W field is a close initial guess but not the discrete
 * fixed point -- it differs from it by the near-source discrete-vs-continuum
 * mismatch (5-point Laplacian of the analytic log vs the deposited TSC
 * blob) plus the residual boundary spatial variation, ~4% peak.  Running
 * the field with the particles FROZEN and the derivative-friction absorber
 * on relaxes that residual to round-off in a few hundred steps (the friction
 * now damps every mode, including the lowest standing mode that the old
 * multiplicative ring missed -- this is what makes the v38 sec:15.9
 * convergence iteration finally work).  Releasing the particles then
 * launches no startup transient, so a symmetric configuration receives no
 * spurious center-of-mass impulse.
 *
 * Measured (stage62, 256x256 pic_binary): the quarter-orbit total-momentum
 * impulse drops from 2.3e-3 of m*v_orb with no settle to ~3e-5 (round-off)
 * with as few as 100 settle steps.  Requires friction configured
 * (gr_sim_set_volume_friction_taper) and particle source deposition on.
 *
 * step_count is reset to 0 after the settle so the released simulation
 * begins at t=0.  n_settle <= 0 reduces to a plain L-W init. */
void gr_sim_init_potentials_settled(gr_sim_t* sim, int n_settle) {
    if (!sim) return;
    gr_sim_init_potentials_lienard_wiechert(sim);
    if (n_settle > 0) {
        const int was_frozen = sim->particles_frozen;
        sim->particles_frozen = 1;
        for (int s = 0; s < n_settle; s++) gr_sim_step(sim);
        sim->particles_frozen = was_frozen;
    }
    sim->step_count = 0;
}

/* Stage 8 — skip the per-step leapfrog when no perturbation dynamics are
 * active. The perturbation fields stay at zero (as initialized) and the
 * particle pusher reads only the background array. */
void gr_sim_set_field_evolution(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->field_evolution_enabled = enabled ? 1 : 0;
}
int gr_sim_get_field_evolution(const gr_sim_t* sim) {
    return sim ? sim->field_evolution_enabled : 1;
}

void gr_sim_set_particle_source_deposition(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->particle_source_deposition = enabled ? 1 : 0;
}
int gr_sim_get_particle_source_deposition(const gr_sim_t* sim) {
    return sim ? sim->particle_source_deposition : 0;
}

void gr_sim_set_esirkepov_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->esirkepov_enabled = enabled ? 1 : 0;
}
int gr_sim_get_esirkepov_enabled(const gr_sim_t* sim) {
    return sim ? sim->esirkepov_enabled : 0;
}
int gr_sim_esirkepov_violations(const gr_sim_t* sim) {
    return sim ? sim->esirkepov_violations : 0;
}

void gr_sim_set_rho_smooth_passes(gr_sim_t* sim, int passes) {
    if (!sim) return;
    sim->rho_smooth_passes = (passes < 0) ? 0 : passes;
}
int gr_sim_get_rho_smooth_passes(const gr_sim_t* sim) {
    return sim ? sim->rho_smooth_passes : 0;
}

void gr_sim_set_j_smooth_passes(gr_sim_t* sim, int passes) {
    if (!sim) return;
    sim->j_smooth_passes = (passes < 0) ? 0 : passes;
}
int gr_sim_get_j_smooth_passes(const gr_sim_t* sim) {
    return sim ? sim->j_smooth_passes : 0;
}

void gr_sim_set_shape_function(gr_sim_t* sim, gr_shape_function_t s) {
    if (!sim) return;
    sim->shape_function = s;
}
gr_shape_function_t gr_sim_get_shape_function(const gr_sim_t* sim) {
    return sim ? sim->shape_function : GR_SHAPE_CIC;
}

void gr_sim_set_force_interp(gr_sim_t* sim, gr_force_interp_t scheme) {
    if (!sim) return;
    sim->force_interp = scheme;
}
gr_force_interp_t gr_sim_get_force_interp(const gr_sim_t* sim) {
    return sim ? sim->force_interp : GR_FORCE_INTERP_LEGACY;
}

void gr_sim_set_periodic_bc(gr_sim_t* sim, int periodic) {
    if (!sim) return;
    sim->periodic_bc = periodic ? 1 : 0;
}
int gr_sim_get_periodic_bc(const gr_sim_t* sim) {
    return sim ? sim->periodic_bc : 0;
}

/* Background evaluation mode — runtime switch between sampled-grid and
 * closed-form analytic paths for the installed background generator. */
void gr_sim_set_bg_mode(gr_sim_t* sim, gr_bg_mode_t mode) {
    if (!sim) return;
    sim->bg_mode = mode;
}
gr_bg_mode_t gr_sim_get_bg_mode(const gr_sim_t* sim) {
    return sim ? sim->bg_mode : GR_BG_MODE_SAMPLED;
}
gr_bg_kind_t gr_sim_get_bg_kind(const gr_sim_t* sim) {
    return sim ? sim->bg_kind : GR_BG_KIND_NONE;
}

void gr_sim_step_n(gr_sim_t* sim, int n) {
    if (!sim || n <= 0) return;
    for (int i = 0; i < n; i++) gr_sim_step(sim);
}

int   gr_sim_step_count(const gr_sim_t* sim) { return sim ? sim->step_count : 0; }
float gr_sim_time(const gr_sim_t* sim)       { return sim ? sim->dt * (float) sim->step_count : 0.0f; }
void  gr_sim_reset_time(gr_sim_t* sim)       { if (sim) sim->step_count = 0; }
float gr_sim_dt(const gr_sim_t* sim)         { return sim ? sim->dt : 0.0f; }
float gr_sim_dx(const gr_sim_t* sim)         { return sim ? sim->dx : 0.0f; }
int   gr_sim_width(const gr_sim_t* sim)      { return sim ? sim->width : 0; }
int   gr_sim_height(const gr_sim_t* sim)     { return sim ? sim->height : 0; }

float* gr_sim_field_ptr(gr_sim_t* sim, gr_field_id_t which) {
    if (!sim || which < 0 || which >= GR_FIELD_COUNT) return NULL;
    return sim->fields[which].curr;
}

/* v35 sublattice classification.  In v35 the implementation has not yet
 * been migrated; storage is still cell-centered.  But the classification
 * table is the eventual one (per §9), and callers can begin using it now
 * without behavior change. */
gr_lattice_t gr_array_lattice(gr_array_id_t which) {
    switch (which) {
    case GR_ARR_PHI_GRAV:
    case GR_ARR_PHI_EM:
    case GR_ARR_RHO_MATTER:
    case GR_ARR_RHO_Q:        return GR_LATTICE_CORNER;
    case GR_ARR_A_GX:
    case GR_ARR_A_X:
    case GR_ARR_J_MX:
    case GR_ARR_J_QX:         return GR_LATTICE_X_EDGE;
    case GR_ARR_A_GY:
    case GR_ARR_A_Y:
    case GR_ARR_J_MY:
    case GR_ARR_J_QY:         return GR_LATTICE_Y_EDGE;
    default:                  return GR_LATTICE_CORNER;
    }
}

void gr_lattice_offset(gr_lattice_t lat, float* dx_out, float* dy_out) {
    if (!dx_out || !dy_out) return;
    switch (lat) {
    case GR_LATTICE_CORNER: *dx_out = 0.0f; *dy_out = 0.0f; break;
    case GR_LATTICE_X_EDGE: *dx_out = 0.5f; *dy_out = 0.0f; break;
    case GR_LATTICE_Y_EDGE: *dx_out = 0.0f; *dy_out = 0.5f; break;
    default:                *dx_out = 0.0f; *dy_out = 0.0f; break;
    }
}

float* gr_sim_array_ptr(gr_sim_t* sim, gr_array_id_t which) {
    if (!sim) return NULL;
    switch (which) {
    /* Potentials: route to gr_field_id_t's current-time-slice buffer.
     * The cast is safe because GR_ARR_PHI_GRAV..GR_ARR_A_Y == 0..5. */
    case GR_ARR_PHI_GRAV:
    case GR_ARR_A_GX:
    case GR_ARR_A_GY:
    case GR_ARR_PHI_EM:
    case GR_ARR_A_X:
    case GR_ARR_A_Y:          return sim->fields[(int) which].curr;
    /* Sources: the existing struct fields. */
    case GR_ARR_RHO_MATTER:   return sim->rho_matter;
    case GR_ARR_J_MX:         return sim->J_mx;
    case GR_ARR_J_MY:         return sim->J_my;
    case GR_ARR_RHO_Q:        return sim->rho_q;
    case GR_ARR_J_QX:         return sim->J_qx;
    case GR_ARR_J_QY:         return sim->J_qy;
    default:                  return NULL;
    }
}

void  gr_sim_set_G_eff(gr_sim_t* sim, float G_eff) {
    if (!sim) return;
    sim->G_eff = G_eff;
    gr_sim_recompute_source_coeffs(sim);
}
float gr_sim_get_G_eff(const gr_sim_t* sim) { return sim ? sim->G_eff : 0.0f; }

void  gr_sim_set_k_e(gr_sim_t* sim, float k_e) {
    if (!sim) return;
    sim->k_e = k_e;
    gr_sim_recompute_source_coeffs(sim);
}
float gr_sim_get_k_e(const gr_sim_t* sim) { return sim ? sim->k_e : 0.0f; }

/* Stage 8 — force-tier selector. calloc-zero gives GR_FORCE_NEWTONIAN by
 * default, matching Stage 7 behavior. */
void gr_sim_set_force_tier(gr_sim_t* sim, gr_force_tier_t tier) {
    if (!sim) return;
    sim->force_tier = tier;
}
gr_force_tier_t gr_sim_get_force_tier(const gr_sim_t* sim) {
    return sim ? sim->force_tier : GR_FORCE_NEWTONIAN;
}

void gr_sim_set_gravitomagnetic_force_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->gravitomagnetic_force_enabled = enabled ? 1 : 0;
}
int gr_sim_get_gravitomagnetic_force_enabled(const gr_sim_t* sim) {
    return sim ? sim->gravitomagnetic_force_enabled : 0;
}

void gr_sim_set_gravitomagnetic_inductive_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->gravitomagnetic_inductive_enabled = enabled ? 1 : 0;
}
int gr_sim_get_gravitomagnetic_inductive_enabled(const gr_sim_t* sim) {
    return sim ? sim->gravitomagnetic_inductive_enabled : 0;
}

void gr_sim_set_em_lorentz_force_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->em_lorentz_force_enabled = enabled ? 1 : 0;
}
int gr_sim_get_em_lorentz_force_enabled(const gr_sim_t* sim) {
    return sim ? sim->em_lorentz_force_enabled : 0;
}

void gr_sim_set_em_inductive_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->em_inductive_enabled = enabled ? 1 : 0;
}
int gr_sim_get_em_inductive_enabled(const gr_sim_t* sim) {
    return sim ? sim->em_inductive_enabled : 0;
}

void gr_sim_set_em_electrostatic_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->em_electrostatic_enabled = enabled ? 1 : 0;
}
int gr_sim_get_em_electrostatic_enabled(const gr_sim_t* sim) {
    return sim ? sim->em_electrostatic_enabled : 0;
}

void gr_sim_set_em_magnetic_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->em_magnetic_enabled = enabled ? 1 : 0;
}
int gr_sim_get_em_magnetic_enabled(const gr_sim_t* sim) {
    return sim ? sim->em_magnetic_enabled : 0;
}

void gr_sim_set_parabolic_gauge_cleaning_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->parabolic_gauge_cleaning_enabled = enabled ? 1 : 0;
}
int gr_sim_get_parabolic_gauge_cleaning_enabled(const gr_sim_t* sim) {
    return sim ? sim->parabolic_gauge_cleaning_enabled : 0;
}
void gr_sim_set_parabolic_gauge_R(gr_sim_t* sim, float R) {
    if (!sim) return;
    sim->parabolic_gauge_R = (R > 0.0f) ? R : 0.0f;
}
float gr_sim_get_parabolic_gauge_R(const gr_sim_t* sim) {
    return sim ? sim->parabolic_gauge_R : 0.0f;
}

void gr_sim_set_em_stress_energy_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->em_stress_energy_enabled = enabled ? 1 : 0;
}
int gr_sim_get_em_stress_energy_enabled(const gr_sim_t* sim) {
    return sim ? sim->em_stress_energy_enabled : 0;
}

/* v41 -- EM stress-energy contribution to gravity sources
 * (gr_sandbox_v38.tex sec:alg_2d step 3).
 *
 * For each interior corner cell, compute the local EM perturbation
 * energy density and Poynting flux from the perturbation potentials
 * and add to rho_matter / J_mx / J_my.  We work in our unit convention
 * eps_0 = 1/(4 pi k_e), mu_0 = 4 pi k_e / c^2:
 *
 *   eps_0 / (2 c^2)  = 1 / (8 pi k_e c^2)
 *   1 / (mu_0 c^2)   = 1 / (4 pi k_e)
 *
 * Approximation: all quantities are evaluated at the corner sublattice
 * via the natural finite-difference / averaging stencils.  J_EM_x is
 * deposited into J_mx at the same storage index (effective x-edge
 * position is dx/2 off the corner; tolerable for smooth field
 * configurations).  Same for J_EM_y. */
void gr_sim_apply_em_stress_energy_sources(gr_sim_t* sim) {
    if (!sim || !sim->em_stress_energy_enabled) return;
    if (sim->k_e <= 0.0f) return;          /* coupling undefined */
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx  = sim->dx;
    const float dt  = sim->dt;
    const float c   = sim->c_eff;
    const float c2  = c * c;
    const float pi  = 3.14159265358979323846f;
    const float eps0_over_2c2 = 1.0f / (8.0f * pi * sim->k_e * c2);
    const float inv_mu0c2     = 1.0f / (4.0f * pi * sim->k_e);
    const float inv_dx        = 1.0f / dx;
    const float inv_2dx       = 1.0f / (2.0f * dx);
    const float inv_dt        = 1.0f / dt;

    const float* phi  = sim->fields[GR_FIELD_PHI_EM].curr;
    const float* Ax_c = sim->fields[GR_FIELD_A_X].curr;
    const float* Ax_p = sim->fields[GR_FIELD_A_X].prev;
    const float* Ay_c = sim->fields[GR_FIELD_A_Y].curr;
    const float* Ay_p = sim->fields[GR_FIELD_A_Y].prev;

    for (int j = 1; j < H - 1; j++) {
        for (int i = 1; i < W - 1; i++) {
            const int k = j * W + i;
            /* grad phi at corner (i,j): centered differences. */
            const float dphi_dx = (phi[k + 1] - phi[k - 1]) * inv_2dx;
            const float dphi_dy = (phi[k + W] - phi[k - W]) * inv_2dx;
            /* d_t A averaged to corner from edge values.
             * Ax storage at [j*W+i] is x-edge at (i+0.5, j); the two
             * x-edges flanking corner (i, j) are at storage [k-1] and [k]. */
            const float Ax_corner_c = 0.5f * (Ax_c[k - 1] + Ax_c[k]);
            const float Ax_corner_p = 0.5f * (Ax_p[k - 1] + Ax_p[k]);
            const float Ay_corner_c = 0.5f * (Ay_c[k - W] + Ay_c[k]);
            const float Ay_corner_p = 0.5f * (Ay_p[k - W] + Ay_p[k]);
            const float dtAx = (Ax_corner_c - Ax_corner_p) * inv_dt;
            const float dtAy = (Ay_corner_c - Ay_corner_p) * inv_dt;
            /* (grad phi + d_t A) -- E = -(this) but we square or take J */
            const float ex = dphi_dx + dtAx;
            const float ey = dphi_dy + dtAy;
            const float e2 = ex * ex + ey * ey;
            /* B_z = curl A at the cell center (i+0.5, j+0.5).  Standard
             * Yee curl: A_y[j, i+1] - A_y[j, i] minus A_x[j+1, i] - A_x[j, i],
             * divided by dx.  Used as the B_z value at corner (i, j) --
             * half-cell offset is a small approximation for smooth fields. */
            const float Bz = (Ay_c[k + 1] - Ay_c[k]) * inv_dx
                           - (Ax_c[k + W] - Ax_c[k]) * inv_dx;
            /* rho_EM = eps_0 / (2 c^2) * (E^2 + c^2 B_z^2) */
            const float rho_em = eps0_over_2c2 * (e2 + c2 * Bz * Bz);
            sim->rho_matter[k] += rho_em;
            /* J_EM = -F_12/(mu_0 c^2) (d_y phi + d_t A_y, -(d_x phi + d_t A_x))
             * which in our (ex, ey) variables is -B_z * (ey, -ex) / (mu_0 c^2). */
            sim->J_mx[k] += -Bz * ey * inv_mu0c2;
            sim->J_my[k] += +Bz * ex * inv_mu0c2;
        }
    }
}

/* Apply the recommended pedagogical defaults bundle.
 * See gr_sandbox_v38 sec kernel_design for empirical justification:
 *   - GR_SHAPE_BUMP at R=6 gives ~5-15x better self-force suppression
 *     than the bare CIC/TSC defaults across orbital v.
 *   - LEWIS_BIRDSALL force interp is uniformly better than LEGACY
 *     (Stage 18 result).
 *   - PML width 16 cells absorbs outgoing waves with O(0.1%) reflection.
 *   - particle_source_deposition = 1 is what closed-loop PIC needs;
 *     leaving it 0 is a footgun (particles wouldn't deposit at all). */
void gr_sim_use_pedagogical_defaults(gr_sim_t* sim) {
    if (!sim) return;
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 6.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_damping(sim, 16);
}

void gr_sim_set_kernel_radius(gr_sim_t* sim, float r) {
    if (!sim) return;
    if (r < 0.5f) r = 0.5f;       /* below 0.5: degenerate / single-cell */
    sim->kernel_radius = r;
}
float gr_sim_get_kernel_radius(const gr_sim_t* sim) {
    return sim ? sim->kernel_radius : 1.5f;
}

void gr_sim_set_em_shapiro_enabled(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    const int want = enabled ? 1 : 0;
    sim->em_shapiro_enabled = want;
    if (want) {
        gr_em_shapiro_recompute_c_local2(sim);
    }
    /* Note: we leave the c_local2 arrays allocated on disable so that
     * re-enable is cheap.  Memory is freed in gr_sim_destroy. */
}
int gr_sim_get_em_shapiro_enabled(const gr_sim_t* sim) {
    return sim ? sim->em_shapiro_enabled : 0;
}

void gr_sim_set_em_shapiro_dynamic(gr_sim_t* sim, int enabled) {
    if (!sim) return;
    sim->em_shapiro_dynamic = enabled ? 1 : 0;
    /* Make sure the arrays exist + reflect the current total Phi_g now. */
    if (sim->em_shapiro_dynamic && sim->em_shapiro_enabled) {
        gr_em_shapiro_recompute_c_local2(sim);
    }
}
int gr_sim_get_em_shapiro_dynamic(const gr_sim_t* sim) {
    return sim ? sim->em_shapiro_dynamic : 0;
}

const float* gr_sim_rho_matter_ptr(const gr_sim_t* sim) {
    return sim ? sim->rho_matter : NULL;
}
const float* gr_sim_rho_q_ptr(const gr_sim_t* sim) {
    return sim ? sim->rho_q : NULL;
}

void gr_sim_clear_sources(gr_sim_t* sim) {
    if (!sim) return;
    const size_t n = (size_t) sim->width * (size_t) sim->height;
    if (sim->rho_matter) memset(sim->rho_matter, 0, n * sizeof(float));
    if (sim->J_mx)       memset(sim->J_mx,       0, n * sizeof(float));
    if (sim->J_my)       memset(sim->J_my,       0, n * sizeof(float));
    if (sim->rho_q)      memset(sim->rho_q,      0, n * sizeof(float));
    if (sim->J_qx)       memset(sim->J_qx,       0, n * sizeof(float));
    if (sim->J_qy)       memset(sim->J_qy,       0, n * sizeof(float));
}

/* Zero every perturbation field buffer (prev/curr/next for all six potentials)
 * plus the current-history buffers.  Needed so a rebuild/reset starts from a
 * quiet field regardless of the init method -- otherwise an 'init=zero' (or
 * 'none') scene keeps the previous run's wave because nothing overwrites the
 * leapfrog buffers.  The L-W / settled inits overwrite the field anyway, so
 * calling this first is always safe. */
void gr_sim_clear_fields(gr_sim_t* sim) {
    if (!sim) return;
    const size_t n = (size_t) sim->width * (size_t) sim->height;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        if (sim->fields[f].prev) memset(sim->fields[f].prev, 0, n * sizeof(float));
        if (sim->fields[f].curr) memset(sim->fields[f].curr, 0, n * sizeof(float));
        if (sim->fields[f].next) memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    if (sim->J_mx_prev) memset(sim->J_mx_prev, 0, n * sizeof(float));
    if (sim->J_my_prev) memset(sim->J_my_prev, 0, n * sizeof(float));
    if (sim->J_qx_prev) memset(sim->J_qx_prev, 0, n * sizeof(float));
    if (sim->J_qy_prev) memset(sim->J_qy_prev, 0, n * sizeof(float));
}

void gr_sim_deposit_point_mass(gr_sim_t* sim, float x, float y, float mass) {
    if (!sim || !sim->rho_matter) return;
    /* rho_matter lives on the CORNER sublattice (§9, v35). */
    gr_cic_deposit_corner(sim->rho_matter, sim->width, sim->height, sim->dx, x, y, mass);
}
void gr_sim_deposit_point_charge(gr_sim_t* sim, float x, float y, float charge) {
    if (!sim || !sim->rho_q) return;
    gr_cic_deposit_corner(sim->rho_q, sim->width, sim->height, sim->dx, x, y, charge);
}

/* Stage 5 — composite deposit for a moving particle. Uses a single CIC kernel
 * (W_2 / bilinear) for all six contributions, satisfying the §9.5 adjoint
 * condition automatically since deposit and interpolation share the same
 * weights. Note: holding the particle position fixed across timesteps while
 * J = rho*v is nonzero violates the discrete continuity equation by exactly
 * v . grad(rho); Stage 5's test measures this directly. */
void gr_sim_deposit_point_particle(gr_sim_t* sim, float x, float y,
                                   float mass, float charge,
                                   float vx, float vy) {
    if (!sim) return;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    /* Per §9, rho deposits to CORNER, J_x to X_EDGE, J_y to Y_EDGE. */
    if (mass   != 0.0f) gr_cic_deposit_corner(sim->rho_matter, W, H, dx, x, y, mass);
    if (charge != 0.0f) gr_cic_deposit_corner(sim->rho_q,      W, H, dx, x, y, charge);
    if (mass   != 0.0f && vx != 0.0f) gr_cic_deposit_xedge(sim->J_mx, W, H, dx, x, y, mass   * vx);
    if (mass   != 0.0f && vy != 0.0f) gr_cic_deposit_yedge(sim->J_my, W, H, dx, x, y, mass   * vy);
    if (charge != 0.0f && vx != 0.0f) gr_cic_deposit_xedge(sim->J_qx, W, H, dx, x, y, charge * vx);
    if (charge != 0.0f && vy != 0.0f) gr_cic_deposit_yedge(sim->J_qy, W, H, dx, x, y, charge * vy);
}

const float* gr_sim_J_mx_ptr(const gr_sim_t* sim) { return sim ? sim->J_mx : NULL; }
const float* gr_sim_J_my_ptr(const gr_sim_t* sim) { return sim ? sim->J_my : NULL; }
const float* gr_sim_J_qx_ptr(const gr_sim_t* sim) { return sim ? sim->J_qx : NULL; }
const float* gr_sim_J_qy_ptr(const gr_sim_t* sim) { return sim ? sim->J_qy : NULL; }

int gr_sim_load_scenario(gr_sim_t* sim, const char* name, const float* params, int n_params) {
    if (!sim || !name) return -1;
    gr_scenarios_init();
    const gr_scenario_t* s = gr_scenario_find(name);
    if (!s) return -2;
    return s->build(sim, params, n_params);
}

int gr_sim_damping_layers(const gr_sim_t* sim) { return sim ? sim->n_damping : 0; }

/* Profile envelope at normalized depth u = d / L (u in [0, 1]).
 * Both kinds return 0 at u=0 (interior side of the layer) and 1 at u=1
 * (outer wall). */
static float damp_profile_envelope(gr_damp_profile_kind_t kind,
                                   float poly_order, float exp_beta,
                                   float u) {
    if (u <= 0.0f) return 0.0f;
    if (u >= 1.0f) return 1.0f;
    switch (kind) {
    case GR_DAMP_POLYNOMIAL:
        return powf(u, poly_order);
    case GR_DAMP_EXPONENTIAL: {
        /* (e^(beta u) - 1) / (e^beta - 1).  Stable for moderate beta. */
        const float num = expf(exp_beta * u) - 1.0f;
        const float den = expf(exp_beta)     - 1.0f;
        return num / den;
    }
    default: return powf(u, poly_order);
    }
}

/* Derived sigma_max from a target round-trip reflection R via the
 * integral identity 2 (sigma_max / c) integral_0^L f(d/L) dd = -ln(R).
 *
 *   integral_0^1 (u^m)            du = 1 / (m+1)
 *   integral_0^1 (e^(beta u) - 1)/(e^beta - 1) du = 1/beta - 1/(e^beta - 1)
 *
 * Result: sigma_max = -c ln(R) / (2 L I), where I is the integral above. */
static float damp_sigma_max_from_R(gr_damp_profile_kind_t kind,
                                   float poly_order, float exp_beta,
                                   float c, float L, float R) {
    if (!(R > 0.0f) || !(R < 1.0f) || !(L > 0.0f)) return 0.0f;
    float I = 1.0f;
    switch (kind) {
    case GR_DAMP_POLYNOMIAL:
        I = 1.0f / (poly_order + 1.0f);
        break;
    case GR_DAMP_EXPONENTIAL: {
        const float eb = expf(exp_beta);
        I = (1.0f / exp_beta) - (1.0f / (eb - 1.0f));
        if (I <= 0.0f) I = 1.0f / (poly_order + 1.0f);  /* defensive */
        break;
    }
    }
    return -c * logf(R) / (2.0f * L * I);
}

void gr_sim_set_damping_config(gr_sim_t* sim, const gr_damp_config_t* cfg) {
    if (!sim || !cfg) return;
    free(sim->damping_d);
    sim->damping_d = NULL;
    sim->n_damping = 0;
    sim->damp_kind = GR_DAMP_POLYNOMIAL;
    sim->damp_poly_order = 0.0f;
    sim->damp_exp_beta = 0.0f;
    sim->damp_target_reflection = 0.0f;
    sim->damp_sigma_max_used = 0.0f;

    const int n_damping = cfg->n_damping;
    if (n_damping <= 0) return;
    if (2 * n_damping >= sim->width || 2 * n_damping >= sim->height) return;

    const int W = sim->width;
    const int H = sim->height;
    if (W > 16384 || H > 16384) return;
    sim->damping_d = (float*) calloc((size_t) W * (size_t) H, sizeof(float));
    if (!sim->damping_d) return;

    /* Resolve config defaults — zero fields mean "use canonical defaults". */
    gr_damp_profile_kind_t kind = cfg->kind;
    float poly_order = cfg->poly_order > 0.0f ? cfg->poly_order : 2.0f;
    float exp_beta   = cfg->exp_beta   > 0.0f ? cfg->exp_beta   : 4.0f;
    float target_R   = cfg->target_reflection > 0.0f ? cfg->target_reflection : 1.0e-3f;

    const float L = (float) n_damping * sim->dx;
    float sigma_max = (cfg->sigma_max_override > 0.0f)
                          ? cfg->sigma_max_override
                          : damp_sigma_max_from_R(kind, poly_order, exp_beta,
                                                  sim->c_eff, L, target_R);

    const float inv_Nd = 1.0f / (float) n_damping;
    const float dt     = sim->dt;

    /* CFL stability enforcement.  See gr_sim_damping_max_stable_sigma_dt
     * comment in grlite.h.  For CRITICAL, the Nyquist mode requires
     *   4 * gamma >= 8 * CFL^2  =>  sigma*dt <= 1 - 2*CFL^2 .
     * Clamp sigma_max if needed.  For MULTIPLICATIVE the formal bound is
     * sigma*dt < 1 (so the (1-sigma*dt) prefactor remains positive). */
    {
        const float cfl_now = sim->cfl;
        float       max_sigma_dt;
        if (cfg->time_form == GR_DAMP_TIME_CRITICAL) {
            max_sigma_dt = 1.0f - 2.0f * cfl_now * cfl_now;
            if (max_sigma_dt < 0.0f) max_sigma_dt = 0.0f;
        } else {
            max_sigma_dt = 0.999f;
        }
        if (sigma_max * dt > max_sigma_dt) {
            sigma_max = (max_sigma_dt > 0.0f) ? (max_sigma_dt / dt) : 0.0f;
        }
    }

    /* Heap-allocate fx/fy: at W=H=16384 these would be 128KB on the
     * stack, which overflows Emscripten's default 64KB stack. */
    float* fx = (float*) malloc((size_t) W * sizeof(float));
    float* fy = (float*) malloc((size_t) H * sizeof(float));
    if (!fx || !fy) { free(fx); free(fy); return; }
    for (int i = 0; i < W; i++) {
        int depth = 0;
        if (i < n_damping)            depth = n_damping - i;
        else if (i >= W - n_damping)  depth = i - (W - n_damping) + 1;
        const float u = (float) depth * inv_Nd;
        fx[i] = damp_profile_envelope(kind, poly_order, exp_beta, u);
    }
    for (int j = 0; j < H; j++) {
        int depth = 0;
        if (j < n_damping)            depth = n_damping - j;
        else if (j >= H - n_damping)  depth = j - (H - n_damping) + 1;
        const float u = (float) depth * inv_Nd;
        fy[j] = damp_profile_envelope(kind, poly_order, exp_beta, u);
    }
    const float sigma_max_dt = sigma_max * dt;
    for (int j = 0; j < H; j++) {
        const int   row = j * W;
        const float fy_j = fy[j];
        for (int i = 0; i < W; i++) {
            const float f_max = (fx[i] > fy_j) ? fx[i] : fy_j;
            sim->damping_d[row + i] = sigma_max_dt * f_max;
        }
    }
    free(fx);
    free(fy);

    sim->n_damping              = n_damping;
    sim->damp_kind              = kind;
    sim->damp_poly_order        = poly_order;
    sim->damp_exp_beta          = exp_beta;
    sim->damp_target_reflection = target_R;
    sim->damp_sigma_max_used    = sigma_max;
    sim->damp_time_form         = cfg->time_form;
}

float gr_sim_damping_max_stable_sigma_dt(const gr_sim_t* sim, gr_damp_time_form_t form) {
    if (!sim) return 0.0f;
    const float cfl = sim->cfl;
    if (form == GR_DAMP_TIME_CRITICAL) {
        const float v = 1.0f - 2.0f * cfl * cfl;
        return v > 0.0f ? v : 0.0f;
    }
    return 0.999f;
}

gr_damp_config_t gr_sim_get_damping_config(const gr_sim_t* sim) {
    gr_damp_config_t out = { 0, GR_DAMP_POLYNOMIAL, 2.0f, 4.0f, 1.0e-3f, 0.0f, GR_DAMP_TIME_MULTIPLICATIVE };
    if (!sim || sim->n_damping <= 0) return out;
    out.n_damping          = sim->n_damping;
    out.kind               = sim->damp_kind;
    out.poly_order         = sim->damp_poly_order;
    out.exp_beta           = sim->damp_exp_beta;
    out.target_reflection  = sim->damp_target_reflection;
    out.sigma_max_override = sim->damp_sigma_max_used;
    out.time_form          = sim->damp_time_form;
    return out;
}

/* Legacy wrapper — uses the §9.6 spec default (polynomial m=2) with the
 * historical sigma_max = 21 c / (2 L) override.  The textbook formula
 * with target_reflection=1e-3 gives sigma_max = 20.72 c / (2 L) (about
 * 1.3% lower), so we override to keep bit-exact backward compat with the
 * stage02_damping baseline. */
void gr_sim_set_damping(gr_sim_t* sim, int n_damping) {
    if (!sim || n_damping <= 0) {
        const gr_damp_config_t off = {0, GR_DAMP_POLYNOMIAL, 2.0f, 0.0f, 1.0e-3f, 0.0f,
                                      GR_DAMP_TIME_MULTIPLICATIVE};
        gr_sim_set_damping_config(sim, &off);
        return;
    }
    const float L = (float) n_damping * sim->dx;
    const gr_damp_config_t cfg = {
        .n_damping          = n_damping,
        .kind               = GR_DAMP_POLYNOMIAL,
        .poly_order         = 2.0f,
        .exp_beta           = 0.0f,
        .target_reflection  = 1.0e-3f,
        .sigma_max_override = 21.0f * sim->c_eff / (2.0f * L),
        .time_form          = GR_DAMP_TIME_MULTIPLICATIVE,
    };
    gr_sim_set_damping_config(sim, &cfg);
}
