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
    const float c_grav  = -4.0f * (float) M_PI * sim->G_eff;
    const float c_em    = -4.0f * (float) M_PI * sim->k_e;
    sim->fields[GR_FIELD_PHI_GRAV].source_coeff = c_grav;
    sim->fields[GR_FIELD_A_GX    ].source_coeff = c_grav * inv_c2;
    sim->fields[GR_FIELD_A_GY    ].source_coeff = c_grav * inv_c2;
    sim->fields[GR_FIELD_PHI_EM  ].source_coeff = c_em;
    sim->fields[GR_FIELD_A_X     ].source_coeff = c_em   * inv_c2;
    sim->fields[GR_FIELD_A_Y     ].source_coeff = c_em   * inv_c2;
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
    free(sim->damping_d);
    free(sim->phi_g_bg);
    free(sim->Agx_bg);
    free(sim->Agy_bg);
    free(sim->phi_bg);
    free(sim->Ax_bg);
    free(sim->Ay_bg);
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
        for (int i = 0; i < sim->n_particles; i++) {
            const gr_particle_t* p = &sim->particles[i];
            const float pmag2 = p->px * p->px + p->py * p->py;
            const float gamma = sqrtf(1.0f + pmag2 / (p->mass * p->mass * c2));
            const float vx    = p->px / (gamma * p->mass);
            const float vy    = p->py / (gamma * p->mass);
            gr_sim_deposit_point_particle(sim, p->x, p->y,
                                          p->mass, p->charge, vx, vy);
        }
    }

    if (sim->field_evolution_enabled) {
        gr_field_leapfrog_step_all(sim);
        /* Three-pointer rotation per field. */
        for (int f = 0; f < GR_FIELD_COUNT; f++) {
            float* tmp           = sim->fields[f].prev;
            sim->fields[f].prev  = sim->fields[f].curr;
            sim->fields[f].curr  = sim->fields[f].next;
            sim->fields[f].next  = tmp;
        }
    }
    /* Stage 7+: push particles each step (Boris-leapfrog kick-drift). */
    if (sim->n_particles > 0) gr_particle_push_all(sim);
    sim->step_count++;
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

void gr_sim_deposit_point_mass(gr_sim_t* sim, float x, float y, float mass) {
    if (!sim || !sim->rho_matter) return;
    gr_cic_deposit_scalar(sim->rho_matter, sim->width, sim->height, sim->dx, x, y, mass);
}
void gr_sim_deposit_point_charge(gr_sim_t* sim, float x, float y, float charge) {
    if (!sim || !sim->rho_q) return;
    gr_cic_deposit_scalar(sim->rho_q, sim->width, sim->height, sim->dx, x, y, charge);
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
    if (mass   != 0.0f) gr_cic_deposit_scalar(sim->rho_matter, W, H, dx, x, y, mass);
    if (charge != 0.0f) gr_cic_deposit_scalar(sim->rho_q,      W, H, dx, x, y, charge);
    if (mass   != 0.0f && vx != 0.0f) gr_cic_deposit_scalar(sim->J_mx, W, H, dx, x, y, mass   * vx);
    if (mass   != 0.0f && vy != 0.0f) gr_cic_deposit_scalar(sim->J_my, W, H, dx, x, y, mass   * vy);
    if (charge != 0.0f && vx != 0.0f) gr_cic_deposit_scalar(sim->J_qx, W, H, dx, x, y, charge * vx);
    if (charge != 0.0f && vy != 0.0f) gr_cic_deposit_scalar(sim->J_qy, W, H, dx, x, y, charge * vy);
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

/* Eq. (eq:damp_profile) — §9.6 — same precomputation as before, unchanged by
 * the Stage 4 field-state refactor (damping is per-cell, not per-field). */
void gr_sim_set_damping(gr_sim_t* sim, int n_damping) {
    if (!sim) return;
    free(sim->damping_d);
    sim->damping_d = NULL;
    sim->n_damping = 0;
    if (n_damping <= 0) return;
    if (2 * n_damping >= sim->width || 2 * n_damping >= sim->height) return;

    const int W = sim->width;
    const int H = sim->height;
    sim->damping_d = (float*) calloc((size_t) W * (size_t) H, sizeof(float));
    if (!sim->damping_d) return;

    sim->n_damping = n_damping;
    const float L         = (float) n_damping * sim->dx;
    const float sigma_max = 21.0f * sim->c_eff / (2.0f * L);
    const float inv_Nd    = 1.0f / (float) n_damping;
    const float dt        = sim->dt;

    float fx2[16384], fy2[16384];
    if (W > 16384 || H > 16384) {
        free(sim->damping_d);
        sim->damping_d = NULL;
        sim->n_damping = 0;
        return;
    }
    for (int i = 0; i < W; i++) {
        int depth = 0;
        if (i < n_damping)            depth = n_damping - i;
        else if (i >= W - n_damping)  depth = i - (W - n_damping) + 1;
        const float f = (float) depth * inv_Nd;
        fx2[i] = f * f;
    }
    for (int j = 0; j < H; j++) {
        int depth = 0;
        if (j < n_damping)            depth = n_damping - j;
        else if (j >= H - n_damping)  depth = j - (H - n_damping) + 1;
        const float f = (float) depth * inv_Nd;
        fy2[j] = f * f;
    }
    const float sigma_max_dt = sigma_max * dt;
    for (int j = 0; j < H; j++) {
        const int   row = j * W;
        const float fy2_j = fy2[j];
        for (int i = 0; i < W; i++) {
            const float f2_max = (fx2[i] > fy2_j) ? fx2[i] : fy2_j;
            sim->damping_d[row + i] = sigma_max_dt * f2_max;
        }
    }
}
