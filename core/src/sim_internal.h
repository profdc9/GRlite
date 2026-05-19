/* Internal simulation struct — visible only to .c files inside core/. */
#ifndef GRLITE_SIM_INTERNAL_H
#define GRLITE_SIM_INTERNAL_H

#include "grlite.h"

/* Per-field state. Each of the six potentials owns three time levels
 * (prev = Phi^{n-1}, curr = Phi^n, next = Phi^{n+1}) plus a pointer to its
 * source array (rho_matter, J_mx, ..., owned by struct gr_sim) and a scalar
 * source coupling coefficient.
 *
 * The leapfrog update is (gr_sandbox_v32.tex §9.2 eq:leapfrog_field, with
 * the §9.7 source term):
 *    Phi^{n+1} = 2 Phi^n - Phi^{n-1} + (c dt)^2 (Lap Phi^n + source_coeff * source^n)
 * applied to all six fields in parallel each step; pointer rotation between
 * prev/curr/next is then done on each independently. */
typedef struct {
    float* prev;
    float* curr;
    float* next;
    const float* source;     /* aliased to a source array owned by gr_sim_t */
    float        source_coeff;
} gr_field_state_t;

struct gr_sim {
    int   width, height;
    float dx;
    float c_eff;
    float dt;
    float cfl;
    float G_eff;        /* gravitational coupling (Stage 3+); default 1.0 */
    float k_e;          /* electric coupling (Stage 4+); default 1.0 */
    int   step_count;

    /* The six potentials. Indexed by gr_field_id_t. */
    gr_field_state_t fields[GR_FIELD_COUNT];

    /* Source arrays — same staggering as fields (cell-centered).
     * Always allocated by gr_sim_create; zero-filled by default. */
    float* rho_matter;  /* sources Phi_g       */
    float* J_mx;        /* sources A_{g,x}     */
    float* J_my;        /* sources A_{g,y}     */
    float* rho_q;       /* sources phi (EM)    */
    float* J_qx;        /* sources A_x         */
    float* J_qy;        /* sources A_y         */

    /* Damping layer (Stage 2). Applied identically to all six fields, since
     * they all satisfy the same wave equation with the same c. */
    int    n_damping;
    float* damping_d;

    /* Background field arrays (Stage 6). Lazily allocated by
     * gr_sim_set_background_*. */
    float* phi_g_bg;
    float* Agx_bg;
    float* Agy_bg;
    float* phi_bg;
    float* Ax_bg;
    float* Ay_bg;

    /* Particles (Stage 7). Dynamic array of gr_particle_t. */
    gr_particle_t* particles;
    int            n_particles;
    int            particles_capacity;

    /* Stage 8 — selects Newtonian vs Tier-2 relativistic gravity force. */
    gr_force_tier_t force_tier;

    /* If 0, the per-step wave-equation leapfrog on the six perturbation
     * fields is skipped (along with the buffer rotation).  Used by Stage 7/8
     * static-background-only tests to avoid the cost of evolving zero fields
     * for tens of thousands of steps.  Default 1 (full step). */
    int field_evolution_enabled;

    /* Background generator parameters — kept alongside the sampled phi_g_bg
     * array so the user can switch between SAMPLED and ANALYTIC at runtime.
     * Only the currently installed kind's fields are meaningful. */
    gr_bg_kind_t bg_kind;
    gr_bg_mode_t bg_mode;
    float        bg_x0, bg_y0;
    float        bg_GM;
    float        bg_eps;
    /* Future-use slots for spinning / charged variants — kept as zeros until
     * the corresponding generator is implemented. */
    float        bg_charge;
    float        bg_spin;
};

/* Defined in field.c — steps all six fields in parallel. */
void gr_field_leapfrog_step_all(struct gr_sim* sim);

/* Defined in particle.c — pushes all particles one timestep (kick-drift).
 * Reads (background + perturbation) Phi_g at each particle position to
 * compute the gravitational force. Subsequent stages will add the EM/GEM
 * vector-potential contributions. */
void gr_particle_push_all(struct gr_sim* sim);

/* Interpolated total Phi_g at (x, y) — bg + curr perturbation, via the
 * cell-centered CIC kernel W_2. Used by both the force pusher and the
 * energy diagnostic. */
float gr_phi_g_total_at(const struct gr_sim* sim, float x, float y);

/* Defined in deposit.c — CIC deposition of a scalar value at sub-cell position.
 * Used by scenarios and by gr_sim_deposit_point_{mass,charge}. */
void gr_cic_deposit_scalar(float* rho, int W, int H, float dx,
                           float x_p, float y_p, float value);

/* Defined in scenarios/registry.c. */
const gr_scenario_t* gr_scenario_find(const char* name);

/* Defined in sim.c — refresh the per-field source coefficients after G_eff
 * or k_e changes. */
void gr_sim_recompute_source_coeffs(struct gr_sim* sim);

/* Defined in background.c — analytic-mode evaluation of the installed
 * background generator at an exact spatial position (x, y).  Returns 1 if a
 * background is installed and was evaluated, 0 otherwise (caller leaves
 * *phi_out / *gx_out / *gy_out as the caller initialized them).  Output:
 *   *phi_out  : value of Phi_g^{bg}(x, y)
 *   *gx_out   : d/dx Phi_g^{bg}(x, y)
 *   *gy_out   : d/dy Phi_g^{bg}(x, y) */
int gr_bg_eval_analytic(const struct gr_sim* sim, float x, float y,
                        float* phi_out, float* gx_out, float* gy_out);

#endif /* GRLITE_SIM_INTERNAL_H */
