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

    /* Stage 10: when nonzero, gr_sim_step auto-deposits every particle's
     * mass, charge, and currents onto the source arrays before each field
     * leapfrog runs.  Default 0 (no automatic deposition; scenarios deposit
     * manually if needed, e.g., Stage 5's static moving_source). */
    int particle_source_deposition;

    /* Stage 11+ (v35 §sec:yee_pivot, answer B1): when nonzero, current
     * deposition (J_mx, J_my, J_qx, J_qy) uses Esirkepov decomposition
     * (Esirkepov 2001) instead of direct CIC, guaranteeing exact discrete
     * continuity (rho^{n+1} - rho^n) / dt + div(J^{n-1/2}) = 0 cell-by-cell.
     * Default 1 (on).  Off-flag (set via gr_sim_set_esirkepov_enabled) is for
     * regression testing only — disabling it brings back PIC grid-heating
     * artifacts on moving sources. */
    int esirkepov_enabled;
    /* Count of timesteps in which a particle's motion exceeded 1 cell in x
     * or y (the 2-cell assumption is violated under CFL but a stiff force
     * impulse could still violate it).  On violation, Esirkepov falls back
     * to direct CIC for that particle that step, and continuity holds only
     * to truncation order rather than exactly.  Public via
     * gr_sim_esirkepov_violations(). */
    int esirkepov_violations;

    /* Background generator parameters — kept alongside the sampled phi_g_bg
     * array so the user can switch between SAMPLED and ANALYTIC at runtime.
     * Only the currently installed kind's fields are meaningful. */
    gr_bg_kind_t bg_kind;
    gr_bg_mode_t bg_mode;
    float        bg_x0, bg_y0;
    float        bg_GM;
    float        bg_eps;
    /* Spin angular momentum of the spinning point-mass generator (z-component
     * only in 2D — the spin axis is perpendicular to the simulation plane). */
    float        bg_Jz;
    /* Reserved slot for the charged variants (Stage 11+). */
    float        bg_charge;
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

/* Defined in deposit.c — per-sublattice CIC deposition (v35 Yee layout).
 * Each variant deposits `value` into the four cells surrounding the
 * sub-cell position (x_p, y_p), with bilinear weights computed against
 * the sublattice's own node offsets. */
void gr_cic_deposit_corner(float* arr, int W, int H, float dx,
                           float x_p, float y_p, float value);
void gr_cic_deposit_xedge (float* arr, int W, int H, float dx,
                           float x_p, float y_p, float value);
void gr_cic_deposit_yedge (float* arr, int W, int H, float dx,
                           float x_p, float y_p, float value);

/* Defined in deposit.c — Esirkepov 2D current deposition for a particle
 * moving from (x0, y0) to (x1, y1) over a timestep of length dt with
 * source coupling `source` (mass for J_m*, charge for J_q*).  Deposits J_x
 * onto the X_EDGE sublattice and J_y onto Y_EDGE such that, paired with a
 * corner-CIC rho deposit at the same endpoints, the discrete continuity
 * (rho^{n+1} - rho^n)/dt + div(J^{n-1/2}) = 0 holds exactly.
 *
 * Returns 1 on successful Esirkepov deposit (2-cell motion case);
 * 0 if either |x1 - x0| > dx or |y1 - y0| > dx (multi-cell crossing) — in
 * which case the caller falls back to gr_cic_deposit_x/yedge and bumps the
 * sim->esirkepov_violations counter. */
int gr_esirkepov_deposit_jxy(float* Jx, float* Jy,
                             int W, int H, float dx, float dt,
                             float x0, float y0, float x1, float y1,
                             float source);

/* Defined in scenarios/registry.c. */
const gr_scenario_t* gr_scenario_find(const char* name);

/* Defined in sim.c — refresh the per-field source coefficients after G_eff
 * or k_e changes. */
void gr_sim_recompute_source_coeffs(struct gr_sim* sim);

/* Defined in background.c — analytic-mode evaluation of the installed
 * background generator at an exact spatial position (x, y).  Returns 1 if a
 * background is installed and was evaluated, 0 otherwise.  Output:
 *   *phi_out  : value of Phi_g^{bg}(x, y)
 *   *gx_out   : d/dx Phi_g^{bg}(x, y)
 *   *gy_out   : d/dy Phi_g^{bg}(x, y) */
int gr_bg_eval_analytic(const struct gr_sim* sim, float x, float y,
                        float* phi_out, float* gx_out, float* gy_out);

/* Defined in background.c — analytic-mode evaluation of the gravitomagnetic
 * vector potential A_g(x, y) for the installed background generator.  Returns
 * 1 if the installed kind supplies a nonzero A_g (i.e., SPINNING_POINT_MASS),
 * 0 otherwise.  When 0 is returned the caller should treat A_g as zero. */
int gr_bg_eval_A_g(const struct gr_sim* sim, float x, float y,
                   float* Ax_out, float* Ay_out);

#endif /* GRLITE_SIM_INTERNAL_H */
