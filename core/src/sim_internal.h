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
    /* Profile parameters used for the current damping_d (so callers can
     * inspect / round-trip via gr_sim_get_damping_config).  When n_damping
     * is 0, these are zero-initialized and ignored. */
    gr_damp_profile_kind_t damp_kind;
    float                  damp_poly_order;
    float                  damp_exp_beta;
    float                  damp_target_reflection;
    float                  damp_sigma_max_used;   /* the sigma_max actually
                                                   * baked into damping_d, in
                                                   * physical units (NOT *dt) */
    gr_damp_time_form_t    damp_time_form;        /* multiplicative vs critical */

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

    /* Tier-1 gravitomagnetic Lorentz force gate (Stage 20+).
     * Default 1 (enabled); 0 disables the +4 m v x B_g piece, useful for
     * isolating the proper-time clock effect from the orbit-shape Lense-
     * Thirring response. */
    int gravitomagnetic_force_enabled;

    /* EM Lorentz force gate (Stage 23+).  Default 1 (enabled); 0 disables
     * the q*(v x B_em) piece (and any future -q*grad phi_em / -q*d_t A
     * pieces).  Mirror of gravitomagnetic_force_enabled. */
    int em_lorentz_force_enabled;

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
    /* Number of binomial-smoothing passes to apply to rho_matter and rho_q
     * after deposit, before the field leapfrog reads them.  Each pass is
     * the [[1,2,1],[2,4,2],[1,2,1]]/16 3x3 stencil.  Default 0 (no smoothing). */
    int rho_smooth_passes;
    /* Shape function for rho deposit + force interp.  Default CIC. */
    gr_shape_function_t shape_function;
    /* Force-interpolation scheme.  LEGACY = FD-then-interp; LEWIS_BIRDSALL
     * = analytic-grad-of-shape interpolation (energy-conserving force
     * pairing per the discrete Lagrangian).  Default LEGACY. */
    gr_force_interp_t   force_interp;
    /* If nonzero, the field leapfrog uses periodic BC instead of zero-
     * Dirichlet at the box edges.  Periodic BC restores translation
     * invariance of the discrete Laplacian — required for HE self-force
     * cancellation to hold at any particle position, not just the box
     * symmetry center.  Outgoing waves wrap around (so damping is
     * essential here; with periodic+damping, wrap energy is absorbed by
     * the damping ring before it returns).  Default 0 (zero-Dirichlet). */
    int periodic_bc;
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
    /* Uniform gravitomagnetic field B_g_z, used by
     * GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC.  Symmetric-gauge potentials:
     *   A_{g,x} = -0.5 B0 (y - bg_y0),  A_{g,y} = +0.5 B0 (x - bg_x0).
     * Stage 20 unit-isolation test for the v x B_g gravitomagnetic Lorentz
     * force (gr_sandbox_v35.tex eq:geodesic_expansion). */
    float        bg_B0;
    /* Uniform EM magnetic field B_z, used by GR_BG_KIND_UNIFORM_MAGNETIC.
     * Symmetric-gauge potentials:
     *   A_x = -0.5 B0_em (y - bg_y0),  A_y = +0.5 B0_em (x - bg_x0).
     * Stage 23 unit-isolation test for the q v x B EM Lorentz force. */
    float        bg_B0_em;
    /* Uniform EM electric field (E_x, E_y), used by GR_BG_KIND_UNIFORM_ELECTRIC.
     *   phi^{bg}(x, y) = -( bg_Ex_em (x - bg_x0) + bg_Ey_em (y - bg_y0) )
     * Stage 24 unit-isolation test for the q E EM Lorentz force. */
    float        bg_Ex_em;
    float        bg_Ey_em;
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

/* TSC (W_3 / quadratic B-spline) deposit + interp on the CORNER
 * sublattice.  3x3 cell support per particle, smoother than CIC,
 * matched deposit/interp kernels preserve HE self-force adjoint at any
 * sub-cell position.  Used for rho_matter / rho_q in Tier-0 orbits to
 * suppress moving-particle deposit aliasing (the dominant PIC heating
 * mode) without resorting to ad-hoc binomial smoothing. */
void  gr_tsc_deposit_corner(float* arr, int W, int H, float dx,
                            float x_p, float y_p, float value);
float gr_tsc_interp_corner (const float* arr, int W, int H, float dx,
                            float x_p, float y_p);

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
 * 1 if the installed kind supplies a nonzero A_g (SPINNING_POINT_MASS or
 * UNIFORM_GRAVITOMAGNETIC), 0 otherwise.  When 0 is returned the caller
 * should treat A_g as zero. */
int gr_bg_eval_A_g(const struct gr_sim* sim, float x, float y,
                   float* Ax_out, float* Ay_out);

/* Defined in background.c — analytic-mode evaluation of the z-component
 * of the gravitomagnetic field B_g = curl(A_g) at (x, y) for the installed
 * background generator.  Returns 1 if a nonzero B_g_z is available, 0
 * otherwise.  Output:
 *   *Bgz_out : B_g_z(x, y) = d/dx A_{g,y} - d/dy A_{g,x}. */
int gr_bg_eval_B_g(const struct gr_sim* sim, float x, float y,
                   float* Bgz_out);

/* Defined in background.c — analytic-mode evaluation of the EM
 * vector potential A_em(x, y) for the installed background generator.
 * Returns 1 if the installed kind supplies a nonzero A_em (currently
 * only UNIFORM_MAGNETIC), 0 otherwise. */
int gr_bg_eval_A_em(const struct gr_sim* sim, float x, float y,
                    float* Ax_out, float* Ay_out);

/* Defined in background.c — analytic-mode evaluation of the z-component
 * of the EM magnetic field B = curl(A) at (x, y).  Returns 1 if a nonzero
 * B_z is available, 0 otherwise.  Output:
 *   *Bz_out : B_z(x, y) = d/dx A_y - d/dy A_x. */
int gr_bg_eval_B_em(const struct gr_sim* sim, float x, float y,
                    float* Bz_out);

/* Defined in background.c — analytic-mode evaluation of the EM scalar
 * potential phi(x, y) and its spatial gradient, for the installed
 * background generator.  Returns 1 if a nonzero phi_em is available,
 * 0 otherwise.  Output:
 *   *phi_out : phi^{bg}(x, y)
 *   *gx_out  : d/dx phi^{bg}(x, y)
 *   *gy_out  : d/dy phi^{bg}(x, y)
 * For UNIFORM_ELECTRIC: phi = -( Ex (x-x0) + Ey (y-y0) ), grad = (-Ex, -Ey). */
int gr_bg_eval_phi_em(const struct gr_sim* sim, float x, float y,
                      float* phi_out, float* gx_out, float* gy_out);

#endif /* GRLITE_SIM_INTERNAL_H */
