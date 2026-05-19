/* GRlite public C API.
 *
 * This header is the surface that both the native test binary and the JS/WASM
 * frontend bind against. Internal struct definitions live in src/sim_internal.h
 * and are not visible here.
 *
 * Spec reference: gr_sandbox_v32.tex §12 "Implementation Plan" — stages 1+.
 */
#ifndef GRLITE_H
#define GRLITE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque simulation handle. */
typedef struct gr_sim gr_sim_t;

/* Field identifiers — six potentials per v33 §12.4 / v32 §9.1. The Stage 4
 * refactor expanded this from {Phi_g} to all six. Numeric values are stable;
 * code may rely on PHI_GRAV == 0. */
typedef enum {
    GR_FIELD_PHI_GRAV = 0,  /* gravitational scalar  Phi_g     */
    GR_FIELD_A_GX     = 1,  /* gravitomagnetic vec.  A_{g,x}   */
    GR_FIELD_A_GY     = 2,  /* gravitomagnetic vec.  A_{g,y}   */
    GR_FIELD_PHI_EM   = 3,  /* electric scalar       phi       */
    GR_FIELD_A_X      = 4,  /* magnetic vector       A_x       */
    GR_FIELD_A_Y      = 5,  /* magnetic vector       A_y       */
    GR_FIELD_COUNT    = 6
} gr_field_id_t;

/* ----------------------------------------------------------------------------
 * Simulation lifecycle
 * --------------------------------------------------------------------------*/

/* Create a simulation on a width-by-height grid with cell size dx (simulation
 * units) and effective speed-of-light c_eff. The timestep is dt = cfl * dx / c_eff.
 *
 * CFL stability: gr_sandbox_v32.tex §9.2 eq:cfl requires cfl <= 1/sqrt(d), so
 * cfl <= 1/sqrt(2) ~= 0.7071 in 2D. Values above will diverge.
 *
 * Returns NULL on invalid parameters or allocation failure. G_eff defaults to
 * 1.0 (set via gr_sim_set_G_eff). */
gr_sim_t* gr_sim_create(int width, int height, float dx, float c_eff, float cfl);

/* Set effective coupling constants (Stage 3 / Stage 4).
 *
 *   Phi_g source coeff:  -4*pi*G_eff
 *   A_g_{x,y} source coeff: -4*pi*G_eff / c_eff^2
 *   phi_em source coeff: -4*pi*k_e
 *   A_{x,y} source coeff: -4*pi*k_e / c_eff^2
 *
 * Both default to 1.0 at gr_sim_create. Setting one updates the corresponding
 * field source coefficients in-place. */
void  gr_sim_set_G_eff(gr_sim_t* sim, float G_eff);
float gr_sim_get_G_eff(const gr_sim_t* sim);
void  gr_sim_set_k_e(gr_sim_t* sim, float k_e);
float gr_sim_get_k_e(const gr_sim_t* sim);

/* Free the simulation. Safe to pass NULL. */
void gr_sim_destroy(gr_sim_t* sim);

/* Advance one timestep. */
void gr_sim_step(gr_sim_t* sim);

/* Advance n timesteps (convenience for tight inner-loops invoked from JS). */
void gr_sim_step_n(gr_sim_t* sim, int n);

/* ----------------------------------------------------------------------------
 * State queries
 * --------------------------------------------------------------------------*/

int   gr_sim_step_count(const gr_sim_t* sim);
float gr_sim_time(const gr_sim_t* sim);  /* simulation time = step_count * dt */
float gr_sim_dt(const gr_sim_t* sim);
float gr_sim_dx(const gr_sim_t* sim);
int   gr_sim_width(const gr_sim_t* sim);
int   gr_sim_height(const gr_sim_t* sim);

/* Pointer to the current (time-n) field array, row-major, size width*height.
 * The pointer remains valid until gr_sim_destroy or the next gr_sim_step (the
 * leapfrog rotates three internal buffers — see src/sim_internal.h). */
float* gr_sim_field_ptr(gr_sim_t* sim, gr_field_id_t which);

/* ----------------------------------------------------------------------------
 * Static source deposition (Stage 3)
 *
 * Sources for the wave-equation RHS are stored on the grid (same staggering
 * as Phi). A scenario builds the source distribution at setup; the leapfrog
 * reads it but never modifies it (Stage 3 == static sources). The discrete
 * Poisson convergence condition at static equilibrium is, for the Phi_g
 * equation,
 *    Lap Phi_g = 4*pi*G_eff * rho_matter
 * which is the v32 §9.7 leapfrog form (-4*pi*G_eff*rho) evaluated when the
 * field has stopped changing.
 *
 * gr_sim_deposit_point_mass implements eq:cic_deposit (§9.5) — bilinear
 * sub-cell weights, four-cell support, total deposited integral = mass.
 * --------------------------------------------------------------------------*/

void gr_sim_clear_sources(gr_sim_t* sim);
void gr_sim_deposit_point_mass(gr_sim_t* sim, float x, float y, float mass);
void gr_sim_deposit_point_charge(gr_sim_t* sim, float x, float y, float charge);

/* Stage 5 — single-call deposit of a point "particle" with mass, charge, and
 * velocity. Lays down all six source contributions at the same CIC stencil:
 *   rho_matter += mass        rho_q  += charge
 *   J_mx       += mass*vx     J_qx   += charge*vx
 *   J_my       += mass*vy     J_qy   += charge*vy */
void gr_sim_deposit_point_particle(gr_sim_t* sim, float x, float y,
                                   float mass, float charge,
                                   float vx, float vy);

/* Read access to all six source arrays. Always non-NULL after gr_sim_create. */
const float* gr_sim_rho_matter_ptr(const gr_sim_t* sim);
const float* gr_sim_rho_q_ptr(const gr_sim_t* sim);
const float* gr_sim_J_mx_ptr(const gr_sim_t* sim);
const float* gr_sim_J_my_ptr(const gr_sim_t* sim);
const float* gr_sim_J_qx_ptr(const gr_sim_t* sim);
const float* gr_sim_J_qy_ptr(const gr_sim_t* sim);

/* ----------------------------------------------------------------------------
 * Lorenz gauge monitoring (Stage 4)
 *
 * Discrete Lorenz residuals (gr_sandbox_v32.tex §9.5 eq:lorenz_gem and
 * eq:lorenz_em):
 *    G_grav = (1/c^2) d_t Phi_g + div A_g          ; ideally 0
 *    G_em   = (1/c^2) d_t phi    + div A           ; ideally 0
 * In the discrete scheme the residual sits at the truncation-error level
 * (k Delta x)^2 (§9.5 monitoring paragraph) once the field has reached its
 * static or quasi-static state. These getters return the RMS over interior
 * cells. Time derivatives use a one-sided backward difference of the curr
 * vs prev rotation slot — call after gr_sim_step.
 * --------------------------------------------------------------------------*/

float gr_sim_gauge_residual_grav(const gr_sim_t* sim);
float gr_sim_gauge_residual_em(const gr_sim_t* sim);

/* ----------------------------------------------------------------------------
 * Particles (Stage 7)
 *
 * A particle carries a 2D position (x, y), 2D relativistic 3-momentum (px, py)
 * = gamma*m*v, a rest mass, and an electric charge. The integrator is the
 * relativistic Boris-leapfrog (kick-drift) form prescribed by §9.2 of v32:
 *
 *   p_{n+1/2}  = p_{n-1/2}  + F^n * dt
 *   gamma     = sqrt(1 + |p|^2 / (m c)^2)
 *   v_{n+1/2}  = p_{n+1/2} / (gamma * m)
 *   x_{n+1}   = x_n + v_{n+1/2} * dt
 *
 * Stages 7-8 (Boris in a background) use only the gravitational gradient
 * (eq:force_gem_pot reduced to -m*grad(Phi_g) when v=0 and A=0). Stage 10+
 * will add the velocity-dependent gravitomagnetic and EM Lorentz pieces.
 * --------------------------------------------------------------------------*/

typedef struct {
    float x, y;
    float px, py;
    float mass;
    float charge;
} gr_particle_t;

/* Force tier — selects which terms enter the gravitational force.
 * (Spec: gr_sandbox_v33.tex §"Practical implementation tiers", line 685.)
 *
 *   NEWTONIAN   — F = m*g, with g = -grad(Phi_g_total). Relativistic momentum
 *                 is still tracked. Valid v/c < 0.01. This is the Stage 7
 *                 behavior.
 *   RELATIVISTIC — Adds the Tier-2 O(v^2/c^2) terms from the geodesic
 *                  expansion (eq:geodesic_expansion §"Expansion of the
 *                  geodesic equation"):
 *                    F = m * [ g + (v^2/c^2)*g + 4*(v.g)*v/c^2 ]
 *                  Required for Stage 8 (perihelion precession). Valid up to
 *                  v/c ~ 0.5; full Tier-3 (with A_g terms) arrives at Stage 10.
 *
 * Default at gr_sim_create: NEWTONIAN. */
typedef enum {
    GR_FORCE_NEWTONIAN   = 0,
    GR_FORCE_RELATIVISTIC = 1
} gr_force_tier_t;

void            gr_sim_set_force_tier(gr_sim_t* sim, gr_force_tier_t tier);
gr_force_tier_t gr_sim_get_force_tier(const gr_sim_t* sim);

int  gr_sim_add_particle(gr_sim_t* sim, float x, float y,
                         float mass, float charge,
                         float vx, float vy);
int  gr_sim_particle_count(const gr_sim_t* sim);
const gr_particle_t* gr_sim_get_particle(const gr_sim_t* sim, int idx);
void gr_sim_clear_particles(gr_sim_t* sim);

/* Total energy of particle `idx`: gamma*m*c^2 + m*Phi_g_total(x_p).
 * Used as the conservation diagnostic in Stage 7+. Returns 0 if idx is out
 * of range. */
float gr_sim_particle_energy(const gr_sim_t* sim, int idx);

/* Enable/disable the per-step wave-equation leapfrog on the six perturbation
 * fields.  When all source arrays are zero and the perturbation fields are
 * zero (e.g., Stage 7/8 — a single test particle in a fixed background, no
 * deposition), the leapfrog step does no useful work and dominates the cost
 * of long-running orbital tests.  Default at gr_sim_create: enabled. */
void gr_sim_set_field_evolution(gr_sim_t* sim, int enabled);
int  gr_sim_get_field_evolution(const gr_sim_t* sim);

/* ----------------------------------------------------------------------------
 * Sampled background field arrays (Stage 6)
 *
 * Each of the six potentials has an optional companion "background" grid array,
 * filled once at scenario setup by a generator and never touched by the
 * leapfrog or damping layer. Force evaluations (Stage 7+) read the sum
 * Phi_total = Phi_bg + Phi_pert via a single CIC interpolation.
 *
 * The arrays are lazily allocated on first set_background_* call. Pass
 * gr_sim_clear_background to free them all (e.g. when switching scenarios).
 * The pointer returned by gr_sim_background_ptr is NULL until the
 * corresponding background has been set.
 *
 * See gr_sandbox_v33.tex §12.6 (sec:stage_bg) for the architecture and the
 * eq:bg_softened_point_mass generator below.
 * --------------------------------------------------------------------------*/

float* gr_sim_background_ptr(gr_sim_t* sim, gr_field_id_t which);
void   gr_sim_clear_background(gr_sim_t* sim);

/* Softened point mass — fills Phi_g^{bg} with
 *   -G*M / sqrt(|x - x0|^2 + epsilon^2)
 * sampled at cell centers. epsilon (in simulation length units) is a smoothing
 * length, recommended ~ few cells, that avoids the 1/r singularity at the
 * source. Calls clear of any previous Phi_g^{bg} and allocates fresh.
 *
 * Also stores the analytic parameters (kind = POINT_MASS, x0, y0, GM, eps)
 * so the background may be evaluated either through the cell-centered sample
 * (GR_BG_MODE_SAMPLED, default) or directly from the closed-form expression
 * at the particle's exact position (GR_BG_MODE_ANALYTIC). See
 * gr_sim_set_bg_mode below. */
void gr_sim_set_background_point_mass(gr_sim_t* sim,
                                      float x0, float y0,
                                      float GM, float epsilon);

/* ----------------------------------------------------------------------------
 * Background evaluation mode
 *
 * The sampled-grid path (SAMPLED, default) evaluates Phi_g^{bg} and its
 * gradient via CIC interpolation of the cell-centered array plus a centered
 * finite-difference. This is the path used in Stage 6's tests and the only
 * path that supports arbitrary user-supplied background fields, but it
 * introduces an O((dx/r)^2) discretization error: at a generic non-grid-
 * aligned particle position the FD gradient picks up a small tangential
 * component that breaks the exact spherical symmetry of the underlying
 * analytic field, producing a small spurious precession even in pure
 * Keplerian dynamics.
 *
 * The analytic path (ANALYTIC) evaluates the background generator (point
 * mass, eventually spinning/charged variants) at the particle's exact
 * position with no grid involvement. The perturbation field — when active
 * in Stage 10+ — still flows through the sampled CIC+FD path; the analytic
 * branch only replaces the background contribution.
 *
 * Default: SAMPLED. Switch with gr_sim_set_bg_mode. The choice is independent
 * of which background generator was last installed.
 * --------------------------------------------------------------------------*/

typedef enum {
    GR_BG_KIND_NONE       = 0,
    GR_BG_KIND_POINT_MASS = 1
    /* Future: SPINNING_POINT_MASS, CHARGED_POINT_MASS, KERR_NEWMAN, ... */
} gr_bg_kind_t;

typedef enum {
    GR_BG_MODE_SAMPLED  = 0,
    GR_BG_MODE_ANALYTIC = 1
} gr_bg_mode_t;

void         gr_sim_set_bg_mode(gr_sim_t* sim, gr_bg_mode_t mode);
gr_bg_mode_t gr_sim_get_bg_mode(const gr_sim_t* sim);
gr_bg_kind_t gr_sim_get_bg_kind(const gr_sim_t* sim);

/* ----------------------------------------------------------------------------
 * Absorbing damping layer (Stage 2)
 *
 * Install (or remove) a quadratic-profile absorbing layer of n_damping cells
 * on each grid edge. Per gr_sandbox_v32.tex §9.6 sec:abc eq:damp_profile:
 *   sigma(x) = sigma_max * (x/L)^2,  L = n_damping * dx,
 *   sigma_max = 21 * c_eff / (2 * L)  =>  round-trip reflection R ~ 1e-3.
 * Pass n_damping = 0 to disable (no per-cell multiply will occur).
 * Default after gr_sim_create: no damping (Stage 1 behavior).
 * --------------------------------------------------------------------------*/
void gr_sim_set_damping(gr_sim_t* sim, int n_damping);
int  gr_sim_damping_layers(const gr_sim_t* sim);

/* ----------------------------------------------------------------------------
 * Scenario registry
 *
 * Scenarios are the single source of truth shared between tests and the web
 * frontend (see memory: grlite-repo-layout). Each scenario lives in
 * core/scenarios/<name>.c and registers itself at startup via
 * gr_scenario_register, called from gr_scenarios_init.
 * --------------------------------------------------------------------------*/

typedef int (*gr_scenario_build_fn)(gr_sim_t* sim, const float* params, int n_params);

typedef struct {
    const char* name;
    gr_scenario_build_fn build;
} gr_scenario_t;

/* Load and apply a registered scenario by name. params is an array of float
 * parameters whose meaning is scenario-specific. Returns 0 on success,
 * negative on error (unknown name, NULL sim, bad params). */
int gr_sim_load_scenario(gr_sim_t* sim, const char* name, const float* params, int n_params);

/* Registry management — used by scenario .c files. */
void gr_scenario_register(const gr_scenario_t* s);
void gr_scenarios_init(void);  /* idempotent; called automatically by load_scenario */

#ifdef __cplusplus
}
#endif

#endif /* GRLITE_H */
