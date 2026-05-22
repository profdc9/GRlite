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
 * code may rely on PHI_GRAV == 0.  These six numeric values also coincide
 * with the first six values of gr_array_id_t below, so a cast is safe. */
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
 * Yee staggered grid layout (v35; specified in §9 of gr_sandbox_vNN.tex but
 * not implemented until v35).
 *
 * Each of the 12 stored grid arrays — 6 potentials and 6 sources — lives at
 * one of three sub-cell positions within each cell-(i,j) box:
 *
 *   GR_LATTICE_CORNER  : sample at (i, j) * dx
 *     Used by scalars: Phi_g, phi_em, rho_matter, rho_q.
 *
 *   GR_LATTICE_X_EDGE  : sample at (i+0.5, j) * dx
 *     Used by x-components: A_{g,x}, A_x, J_{m,x}, J_{q,x}, and the derived
 *     gradient dPhi/dx.
 *
 *   GR_LATTICE_Y_EDGE  : sample at (i, j+0.5) * dx
 *     Used by y-components: A_{g,y}, A_y, J_{m,y}, J_{q,y}, and the derived
 *     gradient dPhi/dy.
 *
 * The cell-center sublattice (i+0.5, j+0.5) * dx is the natural home for
 * derived curl quantities (B_{g,z}, B_z) but they are not stored — computed
 * on demand at particle positions or full-grid on request for visualization.
 *
 * All 12 arrays are allocated identically as N_x * N_y * sizeof(float).  For
 * X_EDGE fields the column i = N_x-1 is a "ghost" column kept at zero (and
 * never written by the leapfrog); for Y_EDGE fields the row j = N_y-1 is the
 * ghost.  CORNER fields use all cells.  This uniform allocation (rather than
 * shape-per-sublattice) preserves cache locality and is GPU/shader friendly.
 *
 * In v35, the implementation is being migrated from a single cell-centered
 * layout to this Yee specification.  Helpers gr_array_lattice() and
 * gr_sim_array_ptr() let callers discover and access arrays generically.
 * --------------------------------------------------------------------------*/

typedef enum {
    GR_LATTICE_CORNER = 0,  /* sample at (i,     j  ) * dx  — scalars     */
    GR_LATTICE_X_EDGE = 1,  /* sample at (i+0.5, j  ) * dx  — x-vectors    */
    GR_LATTICE_Y_EDGE = 2   /* sample at (i,     j+0.5) * dx — y-vectors  */
} gr_lattice_t;

/* Unified identifier for all 12 stored grid arrays — six potentials + six
 * sources.  The first six values match gr_field_id_t numerically (so e.g.
 * GR_ARR_PHI_GRAV == GR_FIELD_PHI_GRAV == 0), allowing safe casting when
 * code that knew only about potentials needs to operate on the broader
 * array list. */
typedef enum {
    /* Potentials (same numeric order as gr_field_id_t) */
    GR_ARR_PHI_GRAV   = 0,   /* corner   */
    GR_ARR_A_GX       = 1,   /* x-edge   */
    GR_ARR_A_GY       = 2,   /* y-edge   */
    GR_ARR_PHI_EM     = 3,   /* corner   */
    GR_ARR_A_X        = 4,   /* x-edge   */
    GR_ARR_A_Y        = 5,   /* y-edge   */
    /* Sources */
    GR_ARR_RHO_MATTER = 6,   /* corner   */
    GR_ARR_J_MX       = 7,   /* x-edge   */
    GR_ARR_J_MY       = 8,   /* y-edge   */
    GR_ARR_RHO_Q      = 9,   /* corner   */
    GR_ARR_J_QX       = 10,  /* x-edge   */
    GR_ARR_J_QY       = 11,  /* y-edge   */
    GR_ARR_COUNT      = 12
} gr_array_id_t;

/* Lattice classification for a given array id.  Constant table lookup. */
gr_lattice_t gr_array_lattice(gr_array_id_t which);

/* Sub-cell offset (in units of dx) of a sublattice's nodes relative to the
 * cell's lower-left corner.  For example, X_EDGE has offset (0.5, 0.0).
 * Useful for translating between array index and physical position. */
void gr_lattice_offset(gr_lattice_t lat, float* dx_out, float* dy_out);

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

/* Unified pointer accessor (v35) — returns the storage array for any of the
 * 12 stored grids (6 potentials + 6 sources).  Use gr_array_lattice(which)
 * to determine the sublattice the returned array's indices map to.  This is
 * the recommended generic access path; the dedicated gr_sim_rho_matter_ptr,
 * gr_sim_J_mx_ptr, ... functions kept for backward compatibility wrap this. */
float* gr_sim_array_ptr(gr_sim_t* sim, gr_array_id_t which);

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
    /* Stage 9: accumulated proper time along the particle's worldline.
     * Updated each gr_sim_step from
     *   d_tau = dt * sqrt(1 + 2 Phi/c^2 - (1 - 2 Phi/c^2) v^2/c^2
     *                       - 8 (v . A_g) / c^2)
     * with Phi, A_g evaluated at the particle's current position.
     * Initialized to 0 in gr_sim_add_particle. */
    float proper_time;
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

/* Gravitomagnetic Lorentz force gate.  When enabled (default), the particle
 * pusher includes the Tier-1 gravitomagnetic Lorentz piece
 *     F_gm = +4 m (v x B_g)
 * per gr_sandbox_v35.tex eq:geodesic_expansion (line 938).  When disabled,
 * the force reduces to the scalar Tier-0 / Tier-2 EIH form only — useful
 * for isolation tests of the proper-time clock effect (stage09), which is
 * derived under the assumption that prograde and retrograde orbits remain
 * on the same Keplerian circle.  With the gravitomagnetic force enabled,
 * the two orbits dynamically differ (Lense-Thirring frame-dragging on the
 * orbit), so the clock-effect formula picks up an additional orbit-shape
 * contribution.  Tests that want the orbit-shape effect (stage20+) leave
 * the flag at its default (enabled); tests that want clock-only isolate
 * by disabling it. */
void gr_sim_set_gravitomagnetic_force_enabled(gr_sim_t* sim, int enabled);
int  gr_sim_get_gravitomagnetic_force_enabled(const gr_sim_t* sim);

/* Gravitomagnetic INDUCTIVE force gate -- the -m*d_t A_g piece per the
 * doc's Tier-3 force law (gr_sandbox_v35.tex eq:eih_full / sec:alg_rel
 * eqbox line 1040).  Default OFF for backward compatibility with the
 * existing gravity-side test stages (which were calibrated without this
 * piece).  Symmetric counterpart to gr_sim_set_em_inductive_enabled.
 *
 * Provided primarily as a diagnostic: Stage 28 turns this on to test
 * whether the inductive heating observed for EM in Stage 27 also
 * appears for gravity.  Uses the same centered-time-difference
 * (A_g^{n+1} - A_g^{n-1})/(2 dt) at t^n. */
void gr_sim_set_gravitomagnetic_inductive_enabled(gr_sim_t* sim, int enabled);
int  gr_sim_get_gravitomagnetic_inductive_enabled(const gr_sim_t* sim);

/* EM Lorentz force gate (analog of the gravitomagnetic gate above).  When
 * enabled (default), the particle pusher includes the EM Lorentz force
 *     F_em = q ( -grad phi - d_t A + v x B )                  (B = curl A)
 * on each charged particle, with phi, A read from the background (analytic
 * or sampled) plus the perturbation fields (always sampled).  The factor
 * on v x B is +1 (no spin-2 enhancement, unlike the GEM coefficient +4).
 * See gr_sandbox_v35.tex eq:eih_full / sec:alg_rel Tier-3 eqbox.
 *
 * Currently implemented pieces (Stage 23+):
 *   v x B   - via the Yee curl of A_x, A_y arrays (mirror of the GM path).
 *
 * Not yet wired (planned):
 *   -grad phi_em  - the scalar Coulomb force (next stage).
 *   -d_t A        - the inductive electric field from time-varying A.
 *
 * The flag exists so future tests can isolate clock effects analogous to
 * Stage 9 Phase B (where the GM Lorentz force is gated off). */
void gr_sim_set_em_lorentz_force_enabled(gr_sim_t* sim, int enabled);
int  gr_sim_get_em_lorentz_force_enabled(const gr_sim_t* sim);

/* Fine-grained gates on the three EM Lorentz pieces — useful for
 * isolating which piece (if any) drives PIC self-heating in
 * closed-loop dynamics tests (Stage 27+).  All default ON. */
void gr_sim_set_em_inductive_enabled(gr_sim_t* sim, int enabled);
int  gr_sim_get_em_inductive_enabled(const gr_sim_t* sim);

/* Discretization for the inductive piece -q d_t A.  Two choices, both
 * exposed for diagnostic comparison:
 *
 *   GR_INDUCTIVE_CENTERED (default; 2nd-order accurate):
 *     (curr - next) / (2 dt) = (A^{n+1} - A^{n-1}) / (2 dt)
 *     Time-CENTER at t^n, but reads A^{n+1} which contains the particle's
 *     OWN current-step deposit.  This makes the self-deposit "visibility"
 *     inconsistent with the spatial-derivative reads (which use .prev =
 *     A^n and so don't contain the current-step deposit), and empirically
 *     causes PIC self-heating in closed-loop dynamics (Stage 27).
 *
 *   GR_INDUCTIVE_BACKWARD (1st-order accurate; closed-loop stable):
 *     (prev - next) / dt = (A^n - A^{n-1}) / dt
 *     Time-CENTER at t^{n-1/2}, but BOTH slices are pre-current-step --
 *     neither A^n nor A^{n-1} contains this step's deposit.  This matches
 *     the .prev convention of the spatial-derivative reads.  The
 *     half-step time mismatch is the standard read-prev-style trade
 *     (cf. Stage 18 read-prev fix for the gravity scalar gradient). */
typedef enum {
    GR_INDUCTIVE_CENTERED = 0,
    GR_INDUCTIVE_BACKWARD = 1
} gr_inductive_disc_t;
void                gr_sim_set_em_inductive_disc(gr_sim_t* sim, gr_inductive_disc_t kind);
gr_inductive_disc_t gr_sim_get_em_inductive_disc(const gr_sim_t* sim);

/* Sign of the inductive piece in the force law -- diagnostic.  Default
 * is +1.0 (the variationally-derived sign of -q d_t A).  Setting to -1.0
 * flips the sign as a probe for whether the closed-loop PIC heating
 * observed in Stage 27/28 is consistent with a sign convention error.
 * Symmetric pair for EM (-q d_t A) and gravity (-m d_t A_g). */
void  gr_sim_set_em_inductive_sign(gr_sim_t* sim, float sign);
float gr_sim_get_em_inductive_sign(const gr_sim_t* sim);
void  gr_sim_set_gravitomagnetic_inductive_sign(gr_sim_t* sim, float sign);
float gr_sim_get_gravitomagnetic_inductive_sign(const gr_sim_t* sim);

/* J time-centering correction (Stage 27/28 diagnostic).
 *
 * Esirkepov produces J^{n-1/2} satisfying discrete continuity with rho at
 * integer steps.  But the wave equation update
 *   A^{n+1} = 2 A^n - A^{n-1} + (c dt)^2 (Lap A^n + sc * src)
 * has its operator centered at integer step n, so the natural source is
 * J^n -- not the half-step-offset J^{n-1/2} we actually feed it.  This
 * is a half-step time-staggering inconsistency baked into the
 * potential-based scheme (standard Yee+Boris avoids it by storing E and
 * B independently with half-step interleaving).
 *
 * When this flag is on, the wave equation source is replaced by a
 * linear-extrapolation estimate of J at integer time:
 *   J^n_approx = 1.5 J^{n-1/2} - 0.5 J^{n-3/2}
 * which is 2nd-order accurate in dt.  Default OFF -- existing tests are
 * calibrated to the raw half-step-staggered behavior. */
void gr_sim_set_j_time_correction_enabled(gr_sim_t* sim, int enabled);
int  gr_sim_get_j_time_correction_enabled(const gr_sim_t* sim);

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

/* Stage 10 — enable automatic per-step deposition of every particle's mass
 * and current (rho_matter, J_mx, J_my) and, if charged, charge and current
 * (rho_q, J_qx, J_qy) onto the grid before each gr_sim_step's leapfrog runs.
 * When enabled, gr_sim_step performs the equivalent of:
 *
 *   gr_sim_clear_sources(sim);
 *   for each particle p:
 *     gr_sim_deposit_point_particle(sim, p->x, p->y, p->mass, p->charge,
 *                                   v_x(p), v_y(p));
 *   ... then the field leapfrog runs as usual.
 *
 * Velocity is read from p^{n-1/2} (lagged half-step), the same convention
 * as the force evaluation.  Default: disabled.  Stage 10+ scenarios turn
 * this on; earlier stages either deposit manually (Stage 5) or skip
 * deposition entirely (Stage 7/8/9, fixed background). */
void gr_sim_set_particle_source_deposition(gr_sim_t* sim, int enabled);
int  gr_sim_get_particle_source_deposition(const gr_sim_t* sim);

/* Esirkepov current deposition (v35 §sec:yee_pivot, answer B1).
 * Default-on; the off-flag is for regression testing only — disabling it
 * reintroduces PIC grid-heating artifacts on moving sources because direct
 * CIC violates the discrete continuity equation
 *   (rho^{n+1} - rho^n)/dt + div(J^{n-1/2}) = 0
 * by exactly v.grad(rho).  gr_sim_esirkepov_violations() reports the count
 * of timesteps in which a particle's motion exceeded 1 cell in x or y and
 * the deposit fell back to direct CIC; a nonzero value warrants caution. */
void gr_sim_set_esirkepov_enabled(gr_sim_t* sim, int enabled);
int  gr_sim_get_esirkepov_enabled(const gr_sim_t* sim);
int  gr_sim_esirkepov_violations(const gr_sim_t* sim);

/* Binomial-smoothing passes applied to rho_matter and rho_q after deposit
 * and before the field leapfrog reads them.  Each pass is the canonical
 * 3x3 [[1,2,1],[2,4,2],[1,2,1]]/16 PIC noise-reduction filter.  Reduces
 * high-frequency aliasing in the moving-particle deposit that the
 * wave-equation leapfrog would otherwise convert into a self-force wake
 * (the dominant PIC heating mode for Tier-0 gravitational orbits).
 * Default 0 (off).  Typical values 1-4. */
void gr_sim_set_rho_smooth_passes(gr_sim_t* sim, int passes);
int  gr_sim_get_rho_smooth_passes(const gr_sim_t* sim);

/* Shape function for rho deposit + force interp (matched pair preserves
 * Hockney-Eastwood adjoint condition).
 *   GR_SHAPE_CIC: linear bilinear (W_2), 2x2 cell footprint. Default.
 *   GR_SHAPE_TSC: quadratic B-spline (W_3), 3x3 cell footprint.  Smoother;
 *     reduces moving-particle deposit aliasing further than CIC + binomial
 *     smoothing.  Cost: ~2.25x more cells per deposit/interp. */
typedef enum {
    GR_SHAPE_CIC = 0,
    GR_SHAPE_TSC = 1
} gr_shape_function_t;
void                gr_sim_set_shape_function(gr_sim_t* sim, gr_shape_function_t s);
gr_shape_function_t gr_sim_get_shape_function(const gr_sim_t* sim);

/* Force-interpolation scheme — selects how the gradient of Phi at the
 * particle is computed from the grid Phi.  Both schemes satisfy the
 * Hockney-Eastwood adjoint condition (F_self = 0 at any stationary
 * sub-cell position).  They differ in their behavior on MOVING particles:
 *
 *   GR_FORCE_INTERP_LEGACY (default):
 *     Compute the FD gradient at each grid corner in the deposit kernel's
 *     footprint, then interpolate that grid gradient to the particle
 *     using the same kernel (W_2 for CIC, W_3 for TSC).
 *     This is "FD-then-interp".  Stage 10 Phase C validates F_self = 0
 *     under this scheme at any sub-cell position (memory:
 *     [[grlite-he-adjoint-translation-invariance]]).
 *
 *   GR_FORCE_INTERP_LEWIS_BIRDSALL:
 *     Compute the gradient at the particle directly as
 *         dPhi/dx_p = sum_g (dW/dx)(x_p - x_g) * Phi_g
 *     using the analytic gradient of the deposit kernel.  This is the
 *     variationally-consistent force from the discrete Lagrangian
 *     (Lewis 1970; Birdsall-Langdon Ch. 14; Brackbill-Forslund).
 *     The discrete work-energy theorem holds exactly, so total energy
 *     (kinetic + interaction + field) is conserved up to round-off,
 *     dramatically reducing moving-particle "PIC self-heating".  Targets
 *     the dominant Phase E heating mechanism that the absorbing boundary
 *     does NOT touch (see [[grlite-damping-sweep-result]]).
 *
 * Default GR_FORCE_INTERP_LEGACY.  All existing scenarios behave
 * identically under the default. */
typedef enum {
    GR_FORCE_INTERP_LEGACY         = 0,
    GR_FORCE_INTERP_LEWIS_BIRDSALL = 1
} gr_force_interp_t;
void              gr_sim_set_force_interp(gr_sim_t* sim, gr_force_interp_t scheme);
gr_force_interp_t gr_sim_get_force_interp(const gr_sim_t* sim);

/* Periodic BC for the field leapfrog: 1 = on, 0 = zero-Dirichlet (default).
 * Periodic BC restores translation invariance of the discrete Laplacian
 * Green's function, which is the precondition for HE self-force = 0 at
 * any particle position (not just the box symmetry center).  Note that
 * outgoing waves wrap around to the opposite edge under periodic BC,
 * which without damping pollutes the simulation; combine with the
 * damping ring or PML when running production simulations. */
void gr_sim_set_periodic_bc(gr_sim_t* sim, int periodic);
int  gr_sim_get_periodic_bc(const gr_sim_t* sim);

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

/* Softened spinning point mass — fills Phi_g^{bg} as in the non-spinning
 * case and also fills A_{g,x}^{bg}, A_{g,y}^{bg} with the gravitomagnetic
 * dipole field for a point spin J_z along +z:
 *
 *   A_g(x) = (G_eff / (2 c^2)) * J × r / (r^2 + epsilon^2)^{3/2}
 *
 * which in 2D with the spin axis perpendicular to the simulation plane
 * reduces to
 *
 *   A_{g,x}(x,y) = -(G_eff J_z / (2 c^2)) * (y - y0) / s^{3/2},
 *   A_{g,y}(x,y) = +(G_eff J_z / (2 c^2)) * (x - x0) / s^{3/2},
 *   s = (x-x0)^2 + (y-y0)^2 + epsilon^2.
 *
 * Stores the kind = SPINNING_POINT_MASS along with (x0, y0, GM, epsilon, J_z)
 * for the analytic-mode evaluator.  Calling this with J_z = 0 is equivalent
 * to gr_sim_set_background_point_mass.  Overrides any previous background. */
void gr_sim_set_background_spinning_point_mass(gr_sim_t* sim,
                                               float x0, float y0,
                                               float GM, float epsilon,
                                               float Jz);

/* Uniform gravitomagnetic field — fills A_{g,x}^{bg}, A_{g,y}^{bg} with the
 * symmetric-gauge potentials that produce a spatially constant B_g_z = B_0:
 *
 *   Phi_g^{bg} = 0                  (no scalar gravity)
 *   A_{g,x}    = -0.5 * B_0 * (y - y_0)
 *   A_{g,y}    = +0.5 * B_0 * (x - x_0)
 *   B_g_z      = d/dx A_{g,y} - d/dy A_{g,x} = B_0   (uniform)
 *
 * The (x_0, y_0) origin only sets the gauge zero of A_g and is otherwise
 * physically irrelevant.  Intended for the Stage 20 unit-isolation test of
 * the gravitomagnetic Lorentz force (cf. gr_sandbox_v35.tex
 * §sec:geodesic_expansion eq:geodesic_expansion): with Phi_g = 0 the force
 * on a particle reduces to F = m * 4 * v x B_g, so the particle traces a
 * circle at the gyrofrequency omega_gm = 4 |B_g| / gamma. */
void gr_sim_set_background_uniform_gravitomagnetic(gr_sim_t* sim,
                                                   float x0, float y0,
                                                   float B0);

/* Uniform electromagnetic magnetic field — EM analog of the gravitomagnetic
 * variant above.  Fills A_x^{bg}, A_y^{bg} with the symmetric-gauge
 * potentials producing a spatially constant magnetic field B_z = B_0:
 *
 *   phi^{bg}   = 0                  (no scalar EM)
 *   A_x        = -0.5 * B_0 * (y - y_0)
 *   A_y        = +0.5 * B_0 * (x - x_0)
 *   B_z        = d/dx A_y - d/dy A_x = B_0  (uniform, by construction)
 *
 * Stage 23 cyclotron test for the EM Lorentz force F = q (v x B); a
 * charged particle traces a circle at omega_c = q |B_0| / (gamma m).
 * Coefficient is +1 on v x B (no spin-2 factor, unlike gravity).  See
 * gr_sandbox_v35.tex eq:eih_full / sec:alg_rel Tier-3 eqbox line 1045. */
void gr_sim_set_background_uniform_magnetic(gr_sim_t* sim,
                                            float x0, float y0,
                                            float B0);

/* Uniform electric field — fills phi^{bg} (EM scalar potential) so the
 * gradient gives a spatially constant (E_x, E_y):
 *
 *   phi^{bg}(x, y)   = -( E_x (x - x_0) + E_y (y - y_0) )
 *   -grad phi^{bg}   = ( E_x, E_y )                       (uniform)
 *   A^{bg}           = 0                                  (no magnetic part)
 *
 * Stage 24 unit-isolation test for the q E piece of the EM Lorentz force.
 * A charged particle at rest accelerates uniformly along the E direction:
 * p(t) = q E t  (relativistic momentum form). */
void gr_sim_set_background_uniform_electric(gr_sim_t* sim,
                                            float x0, float y0,
                                            float Ex, float Ey);

/* Softened point charge — EM analog of gr_sim_set_background_point_mass.
 * Fills phi^{bg} (EM scalar) with the Coulomb potential
 *
 *   phi^{bg}(x) = +k_e * Q / sqrt(|x - x0|^2 + epsilon^2)
 *
 * For a test particle of charge q_test, the resulting force is
 *   F = -q_test * grad phi^{bg} = +q_test * k_e * Q * (r - r0) / s^{3/2}.
 * Like-sign Q and q_test repel; opposite-sign attract.  Stage 26 EM
 * Kepler analog uses Q > 0 and q_test < 0 to get an attractive orbit
 * structurally identical to gravity (gravitational coupling
 * G_eff M_mass <-> EM coupling |q_test Q| k_e / m_test). */
void gr_sim_set_background_point_charge(gr_sim_t* sim,
                                        float x0, float y0,
                                        float Q, float epsilon);

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
    GR_BG_KIND_NONE                    = 0,
    GR_BG_KIND_POINT_MASS              = 1,
    GR_BG_KIND_SPINNING_POINT_MASS     = 2,
    GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC = 3,
    GR_BG_KIND_UNIFORM_MAGNETIC        = 4,
    GR_BG_KIND_UNIFORM_ELECTRIC        = 5,
    GR_BG_KIND_POINT_CHARGE            = 6
    /* Future: CHARGED_POINT_MASS, KERR_NEWMAN, ... */
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
 * Multiplicative (lossy-material) absorber of n_damping cells on each edge.
 * The original spec — gr_sandbox_v32.tex §9.6 sec:abc eq:damp_profile —
 * prescribes a quadratic profile sigma(d) = sigma_max * (d/L)^2 with
 * sigma_max = 21 c / (2 L), targeting round-trip reflection R ~ 1e-3.
 * That is the literature-default polynomial PML profile evaluated at
 * m=2, R=1e-3, and remains the default of gr_sim_set_damping().
 *
 * The gr_sim_set_damping_config() entry point exposes the literature's
 * parametrized profile family — polynomial m=1..4+, exponential beta=...
 * — for systematic experimentation.  See Berenger (1994/1996),
 * Roden & Gedney (2000) for the canonical derivations; the textbook
 * formulas are
 *   POLYNOMIAL:  sigma(d/L) = sigma_max * (d/L)^m,
 *                sigma_max  = -(m+1) c ln(R) / (2 L)
 *   EXPONENTIAL: sigma(d/L) = sigma_max * (e^(beta d/L) - 1) / (e^beta - 1),
 *                sigma_max  = -c ln(R) / (2 L * J(beta))
 *                where J(beta) = 1/beta - 1/(e^beta - 1).
 *
 * Pass n_damping = 0 to disable.  Mutually exclusive with PML (when added). */

typedef enum {
    GR_DAMP_POLYNOMIAL  = 0,  /* sigma(d/L) = sigma_max * (d/L)^poly_order */
    GR_DAMP_EXPONENTIAL = 1   /* sigma(d/L) = sigma_max * (e^(beta*d/L)-1)/(e^beta-1) */
} gr_damp_profile_kind_t;

/* Time-discretization form for the damping kernel.  Treating the
 * leapfrog as a 2nd-order recurrence with characteristic polynomial
 * w^2 + a w + b = 0, each form chooses (a, b) differently:
 *
 *   MULTIPLICATIVE (default): Phi^{n+1} = (1 - sigma dt) * [2 Phi^n -
 *     Phi^{n-1} + (c dt)^2 (Lap Phi + S)].  Roots are complex
 *     conjugates with |w| = sqrt(1 - sigma dt).  Per-step decay is
 *     sqrt(1 - sigma dt) — about half what the prefactor suggests.
 *     Stability bound sigma dt < 1.
 *
 *   CRITICAL: Phi^{n+1} = 2 gamma Phi^n - gamma^2 Phi^{n-1} +
 *     (c dt)^2 (Lap Phi + S), with gamma = 1 - sigma dt.  Both roots
 *     of the characteristic polynomial sit at gamma (double real
 *     root — critically damped).  Per-step decay is exactly gamma.
 *     About 2x the damping rate of MULTIPLICATIVE for the same nominal
 *     sigma.  Continuum PDE is the damped Klein-Gordon equation with
 *     a (sigma)^2 mass term — fine inside an absorbing ring (which
 *     wants exponential screening), wrong for the interior (would
 *     screen the long-range gravity profile).  Use only when the
 *     damping array is nonzero only at the boundary. */
typedef enum {
    GR_DAMP_TIME_MULTIPLICATIVE = 0,
    GR_DAMP_TIME_CRITICAL       = 1
} gr_damp_time_form_t;

typedef struct {
    int                    n_damping;          /* layer thickness (cells); 0 = off */
    gr_damp_profile_kind_t kind;               /* default GR_DAMP_POLYNOMIAL */
    float                  poly_order;         /* m; default 2.0 (current/spec) */
    float                  exp_beta;           /* beta; default 4.0 */
    float                  target_reflection;  /* R for sigma_max derivation; default 1e-3 */
    float                  sigma_max_override; /* > 0 to bypass formula; 0 = use formula */
    gr_damp_time_form_t    time_form;          /* default GR_DAMP_TIME_MULTIPLICATIVE */
} gr_damp_config_t;

/* Legacy entry point — equivalent to set_damping_config with kind=POLYNOMIAL,
 * poly_order=2, target_reflection=1e-3.  Bit-exact backward compat with the
 * §9.6 spec default. */
void gr_sim_set_damping(gr_sim_t* sim, int n_damping);
int  gr_sim_damping_layers(const gr_sim_t* sim);

/* Parametrized damping setup.  See gr_damp_config_t comments above for
 * the supported profile family.  Pass cfg->n_damping = 0 to disable.
 *
 * SAFETY: when cfg->time_form is GR_DAMP_TIME_CRITICAL, the kernel has a
 * tighter CFL bound than the undamped leapfrog (the Nyquist mode goes
 * unstable for any gamma < 1 at CFL = 1/sqrt(2) in 2D).  The
 * implementation enforces the stability bound by clamping sigma_max so
 * that gamma_min = 1 - sigma_max*dt >= 2*CFL^2.  If the requested
 * sigma_max would exceed that, it is clamped silently; readback via
 * gr_sim_get_damping_config returns the clamped value.  Use
 * gr_sim_damping_max_stable_sigma_dt() to query the ceiling before
 * calling. */
void             gr_sim_set_damping_config(gr_sim_t* sim, const gr_damp_config_t* cfg);
gr_damp_config_t gr_sim_get_damping_config(const gr_sim_t* sim);

/* Returns the largest sigma*dt that critical damping can use without
 * destabilizing the Nyquist mode at the current sim's CFL.  Formula:
 *   sigma*dt_max = 1 - 2*CFL^2     (2D, 5-point Laplacian)
 * For the multiplicative form this returns 1.0 (the (1-sigma*dt) > 0
 * bound), since multiplicative is more stable than the undamped leapfrog
 * at the Nyquist mode. */
float gr_sim_damping_max_stable_sigma_dt(const gr_sim_t* sim, gr_damp_time_form_t form);

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
