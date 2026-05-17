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

/* Field identifiers (used by gr_sim_field_ptr). */
typedef enum {
    GR_FIELD_PHI_GRAV = 0,  /* gravitational scalar potential Phi_g */
    GR_FIELD_COUNT
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
 * Returns NULL on invalid parameters or allocation failure.
 */
gr_sim_t* gr_sim_create(int width, int height, float dx, float c_eff, float cfl);

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
