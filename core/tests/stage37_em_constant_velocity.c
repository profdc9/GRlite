/* Stage 37 -- DIAGNOSTIC: constant-velocity charge in a long asymmetric domain.
 *
 * Physics statement: a free charge moving at constant velocity in vacuum does
 * NOT radiate (Larmor radiation requires acceleration).  Equivalently, the
 * self-force on a uniformly-moving charge in Maxwell's equations is zero.
 * In a faithful PIC scheme, the kinetic energy of such a particle must be
 * conserved up to grid-aliasing noise -- there is no physically-allowed
 * channel for KE to leak into the field.
 *
 * Diagnosis target: Stage 30 showed that turning on the -q d_t A inductive
 * piece in the EM Kepler orbit produces an OUTSPIRAL (energy GAIN), whereas
 * physical radiation reaction demands inspiral (energy LOSS).  Stage 36
 * confirmed that flipping the sign of -q d_t A restores inspiral.  This is
 * suspicious: a sign flip on a Lagrangian-derived term should NOT be
 * physically required.  The leading hypothesis is that the inductive piece
 * reads A at a TIME SLICE inconsistent with how -grad phi reads phi (or how
 * J was deposited to source A), so the discrete force law on a moving
 * particle has a spurious longitudinal component that does not vanish for
 * uniform motion.
 *
 * Test design: 128 (x) by 1024 (y) ASYMMETRIC domain.  Launch a single
 * charged particle at low y with constant +y velocity.  No background; the
 * only fields acting on the particle are its own self-deposited
 * (phi_em, A_y).  The asymmetric long axis lets the particle traverse many
 * Liénard-Wiechert wavelengths of its own wake without interacting with the
 * boundary; the damping ring absorbs the wake on each transverse edge.
 *
 * Configurations compared (at v_drift = 0.30 c):
 *   (a) full EM Lorentz, inductive sign +1   (production default)
 *   (b) full EM Lorentz, inductive sign -1   (the Stage 36 sign-flip)
 *   (c) inductive OFF                        (phi + v x B only)
 *   (d) all EM force OFF                     (pure control: 0% drift)
 *
 * Plus a velocity scan with the production config to expose the scaling of
 * the spurious self-force with v -- friction-like (~v), radiation-like
 * (~v^3), or some other power.
 *
 * Physics: in continuous Maxwell, the self-force on a uniformly-moving
 * charge is identically zero -- the -grad phi piece (longitudinal field of
 * the boosted Coulomb potential) is EXACTLY cancelled by the -d_t A piece
 * (from the time-changing A in the boosted frame).  In a PIC discretization
 * the two channels rarely cancel exactly; the residual mismatch becomes a
 * spurious self-force on uniform motion.  If (a) leaks energy but (b) or
 * (c) leaks less, that points to which discretization of -d_t A better
 * cancels -grad phi.  Config (d) measures the "uncancelled phi self-force"
 * baseline; the inductive piece is supposed to cancel a large chunk of it. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); \
        return 1; \
    } \
} while (0)

typedef struct {
    float KE_initial;
    float KE_final;
    float KE_max_excursion;   /* max |KE - KE_initial| over trajectory */
    float vy_initial;
    float vy_final;
    float y_initial;
    float y_final;
    int   n_steps_run;
    int   bailed_out;          /* nonzero if particle left the safe interior */
    /* Per-step spurious-force measurement.  After warmup, vy is sampled
     * over a short window (force_window_steps); accel_per_step is the
     * mean (vy_after_window - vy_after_warmup) / force_window_steps.
     * Reported in code units of (cell / dt^2) since dt = cfl*dx/c_eff is
     * dimensionless in this test.  Acceleration is force/mass.  At
     * non-zero mass it can be converted to force = mass * accel. */
    float accel_per_step;       /* dvy / step over post-warmup window */
    int   force_window_steps;   /* number of steps the window averaged over */
    float dt_used;              /* dt of this run, for converting to dvy/dt */
    float mass_used;             /* mass of the particle (for absolute-force conversion) */
} run_result_t;

static float kinetic_energy(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return (gamma - 1.0f) * m * c_eff * c_eff;
}

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static run_result_t run_config_full_mass(int electrostatic_on, int inductive_on, int magnetic_on,
                                    int lorentz_on, float ind_sign, float v_drift,
                                    float c_eff, float j_shift,
                                    int rho_smooth, int j_smooth,
                                    float cfl, float dx, float mass,
                                    const char* label) {
    /* Physical domain held fixed (128 x 1024 in physical units).  Cells
     * scale as 1/dx so the physical setup is identical at any dx.  This
     * makes the test a true joint spatial+temporal convergence study. */
    const float W_phys  = 128.0f;
    const float H_phys  = 1024.0f;
    const int   W       = (int)(W_phys / dx + 0.5f);
    const int   H       = (int)(H_phys / dx + 0.5f);
    const int   n_damp  = (int)(16.0f  / dx + 0.5f);
    const float Q       = 0.01f;
    /* y_start in PHYSICAL units = 16 (absorber) + 32 (buffer) = 48. */
    const float y_start = 16.0f + 32.0f;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;
    const float dt      = cfl * dx / c_eff;

    /* Drive the particle to y_end and stop before the trailing wake reaches
     * the absorber edge on the leading side.  Buffer = 64 cells at baseline
     * = 64 physical units.  Compute in physical units. */
    const float y_max_safe = H_phys - 16.0f - 64.0f;
    const float distance   = y_max_safe - y_start;
    /* At low v the traversal step count balloons (~1/v).  Cap at 20000 so
     * the v-scan stays bounded; the bail-out check in the loop handles
     * high v where the particle exits the safe region early. */
    int n_steps = (int)(distance / (v_drift * dt));
    if (n_steps > 20000) n_steps = 20000;
    if (n_steps < 100)   n_steps = 100;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) {
        run_result_t r = {0};
        r.bailed_out = 1;
        return r;
    }
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* TSC + LB throughout, with the new TSC X_EDGE / Y_EDGE interp now in
     * place for the inductive -q d_t A channel (Stage 37 follow-up:
     * matched-shape adjoint pairing on the A side). */
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, rho_smooth);
    gr_sim_set_j_smooth_passes(sim, j_smooth);

    /* Turn off ALL gravity so the only physics on the particle is EM.  */
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

    /* EM channel selection per config. */
    gr_sim_set_em_lorentz_force_enabled(sim, lorentz_on);
    gr_sim_set_em_electrostatic_enabled(sim, electrostatic_on);
    gr_sim_set_em_inductive_enabled(sim, inductive_on);
    gr_sim_set_em_magnetic_enabled(sim, magnetic_on);
    gr_sim_set_em_inductive_sign(sim, ind_sign);
    gr_sim_set_j_deposit_shift(sim, j_shift);

    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const float KE0 = kinetic_energy(p0, c_eff);
    const float vy0 = vy_of(p0, c_eff);

    run_result_t res = {0};
    res.KE_initial = KE0;
    res.vy_initial = vy0;
    res.y_initial  = p0->y;
    res.KE_max_excursion = 0.0f;

    /* Let the field "spin up" around the particle for a few cell-traversal
     * times before we start tracking, so we are not measuring start-up
     * transients.  Use ~8 cell traversal times, capped at 2000 steps
     * (2000 * dt * c = 1414 cells of light-travel -- plenty to settle
     * the wake even at very low v). */
    int n_warmup = (int)(8.0f / (v_drift * dt));
    if (n_warmup > 2000) n_warmup = 2000;
    if (n_warmup < 50)   n_warmup = 50;
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    const gr_particle_t* p_after_warmup = gr_sim_get_particle(sim, 0);
    const float KE_post_warmup = kinetic_energy(p_after_warmup, c_eff);
    const float vy_post_warmup = vy_of(p_after_warmup, c_eff);

    /* Per-step force-measurement window.  Sample dvy over the first
     * force_window_steps of the main loop, after warmup.  Short
     * enough that the particle has not yet been captured to the
     * spurious-force fixed point (at higher mass) and gives a clean
     * F/m = dvy/dt reading at the initial v.
     *
     * Choice: 200 steps.  At dt = 0.7 this is 140 units of light-travel
     * time, which is many wake-formation times but still tiny compared
     * to the run length. */
    const int force_window_steps = 200;
    float vy_at_window_end = vy_post_warmup;
    int window_completed = 0;

    /* Now run the main trajectory.  Sample KE at intervals and track the
     * max excursion from the post-warmup KE. */
    const int sample_every = 50;
    int actual_steps = 0;
    /* Bail-out bounds: stay inside the well-formed PML-buffered region
     * on BOTH sides.  At low v the spurious force can drive vy negative
     * and push the particle into the bottom absorber. */
    const float y_min_safe = 16.0f + 16.0f;  /* absorber + small buffer */
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        actual_steps++;
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        if (s == force_window_steps - 1) {
            vy_at_window_end = vy_of(p, c_eff);
            window_completed = 1;
        }
        if (p->y > y_max_safe) break;
        if (p->y < y_min_safe) break;
        if (!isfinite(p->x) || !isfinite(p->y) || !isfinite(p->px) || !isfinite(p->py)) {
            res.bailed_out = 1;
            break;
        }
        if ((s % sample_every) == 0) {
            const float KE = kinetic_energy(p, c_eff);
            const float ex = fabsf(KE - KE_post_warmup);
            if (ex > res.KE_max_excursion) res.KE_max_excursion = ex;
        }
    }
    /* If we exited the loop before completing the force window (e.g.
     * tiny n_steps or early bail-out), use whatever we have, scaled. */
    if (!window_completed && actual_steps > 0) {
        const gr_particle_t* p_now = gr_sim_get_particle(sim, 0);
        vy_at_window_end = vy_of(p_now, c_eff);
        res.accel_per_step = (vy_at_window_end - vy_post_warmup) / (float) actual_steps;
        res.force_window_steps = actual_steps;
    } else {
        res.accel_per_step    = (vy_at_window_end - vy_post_warmup) / (float) force_window_steps;
        res.force_window_steps = force_window_steps;
    }
    res.dt_used   = dt;
    res.mass_used = mass;

    const gr_particle_t* pf = gr_sim_get_particle(sim, 0);
    res.KE_final  = kinetic_energy(pf, c_eff);
    res.vy_final  = vy_of(pf, c_eff);
    res.y_final   = pf->y;
    res.n_steps_run = actual_steps;

    (void) label;
    gr_sim_destroy(sim);
    return res;
}

/* Wrapper layers.  Defaults moved BARE: rho_smooth = 0, j_smooth = 0,
 * mass = 0.01, dx = 1.0, cfl = 1/sqrt(2), c_eff = 1.0, j_shift = 0.0.
 * The early "core configurations" block at the top of main() now runs
 * at rho_smooth = 0 (the bare-deposit honest baseline), which makes
 * the section's numbers DIFFERENT from prior runs that used
 * rho_smooth = 4 implicitly.  See gr_sandbox_v36 sec. "Stage 37" for
 * what the rho_smooth scan revealed: the +1587% bare phi-only drift
 * was being hidden by 4-pass smoothing in earlier work. */
static run_result_t run_config_full_cfl_dx(int electrostatic_on, int inductive_on, int magnetic_on,
                                    int lorentz_on, float ind_sign, float v_drift,
                                    float c_eff, float j_shift,
                                    int rho_smooth, int j_smooth,
                                    float cfl, float dx,
                                    const char* label) {
    return run_config_full_mass(electrostatic_on, inductive_on, magnetic_on, lorentz_on,
                                ind_sign, v_drift, c_eff, j_shift,
                                rho_smooth, j_smooth, cfl, dx, 0.01f, label);
}
static run_result_t run_config_full_cfl(int electrostatic_on, int inductive_on, int magnetic_on,
                                        int lorentz_on, float ind_sign, float v_drift,
                                        float c_eff, float j_shift,
                                        int rho_smooth, int j_smooth,
                                        float cfl,
                                        const char* label) {
    return run_config_full_cfl_dx(electrostatic_on, inductive_on, magnetic_on, lorentz_on,
                                  ind_sign, v_drift, c_eff, j_shift,
                                  rho_smooth, j_smooth, cfl, 1.0f, label);
}
static run_result_t run_config_full(int electrostatic_on, int inductive_on, int magnetic_on,
                                    int lorentz_on, float ind_sign, float v_drift,
                                    float c_eff, float j_shift,
                                    int rho_smooth, int j_smooth,
                                    const char* label) {
    return run_config_full_cfl(electrostatic_on, inductive_on, magnetic_on, lorentz_on,
                               ind_sign, v_drift, c_eff, j_shift,
                               rho_smooth, j_smooth, 1.0f / sqrtf(2.0f), label);
}
static run_result_t run_config_c(int electrostatic_on, int inductive_on, int magnetic_on,
                                 int lorentz_on, float ind_sign, float v_drift,
                                 float c_eff, const char* label) {
    return run_config_full(electrostatic_on, inductive_on, magnetic_on, lorentz_on,
                           ind_sign, v_drift, c_eff, 0.0f, /*rho_smooth*/0, /*j_smooth*/0, label);
}
static run_result_t run_config(int electrostatic_on, int inductive_on, int magnetic_on,
                               int lorentz_on, float ind_sign, float v_drift,
                               const char* label) {
    return run_config_full(electrostatic_on, inductive_on, magnetic_on, lorentz_on,
                           ind_sign, v_drift, 1.0f, 0.0f, /*rho_smooth*/0, /*j_smooth*/0, label);
}

static void print_result(const char* label, const run_result_t* r) {
    const double KE0 = (double) r->KE_initial;
    const double KEf = (double) r->KE_final;
    const double drift_pct = (KE0 > 0.0) ? 100.0 * (KEf - KE0) / KE0 : 0.0;
    const double v_drift_pct = (r->vy_initial != 0.0f)
        ? 100.0 * ((double) r->vy_final - (double) r->vy_initial) / (double) r->vy_initial
        : 0.0;
    printf("%-50s  KE: %.6e -> %.6e  (%+7.3f%%)   vy: %.5f -> %.5f  (%+7.3f%%)   "
           "y: %.1f -> %.1f   steps=%d%s\n",
           label,
           KE0, KEf, drift_pct,
           (double) r->vy_initial, (double) r->vy_final, v_drift_pct,
           (double) r->y_initial, (double) r->y_final,
           r->n_steps_run,
           r->bailed_out ? "  [BAILED]" : "");
}

int main(void) {
    printf("=== stage37_em_constant_velocity ===\n");
    printf("128 x 1024 asymmetric domain.  Single +Q charge launched at low y\n");
    printf("with constant +y velocity.  No background, no second particle.\n");
    printf("A constant-velocity charge in vacuum MUST NOT radiate.  Any KE\n");
    printf("drift here is a discrete PIC artifact -- and any drift that\n");
    printf("DISAPPEARS when the inductive piece is gated off is direct\n");
    printf("evidence that the -q d_t A channel is the source.\n\n");

    const float v_main = 0.30f;

    /* Config (a): full EM Lorentz, sign +1 (production default). */
    run_result_t r_full = run_config(/*es*/1, /*ind*/1, /*mag*/1, /*lor*/1, +1.0f, v_main,
                                     "full");
    /* Config (b): no inductive -- phi + v x B. */
    run_result_t r_noind = run_config(/*es*/1, /*ind*/0, /*mag*/1, /*lor*/1, +1.0f, v_main,
                                      "no-inductive");
    /* Config (c): no magnetic -- phi + inductive (no v x B). */
    run_result_t r_nomag = run_config(/*es*/1, /*ind*/1, /*mag*/0, /*lor*/1, +1.0f, v_main,
                                      "no-magnetic");
    /* Config (d): phi only -- inductive AND magnetic off. */
    run_result_t r_phionly = run_config(/*es*/1, /*ind*/0, /*mag*/0, /*lor*/1, +1.0f, v_main,
                                        "phi-only");
    /* Config (e): pure inertial control -- no EM force at all. */
    run_result_t r_noforce = run_config(/*es*/0, /*ind*/0, /*mag*/0, /*lor*/0, +1.0f, v_main,
                                        "no-force");

    printf("--- core configurations (v_drift = %.2f c) ---\n", (double) v_main);
    print_result("(a) full        (-grad phi - d_t A + v x B)     ", &r_full);
    print_result("(b) no inductive (-grad phi + v x B)            ", &r_noind);
    print_result("(c) no magnetic (-grad phi - d_t A)             ", &r_nomag);
    print_result("(d) phi only    (-grad phi)                     ", &r_phionly);
    print_result("(e) no force    (pure inertial control)         ", &r_noforce);

    printf("\n");
    {
        const double KE0    = (double) r_full.KE_initial;
        const double a_pct  = 100.0 * ((double) r_full.KE_final     - KE0) / KE0;
        const double b_pct  = 100.0 * ((double) r_noind.KE_final    - KE0) / KE0;
        const double c_pct  = 100.0 * ((double) r_nomag.KE_final    - KE0) / KE0;
        const double d_pct  = 100.0 * ((double) r_phionly.KE_final  - KE0) / KE0;
        const double e_pct  = 100.0 * ((double) r_noforce.KE_final  - KE0) / KE0;

        printf("KE drift summary (constant-velocity charge, target: 0%%):\n");
        printf("  (a) full:                  %+12.4f%%\n", a_pct);
        printf("  (b) no inductive (phi+vxB):%+12.4f%%\n", b_pct);
        printf("  (c) no magnetic (phi+ind): %+12.4f%%\n", c_pct);
        printf("  (d) phi only:              %+12.4f%%   <-- isolated phi self-force\n", d_pct);
        printf("  (e) no force:              %+12.4f%%   <-- integrator-noise floor\n", e_pct);
        printf("\n");
        printf("Channel contributions (at v=%.2f c):\n", (double) v_main);
        printf("  phi alone (d):              %+10.4f%%\n", d_pct);
        printf("  v x B added (b - d):        %+10.4f%%\n", b_pct - d_pct);
        printf("  inductive added (c - d):    %+10.4f%%\n", c_pct - d_pct);
        printf("  all three together (a):     %+10.4f%%\n", a_pct);
    }

    /* Velocity scan: side-by-side full force vs no-inductive baseline.
     * - "full" includes -grad phi - d_t A + v x B.
     * - "noind" disables -d_t A; what remains (-grad phi + v x B) is the
     *   "uncancelled phi self-force" baseline.
     * The difference (full - noind) isolates the inductive piece's
     * contribution to the spurious self-force on uniform motion. */
    /* Velocity scan disabled to keep runtime short -- see prior stage37 runs
     * for the v sweep at c_eff=1.  This test now focuses on the c_eff scan. */

    /* c_eff scan disabled to keep runtime short -- see prior stage37 runs
     * for the c_eff sweep result.  Focus this run on the smoothing test. */
    if (0) {
    printf("\n--- c_eff scan (phi-only, v_drift = 0.30) ---\n");
    printf("As v/c shrinks the field response approaches the static Poisson\n");
    printf("Green's function; any residual drift is sub-cell deposit+gradient.\n\n");
    const float c_vals[] = { 1.0f, 2.0f, 4.0f };
    const int   nc       = (int)(sizeof(c_vals) / sizeof(c_vals[0]));
    printf("%-8s %-8s %-15s\n", "c_eff", "v/c", "KE drift (phi only)");
    printf("--------------------------------------------\n");
    for (int i = 0; i < nc; i++) {
        run_result_t r = run_config_c(/*es*/1, /*ind*/0, /*mag*/0, /*lor*/1,
                                       +1.0f, 0.30f, c_vals[i], "c-scan");
        const double KE0 = (double) r.KE_initial;
        const double pct = (KE0 > 0.0)
                           ? 100.0 * ((double) r.KE_final - KE0) / KE0
                           : 0.0;
        printf("%-8.2f %-8.4f %+11.4f%%%s\n",
               (double) c_vals[i], 0.30 / (double) c_vals[i], pct,
               r.bailed_out ? "  [BAIL]" : "");
    }
    printf("\n");
    printf("Interpretation:\n");
    printf("  If drift -> 0 as c_eff grows: artifact is in the Lienard-Wiechert\n");
    printf("                                wake / time-propagation channel.\n");
    printf("  If drift stays finite:        artifact is in the sub-cell deposit+\n");
    printf("                                gradient self-force on a moving source.\n");
    printf("                                The phi gradient at a sub-cell-displaced\n");
    printf("                                particle position simply isn't zero on\n");
    printf("                                this grid for a moving deposit.\n");
    }  /* end disabled c_eff scan */

    /* J-shift scan also disabled for this run. */
    if (0) {
    /* J-deposit-shift scan: ALL THREE force pieces enabled (electrostatic,
     * magnetic, inductive), with the Esirkepov J trajectory shifted by
     * shift * v * dt relative to the rho deposit at x^n.  Hypothesis:
     * the rho/J spatial mismatch (J midpoint at x^n - v dt/2, rho at x^n)
     * is what prevents -grad phi from cancelling -d_t A on a uniformly-
     * moving charge.  Predicted optimum: shift = +0.5 (J midpoint at x^n,
     * colocated with rho). */
    printf("\n--- J-shift scan (full force, v_drift = %.2f c) ---\n", (double) v_main);
    printf("All three EM force pieces enabled.  rho stays at x^n; J trajectory\n");
    printf("shifted by shift * v * dt.  Continuity is broken when shift != 0.\n");
    printf("Predicted optimum: shift = +0.5 (J midpoint colocated with rho).\n\n");
    const float shifts[] = { -1.0f, -0.5f, -0.25f, 0.0f, +0.25f, +0.5f, +0.75f, +1.0f };
    const int   ns       = (int)(sizeof(shifts) / sizeof(shifts[0]));
    printf("%-9s %-15s %-10s\n", "j_shift", "KE drift (full)", "J@x_eff");
    printf("------------------------------------------------\n");
    for (int i = 0; i < ns; i++) {
        run_result_t r = run_config_full(/*es*/1, /*ind*/1, /*mag*/1, /*lor*/1,
                                          +1.0f, v_main, 1.0f, shifts[i],
                                          /*rho_smooth*/4, /*j_smooth*/0, "j-scan");
        const double KE0 = (double) r.KE_initial;
        const double pct = (KE0 > 0.0)
                           ? 100.0 * ((double) r.KE_final - KE0) / KE0
                           : 0.0;
        const double xeff_offset = ((double) shifts[i] - 0.5) * (double) v_main;
        /* xeff_offset is the J-midpoint position relative to x^n in units of dt.
         * 0 means J center is at x^n.  Negative = J behind rho.  Positive = J ahead. */
        printf("%-9.3f %+11.4f%%     %+7.4f%s\n",
               (double) shifts[i], pct, xeff_offset,
               r.bailed_out ? "  [BAIL]" : "");
    }
    printf("\n");
    printf("  J@x_eff = (shift - 0.5) * v_drift: position of J's effective\n");
    printf("           midpoint relative to rho^n, in length units of v*dt.\n");
    printf("           0 = J colocated with rho.\n");
    }  /* end disabled j-shift scan */

    /* Smoothing scan disabled -- hypothesis already confirmed in prior run.
     * At passes=16, TSC-full residual drops to -5.7% (comparable to CIC's
     * +7%), and phi-only collapses from +1587% (passes=0) to +7% (passes=16).
     * See grlite-stage37-uniform-motion-heating.md for the table. */
    if (0) {
    printf("\n--- smoothing scan: does damping kill the artifact? ---\n");
    printf("Tests the hypothesis that TSC's larger -d_t A artifact comes from\n");
    printf("better-preserved (less-damped) wake content vs CIC.\n\n");
    const int passes[] = { 0, 2, 8 };
    const int np      = (int)(sizeof(passes) / sizeof(passes[0]));
    printf("%-10s %-15s %-15s %-15s\n",
           "passes", "phi-only", "no-inductive", "full");
    printf("(rho_sm = j_sm = passes)\n");
    printf("---------------------------------------------------------------\n");
    for (int i = 0; i < np; i++) {
        run_result_t rd = run_config_full(1, 0, 0, 1, +1.0f, v_main, 1.0f, 0.0f,
                                          passes[i], passes[i], "sm-phi");
        run_result_t rb = run_config_full(1, 0, 1, 1, +1.0f, v_main, 1.0f, 0.0f,
                                          passes[i], passes[i], "sm-noind");
        run_result_t ra = run_config_full(1, 1, 1, 1, +1.0f, v_main, 1.0f, 0.0f,
                                          passes[i], passes[i], "sm-full");
        const double KE0 = (double) ra.KE_initial;
        const double pd  = 100.0 * ((double) rd.KE_final - KE0) / KE0;
        const double pb  = 100.0 * ((double) rb.KE_final - KE0) / KE0;
        const double pa  = 100.0 * ((double) ra.KE_final - KE0) / KE0;
        printf("%-10d %+12.4f%%   %+12.4f%%   %+12.4f%%\n",
               passes[i], pd, pb, pa);
    }
    printf("\n");
    printf("Hypothesis: TSC-larger-artifact comes from dispersion-preserved wake.\n");
    printf("  If TSC artifact -> CIC-like as passes grows: confirms hypothesis.\n");
    printf("  If TSC artifact stays large:                 hypothesis refuted.\n");
    }  /* end disabled smoothing scan */

    /* CFL scan disabled -- already established (this session) that drift is
     * essentially flat as dt shrinks 5.66x.  See memory file for table. */
    if (0) {
    /* CFL / dt scan: hold v_drift, c_eff, dx fixed; vary CFL (and therefore
     * dt = CFL * dx / c_eff).  Same physical trajectory, more time-steps to
     * cover it at lower CFL.
     *
     * Predictions if artifact is dominantly...
     *   ...time-truncation (leapfrog is O(dt^2)):
     *        |drift| should fall as ~CFL^2.  Halving CFL -> ~1/4 the drift.
     *   ...spatial discretization (sub-cell self-force on moving deposit):
     *        |drift| should be roughly insensitive to CFL.
     *   ...wave-equation numerical dispersion:
     *        |drift| may GROW as CFL drops, because the 2D Yee wave-eqn
     *        phase error worsens off the magic CFL=1/sqrt(2).
     *
     * The previous tests fix rho_smooth=4, j_smooth=0; we hold those here so
     * results are directly comparable to the (a)/(b)/(c)/(d) numbers above. */
    printf("\n--- CFL / dt scan (full + phi-only + no-inductive, v_drift=0.30 c) ---\n");
    printf("dx and c_eff held fixed; dt = CFL * dx / c_eff.  Same trajectory,\n");
    printf("more steps per traversal at lower CFL.  rho_smooth=4, j_smooth=0\n");
    printf("(production baseline).  If artifact is leapfrog time-truncation,\n");
    printf("|drift| should fall as ~CFL^2.\n\n");
    const float cfl_vals[] = { 0.7071f, 0.5f, 0.25f, 0.125f };
    const int   ncfl       = (int)(sizeof(cfl_vals) / sizeof(cfl_vals[0]));
    printf("%-8s %-10s %-15s %-15s %-15s\n",
           "CFL", "dt", "phi-only", "no-inductive", "full");
    printf("--------------------------------------------------------------------\n");
    for (int i = 0; i < ncfl; i++) {
        const float c_eff_here = 1.0f;
        const float dx_here    = 1.0f;
        const float dt_here    = cfl_vals[i] * dx_here / c_eff_here;
        run_result_t rd = run_config_full_cfl(1, 0, 0, 1, +1.0f, v_main, c_eff_here,
                                              0.0f, 4, 0, cfl_vals[i], "cfl-phi");
        run_result_t rb = run_config_full_cfl(1, 0, 1, 1, +1.0f, v_main, c_eff_here,
                                              0.0f, 4, 0, cfl_vals[i], "cfl-noind");
        run_result_t ra = run_config_full_cfl(1, 1, 1, 1, +1.0f, v_main, c_eff_here,
                                              0.0f, 4, 0, cfl_vals[i], "cfl-full");
        const double KE0 = (double) ra.KE_initial;
        const double pd  = 100.0 * ((double) rd.KE_final - KE0) / KE0;
        const double pb  = 100.0 * ((double) rb.KE_final - KE0) / KE0;
        const double pa  = 100.0 * ((double) ra.KE_final - KE0) / KE0;
        printf("%-8.4f %-10.4f %+12.4f%%   %+12.4f%%   %+12.4f%%\n",
               (double) cfl_vals[i], (double) dt_here, pd, pb, pa);
    }
    printf("\n");
    printf("Interpretation:\n");
    printf("  phi-only scaling ~ CFL^p:  p ~ 2 means leapfrog time-truncation.\n");
    printf("                              p ~ 0 means spatial-discretization-limited.\n");
    printf("                              p < 0 means dispersion-dominated (gets worse).\n");
    printf("  full scaling ~ CFL^p:      same interpretation for the cancelled total.\n");
    }  /* end disabled CFL scan */

    /* dx scan disabled -- data captured in memory file.  Halving dx grows
     * phi-only drift by 1.53x (UV-divergent). */
    if (0) {
    /* dx convergence scan: hold the PHYSICAL setup fixed.  Cells (W, H, n_damp)
     * scale as 1/dx; dt = CFL * dx / c also scales with dx, so n_steps scales
     * as 1/dx.  Halving dx costs ~8x per run (4x cells * 2x steps).  Keep to
     * 2-3 dx values, phi-only and full only.
     *
     * Predictions if artifact is dominantly...
     *   ...spatial-discretization that CONVERGES (e.g. O(dx) or O(dx^2)):
     *       |drift| should fall as some power of dx as the scheme converges.
     *   ...structural cancellation failure that survives the continuum limit:
     *       |drift| roughly insensitive to dx.  The continuous PIC also
     *       fails to cancel; refinement does not help.
     *   ...spatial residual that GROWS in the continuum limit (e.g. an
     *       infrared-divergent self-force from a regularized delta):
     *       |drift| grows as dx -> 0.
     *
     * This is the test that should distinguish "fixable with a finer grid"
     * from "structurally broken cancellation". */
    printf("\n--- dx convergence scan (v_drift=0.30 c) ---\n");
    printf("Physical setup held fixed (128 x 1024 phys; absorbers and buffers\n");
    printf("scale with cells).  dt scales with dx; n_steps scales as 1/dx.\n");
    printf("Cost per dx halving is ~8x.\n\n");
    const float dx_vals[] = { 1.0f, 0.5f };
    const int   ndx       = (int)(sizeof(dx_vals) / sizeof(dx_vals[0]));
    printf("%-6s %-10s %-15s %-15s %-15s\n",
           "dx", "n_steps", "phi-only", "no-inductive", "full");
    printf("--------------------------------------------------------------------\n");
    double prev_phi  = 0.0;
    double prev_full = 0.0;
    for (int i = 0; i < ndx; i++) {
        run_result_t rd = run_config_full_cfl_dx(1, 0, 0, 1, +1.0f, v_main, 1.0f, 0.0f,
                                                 4, 0, 1.0f/sqrtf(2.0f), dx_vals[i], "dx-phi");
        run_result_t rb = run_config_full_cfl_dx(1, 0, 1, 1, +1.0f, v_main, 1.0f, 0.0f,
                                                 4, 0, 1.0f/sqrtf(2.0f), dx_vals[i], "dx-noind");
        run_result_t ra = run_config_full_cfl_dx(1, 1, 1, 1, +1.0f, v_main, 1.0f, 0.0f,
                                                 4, 0, 1.0f/sqrtf(2.0f), dx_vals[i], "dx-full");
        const double KE0 = (double) ra.KE_initial;
        const double pd  = 100.0 * ((double) rd.KE_final - KE0) / KE0;
        const double pb  = 100.0 * ((double) rb.KE_final - KE0) / KE0;
        const double pa  = 100.0 * ((double) ra.KE_final - KE0) / KE0;
        printf("%-6.3f %-10d %+12.4f%%   %+12.4f%%   %+12.4f%%",
               (double) dx_vals[i], ra.n_steps_run, pd, pb, pa);
        if (i > 0) {
            const double r_phi  = (prev_phi  != 0.0) ? pd / prev_phi  : 0.0;
            const double r_full = (prev_full != 0.0) ? pa / prev_full : 0.0;
            printf("   ratio phi=%.3f full=%.3f", r_phi, r_full);
        }
        printf("\n");
        prev_phi  = pd;
        prev_full = pa;
    }
    printf("\n");
    printf("Interpretation of dx halving ratio:\n");
    printf("  ratio ~ 0.25 => O(dx^2) convergence (refinement helps fast)\n");
    printf("  ratio ~ 0.50 => O(dx)   convergence (refinement helps linearly)\n");
    printf("  ratio ~ 1.00 => structural; refinement does NOT help.\n");
    printf("  ratio > 1    => artifact gets WORSE with refinement (UV-divergent).\n");
    }  /* end disabled dx scan */

    /* v_drift scan: does the artifact persist at non-relativistic speeds?
     *
     * Two challenges this scan addresses:
     *  (a) Low v has tiny initial KE (0.5 m v^2), so % KE drift explodes
     *      for any small absolute force.  Report dvy/dt in absolute
     *      units as the primary force-density measure.
     *  (b) At low v with mass = 0.01, the spurious force quickly
     *      captures the particle to the fixed point at vy ~ -0.28c,
     *      so any measurement integrated over the run is dominated by
     *      the captured dynamics, not the force at v_init.  Solution:
     *      use a HEAVIER particle (mass = 1.0, 100x current default)
     *      so the velocity barely changes during the 200-step force
     *      window.  Absolute force is unchanged (same q, same field
     *      operator); only velocity response per step shrinks.
     *
     * Hypothesis being tested: the spurious self-force is structural
     * to the discrete -grad phi / -d_t A pairing on a moving deposit.
     * It should be present at ALL non-zero v including v/c = 0.001. */
    printf("\n--- v_drift scan (BARE deposit, mass=1.0 for clean window measurement) ---\n");
    printf("rho_smooth=0, j_smooth=0, q=0.01, m=1.0.  Particle stays near v_init\n");
    printf("over the 200-step force-measurement window.  accel column reports\n");
    printf("dvy/dt (per unit time) averaged over that window AT THE INITIAL v.\n");
    printf("In a working scheme the accel column would be at the no-force\n");
    printf("integrator-noise floor at every v.\n\n");
    const float v_vals[] = { 0.001f, 0.01f, 0.05f, 0.1f, 0.2f, 0.3f, 0.5f, 0.7f };
    const int   nv       = (int)(sizeof(v_vals) / sizeof(v_vals[0]));
    /* Also a "no force" control at one v for the floor. */
    run_result_t noforce_floor = run_config_full_mass(0, 0, 0, 0, +1.0f, 0.1f, 1.0f, 0.0f,
                                                     0, 0, 1.0f/sqrtf(2.0f), 1.0f, 1.0f, "v-noforce");
    const double a_noforce = (double) noforce_floor.accel_per_step / (double) noforce_floor.dt_used;
    printf("Integrator-noise floor at v=0.1, no-force control:\n");
    printf("  accel = %+10.3e  (this is the threshold a fix must beat)\n\n",
           a_noforce);
    printf("%-7s %-13s %-13s %-13s %-11s %-11s %-11s\n",
           "v/c", "accel(phi)", "accel(noind)", "accel(full)", "dvy_phi", "dvy_full", "KE%(phi)");
    printf("-------------------------------------------------------------------------------------\n");
    for (int i = 0; i < nv; i++) {
        run_result_t rd = run_config_full_mass(1, 0, 0, 1, +1.0f, v_vals[i], 1.0f, 0.0f,
                                              0, 0, 1.0f/sqrtf(2.0f), 1.0f, 1.0f, "v-phi");
        run_result_t rb = run_config_full_mass(1, 0, 1, 1, +1.0f, v_vals[i], 1.0f, 0.0f,
                                              0, 0, 1.0f/sqrtf(2.0f), 1.0f, 1.0f, "v-noind");
        run_result_t ra = run_config_full_mass(1, 1, 1, 1, +1.0f, v_vals[i], 1.0f, 0.0f,
                                              0, 0, 1.0f/sqrtf(2.0f), 1.0f, 1.0f, "v-full");
        const double KE0d = (double) rd.KE_initial;
        const double pd = (KE0d > 0.0) ? 100.0 * ((double) rd.KE_final - KE0d) / KE0d : 0.0;
        const double dvy_phi  = (double) (rd.vy_final - rd.vy_initial);
        const double dvy_full = (double) (ra.vy_final - ra.vy_initial);
        const double accel_d = (double) rd.accel_per_step / (double) rd.dt_used;
        const double accel_b = (double) rb.accel_per_step / (double) rb.dt_used;
        const double accel_a = (double) ra.accel_per_step / (double) ra.dt_used;
        printf("%-7.4f %+12.3e %+12.3e %+12.3e %+10.3e %+10.3e %+11.4f%%\n",
               (double) v_vals[i], accel_d, accel_b, accel_a,
               dvy_phi, dvy_full, pd);
    }
    printf("\n");
    printf("Interpretation:\n");
    printf("  accel column: average dvy/dt over a 200-step window after warmup,\n");
    printf("    at the INITIAL v (particle didn't drift far in those steps because\n");
    printf("    mass=1.0).  This is the spurious-acceleration-at-v measurement.\n");
    printf("  A working scheme should have all accel entries at the no-force floor\n");
    printf("    (~%+9.2e) at every v.\n", a_noforce);
    printf("  Empirical v-scaling of accel: if accel ~ v^k, then on a log-log\n");
    printf("    plot the slope is k.  At low v, k=0 means constant force\n");
    printf("    (independent of v); k=1 friction-like; k=2 radiation-like.\n");

    printf("\nBASELINE COMPLETE (current scheme; numbers are the failure pattern\n");
    printf("GEMPIC must beat per gr_sandbox_v36 sec gempic_benchmark).\n");
    return 0;
}
