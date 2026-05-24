/* Stage 38 -- BENCHMARK B-2: Neutral pair co-moving (gr_sandbox_v36 sec gempic_benchmark).
 *
 * Two equal-and-opposite charges (+Q, -Q) at the same velocity (0, +v).
 * They are separated transversely by `sep` cells along the x-axis.
 * Net rho integrated over the grid is zero (modulo finite deposit
 * shape); net J integrated over the grid is zero (J_y from +Q is
 * cancelled point-for-point by J_y from -Q at the joint x-midpoint
 * line; the two deposits are offset transversely so they do not
 * cancel locally, but the LINE INTEGRAL of J across any plane
 * perpendicular to motion is zero).  At long range the joint field is
 * dipole.
 *
 * Physics expectation: each particle's KE conserved up to dipole-wake
 * radiation (very small at modest v).  No spurious self-force.
 *
 * Why this probe (per v36 sec gempic_benchmark): isolates self-force
 * from "interaction with own wake."  The wake the two particles
 * jointly create is structurally smaller than either's individual
 * wake.  The current scheme's UV-divergent self-force will still hit
 * each particle independently and the test will fail.
 *
 * Acceptance for a fix (per v36):
 *   |Delta v_y| < 1e-6 for both particles over 5000 steps at v=0.001
 *   |Delta v_y| < 1e-4 for both particles over 5000 steps at v=0.1
 *   Total momentum P_+ + P_- conserved to machine precision modulo PML.
 *
 * This file is the BASELINE measurement for the current scheme.
 * Numbers recorded here are the failure pattern GEMPIC must beat. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float vy_init;
    float vy_final_p;       /* +Q particle */
    float vy_final_n;       /* -Q particle */
    float vx_final_p;       /* +Q particle lateral velocity (mutual attraction) */
    float vx_final_n;
    float accel_per_step_p; /* dvy / step over 200-step post-warmup window */
    float accel_per_step_n;
    int   force_window_steps;
    float dt_used;
    float mass_used;
    float Q_used;
    int   sep_cells;
    int   n_steps_run;
    int   bailed_out;
    /* For total-momentum-conservation report */
    float total_py_initial;
    float total_py_final;
    float total_px_initial;
    float total_px_final;
} pair_result_t;

static float vx_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->px / (gamma * m);
}
static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static pair_result_t run_neutral_pair(float v_drift, int sep_cells, float Q, float mass,
                                       int n_steps_cap) {
    const int   W       = 128;
    const int   H       = 256;
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 16;
    const float dt      = cfl * dx / c_eff;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;
    const float cy      = ((float) (H - 1) * 0.5f) * dx;
    const float xp      = cx - 0.5f * (float) sep_cells * dx;
    const float xn      = cx + 0.5f * (float) sep_cells * dx;

    pair_result_t res = {0};
    res.dt_used   = dt;
    res.mass_used = mass;
    res.Q_used    = Q;
    res.sep_cells = sep_cells;
    res.vy_init   = v_drift;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { res.bailed_out = 1; return res; }
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);   /* BARE deposit, per v36 sec gempic_benchmark */
    gr_sim_set_j_smooth_passes(sim, 0);

    /* Gravity off; only EM acts. */
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, +1.0f);   /* variationally correct */
    gr_sim_set_j_deposit_shift(sim, 0.0f);

    gr_sim_add_particle(sim, xp, cy, mass, +Q, 0.0f, v_drift);
    gr_sim_add_particle(sim, xn, cy, mass, -Q, 0.0f, v_drift);

    const gr_particle_t* pp0 = gr_sim_get_particle(sim, 0);
    const gr_particle_t* pn0 = gr_sim_get_particle(sim, 1);
    res.total_py_initial = pp0->py + pn0->py;
    res.total_px_initial = pp0->px + pn0->px;

    /* Warm-up so the dipole wake forms.  Cap at 2000 steps; for slow
     * particles this is plenty of light-crossing times. */
    int n_warmup = 2000;
    if (n_warmup < 50) n_warmup = 50;
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    const gr_particle_t* pp_post = gr_sim_get_particle(sim, 0);
    const gr_particle_t* pn_post = gr_sim_get_particle(sim, 1);
    const float vyp_post = vy_of(pp_post, c_eff);
    const float vyn_post = vy_of(pn_post, c_eff);

    const int force_window_steps = 200;
    float vyp_at_end = vyp_post;
    float vyn_at_end = vyn_post;
    int window_completed = 0;

    /* Bail bounds. */
    const float y_min_safe = (float)(n_damp + 4);
    const float y_max_safe = (float)(H - 1 - n_damp - 4);
    const float x_min_safe = (float)(n_damp + 4);
    const float x_max_safe = (float)(W - 1 - n_damp - 4);

    int actual = 0;
    for (int s = 0; s < n_steps_cap; s++) {
        gr_sim_step(sim);
        actual++;
        const gr_particle_t* pp = gr_sim_get_particle(sim, 0);
        const gr_particle_t* pn = gr_sim_get_particle(sim, 1);
        if (s == force_window_steps - 1) {
            vyp_at_end = vy_of(pp, c_eff);
            vyn_at_end = vy_of(pn, c_eff);
            window_completed = 1;
        }
        if (pp->y < y_min_safe || pp->y > y_max_safe ||
            pn->y < y_min_safe || pn->y > y_max_safe ||
            pp->x < x_min_safe || pp->x > x_max_safe ||
            pn->x < x_min_safe || pn->x > x_max_safe) break;
        if (!isfinite(pp->x) || !isfinite(pp->y) || !isfinite(pn->x) || !isfinite(pn->y)) {
            res.bailed_out = 1;
            break;
        }
    }
    if (!window_completed && actual > 0) {
        const gr_particle_t* pp_now = gr_sim_get_particle(sim, 0);
        const gr_particle_t* pn_now = gr_sim_get_particle(sim, 1);
        vyp_at_end = vy_of(pp_now, c_eff);
        vyn_at_end = vy_of(pn_now, c_eff);
        res.accel_per_step_p   = (vyp_at_end - vyp_post) / (float) actual;
        res.accel_per_step_n   = (vyn_at_end - vyn_post) / (float) actual;
        res.force_window_steps = actual;
    } else {
        res.accel_per_step_p   = (vyp_at_end - vyp_post) / (float) force_window_steps;
        res.accel_per_step_n   = (vyn_at_end - vyn_post) / (float) force_window_steps;
        res.force_window_steps = force_window_steps;
    }

    const gr_particle_t* ppf = gr_sim_get_particle(sim, 0);
    const gr_particle_t* pnf = gr_sim_get_particle(sim, 1);
    res.vy_final_p = vy_of(ppf, c_eff);
    res.vy_final_n = vy_of(pnf, c_eff);
    res.vx_final_p = vx_of(ppf, c_eff);
    res.vx_final_n = vx_of(pnf, c_eff);
    res.total_py_final = ppf->py + pnf->py;
    res.total_px_final = ppf->px + pnf->px;
    res.n_steps_run = actual;

    gr_sim_destroy(sim);
    return res;
}

int main(void) {
    printf("=== stage38_em_neutral_pair (BENCHMARK B-2) ===\n");
    printf("Two opposite charges (+Q, -Q) at velocity (0, +v), separated by\n");
    printf("`sep` cells transversely.  Net rho and net J both zero on integration;\n");
    printf("joint far field is dipole.  Physical answer: each particle's vy\n");
    printf("must be conserved to integrator-noise level.  Any drift = bug.\n\n");
    printf("Setup: 128x256 grid, BARE deposit (rho_smooth=0, j_smooth=0),\n");
    printf("mass=1.0, Q=0.01, sep=8 cells, TSC+LB.  Warmup 2000 steps,\n");
    printf("force-window 200 steps.\n\n");

    const float Q       = 0.01f;
    const float mass    = 1.0f;
    const int   sep     = 8;
    const int   n_cap   = 5000;

    /* Acceptance thresholds (from v36 sec gempic_benchmark B-2). */
    const double thresh_v0001 = 1.0e-6;
    const double thresh_v01   = 1.0e-4;

    const float v_vals[] = { 0.001f, 0.01f, 0.05f, 0.1f, 0.2f, 0.3f };
    const int   nv       = (int)(sizeof(v_vals) / sizeof(v_vals[0]));

    printf("%-7s %-12s %-12s %-12s %-12s %-11s %-11s %-12s\n",
           "v/c", "accel(+Q)", "accel(-Q)", "dvy(+Q)", "dvy(-Q)",
           "vx(+Q)", "vx(-Q)", "dPy_total");
    printf("--------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < nv; i++) {
        pair_result_t r = run_neutral_pair(v_vals[i], sep, Q, mass, n_cap);
        const double accel_p = (double) r.accel_per_step_p / (double) r.dt_used;
        const double accel_n = (double) r.accel_per_step_n / (double) r.dt_used;
        const double dvyp = (double) r.vy_final_p - (double) r.vy_init;
        const double dvyn = (double) r.vy_final_n - (double) r.vy_init;
        const double dpy = (double) r.total_py_final - (double) r.total_py_initial;
        printf("%-7.4f %+12.3e %+12.3e %+12.3e %+12.3e %+10.3e %+10.3e %+11.3e%s\n",
               (double) v_vals[i], accel_p, accel_n, dvyp, dvyn,
               (double) r.vx_final_p, (double) r.vx_final_n,
               dpy,
               r.bailed_out ? " [BAIL]" : "");
    }

    printf("\nInterpretation:\n");
    printf("  accel(+Q), accel(-Q): spurious y-acceleration on each particle.\n");
    printf("    Both should be ZERO for a working scheme (the joint field is\n");
    printf("    dipole at long range and each particle's local field is\n");
    printf("    structurally cancelled by its partner's far-field tail).\n");
    printf("  dvy: integrated change.  Acceptance for a fix:\n");
    printf("    |dvy| < %.0e at v=0.001\n", thresh_v0001);
    printf("    |dvy| < %.0e at v=0.1\n", thresh_v01);
    printf("  vx: lateral approach velocity.  Opposite charges attract;\n");
    printf("    these are NEGATIVE for +Q (moving inward toward -Q) and\n");
    printf("    POSITIVE for -Q.  Physical, not a bug.\n");
    printf("  dPy_total: total y-momentum drift.  Should be << dvy of each\n");
    printf("    particle (the wake carries some momentum via PML absorption).\n");
    printf("    Asymmetry between +Q and -Q dvy indicates the operator pair\n");
    printf("    treats the two signs differently.\n");

    printf("\nBASELINE COMPLETE.\n");
    return 0;
}
