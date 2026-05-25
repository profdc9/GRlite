/* Stage 41 -- empirical check: does GEMPIC pusher actually reduce the
 * Stage 37 uniform-motion artifact?
 *
 * Per gr_sandbox_v37 sec gempic_derivation, the Fourier analysis suggests
 * GEMPIC's residual on a uniformly-moving charge is approximately
 *     F^{GEMPIC} ~ -q/gamma^2 * (Stage 37 phi-only force)
 * which at non-relativistic v is essentially the same magnitude as
 * the current scheme's phi-only force.  This test is the empirical
 * check of that prediction.
 *
 * Setup: identical to Stage 37's v-scan section.  Same domain, same
 * bare deposit, mass=1.0 for clean window measurement.  For each v
 * value, run the test twice: once with the default Boris pusher,
 * once with GR_PUSHER_GEMPIC.  Report per-step acceleration
 * (dvy/dt) for both, plus the ratio.
 *
 * Acceptance criteria for GEMPIC to be considered a fix:
 *   |accel(GEMPIC)| / |accel(Boris)| < 0.1 at every v.
 * If the Fourier analysis is right, the ratio will be ~1 at low v and
 * ~(1 - v^2/c^2) at higher v -- in which case GEMPIC is NOT a fix and
 * we should not invest in the implementation polish or feature flag
 * promotion. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float accel_per_step;
    float vy_init;
    float vy_after_window;
    float dt_used;
    int   bailed;
} v_result_t;

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static v_result_t run_one(float v_drift, gr_pusher_kind_t pusher) {
    const int   W       = 128;
    const int   H       = 1024;
    const float dx      = 1.0f;
    const float c_eff   = 1.0f;
    const float cfl     = 1.0f / sqrtf(2.0f);
    const int   n_damp  = 16;
    const float Q       = 0.01f;
    const float mass    = 1.0f;
    const float y_start = 16.0f + 32.0f;
    const float dt      = cfl * dx / c_eff;
    const float cx      = ((float) (W - 1) * 0.5f) * dx;

    v_result_t r = {0};
    r.dt_used = dt;
    r.vy_init = v_drift;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { r.bailed = 1; return r; }
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, +1.0f);
    gr_sim_set_pusher(sim, pusher);

    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    /* Warm-up to let the wake form. */
    const int n_warmup = 2000;
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    const gr_particle_t* p_post = gr_sim_get_particle(sim, 0);
    const float vy_post_warmup = vy_of(p_post, c_eff);

    /* Force-measurement window. */
    const int force_window = 200;
    for (int s = 0; s < force_window; s++) gr_sim_step(sim);

    const gr_particle_t* p_after = gr_sim_get_particle(sim, 0);
    r.vy_after_window = vy_of(p_after, c_eff);
    r.accel_per_step = (r.vy_after_window - vy_post_warmup) / (float) force_window;

    gr_sim_destroy(sim);
    return r;
}

int main(void) {
    printf("=== stage41_gempic_check ===\n");
    printf("Empirical comparison of the Stage 37 uniform-motion artifact\n");
    printf("under the Boris pusher (default) and the GEMPIC pusher (v37).\n");
    printf("Same setup as Stage 37 v-scan, mass=1.0, bare deposit.\n");
    printf("accel = dvy/dt averaged over 200-step post-warmup window.\n\n");

    const float v_vals[] = { 0.001f, 0.01f, 0.05f, 0.1f, 0.2f, 0.3f };
    const int   nv       = (int)(sizeof(v_vals) / sizeof(v_vals[0]));

    printf("%-7s %-14s %-14s %-12s %-12s\n",
           "v/c", "accel(Boris)", "accel(GEMPIC)", "ratio", "GEMPIC/Boris");
    printf("---------------------------------------------------------------\n");

    for (int i = 0; i < nv; i++) {
        v_result_t boris  = run_one(v_vals[i], GR_PUSHER_BORIS);
        v_result_t gempic = run_one(v_vals[i], GR_PUSHER_GEMPIC);
        const double aB = (double) boris.accel_per_step  / (double) boris.dt_used;
        const double aG = (double) gempic.accel_per_step / (double) gempic.dt_used;
        const double ratio = (aB != 0.0) ? aG / aB : 0.0;
        const double mag_ratio = (fabs(aB) > 0.0)
                                 ? fabs(aG) / fabs(aB) : 0.0;
        printf("%-7.4f %+13.3e %+13.3e %+11.3f %11.3f%s%s\n",
               (double) v_vals[i], aB, aG, ratio, mag_ratio,
               boris.bailed  ? " B-BAIL"  : "",
               gempic.bailed ? " G-BAIL" : "");
    }

    printf("\nInterpretation:\n");
    printf("  Fourier-predicted GEMPIC/Boris ratio (magnitude):\n");
    printf("    v=0.001: ~1.000\n");
    printf("    v=0.1:   ~0.990  (1 - v^2/c^2)\n");
    printf("    v=0.3:   ~0.910\n");
    printf("  If measured ratio matches prediction => GEMPIC is NOT a fix.\n");
    printf("  If measured ratio is << 1 at every v => GEMPIC IS a fix and\n");
    printf("    we should pursue feature-flag promotion.\n");
    printf("  Acceptance for fix: ratio < 0.1 at every v.\n");

    return 0;
}
