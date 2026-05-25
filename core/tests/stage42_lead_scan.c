/* Stage 42 -- diagnostic experiment: shift the EM field-gather position
 * ahead of the particle by lead * v * dt and scan lead.
 *
 * Hypothesis (Dan, 2026-05-24): the discrete wake's effective center
 * of symmetry is offset from the particle's actual position by some
 * fraction of a step (the wake "trails" the particle because the
 * wave equation propagates at finite c on the grid).  If we gather
 * the field AHEAD of the particle in the direction of motion, we
 * might land on a point where the discrete grad phi vanishes.
 *
 * Method: same Stage 37 setup; scan lead in {-1, -0.5, -0.25, 0,
 * +0.25, +0.5, +1.0}; for each lead value report per-step spurious
 * acceleration at a representative v.  Run for both pushers (Boris
 * default and GEMPIC) so we can see whether the same lead value is
 * optimal in both.
 *
 * Outcomes:
 *   - If an optimal lead exists where accel ~ noise floor: deposit-
 *     side mitigation found.
 *   - If accel is monotone in lead at all v: indicates direction but
 *     no clean cancellation point.
 *   - If optimum varies wildly with v: the offset isn't a constant
 *     fraction of v*dt and we need a more sophisticated correction. */

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
    float dt_used;
    int   bailed;
} lead_result_t;

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m * m * c_eff * c_eff));
    return p->py / (gamma * m);
}

static lead_result_t run_one(float v_drift, float lead, gr_pusher_kind_t pusher) {
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

    lead_result_t r = {0};
    r.dt_used = dt;

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
    gr_sim_set_phi_gather_lead(sim, lead);

    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    const int n_warmup = 2000;
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    const gr_particle_t* p_post = gr_sim_get_particle(sim, 0);
    const float vy_post = vy_of(p_post, c_eff);

    const int force_window = 200;
    for (int s = 0; s < force_window; s++) gr_sim_step(sim);

    const gr_particle_t* p_after = gr_sim_get_particle(sim, 0);
    const float vy_end = vy_of(p_after, c_eff);
    r.accel_per_step = (vy_end - vy_post) / (float) force_window;

    gr_sim_destroy(sim);
    return r;
}

static void scan_at(float v, gr_pusher_kind_t pusher, const char* pusher_name) {
    const float leads[] = { -1.0f, -0.5f, -0.25f, 0.0f, +0.25f, +0.5f, +1.0f };
    const int   nl      = (int)(sizeof(leads) / sizeof(leads[0]));
    printf("--- v_drift = %.4f, pusher = %s ---\n", (double) v, pusher_name);
    printf("%-7s %-13s\n", "lead", "accel/dt");
    printf("---------------------------\n");
    double best_abs = 1.0e30;
    float best_lead = 0.0f;
    for (int i = 0; i < nl; i++) {
        lead_result_t r = run_one(v, leads[i], pusher);
        const double a = (double) r.accel_per_step / (double) r.dt_used;
        printf("%+5.2f  %+12.3e%s\n",
               (double) leads[i], a, r.bailed ? "  [BAIL]" : "");
        if (fabs(a) < best_abs) { best_abs = fabs(a); best_lead = leads[i]; }
    }
    printf("min |accel| at lead = %+.3f  (|accel| = %.3e)\n\n",
           (double) best_lead, best_abs);
}

int main(void) {
    printf("=== stage42_lead_scan ===\n");
    printf("Scan EM field-gather position ahead-of-particle by lead*v*dt.\n");
    printf("If discrete wake center-of-symmetry is offset from x_p, some\n");
    printf("non-zero lead should give |accel| ~ noise floor.\n\n");

    /* Scan a few representative v values for both pushers. */
    const float v_vals[] = { 0.001f, 0.1f, 0.3f };
    const int   nv       = (int)(sizeof(v_vals) / sizeof(v_vals[0]));
    for (int i = 0; i < nv; i++) {
        scan_at(v_vals[i], GR_PUSHER_BORIS,  "Boris");
    }
    for (int i = 0; i < nv; i++) {
        scan_at(v_vals[i], GR_PUSHER_GEMPIC, "GEMPIC");
    }

    return 0;
}
