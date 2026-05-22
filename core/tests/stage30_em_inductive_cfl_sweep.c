/* Stage 30 — CFL convergence sweep for the EM inductive PIC inspiral.
 *
 * Mirror of Stage 29, but for EM closed-loop PIC (Stage 27 setup) instead
 * of gravity (Stage 28).  Confirms that the EM-side inspiral observed
 * after the half-step A fix is physical (CFL-independent plateau) rather
 * than residual numerical truncation.
 *
 * Setup: charged particle (q=-1e-3, m=1e-3) orbiting a softened POINT_CHARGE
 * background (Q=+1) at r=20.  Field evolution + auto-deposition both on;
 * full EM Lorentz force (electrostatic + magnetic + inductive) all
 * enabled.  Sweep CFL by halving from 1/sqrt(2) down to 1/(8*sqrt(2)).
 *
 * Predictions (per Stage 29 result):
 *   - No-inductive drift should CONVERGE TOWARD ZERO as CFL -> 0
 *     (spatial pieces alone produce no radiation; residual is numerical).
 *   - Inductive-on drift should CONVERGE TO A PLATEAU as CFL -> 0
 *     (the genuine 2D radiation-reaction inspiral).
 *   - Consecutive-CFL contribution ratio should approach 1 (physical),
 *     not 4 (O(dt^2) numerical). */

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
        exit(1); \
    } \
} while (0)

typedef struct {
    float cfl;
    float dt;
    int   n_steps_total;
    float drift_noind;
    float drift_indon;
    float drift_indcontrib;
    int   n_completed_noind;
    int   n_completed_indon;
} sweep_result_t;

static float run(float Q, float q_test, float m_test, float r_orb, float cfl,
                 int n_orbits, int em_inductive_enabled,
                 int* n_completed_out, int* n_steps_total_out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float k_e    = 1.0f;

    const float g_mag  = fabsf(q_test * Q) * k_e * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f) / m_test;
    const float rg     = r_orb * g_mag;
    const float rg2_c2 = rg * rg / (c_eff * c_eff);
    const float u_v2   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ = sqrtf(u_v2);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_circ;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return NAN;
    gr_sim_set_damping(sim, 16);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    /* CIC + LEGACY: EM gradient path doesn't yet have TSC/LB variants.
     * rho smoothing 4 (mirrors gravity production). */
    gr_sim_set_shape_function(sim, GR_SHAPE_CIC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEGACY);
    gr_sim_set_rho_smooth_passes(sim, 4);
    gr_sim_set_em_inductive_enabled(sim, em_inductive_enabled);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test,
                        /*vx=*/0.0f, /*vy=*/v_circ);

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.2f * (float) n_orbits * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    float r_at_end = r_orb;
    int   s_used = 0;
    for (int s = 0; s < n_max && wraps < n_orbits; s++) {
        gr_sim_step(sim);
        s_used = s + 1;
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) { gr_sim_destroy(sim); return NAN; }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            r_at_end = r_now;
        }
        th_prev = th;
    }
    if (n_completed_out) *n_completed_out = wraps;
    if (n_steps_total_out) *n_steps_total_out = s_used;
    gr_sim_destroy(sim);
    if (wraps < n_orbits) return NAN;
    return (r_at_end - r_orb) / r_orb;
}

int main(void) {
    printf("=== stage30_em_inductive_cfl_sweep ===\n");
    printf("Does the EM inductive-piece inspiral rate scale with CFL?\n");
    printf("Setup: EM closed-loop PIC, Q=+1, q=-1e-3, m=1e-3, r=20, 4 orbits.\n\n");

    const float Q       = +1.0f;
    const float q_test  = -1.0e-3f;
    const float m_test  = +1.0e-3f;
    const float r_orb   = 20.0f;
    const int   N       = 4;

    const float cfls[] = {
        1.0f / sqrtf(2.0f),
        1.0f / (2.0f * sqrtf(2.0f)),
        1.0f / (4.0f * sqrtf(2.0f)),
        1.0f / (8.0f * sqrtf(2.0f)),
    };
    const int n_cfl = (int) (sizeof(cfls) / sizeof(cfls[0]));
    sweep_result_t results[8];

    printf("%-6s %-10s %-10s %-12s %-12s %-12s %-12s\n",
           "i", "CFL", "dt", "steps/4-orb", "no induct", "induct +1", "induct contrib");
    printf("---------------------------------------------------------------------------\n");

    for (int i = 0; i < n_cfl; i++) {
        const float cfl = cfls[i];
        const float dt  = cfl;
        int n_noind = 0, n_indon = 0;
        int steps_noind = 0, steps_indon = 0;
        const float d_noind = run(Q, q_test, m_test, r_orb, cfl, N, 0, &n_noind, &steps_noind);
        const float d_indon = run(Q, q_test, m_test, r_orb, cfl, N, 1, &n_indon, &steps_indon);

        results[i].cfl = cfl;
        results[i].dt  = dt;
        results[i].n_steps_total = steps_indon;
        results[i].drift_noind   = d_noind;
        results[i].drift_indon   = d_indon;
        results[i].drift_indcontrib = (isfinite(d_noind) && isfinite(d_indon))
                                          ? (d_indon - d_noind) : NAN;
        results[i].n_completed_noind = n_noind;
        results[i].n_completed_indon = n_indon;

        char no_s[20], in_s[20], dc_s[20];
        if (!isfinite(d_noind)) snprintf(no_s, sizeof(no_s), "<unbound>");
        else snprintf(no_s, sizeof(no_s), "%+7.4f%%", 100.0 * (double) d_noind);
        if (!isfinite(d_indon)) snprintf(in_s, sizeof(in_s), "<unbound>");
        else snprintf(in_s, sizeof(in_s), "%+7.4f%%", 100.0 * (double) d_indon);
        if (!isfinite(results[i].drift_indcontrib)) snprintf(dc_s, sizeof(dc_s), "—");
        else snprintf(dc_s, sizeof(dc_s), "%+7.4f%%", 100.0 * (double) results[i].drift_indcontrib);

        printf("%-6d %-10.5f %-10.5f %-12d %-12s %-12s %-12s\n",
               i, (double) cfl, (double) dt, steps_indon, no_s, in_s, dc_s);
        fflush(stdout);
    }
    printf("\n");

    printf("Inductive contribution ratio (consecutive CFL halvings):\n");
    for (int i = 1; i < n_cfl; i++) {
        const float dprev = results[i - 1].drift_indcontrib;
        const float dcur  = results[i].drift_indcontrib;
        if (isfinite(dprev) && isfinite(dcur) && fabsf(dcur) > 1e-8f) {
            const float ratio = dprev / dcur;
            printf("  CFL %.5f -> %.5f:   |contrib_prev / contrib_cur| = %.3f\n",
                   (double) results[i - 1].cfl, (double) results[i].cfl,
                   fabs((double) ratio));
        }
    }
    printf("\n");
    printf("Interpretation (per Stage 29 gravity result):\n");
    printf("  ratio ~ 1.0 : CFL-independent (physical radiation reaction).\n");
    printf("  ratio ~ 4.0 : O(dt^2) numerical truncation.\n");

    for (int i = 0; i < n_cfl; i++) {
        TEST_ASSERT(results[i].n_completed_noind == N,
                    "no-inductive run at CFL=%.4f didn't complete %d orbits",
                    (double) cfls[i], N);
        TEST_ASSERT(results[i].n_completed_indon == N,
                    "inductive-on run at CFL=%.4f didn't complete %d orbits",
                    (double) cfls[i], N);
    }

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
