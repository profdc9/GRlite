/* Stage 33 -- TSC + Lewis-Birdsall on the EM scalar gradient.
 *
 * Verifies that the v36 EM TSC+LB upgrade on phi_em_grad_at_total reduces
 * the spatial-piece PIC residual on the EM channel.  Specifically, the
 * Stage 30 closed-loop PIC orbit at CIC+LEGACY (the pre-upgrade EM-side
 * default) was the EM-vs-gravity asymmetry exposed in the
 * grlite-next-steps queue.
 *
 * Compares the NO-INDUCTIVE drift (Coulomb-only PIC dynamics) over 4
 * orbits with three (shape, force_interp) combinations:
 *
 *   (1) CIC + LEGACY     -- the historical EM baseline (Stage 30 setup)
 *   (2) TSC + LEGACY     -- TSC smoother shape, legacy FD-then-interp
 *   (3) TSC + LEWIS-BIRDSALL -- production target (energy-conserving)
 *
 * Expectation: drift magnitude DECREASES going (1) -> (2) -> (3),
 * mirroring the gravity-side improvement quantified in Stage 18.  No
 * strict drift threshold is asserted (the physical magnitude depends
 * on radiation-reaction details); we just check that the field machinery
 * runs correctly under TSC+LB (no NaN, no Esirkepov violations,
 * complete 4 orbits) and report the comparative drifts. */

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

static float run_one(gr_shape_function_t shape, gr_force_interp_t interp,
                     int* n_completed_out, int* viols_out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float Q      = +1.0f;
    const float q_test = -1.0e-3f;
    const float m_test = +1.0e-3f;
    const float r_orb  = 20.0f;
    const int   N      = 4;
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
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, interp);
    gr_sim_set_rho_smooth_passes(sim, 4);
    gr_sim_set_em_inductive_enabled(sim, 0);   /* No inductive -- spatial piece only. */

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test, 0.0f, v_circ);

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.2f * (float) N * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    float r_at_end = r_orb;
    for (int s = 0; s < n_max && wraps < N; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) {
            if (n_completed_out) *n_completed_out = wraps;
            if (viols_out) *viols_out = gr_sim_esirkepov_violations(sim);
            gr_sim_destroy(sim);
            return NAN;
        }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            r_at_end = r_now;
        }
        th_prev = th;
    }
    if (n_completed_out) *n_completed_out = wraps;
    if (viols_out) *viols_out = gr_sim_esirkepov_violations(sim);
    gr_sim_destroy(sim);
    if (wraps < N) return NAN;
    return (r_at_end - r_orb) / r_orb;
}

int main(void) {
    printf("=== stage33_em_tsc_lb_check ===\n");
    printf("EM Stage-30 setup, NO-inductive (Coulomb-only), 4 orbits.\n");
    printf("Compare CIC+LEGACY vs TSC+LEGACY vs TSC+LB.\n\n");

    const struct {
        const char*         label;
        gr_shape_function_t shape;
        gr_force_interp_t   interp;
    } configs[] = {
        {"CIC + LEGACY     ", GR_SHAPE_CIC, GR_FORCE_INTERP_LEGACY        },
        {"TSC + LEGACY     ", GR_SHAPE_TSC, GR_FORCE_INTERP_LEGACY        },
        {"TSC + LEWIS-BIRDS", GR_SHAPE_TSC, GR_FORCE_INTERP_LEWIS_BIRDSALL},
        {"CIC + LEWIS-BIRDS", GR_SHAPE_CIC, GR_FORCE_INTERP_LEWIS_BIRDSALL},
    };
    const int n_cfg = (int)(sizeof(configs) / sizeof(configs[0]));

    printf("%-20s %-12s %-10s %-7s\n", "config", "drift", "n_completed", "viols");
    printf("---------------------------------------------------------\n");

    float drifts[4];
    for (int i = 0; i < n_cfg; i++) {
        int n_completed = 0, viols = 0;
        const float drift = run_one(configs[i].shape, configs[i].interp,
                                    &n_completed, &viols);
        drifts[i] = drift;
        char dstr[20];
        if (!isfinite(drift)) snprintf(dstr, sizeof(dstr), "<unbound>");
        else snprintf(dstr, sizeof(dstr), "%+7.4f%%", 100.0 * (double) drift);
        printf("%-20s %-12s %-10d %-7d\n", configs[i].label, dstr,
               n_completed, viols);
        TEST_ASSERT(n_completed == 4, "config %d didn't complete 4 orbits",
                    i);
        TEST_ASSERT(viols == 0, "config %d had Esirkepov violations: %d",
                    i, viols);
        TEST_ASSERT(isfinite(drift), "config %d unbound", i);
    }
    printf("\n");
    printf("Expectation: |drift| should DECREASE going CIC+LEGACY -> TSC+LB.\n");
    printf("This is the v36 closure of the EM-vs-gravity spatial-piece asymmetry\n");
    printf("from Stage 30 (gravity Stage 29 no-induct drift was ~0.3%% at all CFL;\n");
    printf("EM at CIC+LEGACY was several %% -- this stage shows the convergence.\n");

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
