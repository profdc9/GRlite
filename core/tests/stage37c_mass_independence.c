/* Stage 37c -- is the spurious self-force mass-independent?
 *
 * Dan, 2026-05-27: "It's interesting the spurious force varies with
 * mass.  It must mean it's not Lorentz invariant."
 *
 * Claim under test: spurious force on a uniformly-moving charge depends
 * only on v (and q, grid state) -- NOT on the particle's mass.  If true,
 * the apparent mass-dependence of dKE / dvy in long runs is downstream:
 * smaller mass -> larger a = F/m -> particle traverses more of the F(v)
 * curve during the run -> different averaged drag.
 *
 * Test: measure F = m * a directly over a SHORT distance window, before
 * the particle's velocity has changed appreciably.  Sweep mass at fixed
 * initial v.  If F is constant across mass at each v, hypothesis confirmed.
 * If F varies with mass, there's a deeper bug.
 *
 * Equal-distance methodology throughout: n_steps = ceil(D / (v * dt)),
 * no caps.  Warmup distance = 8 cells; measurement distance = 2 cells. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float vy_of(const gr_particle_t* p, float c_eff) {
    const float m = p->mass;
    const float gamma = sqrtf(1.0f + (p->px*p->px + p->py*p->py) / (m*m*c_eff*c_eff));
    return p->py / (gamma * m);
}

typedef struct {
    float F_y;           /* spurious force = m * (dvy/dt), measured */
    float a_y;           /* dvy/dt averaged over the measurement window */
    float v_initial;     /* vy at start of measurement window */
    float v_final;       /* vy at end of measurement window */
    int   n_meas_steps;
} result_t;

static result_t measure(gr_shape_function_t shape, float Rk, int smooth,
                         float mass, float v_drift) {
    const int W = 128, H = 1024, n_damp = 16;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float Q = 0.01f;
    const float y_start = (float)(n_damp + 32) * dx;
    const float dt = cfl * dx / c_eff;

    /* Equal-distance windows.  warmup = 8 cells of travel, measurement
     * = 2 cells.  At v=0.05 measurement window = 57 steps; at v=0.1 = 29
     * steps; at v=0.2 = 15 steps.  Plenty for averaging over sub-cell
     * offsets in either x or y. */
    const int n_warmup = (int) ceilf(8.0f / (v_drift * dt));
    const int n_measure = (int) ceilf(2.0f / (v_drift * dt));

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, n_damp);
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, smooth);
    gr_sim_set_j_smooth_passes(sim, smooth);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);

    const float cx = ((float)(W-1) * 0.5f) * dx;
    gr_sim_add_particle(sim, cx, y_start, mass, Q, 0.0f, v_drift);

    /* Warmup -- let the wake form (8 cells of travel). */
    for (int s = 0; s < n_warmup; s++) gr_sim_step(sim);

    const gr_particle_t* p_w = gr_sim_get_particle(sim, 0);
    const float vy_at_warmup_end = vy_of(p_w, c_eff);

    /* Measure dvy over 2 cells of travel. */
    for (int s = 0; s < n_measure; s++) gr_sim_step(sim);
    const gr_particle_t* p_m = gr_sim_get_particle(sim, 0);
    const float vy_at_meas_end = vy_of(p_m, c_eff);

    result_t r;
    r.v_initial = vy_at_warmup_end;
    r.v_final   = vy_at_meas_end;
    r.a_y       = (vy_at_meas_end - vy_at_warmup_end) / ((float)n_measure * dt);
    r.F_y       = mass * r.a_y;
    r.n_meas_steps = n_measure;
    gr_sim_destroy(sim);
    return r;
}

static void scan(const char* kernel_label, gr_shape_function_t shape, float Rk, int smooth) {
    printf("=== %s ===\n", kernel_label);
    const float masses[] = {1.0f, 0.1f, 0.01f, 0.001f};
    const float vs[]     = {0.05f, 0.1f, 0.2f};
    const int Nm = 4, Nv = 3;

    printf("  v/c    mass=1.0       mass=0.1       mass=0.01      mass=0.001     (F = m*a, want same across mass)\n");
    printf("  ------------------------------------------------------------------------\n");
    for (int iv = 0; iv < Nv; iv++) {
        printf("  %.3f", (double)vs[iv]);
        float F_vals[4];
        for (int im = 0; im < Nm; im++) {
            result_t r = measure(shape, Rk, smooth, masses[im], vs[iv]);
            F_vals[im] = r.F_y;
            printf("  %+11.4e", (double)r.F_y);
        }
        /* Coefficient-of-variation: how much does F vary across mass for this v? */
        float mean = 0.0f;
        for (int im = 0; im < Nm; im++) mean += F_vals[im];
        mean /= (float)Nm;
        float var = 0.0f;
        for (int im = 0; im < Nm; im++) var += (F_vals[im] - mean)*(F_vals[im] - mean);
        var /= (float)Nm;
        const float cv = (fabsf(mean) > 0.0f) ? sqrtf(var) / fabsf(mean) : 0.0f;
        printf("    CoV=%.2f%%\n", (double)(100.0f * cv));
    }
    printf("\n");
}

int main(void) {
    printf("=== stage37c_mass_independence ===\n");
    printf("Spurious force F = m * a measured over 2-cell equal-distance window\n");
    printf("after 8-cell equal-distance warmup.  If F(v) is mass-independent,\n");
    printf("the four mass columns should agree to within float32 noise (CoV ~ few %%).\n\n");

    scan("TSC bare",   GR_SHAPE_TSC,  0.0f, 0);
    scan("TSC + N=16", GR_SHAPE_TSC,  0.0f, 16);
    scan("BUMP R=8",   GR_SHAPE_BUMP, 8.0f, 0);
    return 0;
}
