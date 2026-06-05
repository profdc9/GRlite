/* Stage 62 -- Liénard-Wiechert field-init quality.
 *
 * Two clean, decoupled measurements of how good the v42 L-W initialization
 * (with per-particle boundary-mean subtraction, 2026-06-05) is:
 *
 *   Part A -- frozen-particle field transient (pure init quality).
 *     Load pic_static (single mass at box center), FREEZE the particle so
 *     no dynamics occur, init the field, and run.  If the IC equals the
 *     discrete fixed point the field never moves; any motion is the
 *     init mismatch (near-source discrete-vs-continuum + boundary residual)
 *     being shed as a transient that friction then damps.  Metric: peak
 *     |Phi_g - Phi_g(0)| over the grid, relative to |Phi_g| at the source.
 *
 *   Part B -- total-momentum (COM-impulse) drift for a symmetric binary.
 *     Load pic_binary.  The setup is symmetric under point inversion through
 *     the box center, so the total particle momentum P = p0 + p1 is zero by
 *     symmetry and MUST stay zero under internal forces -- any nonzero P is
 *     a symmetry-breaking numerical artifact (the COM-drift cause).  Unlike
 *     COM position, P is not confounded by orbital eccentricity/instability
 *     (internal forces never change it).  Metric: max |P| over a couple of
 *     orbits, relative to a single particle's momentum scale m*v_orb. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void configure_absorber(gr_sim_t* sim, int n_ring) {
    /* Mirrors web/src/main.ts bindApi: Dirichlet outer BC, multiplicative
     * ring OFF, centered-implicit derivative friction ON. */
    gr_sim_set_outer_bc_neumann(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_volume_friction_taper(sim, 0.001f, 0.02f, n_ring);
}

static int part_a_frozen_transient(void) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const int   n_ring = 32;
    const float mass = 0.01f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_G_eff(sim, 1.0f);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    configure_absorber(sim, n_ring);

    const float params[1] = { mass };
    if (gr_sim_load_scenario(sim, "pic_static", params, 1) != 0) {
        printf("  (pic_static load failed)\n"); gr_sim_destroy(sim); return 1;
    }
    gr_sim_init_potentials_lienard_wiechert(sim);

    /* Snapshot the initial field. */
    const size_t ncells = (size_t) W * (size_t) H;
    float* phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
    float phi0_min = 1e30f, phi0_max = -1e30f;
    for (size_t k = 0; k < ncells; k++) {
        if (phi[k] < phi0_min) phi0_min = phi[k];
        if (phi[k] > phi0_max) phi0_max = phi[k];
    }
    const float phi_scale = phi0_max - phi0_min;
    static float phi0[256 * 256];
    for (size_t k = 0; k < ncells; k++) phi0[k] = phi[k];

    /* Freeze the particle: only the field evolves. */
    gr_sim_set_particles_frozen(sim, 1);

    printf("=== Part A: frozen-particle field transient (pic_static) ===\n");
    printf("  Phi_g range at init = [%.4e, %.4e]  (scale %.4e)\n",
           (double) phi0_min, (double) phi0_max, (double) phi_scale);
    printf("  step    max|dPhi_g| / scale\n");

    float worst = 0.0f;
    const int n_steps = 3000;
    for (int s = 1; s <= n_steps; s++) {
        gr_sim_step(sim);
        if (s % 300 == 0) {
            phi = gr_sim_field_ptr(sim, GR_FIELD_PHI_GRAV);
            float md = 0.0f;
            for (size_t k = 0; k < ncells; k++) {
                const float d = fabsf(phi[k] - phi0[k]);
                if (d > md) md = d;
            }
            if (md > worst) worst = md;
            printf("  %5d    %.4e\n", s, (double)(md / phi_scale));
        }
    }
    printf("  worst transient / scale = %.4e (%.2f%%)\n\n",
           (double)(worst / phi_scale), (double)(100.0f * worst / phi_scale));
    gr_sim_destroy(sim);
    return 0;
}

static float part_b_momentum_drift(float v_factor, int presettle_steps,
                                    const char* label) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const int   n_ring = 32;
    const float mass = 0.01f, r_orb = 15.0f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_G_eff(sim, 1.0f);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    configure_absorber(sim, n_ring);

    const float params[3] = { mass, r_orb, v_factor };
    if (gr_sim_load_scenario(sim, "pic_binary", params, 3) != 0) {
        printf("  (pic_binary load failed)\n"); gr_sim_destroy(sim); return -1.0f;
    }
    /* Use the shipped helper: boundary-mean L-W init + n_settle frozen
     * friction settle.  presettle_steps == 0 reduces to a plain init. */
    gr_sim_init_potentials_settled(sim, presettle_steps);

    const float v_orb = v_factor * sqrtf(1.0f * mass);
    const float dt    = cfl * dx / c_eff;
    const int   spo   = (int)((2.0 * M_PI * r_orb / (sqrtf(mass))) / dt);
    const float p_scale = mass * (v_orb > 0.0f ? v_orb : sqrtf(mass));

    float worst = 0.0f, early = 0.0f;
    const int n_steps = 2 * spo;
    for (int s = 1; s <= n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        const float Px = p0->px + p1->px, Py = p0->py + p1->py;
        const float P  = sqrtf(Px * Px + Py * Py);
        if (P > worst) worst = P;
        if (s == spo / 4 && early == 0.0f) early = P;   /* ~quarter orbit */
    }
    printf("  %-32s worst|P|/ps=%.3e   quarter-orbit|P|/ps=%.3e\n",
           label, (double)(worst / p_scale), (double)(early / p_scale));
    gr_sim_destroy(sim);
    return worst / p_scale;
}

int main(void) {
    printf("=== stage62_lw_init_com_drift ===\n\n");
    part_a_frozen_transient();

    printf("=== Part B: total-momentum (COM-impulse) drift ===\n");
    printf("  P=0 by point-inversion symmetry; any |P|>0 is symmetry-breaking.\n");
    printf("  Discriminator: presettle pins the field at its fixed point before\n");
    printf("  release, so a drop with presettle => init-transient cause;\n");
    printf("  no drop => orbital-dynamics cause (separate v_factor issue).\n\n");
    part_b_momentum_drift(1.00f, 0,    "v=1.00, no presettle");
    part_b_momentum_drift(1.00f, 2000, "v=1.00, presettle 2000");
    part_b_momentum_drift(0.98f, 0,    "v=0.98, no presettle");
    part_b_momentum_drift(0.98f, 2000, "v=0.98, presettle 2000");

    printf("\n=== Part C: presettle-count sweep (quarter-orbit init impulse) ===\n");
    printf("  Pick the smallest settle that drives the init impulse to ~round-off.\n");
    part_b_momentum_drift(1.00f, 0,   "v=1.00, settle    0");
    part_b_momentum_drift(1.00f, 100, "v=1.00, settle  100");
    part_b_momentum_drift(1.00f, 200, "v=1.00, settle  200");
    part_b_momentum_drift(1.00f, 400, "v=1.00, settle  400");
    part_b_momentum_drift(1.00f, 800, "v=1.00, settle  800");

    /* Measurement, not a pass/fail gate. */
    return 0;
}
