/* Stage 63 -- pic_binary orbital stability across deposition/force knobs.
 *
 * With the v42 field init fixed (boundary-mean L-W + friction settle, so no
 * startup transient) and the Dirichlet+friction absorber, this isolates the
 * orbital instability itself: does the bound binary keep a steady separation,
 * and does the total momentum (COM) stay near zero?
 *
 * Sweeps the knobs the prior sessions identified as load-bearing for stable
 * PIC orbits (memory grlite-self-vs-mutual-force-direction, grlite-v39-self-
 * field, grlite-lewis-birdsall-result): shape function (TSC vs BUMP R), the
 * per-particle self-field subtraction, and v_factor.
 *
 * Metrics over n_orbits:
 *   sep_range  = (max_sep - min_sep) / mean_sep   (eccentricity + breathing)
 *   sep_growth = sep_final / sep_initial          (secular unbinding/infall)
 *   maxP       = max |p0 + p1| / (m v_orb)         (COM-drift / symmetry break)
 */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run(const char* label, gr_shape_function_t shape, float kernel_R,
                int self_field, float v_factor, int n_orbits) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const int   n_ring = 32;
    const float mass = 0.01f, r_orb = 15.0f;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_set_G_eff(sim, 1.0f);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_outer_bc_neumann(sim, 0);
    gr_sim_set_damping(sim, 0);
    gr_sim_set_volume_friction_taper(sim, 0.001f, 0.02f, n_ring);

    const float params[3] = { mass, r_orb, v_factor };
    if (gr_sim_load_scenario(sim, "pic_binary", params, 3) != 0) {
        printf("  %-34s LOAD FAILED\n", label); gr_sim_destroy(sim); return;
    }
    /* Override deposition/force knobs AFTER load (the scenario sets TSC+LB). */
    gr_sim_set_shape_function(sim, shape);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, kernel_R);
    if (self_field) {
        gr_sim_particle_enable_self_field(sim, 0);
        gr_sim_particle_enable_self_field(sim, 1);
    }

    gr_sim_init_potentials_settled(sim, W);   /* fixed init */

    const float v_orb   = sqrtf(1.0f * mass);
    const float dt      = cfl * dx / c_eff;
    const int   spo     = (int)((2.0 * M_PI * r_orb / v_orb) / dt);
    const float p_scale = mass * v_orb;

    const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
    const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
    float sx = p0->x - p1->x, sy = p0->y - p1->y;
    const float sep0 = sqrtf(sx * sx + sy * sy);

    float min_sep = sep0, max_sep = sep0, sum_sep = 0.0f, maxP = 0.0f, sep_final = sep0;
    long  n_samp = 0;
    const int n_steps = n_orbits * spo;
    for (int s = 1; s <= n_steps; s++) {
        gr_sim_step(sim);
        if (s % 50 == 0) {
            p0 = gr_sim_get_particle(sim, 0);
            p1 = gr_sim_get_particle(sim, 1);
            sx = p0->x - p1->x; sy = p0->y - p1->y;
            const float sep = sqrtf(sx * sx + sy * sy);
            if (sep < min_sep) min_sep = sep;
            if (sep > max_sep) max_sep = sep;
            sum_sep += sep; n_samp++; sep_final = sep;
            const float Px = p0->px + p1->px, Py = p0->py + p1->py;
            const float P = sqrtf(Px * Px + Py * Py);
            if (P > maxP) maxP = P;
        }
    }
    const float mean_sep = (n_samp > 0) ? sum_sep / (float) n_samp : sep0;
    printf("  %-34s range=%6.1f%%  growth=%5.2f  maxP/ps=%.3e\n",
           label,
           (double)(100.0f * (max_sep - min_sep) / mean_sep),
           (double)(sep_final / sep0),
           (double)(maxP / p_scale));
    gr_sim_destroy(sim);
}

int main(void) {
    const int n_orbits = 4;
    printf("=== stage63_orbit_stability ===\n");
    printf("pic_binary, fixed settled init, Dirichlet+friction, %d orbits\n", n_orbits);
    printf("  range = (max-min)/mean separation; growth = sep_final/sep_0;\n");
    printf("  maxP  = max |p0+p1| / (m v_orb)  (0 = COM stays put)\n\n");

    run("TSC,        no-self, v=1.00", GR_SHAPE_TSC,  0.0f, 0, 1.00f, n_orbits);
    run("BUMP R=8,   no-self, v=1.00", GR_SHAPE_BUMP, 8.0f, 0, 1.00f, n_orbits);
    run("BUMP R=8,   self,    v=1.00", GR_SHAPE_BUMP, 8.0f, 1, 1.00f, n_orbits);
    run("BUMP R=8,   self,    v=0.98", GR_SHAPE_BUMP, 8.0f, 1, 0.98f, n_orbits);
    run("BUMP R=8,   self,    v=0.95", GR_SHAPE_BUMP, 8.0f, 1, 0.95f, n_orbits);
    return 0;
}
