/* Stage 49 -- GR analog of Stage 48 Case C.
 *
 * Verifies that two-body GR binary (scenario pic_binary, gravity only)
 * at proper non-rel circular conditions conserves E and L_z at the
 * same level as the EM case (Stage 48 Case C).
 *
 * Setup parallels Stage 48 Case C:
 *   m = 1.0 (heavy enough to make orbit non-rel),
 *   G_eff = 1e-4 (so v_orb = sqrt(G_eff * m) = 0.01 = same as EM Case C),
 *   r = 20, v_factor = 1 (circular).
 *
 * V_eff inner turning point d_min = sqrt(L^2 / (2 m_red G_eff m^2))
 *                                 = sqrt((0.4)^2 / (2 * 0.5 * 1e-4 * 1))
 *                                 = sqrt(0.16 / 1e-4) = 40 cells = d_0.
 * So circular at d_0 = 40. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run(const char* tag, gr_shape_function_t shape, float Rk, int smooth) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass = 1.0f, r_orb = 20.0f * dx;
    const float G_eff_val = 1e-4f;
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float v_orb = sqrtf(G_eff_val * mass);
    const float T_ana = 2.0f * (float)M_PI * r_orb / v_orb;
    const int   n_orbits = 4;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    /* Set G_eff BEFORE scenario load: scenario uses it for v_orb. */
    gr_sim_set_G_eff(sim, G_eff_val);
    const float params[5] = {mass, r_orb, 1.0f, cx, cy};
    gr_sim_load_scenario(sim, "pic_binary", params, 5);
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, smooth);
    gr_sim_set_j_smooth_passes(sim, smooth);
    /* pic_binary now enables ALL gravity force pieces (grad Phi_g,
     * 4 m v x B_g, d_t A_g) to mirror pic_binary_em's setup.  Nothing
     * to toggle here. */

    const float dt = gr_sim_dt(sim);
    const int n_steps = (int)((float)n_orbits * T_ana / dt);
    float Lz0 = 0.0f, E0 = 0.0f;
    {
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(mass*mass));
        const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(mass*mass));
        const float dx0 = p1->x - p0->x, dy0 = p1->y - p0->y;
        const float d = sqrtf(dx0*dx0 + dy0*dy0);
        /* 2D log gravitational PE (per analog of EM): U = +2 G_eff m^2 ln(d) */
        E0 = (g0 - 1.0f)*mass + (g1 - 1.0f)*mass + 2.0f*G_eff_val*mass*mass*logf(d);
        Lz0 = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
            + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
    }
    float Px_max = 0.0f, Py_max = 0.0f;
    float E_min = E0, E_max = E0, Lz_min = Lz0, Lz_max = Lz0;
    float d_min_seen = 40.0f, d_max_seen = 40.0f;
    int unbound = 0;
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        if (!isfinite(p0->x) || !isfinite(p1->x)) { unbound = 1; break; }
        const float Px = p0->px + p1->px;
        const float Py = p0->py + p1->py;
        if (fabsf(Px) > Px_max) Px_max = fabsf(Px);
        if (fabsf(Py) > Py_max) Py_max = fabsf(Py);
        if ((s % 256) == 0) {
            const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(mass*mass));
            const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(mass*mass));
            const float dx0 = p1->x - p0->x, dy0 = p1->y - p0->y;
            const float d = sqrtf(dx0*dx0 + dy0*dy0);
            if (d < d_min_seen) d_min_seen = d;
            if (d > d_max_seen) d_max_seen = d;
            const float E = (g0 - 1.0f)*mass + (g1 - 1.0f)*mass + 2.0f*G_eff_val*mass*mass*logf(d);
            const float Lz = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
                           + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
            if (E < E_min) E_min = E;  if (E > E_max) E_max = E;
            if (Lz < Lz_min) Lz_min = Lz;  if (Lz > Lz_max) Lz_max = Lz;
        }
    }
    if (unbound) {
        printf("  %-16s UNBOUND (NaN)\n", tag);
    } else {
        printf("  %-16s d=[%.2f, %.2f]  |Px|max=%.2e  dE/|E0|=[%+.3f%%, %+.3f%%]  dL/|L0|=[%+.3f%%, %+.3f%%]\n",
               tag, (double)d_min_seen, (double)d_max_seen, (double)Px_max,
               (double)(100.0f*(E_min-E0)/fabsf(E0)),
               (double)(100.0f*(E_max-E0)/fabsf(E0)),
               (double)(100.0f*(Lz_min-Lz0)/fabsf(Lz0)),
               (double)(100.0f*(Lz_max-Lz0)/fabsf(Lz0)));
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage49_gr_binary_conservation (parallel of Stage 48 Case C) ===\n");
    printf("m=1.0, G_eff=1e-4, v_orb=0.01, r=20.  4 orbits.  PURE GRAVITY.\n");
    printf("All gravity force terms ON (grad Phi_g, 4 m v x B_g, d_t A_g).\n\n");
    run("TSC bare",     GR_SHAPE_TSC,  0.0f, 0);
    run("TSC + N=16",   GR_SHAPE_TSC,  0.0f, 16);
    run("BUMP R=6",     GR_SHAPE_BUMP, 6.0f, 0);
    run("BUMP R=8",     GR_SHAPE_BUMP, 8.0f, 0);
    return 0;
}
