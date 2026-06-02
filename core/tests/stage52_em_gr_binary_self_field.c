/* Stage 52 -- binary orbit with v39 per-particle self-field subtraction.
 *
 * Setup mirrors Stage 50 (EM and GR binary at matched coupling
 * magnitudes, BUMP R=8, mass_inertia=1, v_orb=0.01).  For each system
 * we run two configurations:
 *   (1) self-field DISABLED -- baseline = current PIC behavior
 *   (2) self-field ENABLED on BOTH particles, eps=0
 *
 * Predicted: with self-field on, each particle sees only the field
 * from the OTHER particle.  The +13% EM rosette in Stage 50 should
 * collapse to ~0 (because the rosette came from float-roundoff in the
 * bipolar rho_q's dipole cancellation + the spurious self-force).
 * GR's already-tight orbit should stay clean. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { int is_em; int self_on; } cfg_t;

static void run(cfg_t c) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass_inertia = 1.0f, Q = 0.01f, r_orb = 20.0f * dx;
    const float k_e_or_G_eff = c.is_em ? 1.0f : 1e-4f;
    const float per_particle = c.is_em ? Q : mass_inertia;
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float v_orb = sqrtf(k_e_or_G_eff * per_particle * per_particle / mass_inertia);
    const float T_ana = 2.0f * (float)M_PI * r_orb / v_orb;
    const int n_orbits = 4;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    if (c.is_em) {
        gr_sim_set_k_e(sim, k_e_or_G_eff);
        const float params[6] = {mass_inertia, r_orb, 1.0f, cx, cy, Q};
        gr_sim_load_scenario(sim, "pic_binary_em", params, 6);
        gr_sim_set_G_eff(sim, 0.0f);
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    } else {
        gr_sim_set_G_eff(sim, k_e_or_G_eff);
        const float params[5] = {mass_inertia, r_orb, 1.0f, cx, cy};
        gr_sim_load_scenario(sim, "pic_binary", params, 5);
        gr_sim_set_k_e(sim, 0.0f);
        gr_sim_set_em_lorentz_force_enabled(sim, 0);
        gr_sim_set_em_inductive_enabled(sim, 0);
    }
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 8.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);

    if (c.self_on) {
        gr_sim_particle_enable_self_field(sim, 0);
        gr_sim_particle_enable_self_field(sim, 1);
    }

    const float dt = gr_sim_dt(sim);
    const int n_steps = (int)((float)n_orbits * T_ana / dt);
    float Lz0 = 0.0f, E0 = 0.0f;
    {
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(mass_inertia*mass_inertia));
        const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(mass_inertia*mass_inertia));
        const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
        const float d = sqrtf(dxs*dxs + dys*dys);
        E0 = (g0 - 1.0f)*mass_inertia + (g1 - 1.0f)*mass_inertia
           + 2.0f * k_e_or_G_eff * per_particle * per_particle * logf(d);
        Lz0 = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
            + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
    }
    float Px_max = 0.0f, E_min = E0, E_max = E0, Lz_min = Lz0, Lz_max = Lz0;
    float d_min = 40.0f, d_max = 40.0f;
    int unbound = 0;
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        if (!isfinite(p0->x) || !isfinite(p1->x)) { unbound = 1; break; }
        const float Px = p0->px + p1->px;
        if (fabsf(Px) > Px_max) Px_max = fabsf(Px);
        if ((s % 256) == 0) {
            const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(mass_inertia*mass_inertia));
            const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(mass_inertia*mass_inertia));
            const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
            const float d = sqrtf(dxs*dxs + dys*dys);
            if (d < d_min) d_min = d;
            if (d > d_max) d_max = d;
            const float E = (g0 - 1.0f)*mass_inertia + (g1 - 1.0f)*mass_inertia
                          + 2.0f * k_e_or_G_eff * per_particle * per_particle * logf(d);
            const float Lz = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
                           + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
            if (E < E_min) E_min = E;  if (E > E_max) E_max = E;
            if (Lz < Lz_min) Lz_min = Lz;
            if (Lz > Lz_max) Lz_max = Lz;
        }
    }
    if (unbound) {
        printf("  %-3s self=%s    UNBOUND\n", c.is_em ? "EM" : "GR", c.self_on ? "on " : "off");
    } else {
        printf("  %-3s self=%s   d=[%6.2f,%6.2f]  |Px|max=%.2e  dE=[%+.3f%%,%+.3f%%]  dL=[%+.3f%%,%+.3f%%]\n",
               c.is_em ? "EM" : "GR", c.self_on ? "on " : "off",
               (double)d_min, (double)d_max, (double)Px_max,
               (double)(100.0f*(E_min-E0)/fabsf(E0)),
               (double)(100.0f*(E_max-E0)/fabsf(E0)),
               (double)(100.0f*(Lz_min-Lz0)/fabsf(Lz0)),
               (double)(100.0f*(Lz_max-Lz0)/fabsf(Lz0)));
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage52_em_gr_binary_self_field ===\n");
    printf("Binary orbit at matched configs (Stage 50 setup), BUMP R=8.\n");
    printf("self=off: baseline.  self=on: per-particle field-set enabled, eps=0.\n");
    printf("Prediction: self=on collapses EM's +13%% rosette to GR-like cleanliness.\n\n");
    cfg_t cfgs[] = {{1, 0}, {1, 1}, {0, 0}, {0, 1}};
    for (int i = 0; i < 4; i++) run(cfgs[i]);
    return 0;
}
