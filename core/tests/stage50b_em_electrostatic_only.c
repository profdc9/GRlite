/* Stage 50b -- isolate the +13% EM outward drift from Stage 50.
 *
 * Dan, 2026-05-26: "Let's verify your hypothesis by using only
 * electrostatic force in the simulation."
 *
 * Stage 50 ran EM with all three force terms ON (-q grad phi, -q d_t A,
 * +q v x B).  Saw a +13% outward drift in d_max while GR (with the
 * analog +4 m v x B_g and -m d_t A_g and -m grad Phi_g all ON) stayed
 * at d=40 +/- 0.1%.  Hypothesis: the J^{n-1/2} -> J^n half-step timing
 * mismatch in Esirkepov creates a small numerical "anti-radiation
 * reaction" in the -q d_t A channel, slowly adding energy.
 *
 * This test isolates the channels.  For BUMP R=8 (the cleanest config):
 *   (A) full EM            : grad phi + d_t A + v x B   (Stage 50 baseline)
 *   (B) electrostatic only : grad phi                   (NO d_t A, NO v x B)
 *   (C) no inductive       : grad phi + v x B           (NO d_t A)
 *   (D) no magnetic        : grad phi + d_t A           (NO v x B)
 *
 * If hypothesis is right: (B) and (D) should be clean (no outward drift),
 * (C) might still drift if v x B is contributing too, (A) baseline matches
 * what Stage 50 saw. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run(const char* tag, int electrostatic, int magnetic, int inductive) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass_inertia = 1.0f, Q = 0.01f, r_orb = 20.0f * dx, k_e = 1.0f;
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float v_orb = Q * sqrtf(k_e / mass_inertia);
    const float T_ana = 2.0f * (float)M_PI * r_orb / v_orb;
    const int n_orbits = 4;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    const float params[6] = {mass_inertia, r_orb, 1.0f, cx, cy, Q};
    gr_sim_load_scenario(sim, "pic_binary_em", params, 6);
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 8.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    /* The dial knobs for this experiment: */
    gr_sim_set_em_electrostatic_enabled(sim, electrostatic);
    gr_sim_set_em_magnetic_enabled     (sim, magnetic);
    gr_sim_set_em_inductive_enabled    (sim, inductive);

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
        E0 = (g0 - 1.0f)*mass_inertia + (g1 - 1.0f)*mass_inertia + 2.0f*k_e*Q*Q*logf(d);
        Lz0 = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
            + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
    }
    float Px_max = 0.0f, E_min = E0, E_max = E0, Lz_min = Lz0, Lz_max = Lz0;
    float d_min = 40.0f, d_max = 40.0f;
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        if (!isfinite(p0->x)) { printf("  %-25s NaN\n", tag); gr_sim_destroy(sim); return; }
        const float Px = p0->px + p1->px;
        if (fabsf(Px) > Px_max) Px_max = fabsf(Px);
        if ((s % 256) == 0) {
            const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(mass_inertia*mass_inertia));
            const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(mass_inertia*mass_inertia));
            const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
            const float d = sqrtf(dxs*dxs + dys*dys);
            if (d < d_min) d_min = d;
            if (d > d_max) d_max = d;
            const float E = (g0 - 1.0f)*mass_inertia + (g1 - 1.0f)*mass_inertia + 2.0f*k_e*Q*Q*logf(d);
            const float Lz = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
                           + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
            if (E < E_min) E_min = E;  if (E > E_max) E_max = E;
            if (Lz < Lz_min) Lz_min = Lz;
            if (Lz > Lz_max) Lz_max = Lz;
        }
    }
    printf("  %-25s d=[%6.2f,%6.2f]  Px=%.2e  dE=[%+.3f%%,%+.3f%%]  dL=[%+.3f%%,%+.3f%%]\n",
           tag, (double)d_min, (double)d_max, (double)Px_max,
           (double)(100.0f*(E_min-E0)/fabsf(E0)),
           (double)(100.0f*(E_max-E0)/fabsf(E0)),
           (double)(100.0f*(Lz_min-Lz0)/fabsf(Lz0)),
           (double)(100.0f*(Lz_max-Lz0)/fabsf(Lz0)));
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage50b_em_electrostatic_only ===\n");
    printf("EM binary (k_e=1, Q=0.01, m_inertia=1, v_orb=0.01, r=20).  BUMP R=8.\n");
    printf("Isolate which EM force term causes the +13%% outward drift.\n\n");
    /*    tag                       elec  mag  ind */
    run("(A) full EM",              1, 1, 1);
    run("(B) electrostatic only",   1, 0, 0);
    run("(C) no inductive",         1, 1, 0);
    run("(D) no magnetic",          1, 0, 1);
    return 0;
}
