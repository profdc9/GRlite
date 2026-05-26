/* Stage 50 -- EM and GR binary side-by-side at maximally-parallel settings.
 *
 * Dan, 2026-05-26: "I want to make the EM and GR paths as close as
 * possible, so the main difference is that in EM you have opposite
 * charges and the Shapiro terms in FDTD, and in GR you have like
 * charges (mass).  This includes all force terms for both: phi,
 * dA/dt, magnetic, etc."
 *
 * Configuration mirroring:
 *
 *                          EM (pic_binary_em)           GR (pic_binary, updated)
 *   shape                  TSC                          TSC
 *   force interp           Lewis-Birdsall               Lewis-Birdsall
 *   bg_mode                SAMPLED                      SAMPLED
 *   field evolution        ON                           ON
 *   particle source dep    ON                           ON
 *   force pieces           -q grad phi  -q d_t A        -m grad Phi_g  -m d_t A_g
 *                          +q v x B                     +4 m v x B_g
 *   coupling const         k_e = 1                      G_eff = 1e-4
 *   per-particle coupling  q = +/-0.01                  m = +1.0   (same sign!)
 *   v_orb                  Q sqrt(k_e / m_inertia)      sqrt(G_eff * m)
 *   m_inertia              1.0                          1.0 (= coupling)
 *
 * The factor of 4 on v x B_g (spin-2 GEM enhancement) is genuine physics
 * and stays.  The Shapiro c_local^2 modification to the EM FDTD operator
 * is EM-only; gravity FDTD uses uniform c.  These are the two physics
 * differences Dan called out.
 *
 * Parameters chosen to make EM and GR forces equal in magnitude:
 *   F_em = k_e Q^2 / r = 1 * (0.01)^2 / 20 = 5e-6
 *   F_gr = G_eff m^2 / r = 1e-4 * 1^2 / 20 = 5e-6
 * v_orb = sqrt(F * r / m_inertia) = sqrt(5e-6 * 20 / 1) = 0.01 in both. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    int   is_em;
    float coupling;        /* k_e for EM, G_eff for GR */
    float per_particle;    /* Q for EM, m for GR */
    float mass_inertia;
    const char* shape_tag;
    gr_shape_function_t shape;
    float bump_R;
    int   smooth;
} cfg_t;

static void run(const cfg_t* c) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, c_eff = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float r_orb = 20.0f * dx;
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    /* Both: v_orb depends only on (coupling * per_particle^2 / mass_inertia) */
    const float v_orb = sqrtf(c->coupling * c->per_particle * c->per_particle / c->mass_inertia);
    const float T_ana = 2.0f * (float)M_PI * r_orb / v_orb;
    const int n_orbits = 4;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    if (c->is_em) {
        gr_sim_set_k_e(sim, c->coupling);
        const float params[6] = {c->mass_inertia, r_orb, 1.0f, cx, cy, c->per_particle};
        gr_sim_load_scenario(sim, "pic_binary_em", params, 6);
        /* EM scenario sets em_lorentz + em_inductive ON.  Disable gravity
         * to isolate pure EM physics. */
        gr_sim_set_G_eff(sim, 0.0f);
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    } else {
        gr_sim_set_G_eff(sim, c->coupling);
        const float params[5] = {c->per_particle, r_orb, 1.0f, cx, cy};
        gr_sim_load_scenario(sim, "pic_binary", params, 5);
        /* GR scenario (updated 2026-05-26) sets all gravitomagnetic pieces
         * ON to mirror EM.  Disable EM to isolate pure gravity. */
        gr_sim_set_k_e(sim, 0.0f);
        gr_sim_set_em_lorentz_force_enabled(sim, 0);
        gr_sim_set_em_inductive_enabled(sim, 0);
    }
    gr_sim_set_shape_function(sim, c->shape);
    if (c->shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, c->bump_R);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, c->smooth);
    gr_sim_set_j_smooth_passes(sim, c->smooth);

    const float dt = gr_sim_dt(sim);
    const int n_steps = (int)((float)n_orbits * T_ana / dt);

    float Lz0 = 0.0f, E0 = 0.0f;
    {
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(c->mass_inertia*c->mass_inertia));
        const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(c->mass_inertia*c->mass_inertia));
        const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
        const float d = sqrtf(dxs*dxs + dys*dys);
        /* Analytic 2D-log interaction PE:
         *   EM (opposite charges):  U = +2 k_e Q^2 ln(d)
         *   GR (same masses):       U = +2 G_eff m^2 ln(d)
         * Same form, just with the relevant coupling and per-particle value. */
        E0 = (g0 - 1.0f)*c->mass_inertia + (g1 - 1.0f)*c->mass_inertia
           + 2.0f * c->coupling * c->per_particle * c->per_particle * logf(d);
        Lz0 = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
            + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
    }
    float Px_max = 0.0f;
    float E_min = E0, E_max = E0, Lz_min = Lz0, Lz_max = Lz0;
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
            const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(c->mass_inertia*c->mass_inertia));
            const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(c->mass_inertia*c->mass_inertia));
            const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
            const float d = sqrtf(dxs*dxs + dys*dys);
            if (d < d_min) d_min = d;
            if (d > d_max) d_max = d;
            const float E = (g0 - 1.0f)*c->mass_inertia + (g1 - 1.0f)*c->mass_inertia
                          + 2.0f * c->coupling * c->per_particle * c->per_particle * logf(d);
            const float Lz = ((p0->x - cx)*p0->py - (p0->y - cy)*p0->px)
                           + ((p1->x - cx)*p1->py - (p1->y - cy)*p1->px);
            if (E < E_min) E_min = E;  if (E > E_max) E_max = E;
            if (Lz < Lz_min) Lz_min = Lz;
            if (Lz > Lz_max) Lz_max = Lz;
        }
    }
    if (unbound) {
        printf("  %-3s  %-12s   UNBOUND\n", c->is_em ? "EM" : "GR", c->shape_tag);
    } else {
        printf("  %-3s  %-12s  d=[%6.2f,%6.2f]  |Px|max=%.2e  dE=[%+.3f%%,%+.3f%%]  dL=[%+.3f%%,%+.3f%%]\n",
               c->is_em ? "EM" : "GR", c->shape_tag, (double)d_min, (double)d_max,
               (double)Px_max,
               (double)(100.0f*(E_min-E0)/fabsf(E0)),
               (double)(100.0f*(E_max-E0)/fabsf(E0)),
               (double)(100.0f*(Lz_min-Lz0)/fabsf(Lz0)),
               (double)(100.0f*(Lz_max-Lz0)/fabsf(Lz0)));
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage50_em_gr_parallel ===\n");
    printf("EM (k_e=1, Q=0.01, m_inertia=1.0, opposite charges)\n");
    printf("GR (G_eff=1e-4, m=1.0, m_inertia=1.0, same masses)\n");
    printf("Force F = 5e-6 at d=40 in both; v_orb=0.01.  All force terms ON in both.\n\n");

    const struct { const char* tag; gr_shape_function_t shape; float Rk; int smooth; } configs[4] = {
        {"TSC bare",   GR_SHAPE_TSC,  0.0f, 0},
        {"TSC + N=16", GR_SHAPE_TSC,  0.0f, 16},
        {"BUMP R=6",   GR_SHAPE_BUMP, 6.0f, 0},
        {"BUMP R=8",   GR_SHAPE_BUMP, 8.0f, 0},
    };
    for (int i = 0; i < 4; i++) {
        cfg_t cem = { 1, 1.0f,  0.01f, 1.0f, configs[i].tag, configs[i].shape, configs[i].Rk, configs[i].smooth };
        cfg_t cgr = { 0, 1e-4f, 1.0f,  1.0f, configs[i].tag, configs[i].shape, configs[i].Rk, configs[i].smooth };
        run(&cem);
        run(&cgr);
    }
    return 0;
}
