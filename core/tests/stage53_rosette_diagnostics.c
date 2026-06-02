/* Stage 53 -- isolate the residual 7% rosette in Stage 52.
 *
 * Dan, 2026-06-01: "Investigate the rosette... it is a likely
 * combination of Bertrand's theorem artifacts and each mass/charge
 * responding to the delayed potential of the other."
 *
 * Two-pronged hypothesis to test:
 *   (a) Bertrand: 2D log Coulomb is not Bertrand-closing, so any IC
 *       perturbation gives a non-closed precessing rosette.
 *   (b) Delayed potential: at t=0 the field is zero (no mutual force).
 *       It takes ~ d/c sim units for the field to build up.  During
 *       this time each particle drifts radially outward (no centripetal
 *       force), leaving the orbit off-circular by the time the field
 *       has equilibrated.  This off-circular IC then drives the
 *       Bertrand rosette via mechanism (a).
 *
 * Test: vary c_eff.  v_orb does NOT depend on c, so doubling c keeps
 * the orbit at the same radius and period but halves tau = d/c.
 * Smaller tau/T means smaller field-buildup IC perturbation.
 *
 * If the rosette amplitude scales as 1/c, mechanism (b) dominates.
 * If c-independent, mechanism (a) acting on some other source dominates.
 *
 * Self-field ENABLED (spurious self-force removed).  BUMP R=8.  2 orbits. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run(const char* tag, int is_em, float c_eff) {
    const int   W = 256, H = 256;
    const float dx = 1.0f, cfl = 1.0f / sqrtf(2.0f);
    const float mass_inertia = 1.0f, Q = 0.01f, r_orb = 20.0f * dx;
    const float coupling = is_em ? 1.0f : 1e-4f;
    const float per_particle = is_em ? Q : mass_inertia;
    const float cx = ((float)(W - 1) * 0.5f) * dx;
    const float cy = ((float)(H - 1) * 0.5f) * dx;
    const float v_orb = sqrtf(coupling * per_particle * per_particle / mass_inertia);
    const float T_ana = 2.0f * (float)M_PI * r_orb / v_orb;
    const float tau   = 2.0f * r_orb / c_eff;
    const int n_orbits = 2;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);
    if (is_em) {
        gr_sim_set_k_e(sim, coupling);
        const float params[6] = {mass_inertia, r_orb, 1.0f, cx, cy, Q};
        gr_sim_load_scenario(sim, "pic_binary_em", params, 6);
        gr_sim_set_G_eff(sim, 0.0f);
        gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
        gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);
    } else {
        gr_sim_set_G_eff(sim, coupling);
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
    gr_sim_particle_enable_self_field(sim, 0);
    gr_sim_particle_enable_self_field(sim, 1);

    const float dt = gr_sim_dt(sim);
    const int n_steps = (int)((float)n_orbits * T_ana / dt);
    float d_min = 40.0f, d_max = 40.0f;
    int unbound = 0;
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        if (!isfinite(p0->x)) { unbound = 1; break; }
        if ((s % 64) == 0) {
            const float dxs = p1->x - p0->x, dys = p1->y - p0->y;
            const float d = sqrtf(dxs*dxs + dys*dys);
            if (d < d_min) d_min = d;
            if (d > d_max) d_max = d;
        }
    }
    if (unbound) {
        printf("  %-3s c=%4.1f  tau/T=%.5f   UNBOUND\n",
               tag, (double)c_eff, (double)(tau/T_ana));
    } else {
        const float dlow  = 100.0f * (d_min - 40.0f) / 40.0f;
        const float dhigh = 100.0f * (d_max - 40.0f) / 40.0f;
        printf("  %-3s c=%4.1f  tau/T=%.5f   d_min=%6.2f (%+5.2f%%)  d_max=%6.2f (%+5.2f%%)  range=%5.2f%%\n",
               tag, (double)c_eff, (double)(tau/T_ana),
               (double)d_min, (double)dlow,
               (double)d_max, (double)dhigh,
               (double)(dhigh - dlow));
    }
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage53_rosette_diagnostics ===\n");
    printf("Binary orbit, self-field ON (eps=0), BUMP R=8.  2 orbits each.\n");
    printf("Vary c_eff to test delayed-potential hypothesis.\n\n");

    run("EM", 1, 1.0f);
    run("EM", 1, 2.0f);
    run("EM", 1, 4.0f);
    run("GR", 0, 1.0f);
    run("GR", 0, 2.0f);
    run("GR", 0, 4.0f);
    return 0;
}
