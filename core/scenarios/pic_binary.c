/* Scenario "pic_binary" — two equal-mass particles in a circular orbit
 * around their common center, with ALL forces coming from the FDTD
 * perturbation field (no analytic background).  This is the new-physics
 * scenario the v35 Yee + Esirkepov refactor was built to enable: PIC
 * grid-heating prevented bound mutual orbits with the v34 cell-centered
 * deposit, and Stage 11 (Esirkepov continuity) confirms the discrete
 * charge conservation that makes self-consistent FDTD dynamics tractable.
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: mass m of each particle  (default 0.01)
 *   params[1]: orbital radius r from COM (default 15.0 * dx)
 *   params[2]: v_factor                  (default 1.0; 0 = free-fall start)
 *   params[3]: cx of COM                 (default (W/2) * dx)
 *   params[4]: cy of COM                 (default (H/2) * dx)
 *
 * Particle velocities are v_factor * v_orb where v_orb = sqrt(G_eff m) is
 * the analytic circular-orbit speed for the 2D log potential.  v_factor=0
 * gives a free-fall release (mutual attraction pulls both particles in);
 * v_factor=1 attempts a circular orbit; intermediate values give elliptic
 * orbits.  Note that FDTD circular orbits are tricky to stabilize because
 * the field has to propagate between particles before the mutual force
 * appears (~separation/c wave-crossing time), during which the orbital
 * momentum carries the particles outward without an inward force; the
 * resulting elliptic-perturbed orbit may unbind for tight orbits.  Free-
 * fall avoids this issue and directly demonstrates the mutual attraction.
 *
 * For two equal masses m at separation d = 2r in the 2D linearized GR
 * potential Phi_g(r) = 2 G_eff m log(r) + const, Newton's law gives
 * mutual force F = 2 G_eff m^2 / d = G_eff m^2 / r.  Centripetal balance
 * for each particle at radius r from COM:
 *   m v^2 / r = G_eff m^2 / r   =>   v = sqrt(G_eff m)
 * Orbital period (if bound): T = 2 pi r / v = 2 pi r / sqrt(G_eff m). */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_pic_binary(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float mass     = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 0.01f;
    const float r_orb    = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 15.0f * dx;
    /* v_factor default 1.0 = the continuum analytic circular speed
     * sqrt(G m).  NOTE: with the BUMP R=8 kernel the discrete force during
     * motion is ~8% stronger, so the matched (near-circular) value is
     * v_factor ~ 1.05 -- the web demo passes that explicitly (stage63).
     * The default is left at 1.0 so existing tests that omit params[2]
     * (e.g. stage12) keep the documented continuum reference. */
    const float v_factor = (n_params >= 3) ? params[2] : 1.0f;
    /* True box center (W-1)/2 so the absorbing damping is symmetric
     * around the orbit / COM. */
    const float cx       = (n_params >= 4) ? params[3] : ((float) (W - 1) * 0.5f) * dx;
    const float cy       = (n_params >= 5) ? params[4] : ((float) (H - 1) * 0.5f) * dx;

    const float G_eff = gr_sim_get_G_eff(sim);
    const float v_orb = v_factor * sqrtf(G_eff * mass);

    /* Zero perturbation fields. */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    gr_sim_clear_sources(sim);
    gr_sim_clear_background(sim);
    gr_sim_clear_particles(sim);

    /* Vacuum perturbation-only setup. */
    gr_sim_set_bg_mode(sim, GR_BG_MODE_SAMPLED);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* Production defaults: BUMP kernel (R=8) + Lewis-Birdsall.  stage63
     * (2026-06-05) shows TSC deposition produces a catastrophic center-of-
     * mass drift for the symmetric binary (max |p0+p1| ~ 1.6 m*v_orb over
     * 4 orbits) and 130% separation breathing; switching to the wider,
     * smoother BUMP R=8 kernel pins the COM (max |p0+p1| ~ 4e-4 m*v_orb,
     * a ~3700x reduction) and cuts the breathing to ~20%.  See memory
     * [[grlite-self-vs-mutual-force-direction]] (BUMP R=6-8 is production)
     * and [[grlite-v42-orbit-stability]]. */
    gr_sim_set_shape_function(sim, GR_SHAPE_BUMP);
    gr_sim_set_kernel_radius(sim, 8.0f);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    /* Enable ALL gravity force pieces: -m grad Phi_g, +4 m v x B_g, and
     * -m d_t A_g.  This mirrors pic_binary_em's setup, where -q grad phi,
     * +q v x B, and -q d_t A are all on.  The first two are on by default
     * via gravitomagnetic_force_enabled = 1; the d_t A_g piece defaults
     * OFF (since most pre-v37 tests didn't want it), so we enable it
     * explicitly here.  Stage 28 verified that d_t A_g works correctly. */
    gr_sim_set_gravitomagnetic_force_enabled(sim, 1);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 1);

    /* Two particles counter-orbiting around (cx, cy).
     *   p0 at (cx - r, cy), velocity (0, +v)
     *   p1 at (cx + r, cy), velocity (0, -v)
     * gives a counter-clockwise orbit viewed from +z. */
    gr_sim_add_particle(sim, cx - r_orb, cy, mass, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/+v_orb);
    gr_sim_add_particle(sim, cx + r_orb, cy, mass, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/-v_orb);

    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_PIC_BINARY = {
    .name  = "pic_binary",
    .build = build_pic_binary,
};

void gr_scenario_register_pic_binary(void) {
    gr_scenario_register(&SCENARIO_PIC_BINARY);
}
