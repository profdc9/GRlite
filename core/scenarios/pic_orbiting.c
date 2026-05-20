/* Scenario "pic_orbiting" — single particle in a Stage 7-style circular
 * orbit around an analytic softened point-mass background, but with the
 * perturbation FDTD running and auto-deposition enabled.  Stage 10
 * orbit-regression test (test E) uses this to verify that the additional
 * perturbation forces from the particle's own retarded field do not
 * significantly alter the orbit over a few periods.
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: GM (central mass * G_eff)  (default 1.0)
 *   params[1]: orbital radius r           (default 20.0)
 *   params[2]: smoothing length eps       (default 1.0)
 *   params[3]: orbiting particle mass     (default 0.001 — weak coupling)
 *   params[4]: x0  (default width *dx / 2)
 *   params[5]: y0  (default height*dx / 2)
 *
 * Note on the default mass: the orbiting particle's mass enters the test
 * only through (a) its auto-deposited perturbation source and (b) the
 * radiation reaction back on itself.  Both scale with the deposit; the
 * test-particle orbit dynamics around the BACKGROUND mass don't depend on
 * the test particle's mass at all.  At m_test = 1 and GM_bg = 1, the
 * deposit is strong enough that grid-heating / self-force artifacts
 * overwhelm the orbit on timescales much shorter than the orbital period.
 * m_test = 0.001 puts us in the weak-coupling regime where the doc's
 * Stage 10 assumptions actually hold.
 *
 * Identical initial conditions to kepler_orbit; only the post-load
 * configuration differs (perturbation FDTD + auto-deposition both on,
 * analytic background still on for the central mass). */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_pic_orbiting(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float GM     = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 1.0f;
    const float r      = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 20.0f;
    const float eps    = (n_params >= 3 && params[2] > 0.0f) ? params[2] : 1.0f;
    const float m_test = (n_params >= 4 && params[3] > 0.0f) ? params[3] : 0.001f;
    const float x0     = (n_params >= 5) ? params[4] : ((float) W * 0.5f) * dx;
    const float y0     = (n_params >= 6) ? params[5] : ((float) H * 0.5f) * dx;

    /* Zero perturbation fields and sources. */
    const size_t n = (size_t) W * (size_t) H;
    for (int f = 0; f < GR_FIELD_COUNT; f++) {
        memset(sim->fields[f].prev, 0, n * sizeof(float));
        memset(sim->fields[f].curr, 0, n * sizeof(float));
        memset(sim->fields[f].next, 0, n * sizeof(float));
    }
    gr_sim_clear_sources(sim);
    gr_sim_clear_background(sim);
    gr_sim_set_background_point_mass(sim, x0, y0, GM, eps);
    gr_sim_clear_particles(sim);

    /* Stage 10 configuration: analytic background for the central mass,
     * perturbation FDTD on for the orbit's radiated field, auto-deposition
     * on so the particle sources the perturbation. */
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    /* Apply binomial smoothing on rho to suppress moving-particle deposit
     * aliasing — the dominant Tier-0 PIC heating mode.  Empirically this
     * lifts the bound-orbit coupling threshold from m_test ~ 1e-6 (Yee+
     * Esirkepov alone) to m_test ~ 1e-4 holding 4+ orbits cleanly at 4
     * passes.  Cost: 4 light low-pass-filter sweeps per step, negligible. */
    gr_sim_set_rho_smooth_passes(sim, 4);

    /* Relativistic circular velocity for the softened force law (identical
     * to kepler_orbit's derivation). */
    const float r2_e2_15 = powf(r * r + eps * eps, 1.5f);
    const float g_mag    = GM * r / r2_e2_15;
    const float rg       = r * g_mag;
    const float rg2_c2   = rg * rg / (sim->c_eff * sim->c_eff);
    const float u        = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ   = sqrtf(u);

    gr_sim_add_particle(sim, x0 + r, y0,
                        m_test, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/v_circ);
    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_PIC_ORBITING = {
    .name  = "pic_orbiting",
    .build = build_pic_orbiting,
};

void gr_scenario_register_pic_orbiting(void) {
    gr_scenario_register(&SCENARIO_PIC_ORBITING);
}
