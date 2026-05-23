/* Scenario "pic_binary_em" -- EM analog of pic_binary.
 *
 * Two equal-mass, opposite-charge (+Q, -Q) particles in circular orbit
 * around their common center, with ALL EM forces coming from the FDTD
 * perturbation field (no analytic background).  Each particle deposits
 * (rho_q, J_qx, J_qy) via Esirkepov; each feels the other's perturbation
 * phi_em + A_em via the full EM Lorentz force F_em = q(-grad phi - d_t A
 * + v x B).  The inductive piece is production-ON (Stage 27/30 resolution).
 *
 * SIGN CONVENTION (standard Maxwell after the v36 c_em flip):
 *   PDE static limit : Lap phi = -4 pi k_e rho_q
 *   2D log Green's   : phi_pert(r) = -2 k_e Q ln(r) + const
 *   F on q'          : -q' grad phi = +2 k_e q q' / r * r_hat
 *                      => REPULSIVE if q q' > 0 (like charges),
 *                         ATTRACTIVE if q q' < 0 (opposite charges).
 * The optional analytic point-charge background
 * (gr_sim_set_background_point_charge) uses a 3D-1/r SHAPE for
 * pedagogical familiarity (Stages 26/27); its sign convention matches
 * standard Maxwell now that the PIC source coefficient is also
 * standard.  Analytic and PIC paths are sign-consistent (both
 * "opposite attracts"); they differ only in spatial shape (1/r vs
 * 2D-log), so mixing them is meaningful, not contradictory.
 *
 * Pedagogical purpose: the classical-atom-collapse demo.  Opposite
 * charges attract via the 2D-log perturbation; the inductive back-
 * reaction (-q d_t A) slowly spirals the orbit in -- a self-consistent
 * 2D radiation-reaction inspiral.  Stage 27 already showed this for a
 * single test particle around a softened-charge background; this
 * scenario promotes both bodies to active sources.
 *
 * Parameters (defaults applied when not provided):
 *   params[0]: mass m of each particle      (default 0.01)
 *   params[1]: orbital radius r from COM    (default 10.0 * dx)
 *   params[2]: v_factor                     (default 1.0; 0 = free-fall)
 *   params[3]: cx of COM                    (default (W-1)/2 * dx)
 *   params[4]: cy of COM                    (default (H-1)/2 * dx)
 *   params[5]: charge magnitude Q           (default 0.01)
 *
 * Particle velocities are v_factor * v_orb with v_orb = Q * sqrt(k_e / m)
 * from 2D-log Coulomb centripetal balance:
 *
 *   Each particle's charge q deposits source rho_q = q delta(x-x_p).  The
 *   2D wave equation's static Green's function (-1/(2 pi)) ln(r) gives the
 *   perturbation phi_em^{pert}(r) sourced by a unit charge.  With source
 *   coupling -4 pi k_e, a charge q produces
 *       phi_em^{pert}(r) = 2 k_e q ln(r) + const.
 *   So the EM force on a second charge q' at separation d from charge q is
 *       F = -q' grad phi = -2 k_e q q' / d   (radial)
 *       sign: attractive if q q' < 0.
 *
 *   For opposite charges (q1 = +Q, q2 = -Q) at d = 2 r:
 *       |F| = 2 k_e Q^2 / (2 r) = k_e Q^2 / r          (attractive)
 *   Centripetal balance for each particle at radius r from COM:
 *       m v^2 / r = k_e Q^2 / r   =>   v = Q sqrt(k_e / m)
 *   (Note v_orb is INDEPENDENT of r, the 2D-log signature -- same form
 *   as the gravity binary's v_orb = sqrt(G_eff m), with G_eff m -> k_e Q^2.)
 *   Period (if bound): T = 2 pi r / v.
 *
 * Stability caveats inherited from pic_binary (gravity): PIC heating
 * at sub-cell scales injects kinetic energy near closest approach; v34
 * cell-centered unbound even short orbits, the v35 Yee+Esirkepov path
 * holds for the early-orbit phase.  Plus, with the inductive piece ON,
 * there is a real radiation-reaction inspiral on top of any residual
 * numerical heating -- Stage 32 measures the early-phase behavior to
 * isolate the validity regime. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <string.h>

static int build_pic_binary_em(gr_sim_t* sim, const float* params, int n_params) {
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    const float mass     = (n_params >= 1 && params[0] > 0.0f) ? params[0] : 0.01f;
    const float r_orb    = (n_params >= 2 && params[1] > 0.0f) ? params[1] : 10.0f * dx;
    const float v_factor = (n_params >= 3) ? params[2] : 1.0f;
    const float cx       = (n_params >= 4) ? params[3] : ((float) (W - 1) * 0.5f) * dx;
    const float cy       = (n_params >= 5) ? params[4] : ((float) (H - 1) * 0.5f) * dx;
    const float Q        = (n_params >= 6 && params[5] > 0.0f) ? params[5] : 0.01f;

    const float k_e   = gr_sim_get_k_e(sim);
    const float v_orb = v_factor * Q * sqrtf(k_e / mass);

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

    /* Vacuum perturbation-only EM setup. */
    gr_sim_set_bg_mode(sim, GR_BG_MODE_SAMPLED);
    gr_sim_set_field_evolution(sim, 1);
    gr_sim_set_particle_source_deposition(sim, 1);
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    /* EM Lorentz force with inductive piece ON (production default after
     * the half-step-A fix landed in 9769217).  Gravity force is irrelevant
     * here -- the particles have m for inertia only; no rho_matter sourcing
     * matters because we're isolating EM physics. */
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);

    /* Two OPPOSITE-SIGN-charge particles counter-orbiting around (cx, cy):
     *   p0 (charge +Q) at (cx - r, cy), velocity (0, +v)
     *   p1 (charge -Q) at (cx + r, cy), velocity (0, -v)
     * gives an attractive mutual orbit -- standard Maxwell sign
     * convention (opposite attracts).  Velocities are tangential to the
     * orbital plane; v_orb's magnitude is the same as for the
     * gravity binary's centripetal balance (formula above). */
    gr_sim_add_particle(sim, cx - r_orb, cy, mass, /*charge=*/+Q,
                        /*vx=*/0.0f, /*vy=*/+v_orb);
    gr_sim_add_particle(sim, cx + r_orb, cy, mass, /*charge=*/-Q,
                        /*vx=*/0.0f, /*vy=*/-v_orb);

    sim->step_count = 0;
    return 0;
}

static const gr_scenario_t SCENARIO_PIC_BINARY_EM = {
    .name  = "pic_binary_em",
    .build = build_pic_binary_em,
};

void gr_scenario_register_pic_binary_em(void) {
    gr_scenario_register(&SCENARIO_PIC_BINARY_EM);
}
