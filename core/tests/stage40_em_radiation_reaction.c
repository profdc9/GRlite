/* Stage 40 -- BENCHMARK B-6: Radiation reaction sanity check (gr_sandbox_v36 sec gempic_benchmark).
 *
 * EXPLORATORY ONLY -- NOT an acceptance gate.  See v36 sec
 * gempic_benchmark for the reasoning: in 2D retarded fields decay as
 * 1/r, not 1/r^2, so the textbook 3D Larmor formula does not apply.
 * The correct analytic benchmark is the 2D Lienard-Wiechert potential
 * which is not yet implemented.  This test produces a QUALITATIVE
 * report (orbit inspirals, loss rates exist, empirical v^k exponent)
 * not a hard threshold.
 *
 * Setup: charged particle orbits a softened point-charge background
 * (same setup as Stage 27).  Vary v_orbital by varying m_test at
 * fixed Q, q_test, r_orb.  Run N orbits, sample E and L per orbit.
 * Linear-fit the loss rates.
 *
 * Why this probe (per v36): an "accelerating charge MUST radiate"
 * check.  A scheme that conserves orbital KE exactly would be wrong
 * (back-reaction-free).  A working GEMPIC should produce monotone
 * inspiral with a non-zero loss rate.  Quantitative matching to a
 * 2D analytic is a separate follow-up scope item. */

#define _USE_MATH_DEFINES
#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float v_orb_target;      /* requested v_circ */
    float v_orb_actual;      /* actual v at start (relativistic) */
    int   n_orbits_done;
    int   nan;
    float E_at_orbit[8];     /* kinetic+potential energy at each orbit completion */
    float L_at_orbit[8];     /* z-component of angular momentum at each orbit */
    float r_at_orbit[8];     /* radius at each orbit */
    /* Linear-fit slopes: per-orbit rates. */
    float dE_per_orbit;       /* (E[N] - E[0]) / N */
    float dL_per_orbit;
    float dr_per_orbit;
    /* Same but with field evolution OFF (control: no radiation reaction). */
    float E_at_orbit_noff[8];
    float L_at_orbit_noff[8];
    float r_at_orbit_noff[8];
    float dE_per_orbit_noff;
    float dL_per_orbit_noff;
    float dr_per_orbit_noff;
} radiation_result_t;

static int run_orbit(float Q, float q_test, float m_test, float r_orb,
                     int n_orbits, int field_evolution,
                     float* E_arr, float* L_arr, float* r_arr,
                     float* v_actual_out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float k_e    = 1.0f;

    const float g_mag  = fabsf(q_test * Q) * k_e * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f) / m_test;
    const float rg     = r_orb * g_mag;
    const float rg2_c2 = rg * rg / (c_eff * c_eff);
    const float u_v2   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ = sqrtf(u_v2);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_circ;
    *v_actual_out = v_circ;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) return -1;
    gr_sim_set_damping(sim, 16);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, field_evolution);
    gr_sim_set_particle_source_deposition(sim, field_evolution);
    gr_sim_set_em_lorentz_force_enabled(sim, 1);
    gr_sim_set_em_electrostatic_enabled(sim, 1);
    gr_sim_set_em_inductive_enabled(sim, 1);
    gr_sim_set_em_magnetic_enabled(sim, 1);
    gr_sim_set_em_inductive_sign(sim, +1.0f);
    gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    /* TSC + LB everywhere (current production for EM).  rho_smooth = 0
     * per v36 sec gempic_benchmark "design principles" (bare deposit). */
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 0);
    gr_sim_set_j_smooth_passes(sim, 0);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test, 0.0f, v_circ);

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.2f * (float) n_orbits * T_ana / dt);

    /* Sample E, L, r at orbit 0 (start). */
    {
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r  = sqrtf(rx * rx + ry * ry);
        const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m_test * m_test * c_eff * c_eff));
        const float KE = (gamma - 1.0f) * m_test * c_eff * c_eff;
        /* 2D Coulomb potential: U = q*Q*k_e * ln(sqrt(r^2 + eps^2)).
         * For our purposes use the analytic background. */
        const float U  = q_test * Q * k_e * 0.5f * logf(r * r + eps * eps);
        const float vx = p->px / (gamma * m_test);
        const float vy = p->py / (gamma * m_test);
        const float Lz = m_test * (rx * vy - ry * vx);
        E_arr[0] = KE + U;
        L_arr[0] = Lz;
        r_arr[0] = r;
    }

    int wraps = 0;
    float th_prev = 0.0f;
    for (int s = 0; s < n_max && wraps < n_orbits; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) { gr_sim_destroy(sim); return -1; }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            if (wraps < 8) {
                const float gamma = sqrtf(1.0f + (p->px * p->px + p->py * p->py) / (m_test * m_test * c_eff * c_eff));
                const float KE = (gamma - 1.0f) * m_test * c_eff * c_eff;
                const float U  = q_test * Q * k_e * 0.5f * logf(r_now * r_now + eps * eps);
                const float vx = p->px / (gamma * m_test);
                const float vy = p->py / (gamma * m_test);
                E_arr[wraps] = KE + U;
                L_arr[wraps] = m_test * (rx * vy - ry * vx);
                r_arr[wraps] = r_now;
            }
        }
        th_prev = th;
    }
    (void) dt;
    gr_sim_destroy(sim);
    return wraps;
}

static radiation_result_t measure_radiation(float Q, float q_test, float m_test,
                                              float r_orb, int n_orbits) {
    radiation_result_t r = {0};
    r.v_orb_target = sqrtf(fabsf(q_test * Q) * r_orb
                           / powf(r_orb * r_orb + 1.0f, 1.5f) / m_test);
    float v_actual;
    int wraps_field = run_orbit(Q, q_test, m_test, r_orb, n_orbits, /*field_evol=*/1,
                                 r.E_at_orbit, r.L_at_orbit, r.r_at_orbit, &v_actual);
    r.v_orb_actual = v_actual;
    if (wraps_field < 0) { r.nan = 1; return r; }
    r.n_orbits_done = wraps_field;
    int wraps_noff = run_orbit(Q, q_test, m_test, r_orb, n_orbits, /*field_evol=*/0,
                                r.E_at_orbit_noff, r.L_at_orbit_noff,
                                r.r_at_orbit_noff, &v_actual);
    if (wraps_noff < 0) { r.nan = 1; return r; }

    /* Simple slope: (last - first) / n. */
    if (wraps_field > 0) {
        r.dE_per_orbit = (r.E_at_orbit[wraps_field] - r.E_at_orbit[0]) / (float) wraps_field;
        r.dL_per_orbit = (r.L_at_orbit[wraps_field] - r.L_at_orbit[0]) / (float) wraps_field;
        r.dr_per_orbit = (r.r_at_orbit[wraps_field] - r.r_at_orbit[0]) / (float) wraps_field;
    }
    if (wraps_noff > 0) {
        r.dE_per_orbit_noff = (r.E_at_orbit_noff[wraps_noff] - r.E_at_orbit_noff[0]) / (float) wraps_noff;
        r.dL_per_orbit_noff = (r.L_at_orbit_noff[wraps_noff] - r.L_at_orbit_noff[0]) / (float) wraps_noff;
        r.dr_per_orbit_noff = (r.r_at_orbit_noff[wraps_noff] - r.r_at_orbit_noff[0]) / (float) wraps_noff;
    }
    return r;
}

int main(void) {
    printf("=== stage40_em_radiation_reaction (BENCHMARK B-6, EXPLORATORY) ===\n");
    printf("Charged particle in PIC orbit around an analytic point-charge background.\n");
    printf("Report E, L, r per orbit and the per-orbit loss rates.  Field-evolution-\n");
    printf("OFF rows are the no-radiation control.  Field-evolution-ON rows include\n");
    printf("PIC self-field and should show monotone inspiral if the scheme is\n");
    printf("producing back-reaction.\n\n");
    printf("NOT a hard gate.  See v36 sec gempic_benchmark on why 2D physics\n");
    printf("makes quantitative matching to textbook Larmor inappropriate; this\n");
    printf("is a qualitative report.\n\n");

    const float Q       = +1.0f;
    const float q_test  = -1.0e-3f;
    const float r_orb   = 20.0f;
    const int   N_orbits = 4;

    /* Vary m_test to get different v_circ. */
    const float m_vals[]   = { 4.0e-3f, 1.0e-3f, 2.5e-4f };
    const int   nm         = (int)(sizeof(m_vals) / sizeof(m_vals[0]));

    for (int i = 0; i < nm; i++) {
        const float m_test = m_vals[i];
        radiation_result_t rr = measure_radiation(Q, q_test, m_test, r_orb, N_orbits);
        const float g_eff = fabsf(q_test * Q) / m_test;
        printf("--- m_test=%g  g_eff=%g  v_circ=%.4f  (%d orbits completed of %d) ---\n",
               (double) m_test, (double) g_eff, (double) rr.v_orb_actual,
               rr.n_orbits_done, N_orbits);
        if (rr.nan) {
            printf("  NaN -- run diverged or particle left domain.\n\n");
            continue;
        }
        printf("  %-6s %-13s %-13s %-13s %-13s %-13s %-13s\n",
               "orbit", "E(field)", "L(field)", "r(field)", "E(noff)", "L(noff)", "r(noff)");
        for (int k = 0; k <= rr.n_orbits_done; k++) {
            printf("  %-6d %+12.4e %+12.4e %12.4e %+12.4e %+12.4e %12.4e\n",
                   k,
                   (double) rr.E_at_orbit[k], (double) rr.L_at_orbit[k], (double) rr.r_at_orbit[k],
                   (double) rr.E_at_orbit_noff[k], (double) rr.L_at_orbit_noff[k], (double) rr.r_at_orbit_noff[k]);
        }
        printf("  dE/orbit (field) = %+.4e   dL/orbit (field) = %+.4e   dr/orbit (field) = %+.4e\n",
               (double) rr.dE_per_orbit, (double) rr.dL_per_orbit, (double) rr.dr_per_orbit);
        printf("  dE/orbit (noff)  = %+.4e   dL/orbit (noff)  = %+.4e   dr/orbit (noff)  = %+.4e\n",
               (double) rr.dE_per_orbit_noff, (double) rr.dL_per_orbit_noff, (double) rr.dr_per_orbit_noff);
        printf("\n");
    }

    printf("Interpretation:\n");
    printf("  field-evolution-OFF rows should show near-zero per-orbit drift\n");
    printf("    (this is the no-radiation control; remaining drift is\n");
    printf("    Boris-pusher truncation error).\n");
    printf("  field-evolution-ON rows show the PIC dynamics.  In a working\n");
    printf("    scheme, dE/orbit and dL/orbit should be NEGATIVE (energy\n");
    printf("    radiating away, monotone inspiral, dr/orbit < 0).\n");
    printf("  The current scheme is known to give the WRONG sign at the\n");
    printf("    operator level (Stage 37 diagnosis); orbit inspiral vs\n");
    printf("    outspiral here is downstream of that and not interpretable\n");
    printf("    as physical radiation reaction.  Numbers recorded as the\n");
    printf("    pre-GEMPIC baseline.\n");

    printf("\nBASELINE COMPLETE.\n");
    return 0;
}
