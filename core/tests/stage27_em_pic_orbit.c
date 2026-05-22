/* Stage 27 — EM analog of Stage 10/18b: charged particle orbiting a
 * central charge, with the EM perturbation FDTD active.
 *
 * This is the EM closed-loop PIC test.  A test charge orbits a softened
 * analytic POINT_CHARGE background.  The PERTURBATION EM fields
 * (phi_em, A_x, A_y) evolve via the wave equation under sources
 * deposited by the orbiting charge itself.  The test charge then
 * experiences the full EM Lorentz force from the perturbation,
 *
 *     F_em = q ( -grad phi_em - d_t A + v x B )
 *
 * including the inductive piece -q d_t A.  This is the EM analog of the
 * Phase-E gravity problem (memory [[grlite-lewis-birdsall-result]]):
 * does the moving-particle self-deposit feed back into a stable orbit,
 * or does numerical PIC self-heating (or its EM-side analog) unbind it?
 *
 * Physical expectation: a classical orbiting charge radiates EM waves
 * (Larmor formula), loses energy, and spirals inward.  This is the
 * 19th-century "atom collapsing" picture that quantum mechanics famously
 * fixed.  Our 2D-softened toy should reproduce qualitatively:
 *     P_rad ~ q^2 a^2 / c^3   (Larmor),   tau_decay ~ E_orbit / P_rad.
 *
 * Setup mirrors gravity's pic_orbiting at m_test = 1e-3 (the
 * weak-coupling regime where Phase E is stable post-LB-fix):
 *
 *   Q_central = +1, eps = 1, r = 20
 *   q_test    = -1e-3 (opposite sign for attractive orbit)
 *   m_test    = +1e-3 (same magnitude -> g_eff = |q*Q|/m * k_e = 1, same
 *                       orbital dynamics as gravity at m=1e-3)
 *   c = 1, k_e = 1, G_eff = 1
 *
 *   field_evolution ON, particle_source_deposition ON
 *   damping ring 16 cells (absorb outgoing radiation, suppress wraparound)
 *   ANALYTIC bg for the central charge (no sampled-bg discretization)
 *   CIC shape + LEGACY force interp (TSC and LB not yet ported to EM
 *      gradient path; this is the EM-CIC-legacy baseline analogous to
 *      pre-Stage-18 gravity)
 *
 * What we report:
 *   - radial drift over 4 orbital periods (Phase-E-style observable)
 *   - is the orbit (a) stable, (b) inspiral (radiation reaction wins),
 *     or (c) unbound (PIC heating wins)
 *   - sign and magnitude of the drift identifies which mechanism dominates */

#include "grlite.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, fmt, ...) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL %s:%d: " fmt "\n", __FILE__, __LINE__, ##__VA_ARGS__); \
        exit(1); \
    } \
} while (0)

typedef struct {
    float r_at_orbit[5];   /* r at orbits 0, 1, 2, 3, 4 (probed at theta=pi crossings) */
    int   n_completed;
    int   nan;
} result_t;

static void run(float Q, float q_test, float m_test, float r_orb,
                int n_orbits, int field_evolution, int inductive_enabled,
                gr_inductive_disc_t disc, float inductive_sign,
                result_t* out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float k_e    = 1.0f;

    /* Relativistic circular velocity: gamma v^2 = g r with
     * g = |q Q| k_e r / (r^2 + eps^2)^{3/2} / m_test. */
    const float g_mag  = fabsf(q_test * Q) * k_e * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f) / m_test;
    const float rg     = r_orb * g_mag;
    const float rg2_c2 = rg * rg / (c_eff * c_eff);
    const float u_v2   = (sqrtf(rg2_c2 * rg2_c2 + 4.0f * rg * rg) - rg2_c2) * 0.5f;
    const float v_circ = sqrtf(u_v2);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_circ;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { out->nan = 1; return; }
    gr_sim_set_damping(sim, 16);
    gr_sim_set_force_tier(sim, GR_FORCE_NEWTONIAN);
    gr_sim_set_field_evolution(sim, field_evolution);
    gr_sim_set_particle_source_deposition(sim, field_evolution);
    gr_sim_set_em_inductive_enabled(sim, inductive_enabled);
    gr_sim_set_em_inductive_disc(sim, disc);
    gr_sim_set_em_inductive_sign(sim, inductive_sign);
    gr_sim_set_background_point_charge(sim, cx, cy, Q, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    /* CIC + LEGACY: EM gradient path doesn't yet have TSC/LB variants.
     * rho_smooth_passes = 4 matches gravity's pic_orbiting production
     * defaults — empirically reduces moving-particle deposit aliasing. */
    gr_sim_set_shape_function(sim, GR_SHAPE_CIC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEGACY);
    gr_sim_set_rho_smooth_passes(sim, 4);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, q_test,
                        /*vx=*/0.0f, /*vy=*/v_circ);

    /* Probe r at each pi-wrap (i.e., once per orbit, at the diametrically-
     * opposite side of the orbit from the start).  The radius at the
     * starting side may differ from r at the opposite side by O(orbit
     * eccentricity); both are observable. */
    out->n_completed = 0;
    for (int k = 0; k <= n_orbits; k++) out->r_at_orbit[k] = 0.0f;

    /* Initial r at start (orbit "0" — t=0 reading). */
    {
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        out->r_at_orbit[0] = sqrtf(rx * rx + ry * ry);
    }

    const float dt = gr_sim_dt(sim);
    const int   n_max = (int) (1.2f * (float) n_orbits * T_ana / dt);
    float th_prev = 0.0f;
    int   wraps = 0;
    for (int s = 0; s < n_max && wraps < n_orbits; s++) {
        gr_sim_step(sim);
        const gr_particle_t* p = gr_sim_get_particle(sim, 0);
        const float rx = p->x - cx;
        const float ry = p->y - cy;
        const float r_now = sqrtf(rx * rx + ry * ry);
        if (!isfinite(r_now)) { out->nan = 1; gr_sim_destroy(sim); return; }
        const float th = atan2f(ry, rx);
        if (th_prev > 0.9f * (float) M_PI && th < -0.9f * (float) M_PI) {
            wraps++;
            if (wraps <= n_orbits) {
                out->r_at_orbit[wraps] = r_now;
                out->n_completed = wraps;
            }
        }
        th_prev = th;
    }
    out->nan = 0;
    gr_sim_destroy(sim);
}

int main(void) {
    printf("=== stage27_em_pic_orbit ===\n");
    printf("Closed-loop EM PIC: charged orbit + self-consistent EM perturbation.\n");
    printf("Spec: gr_sandbox_v35.tex sec:alg_rel; EM analog of Stage 10/18b.\n\n");

    const float Q       = +1.0f;
    const float q_test  = -1.0e-3f;
    const float m_test  = +1.0e-3f;
    const float r_orb   = 20.0f;
    const int   N       = 4;

    printf("Parameters: Q=%+g, q_test=%+g, m_test=%+g, r=%g\n",
           (double) Q, (double) q_test, (double) m_test, (double) r_orb);
    printf("            Coupling g_eff = |q*Q|*k_e/m = %g (matches gravity at m=1e-3).\n\n",
           (double) (fabsf(q_test * Q) / m_test));

    /* Baseline: field evolution OFF.  Analytic-bg-only orbit. */
    result_t base;
    run(Q, q_test, m_test, r_orb, N, /*field_evolution=*/0, /*induct=*/1,
        GR_INDUCTIVE_CENTERED, +1.0f, &base);
    TEST_ASSERT(!base.nan, "baseline went NaN");
    TEST_ASSERT(base.n_completed == N,
                "baseline didn't complete %d orbits (got %d)", N, base.n_completed);

    /* PIC variant A: all pieces, CENTERED inductive (default). */
    result_t pic_centered;
    run(Q, q_test, m_test, r_orb, N, /*field_evolution=*/1, /*induct=*/1,
        GR_INDUCTIVE_CENTERED, +1.0f, &pic_centered);

    /* PIC variant B: inductive disabled entirely. */
    result_t pic_noind;
    run(Q, q_test, m_test, r_orb, N, /*field_evolution=*/1, /*induct=*/0,
        GR_INDUCTIVE_CENTERED, +1.0f, &pic_noind);

    /* PIC variant C: inductive ENABLED, BACKWARD difference. */
    result_t pic_back;
    run(Q, q_test, m_test, r_orb, N, /*field_evolution=*/1, /*induct=*/1,
        GR_INDUCTIVE_BACKWARD, +1.0f, &pic_back);

    /* PIC variant D: inductive ENABLED, CENTERED diff, SIGN FLIPPED (-1).
     * Diagnostic: if the heating is purely a sign-convention error, this
     * should turn outward heating into inward inspiral. */
    result_t pic_flip;
    run(Q, q_test, m_test, r_orb, N, /*field_evolution=*/1, /*induct=*/1,
        GR_INDUCTIVE_CENTERED, -1.0f, &pic_flip);

    printf("Radius drift at each pi-wrap:\n");
    printf("  %-7s %-10s %-14s %-14s %-14s %-14s\n",
           "orbit", "baseline", "centered (+1)", "no inductive", "backward (+1)", "centered (-1)");
    for (int k = 1; k <= N; k++) {
        const float bd = (base.r_at_orbit[k]          - r_orb) / r_orb;
        const float cd = (pic_centered.r_at_orbit[k]  - r_orb) / r_orb;
        const float nd = (pic_noind.r_at_orbit[k]     - r_orb) / r_orb;
        const float xd = (pic_back.r_at_orbit[k]      - r_orb) / r_orb;
        const float fd = (pic_flip.r_at_orbit[k]      - r_orb) / r_orb;
        char ccell[32], ncell[32], xcell[32], fcell[32];
        if (k > pic_centered.n_completed) snprintf(ccell, sizeof(ccell), "<unbound>"); else snprintf(ccell, sizeof(ccell), "%+7.3f%%", 100.0 * (double) cd);
        if (k > pic_noind.n_completed)    snprintf(ncell, sizeof(ncell), "<unbound>"); else snprintf(ncell, sizeof(ncell), "%+7.3f%%", 100.0 * (double) nd);
        if (k > pic_back.n_completed)     snprintf(xcell, sizeof(xcell), "<unbound>"); else snprintf(xcell, sizeof(xcell), "%+7.3f%%", 100.0 * (double) xd);
        if (k > pic_flip.n_completed)     snprintf(fcell, sizeof(fcell), "<unbound>"); else snprintf(fcell, sizeof(fcell), "%+7.3f%%", 100.0 * (double) fd);
        printf("  %-7d %+7.3f%%   %-14s %-14s %-14s %-14s\n",
               k, 100.0 * (double) bd, ccell, ncell, xcell, fcell);
    }
    printf("\n");

    /* Verdict. */
    const float drift_cen   = (pic_centered.n_completed == N) ? (pic_centered.r_at_orbit[N] - r_orb) / r_orb : NAN;
    const float drift_noind = (pic_noind.n_completed   == N) ? (pic_noind.r_at_orbit[N]    - r_orb) / r_orb : NAN;
    const float drift_back  = (pic_back.n_completed    == N) ? (pic_back.r_at_orbit[N]     - r_orb) / r_orb : NAN;
    const float drift_flip  = (pic_flip.n_completed    == N) ? (pic_flip.r_at_orbit[N]     - r_orb) / r_orb : NAN;

    printf("After %d orbits:\n", N);
    if (isfinite(drift_cen))   printf("  Centered inductive +1 (default):       drift = %+.3f%%\n", 100.0 * (double) drift_cen);
    else                       printf("  Centered inductive +1: UNBOUND\n");
    if (isfinite(drift_noind)) printf("  Inductive disabled:                    drift = %+.3f%%\n", 100.0 * (double) drift_noind);
    else                       printf("  Inductive disabled: UNBOUND\n");
    if (isfinite(drift_back))  printf("  Backward-diff inductive +1:            drift = %+.3f%%\n", 100.0 * (double) drift_back);
    else                       printf("  Backward-diff inductive: UNBOUND\n");
    if (isfinite(drift_flip))  printf("  Centered inductive -1 (sign flipped):  drift = %+.3f%%\n", 100.0 * (double) drift_flip);
    else                       printf("  Centered inductive -1: UNBOUND\n");

    /* Soft assertion. */
    int stable_count = 0;
    if (isfinite(drift_cen)   && fabsf(drift_cen)   < 0.5f) stable_count++;
    if (isfinite(drift_noind) && fabsf(drift_noind) < 0.5f) stable_count++;
    if (isfinite(drift_back)  && fabsf(drift_back)  < 0.5f) stable_count++;
    if (isfinite(drift_flip)  && fabsf(drift_flip)  < 0.5f) stable_count++;
    TEST_ASSERT(stable_count >= 1, "All PIC variants unbound");

    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
