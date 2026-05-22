/* Stage 28 — does gravity exhibit the same inductive PIC heating as
 * EM (Stage 27) when the -m d_t A_g piece is enabled?
 *
 * Diagnostic to test the user's hypothesis: gravity at m=1e-3 appears
 * stable (Stage 18b: -2.4%/orbit inward, consistent with radiation
 * reaction), but only because the gravity implementation DOES NOT
 * include the -m d_t A_g term that the doc's Tier-3 force
 * (gr_sandbox_v35.tex eq:eih_full / sec:alg_rel eqbox line 1040)
 * actually prescribes.  If we wire that piece in symmetrically with
 * EM's -q d_t A piece, does gravity also heat?
 *
 * If YES: the inductive heating is a STRUCTURAL discretization issue
 * with the simple centered-time-difference (not EM-specific), and the
 * gravity result was an accidental stability win from omitting the
 * piece.  Both sectors would need an implicit/symplectic time
 * integration fix.
 *
 * If NO: there is something specific to EM that drives the heating
 * (perhaps the coefficient +1 on v x B vs gravity's +4 spin-2
 * enhancement, or the source-coefficient sign convention), worth
 * deeper investigation. */

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
    float r_at_orbit[5];
    int   n_completed;
    int   nan;
} result_t;

static void run(float GM, float m_test, float r_orb, int n_orbits,
                int field_evolution, int grav_inductive_enabled,
                float inductive_sign, int j_time_correction,
                result_t* out) {
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float eps    = 1.0f;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;

    const float g_mag  = GM * r_orb
                       / powf(r_orb * r_orb + eps * eps, 1.5f);
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
    gr_sim_set_background_point_mass(sim, cx, cy, GM, eps);
    gr_sim_set_bg_mode(sim, GR_BG_MODE_ANALYTIC);
    /* Match Stage 18b/pic_orbiting production settings. */
    gr_sim_set_shape_function(sim, GR_SHAPE_TSC);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, 4);
    /* The diagnostic switches. */
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, grav_inductive_enabled);
    gr_sim_set_gravitomagnetic_inductive_sign(sim, inductive_sign);
    gr_sim_set_j_time_correction_enabled(sim, j_time_correction);

    gr_sim_add_particle(sim, cx + r_orb, cy, m_test, /*charge=*/0.0f,
                        /*vx=*/0.0f, /*vy=*/v_circ);

    out->n_completed = 0;
    for (int k = 0; k <= n_orbits; k++) out->r_at_orbit[k] = 0.0f;
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
    printf("=== stage28_grav_inductive_heating ===\n");
    printf("Same setup as Stage 18b (gravity PIC orbit at m_test = 1e-3),\n");
    printf("but with the -m d_t A_g inductive piece toggled on/off.\n\n");

    const float GM     = 1.0f;
    const float m_test = 1.0e-3f;
    const float r_orb  = 20.0f;
    const int   N      = 4;

    /* Reference: field_evolution OFF (analytic bg only). */
    result_t base;
    run(GM, m_test, r_orb, N, /*field_ev=*/0, /*grav_ind=*/0, /*sign=*/+1.0f, /*jtc=*/0, &base);
    TEST_ASSERT(!base.nan, "baseline NaN");

    /* Stage 18b current default: PIC with NO inductive. */
    result_t pic_noind;
    run(GM, m_test, r_orb, N, /*field_ev=*/1, /*grav_ind=*/0, /*sign=*/+1.0f, /*jtc=*/0, &pic_noind);

    /* PIC with inductive -m d_t A_g enabled, NORMAL (+1) sign. */
    result_t pic_ind;
    run(GM, m_test, r_orb, N, /*field_ev=*/1, /*grav_ind=*/1, /*sign=*/+1.0f, /*jtc=*/0, &pic_ind);

    /* PIC with inductive ENABLED, sign FLIPPED (-1). */
    result_t pic_ind_flip;
    run(GM, m_test, r_orb, N, /*field_ev=*/1, /*grav_ind=*/1, /*sign=*/-1.0f, /*jtc=*/0, &pic_ind_flip);

    /* PIC with inductive ENABLED + J time-correction.  Tests whether the
     * source-time-staggering is the main culprit. */
    result_t pic_ind_jtc;
    run(GM, m_test, r_orb, N, /*field_ev=*/1, /*grav_ind=*/1, /*sign=*/+1.0f, /*jtc=*/1, &pic_ind_jtc);

    printf("Parameters: GM=%g, m_test=%g, r=%g (matches Stage 18b regime).\n\n",
           (double) GM, (double) m_test, (double) r_orb);

    printf("Radius drift at each pi-wrap:\n");
    printf("  %-8s %-13s %-13s %-13s %-13s %-15s\n",
           "orbit", "baseline", "no induct", "ind +1", "ind -1", "ind +1 + J-tc");
    for (int k = 1; k <= N; k++) {
        const float bd = (base.r_at_orbit[k] - r_orb) / r_orb;
        char ncell[32], icell[32], fcell[32], jcell[32];
        if (k > pic_noind.n_completed)    snprintf(ncell, sizeof(ncell), "<unbound>");
        else snprintf(ncell, sizeof(ncell), "%+7.3f%%", 100.0 * (double) ((pic_noind.r_at_orbit[k] - r_orb) / r_orb));
        if (k > pic_ind.n_completed)      snprintf(icell, sizeof(icell), "<unbound>");
        else snprintf(icell, sizeof(icell), "%+7.3f%%", 100.0 * (double) ((pic_ind.r_at_orbit[k] - r_orb) / r_orb));
        if (k > pic_ind_flip.n_completed) snprintf(fcell, sizeof(fcell), "<unbound>");
        else snprintf(fcell, sizeof(fcell), "%+7.3f%%", 100.0 * (double) ((pic_ind_flip.r_at_orbit[k] - r_orb) / r_orb));
        if (k > pic_ind_jtc.n_completed)  snprintf(jcell, sizeof(jcell), "<unbound>");
        else snprintf(jcell, sizeof(jcell), "%+7.3f%%", 100.0 * (double) ((pic_ind_jtc.r_at_orbit[k] - r_orb) / r_orb));
        printf("  %-8d  %+7.3f%%      %-13s %-13s %-13s %-15s\n",
               k, 100.0 * (double) bd, ncell, icell, fcell, jcell);
    }
    printf("\n");

    const float drift_noind   = (pic_noind.n_completed    == N) ? (pic_noind.r_at_orbit[N]    - r_orb) / r_orb : NAN;
    const float drift_ind     = (pic_ind.n_completed      == N) ? (pic_ind.r_at_orbit[N]      - r_orb) / r_orb : NAN;
    const float drift_flip    = (pic_ind_flip.n_completed == N) ? (pic_ind_flip.r_at_orbit[N] - r_orb) / r_orb : NAN;

    printf("After %d orbits:\n", N);
    if (isfinite(drift_noind)) printf("  Gravity PIC, no inductive (Stage 18b baseline): drift = %+.3f%%\n", 100.0 * (double) drift_noind);
    else                       printf("  Gravity PIC, no inductive: UNBOUND\n");
    if (isfinite(drift_ind))   printf("  Gravity PIC, inductive +1 (normal sign):        drift = %+.3f%%\n", 100.0 * (double) drift_ind);
    else                       printf("  Gravity PIC, inductive +1: UNBOUND\n");
    if (isfinite(drift_flip))  printf("  Gravity PIC, inductive -1 (sign flipped):       drift = %+.3f%%\n", 100.0 * (double) drift_flip);
    else                       printf("  Gravity PIC, inductive -1: UNBOUND\n");
    const float drift_jtc = (pic_ind_jtc.n_completed == N) ? (pic_ind_jtc.r_at_orbit[N] - r_orb) / r_orb : NAN;
    if (isfinite(drift_jtc))  printf("  Gravity PIC, inductive +1 + J time-correction:  drift = %+.3f%%\n", 100.0 * (double) drift_jtc);
    else                      printf("  Gravity PIC, inductive +1 + J-tc: UNBOUND\n");
    printf("\n");
    printf("Comparison with EM Stage 27 (q=m=1e-3, same coupling):\n");
    printf("  EM PIC, no inductive:            drift = -7.0%%   (inspiral)\n");
    printf("  EM PIC, with -q d_t A centered:  drift = +18.6%%  (heating)\n");

    /* No assertion — this is purely diagnostic. */
    printf("\nALL CHECKS PASSED.\n");
    return 0;
}
