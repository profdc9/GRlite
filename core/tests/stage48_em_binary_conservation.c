/* Stage 48 -- equal-coupling EM binary CONSERVATION diagnostic.
 *
 * Motivation (Dan, 2026-05-26): The 2D-log Coulomb force scales as 1/d
 * (not 1/d^2 like 3D Newton/Kepler).  By Bertrand's theorem only 1/r^2
 * and r^2 central forces give CLOSED orbits in 3D; the 2D log potential
 * is neither.  So even with perfect numerics, any deviation from an
 * exactly circular initial condition produces a non-closed rosette with
 * radial oscillation between r_min and r_max while angle advances.
 *
 * Stage 47's "% drift in radius at each theta-wrap" metric conflates
 * physical rosette behavior with numerical heating.  Stage 48 instead
 * uses CONSERVATION LAWS, which must hold regardless of orbit shape:
 *
 *   E_total      = E_kin + U_int (analytic 2D Coulomb)
 *   L_z (total)  about COM
 *   P_x, P_y     (linear momentum, ~0 by symmetry)
 *
 * Radiation reaction can drain E and |L| slowly (~ Larmor-2D), but over
 * a few orbits at non-rel v this should be a tiny effect.  Anything
 * larger is numerical.  We ALSO log the max |phi_em| at the PML inner
 * boundary to detect reflection.
 *
 * Same configurations as Stage 47: TSC bare, TSC+N=16, BUMP R=6/R=8.
 * Also a non-relativistic v=0.01 sweep (10x slower) per Dan's third
 * check item. */

#define _USE_MATH_DEFINES
#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NMAX 16
typedef struct {
    /* Snapshots at fixed-time intervals (t_k = k * (run_total/NMAX-1)). */
    float t[NMAX];
    float E_kin[NMAX];        /* relativistic kinetic energy sum */
    float U_int[NMAX];        /* analytic 2D Coulomb interaction PE (attractive) */
    float L_z[NMAX];          /* total z-angular momentum about COM */
    float Px[NMAX], Py[NMAX];
    float sep[NMAX];          /* separation |r1 - r2| */
    /* Run-wide maxima. */
    float max_phi_em_pml_inner;     /* max |phi_em| at PML inner edge cells */
    float max_phi_em_bulk;          /* max |phi_em| anywhere outside PML */
    int   n_samples;
    int   nan;
    const char* tag;
} diag_t;

static void snapshot(diag_t* out, int idx, float t,
                      const gr_particle_t* p0, const gr_particle_t* p1,
                      float mass, float k_e, float Q, float cx, float cy) {
    const float g0 = sqrtf(1.0f + (p0->px*p0->px + p0->py*p0->py)/(mass*mass));
    const float g1 = sqrtf(1.0f + (p1->px*p1->px + p1->py*p1->py)/(mass*mass));
    const float dxs = p1->x - p0->x;
    const float dys = p1->y - p0->y;
    const float d   = sqrtf(dxs*dxs + dys*dys);
    const float Rx0 = p0->x - cx, Ry0 = p0->y - cy;
    const float Rx1 = p1->x - cx, Ry1 = p1->y - cy;
    out->t[idx]     = t;
    out->E_kin[idx] = (g0 - 1.0f)*mass + (g1 - 1.0f)*mass;
    /* Attractive 2D Coulomb PE: U(d) = +2 k_e Q^2 ln(d) (rises with d). */
    out->U_int[idx] = +2.0f*k_e*Q*Q*logf(d);
    out->L_z[idx]   = (Rx0*p0->py - Ry0*p0->px)
                    + (Rx1*p1->py - Ry1*p1->px);
    out->Px[idx]    = p0->px + p1->px;
    out->Py[idx]    = p0->py + p1->py;
    out->sep[idx]   = d;
}

static void run_diagnostic(int n_orbits, float v_factor, float mass_in,
                            gr_shape_function_t shape, float Rk, int smooth,
                            diag_t* out) {
    out->nan = 0;     /* prior version left this uninitialized */
    const int   W      = 256, H = 256;
    const float dx     = 1.0f;
    const float c_eff  = 1.0f;
    const float cfl    = 1.0f / sqrtf(2.0f);
    const float mass   = mass_in;
    const float Q      = 0.01f;
    const float r_orb  = 20.0f * dx;
    const float cx     = ((float) (W - 1) * 0.5f) * dx;
    const float cy     = ((float) (H - 1) * 0.5f) * dx;
    const float k_e    = 1.0f;
    const float v_orb  = v_factor * Q * sqrtf(k_e / mass);
    const float T_ana  = 2.0f * (float) M_PI * r_orb / v_orb;

    gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
    if (!sim) { out->nan = 1; return; }
    gr_sim_use_pedagogical_defaults(sim);
    gr_sim_set_damping(sim, 16);

    const float params[6] = {mass, r_orb, v_factor, cx, cy, Q};
    if (gr_sim_load_scenario(sim, "pic_binary_em", params, 6) != 0) {
        out->nan = 1; gr_sim_destroy(sim); return;
    }
    /* Override comparison config AFTER scenario load. */
    gr_sim_set_shape_function(sim, shape);
    if (shape == GR_SHAPE_BUMP) gr_sim_set_kernel_radius(sim, Rk);
    gr_sim_set_force_interp(sim, GR_FORCE_INTERP_LEWIS_BIRDSALL);
    gr_sim_set_rho_smooth_passes(sim, smooth);
    gr_sim_set_j_smooth_passes(sim, smooth);
    /* PURE EM (gravity disabled, per Stage 47 correction). */
    gr_sim_set_G_eff(sim, 0.0f);
    gr_sim_set_gravitomagnetic_force_enabled(sim, 0);
    gr_sim_set_gravitomagnetic_inductive_enabled(sim, 0);

    const float dt = gr_sim_dt(sim);
    const int   n_steps = (int) ((float) n_orbits * T_ana / dt);
    /* Take NMAX evenly-spaced samples in time over the full n_orbits run.
     * Time-uniform sampling decouples the diagnostic from orbit-counting
     * (which is ambiguous when the orbit is a precessing rosette). */
    const int sample_every = n_steps / (NMAX - 1);
    out->n_samples = 0;
    out->max_phi_em_pml_inner = 0.0f;
    out->max_phi_em_bulk = 0.0f;
    const int pml = 16;

    /* Initial sample at t=0. */
    {
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        snapshot(out, 0, 0.0f, p0, p1, mass, k_e, Q, cx, cy);
        out->n_samples = 1;
    }

    float t = 0.0f;
    for (int s = 0; s < n_steps; s++) {
        gr_sim_step(sim);
        t += dt;
        const gr_particle_t* p0 = gr_sim_get_particle(sim, 0);
        const gr_particle_t* p1 = gr_sim_get_particle(sim, 1);
        if (!isfinite(p0->x) || !isfinite(p1->x)) { out->nan = 1; break; }

        /* Time-uniform snapshot. */
        if ((s+1) % sample_every == 0 && out->n_samples < NMAX) {
            snapshot(out, out->n_samples, t, p0, p1, mass, k_e, Q, cx, cy);
            out->n_samples++;
        }

        /* PML inner-edge probe + bulk amplitude. */
        if ((s % 64) == 0) {
            const float* phi_em = sim->fields[GR_FIELD_PHI_EM].curr;
            float m_pml = 0.0f;
            for (int i = pml; i < W - pml; i++) {
                float v = fabsf(phi_em[pml*W + i]);
                if (v > m_pml) m_pml = v;
                v = fabsf(phi_em[(H-pml-1)*W + i]);
                if (v > m_pml) m_pml = v;
            }
            for (int j = pml + 1; j < H - pml - 1; j++) {
                float v = fabsf(phi_em[j*W + pml]);
                if (v > m_pml) m_pml = v;
                v = fabsf(phi_em[j*W + (W-pml-1)]);
                if (v > m_pml) m_pml = v;
            }
            if (m_pml > out->max_phi_em_pml_inner) out->max_phi_em_pml_inner = m_pml;
            /* Bulk peak: max within interior square 2*pml in from edges. */
            float m_bulk = 0.0f;
            for (int j = 2*pml; j < H - 2*pml; j += 4) {
                for (int i = 2*pml; i < W - 2*pml; i += 4) {
                    const float v = fabsf(phi_em[j*W + i]);
                    if (v > m_bulk) m_bulk = v;
                }
            }
            if (m_bulk > out->max_phi_em_bulk) out->max_phi_em_bulk = m_bulk;
        }
    }
    out->nan |= 0;
    gr_sim_destroy(sim);
}

static void print_table(const diag_t* d) {
    if (d->n_samples < 2) {
        printf("  %-15s : NaN before any samples\n", d->tag); return;
    }
    printf("  %-15s (%d samples, NaN=%d):\n", d->tag, d->n_samples, d->nan);
    const float E0 = d->E_kin[0] + d->U_int[0];
    const float L0 = d->L_z[0];
    const float S0 = d->sep[0];
    printf("    sample t          d/d0       L/L0      E_K          U_int        E_tot       dE/|E0|\n");
    for (int k = 0; k < d->n_samples; k++) {
        const float Et = d->E_kin[k] + d->U_int[k];
        printf("    %2d   %8.1f  %+7.2f%%  %+7.2f%%  %+11.5e %+11.5e %+11.5e %+8.3f%%\n",
               k, (double)d->t[k],
               (double)(100.0f*(d->sep[k] - S0)/S0),
               (double)(100.0f*(d->L_z[k] - L0)/L0),
               (double)d->E_kin[k], (double)d->U_int[k], (double)Et,
               (double)(100.0f*(Et - E0)/fabsf(E0)));
    }
    printf("    max |phi_em| PML edge: %.3e   bulk: %.3e   ratio: %.2f%%\n",
           (double)d->max_phi_em_pml_inner, (double)d->max_phi_em_bulk,
           (double)(100.0f * d->max_phi_em_pml_inner / fmaxf(d->max_phi_em_bulk, 1e-30f)));
    printf("    Px range: [%+.3e .. %+.3e]   Py range: [%+.3e .. %+.3e]\n",
           (double)d->Px[0], (double)d->Px[d->n_samples-1],
           (double)d->Py[0], (double)d->Py[d->n_samples-1]);
}

int main(void) {
    printf("=== stage48_em_binary_conservation ===\n");
    printf("Equal-mass m=0.01, opposite charge q=+/-0.01.  r=20.\n");
    printf("Conservation diagnostic (E, L_z, P) -- decouples physical rosette\n");
    printf("precession (Bertrand's theorem on 1/r force) from numerical heating.\n\n");

    const int N = 4;

    printf("--- (A) m=0.01, v_factor=1: v_orb=0.1 (mildly relativistic, CIRCULAR) ---\n");
    {
        diag_t d_tsc = {0};    d_tsc.tag    = "TSC bare";
        diag_t d_tsc16 = {0};  d_tsc16.tag  = "TSC + N=16";
        diag_t d_r6 = {0};     d_r6.tag     = "BUMP R=6";
        diag_t d_r8 = {0};     d_r8.tag     = "BUMP R=8";
        run_diagnostic(N, 1.0f, 0.01f, GR_SHAPE_TSC,  0.0f, 0,  &d_tsc);
        run_diagnostic(N, 1.0f, 0.01f, GR_SHAPE_TSC,  0.0f, 16, &d_tsc16);
        run_diagnostic(N, 1.0f, 0.01f, GR_SHAPE_BUMP, 6.0f, 0,  &d_r6);
        run_diagnostic(N, 1.0f, 0.01f, GR_SHAPE_BUMP, 8.0f, 0,  &d_r8);
        print_table(&d_tsc);  print_table(&d_tsc16);
        print_table(&d_r6);   print_table(&d_r8);
    }
    printf("\n");

    printf("--- (B) m=0.01, v_factor=0.1: v_orb=0.01 (LOW-L, orbit dives to d~1.55) ---\n");
    printf("    Inner turning point d_in approximate 1.55 cells (below kernel width):\n");
    printf("    orbit physically collapses; all configs fail when d < kernel scale.\n");
    {
        diag_t d_tsc = {0};    d_tsc.tag    = "TSC bare";
        diag_t d_tsc16 = {0};  d_tsc16.tag  = "TSC + N=16";
        diag_t d_r6 = {0};     d_r6.tag     = "BUMP R=6";
        diag_t d_r8 = {0};     d_r8.tag     = "BUMP R=8";
        run_diagnostic(N, 0.1f, 0.01f, GR_SHAPE_TSC,  0.0f, 0,  &d_tsc);
        run_diagnostic(N, 0.1f, 0.01f, GR_SHAPE_TSC,  0.0f, 16, &d_tsc16);
        run_diagnostic(N, 0.1f, 0.01f, GR_SHAPE_BUMP, 6.0f, 0,  &d_r6);
        run_diagnostic(N, 0.1f, 0.01f, GR_SHAPE_BUMP, 8.0f, 0,  &d_r8);
        print_table(&d_tsc);  print_table(&d_tsc16);
        print_table(&d_r6);   print_table(&d_r8);
    }
    printf("\n");

    printf("--- (C) m=1.0,  v_factor=1: v_orb=0.01 (PROPER NON-REL CIRCULAR at d=40) ---\n");
    {
        diag_t d_tsc = {0};    d_tsc.tag    = "TSC bare";
        diag_t d_tsc16 = {0};  d_tsc16.tag  = "TSC + N=16";
        diag_t d_r6 = {0};     d_r6.tag     = "BUMP R=6";
        diag_t d_r8 = {0};     d_r8.tag     = "BUMP R=8";
        run_diagnostic(N, 1.0f, 1.0f, GR_SHAPE_TSC,  0.0f, 0,  &d_tsc);
        run_diagnostic(N, 1.0f, 1.0f, GR_SHAPE_TSC,  0.0f, 16, &d_tsc16);
        run_diagnostic(N, 1.0f, 1.0f, GR_SHAPE_BUMP, 6.0f, 0,  &d_r6);
        run_diagnostic(N, 1.0f, 1.0f, GR_SHAPE_BUMP, 8.0f, 0,  &d_r8);
        print_table(&d_tsc);  print_table(&d_tsc16);
        print_table(&d_r6);   print_table(&d_r8);
    }
    return 0;
}
