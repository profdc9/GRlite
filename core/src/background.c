/* Background field array management and generators.
 * Spec reference: gr_sandbox_vNN.tex §12.6 (sec:stage_bg)
 *                 "Stage 6: Sampled background field arrays" +
 *                 §sec:yee_pivot (v35 Yee sublattice layout).
 *
 * Per §9 each background array lives on its own Yee sublattice:
 *   phi_g_bg, phi_bg : GR_LATTICE_CORNER  - nodes at (i,   j  ) * dx
 *   Agx_bg,   Ax_bg  : GR_LATTICE_X_EDGE  - nodes at (i+0.5, j ) * dx
 *   Agy_bg,   Ay_bg  : GR_LATTICE_Y_EDGE  - nodes at (i,   j+0.5) * dx
 *
 * The installers below sample the analytic generator at each array's own
 * sublattice node positions.  The analytic-mode evaluators
 * (gr_bg_eval_analytic, gr_bg_eval_A_g) read the closed form directly at
 * an arbitrary (x, y), so they are sublattice-agnostic. */

#include "grlite.h"
#include "sim_internal.h"

#include <math.h>
#include <stdlib.h>

/* Lookup the storage slot for a given field id. Returns a pointer-to-pointer
 * so callers can both read and (re)assign the array. NULL for unknown ids. */
static float** background_slot(gr_sim_t* sim, gr_field_id_t which) {
    if (!sim) return NULL;
    switch (which) {
    case GR_FIELD_PHI_GRAV: return &sim->phi_g_bg;
    case GR_FIELD_A_GX:     return &sim->Agx_bg;
    case GR_FIELD_A_GY:     return &sim->Agy_bg;
    case GR_FIELD_PHI_EM:   return &sim->phi_bg;
    case GR_FIELD_A_X:      return &sim->Ax_bg;
    case GR_FIELD_A_Y:      return &sim->Ay_bg;
    default: return NULL;
    }
}

float* gr_sim_background_ptr(gr_sim_t* sim, gr_field_id_t which) {
    float** slot = background_slot(sim, which);
    return slot ? *slot : NULL;
}

/* Free the six sampled background arrays (but keep the body list / params). */
static void free_bg_arrays(gr_sim_t* sim) {
    free(sim->phi_g_bg); sim->phi_g_bg = NULL;
    free(sim->Agx_bg);   sim->Agx_bg   = NULL;
    free(sim->Agy_bg);   sim->Agy_bg   = NULL;
    free(sim->phi_bg);   sim->phi_bg   = NULL;
    free(sim->Ax_bg);    sim->Ax_bg    = NULL;
    free(sim->Ay_bg);    sim->Ay_bg    = NULL;
}

void gr_sim_clear_background(gr_sim_t* sim) {
    if (!sim) return;
    free_bg_arrays(sim);
    /* Drop the compact-body list too. */
    sim->bg_nbody = 0;
    /* Clear analytic generator parameters too — they're paired with the
     * sampled array.  Mode is left alone (user-set). */
    sim->bg_kind   = GR_BG_KIND_NONE;
    sim->bg_x0     = 0.0f;
    sim->bg_y0     = 0.0f;
    sim->bg_GM         = 0.0f;
    sim->bg_eps        = 0.0f;
    sim->bg_charge     = 0.0f;
    sim->bg_Jz         = 0.0f;
    sim->bg_B0         = 0.0f;
    sim->bg_B0_em      = 0.0f;
    sim->bg_B_prime_em = 0.0f;
    sim->bg_Ex_em      = 0.0f;
    sim->bg_Ey_em  = 0.0f;
    sim->bg_g_em   = 0.0f;
}

/* Lazily allocate and zero a background array if not already present. */
static float* ensure_bg_alloc(gr_sim_t* sim, gr_field_id_t which) {
    float** slot = background_slot(sim, which);
    if (!slot) return NULL;
    if (*slot) return *slot;
    const size_t n = (size_t) sim->width * (size_t) sim->height;
    *slot = (float*) calloc(n, sizeof(float));
    return *slot;
}

/* ---- Shared field generators (fill arrays; no clear, no param store) -------
 * Each writes one component's analytic field onto its Yee sublattice.  The
 * public setters below clear + call these + store the generator params; the
 * unified gr_sim_set_background_body clears ONCE and calls all three so the
 * components superpose. */

/* Phi_g^{bg}(x) = -GM / sqrt(r^2 + eps^2) on the CORNER sublattice. */
static void fill_bg_grav_scalar(gr_sim_t* sim, float x0, float y0,
                                float GM, float epsilon) {
    float* phi_bg = ensure_bg_alloc(sim, GR_FIELD_PHI_GRAV);
    if (!phi_bg) return;
    const int W = sim->width, H = sim->height;
    const float dx = sim->dx, eps2 = epsilon * epsilon;
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx, dy = y - y0;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xx = (float) i * dx, dxi = xx - x0;
            phi_bg[row + i] += -GM / sqrtf(dxi * dxi + dy * dy + eps2);
        }
    }
}

/* Gravitomagnetic dipole A_g = (G Jz / 2c^2) (J x r)/s^{3/2}.  Agx on X_EDGE
 * (i+0.5, j), Agy on Y_EDGE (i, j+0.5). */
static void fill_bg_gravomag_dipole(gr_sim_t* sim, float x0, float y0,
                                    float Jz, float epsilon) {
    float* Agx = ensure_bg_alloc(sim, GR_FIELD_A_GX);
    float* Agy = ensure_bg_alloc(sim, GR_FIELD_A_GY);
    if (!Agx || !Agy) return;
    const int W = sim->width, H = sim->height;
    const float dx = sim->dx, eps2 = epsilon * epsilon;
    const float coeff = sim->G_eff * Jz / (2.0f * sim->c_eff * sim->c_eff);
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx, dyc = y - y0;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xc = ((float) i + 0.5f) * dx, dxc = xc - x0;
            const float s2 = dxc * dxc + dyc * dyc + eps2;
            Agx[row + i] += -coeff * dyc / (s2 * sqrtf(s2));
        }
    }
    for (int j = 0; j < H; j++) {
        const float y = ((float) j + 0.5f) * dx, dyc = y - y0;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xc = (float) i * dx, dxc = xc - x0;
            const float s2 = dxc * dxc + dyc * dyc + eps2;
            Agy[row + i] += +coeff * dxc / (s2 * sqrtf(s2));
        }
    }
}

/* EM magnetic dipole A_em = (k_e mu_z / c^2) (z_hat x r)/s^{3/2} from a
 * magnetic moment mu_z (the EM analog of fill_bg_gravomag_dipole).  Ax on
 * X_EDGE (i+0.5, j), Ay on Y_EDGE (i, j+0.5). */
static void fill_bg_magnetic_dipole(gr_sim_t* sim, float x0, float y0,
                                    float mu_z, float epsilon) {
    float* Ax = ensure_bg_alloc(sim, GR_FIELD_A_X);
    float* Ay = ensure_bg_alloc(sim, GR_FIELD_A_Y);
    if (!Ax || !Ay) return;
    const int W = sim->width, H = sim->height;
    const float dx = sim->dx, eps2 = epsilon * epsilon;
    const float coeff = sim->k_e * mu_z / (sim->c_eff * sim->c_eff);
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx, dyc = y - y0;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xc = ((float) i + 0.5f) * dx, dxc = xc - x0;
            const float s2 = dxc * dxc + dyc * dyc + eps2;
            Ax[row + i] += -coeff * dyc / (s2 * sqrtf(s2));
        }
    }
    for (int j = 0; j < H; j++) {
        const float y = ((float) j + 0.5f) * dx, dyc = y - y0;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xc = (float) i * dx, dxc = xc - x0;
            const float s2 = dxc * dxc + dyc * dyc + eps2;
            Ay[row + i] += +coeff * dxc / (s2 * sqrtf(s2));
        }
    }
}

/* EM magnetic moment of a spinning charged body: mu_z = g (Q/2M) Jz, with the
 * parameter GM = G_eff*M => mu_z = g*Q*Jz*G_eff/(2*GM).  Zero if GM==0. */
static float compact_body_mu_z(const gr_sim_t* sim,
                               float GM, float Q, float Jz, float g_em) {
    if (GM == 0.0f) return 0.0f;
    return g_em * Q * Jz * sim->G_eff / (2.0f * GM);
}

/* Coulomb phi^{bg}(x) = +k_e Q / sqrt(r^2 + eps^2) on the CORNER sublattice. */
static void fill_bg_coulomb(gr_sim_t* sim, float x0, float y0,
                            float Q, float epsilon) {
    float* phi = ensure_bg_alloc(sim, GR_FIELD_PHI_EM);
    if (!phi) return;
    const int W = sim->width, H = sim->height;
    const float dx = sim->dx, eps2 = epsilon * epsilon;
    const float coeff = sim->k_e * Q;
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx, dy = y - y0;
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xx = (float) i * dx, dxi = xx - x0;
            phi[row + i] += coeff / sqrtf(dxi * dxi + dy * dy + eps2);
        }
    }
}

/* Eq. (eq:bg_softened_point_mass) — gr_sandbox_vNN.tex §12.6.
 *
 *   Phi_g^{bg}(x) = -G_eff * M / sqrt(|x - x0|^2 + epsilon^2)
 *
 * Sampled onto the CORNER sublattice: phi_g_bg[j*W + i] holds the value at
 * position (i, j) * dx (v35 §sec:yee_pivot).  The smoothing length epsilon
 * avoids the 1/r singularity; ~few cells is typical. */
void gr_sim_set_background_point_mass(gr_sim_t* sim,
                                      float x0, float y0,
                                      float GM, float epsilon) {
    if (!sim) return;
    gr_sim_clear_background(sim);   /* generators are additive — start from zero */
    fill_bg_grav_scalar(sim, x0, y0, GM, epsilon);
    /* Store generator parameters for the analytic-mode evaluator. */
    sim->bg_kind = GR_BG_KIND_POINT_MASS;
    sim->bg_x0   = x0;
    sim->bg_y0   = y0;
    sim->bg_GM   = GM;
    sim->bg_eps  = epsilon;
    sim->bg_Jz   = 0.0f;
}

/* Spinning softened point mass — fills Phi_g^{bg} via the same generator as
 * gr_sim_set_background_point_mass, then additionally fills the A_{g,x} and
 * A_{g,y} cell-centered arrays with the gravitomagnetic dipole field
 *
 *   A_g(x) = (G_eff/(2 c^2)) * J × r / (r^2 + epsilon^2)^{3/2}
 *
 * For a spin axis along +z and 2D in-plane positions this is:
 *
 *   A_{g,x}(x,y) = -(G_eff J_z / (2 c^2)) * dy / s^{3/2}
 *   A_{g,y}(x,y) = +(G_eff J_z / (2 c^2)) * dx / s^{3/2}
 *
 * with dx = x - x0, dy = y - y0, s = dx^2 + dy^2 + eps^2.  The 1/(2 c^2)
 * factor is the dipole-moment coefficient in the simulation's GEM
 * convention (vector-potential source equation
 * grad^2 A_g = -(4 pi G / c^2) j_mass, identical to EM up to sign;
 * the doc's spin-2 factor of 4 sits in the FORCE law, not in A_g itself). */
void gr_sim_set_background_spinning_point_mass(gr_sim_t* sim,
                                               float x0, float y0,
                                               float GM, float epsilon,
                                               float Jz) {
    if (!sim) return;
    gr_sim_clear_background(sim);   /* generators are additive — start from zero */
    /* Scalar (grav) piece + gravitomagnetic dipole, via the shared generators. */
    fill_bg_grav_scalar(sim, x0, y0, GM, epsilon);
    fill_bg_gravomag_dipole(sim, x0, y0, Jz, epsilon);
    sim->bg_kind = GR_BG_KIND_SPINNING_POINT_MASS;
    sim->bg_x0   = x0;
    sim->bg_y0   = y0;
    sim->bg_GM   = GM;
    sim->bg_eps  = epsilon;
    sim->bg_Jz   = Jz;
}

/* Uniform gravitomagnetic background.  Symmetric-gauge potentials produce
 * a spatially constant B_g_z = B_0:
 *   A_{g,x}(x,y) = -0.5 B_0 (y - y_0)   ;  on X_EDGE sublattice
 *   A_{g,y}(x,y) = +0.5 B_0 (x - x_0)   ;  on Y_EDGE sublattice
 *   B_g_z        = d/dx A_{g,y} - d/dy A_{g,x} = B_0    (uniform, by construction)
 *   Phi_g        = 0                                    (no scalar gravity)
 *
 * Stage 20 unit-isolation test for the gravitomagnetic Lorentz force; see
 * gr_sandbox_v35.tex §sec:geodesic_expansion eq:geodesic_expansion (line 938)
 * for the +4 v x B_g coefficient in the doc's GEM-with-spin-2-factor
 * convention. */
void gr_sim_set_background_uniform_gravitomagnetic(gr_sim_t* sim,
                                                   float x0, float y0,
                                                   float B0) {
    if (!sim) return;
    /* Clear any prior background fields, then re-allocate A_g arrays.
     * Phi_g_bg is intentionally left NULL — the SAMPLED-mode path treats
     * NULL as zero, which is exactly what we want here. */
    gr_sim_clear_background(sim);
    float* Agx = ensure_bg_alloc(sim, GR_FIELD_A_GX);
    float* Agy = ensure_bg_alloc(sim, GR_FIELD_A_GY);
    if (!Agx || !Agy) return;

    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    /* X_EDGE sublattice for Agx: nodes at (i + 0.5, j) * dx.
     *   A_{g,x}(x, y) = -0.5 B_0 (y - y_0). */
    for (int j = 0; j < H; j++) {
        const float y  = (float) j * dx;
        const float dy = y - y0;
        const int   row = j * W;
        const float val = -0.5f * B0 * dy;
        for (int i = 0; i < W; i++) {
            Agx[row + i] = val;     /* independent of x for this gauge */
        }
    }
    /* Y_EDGE sublattice for Agy: nodes at (i, j + 0.5) * dx.
     *   A_{g,y}(x, y) = +0.5 B_0 (x - x_0). */
    for (int j = 0; j < H; j++) {
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = (float) i * dx;
            const float dxc = x - x0;
            Agy[row + i] = 0.5f * B0 * dxc;     /* independent of y */
        }
    }

    sim->bg_kind = GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC;
    sim->bg_x0   = x0;
    sim->bg_y0   = y0;
    sim->bg_GM   = 0.0f;
    sim->bg_eps  = 0.0f;
    sim->bg_Jz   = 0.0f;
    sim->bg_B0   = B0;
}

/* Uniform EM magnetic background — EM analog of
 * gr_sim_set_background_uniform_gravitomagnetic.  Symmetric-gauge
 * potentials fill A_x_bg, A_y_bg arrays such that
 *   B_z = d/dx A_y - d/dy A_x = B_0   (uniform).
 * No scalar EM is installed; phi_bg is left NULL (sampled mode treats
 * NULL as zero). */
void gr_sim_set_background_uniform_magnetic(gr_sim_t* sim,
                                            float x0, float y0,
                                            float B0) {
    if (!sim) return;
    gr_sim_clear_background(sim);
    float* Ax = ensure_bg_alloc(sim, GR_FIELD_A_X);
    float* Ay = ensure_bg_alloc(sim, GR_FIELD_A_Y);
    if (!Ax || !Ay) return;

    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    /* X_EDGE sublattice for A_x: nodes at (i+0.5, j) * dx.
     *   A_x(x, y) = -0.5 B_0 (y - y_0). */
    for (int j = 0; j < H; j++) {
        const float y  = (float) j * dx;
        const float dy = y - y0;
        const int   row = j * W;
        const float val = -0.5f * B0 * dy;
        for (int i = 0; i < W; i++) {
            Ax[row + i] = val;
        }
    }
    /* Y_EDGE sublattice for A_y: nodes at (i, j+0.5) * dx.
     *   A_y(x, y) = +0.5 B_0 (x - x_0). */
    for (int j = 0; j < H; j++) {
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = (float) i * dx;
            const float dxc = x - x0;
            Ay[row + i] = 0.5f * B0 * dxc;
        }
    }

    sim->bg_kind  = GR_BG_KIND_UNIFORM_MAGNETIC;
    sim->bg_x0    = x0;
    sim->bg_y0    = y0;
    sim->bg_GM    = 0.0f;
    sim->bg_eps   = 0.0f;
    sim->bg_Jz    = 0.0f;
    sim->bg_B0    = 0.0f;
    sim->bg_B0_em = B0;
}

/* v40 Linear-B background.  A_x = 0; A_y(x) = B0 (x-x0) + 0.5 B' (x-x0)^2
 * gives B_z = d_x A_y = B0 + B' (x - x0). */
void gr_sim_set_background_linear_magnetic(gr_sim_t* sim,
                                            float x0, float y0,
                                            float B0, float B_prime) {
    if (!sim) return;
    gr_sim_clear_background(sim);
    float* Ax = ensure_bg_alloc(sim, GR_FIELD_A_X);
    float* Ay = ensure_bg_alloc(sim, GR_FIELD_A_Y);
    if (!Ax || !Ay) return;

    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    /* A_x at X_EDGE = 0 everywhere. */
    for (int j = 0; j < H; j++) {
        const int row = j * W;
        for (int i = 0; i < W; i++) Ax[row + i] = 0.0f;
    }
    /* A_y at Y_EDGE sublattice -- node positions (i, j+0.5)*dx. */
    for (int j = 0; j < H; j++) {
        const int row = j * W;
        for (int i = 0; i < W; i++) {
            const float xv  = (float) i * dx;
            const float xrm = xv - x0;
            Ay[row + i] = B0 * xrm + 0.5f * B_prime * xrm * xrm;
        }
    }

    sim->bg_kind        = GR_BG_KIND_LINEAR_MAGNETIC;
    sim->bg_x0          = x0;
    sim->bg_y0          = y0;
    sim->bg_GM          = 0.0f;
    sim->bg_eps         = 0.0f;
    sim->bg_Jz          = 0.0f;
    sim->bg_B0          = 0.0f;
    sim->bg_B0_em       = B0;
    sim->bg_B_prime_em  = B_prime;
}

/* Uniform EM electric background — fills the EM scalar phi_bg with a
 * linear potential whose negative gradient is the desired uniform field:
 *   phi^{bg}(x, y) = -( E_x (x - x_0) + E_y (y - y_0) ),
 *   -grad phi^{bg} = ( E_x, E_y )                       (uniform).
 * No A^{bg} is installed.  Used by Stage 24 to isolate the q E piece of
 * the EM Lorentz force from the v x B piece. */
void gr_sim_set_background_uniform_electric(gr_sim_t* sim,
                                            float x0, float y0,
                                            float Ex, float Ey) {
    if (!sim) return;
    gr_sim_clear_background(sim);
    float* phi = ensure_bg_alloc(sim, GR_FIELD_PHI_EM);
    if (!phi) return;

    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;

    /* CORNER sublattice for phi: nodes at (i, j) * dx. */
    for (int j = 0; j < H; j++) {
        const float y  = (float) j * dx;
        const float dyi = y - y0;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = (float) i * dx;
            const float dxi = x - x0;
            phi[row + i] = -(Ex * dxi + Ey * dyi);
        }
    }

    sim->bg_kind  = GR_BG_KIND_UNIFORM_ELECTRIC;
    sim->bg_x0    = x0;
    sim->bg_y0    = y0;
    sim->bg_GM    = 0.0f;
    sim->bg_eps   = 0.0f;
    sim->bg_Jz    = 0.0f;
    sim->bg_B0    = 0.0f;
    sim->bg_B0_em = 0.0f;
    sim->bg_Ex_em = Ex;
    sim->bg_Ey_em = Ey;
}

/* Softened point-charge background.  EM analog of point-mass:
 *   phi^{bg}(x, y) = +k_e * Q / sqrt(r^2 + epsilon^2)
 * sampled at CORNER sublattice nodes (i, j) * dx.  Like the mass case,
 * the smoothing length epsilon ~ a few cells avoids the 1/r singularity. */
void gr_sim_set_background_point_charge(gr_sim_t* sim,
                                        float x0, float y0,
                                        float Q, float epsilon) {
    if (!sim) return;
    gr_sim_clear_background(sim);
    fill_bg_coulomb(sim, x0, y0, Q, epsilon);

    sim->bg_kind  = GR_BG_KIND_POINT_CHARGE;
    sim->bg_x0    = x0;
    sim->bg_y0    = y0;
    sim->bg_GM    = 0.0f;
    sim->bg_eps   = epsilon;
    sim->bg_Jz    = 0.0f;
    sim->bg_B0    = 0.0f;
    sim->bg_B0_em = 0.0f;
    sim->bg_Ex_em = 0.0f;
    sim->bg_Ey_em = 0.0f;
    sim->bg_charge = Q;
}

/* Unified compact body (M, Q, J_z) -- see header.  Clears once, then
 * superposes the grav scalar, gravitomagnetic dipole, and Coulomb generators.
 * Components with a zero coefficient are skipped (the sampled path treats a
 * NULL array as zero; the analytic evaluators use the stored zero param). */
/* Accumulate one compact body's sampled field into the (already-allocated or
 * lazily-allocated, zero-initialized) background arrays.  Additive, so several
 * bodies superpose.  Skips zero-coefficient components. */
static void accumulate_body_sampled(gr_sim_t* sim, const gr_bg_body_t* b) {
    const float mu_z = compact_body_mu_z(sim, b->GM, b->Q, b->Jz, b->g_em);
    if (b->GM != 0.0f) fill_bg_grav_scalar(sim, b->x, b->y, b->GM, b->eps);
    if (b->Jz != 0.0f) fill_bg_gravomag_dipole(sim, b->x, b->y, b->Jz, b->eps);
    if (b->Q  != 0.0f) fill_bg_coulomb(sim, b->x, b->y, b->Q, b->eps);
    if (mu_z  != 0.0f) fill_bg_magnetic_dipole(sim, b->x, b->y, mu_z, b->eps);
}

/* Re-derive the sampled arrays from the current compact-body list (after an
 * in-place edit / removal): drop the arrays, then re-accumulate every body. */
static void rebuild_sampled_from_bodies(gr_sim_t* sim) {
    free_bg_arrays(sim);
    for (int i = 0; i < sim->bg_nbody; i++)
        accumulate_body_sampled(sim, &sim->bg_bodies[i]);
}

int gr_sim_add_background_body(gr_sim_t* sim,
                               float x0, float y0,
                               float GM, float Q, float Jz,
                               float g_em, float epsilon) {
    if (!sim) return -1;
    if (sim->bg_nbody >= GR_MAX_BG_BODIES) return -1;
    const int idx = sim->bg_nbody;
    gr_bg_body_t* b = &sim->bg_bodies[idx];
    b->x = x0; b->y = y0; b->GM = GM; b->Q = Q; b->Jz = Jz; b->g_em = g_em; b->eps = epsilon;
    accumulate_body_sampled(sim, b);   /* additive into the sampled arrays */
    sim->bg_nbody = idx + 1;
    sim->bg_kind  = GR_BG_KIND_COMPACT_BODY;
    /* Mirror the first body into the legacy scalars (informational; the
     * COMPACT_BODY evaluators read bg_bodies, not these). */
    if (idx == 0) {
        sim->bg_x0 = x0; sim->bg_y0 = y0; sim->bg_GM = GM; sim->bg_eps = epsilon;
        sim->bg_Jz = Jz; sim->bg_charge = Q; sim->bg_g_em = g_em;
        sim->bg_B0 = 0.0f; sim->bg_B0_em = 0.0f; sim->bg_Ex_em = 0.0f; sim->bg_Ey_em = 0.0f;
    }
    return idx;
}

int gr_sim_background_count(const gr_sim_t* sim) {
    return sim ? sim->bg_nbody : 0;
}

float* gr_sim_get_background_body_ptr(gr_sim_t* sim, int i) {
    if (!sim || i < 0 || i >= sim->bg_nbody) return NULL;
    return &sim->bg_bodies[i].x;   /* 7 contiguous floats */
}

void gr_sim_set_background_body_at(gr_sim_t* sim, int i,
                                   float x0, float y0,
                                   float GM, float Q, float Jz,
                                   float g_em, float epsilon) {
    if (!sim || i < 0 || i >= sim->bg_nbody) return;
    gr_bg_body_t* b = &sim->bg_bodies[i];
    b->x = x0; b->y = y0; b->GM = GM; b->Q = Q; b->Jz = Jz; b->g_em = g_em; b->eps = epsilon;
    rebuild_sampled_from_bodies(sim);   /* the edit changes the whole field */
}

/* Single compact body — back-compat convenience: clear, then add one body. */
void gr_sim_set_background_body(gr_sim_t* sim,
                                float x0, float y0,
                                float GM, float Q, float Jz,
                                float g_em, float epsilon) {
    if (!sim) return;
    gr_sim_clear_background(sim);
    gr_sim_add_background_body(sim, x0, y0, GM, Q, Jz, g_em, epsilon);
}

/* Analytic-mode evaluation of the installed background generator at the
 * particle's exact (x, y).  See gr_sandbox §12.6 / §12.8 for the rationale:
 * the sampled CIC+FD path introduces an O((dx/r)^2) tangential force error
 * that breaks angular-momentum conservation in test-particle orbits.  The
 * analytic path replaces that with the closed-form expression, eliminating
 * the grid-induced precession in Stage 7/8 tests entirely. */
int gr_bg_eval_analytic(const struct gr_sim* sim, float x, float y,
                        float* phi_out, float* gx_out, float* gy_out) {
    if (!sim || sim->bg_kind == GR_BG_KIND_NONE) return 0;

    switch (sim->bg_kind) {
    case GR_BG_KIND_POINT_MASS:
    case GR_BG_KIND_SPINNING_POINT_MASS: {
        /* Scalar (grav) piece uses bg_GM for the single-body kinds. */
        /*   Phi(x,y) = -G*M / sqrt(r^2 + eps^2)                    */
        /*   grad     =  G*M * (r - r0) / (r^2 + eps^2)^{3/2}       */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        const float r2  = dxi * dxi + dyi * dyi;
        const float s2  = r2 + sim->bg_eps * sim->bg_eps;
        const float inv_s  = 1.0f / sqrtf(s2);
        const float inv_s3 = inv_s / s2;            /* 1 / (s2)^{3/2} */
        *phi_out = -sim->bg_GM * inv_s;
        *gx_out  =  sim->bg_GM * dxi * inv_s3;
        *gy_out  =  sim->bg_GM * dyi * inv_s3;
        return 1;
    }
    case GR_BG_KIND_COMPACT_BODY: {
        /* Sum the grav scalar + gradient over all compact bodies. */
        float phi = 0.0f, gx = 0.0f, gy = 0.0f;
        for (int k = 0; k < sim->bg_nbody; k++) {
            const gr_bg_body_t* b = &sim->bg_bodies[k];
            const float dxi = x - b->x, dyi = y - b->y;
            const float s2  = dxi * dxi + dyi * dyi + b->eps * b->eps;
            const float inv_s = 1.0f / sqrtf(s2), inv_s3 = inv_s / s2;
            phi += -b->GM * inv_s;
            gx  +=  b->GM * dxi * inv_s3;
            gy  +=  b->GM * dyi * inv_s3;
        }
        *phi_out = phi; *gx_out = gx; *gy_out = gy;
        return 1;
    }
    case GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC:
    case GR_BG_KIND_UNIFORM_MAGNETIC:
    case GR_BG_KIND_UNIFORM_ELECTRIC:
    case GR_BG_KIND_POINT_CHARGE: {
        /* No scalar gravity in these backgrounds. */
        *phi_out = 0.0f;
        *gx_out  = 0.0f;
        *gy_out  = 0.0f;
        (void) x; (void) y;
        return 1;
    }
    default:
        return 0;
    }
}

int gr_bg_eval_A_g(const struct gr_sim* sim, float x, float y,
                   float* Ax_out, float* Ay_out) {
    if (!sim) {
        *Ax_out = 0.0f;
        *Ay_out = 0.0f;
        return 0;
    }
    switch (sim->bg_kind) {
    case GR_BG_KIND_SPINNING_POINT_MASS: {
        /* Same dipole formula used in the sampler (bg_Jz=0 => zero A_g). */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        const float s2  = dxi * dxi + dyi * dyi + sim->bg_eps * sim->bg_eps;
        const float inv_s3 = 1.0f / (s2 * sqrtf(s2));
        const float coeff  = sim->G_eff * sim->bg_Jz
                           / (2.0f * sim->c_eff * sim->c_eff);
        *Ax_out = -coeff * dyi * inv_s3;
        *Ay_out =  coeff * dxi * inv_s3;
        return 1;
    }
    case GR_BG_KIND_COMPACT_BODY: {
        /* Sum the gravitomagnetic dipole over all compact bodies. */
        const float k = sim->G_eff / (2.0f * sim->c_eff * sim->c_eff);
        float ax = 0.0f, ay = 0.0f;
        for (int n = 0; n < sim->bg_nbody; n++) {
            const gr_bg_body_t* b = &sim->bg_bodies[n];
            if (b->Jz == 0.0f) continue;
            const float dxi = x - b->x, dyi = y - b->y;
            const float s2  = dxi * dxi + dyi * dyi + b->eps * b->eps;
            const float inv_s3 = 1.0f / (s2 * sqrtf(s2));
            const float coeff = k * b->Jz;
            ax += -coeff * dyi * inv_s3;
            ay +=  coeff * dxi * inv_s3;
        }
        *Ax_out = ax; *Ay_out = ay;
        return 1;
    }
    case GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC: {
        /* Symmetric-gauge form:
         *   A_{g,x} = -0.5 B_0 (y - y_0),  A_{g,y} = +0.5 B_0 (x - x_0). */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        *Ax_out = -0.5f * sim->bg_B0 * dyi;
        *Ay_out =  0.5f * sim->bg_B0 * dxi;
        return 1;
    }
    default:
        *Ax_out = 0.0f;
        *Ay_out = 0.0f;
        return 0;
    }
}

/* Analytic-mode B_g_z(x, y) = (curl A_g)_z = d/dx A_{g,y} - d/dy A_{g,x}. */
int gr_bg_eval_B_g(const struct gr_sim* sim, float x, float y,
                   float* Bgz_out) {
    if (!sim) {
        *Bgz_out = 0.0f;
        return 0;
    }
    switch (sim->bg_kind) {
    case GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC: {
        /* Constant by construction. */
        (void) x; (void) y;
        *Bgz_out = sim->bg_B0;
        return 1;
    }
    case GR_BG_KIND_SPINNING_POINT_MASS: {
        /* Differentiating A_g = coeff * (J × r) / s^{3/2} gives
         *   A_{g,x} = -k dy / s^{3/2},     A_{g,y} =  k dx / s^{3/2}
         *   s       = dx^2 + dy^2 + eps^2,   k = G J_z / (2 c^2)
         *   B_g_z = d/dx A_{g,y} - d/dy A_{g,x}
         *         = k (2 eps^2 - (dx^2 + dy^2)) / s^{5/2}. */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        const float r2  = dxi * dxi + dyi * dyi;
        const float eps2 = sim->bg_eps * sim->bg_eps;
        const float s2  = r2 + eps2;
        const float inv_s52 = 1.0f / (s2 * s2 * sqrtf(s2));
        const float coeff  = sim->G_eff * sim->bg_Jz
                           / (2.0f * sim->c_eff * sim->c_eff);
        *Bgz_out = coeff * (2.0f * eps2 - r2) * inv_s52;
        return 1;
    }
    case GR_BG_KIND_COMPACT_BODY: {
        /* Sum B_g_z = k (2 eps^2 - r^2) / s^{5/2} over all compact bodies. */
        const float k = sim->G_eff / (2.0f * sim->c_eff * sim->c_eff);
        float bz = 0.0f;
        for (int n = 0; n < sim->bg_nbody; n++) {
            const gr_bg_body_t* b = &sim->bg_bodies[n];
            if (b->Jz == 0.0f) continue;
            const float dxi = x - b->x, dyi = y - b->y;
            const float r2 = dxi * dxi + dyi * dyi;
            const float eps2 = b->eps * b->eps;
            const float s2 = r2 + eps2;
            const float inv_s52 = 1.0f / (s2 * s2 * sqrtf(s2));
            bz += k * b->Jz * (2.0f * eps2 - r2) * inv_s52;
        }
        *Bgz_out = bz;
        return 1;
    }
    default:
        *Bgz_out = 0.0f;
        return 0;
    }
}

/* Analytic-mode evaluation of the EM vector potential A_em(x, y). */
int gr_bg_eval_A_em(const struct gr_sim* sim, float x, float y,
                    float* Ax_out, float* Ay_out) {
    if (!sim) {
        *Ax_out = 0.0f;
        *Ay_out = 0.0f;
        return 0;
    }
    switch (sim->bg_kind) {
    case GR_BG_KIND_UNIFORM_MAGNETIC: {
        /* Symmetric-gauge form:
         *   A_x = -0.5 B_0 (y - y_0),  A_y = +0.5 B_0 (x - x_0). */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        *Ax_out = -0.5f * sim->bg_B0_em * dyi;
        *Ay_out =  0.5f * sim->bg_B0_em * dxi;
        return 1;
    }
    case GR_BG_KIND_LINEAR_MAGNETIC: {
        /* A_x = 0; A_y = B0 (x - x0) + 0.5 B' (x - x0)^2 ->
         * B_z = d_x A_y = B0 + B' (x - x0). */
        const float dxi = x - sim->bg_x0;
        (void) y;
        *Ax_out = 0.0f;
        *Ay_out = sim->bg_B0_em * dxi + 0.5f * sim->bg_B_prime_em * dxi * dxi;
        return 1;
    }
    case GR_BG_KIND_COMPACT_BODY: {
        /* Sum the magnetic dipole of each spinning charged body:
         * A_em = (k_e mu_z/c^2)(z x r)/s^{3/2}, mirroring A_g (mu_z=0 => zero). */
        const float kc = sim->k_e / (sim->c_eff * sim->c_eff);
        float ax = 0.0f, ay = 0.0f;
        for (int n = 0; n < sim->bg_nbody; n++) {
            const gr_bg_body_t* b = &sim->bg_bodies[n];
            const float mu_z = compact_body_mu_z(sim, b->GM, b->Q, b->Jz, b->g_em);
            if (mu_z == 0.0f) continue;
            const float dxi = x - b->x, dyi = y - b->y;
            const float s2  = dxi * dxi + dyi * dyi + b->eps * b->eps;
            const float inv_s3 = 1.0f / (s2 * sqrtf(s2));
            const float coeff  = kc * mu_z;
            ax += -coeff * dyi * inv_s3;
            ay +=  coeff * dxi * inv_s3;
        }
        *Ax_out = ax; *Ay_out = ay;
        return 1;
    }
    default:
        *Ax_out = 0.0f;
        *Ay_out = 0.0f;
        return 0;
    }
}

/* Analytic-mode B_em_z(x, y) = (curl A)_z = d/dx A_y - d/dy A_x. */
int gr_bg_eval_B_em(const struct gr_sim* sim, float x, float y,
                    float* Bz_out) {
    if (!sim) {
        *Bz_out = 0.0f;
        return 0;
    }
    switch (sim->bg_kind) {
    case GR_BG_KIND_UNIFORM_MAGNETIC: {
        /* Constant by construction. */
        (void) x; (void) y;
        *Bz_out = sim->bg_B0_em;
        return 1;
    }
    case GR_BG_KIND_LINEAR_MAGNETIC: {
        /* B_z = B0 + B' (x - x0). */
        (void) y;
        *Bz_out = sim->bg_B0_em + sim->bg_B_prime_em * (x - sim->bg_x0);
        return 1;
    }
    case GR_BG_KIND_COMPACT_BODY: {
        /* Curl of the magnetic dipole summed over bodies (mu_z=0 => zero). */
        const float kc = sim->k_e / (sim->c_eff * sim->c_eff);
        float bz = 0.0f;
        for (int n = 0; n < sim->bg_nbody; n++) {
            const gr_bg_body_t* b = &sim->bg_bodies[n];
            const float mu_z = compact_body_mu_z(sim, b->GM, b->Q, b->Jz, b->g_em);
            if (mu_z == 0.0f) continue;
            const float dxi = x - b->x, dyi = y - b->y;
            const float r2  = dxi * dxi + dyi * dyi;
            const float eps2 = b->eps * b->eps;
            const float s2  = r2 + eps2;
            const float inv_s52 = 1.0f / (s2 * s2 * sqrtf(s2));
            bz += kc * mu_z * (2.0f * eps2 - r2) * inv_s52;
        }
        *Bz_out = bz;
        return 1;
    }
    default:
        *Bz_out = 0.0f;
        return 0;
    }
}

/* Evaluate Phi_g^{bg}(x, y) at an arbitrary point using whichever path is
 * available: prefer the analytic generator (exact), else CIC-interpolate
 * the sampled phi_g_bg array on the CORNER sublattice, else return 0.
 * Used by the Shapiro c_local^2 recompute below.  Out-of-range positions
 * (outside the grid) snap to the nearest cell-edge sample. */
static float bg_phi_g_at(const struct gr_sim* sim, float xp, float yp) {
    if (!sim) return 0.0f;
    /* Analytic path: only the kinds that set phi_g_bg (POINT_MASS,
     * SPINNING_POINT_MASS) return a nonzero Phi_g. */
    if (sim->bg_kind == GR_BG_KIND_POINT_MASS
        || sim->bg_kind == GR_BG_KIND_SPINNING_POINT_MASS) {
        float phi, gx, gy;
        if (gr_bg_eval_analytic(sim, xp, yp, &phi, &gx, &gy)) {
            return phi;
        }
    }
    /* Fall back to sampled CORNER array with bilinear interpolation. */
    if (!sim->phi_g_bg) return 0.0f;
    const int   W  = sim->width;
    const int   H  = sim->height;
    const float dx = sim->dx;
    /* CORNER nodes at (i, j) * dx. */
    float u = xp / dx;
    float v = yp / dx;
    if (u < 0.0f) u = 0.0f;
    if (v < 0.0f) v = 0.0f;
    if (u > (float) (W - 1)) u = (float) (W - 1);
    if (v > (float) (H - 1)) v = (float) (H - 1);
    int i0 = (int) u;
    int j0 = (int) v;
    if (i0 > W - 2) i0 = W - 2;
    if (j0 > H - 2) j0 = H - 2;
    const float fx = u - (float) i0;
    const float fy = v - (float) j0;
    const float* a = sim->phi_g_bg;
    const float v00 = a[ j0      * W + i0    ];
    const float v10 = a[ j0      * W + i0 + 1];
    const float v01 = a[(j0 + 1) * W + i0    ];
    const float v11 = a[(j0 + 1) * W + i0 + 1];
    const float wy0 = 1.0f - fy;
    const float wy1 = fy;
    const float wx0 = 1.0f - fx;
    const float wx1 = fx;
    return wy0 * (wx0 * v00 + wx1 * v10) + wy1 * (wx0 * v01 + wx1 * v11);
}

/* Shapiro delay (Stage 31+): the EM wave equation propagates at
 *
 *   c_local(x) = c * (1 + 2 Phi_g(x) / c^2)
 *
 * (gr_sandbox_v35.tex sec:shapiro eq:c_local).  We precompute c_local^2
 * sampled at each EM-field sublattice's own node positions so the
 * field-evolution kernel can read it without per-step interpolation.
 *
 * Floor:  for the weak-field regime in which the c_local approximation
 * is valid, (1 + 2 Phi_g/c^2) > 0.  We clamp to 1e-3 as a defensive
 * floor against pathological strong-field probes that would otherwise
 * make c_local^2 vanish or go negative; in that regime the linearized
 * approximation is no longer trustworthy anyway. */
void gr_em_shapiro_recompute_c_local2(struct gr_sim* sim) {
    if (!sim) return;
    const int   W   = sim->width;
    const int   H   = sim->height;
    const float dx  = sim->dx;
    const size_t n  = (size_t) W * (size_t) H;
    const float c   = sim->c_eff;
    const float c2  = c * c;
    const float inv_c2 = 1.0f / c2;

    if (!sim->c_local2_corner) sim->c_local2_corner = (float*) calloc(n, sizeof(float));
    if (!sim->c_local2_xedge)  sim->c_local2_xedge  = (float*) calloc(n, sizeof(float));
    if (!sim->c_local2_yedge)  sim->c_local2_yedge  = (float*) calloc(n, sizeof(float));
    if (!sim->c_local2_corner || !sim->c_local2_xedge || !sim->c_local2_yedge) return;

    /* Include the perturbation Phi_g (deposited masses) on top of the analytic
     * background, so a deposited particle's (2D log r) well also bends the EM
     * wave -- not just the 1/r background.  Phi_g^pert lives on the CORNER grid;
     * the edge sublattices average adjacent corners.  This is recomputed per
     * step when em_shapiro_dynamic is set (a moving lens). */
    const float* phig = sim->fields[GR_FIELD_PHI_GRAV].curr;

    /* CORNER sublattice: nodes at (i, j) * dx -- for phi_em. */
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = (float) i * dx;
            const float phi = bg_phi_g_at(sim, x, y) + (phig ? phig[row + i] : 0.0f);
            float f = 1.0f + 2.0f * phi * inv_c2;
            if (f < 1.0e-3f) f = 1.0e-3f;
            sim->c_local2_corner[row + i] = c2 * f * f;
        }
    }
    /* X_EDGE sublattice: nodes at (i + 0.5, j) * dx -- for A_x. */
    for (int j = 0; j < H; j++) {
        const float y = (float) j * dx;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = ((float) i + 0.5f) * dx;
            const float pp  = phig ? 0.5f * (phig[row + i] + phig[row + (i + 1 < W ? i + 1 : i)]) : 0.0f;
            const float phi = bg_phi_g_at(sim, x, y) + pp;
            float f = 1.0f + 2.0f * phi * inv_c2;
            if (f < 1.0e-3f) f = 1.0e-3f;
            sim->c_local2_xedge[row + i] = c2 * f * f;
        }
    }
    /* Y_EDGE sublattice: nodes at (i, j + 0.5) * dx -- for A_y. */
    for (int j = 0; j < H; j++) {
        const float y = ((float) j + 0.5f) * dx;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = (float) i * dx;
            const float pp  = phig ? 0.5f * (phig[row + i] + phig[(j + 1 < H ? j + 1 : j) * W + i]) : 0.0f;
            const float phi = bg_phi_g_at(sim, x, y) + pp;
            float f = 1.0f + 2.0f * phi * inv_c2;
            if (f < 1.0e-3f) f = 1.0e-3f;
            sim->c_local2_yedge[row + i] = c2 * f * f;
        }
    }
}

/* Analytic-mode evaluation of the EM scalar potential phi^{bg}(x, y)
 * and its gradient. */
int gr_bg_eval_phi_em(const struct gr_sim* sim, float x, float y,
                      float* phi_out, float* gx_out, float* gy_out) {
    if (!sim) {
        *phi_out = 0.0f;
        *gx_out  = 0.0f;
        *gy_out  = 0.0f;
        return 0;
    }
    switch (sim->bg_kind) {
    case GR_BG_KIND_UNIFORM_ELECTRIC: {
        /* phi = -( Ex (x-x0) + Ey (y-y0) ) ; grad = (-Ex, -Ey). */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        *phi_out = -(sim->bg_Ex_em * dxi + sim->bg_Ey_em * dyi);
        *gx_out  = -sim->bg_Ex_em;
        *gy_out  = -sim->bg_Ey_em;
        return 1;
    }
    case GR_BG_KIND_POINT_CHARGE: {
        /* phi(x,y) = +k_e Q / sqrt(r^2 + eps^2)  (bg_charge=0 => zero E),
         * grad     = -k_e Q (r - r0) / (r^2 + eps^2)^{3/2}. */
        const float dxi = x - sim->bg_x0;
        const float dyi = y - sim->bg_y0;
        const float r2  = dxi * dxi + dyi * dyi;
        const float s2  = r2 + sim->bg_eps * sim->bg_eps;
        const float inv_s  = 1.0f / sqrtf(s2);
        const float inv_s3 = inv_s / s2;
        const float coeff  = sim->k_e * sim->bg_charge;
        *phi_out = coeff * inv_s;
        *gx_out  = -coeff * dxi * inv_s3;
        *gy_out  = -coeff * dyi * inv_s3;
        return 1;
    }
    case GR_BG_KIND_COMPACT_BODY: {
        /* Sum the Coulomb scalar + gradient over all charged bodies. */
        float phi = 0.0f, gx = 0.0f, gy = 0.0f;
        for (int n = 0; n < sim->bg_nbody; n++) {
            const gr_bg_body_t* b = &sim->bg_bodies[n];
            if (b->Q == 0.0f) continue;
            const float dxi = x - b->x, dyi = y - b->y;
            const float s2  = dxi * dxi + dyi * dyi + b->eps * b->eps;
            const float inv_s = 1.0f / sqrtf(s2), inv_s3 = inv_s / s2;
            const float coeff = sim->k_e * b->Q;
            phi += coeff * inv_s;
            gx  += -coeff * dxi * inv_s3;
            gy  += -coeff * dyi * inv_s3;
        }
        *phi_out = phi; *gx_out = gx; *gy_out = gy;
        return 1;
    }
    default:
        *phi_out = 0.0f;
        *gx_out  = 0.0f;
        *gy_out  = 0.0f;
        return 0;
    }
}
