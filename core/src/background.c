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

void gr_sim_clear_background(gr_sim_t* sim) {
    if (!sim) return;
    free(sim->phi_g_bg); sim->phi_g_bg = NULL;
    free(sim->Agx_bg);   sim->Agx_bg   = NULL;
    free(sim->Agy_bg);   sim->Agy_bg   = NULL;
    free(sim->phi_bg);   sim->phi_bg   = NULL;
    free(sim->Ax_bg);    sim->Ax_bg    = NULL;
    free(sim->Ay_bg);    sim->Ay_bg    = NULL;
    /* Clear analytic generator parameters too — they're paired with the
     * sampled array.  Mode is left alone (user-set). */
    sim->bg_kind   = GR_BG_KIND_NONE;
    sim->bg_x0     = 0.0f;
    sim->bg_y0     = 0.0f;
    sim->bg_GM     = 0.0f;
    sim->bg_eps    = 0.0f;
    sim->bg_charge = 0.0f;
    sim->bg_Jz     = 0.0f;
    sim->bg_B0     = 0.0f;
    sim->bg_B0_em  = 0.0f;
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
    float* phi_bg = ensure_bg_alloc(sim, GR_FIELD_PHI_GRAV);
    if (!phi_bg) return;

    const int   W    = sim->width;
    const int   H    = sim->height;
    const float dx   = sim->dx;
    const float eps2 = epsilon * epsilon;

    /* CORNER sublattice: nodes at (i, j) * dx — no offset. */
    for (int j = 0; j < H; j++) {
        const float y  = (float) j * dx;
        const float dy = y - y0;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float x   = (float) i * dx;
            const float dxi = x - x0;
            const float r2  = dxi * dxi + dy * dy;
            phi_bg[row + i] = -GM / sqrtf(r2 + eps2);
        }
    }
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
    /* First fill the scalar piece using the same generator as the
     * non-spinning case (and stash kind/params), then overwrite kind and
     * fill the vector potential arrays. */
    gr_sim_set_background_point_mass(sim, x0, y0, GM, epsilon);
    sim->bg_kind = GR_BG_KIND_SPINNING_POINT_MASS;
    sim->bg_Jz   = Jz;

    float* Agx = ensure_bg_alloc(sim, GR_FIELD_A_GX);
    float* Agy = ensure_bg_alloc(sim, GR_FIELD_A_GY);
    if (!Agx || !Agy) return;

    const int   W      = sim->width;
    const int   H      = sim->height;
    const float dx     = sim->dx;
    const float eps2   = epsilon * epsilon;
    /* Dipole coefficient: A_g = (k/(s^{3/2})) * (J × r),
     * with k = G_eff J_z / (2 c^2).  We also include the G_eff factor here
     * for consistency with the rest of the GEM-source convention. */
    const float coeff  = sim->G_eff * Jz
                       / (2.0f * sim->c_eff * sim->c_eff);

    /* X_EDGE sublattice for Agx: nodes at (i + 0.5, j) * dx. */
    for (int j = 0; j < H; j++) {
        const float y  = (float) j * dx;
        const float dyc = y - y0;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float xc   = ((float) i + 0.5f) * dx;
            const float dxc  = xc - x0;
            const float s2   = dxc * dxc + dyc * dyc + eps2;
            const float inv_s3 = 1.0f / (s2 * sqrtf(s2));
            Agx[row + i] = -coeff * dyc * inv_s3;
        }
    }
    /* Y_EDGE sublattice for Agy: nodes at (i, j + 0.5) * dx. */
    for (int j = 0; j < H; j++) {
        const float y  = ((float) j + 0.5f) * dx;
        const float dyc = y - y0;
        const int   row = j * W;
        for (int i = 0; i < W; i++) {
            const float xc   = (float) i * dx;
            const float dxc  = xc - x0;
            const float s2   = dxc * dxc + dyc * dyc + eps2;
            const float inv_s3 = 1.0f / (s2 * sqrtf(s2));
            Agy[row + i] = +coeff * dxc * inv_s3;
        }
    }
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
        /* Scalar piece is identical for both kinds — only A_g differs. */
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
        (void) x; (void) y;
        return 1;
    }
    case GR_BG_KIND_UNIFORM_GRAVITOMAGNETIC:
    case GR_BG_KIND_UNIFORM_MAGNETIC: {
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
        /* Same dipole formula used in the sampler. */
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
         *   B_g_z = coeff * (2 r^2 + 3 eps^2 - r^2) / s^{5/2}  ... let me just
         * recompute from scratch using the explicit form below.
         *
         *   A_{g,x} = -k dy / s^{3/2},     A_{g,y} =  k dx / s^{3/2}
         *   s       = dx^2 + dy^2 + eps^2,   k = G J_z / (2 c^2)
         *   d/dx A_{g,y} = k [ 1/s^{3/2} - 3 dx^2 / s^{5/2} ]
         *   d/dy A_{g,x} = k [ -1/s^{3/2} + 3 dy^2 / s^{5/2} ]
         *   B_g_z = d/dx A_{g,y} - d/dy A_{g,x}
         *         = k [ 2/s^{3/2} - 3 (dx^2 + dy^2) / s^{5/2} ]
         *         = k (2 s - 3 (dx^2 + dy^2)) / s^{5/2}
         *         = k (2 eps^2 - (dx^2 + dy^2)) / s^{5/2}    after grouping
         *  (using 2 s = 2(dx^2 + dy^2) + 2 eps^2). */
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
    default:
        *Bz_out = 0.0f;
        return 0;
    }
}
