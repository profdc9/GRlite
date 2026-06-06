/* Simulation + display configuration constants.
 *
 * These are the current hard-coded defaults carried over from the Stage-1
 * monolith.  Phase 2 moves the per-scenario subset of these into the JSON
 * scenario; the grid / CFL / colormap defaults stay here. */

export const GRID_W = 384;
export const GRID_H = 384;
export const DX = 1.0;
export const C_EFF = 1.0;
export const CFL = 1.0 / Math.sqrt(2);

/* Absorbing-ring depth (cells) used for the friction wall ramp. */
export const N_DAMPING = 48;

export const STEPS_PER_FRAME = 4;

/* Field-init settle length: after the boundary-mean Liénard-Wiechert direct
 * sum, run this many frozen-friction steps to relax to the discrete fixed
 * point (see gr_sandbox_v38_addendum §field init). */
export const N_SETTLE = Math.max(GRID_W, GRID_H);

/* Deposition / force enum values (mirror grlite.h). */
export const GR_SHAPE_CIC = 0;
export const GR_SHAPE_TSC = 1;
export const GR_SHAPE_BUMP = 2;
export const GR_FORCE_INTERP_LEGACY = 0;
export const GR_FORCE_INTERP_LB = 1;

/* GR field ids (mirror gr_field_id_t in grlite.h). */
export const GR_FIELD_PHI_GRAV = 0;
export const GR_FIELD_A_GX = 1;
export const GR_FIELD_A_GY = 2;
export const GR_FIELD_PHI_EM = 3;
export const GR_FIELD_A_X = 4;
export const GR_FIELD_A_Y = 5;

/* Production deposition/force config (stage63: BUMP R=8 pins the binary COM). */
export const KERNEL_RADIUS_CELLS = 8.0;

/* v42 absorber: Dirichlet outer BC + wall-ramp-only derivative friction
 * (interior floor 0 to avoid bleeding orbital near-field energy). */
export const FRICTION_FLOOR = 0.0;
export const FRICTION_WALL = 0.02;

/* Fallback colormap gain when the field is identically zero. */
export const DISPLAY_SCALE = 5.0;

/* Particle struct stride (floats) in gr_particle_t.  7 pre-v40, 10 with spin. */
export const PARTICLE_STRIDE_F32_OLD = 7;
export const PARTICLE_STRIDE_F32_NEW = 10;

/* Per-scenario C-scenario parameter packages (Phase 2 replaces this with
 * JSON-driven building). */
export interface ScenarioSpec {
    name: string;
    params: number[];
}
export const SCENARIOS: { [k: string]: ScenarioSpec } = {
    pic_binary: {
        name: 'pic_binary',
        /* mass, orbital radius, v_factor.  v_factor=1.05 is the BUMP-R8
         * calibrated near-circular speed (stage63). */
        params: [0.01, 15.0 * DX, 1.05],
    },
    pic_static: {
        name: 'pic_static',
        params: [0.01],
    },
};
