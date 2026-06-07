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

export const STEPS_PER_FRAME = 4;

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

/* Fallback colormap gain when the field is identically zero. */
export const DISPLAY_SCALE = 5.0;

/* Particle struct stride (floats) in gr_particle_t.  7 pre-v40, 10 with spin,
 * 18 with the v43 force breakdown (fgrav/fem/fspin/ftot), 23 with the v43 drive
 * (drive_vx/vy, drive_omega, drive_phase, forces_enabled). */
export const PARTICLE_STRIDE_F32_OLD = 7;
export const PARTICLE_STRIDE_F32_NEW = 23;
