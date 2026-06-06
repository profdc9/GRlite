/* Central app state: ONLY genuine runtime/UI state that is not part of the
 * scenario.  All view settings (field/source/vectors/trails/velocity) and the
 * scenario name live on the canonical Scenario (current.view / current.name) --
 * that is their single source of truth; mirroring them here caused drift.  The
 * render loop reads view settings from `current`; UI + bridge mutate `current`. */

export interface AppState {
    paused: boolean;
    singleStep: boolean;
    /* Set when a live edit has diverged the running sim from the canonical
     * scenario -- run is no longer reproducible until reset. */
    liveModified: boolean;
    /* Index of the selected particle for the inspector, or -1. */
    selected: number;
}

export function createState(): AppState {
    return {
        paused: true,
        singleStep: false,
        liveModified: false,
        selected: -1,
    };
}
