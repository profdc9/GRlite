/* Central app state: ONLY genuine runtime/UI state that is not part of the
 * scenario.  All view settings (field/source/vectors/trails/velocity) and the
 * scenario name live on the canonical Scenario (current.view / current.name) --
 * that is their single source of truth; mirroring them here caused drift.  The
 * render loop reads view settings from `current`; UI + bridge mutate `current`. */

/* What the inspector is editing: a particle or a background body (or nothing).
 * One selection at a time -- selecting one kind clears the other. */
export type SelKind = 'particle' | 'body';
export interface Selection { kind: SelKind; index: number; }

export interface AppState {
    paused: boolean;
    singleStep: boolean;
    /* Set when a live edit has diverged the running sim from the canonical
     * scenario -- run is no longer reproducible until reset. */
    liveModified: boolean;
    /* The selected object (particle or background body), or null. */
    selection: Selection | null;
}

export function createState(): AppState {
    return {
        paused: true,
        singleStep: false,
        liveModified: false,
        selection: null,
    };
}

/* The selected index of a given kind, or -1 if a different kind / nothing is
 * selected. */
export function selectedOf(state: AppState, kind: SelKind): number {
    return state.selection && state.selection.kind === kind ? state.selection.index : -1;
}
