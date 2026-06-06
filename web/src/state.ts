/* Central app state.  The render loop reads it; UI and (later) the debug
 * bridge mutate it.  Kept deliberately small and plain. */

import { GR_FIELD_PHI_GRAV } from './sim/config';

export interface AppState {
    scenario: string;
    paused: boolean;
    singleStep: boolean;
    /* Which field the FieldPass displays. */
    viewField: number;
    /* Set when a live edit has diverged the running sim from the canonical
     * scenario -- run is no longer reproducible until reset. */
    liveModified: boolean;
    /* Index of the selected particle for the inspector, or -1. */
    selected: number;
}

export function createState(scenario: string): AppState {
    return {
        scenario,
        paused: true,
        singleStep: false,
        viewField: GR_FIELD_PHI_GRAV,
        liveModified: false,
        selected: -1,
    };
}
