/* Wires the toolbar buttons + scenario selector to callbacks.  Pure DOM,
 * no framework.  The richer inspectors land in Phase 3. */

export interface ControlHandlers {
    onTogglePause: () => void;
    onStep: () => void;
    onReset: () => void;
    onScenarioChange: (name: string) => void;
    onProbe?: () => void;
}

export interface Controls {
    setPaused: (paused: boolean) => void;
    setStatus: (text: string) => void;
}

export function wireControls(h: ControlHandlers): Controls {
    const statusEl = document.getElementById('status') as HTMLDivElement;
    const resetBtn = document.getElementById('reset') as HTMLButtonElement;
    const pauseBtn = document.getElementById('pause') as HTMLButtonElement;
    const stepBtn = document.getElementById('stepFrame') as HTMLButtonElement;
    const scenarioSel = document.getElementById('scenario') as HTMLSelectElement;
    const probeBtn = document.getElementById('probe') as HTMLButtonElement | null;

    pauseBtn.addEventListener('click', h.onTogglePause);
    stepBtn.addEventListener('click', h.onStep);
    resetBtn.addEventListener('click', h.onReset);
    scenarioSel.addEventListener('change', () => h.onScenarioChange(scenarioSel.value));
    if (probeBtn && h.onProbe) probeBtn.addEventListener('click', h.onProbe);

    return {
        setPaused: (paused: boolean) => { pauseBtn.textContent = paused ? 'resume' : 'pause'; },
        setStatus: (text: string) => { statusEl.textContent = text; },
    };
}
