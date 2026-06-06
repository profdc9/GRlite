/* Toolbar + selector wiring.  Pure DOM, no framework. */

import { FIELD_VIEWS } from '../sim/fieldViews';

export interface ControlHandlers {
    onTogglePause: () => void;
    onStep: () => void;
    onReset: () => void;
    onScenarioChange: (name: string) => void;
    onFieldChange: (index: number) => void;
    onToggleTrails: () => void;
    onUndo: () => void;
    onToggleRadiation: () => void;
    onCopyLink?: () => void;
    onSaveJson?: () => void;
    onProbe?: () => void;
}

export interface Controls {
    setPaused: (paused: boolean) => void;
    setStatus: (text: string) => void;
    setTrails: (on: boolean) => void;
    setField: (index: number) => void;
    setScenarios: (items: { value: string; label: string }[], selected?: string) => void;
    setScenario: (value: string) => void;
    setRadiation: (on: boolean) => void;
}

export function wireControls(h: ControlHandlers): Controls {
    const statusEl = document.getElementById('status') as HTMLDivElement;
    const resetBtn = document.getElementById('reset') as HTMLButtonElement;
    const pauseBtn = document.getElementById('pause') as HTMLButtonElement;
    const stepBtn = document.getElementById('stepFrame') as HTMLButtonElement;
    const scenarioSel = document.getElementById('scenario') as HTMLSelectElement;
    const fieldSel = document.getElementById('field') as HTMLSelectElement;
    const trailsBtn = document.getElementById('trails') as HTMLButtonElement;
    const probeBtn = document.getElementById('probe') as HTMLButtonElement | null;

    /* Populate field selector from FIELD_VIEWS. */
    fieldSel.innerHTML = '';
    FIELD_VIEWS.forEach((v, i) => {
        const opt = document.createElement('option');
        opt.value = String(i);
        opt.textContent = v.label;
        fieldSel.appendChild(opt);
    });

    pauseBtn.addEventListener('click', h.onTogglePause);
    stepBtn.addEventListener('click', h.onStep);
    resetBtn.addEventListener('click', h.onReset);
    scenarioSel.addEventListener('change', () => h.onScenarioChange(scenarioSel.value));
    fieldSel.addEventListener('change', () => h.onFieldChange(parseInt(fieldSel.value, 10)));
    trailsBtn.addEventListener('click', h.onToggleTrails);
    const undoBtn = document.getElementById('undo') as HTMLButtonElement;
    const radBtn = document.getElementById('radiation') as HTMLButtonElement;
    undoBtn.addEventListener('click', h.onUndo);
    radBtn.addEventListener('click', h.onToggleRadiation);
    /* Panel show/hide toggles (DOM-only). */
    const togGlobal = document.getElementById('toggleGlobal') as HTMLButtonElement;
    const togParticle = document.getElementById('toggleParticle') as HTMLButtonElement;
    const panelGlobal = document.getElementById('panelGlobal') as HTMLElement;
    const panelParticle = document.getElementById('panelParticle') as HTMLElement;
    togGlobal.addEventListener('click', () => panelGlobal.classList.toggle('hidden'));
    togParticle.addEventListener('click', () => panelParticle.classList.toggle('hidden'));
    const copyBtn = document.getElementById('copylink') as HTMLButtonElement | null;
    const saveBtn = document.getElementById('savejson') as HTMLButtonElement | null;
    if (copyBtn && h.onCopyLink) copyBtn.addEventListener('click', h.onCopyLink);
    if (saveBtn && h.onSaveJson) saveBtn.addEventListener('click', h.onSaveJson);
    if (probeBtn && h.onProbe) probeBtn.addEventListener('click', h.onProbe);

    return {
        setPaused: (paused) => { pauseBtn.textContent = paused ? 'resume' : 'pause'; },
        setStatus: (text) => { statusEl.textContent = text; },
        setTrails: (on) => { trailsBtn.textContent = `trails: ${on ? 'on' : 'off'}`; },
        setField: (index) => { fieldSel.value = String(index); },
        setRadiation: (on) => { radBtn.textContent = `radiation: ${on ? 'on' : 'off'}`; },
        setScenarios: (items, selected) => {
            scenarioSel.innerHTML = '';
            for (const it of items) {
                const o = document.createElement('option');
                o.value = it.value; o.textContent = it.label;
                scenarioSel.appendChild(o);
            }
            if (selected !== undefined) scenarioSel.value = selected;
        },
        setScenario: (value) => { scenarioSel.value = value; },
    };
}
