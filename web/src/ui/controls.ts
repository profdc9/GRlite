/* Toolbar + selector wiring.  Pure DOM, no framework. */

import { FIELD_VIEWS, VECTOR_FIELDS } from '../sim/fieldViews';
import type { FieldSourceName } from '../sim/scenario';

export interface ControlHandlers {
    onTogglePause: () => void;
    onStep: () => void;
    onReset: () => void;
    onScenarioChange: (name: string) => void;
    onFieldChange: (index: number) => void;
    onSourceChange: (source: FieldSourceName) => void;
    onVectorsChange: (index: number) => void;
    onVectorSpacingChange: (cells: number) => void;
    onToggleTrails: () => void;
    onToggleVelocity: () => void;
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
    setVelocity: (on: boolean) => void;
    setField: (index: number) => void;
    setSource: (source: FieldSourceName) => void;
    setVectors: (index: number) => void;
    setVectorSpacing: (cells: number) => void;
    setScenarios: (items: { value: string; label: string }[], selected?: string) => void;
    setScenario: (value: string) => void;
    /* Select a synthetic "(shared)" entry for a scenario with no library file
     * (loaded from a URL hash / edited), so the dropdown reflects what's loaded. */
    setCustomScenario: (label: string) => void;
    setRadiation: (on: boolean) => void;
}

export function wireControls(h: ControlHandlers): Controls {
    const statusEl = document.getElementById('status') as HTMLDivElement;
    const resetBtn = document.getElementById('reset') as HTMLButtonElement;
    const pauseBtn = document.getElementById('pause') as HTMLButtonElement;
    const stepBtn = document.getElementById('stepFrame') as HTMLButtonElement;
    const scenarioSel = document.getElementById('scenario') as HTMLSelectElement;
    const fieldSel = document.getElementById('field') as HTMLSelectElement;
    const sourceSel = document.getElementById('source') as HTMLSelectElement;
    const vectorsSel = document.getElementById('vectors') as HTMLSelectElement;
    const vecSpacingInp = document.getElementById('vecspacing') as HTMLInputElement;
    const trailsBtn = document.getElementById('trails') as HTMLButtonElement;
    const velocityBtn = document.getElementById('velocity') as HTMLButtonElement;
    const probeBtn = document.getElementById('probe') as HTMLButtonElement | null;

    /* Populate field selector from FIELD_VIEWS. */
    fieldSel.innerHTML = '';
    FIELD_VIEWS.forEach((v, i) => {
        const opt = document.createElement('option');
        opt.value = String(i);
        opt.textContent = v.label;
        fieldSel.appendChild(opt);
    });

    /* Populate vector-field selector from VECTOR_FIELDS (index 0 = none). */
    vectorsSel.innerHTML = '';
    VECTOR_FIELDS.forEach((v, i) => {
        const opt = document.createElement('option');
        opt.value = String(i);
        opt.textContent = v.label;
        vectorsSel.appendChild(opt);
    });

    pauseBtn.addEventListener('click', h.onTogglePause);
    stepBtn.addEventListener('click', h.onStep);
    resetBtn.addEventListener('click', h.onReset);
    scenarioSel.addEventListener('change', () => h.onScenarioChange(scenarioSel.value));
    fieldSel.addEventListener('change', () => h.onFieldChange(parseInt(fieldSel.value, 10)));
    sourceSel.addEventListener('change', () => h.onSourceChange(sourceSel.value as FieldSourceName));
    vectorsSel.addEventListener('change', () => h.onVectorsChange(parseInt(vectorsSel.value, 10)));
    vecSpacingInp.addEventListener('change', () => {
        const n = parseInt(vecSpacingInp.value, 10);
        if (Number.isFinite(n) && n >= 2) h.onVectorSpacingChange(n);
    });
    trailsBtn.addEventListener('click', h.onToggleTrails);
    velocityBtn.addEventListener('click', h.onToggleVelocity);
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
        setVelocity: (on) => { velocityBtn.textContent = `velocity: ${on ? 'on' : 'off'}`; },
        setField: (index) => { fieldSel.value = String(index); },
        setSource: (source) => { sourceSel.value = source; },
        setVectors: (index) => { vectorsSel.value = String(index); },
        setVectorSpacing: (cells) => { vecSpacingInp.value = String(cells); },
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
        setCustomScenario: (label) => {
            let opt = scenarioSel.querySelector('option[value="__custom__"]') as HTMLOptionElement | null;
            if (!opt) {
                opt = document.createElement('option');
                opt.value = '__custom__';
                scenarioSel.appendChild(opt);
            }
            opt.textContent = `(shared) ${label}`;
            scenarioSel.value = '__custom__';
        },
    };
}
