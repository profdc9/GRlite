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
    onUndo: () => void;
    onAddParticle?: () => void;
    onCopyLink?: () => void;
    onSaveJson?: () => void;
    onOpenJson?: () => void;
    onSaveScene?: () => void;
    /* Delete the scene identified by the highlighted listbox value (a "local:"
     * key for a saved scene; anything else is rejected by the handler). */
    onDeleteScene?: (value: string) => void;
}

/* A scenario dropdown entry; `group` puts it under an <optgroup>. */
export interface ScenarioItem { value: string; label: string; group?: string; }

export interface Controls {
    setPaused: (paused: boolean) => void;
    setStatus: (text: string) => void;
    /* Reveal the particle inspector (e.g. after adding a particle, so the new
     * one is immediately editable even if the panel was toggled off). */
    showParticlePanel: () => void;
    /* Feedback line inside the scenario modal (save/delete results), since the
     * page status line is hidden behind the modal backdrop. */
    setSceneMsg: (text: string, ok?: boolean) => void;
    setField: (index: number) => void;
    setSource: (source: FieldSourceName) => void;
    setVectors: (index: number) => void;
    setVectorSpacing: (cells: number) => void;
    setScenarios: (items: ScenarioItem[], selected?: string) => void;
    setScenario: (value: string) => void;
    /* Select a synthetic "(shared)" entry for a scenario with no library file
     * (loaded from a URL hash / edited), so the dropdown reflects what's loaded. */
    setCustomScenario: (label: string) => void;
}

export function wireControls(h: ControlHandlers): Controls {
    const statusEl = document.getElementById('status') as HTMLDivElement;
    const resetBtn = document.getElementById('reset') as HTMLButtonElement;
    const pauseBtn = document.getElementById('pause') as HTMLButtonElement;
    const stepBtn = document.getElementById('stepFrame') as HTMLButtonElement;
    /* Scenario picker: a "scenario list" button opens a modal holding a
     * scrollable listbox; selection is applied on Load, not on highlight. */
    const scenarioSel = document.getElementById('scenarioSelect') as HTMLSelectElement;
    const scenarioListBtn = document.getElementById('scenariolist') as HTMLButtonElement;
    const scenarioModal = document.getElementById('scenarioModal') as HTMLElement;
    const scenarioLoadBtn = document.getElementById('scenarioLoad') as HTMLButtonElement;
    const scenarioCancelBtn = document.getElementById('scenarioCancel') as HTMLButtonElement;
    const scenarioMsgEl = document.getElementById('scenarioMsg') as HTMLElement;
    const fieldSel = document.getElementById('field') as HTMLSelectElement;
    const sourceSel = document.getElementById('source') as HTMLSelectElement;
    const vectorsSel = document.getElementById('vectors') as HTMLSelectElement;
    const vecSpacingInp = document.getElementById('vecspacing') as HTMLInputElement;

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
    /* Modal open/close + load.  On open, snapshot the highlighted value so a
     * Cancel restores it (the highlight tracks what's actually loaded, kept in
     * sync by setScenario/setCustomScenario).  Load applies the highlight. */
    let scenarioReopenValue = '';
    const openScenarioModal = (): void => {
        scenarioReopenValue = scenarioSel.value;
        scenarioMsgEl.textContent = '';
        scenarioModal.classList.remove('hidden');
        scenarioSel.focus();
    };
    const closeScenarioModal = (): void => scenarioModal.classList.add('hidden');
    const cancelScenarioModal = (): void => { scenarioSel.value = scenarioReopenValue; closeScenarioModal(); };
    const loadScenarioModal = (): void => {
        if (!scenarioSel.value) return;
        closeScenarioModal();
        h.onScenarioChange(scenarioSel.value);
    };
    scenarioListBtn.addEventListener('click', openScenarioModal);
    scenarioLoadBtn.addEventListener('click', loadScenarioModal);
    scenarioCancelBtn.addEventListener('click', cancelScenarioModal);
    scenarioSel.addEventListener('dblclick', loadScenarioModal);
    scenarioModal.addEventListener('click', (e) => { if (e.target === scenarioModal) cancelScenarioModal(); });
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && !scenarioModal.classList.contains('hidden')) cancelScenarioModal();
    });
    fieldSel.addEventListener('change', () => h.onFieldChange(parseInt(fieldSel.value, 10)));
    sourceSel.addEventListener('change', () => h.onSourceChange(sourceSel.value as FieldSourceName));
    vectorsSel.addEventListener('change', () => h.onVectorsChange(parseInt(vectorsSel.value, 10)));
    vecSpacingInp.addEventListener('change', () => {
        const n = parseInt(vecSpacingInp.value, 10);
        if (Number.isFinite(n) && n >= 2) h.onVectorSpacingChange(n);
    });
    const undoBtn = document.getElementById('undo') as HTMLButtonElement;
    undoBtn.addEventListener('click', h.onUndo);
    /* Panel show/hide toggles (DOM-only). */
    const togGlobal = document.getElementById('toggleGlobal') as HTMLButtonElement;
    const togParticle = document.getElementById('toggleParticle') as HTMLButtonElement;
    const panelGlobal = document.getElementById('panelGlobal') as HTMLElement;
    const panelParticle = document.getElementById('panelParticle') as HTMLElement;
    togGlobal.addEventListener('click', () => panelGlobal.classList.toggle('hidden'));
    togParticle.addEventListener('click', () => panelParticle.classList.toggle('hidden'));
    const addParticleBtn = document.getElementById('addparticle') as HTMLButtonElement | null;
    if (addParticleBtn && h.onAddParticle) addParticleBtn.addEventListener('click', h.onAddParticle);
    const copyBtn = document.getElementById('copylink') as HTMLButtonElement | null;
    const saveBtn = document.getElementById('savejson') as HTMLButtonElement | null;
    const jsonBtn = document.getElementById('jsonbtn') as HTMLButtonElement | null;
    if (copyBtn && h.onCopyLink) copyBtn.addEventListener('click', h.onCopyLink);
    if (saveBtn && h.onSaveJson) saveBtn.addEventListener('click', h.onSaveJson);
    if (jsonBtn && h.onOpenJson) jsonBtn.addEventListener('click', h.onOpenJson);
    const saveSceneBtn = document.getElementById('savescene') as HTMLButtonElement | null;
    const delSceneBtn = document.getElementById('deletescene') as HTMLButtonElement | null;
    if (saveSceneBtn && h.onSaveScene) saveSceneBtn.addEventListener('click', h.onSaveScene);
    if (delSceneBtn && h.onDeleteScene) delSceneBtn.addEventListener('click', () => h.onDeleteScene!(scenarioSel.value));

    return {
        setPaused: (paused) => { pauseBtn.textContent = paused ? 'resume' : 'pause'; },
        setStatus: (text) => { statusEl.textContent = text; },
        showParticlePanel: () => { panelParticle.classList.remove('hidden'); },
        setSceneMsg: (text, ok = false) => { scenarioMsgEl.textContent = text; scenarioMsgEl.style.color = ok ? '#7c9' : '#e88'; },
        setField: (index) => { fieldSel.value = String(index); },
        setSource: (source) => { sourceSel.value = source; },
        setVectors: (index) => { vectorsSel.value = String(index); },
        setVectorSpacing: (cells) => { vecSpacingInp.value = String(cells); },
        setScenarios: (items, selected) => {
            scenarioSel.innerHTML = '';
            /* Group entries under <optgroup> by their `group` field (preserving
             * first-seen order); ungrouped entries go straight in. */
            const groups = new Map<string, ScenarioItem[]>();
            for (const it of items) {
                const g = it.group ?? '';
                if (!groups.has(g)) groups.set(g, []);
                groups.get(g)!.push(it);
            }
            const mkOpt = (it: ScenarioItem): HTMLOptionElement => {
                const o = document.createElement('option');
                o.value = it.value; o.textContent = it.label;
                return o;
            };
            for (const [g, list] of groups) {
                if (g) {
                    const og = document.createElement('optgroup');
                    og.label = g;
                    for (const it of list) og.appendChild(mkOpt(it));
                    scenarioSel.appendChild(og);
                } else {
                    for (const it of list) scenarioSel.appendChild(mkOpt(it));
                }
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
