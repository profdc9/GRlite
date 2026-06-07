/* Scenario JSON import/export modal.
 *
 * Shows the current scenario as editable JSON (copy to clipboard / save), and
 * lets the user paste a scenario in and load it.  Pasted text is fully vetted
 * (syntax + format + value ranges via parseScenarioText) before it is applied;
 * problems are shown in the modal and the load is refused. */

import { toJSON, parseScenarioText } from '../sim/serialization';
import type { Scenario } from '../sim/scenario';

export interface JsonModalHandlers {
    getScenario: () => Scenario;     // current scenario, to populate the textarea
    onLoad: (scn: Scenario) => void; // apply a successfully-parsed scenario
}

export interface JsonModal { open: () => void; }

export function initJsonModal(h: JsonModalHandlers): JsonModal {
    const modal = document.getElementById('jsonModal') as HTMLElement;
    const text = document.getElementById('jsonText') as HTMLTextAreaElement;
    const msg = document.getElementById('jsonMsg') as HTMLElement;

    const setMsg = (s: string, ok = false): void => {
        msg.textContent = s;
        msg.style.color = ok ? '#7c9' : '#e88';
    };
    const fill = (): void => { text.value = toJSON(h.getScenario()); };
    const close = (): void => modal.classList.add('hidden');
    const open = (): void => { fill(); setMsg(''); modal.classList.remove('hidden'); text.focus(); };

    document.getElementById('jsonClose')!.addEventListener('click', close);
    /* Click the backdrop (outside the box) to dismiss. */
    modal.addEventListener('click', (e) => { if (e.target === modal) close(); });
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape' && !modal.classList.contains('hidden')) close();
    });

    document.getElementById('jsonCopy')!.addEventListener('click', () => {
        void navigator.clipboard?.writeText(text.value)
            .then(() => setMsg('copied to clipboard', true))
            .catch(() => { text.select(); setMsg('press Ctrl+C to copy the selected text'); });
    });
    document.getElementById('jsonRevert')!.addEventListener('click', () => {
        fill(); setMsg('reverted to the current scenario', true);
    });
    document.getElementById('jsonLoad')!.addEventListener('click', () => {
        const res = parseScenarioText(text.value);
        if ('errors' in res) { setMsg(res.errors.join('\n')); return; }
        h.onLoad(res.scenario);
        setMsg('loaded', true);
        close();
    });

    return { open };
}
