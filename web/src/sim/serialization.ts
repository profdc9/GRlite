/* Scenario <-> URL hash (and JSON text).  A scenario fully determines a run,
 * so a URL hash is a shareable, reproducible bug report: open the link, get
 * the exact state.  Encoding is base64url of the JSON (scenarios are small,
 * well under URL limits). */

import { validate, type Scenario } from './scenario';

function b64urlEncode(str: string): string {
    /* UTF-8 safe: percent-encode then btoa. */
    const b64 = btoa(encodeURIComponent(str).replace(/%([0-9A-F]{2})/g,
        (_, h) => String.fromCharCode(parseInt(h, 16))));
    return b64.replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
}

function b64urlDecode(s: string): string {
    const b64 = s.replace(/-/g, '+').replace(/_/g, '/');
    const bin = atob(b64);
    let pct = '';
    for (let i = 0; i < bin.length; i++) {
        pct += '%' + ('00' + bin.charCodeAt(i).toString(16)).slice(-2);
    }
    return decodeURIComponent(pct);
}

export function toJSON(scn: Scenario): string {
    return JSON.stringify(scn, null, 2);
}

export function fromJSON(text: string): Scenario {
    return validate(JSON.parse(text));
}

/* The hash carries the full scenario (s=) plus, when the run came from a named
 * library scene, that scene's file key (f=) so a reload can re-select the right
 * dropdown entry instead of leaving it stale. */
export function encodeToHash(scn: Scenario, file?: string | null): string {
    const s = 's=' + b64urlEncode(JSON.stringify(scn));
    return file ? `${s}&f=${encodeURIComponent(file)}` : s;
}

/* Read the scenario (and its library file key, if any) from the current URL
 * hash, or null if absent/invalid. */
export function readFromHash(hash = location.hash): { scenario: Scenario; file: string | null } | null {
    const m = /[#&]s=([^&]+)/.exec(hash);
    if (!m) return null;
    try {
        const scenario = validate(JSON.parse(b64urlDecode(m[1])));
        const fm = /[#&]f=([^&]+)/.exec(hash);
        return { scenario, file: fm ? decodeURIComponent(fm[1]) : null };
    } catch (e) { console.warn('bad scenario in URL hash:', (e as Error).message); return null; }
}

/* Write the scenario into the URL hash without reloading or adding history. */
export function writeToHash(scn: Scenario, file?: string | null): void {
    history.replaceState(null, '', '#' + encodeToHash(scn, file));
}

/* Trigger a .json file download of the scenario. */
export function downloadScenario(scn: Scenario): void {
    const blob = new Blob([toJSON(scn)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${scn.name.replace(/\s+/g, '-') || 'scenario'}.json`;
    a.click();
    URL.revokeObjectURL(url);
}
