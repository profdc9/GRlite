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

export function encodeToHash(scn: Scenario): string {
    return 's=' + b64urlEncode(JSON.stringify(scn));
}

/* Read a scenario from the current URL hash, or null if absent/invalid. */
export function readFromHash(hash = location.hash): Scenario | null {
    const m = /[#&]s=([^&]+)/.exec(hash);
    if (!m) return null;
    try { return validate(JSON.parse(b64urlDecode(m[1]))); }
    catch (e) { console.warn('bad scenario in URL hash:', (e as Error).message); return null; }
}

/* Write the scenario into the URL hash without reloading or adding history. */
export function writeToHash(scn: Scenario): void {
    const h = '#' + encodeToHash(scn);
    history.replaceState(null, '', h);
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
