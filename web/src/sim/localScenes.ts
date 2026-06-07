/* User's personal scenario library, persisted in localStorage.
 *
 * Stored as one JSON object { [name]: scenarioJsonText } under a single key.
 * Scenes are kept as JSON *text* so loading reuses the same validated pipeline
 * (parseScenarioText) as the library/URL/paste paths.  All access is wrapped --
 * localStorage may be unavailable (private mode) or full, in which case saves
 * silently no-op and lists come back empty. */

const KEY = 'grlite.localScenes';

function readAll(): Record<string, string> {
    try {
        const raw = localStorage.getItem(KEY);
        const obj = raw ? JSON.parse(raw) : {};
        return (obj && typeof obj === 'object') ? obj as Record<string, string> : {};
    } catch { return {}; }
}

function writeAll(m: Record<string, string>): boolean {
    try { localStorage.setItem(KEY, JSON.stringify(m)); return true; }
    catch (e) { console.warn('localStorage write failed:', (e as Error).message); return false; }
}

export function listLocalScenes(): string[] {
    return Object.keys(readAll()).sort((a, b) => a.localeCompare(b));
}

export function getLocalSceneText(name: string): string | null {
    const m = readAll();
    return Object.prototype.hasOwnProperty.call(m, name) ? m[name] : null;
}

export function saveLocalScene(name: string, jsonText: string): boolean {
    const m = readAll();
    m[name] = jsonText;
    return writeAll(m);
}

export function deleteLocalScene(name: string): void {
    const m = readAll();
    delete m[name];
    writeAll(m);
}
