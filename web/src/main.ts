/* GRlite web frontend — entry point.
 *
 * JSON-driven: a Scenario (sim/scenario.ts) is the source of truth.  The sim
 * is built from it via applyScenario (sim/build.ts); the URL hash carries
 * the canonical scenario so a link reproduces a run exactly. */

import { loadCore, type GRliteCore } from './wasm/binding';
import { World } from './sim/world';
import { createState, selectedOf } from './state';
import { wireControls } from './ui/controls';
import { Inspector } from './ui/inspector';
import { initJsonModal } from './ui/jsonModal';
import { FieldPass } from './render/fieldPass';
import { Overlay2D, createTrails, VEL_PX_PER_C, type Trails } from './render/overlay2d';
import { computeView, computeVectorField, colormapEnabled } from './sim/fieldViews';
import { fullReport, stabilityProbe } from './dev/diagnostics';
import { STEPS_PER_FRAME } from './sim/config';
import { applyScenario } from './sim/build';
import { emptyScenario, defaultParticle, defaultBackground, type Scenario } from './sim/scenario';
import { readFromHash, writeToHash, downloadScenario, parseScenarioText, toJSON } from './sim/serialization';
import { listLocalScenes, getLocalSceneText, saveLocalScene, deleteLocalScene } from './sim/localScenes';

/* Dropdown value prefix marking a user-saved (localStorage) scene. */
const LOCAL_PREFIX = 'local:';

const DEFAULT_SCENARIO = 'pic_binary';

async function main(): Promise<void> {
    const canvas = document.getElementById('view') as HTMLCanvasElement;
    const overlayCanvas = document.getElementById('overlay') as HTMLCanvasElement;

    const state = createState();

    const gl = canvas.getContext('webgl2',
        { antialias: true, premultipliedAlpha: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL2 not available');

    /* Grid-dependent objects (recreated when a scenario's grid differs). */
    let world!: World;
    let fieldPass!: FieldPass;
    let overlay!: Overlay2D;
    let trails!: Trails;
    let core!: GRliteCore;
    let current: Scenario = emptyScenario('—');   // placeholder until first load
    /* Library file key of the loaded scene ('' = custom / from a shared hash).
     * Persisted in the URL hash (f=) so a reload re-selects the right dropdown
     * entry instead of leaving it stale. */
    let currentFile = '';
    /* Library scene manifest (value=file key, label=dropdown text); filled once
     * the manifest is fetched.  Declared up here so the load helpers can map a
     * file key to a dropdown entry. */
    let scenes: { value: string; label: string }[] = [];
    /* Scenario JSON modal (assigned after controls/inspector exist; the toolbar
     * handler calls jsonModal.open() only at click time). */
    let jsonModal: ReturnType<typeof initJsonModal>;
    /* Force the next frame to recompute the field view even while paused -- set
     * on (re)load/reset so the screen reflects the rebuilt (e.g. zeroed) field
     * immediately instead of showing the cached pre-reset colormap. */
    let viewDirty = false;

    function setupGrid(g: Scenario['grid']): void {
        world = new World(core, g.W, g.H, g.dx, g.cEff, g.cfl);
        fieldPass = new FieldPass(gl!, g.W, g.H);
        overlay = new Overlay2D(overlayCanvas, g.W, g.H);
        trails = createTrails();
    }

    /* Push the canonical view settings (current.view) into the toolbar widgets.
     * current.view is the single source of truth; this just reflects it in the UI. */
    function syncViewControls(): void {
        controls.setField(current.view.field);
        controls.setSource(current.view.source);
        controls.setVectors(current.view.vectorField);
        controls.setVectorSpacing(current.view.vectorSpacing);
    }

    function loadScenario(scn: Scenario): void {
        controls.setStatus(`building ${scn.name}…`);
        if (!world || world.W !== scn.grid.W || world.H !== scn.grid.H) setupGrid(scn.grid);
        applyScenario(world, scn);
        current = scn;
        state.paused = true;
        state.liveModified = false;
        state.selection = null;
        trails.reset(world.particleCount());
        controls.setPaused(true);
        viewDirty = true;                 // redraw the field view for the new state
        syncViewControls();
        inspector.renderEditors();
        writeToHash(scn, currentFile || null);
        console.log(`[${scn.name}] built\n` + fullReport(world));
        controls.setStatus(`paused on ${scn.name} (press resume)`);
    }

    /* Edit model: snapshot for undo, mutate the canonical scenario, rebuild
     * from it (deterministic edit -> restart).  Undo pops + rebuilds. */
    const undoStack: Scenario[] = [];
    function applyEdit(mutate: (s: Scenario) => void): void {
        undoStack.push(structuredClone(current));
        if (undoStack.length > 50) undoStack.shift();
        /* loadScenario clears the selection (a fresh scenario has a different
         * particle/body set); for an in-place edit the set is unchanged, so
         * restore the selection if it's still valid -- otherwise editing a field
         * would deselect the object being edited. */
        const keep = state.selection;
        mutate(current);
        loadScenario(current);
        if (keep) {
            const n = keep.kind === 'particle' ? current.particles.length : current.background.length;
            if (keep.index >= 0 && keep.index < n) {
                state.selection = keep;
                inspector.renderEditors();
            }
        }
    }
    function undo(): void {
        const prev = undoStack.pop();
        if (!prev) { controls.setStatus('nothing to undo'); return; }
        loadScenario(prev);
    }

    /* View-only edit: mutate the canonical scenario but DON'T rebuild the sim,
     * so per-particle ticks/clock and tick Δ can be toggled mid-run.  The hash
     * is updated so the visualization choice is still shareable. */
    function applyViewEdit(mutate: (s: Scenario) => void): void {
        mutate(current);
        writeToHash(current, currentFile || null);
        inspector.renderEditors();
    }

    /* Is `file` a value that exists as a dropdown option (a library file key or
     * a saved local scene)?  Used so the dropdown can re-select what loaded. */
    function isDropdownValue(file: string): boolean {
        if (file.startsWith(LOCAL_PREFIX)) return listLocalScenes().includes(file.slice(LOCAL_PREFIX.length));
        return scenes.some((s) => s.value === file);
    }

    /* Rebuild the scenario dropdown: demo library + the user's saved scenes. */
    function refreshScenarioList(selected?: string): void {
        const items = [
            ...scenes.map((s) => ({ value: s.value, label: s.label, group: 'library' })),
            ...listLocalScenes().map((n) => ({ value: LOCAL_PREFIX + n, label: n, group: 'saved' })),
        ];
        controls.setScenarios(items, selected);
    }

    /* Apply an already-validated scenario: tag it with its dropdown value (library
     * file or local:name) or mark it custom/"(shared)", then build it. */
    function applyLoadedScenario(scn: Scenario, file: string | null): void {
        currentFile = (file && isDropdownValue(file)) ? file : '';
        loadScenario(scn);
        if (currentFile) controls.setScenario(currentFile);
        else controls.setCustomScenario(scn.name);
    }

    /* Load one of the user's saved scenes through the shared validated pipeline. */
    function loadLocalScene(name: string): void {
        const text = getLocalSceneText(name);
        if (text === null) { controls.setStatus(`no saved scene "${name}"`); return; }
        const errs = loadScenarioFromText(text, LOCAL_PREFIX + name);
        if (errs.length) controls.setStatus(`saved scene "${name}" failed validation: ${errs[0]}`);
    }

    /* THE one validated load path: a JSON string -> parseScenarioText (syntax +
     * format + value ranges) -> applyLoadedScenario.  Used by the library
     * fetch, the URL hash, and the paste modal alike.  Returns [] on success or
     * the list of problems (caller decides how to surface them). */
    function loadScenarioFromText(text: string, file: string | null): string[] {
        const res = parseScenarioText(text);
        if ('errors' in res) return res.errors;
        applyLoadedScenario(res.scenario, file);
        return [];
    }

    /* Load a scenario by name from the JSON library (public/scenes/<name>.json).
     * The JSON files are the single source of truth -- there is no in-code copy. */
    async function buildByName(name: string): Promise<void> {
        try {
            const res = await fetch(`/scenes/${name}.json`, { cache: 'no-cache' });
            if (!res.ok) throw new Error(`HTTP ${res.status}`);
            const errs = loadScenarioFromText(await res.text(), name);
            if (errs.length) {
                const msg = `scene "${name}" failed validation: ${errs[0]}`;
                console.error(msg, errs);
                controls.setStatus(msg);
            }
        } catch (e) {
            const msg = `failed to load /scenes/${name}.json: ${(e as Error).message}`;
            console.error(msg);
            controls.setStatus(msg);
        }
    }

    const controls = wireControls({
        onTogglePause: () => { state.paused = !state.paused; controls.setPaused(state.paused); },
        onStep: () => { state.singleStep = true; },
        onReset: () => loadScenario(current),
        onScenarioChange: (value) => {
            if (value === '__custom__') return;                  // synthetic "(shared)" entry
            if (value.startsWith(LOCAL_PREFIX)) loadLocalScene(value.slice(LOCAL_PREFIX.length));
            else void buildByName(value);
        },
        onFieldChange: (i) => { current.view.field = i; },
        onSourceChange: (sourceName) => { current.view.source = sourceName; },
        onVectorsChange: (i) => { current.view.vectorField = i; },
        onVectorSpacingChange: (n) => { current.view.vectorSpacing = n; },
        onUndo: () => undo(),
        /* Add a fresh particle at the grid center and select it for editing.
         * applyEdit rebuilds (and resets selection), so set selection + re-render
         * afterwards.  Undo removes it. */
        onAddParticle: () => {
            applyEdit((sc) => { sc.particles.push(defaultParticle(sc.grid.W, sc.grid.H)); });
            state.selection = { kind: 'particle', index: current.particles.length - 1 };
            controls.showParticlePanel();
            inspector.renderEditors();
            controls.setStatus('added particle — edit it in the selection inspector');
        },
        /* Add a fixed background body at the grid center (GM 0.01 so it's visible)
         * and select it for editing.  Undoable. */
        onAddBody: () => {
            applyEdit((sc) => { sc.background.push({ ...defaultBackground(sc.grid.W, sc.grid.H), GM: 0.01 }); });
            state.selection = { kind: 'body', index: current.background.length - 1 };
            controls.showParticlePanel();
            inspector.renderEditors();
            controls.setStatus('added background body — edit it in the selection inspector');
        },
        onCopyLink: () => {
            writeToHash(current, currentFile || null);
            void navigator.clipboard?.writeText(location.href);
            controls.setStatus('copied shareable URL to clipboard');
        },
        onSaveJson: () => downloadScenario(current),
        onOpenJson: () => jsonModal.open(),
        /* Save the current scenario to this browser under a name; it then appears
         * under "saved" in the dropdown.  The saved JSON carries the chosen name. */
        onSaveScene: () => {
            const name = (window.prompt('Save scenario as:', current.name) || '').trim();
            if (!name) return;
            if (name.startsWith(LOCAL_PREFIX)) { controls.setSceneMsg(`name cannot start with "${LOCAL_PREFIX}"`); return; }
            current.name = name;
            if (!saveLocalScene(name, toJSON(current))) {
                controls.setSceneMsg('save failed (localStorage unavailable or full)');
                return;
            }
            currentFile = LOCAL_PREFIX + name;
            writeToHash(current, currentFile);
            refreshScenarioList(currentFile);
            inspector.renderEditors();
            controls.setSceneMsg(`saved scene "${name}"`, true);
        },
        /* Delete the saved scene highlighted in the modal list (no effect on the
         * running sim, even if it's the one currently loaded). */
        onDeleteScene: (value) => {
            if (!value.startsWith(LOCAL_PREFIX)) { controls.setSceneMsg('select a saved scene to delete'); return; }
            const name = value.slice(LOCAL_PREFIX.length);
            if (!window.confirm(`Delete saved scene "${name}"? This cannot be undone.`)) return;
            deleteLocalScene(name);
            if (currentFile === value) {
                /* Deleted the loaded scene: keep it running, just drop the link. */
                currentFile = '';
                refreshScenarioList();
                controls.setCustomScenario(current.name);
            } else {
                refreshScenarioList(currentFile || undefined);
            }
            controls.setSceneMsg(`deleted scene "${name}"`, true);
        },
    });
    syncViewControls();

    /* Frozen-particle field-stability probe -- a dev diagnostic, not a sandbox
     * feature.  No toolbar button; run `grProbe()` (optionally grProbe(batches,
     * perBatch)) from the browser console.  Logs Φ_g's swing over time. */
    (window as unknown as { grProbe: (b?: number, p?: number) => void }).grProbe =
        (b = 20, p = 200) => { console.log(stabilityProbe(world, b, p)); };

    controls.setStatus('loading WASM…');
    core = await loadCore();

    const inspector = new Inspector({
        getScenario: () => current,
        getState: () => state,
        getWorld: () => world,
        onEdit: applyEdit,
        onViewEdit: applyViewEdit,
        /* Remove the particle; applyEdit rebuilds and clears the selection, so
         * the panel falls back to "click a particle".  Undoable. */
        onDeleteParticle: (index) => {
            if (index < 0 || index >= current.particles.length) return;
            applyEdit((sc) => { sc.particles.splice(index, 1); });
            state.selection = null;       // don't silently jump to a neighbor
            inspector.renderEditors();
            controls.setStatus('deleted particle (undo to restore)');
        },
        onDeleteBody: (index) => {
            if (index < 0 || index >= current.background.length) return;
            applyEdit((sc) => { sc.background.splice(index, 1); });
            state.selection = null;
            inspector.renderEditors();
            controls.setStatus('deleted background body (undo to restore)');
        },
        onSelectBody: (index) => {
            state.selection = { kind: 'body', index };
            controls.showParticlePanel();
            inspector.renderEditors();
        },
    });

    /* Scenario JSON import/export modal.  A pasted scenario is a custom one (no
     * library file), so clear currentFile and reflect it as "(shared)". */
    jsonModal = initJsonModal({
        getScenario: () => current,
        /* Pasted text goes through the same validated pipeline; a pasted scenario
         * has no library file (custom).  Returns problems for the modal to show. */
        onLoadText: (text) => loadScenarioFromText(text, null),
    });

    /* Canvas interaction: click selects the nearest particle; once a particle is
     * selected you can drag it to reposition (grab the marker), or Shift-drag its
     * velocity arrow tip to set its initial velocity.  Live feedback is always
     * written straight to the sim (setParticleFields).  The drop's meaning
     * depends on whether the sim has run yet:
     *   - at the initial state (step 0): bake it into the scenario via applyEdit
     *     -- a deterministic restart, undoable, saved/shared (you're editing the
     *     scenario's initial conditions).
     *   - after the sim has advanced: leave it as a live-only nudge to the
     *     physics-engine state; a reset reverts to the scenario. */
    const GRAB_TOL = 16;                       // px radius to grab a handle
    let drag: { mode: 'pos' | 'vel'; index: number; moved: boolean } | null = null;

    const pointerPx = (ev: PointerEvent): { px: number; py: number } => {
        const r = overlayCanvas.getBoundingClientRect();
        return {
            px: (ev.clientX - r.left) * (overlayCanvas.width / r.width),
            py: (ev.clientY - r.top) * (overlayCanvas.height / r.height),
        };
    };
    const cw = () => overlayCanvas.width, ch = () => overlayCanvas.height;
    /* sim<->canvas mapping, matching Overlay2D.sx/sy (sim y is +up). */
    const toCanvas = (x: number, y: number) => ({ sx: (x / world.W) * cw(), sy: ch() - (y / world.H) * ch() });
    const toWorld = (px: number, py: number) => ({ x: (px / cw()) * world.W, y: ((ch() - py) / ch()) * world.H });
    const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));

    /* Canvas pixel of a particle's velocity-arrow tip (mirrors Overlay2D). */
    const velTipPx = (p: ReturnType<typeof world.particle>) => {
        const { sx, sy } = toCanvas(p.x, p.y);
        const g = p.mass > 0 ? Math.sqrt(1 + (p.px * p.px + p.py * p.py) / (p.mass * p.mass)) : 1;
        const vx = p.mass > 0 ? p.px / (g * p.mass) : 0;
        const vy = p.mass > 0 ? p.py / (g * p.mass) : 0;
        return { sx: sx + vx * VEL_PX_PER_C, sy: sy - vy * VEL_PX_PER_C };
    };

    /* Nearest background body to a canvas pixel (hit-test its circle), or -1.
     * Tolerance is the larger of the marker grab radius and the drawn circle. */
    const pickBody = (px: number, py: number): number => {
        let best = -1, bestD = Infinity;
        current.background.forEach((b, i) => {
            const c = toCanvas(b.x, b.y);
            const rPx = Math.max(6, b.epsilon * (cw() / world.W));
            const d = Math.hypot(px - c.sx, py - c.sy);
            if (d <= Math.max(GRAB_TOL, rPx + 6) && d < bestD) { bestD = d; best = i; }
        });
        return best;
    };

    /* Write live drag values straight to the sim for instant feedback. */
    const applyDrag = (px: number, py: number): void => {
        if (!drag) return;
        const i = drag.index;
        const p = world.particle(i);
        if (drag.mode === 'pos') {
            const w = toWorld(px, py);
            world.setParticleFields(i, { x: clamp(w.x, 0, world.W), y: clamp(w.y, 0, world.H) });
        } else {
            const { sx, sy } = toCanvas(p.x, p.y);
            let vx = (px - sx) / VEL_PX_PER_C;
            let vy = (sy - py) / VEL_PX_PER_C;          // canvas +down -> sim +up
            const cMax = 0.99 * current.grid.cEff;
            const sp = Math.hypot(vx, vy);
            if (sp > cMax) { vx *= cMax / sp; vy *= cMax / sp; }
            const c = current.grid.cEff;
            const gamma = 1 / Math.sqrt(Math.max(1e-9, 1 - (vx * vx + vy * vy) / (c * c)));
            world.setParticleFields(i, { px: gamma * p.mass * vx, py: gamma * p.mass * vy });
        }
        state.liveModified = true;
        drag.moved = true;
    };

    overlayCanvas.addEventListener('pointerdown', (ev) => {
        const { px, py } = pointerPx(ev);
        const sel = selectedOf(state, 'particle');
        if (sel >= 0 && sel < world.particleCount()) {
            const p = world.particle(sel);
            if (ev.shiftKey) {
                /* Shift => velocity: grab the arrow tip (or the marker when v=0). */
                const tip = velTipPx(p);
                const m = toCanvas(p.x, p.y);
                if (Math.hypot(px - tip.sx, py - tip.sy) <= GRAB_TOL ||
                    Math.hypot(px - m.sx, py - m.sy) <= GRAB_TOL) {
                    drag = { mode: 'vel', index: sel, moved: false };
                }
            } else {
                const m = toCanvas(p.x, p.y);
                if (Math.hypot(px - m.sx, py - m.sy) <= GRAB_TOL) drag = { mode: 'pos', index: sel, moved: false };
            }
        }
        if (drag) {
            overlayCanvas.setPointerCapture(ev.pointerId);
            overlayCanvas.style.cursor = 'grabbing';
            ev.preventDefault();
            applyDrag(px, py);
        }
    });

    overlayCanvas.addEventListener('pointermove', (ev) => {
        const { px, py } = pointerPx(ev);
        if (drag) { applyDrag(px, py); return; }
        /* Hover cursor: 'grab' over the selected particle's draggable handle
         * (marker, or arrow tip when Shift is held) to advertise drag-to-place. */
        const sel = selectedOf(state, 'particle');
        let over = false;
        if (sel >= 0 && sel < world.particleCount()) {
            const p = world.particle(sel);
            const h = ev.shiftKey ? velTipPx(p) : toCanvas(p.x, p.y);
            over = Math.hypot(px - h.sx, py - h.sy) <= GRAB_TOL;
        }
        overlayCanvas.style.cursor = over ? 'grab' : 'crosshair';
    });

    const endDrag = (ev: PointerEvent): void => {
        if (!drag) {
            /* A plain click (no drag started): select the nearest particle,
             * else the nearest background body, else clear the selection. */
            const { px, py } = pointerPx(ev);
            const pick = overlay.pick(px, py, world.particles());
            if (pick >= 0) state.selection = { kind: 'particle', index: pick };
            else {
                const bp = pickBody(px, py);
                state.selection = bp >= 0 ? { kind: 'body', index: bp } : null;
            }
            inspector.renderEditors();
            return;
        }
        const { index, mode, moved } = drag;
        drag = null;
        overlayCanvas.style.cursor = 'crosshair';
        if (overlayCanvas.hasPointerCapture(ev.pointerId)) overlayCanvas.releasePointerCapture(ev.pointerId);
        if (!moved) return;                       // a click on the handle: nothing to commit
        /* After the sim has advanced this is a live-only nudge (liveModified is
         * already set by applyDrag); a reset reverts it.  Only at the initial
         * state does the drag edit the scenario's initial conditions. */
        if (world.stepCount() > 0) return;
        /* Bake the live (dragged) value into the scenario; applyEdit rebuilds
         * deterministically and preserves the selection. */
        const p = world.particle(index);
        if (mode === 'pos') {
            applyEdit((sc) => { sc.particles[index].x = p.x; sc.particles[index].y = p.y; });
        } else {
            const g = p.mass > 0 ? Math.sqrt(1 + (p.px * p.px + p.py * p.py) / (p.mass * p.mass)) : 1;
            const vx = p.mass > 0 ? p.px / (g * p.mass) : 0;
            const vy = p.mass > 0 ? p.py / (g * p.mass) : 0;
            applyEdit((sc) => { sc.particles[index].vx = vx; sc.particles[index].vy = vy; });
        }
    };
    overlayCanvas.addEventListener('pointerup', endDrag);
    overlayCanvas.addEventListener('pointercancel', (ev) => {
        if (drag) { drag = null; overlayCanvas.style.cursor = 'crosshair'; }
        if (overlayCanvas.hasPointerCapture(ev.pointerId)) overlayCanvas.releasePointerCapture(ev.pointerId);
    });

    /* Populate the scenario dropdown from the library manifest (the JSON
     * files define what's available; no hard-coded option list). */
    try {
        const res = await fetch('/scenes/index.json', { cache: 'no-cache' });
        if (res.ok) {
            const m = await res.json() as { scenes?: { file: string; label?: string }[] };
            scenes = (m.scenes ?? []).map((s) => ({ value: s.file, label: s.label ?? s.file }));
        }
    } catch (e) { console.warn('scene manifest fetch failed:', (e as Error).message); }
    if (scenes.length === 0) scenes = [{ value: DEFAULT_SCENARIO, label: DEFAULT_SCENARIO }];
    refreshScenarioList(DEFAULT_SCENARIO);   // demo library + the user's saved scenes

    /* Initial scenario: URL hash if present (shareable/reproducible state), else
     * the default library file.  Keep the dropdown in sync with what actually
     * loads: re-select the library file if the hash carries one, otherwise show
     * a "(shared)" entry so the dropdown never disagrees with the view. */
    await new Promise<void>((r) => requestAnimationFrame(() => r()));
    const hashState = readFromHash();
    if (hashState) {
        /* Same validated pipeline as the paste modal: decoded JSON string ->
         * loadScenarioFromText.  On bad/invalid hash, warn and fall back. */
        const errs = loadScenarioFromText(hashState.json, hashState.file);
        if (errs.length) {
            console.warn('ignoring invalid scenario in URL hash:', errs);
            await buildByName(DEFAULT_SCENARIO);
        }
    } else {
        await buildByName(DEFAULT_SCENARIO);
    }

    /* Debug bridge (dev only). */
    if ((import.meta as unknown as { env?: { DEV?: boolean } }).env?.DEV) {
        const { startBridge } = await import('./dev/bridge');
        startBridge({
            getWorld: () => world,
            getScenario: () => current,
            state, canvas, overlayCanvas,
            setStatus: controls.setStatus,
            setPaused: controls.setPaused,
            setViewField: (i) => { current.view.field = i; controls.setField(i); },
            reset: () => loadScenario(current),
            buildByName,
        });
    }

    /* Field-view caches: the colormap + vector arrays only change when the sim
     * advances or the view selection changes, so we skip the W*H CPU passes on
     * static (paused) frames and just re-draw the cached data. */
    let cmCache: Float32Array | null = null;
    let vecCache: ReturnType<typeof computeVectorField> = null;
    let viewKey = '';

    /* Proper-time clock: instantaneous dτ/dt per particle, finite-differenced
     * from the previous advanced frame (engine-integrated τ vs coordinate t). */
    let clockPrevT = -1;
    let clockPrevTau: number[] = [];
    let clockRate: number[] = [];

    function frame(): void {
        const advanced = !state.paused || state.singleStep;
        if (advanced) {
            world.step(state.singleStep ? 1 : STEPS_PER_FRAME);
            state.singleStep = false;
        }

        const view = current.view;
        const key = `${view.field}|${view.source}|${view.vectorField}`;
        const recompute = advanced || viewDirty || key !== viewKey || cmCache === null;
        viewKey = key;
        viewDirty = false;

        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clearColor(0, 0, 0, 1);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
        if (colormapEnabled(view.field)) {
            if (recompute) cmCache = computeView(world, view.field, view.source);
            fieldPass.render(cmCache!);
        } else {
            cmCache = null;
        }

        if (recompute) {
            vecCache = view.vectorField > 0
                ? computeVectorField(world, view.vectorField, view.source)
                : null;
        }
        const vec = vecCache;

        const parts = world.particles();
        const tNow = world.time();
        /* Always record the track (cheap) so ticks/clock work even when the
         * trail LINE is hidden; the line draw itself is gated by showTrails. */
        if (advanced) trails.push(parts, tNow);

        /* Update dτ/dt estimates (only on a forward advance). */
        if (advanced) {
            const dT = tNow - clockPrevT;
            for (let i = 0; i < parts.length; i++) {
                const prev = clockPrevTau[i];
                if (clockPrevT >= 0 && dT > 1e-9 && prev !== undefined) {
                    clockRate[i] = (parts[i].properTime - prev) / dT;
                } else if (clockRate[i] === undefined) {
                    clockRate[i] = 1;
                }
                clockPrevTau[i] = parts[i].properTime;
            }
            clockPrevT = tNow;
        }

        /* Per-particle visualization flags (ticks / clock) from the scenario. */
        const display = parts.map((p) => ({
            ticks: current.particles[p.index]?.ticks ?? 'none',
            clock: current.particles[p.index]?.clock ?? false,
            rate: clockRate[p.index] ?? 1,
        }));

        overlay.render(parts, trails,
            { showTrails: view.showTrails, showVelocity: view.showVelocity,
              selectedParticle: selectedOf(state, 'particle'), selectedBody: selectedOf(state, 'body'),
              vectors: vec ? { data: vec, spacing: view.vectorSpacing } : null,
              display, tickInterval: view.tickInterval, time: tNow,
              forces: view.forceArrows,
              bodies: current.background.map((b) => ({ x: b.x, y: b.y, r: b.epsilon })) });

        inspector.updateLive();
        controls.setStatus(
            `${current.name}  step ${world.stepCount().toString().padStart(6)}  ` +
            `t=${world.time().toFixed(2)}  parts=${parts.length}  ` +
            `${state.liveModified ? '[LIVE-MODIFIED] ' : ''}${state.paused ? '[PAUSED]' : ''}`);
        requestAnimationFrame(frame);
    }
    frame();
}

main().catch((err) => {
    const s = document.getElementById('status');
    if (s) s.textContent = `error: ${err instanceof Error ? err.message : String(err)}`;
    console.error(err);
});
