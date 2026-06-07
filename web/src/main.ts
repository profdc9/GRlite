/* GRlite web frontend — entry point.
 *
 * JSON-driven: a Scenario (sim/scenario.ts) is the source of truth.  The sim
 * is built from it via applyScenario (sim/build.ts); the URL hash carries
 * the canonical scenario so a link reproduces a run exactly. */

import { loadCore, type GRliteCore } from './wasm/binding';
import { World } from './sim/world';
import { createState } from './state';
import { wireControls } from './ui/controls';
import { Inspector } from './ui/inspector';
import { initJsonModal } from './ui/jsonModal';
import { FieldPass } from './render/fieldPass';
import { Overlay2D, createTrails, type Trails } from './render/overlay2d';
import { computeView, computeVectorField, colormapEnabled } from './sim/fieldViews';
import { fullReport, stabilityProbe } from './dev/diagnostics';
import { STEPS_PER_FRAME } from './sim/config';
import { applyScenario } from './sim/build';
import { emptyScenario, type Scenario } from './sim/scenario';
import { readFromHash, writeToHash, downloadScenario, parseScenarioText } from './sim/serialization';

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
        controls.setTrails(current.view.showTrails);
        controls.setVelocity(current.view.showVelocity);
    }

    function loadScenario(scn: Scenario): void {
        controls.setStatus(`building ${scn.name}…`);
        if (!world || world.W !== scn.grid.W || world.H !== scn.grid.H) setupGrid(scn.grid);
        applyScenario(world, scn);
        current = scn;
        state.paused = true;
        state.liveModified = false;
        state.selected = -1;
        trails.reset(world.particleCount());
        controls.setPaused(true);
        viewDirty = true;                 // redraw the field view for the new state
        syncViewControls();
        controls.setRadiation(scn.global.noRadiation);
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
        mutate(current);
        loadScenario(current);
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

    /* Apply an already-validated scenario: tag it with its library file (or mark
     * it custom/"(shared)") and build it.  The single post-validation entry. */
    function applyLoadedScenario(scn: Scenario, file: string | null): void {
        currentFile = (file && scenes.some((s) => s.value === file)) ? file : '';
        loadScenario(scn);
        if (currentFile) controls.setScenario(currentFile);
        else controls.setCustomScenario(scn.name);
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
        onScenarioChange: (name) => { void buildByName(name); },
        onFieldChange: (i) => { current.view.field = i; },
        onSourceChange: (sourceName) => { current.view.source = sourceName; },
        onVectorsChange: (i) => { current.view.vectorField = i; },
        onVectorSpacingChange: (n) => { current.view.vectorSpacing = n; },
        onToggleTrails: () => { current.view.showTrails = !current.view.showTrails; controls.setTrails(current.view.showTrails); },
        onToggleVelocity: () => { current.view.showVelocity = !current.view.showVelocity; controls.setVelocity(current.view.showVelocity); },
        onUndo: () => undo(),
        onToggleRadiation: () => applyEdit((sc) => { sc.global.noRadiation = !sc.global.noRadiation; }),
        onCopyLink: () => {
            writeToHash(current, currentFile || null);
            void navigator.clipboard?.writeText(location.href);
            controls.setStatus('copied shareable URL to clipboard');
        },
        onSaveJson: () => downloadScenario(current),
        onOpenJson: () => jsonModal.open(),
        onProbe: () => {
            controls.setStatus('running stability probe…');
            requestAnimationFrame(() => {
                console.log(stabilityProbe(world, 20, 200));
                controls.setStatus('probe complete — see console');
            });
        },
    });
    syncViewControls();

    controls.setStatus('loading WASM…');
    core = await loadCore();

    const inspector = new Inspector({
        getScenario: () => current,
        getState: () => state,
        getWorld: () => world,
        onEdit: applyEdit,
        onViewEdit: applyViewEdit,
    });

    /* Scenario JSON import/export modal.  A pasted scenario is a custom one (no
     * library file), so clear currentFile and reflect it as "(shared)". */
    jsonModal = initJsonModal({
        getScenario: () => current,
        /* Pasted text goes through the same validated pipeline; a pasted scenario
         * has no library file (custom).  Returns problems for the modal to show. */
        onLoadText: (text) => loadScenarioFromText(text, null),
    });

    /* Click to select the nearest particle. */
    overlayCanvas.addEventListener('click', (ev) => {
        const r = overlayCanvas.getBoundingClientRect();
        const px = (ev.clientX - r.left) * (overlayCanvas.width / r.width);
        const py = (ev.clientY - r.top) * (overlayCanvas.height / r.height);
        state.selected = overlay.pick(px, py, world.particles());
        inspector.renderEditors();
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
    controls.setScenarios(scenes, DEFAULT_SCENARIO);

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
            { showTrails: view.showTrails, showVelocity: view.showVelocity, selected: state.selected,
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
