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
import { FieldPass } from './render/fieldPass';
import { Overlay2D, createTrails, type Trails } from './render/overlay2d';
import { computeView } from './sim/fieldViews';
import { fullReport, stabilityProbe } from './dev/diagnostics';
import { STEPS_PER_FRAME } from './sim/config';
import { applyScenario } from './sim/build';
import { emptyScenario, validate, type Scenario } from './sim/scenario';
import { readFromHash, writeToHash, downloadScenario } from './sim/serialization';

const DEFAULT_SCENARIO = 'pic_binary';

async function main(): Promise<void> {
    const canvas = document.getElementById('view') as HTMLCanvasElement;
    const overlayCanvas = document.getElementById('overlay') as HTMLCanvasElement;

    const state = createState('pic_binary');

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

    function setupGrid(g: Scenario['grid']): void {
        world = new World(core, g.W, g.H, g.dx, g.cEff, g.cfl);
        fieldPass = new FieldPass(gl!, g.W, g.H);
        overlay = new Overlay2D(overlayCanvas, g.W, g.H);
        trails = createTrails();
    }

    function loadScenario(scn: Scenario): void {
        controls.setStatus(`building ${scn.name}…`);
        if (!world || world.W !== scn.grid.W || world.H !== scn.grid.H) setupGrid(scn.grid);
        applyScenario(world, scn);
        current = scn;
        state.scenario = scn.name;
        state.viewField = scn.view.field;
        state.showTrails = scn.view.showTrails;
        state.showVelocity = scn.view.showVelocity;
        state.paused = true;
        state.liveModified = false;
        state.selected = -1;
        trails.reset(world.particleCount());
        controls.setPaused(true);
        controls.setField(state.viewField);
        controls.setTrails(state.showTrails);
        controls.setRadiation(scn.global.noRadiation);
        inspector.renderEditors();
        writeToHash(scn);
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

    /* Load a scenario by name from the JSON library (public/scenes/<name>.json).
     * The JSON files are the single source of truth -- there is no in-code copy. */
    async function buildByName(name: string): Promise<void> {
        try {
            const res = await fetch(`/scenes/${name}.json`, { cache: 'no-cache' });
            if (!res.ok) throw new Error(`HTTP ${res.status}`);
            loadScenario(validate(await res.json()));
            controls.setScenario(name);
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
        onFieldChange: (i) => { state.viewField = i; },
        onToggleTrails: () => { state.showTrails = !state.showTrails; controls.setTrails(state.showTrails); },
        onUndo: () => undo(),
        onToggleRadiation: () => applyEdit((sc) => { sc.global.noRadiation = !sc.global.noRadiation; }),
        onCopyLink: () => {
            writeToHash(current);
            void navigator.clipboard?.writeText(location.href);
            controls.setStatus('copied shareable URL to clipboard');
        },
        onSaveJson: () => downloadScenario(current),
        onProbe: () => {
            controls.setStatus('running stability probe…');
            requestAnimationFrame(() => {
                console.log(stabilityProbe(world, 20, 200));
                controls.setStatus('probe complete — see console');
            });
        },
    });
    controls.setTrails(state.showTrails);

    controls.setStatus('loading WASM…');
    core = await loadCore();

    const inspector = new Inspector({
        getScenario: () => current,
        getState: () => state,
        getWorld: () => world,
        onEdit: applyEdit,
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
    {
        let scenes: { value: string; label: string }[] = [];
        try {
            const res = await fetch('/scenes/index.json', { cache: 'no-cache' });
            if (res.ok) {
                const m = await res.json() as { scenes?: { file: string; label?: string }[] };
                scenes = (m.scenes ?? []).map((s) => ({ value: s.file, label: s.label ?? s.file }));
            }
        } catch (e) { console.warn('scene manifest fetch failed:', (e as Error).message); }
        if (scenes.length === 0) scenes = [{ value: DEFAULT_SCENARIO, label: DEFAULT_SCENARIO }];
        controls.setScenarios(scenes, DEFAULT_SCENARIO);
    }

    /* Initial scenario: URL hash if present, else the default library file. */
    await new Promise<void>((r) => requestAnimationFrame(() => r()));
    const hashScn = readFromHash();
    if (hashScn) loadScenario(hashScn);
    else await buildByName(DEFAULT_SCENARIO);

    /* Debug bridge (dev only). */
    if ((import.meta as unknown as { env?: { DEV?: boolean } }).env?.DEV) {
        const { startBridge } = await import('./dev/bridge');
        startBridge({
            getWorld: () => world,
            state, canvas, overlayCanvas,
            setStatus: controls.setStatus,
            setPaused: controls.setPaused,
            reset: () => loadScenario(current),
            buildByName,
        });
    }

    function frame(): void {
        const advanced = !state.paused || state.singleStep;
        if (advanced) {
            world.step(state.singleStep ? 1 : STEPS_PER_FRAME);
            state.singleStep = false;
        }

        const data = computeView(world, state.viewField);
        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clearColor(0, 0, 0, 1);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
        fieldPass.render(data);

        const parts = world.particles();
        if (advanced && state.showTrails) trails.push(parts);
        overlay.render(parts, trails,
            { showTrails: state.showTrails, showVelocity: state.showVelocity, selected: state.selected });

        inspector.updateLive();
        controls.setStatus(
            `${state.scenario}  step ${world.stepCount().toString().padStart(6)}  ` +
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
