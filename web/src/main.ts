/* GRlite web frontend — entry point.
 *
 * Thin wiring: load the WASM core, build a World, set up the field pass +
 * Canvas2D overlay + inspector, wire the controls, run the frame loop.
 * Substance lives in the modules (wasm/, sim/, render/, dev/, ui/). */

import { loadCore } from './wasm/binding';
import { World } from './sim/world';
import { createState } from './state';
import { wireControls } from './ui/controls';
import { Inspector } from './ui/inspector';
import { FieldPass } from './render/fieldPass';
import { Overlay2D, createTrails } from './render/overlay2d';
import { computeView } from './sim/fieldViews';
import { fullReport, stabilityProbe } from './dev/diagnostics';
import { STEPS_PER_FRAME } from './sim/config';

async function main(): Promise<void> {
    const canvas = document.getElementById('view') as HTMLCanvasElement;
    const overlayCanvas = document.getElementById('overlay') as HTMLCanvasElement;
    const scenarioSel = document.getElementById('scenario') as HTMLSelectElement;

    const state = createState(scenarioSel?.value || 'pic_binary');

    const controls = wireControls({
        onTogglePause: () => { state.paused = !state.paused; controls.setPaused(state.paused); },
        onStep: () => { state.singleStep = true; },
        onReset: () => rebuild(state.scenario),
        onScenarioChange: (name) => { state.scenario = name; rebuild(name); },
        onFieldChange: (i) => { state.viewField = i; },
        onToggleTrails: () => { state.showTrails = !state.showTrails; controls.setTrails(state.showTrails); },
        onProbe: () => {
            controls.setStatus('running stability probe…');
            requestAnimationFrame(() => {
                console.log(stabilityProbe(world, 20, 200));
                controls.setStatus('probe complete — see console');
            });
        },
    });
    controls.setTrails(state.showTrails);
    controls.setField(state.viewField);

    controls.setStatus('loading WASM…');
    const core = await loadCore();
    const world = new World(core);

    const gl = canvas.getContext('webgl2',
        { antialias: true, premultipliedAlpha: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL2 not available');
    const fieldPass = new FieldPass(gl, world.W, world.H);
    const overlay = new Overlay2D(overlayCanvas, world.W, world.H);
    const inspector = new Inspector();
    const trails = createTrails();

    /* Click to select the nearest particle. */
    overlayCanvas.addEventListener('click', (ev) => {
        const r = overlayCanvas.getBoundingClientRect();
        const px = (ev.clientX - r.left) * (overlayCanvas.width / r.width);
        const py = (ev.clientY - r.top) * (overlayCanvas.height / r.height);
        state.selected = overlay.pick(px, py, world.particles());
    });

    function rebuild(name: string): void {
        controls.setStatus(`building ${name}…`);
        const rc = world.buildScenario(name);
        if (rc !== 0) { controls.setStatus(`failed to build ${name} (rc=${rc})`); return; }
        state.paused = true;
        state.liveModified = false;
        state.selected = -1;
        trails.reset(world.particleCount());
        controls.setPaused(true);
        console.log(`[${name}] built\n` + fullReport(world));
        controls.setStatus(`paused on ${name} (press resume)`);
    }

    controls.setStatus(`building ${state.scenario}…`);
    await new Promise<void>((r) => requestAnimationFrame(() => r()));
    rebuild(state.scenario);

    /* Debug bridge (dev only). */
    if ((import.meta as unknown as { env?: { DEV?: boolean } }).env?.DEV) {
        const { startBridge } = await import('./dev/bridge');
        startBridge({
            world, state, canvas, overlayCanvas,
            setStatus: controls.setStatus,
            setPaused: controls.setPaused,
            rebuild,
        });
    }

    function frame(): void {
        const advanced = !state.paused || state.singleStep;
        if (advanced) {
            world.step(state.singleStep ? 1 : STEPS_PER_FRAME);
            state.singleStep = false;
        }

        /* Field pass (selected view, derived on CPU where needed). */
        const data = computeView(world, state.viewField);
        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clearColor(0, 0, 0, 1);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
        fieldPass.render(data);

        /* Overlay pass (Canvas2D). */
        const parts = world.particles();
        if (advanced && state.showTrails) trails.push(parts);
        overlay.render(parts, trails,
            { showTrails: state.showTrails, showVelocity: state.showVelocity, selected: state.selected });

        inspector.update(world, state);
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
