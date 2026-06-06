/* GRlite web frontend — entry point.
 *
 * Thin wiring: load the WASM core, build a World, set up the two render
 * passes, wire the controls, and run the frame loop.  All substance lives
 * in the modules (wasm/, sim/, render/, dev/, ui/). */

import { loadCore } from './wasm/binding';
import { World } from './sim/world';
import { createState } from './state';
import { wireControls } from './ui/controls';
import { FieldPass } from './render/fieldPass';
import { OverlayPass } from './render/overlayPass';
import { fullReport, stabilityProbe } from './dev/diagnostics';
import { STEPS_PER_FRAME } from './sim/config';

async function main(): Promise<void> {
    const canvas = document.getElementById('view') as HTMLCanvasElement;
    const scenarioSel = document.getElementById('scenario') as HTMLSelectElement;

    const state = createState(scenarioSel?.value || 'pic_binary');

    const controls = wireControls({
        onTogglePause: () => { state.paused = !state.paused; controls.setPaused(state.paused); },
        onStep: () => { state.singleStep = true; },
        onReset: () => { rebuild(state.scenario); },
        onScenarioChange: (name) => { state.scenario = name; rebuild(name); },
        onProbe: () => {
            controls.setStatus('running stability probe…');
            requestAnimationFrame(() => {
                console.log(stabilityProbe(world, 20, 200));
                controls.setStatus('probe complete — see console');
            });
        },
    });

    controls.setStatus('loading WASM…');
    const core = await loadCore();
    const world = new World(core);

    const gl = canvas.getContext('webgl2', { antialias: true, premultipliedAlpha: false });
    if (!gl) throw new Error('WebGL2 not available');
    const fieldPass = new FieldPass(gl, world.W, world.H);
    const overlay = new OverlayPass(gl, world.W, world.H,
        Math.min(canvas.width, canvas.height) * 0.018);

    const PARTICLE_MAX = 256;
    const xyBuf = new Float32Array(PARTICLE_MAX * 2);

    function rebuild(name: string): void {
        controls.setStatus(`building ${name}…`);
        const rc = world.buildScenario(name);
        if (rc !== 0) { controls.setStatus(`failed to build ${name} (rc=${rc})`); return; }
        state.paused = true;
        state.liveModified = false;
        controls.setPaused(true);
        console.log(`[${name}] built\n` + fullReport(world));
        controls.setStatus(`paused on ${name} (press resume)`);
    }

    /* Initial build. */
    controls.setStatus(`building ${state.scenario}…`);
    await new Promise<void>((r) => requestAnimationFrame(() => r()));
    rebuild(state.scenario);

    function frame(): void {
        if (!state.paused || state.singleStep) {
            world.step(state.singleStep ? 1 : STEPS_PER_FRAME);
            state.singleStep = false;
        }

        const data = world.fieldView(state.viewField);
        gl!.viewport(0, 0, canvas.width, canvas.height);
        gl!.clearColor(0, 0, 0, 1);
        gl!.clear(gl!.COLOR_BUFFER_BIT);
        fieldPass.render(data);

        const n = Math.min(world.particleCount(), PARTICLE_MAX);
        for (let i = 0; i < n; i++) {
            const v = world.particleView(i);
            xyBuf[2 * i] = v[0];
            xyBuf[2 * i + 1] = v[1];
        }
        overlay.render(xyBuf, n);

        controls.setStatus(
            `${state.scenario}  step ${world.stepCount().toString().padStart(6)}  ` +
            `t=${world.time().toFixed(2)}  parts=${n}  ` +
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
