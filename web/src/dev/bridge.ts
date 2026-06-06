/* Debug-bridge client (dev only).
 *
 * Connects to the grlite-bridge MCP server's WebSocket (ws://127.0.0.1:8787)
 * and answers RPC requests by calling the World / diagnostics API, returning
 * computed TEXT so the agent on the other end reads numbers, not pixels.
 *
 * Gated to dev builds (import.meta.env.DEV).  Never bundled in production;
 * binds only to localhost via the server side.  Auto-reconnects so the
 * long-lived bridge survives app reloads. */

import type { World } from '../sim/world';
import type { AppState } from '../state';
import {
    fullReport, particleReport, conservationReport, fieldStats, fieldProfile,
} from './diagnostics';

export interface BridgeContext {
    world: World;
    state: AppState;
    canvas: HTMLCanvasElement;          // WebGL field canvas
    overlayCanvas?: HTMLCanvasElement;  // Canvas2D overlay (trails/arrows)
    setStatus: (s: string) => void;
    setPaused: (p: boolean) => void;
    rebuild: (scenario: string) => void;
}

/* Composite the WebGL field canvas and the Canvas2D overlay into one PNG. */
function compositeScreenshot(ctx: BridgeContext): string {
    const base = ctx.canvas;
    const tmp = document.createElement('canvas');
    tmp.width = base.width; tmp.height = base.height;
    const c = tmp.getContext('2d');
    if (!c) return base.toDataURL('image/png');
    c.drawImage(base, 0, 0);
    if (ctx.overlayCanvas) c.drawImage(ctx.overlayCanvas, 0, 0, tmp.width, tmp.height);
    return tmp.toDataURL('image/png');
}

type Handler = (p: any) => unknown;

function buildHandlers(ctx: BridgeContext): Record<string, Handler> {
    const { world, state } = ctx;
    return {
        status: () => ({
            connected: true,
            scenario: state.scenario,
            step: world.stepCount(),
            time: world.time(),
            paused: state.paused,
            liveModified: state.liveModified,
            particleCount: world.particleCount(),
            grid: `${world.W}x${world.H}`,
        }),
        report: () => fullReport(world),
        particles: () => particleReport(world),
        conservation: () => conservationReport(world),
        field_stats: (p) => fieldStats(world, p?.field ?? 0),
        field_profile: (p) =>
            fieldProfile(world, p?.field ?? 0, p?.axis ?? 'row', p?.at ?? (world.H >> 1), p?.stride ?? 16),
        step: (p) => { world.step(Math.max(1, p?.n ?? 1)); return `step=${world.stepCount()} t=${world.time().toFixed(3)}`; },
        play: () => { state.paused = false; ctx.setPaused(false); return 'playing'; },
        pause: () => { state.paused = true; ctx.setPaused(true); return 'paused'; },
        reset: () => { ctx.rebuild(state.scenario); return `rebuilt ${state.scenario}`; },
        build: (p) => {
            if (!p?.scenario) throw new Error('build requires { scenario }');
            state.scenario = p.scenario; ctx.rebuild(p.scenario);
            return `built ${p.scenario}`;
        },
        set_particle: (p) => {
            if (p?.index == null) throw new Error('set_particle requires { index }');
            world.setParticleFields(p.index, p);
            state.liveModified = true;
            return `particle ${p.index} updated (run is now live-modified)`;
        },
        set_global: (p) => {
            let live = false;
            if (p?.gEff !== undefined) { world.setGEff(p.gEff); live = true; }
            if (p?.kE !== undefined) { world.setKE(p.kE); live = true; }
            if (p?.viewField !== undefined) state.viewField = p.viewField;
            if (p?.paused !== undefined) { state.paused = !!p.paused; ctx.setPaused(state.paused); }
            if (live) state.liveModified = true;
            return `global updated${live ? ' (run is now live-modified)' : ''}`;
        },
        screenshot: () => compositeScreenshot(ctx),
    };
}

export function startBridge(ctx: BridgeContext, url = 'ws://127.0.0.1:8787'): void {
    const handlers = buildHandlers(ctx);
    let backoff = 1000;

    const connect = () => {
        let ws: WebSocket;
        try { ws = new WebSocket(url); }
        catch { setTimeout(connect, backoff); return; }

        ws.onopen = () => { backoff = 1000; console.info(`[bridge] connected to ${url}`); };
        ws.onclose = () => { setTimeout(connect, backoff); backoff = Math.min(backoff * 2, 15000); };
        ws.onerror = () => { /* onclose will handle reconnect */ };
        ws.onmessage = (ev) => {
            let req: { id: number; method: string; params?: any };
            try { req = JSON.parse(ev.data as string); } catch { return; }
            const h = handlers[req.method];
            if (!h) { ws.send(JSON.stringify({ id: req.id, ok: false, error: `unknown method "${req.method}"` })); return; }
            try {
                const result = h(req.params ?? {});
                ws.send(JSON.stringify({ id: req.id, ok: true, result }));
            } catch (e) {
                ws.send(JSON.stringify({ id: req.id, ok: false, error: (e as Error).message }));
            }
        };
    };
    connect();
}
