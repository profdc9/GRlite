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
import type { Scenario } from '../sim/scenario';
import {
    fullReport, particleReport, conservationReport, fieldStats, fieldProfile,
} from './diagnostics';

export interface BridgeContext {
    getWorld: () => World;              // getter: the World may be recreated on grid change
    getScenario: () => Scenario;       // canonical scenario (name + view live here)
    state: AppState;
    canvas: HTMLCanvasElement;          // WebGL field canvas
    overlayCanvas?: HTMLCanvasElement;  // Canvas2D overlay (trails/arrows)
    setStatus: (s: string) => void;
    setPaused: (p: boolean) => void;
    setViewField: (i: number) => void; // set the colormap field (updates UI too)
    reset: () => void;                  // rebuild the current scenario
    buildByName: (name: string) => void | Promise<void>;
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

type Handler = (p: any) => unknown | Promise<unknown>;

function buildHandlers(ctx: BridgeContext): Record<string, Handler> {
    const { state } = ctx;
    const w = () => ctx.getWorld();
    return {
        status: () => ({
            connected: true,
            scenario: ctx.getScenario().name,
            step: w().stepCount(),
            time: w().time(),
            paused: state.paused,
            liveModified: state.liveModified,
            particleCount: w().particleCount(),
            grid: `${w().W}x${w().H}`,
        }),
        report: () => fullReport(w()),
        particles: () => particleReport(w()),
        conservation: () => conservationReport(w()),
        field_stats: (p) => fieldStats(w(), p?.field ?? 0),
        field_profile: (p) =>
            fieldProfile(w(), p?.field ?? 0, p?.axis ?? 'row', p?.at ?? (w().H >> 1), p?.stride ?? 16),
        step: (p) => { w().step(Math.max(1, p?.n ?? 1)); return `step=${w().stepCount()} t=${w().time().toFixed(3)}`; },
        play: () => { state.paused = false; ctx.setPaused(false); return 'playing'; },
        pause: () => { state.paused = true; ctx.setPaused(true); return 'paused'; },
        reset: () => { ctx.reset(); return `rebuilt ${ctx.getScenario().name}`; },
        build: async (p) => {
            if (!p?.scenario) throw new Error('build requires { scenario }');
            await ctx.buildByName(p.scenario);
            return `built ${p.scenario}`;
        },
        set_particle: (p) => {
            if (p?.index == null) throw new Error('set_particle requires { index }');
            w().setParticleFields(p.index, p);
            state.liveModified = true;
            return `particle ${p.index} updated (run is now live-modified)`;
        },
        /* List the background compact bodies (index + params). */
        background: () => {
            const n = w().backgroundCount();
            const bodies = [];
            for (let i = 0; i < n; i++) bodies.push({ index: i, ...w().backgroundBody(i) });
            return { count: n, bodies };
        },
        /* Live-edit background body i: { index, x?, y?, GM?, Q?, Jz?, gFactor?, eps? }. */
        set_background: (p) => {
            if (p?.index == null) throw new Error('set_background requires { index }');
            const cur = w().backgroundBody(p.index);
            if (!cur) throw new Error(`no background body ${p.index} (count=${w().backgroundCount()})`);
            w().setBackgroundBodyAt(p.index, {
                x: p.x ?? cur.x, y: p.y ?? cur.y, GM: p.GM ?? cur.GM, Q: p.Q ?? cur.Q,
                Jz: p.Jz ?? cur.Jz, gFactor: p.gFactor ?? cur.gFactor, eps: p.eps ?? cur.eps,
            });
            state.liveModified = true;
            return `background body ${p.index} updated (run is now live-modified)`;
        },
        /* Live-add a background body; returns its index or errors when full. */
        add_background: (p) => {
            const idx = w().addBackgroundBody({
                x: p?.x ?? (w().W >> 1), y: p?.y ?? (w().H >> 1), GM: p?.GM ?? 0.01,
                Q: p?.Q ?? 0, Jz: p?.Jz ?? 0, gFactor: p?.gFactor ?? 2, eps: p?.eps ?? 8,
            });
            if (idx < 0) throw new Error('background body list full (max 16)');
            state.liveModified = true;
            return `added background body ${idx} (run is now live-modified)`;
        },
        set_global: (p) => {
            let live = false;
            if (p?.gEff !== undefined) { w().setGEff(p.gEff); live = true; }
            if (p?.kE !== undefined) { w().setKE(p.kE); live = true; }
            if (p?.viewField !== undefined) ctx.setViewField(p.viewField);
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
            /* Supports both sync and async handlers. */
            Promise.resolve().then(() => h(req.params ?? {}))
                .then((result) => ws.send(JSON.stringify({ id: req.id, ok: true, result })))
                .catch((e) => ws.send(JSON.stringify({ id: req.id, ok: false, error: (e as Error).message })));
        };
    };
    connect();
}
