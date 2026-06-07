#!/usr/bin/env node
/* grlite-bridge — a debug RPC bridge.
 *
 *   Claude  ──MCP (stdio)──  this server  ──WebSocket (127.0.0.1:8787)──  GRlite web app
 *
 * The server is long-lived (Claude Code spawns it once).  The web app is a
 * WebSocket *client* that connects/reconnects on each reload, so restarting
 * the app never restarts this bridge.  Each MCP tool forwards a request to
 * the connected app and returns the app's computed TEXT result.
 *
 * IMPORTANT: MCP speaks JSON-RPC over stdout.  Never write logs to stdout --
 * use stderr (console.error) only. */

import { McpServer } from '@modelcontextprotocol/sdk/server/mcp.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import { z } from 'zod';
import { WebSocketServer } from 'ws';
import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const PORT = 8787;
const HOST = '127.0.0.1';
const HERE = dirname(fileURLToPath(import.meta.url));

/* ---- WebSocket side: talk to the app ---------------------------------- */

let appSocket = null;
let nextId = 1;
const pending = new Map();   // id -> { resolve, reject, timer }

const wss = new WebSocketServer({ host: HOST, port: PORT });
wss.on('listening', () => console.error(`[bridge] WS listening on ws://${HOST}:${PORT}`));
wss.on('error', (e) => console.error('[bridge] WS server error:', e.message));
wss.on('connection', (sock) => {
    console.error('[bridge] app connected');
    appSocket = sock;
    sock.on('message', (buf) => {
        let msg;
        try { msg = JSON.parse(buf.toString()); } catch { return; }
        if (msg && msg.id != null && pending.has(msg.id)) {
            const p = pending.get(msg.id);
            pending.delete(msg.id);
            clearTimeout(p.timer);
            if (msg.ok) p.resolve(msg.result);
            else p.reject(new Error(msg.error || 'app error'));
        }
    });
    sock.on('close', () => { if (appSocket === sock) appSocket = null; console.error('[bridge] app disconnected'); });
    sock.on('error', (e) => console.error('[bridge] app socket error:', e.message));
});

function callApp(method, params = {}, timeoutMs = 8000) {
    return new Promise((resolve, reject) => {
        if (!appSocket || appSocket.readyState !== appSocket.OPEN) {
            reject(new Error(`No GRlite app connected on ws://${HOST}:${PORT}. ` +
                `Open the app under the Vite dev server (the bridge client auto-connects in dev).`));
            return;
        }
        const id = nextId++;
        const timer = setTimeout(() => {
            pending.delete(id);
            reject(new Error(`app did not respond to "${method}" within ${timeoutMs}ms`));
        }, timeoutMs);
        pending.set(id, { resolve, reject, timer });
        appSocket.send(JSON.stringify({ id, method, params }));
    });
}

/* ---- MCP side: tools for Claude --------------------------------------- */

const server = new McpServer({ name: 'grlite-bridge', version: '0.1.0' });

const textResult = (text) => ({ content: [{ type: 'text', text: String(text) }] });
const errResult = (e) => ({ content: [{ type: 'text', text: `error: ${e.message || e}` }], isError: true });

/* Register a tool that forwards `method` to the app and returns its text. */
function passthrough(name, description, shape, method, mapParams) {
    server.registerTool(name, { description, inputSchema: shape ?? {} },
        async (args) => {
            try {
                const params = mapParams ? mapParams(args) : (args ?? {});
                const res = await callApp(method, params);
                return textResult(typeof res === 'string' ? res : JSON.stringify(res, null, 2));
            } catch (e) { return errResult(e); }
        });
}

passthrough('grlite_status', 'Connection + sim status (step, time, paused, scenario, particle count).', {}, 'status');
passthrough('grlite_report', 'Full text report: step/time, total momentum, COM, separation, per-particle state, Phi_g field stats.', {}, 'report');
passthrough('grlite_particles', 'Per-particle text report (position, momentum, mass, charge, speed, proper time).', {}, 'particles');
passthrough('grlite_conservation', 'Total momentum, center of mass, and separation (symmetry/COM-drift check).', {}, 'conservation');
passthrough('grlite_field_stats',
    'Min/max (with locations), mean, and RMS of a field. field: 0=Phi_g,1=A_gx,2=A_gy,3=phi_em,4=A_x,5=A_y.',
    { field: z.number().int().min(0).max(5).default(0) }, 'field_stats');
passthrough('grlite_field_profile',
    'A 1D slice of a field as (index,value) pairs. axis row=fixed y, col=fixed x.',
    { field: z.number().int().min(0).max(5).default(0),
      axis: z.enum(['row', 'col']).default('row'),
      at: z.number().int().default(192),
      stride: z.number().int().min(1).default(16) }, 'field_profile');
passthrough('grlite_step', 'Advance the simulation by n steps and report new step/time.',
    { n: z.number().int().min(1).default(1) }, 'step');
passthrough('grlite_play', 'Resume the simulation (unpause).', {}, 'play');
passthrough('grlite_pause', 'Pause the simulation.', {}, 'pause');
passthrough('grlite_reset', 'Rebuild the current scenario from scratch (back to t=0).', {}, 'reset');
passthrough('grlite_build', 'Build a named scenario (e.g. pic_binary, pic_static).',
    { scenario: z.string() }, 'build');
passthrough('grlite_set_particle',
    'Live-edit a particle (marks the run live-modified / non-reproducible). Only provided fields change.',
    { index: z.number().int().min(0),
      mass: z.number().optional(), charge: z.number().optional(),
      x: z.number().optional(), y: z.number().optional(),
      px: z.number().optional(), py: z.number().optional() }, 'set_particle');
passthrough('grlite_background',
    'List the background compact bodies (index + x,y,GM,Q,Jz,gFactor,eps).', {}, 'background');
passthrough('grlite_set_background',
    'Live-edit background body i (marks the run live-modified). Only provided fields change.',
    { index: z.number().int().min(0),
      x: z.number().optional(), y: z.number().optional(),
      GM: z.number().optional(), Q: z.number().optional(), Jz: z.number().optional(),
      gFactor: z.number().optional(), eps: z.number().optional() }, 'set_background');
passthrough('grlite_add_background',
    'Live-add a background compact body (max 16); returns its index. Defaults: center, GM 0.01, eps 8.',
    { x: z.number().optional(), y: z.number().optional(),
      GM: z.number().optional(), Q: z.number().optional(), Jz: z.number().optional(),
      gFactor: z.number().optional(), eps: z.number().optional() }, 'add_background');
passthrough('grlite_set_global',
    'Live-edit global state: gEff, kE (live-modify sim), viewField (display only), paused.',
    { gEff: z.number().optional(), kE: z.number().optional(),
      viewField: z.number().int().min(0).max(5).optional(),
      paused: z.boolean().optional() }, 'set_global');

/* Screenshot: app returns a data URL; save it to a PNG and return the path. */
server.registerTool('grlite_screenshot',
    { description: 'Capture the canvas to a PNG file on disk and return its path.', inputSchema: {} },
    async () => {
        try {
            const dataUrl = await callApp('screenshot', {}, 15000);
            const b64 = String(dataUrl).replace(/^data:image\/png;base64,/, '');
            mkdirSync(join(HERE, 'shots'), { recursive: true });
            const path = join(HERE, 'shots', `shot-${Date.now()}.png`);
            writeFileSync(path, Buffer.from(b64, 'base64'));
            return textResult(`saved screenshot: ${path}`);
        } catch (e) { return errResult(e); }
    });

const transport = new StdioServerTransport();
await server.connect(transport);
console.error('[bridge] MCP server ready on stdio');
