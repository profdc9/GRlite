/* End-to-end smoke test for the bridge.
 *
 *   MCP client (this test) ──stdio── server.mjs ──WS── mock app (this test)
 *
 * Spawns the real server via the MCP stdio client, connects a mock app that
 * answers requests, then calls tools and checks the relayed results. */

import { Client } from '@modelcontextprotocol/sdk/client/index.js';
import { StdioClientTransport } from '@modelcontextprotocol/sdk/client/stdio.js';
import { WebSocket } from 'ws';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const HERE = dirname(fileURLToPath(import.meta.url));
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
let failures = 0;
const check = (cond, msg) => { console.log(`${cond ? 'PASS' : 'FAIL'}: ${msg}`); if (!cond) failures++; };

/* 1. Spawn the server via the MCP stdio client. */
const transport = new StdioClientTransport({ command: 'node', args: [join(HERE, 'server.mjs')] });
const client = new Client({ name: 'grlite-bridge-test', version: '0.0.0' });
await client.connect(transport);
console.log('[test] MCP connected');

/* 2. Tools are listed. */
const { tools } = await client.listTools();
const names = tools.map((t) => t.name);
check(names.includes('grlite_status'), 'grlite_status tool registered');
check(names.includes('grlite_set_particle'), 'grlite_set_particle tool registered');
check(names.length >= 10, `>=10 tools registered (got ${names.length})`);

/* 3. With no app connected, a tool reports the disconnected state. */
const noApp = await client.callTool({ name: 'grlite_status', arguments: {} });
check(/No GRlite app connected/.test(noApp.content[0].text), 'reports "no app connected" before app joins');

/* 4. Connect a mock app over WS and answer requests. */
await sleep(300);
const app = new WebSocket('ws://127.0.0.1:8787');
await new Promise((res, rej) => { app.on('open', res); app.on('error', rej); });
app.on('message', (buf) => {
    const req = JSON.parse(buf.toString());
    const replies = {
        status: { connected: true, step: 42, scenario: 'pic_binary' },
        report: 'MOCK REPORT: step=42 P_total=(0,0)',
        set_particle: `particle ${req.params?.index} updated`,
        step: `step=43 t=0.707`,
    };
    const result = replies[req.method] ?? `mock:${req.method}`;
    app.send(JSON.stringify({ id: req.id, ok: true, result }));
});
console.log('[test] mock app connected');
await sleep(200);

/* 5. Tools now relay to the app. */
const status = await client.callTool({ name: 'grlite_status', arguments: {} });
check(/pic_binary/.test(status.content[0].text), 'grlite_status relays app result');

const report = await client.callTool({ name: 'grlite_report', arguments: {} });
check(/MOCK REPORT/.test(report.content[0].text), 'grlite_report relays text');

const setp = await client.callTool({ name: 'grlite_set_particle', arguments: { index: 1, charge: 0.5 } });
check(/particle 1 updated/.test(setp.content[0].text), 'grlite_set_particle relays with params');

const step = await client.callTool({ name: 'grlite_step', arguments: { n: 5 } });
check(/step=43/.test(step.content[0].text), 'grlite_step relays');

await client.close();
app.close();
console.log(failures === 0 ? '\nALL BRIDGE E2E CHECKS PASSED' : `\n${failures} CHECK(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
