# grlite-bridge

A debug RPC bridge that lets an MCP client (Claude Code) query and drive the
running GRlite web app.

```
  Claude Code ──MCP (stdio)── server.mjs ──WebSocket (127.0.0.1:8787)── web app (dev)
```

- The **server is long-lived** (Claude Code spawns it once per session).
- The **web app is the WebSocket client**; it connects/reconnects on every
  reload, so restarting the app never restarts the bridge.
- The app does the computation and returns **text** (field stats, per-particle
  state, conservation sums, 1D profiles) — the agent reads numbers, not pixels.
- Dev-only: the app-side client (`web/src/dev/bridge.ts`) starts only under the
  Vite dev server (`import.meta.env.DEV`) and the server binds to localhost.

## Setup

```bash
cd tools/grlite-bridge && npm install
```

Register with Claude Code (project `.mcp.json`, already provided):

```json
{ "mcpServers": { "grlite-bridge": {
    "command": "node",
    "args": ["C:/Users/dmarks/Documents/cursor/GRlite/tools/grlite-bridge/server.mjs"] } } }
```

Restart Claude Code so it spawns the server and exposes the `grlite_*` tools.
Then run the web app (`cd web && npm run dev`) and open it; it auto-connects.

## Tools

`grlite_status` · `grlite_report` · `grlite_particles` · `grlite_conservation` ·
`grlite_field_stats{field}` · `grlite_field_profile{field,axis,at,stride}` ·
`grlite_step{n}` · `grlite_play` · `grlite_pause` · `grlite_reset` ·
`grlite_build{scenario}` · `grlite_set_particle{index,...}` ·
`grlite_set_global{gEff,kE,viewField,paused}` · `grlite_screenshot`

## Test

```bash
npm test     # spawns the server, connects a mock app, checks tool relay
```
