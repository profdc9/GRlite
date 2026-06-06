# GRlite web sandbox

Browser front-end for the GRlite 2D linearized-GR + EM PIC core. Built first
and foremost as an **interactive sandbox for debugging the physics engine** —
visualizations surface bugs that raw numbers hide.

```
cd web
npm install
npm run dev        # Vite dev server (http://localhost:5173)
npm run build      # tsc --noEmit && vite build
```

The WASM core lives in `public/grlite/grlite.{js,wasm}`, built from `core/`
via `mingw32-make wasm` (see `core/Makefile`). Rebuild that whenever the C
changes; the dev page cache-busts the import so a fresh build is picked up.

## Architecture

Plain TypeScript + DOM (no UI framework). WebGL2 renders the scalar field;
a stacked Canvas2D overlay renders particles, trails, arrows, and selection.

```
src/
  wasm/binding.ts     loadCore() + the cwrap FFI surface (GRliteCore)
  sim/config.ts       grid/CFL/absorber/enum constants, scenario specs
  sim/world.ts        World: ergonomic wrapper over a sim handle
                      (typed particle reads, field views, buildScenario,
                      step, live setParticleFields via heap write)
  sim/fieldViews.ts   selectable field views; derived B_gz/B_z (curl) and
                      |g| (grad-mag) computed CPU-side per frame
  render/gl.ts        shader compile/link helpers
  render/fieldPass.ts WebGL2 scalar field + diverging colormap (auto-norm)
  render/overlay2d.ts Canvas2D overlay: trails, velocity arrows, markers,
                      selection, click-picking
  ui/controls.ts      toolbar + scenario/field selectors
  ui/inspector.ts     read-only global + selected-particle panels
  dev/diagnostics.ts  text reports (field stats, profiles, conservation,
                      stability probe)
  dev/bridge.ts       debug-bridge WebSocket client (dev only)
  state.ts            central AppState
  main.ts             thin wiring + frame loop
```

## Design decisions

- **JSON scenarios are the source of truth** (Phase 2, in progress). JS reads
  a scenario JSON and builds the sim through the WASM API rather than the C
  `gr_sim_load_scenario` path (C scenarios stay for native tests). Schema:
  `{ program, format, version, grid, global, background, particles, view }`,
  where `global` carries every sim switch + absorber + init config so a
  scenario fully and reproducibly determines a run. Scenarios are shareable
  via a compressed URL hash (a bug = a link) and double as engine tests.
- **Consistency model:** the running sim is a deterministic function of the
  JSON. Editing while paused at t=0 updates the canonical JSON (reproducible).
  Live edits during a run (e.g. dragging mass/charge) are allowed but flip a
  **live-modified** flag — the run is exploratory and not reproducible until
  reset. Reproducibility is guaranteed only if the scenario is unchanged
  after start.
- **Undo** (Phase 3) is a snapshot stack of the JSON; the field state is
  derived, so undo rebuilds rather than storing grid arrays.
- Reference architecture: `profdc9/PhysicsKitchen` (same UI constraints,
  plain TS, world/serialization/snapshot/undo + per-object & world panels).

## Debug bridge

`tools/grlite-bridge/` is an MCP (stdio) server with a WebSocket relay so an
agent can query and drive the running app and get **computed text** back
(the app does the math). See `tools/grlite-bridge/README.md`. The app-side
client (`src/dev/bridge.ts`) is dev-only and auto-reconnects.

## Status (phases)

- **1a** ✅ module refactor (from the Stage-1 monolith).
- **1b** ✅ debug bridge (MCP ↔ app).
- **1c** ✅ field selector, particle trails, velocity arrows, read-only
  inspector, click-to-select.
- **1c.1** ⬜ force-vector overlay (needs a per-particle force readout from C).
- **2** ⬜ JSON source of truth: schema, build-from-JSON, URL-hash
  load/save, scenario library mirroring the C scenarios.
- **3** ⬜ editable inspectors + live edits + undo.
- **4** ⬜ JSON + HTML test harness (reuses the bridge).
