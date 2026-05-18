# GRlite

A 2D linearized General Relativity sandbox for interactive pedagogical use.
The current authoritative spec is [`docs/gr_sandbox_v33.tex`](docs/gr_sandbox_v33.tex)
(v32 is preserved as the original); the staged build plan is in §12.

This README covers how to build and run what's currently implemented.

## Status

Of the 17 stages in v33's §12:

- **Stage 1** — scalar wave equation, free propagation ✓
- **Stage 2** — absorbing damping layer ✓
- **Stage 3** — static CIC source, Poisson convergence ✓
- **Stage 6** — sampled background field arrays ✓ (built ahead of Stages 4–5
  so single-particle dynamics tests can be written against a prescribed background)
- Stages 4, 5, 7–17 — pending

## Repo layout

```
core/      C simulation core (compiled native for tests + WASM for the web)
  include/   public API headers
  src/       primitives — field.c (leapfrog), sim.c (lifecycle)
  scenarios/ each .c file is one scenario; registry.c is the lookup table
  tests/     native test binaries; each test = scenario + analytic assertion
  Makefile
web/       TypeScript + Vite frontend
  src/       main.ts loads the WASM and renders the field via WebGL2
docs/      LaTeX design documents (authoritative spec)
```

## Toolchain

- **C/WASM**: Emscripten **5.0.7** at `C:\Users\dmarks\src\emsdk`
- **Native C**: MinGW-w64 gcc **15.2.0** (`mingw32-make 4.4.1` for the build)
- **Frontend**: Node **25.x**, npm or pnpm; Vite **^6.0**, TypeScript **^5.5**
- **Optional desktop**: Tauri 2.10.1 + Rust 1.94 (not used yet)

Activation per shell (Emscripten doesn't persist to PATH):

```powershell
# PowerShell
& 'C:\Users\dmarks\src\emsdk\emsdk_env.ps1' | Out-Null
```

## Build & run

### Native tests (fast iteration, no browser)

```bash
cd core
mingw32-make test
```

Builds `core/build/native/stage01_wave` and runs it. Stage 1 verifies wavefront
propagation at $r = c\,t$, $1/\sqrt{r}$ amplitude falloff in 2D, and CFL stability
boundary behavior.

### WASM build + web frontend

```powershell
# 1. Activate Emscripten (per shell)
& 'C:\Users\dmarks\src\emsdk\emsdk_env.ps1' | Out-Null

# 2. Build the WASM into web/public/grlite/
cd core
mingw32-make wasm

# 3. Install web deps once, then start the dev server
cd ..\web
npm install
npm run dev
```

Open the URL Vite prints (typically <http://localhost:5173>). You'll see a
Gaussian pulse propagating outward from the center of a 256×256 grid.

### Clean

```bash
cd core
mingw32-make clean
```

## Build flags worth knowing

```bash
mingw32-make test SANITIZE=1     # turn on -fsanitize=address,undefined for native
mingw32-make CC=clang test       # use clang instead of gcc
```

## Code conventions

- **Formula traceability**: every formula or numerical method in C/TS code carries
  a comment referencing its equation label and section in `docs/gr_sandbox_v32.tex`
  (e.g. `Eq. (eq:leapfrog_field) — §9.2 "Time stepping: the leapfrog scheme"`).
  This lets reviewers cross-check code against spec without context-switching.
- **Scenarios are the single source of truth**. Tests and the web frontend both
  load scenarios by name through `gr_sim_load_scenario`; no test or UI inlines
  initial conditions directly.
- **No new dependencies without justification**. The default is to write the ~50
  lines ourselves; a dep must bring significant, necessary capability.

## License

MIT — see [`LICENSE`](LICENSE).
