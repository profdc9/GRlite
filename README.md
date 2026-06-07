# GRlite

[![Deploy to GitHub Pages](https://github.com/profdc9/GRlite/actions/workflows/deploy.yml/badge.svg)](https://github.com/profdc9/GRlite/actions/workflows/deploy.yml)

**▶ [Live sandbox](https://profdc9.github.io/GRlite/)** — runs in the browser
(WebGL2 + WASM), no install. Auto-deployed from `main` to GitHub Pages by
[`.github/workflows/deploy.yml`](.github/workflows/deploy.yml).

A 2D linearized General Relativity sandbox for interactive pedagogical use.
The current authoritative spec is [`docs/gr_sandbox_v38.tex`](docs/gr_sandbox_v38.tex)
(v32 is preserved as the original); the staged build plan is in §12 of v32 and
the v33–v38 revisions layer subsequent work on top. The latest design state is
captured in v38's "Parameterizable Macroparticle Kernel" section.

This README covers how to build and run what's currently implemented.

## Status

Original 17-stage build plan: complete. Beyond that, the EM PIC physics has
been developed through Stages 18–46 across the v33–v38 design revisions.

### Original 17 stages
- **Stage 1** — scalar wave equation, free propagation ✓
- **Stage 2** — absorbing damping layer ✓
- **Stage 3** — static CIC source, Poisson convergence ✓
- **Stage 4** — all six wave equations, Lorenz-gauge monitoring ✓
- **Stage 5** — moving source, vector potentials, Yee curl ✓
- **Stage 6** — sampled background field arrays ✓
- **Stage 7** — single test particle, Boris pusher, Keplerian orbit ✓
- **Stage 8** — relativistic corrections, perihelion precession (EIH 1PN) ✓
- **Stage 9** — proper time accumulation; gravitomagnetic clock effect ✓
- **Stage 10** — perturbation FDTD with moving particle source ✓
- **Stage 11** — Esirkepov continuity (discrete charge conservation) ✓
- **Stage 12** — binary orbit ✓
- **Stages 13–17** — damping/CFL parameter sweeps ✓

### EM extension and self-force investigation (v33–v38)
- **Stages 18–22** — Lewis-Birdsall force interpolation; gravitomagnetic coupling
- **Stages 23–26** — EM force law (cyclotron, uniform E, inductive E, EM Kepler)
- **Stages 27–28** — closed-loop EM PIC orbit; gravity inductive heating diagnostic
- **Stages 29–32** — EM inductive CFL sweeps; binary PIC
- **Stages 33–35** — TSC+LB infrastructure; EM chain sign verification; magnetic-only
- **Stage 37** — uniform-motion self-force diagnosis (started the v37/v38 arc)
- **Stages 38–40** — benchmark suite (neutral pair, reciprocity, radiation reaction)
- **Stages 43–46** — bump kernel validation (Tier-1 and Tier-2)

## EM PIC production setup (v38)

The closed-loop EM/GR PIC scheme uses:
- **Bump kernel** (`GR_SHAPE_BUMP`) at radius R ≈ 6 cells, applied consistently to
  ρ deposit, J Esirkepov, A edge interp, φ-gradient gather, and curl-A
- **Lewis-Birdsall** force interpolation (`GR_FORCE_INTERP_LEWIS_BIRDSALL`)
- **Half-step A** convention (1-step ∂_t A stencil, J source aligned to A's
  half-step time)
- **PML** damping layer (16 cells default)

The recommended way to set this up:

```c
gr_sim_t* sim = gr_sim_create(W, H, dx, c_eff, cfl);
gr_sim_use_pedagogical_defaults(sim);   // bump kernel R=6, LB, deposition, PML
// ... add particles, background, etc.
```

This single call sets the entire production bundle. See v38 §kernel_design for
the empirical validation (Stages 43–46).

## Repo layout

```
core/      C simulation core (compiled native for tests + WASM for the web)
  include/   public API headers
  src/       primitives — field.c (leapfrog), sim.c (lifecycle), deposit.c,
             particle.c, background.c, gauge.c
  scenarios/ each .c file is one scenario; registry.c is the lookup table
  tests/     native test binaries; each test = scenario + analytic assertion
  Makefile
web/       TypeScript + Vite frontend
  src/       main.ts loads the WASM and renders the field via WebGL2
docs/      LaTeX design documents
  gr_sandbox_v32.tex   original spec
  gr_sandbox_v33.tex   ... successive design revisions
  ...
  gr_sandbox_v38.tex   current authoritative spec
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

Builds and runs all stage tests. Stage 1 verifies wavefront propagation at
$r = c\,t$, $1/\sqrt{r}$ amplitude falloff in 2D, and CFL stability boundary
behavior. Subsequent stages build on top.

### EM PIC benchmark suite

```bash
cd core
mingw32-make benchmark-pic
```

Runs Stages 37–40 (the self-force diagnostic and benchmark suite). Recommended
after any change to the EM PIC scheme.

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
  a comment referencing its equation label and section in the authoritative
  design document (currently `docs/gr_sandbox_v38.tex`; v32 references remain
  valid as the labels are preserved across revisions).
- **Scenarios are the single source of truth**. Tests and the web frontend both
  load scenarios by name through `gr_sim_load_scenario`; no test or UI inlines
  initial conditions directly.
- **No new dependencies without justification**. The default is to write the ~50
  lines ourselves; a dep must bring significant, necessary capability.

## License

MIT — see [`LICENSE`](LICENSE).
