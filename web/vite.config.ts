import { defineConfig } from 'vite';

// The WASM module (grlite.js + grlite.wasm) is produced by `mingw32-make wasm`
// in ../core/ and written to web/public/grlite/. Vite serves anything under
// public/ at the URL root, so the runtime loads it via `/grlite/grlite.js`.
//
// The dynamic import for that path is intentionally external — Rollup must not
// try to bundle a runtime asset that doesn't exist at build time. The
// `/* @vite-ignore */` comment handles the dev server; this `external` entry
// handles the production build.
export default defineConfig({
  server: { port: 5173 },
  build: {
    target: 'es2022',
    rollupOptions: {
      external: ['/grlite/grlite.js'],
    },
  },
});
