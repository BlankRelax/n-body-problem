// Builds the simulator's GUI for the web.
//
//   node scripts/build.mjs
//
// 1. compiles the `nbody-gui` crate to WebAssembly (release),
// 2. generates the JavaScript bindings with wasm-bindgen,
// 3. writes dist/ — a static site (index.html, the module, the stylesheet)
//    that any web server can host, and that other sites embed by calling the
//    module's exported `start(container)`,
// 4. writes dist-single/n-body-problem.html — one self-contained file with
//    the WebAssembly inlined, which opens straight from disk in a browser.
//
// Needs: Rust with the wasm32-unknown-unknown target, and wasm-bindgen-cli at
// the exact version pinned in crates/nbody-gui/Cargo.toml
// (cargo install wasm-bindgen-cli --version <v> --locked).

import { execFileSync } from 'node:child_process';
import { cpSync, mkdirSync, readFileSync, rmSync, statSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = fileURLToPath(new URL('..', import.meta.url));
const DIST = join(ROOT, 'dist');
const SINGLE = join(ROOT, 'dist-single');
const WEB = join(ROOT, 'web');
const NAME = 'nbody_gui';

const run = (cmd, args) => execFileSync(cmd, args, { cwd: ROOT, stdio: 'inherit' });

run('cargo', ['build', '--release', '--target', 'wasm32-unknown-unknown', '-p', 'nbody-gui']);
rmSync(DIST, { recursive: true, force: true });
run('wasm-bindgen', [
  '--target', 'web',
  '--no-typescript',
  '--out-dir', join(DIST, 'pkg'),
  join(ROOT, 'target', 'wasm32-unknown-unknown', 'release', `${NAME}.wasm`),
]);
for (const file of ['index.html', 'main.js', 'nbody.css']) {
  cpSync(join(WEB, file), join(DIST, file));
}

// --- the single-file build ------------------------------------------------------
const glue = readFileSync(join(DIST, 'pkg', `${NAME}.js`), 'utf8');
const wasm = readFileSync(join(DIST, 'pkg', `${NAME}_bg.wasm`)).toString('base64');
const css = readFileSync(join(WEB, 'nbody.css'), 'utf8');
if (glue.includes('</script')) throw new Error('bindings contain </script>');
const page = readFileSync(join(WEB, 'index.html'), 'utf8')
  .replace(/\s*<link rel="stylesheet"[^>]*>/g, '')
  .replace('</head>', `    <style>\n${css}\n    </style>\n  </head>`)
  .replace(
    /<script type="module" src="main.js"><\/script>/,
    () => `<script type="module">
${glue}
initSync({ module: Uint8Array.from(atob(${JSON.stringify(wasm)}), (c) => c.charCodeAt(0)) });
start(document.getElementById('game'));
</script>`,
  );
rmSync(SINGLE, { recursive: true, force: true });
mkdirSync(SINGLE);
writeFileSync(join(SINGLE, 'n-body-problem.html'), page);

const size = (p) => `${(statSync(p).size / 1024 / 1024).toFixed(2)} MB`;
console.log(`dist/pkg/${NAME}_bg.wasm  ${size(join(DIST, 'pkg', `${NAME}_bg.wasm`))}`);
console.log(`dist-single/n-body-problem.html  ${size(join(SINGLE, 'n-body-problem.html'))}`);
