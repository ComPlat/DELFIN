/**
 * The window as one executable: `bun build --compile` of src/main.tsx, with
 * Solid's JSX transform, and inside it the native library (the patched
 * build installed by build-opentui-native.sh), the tree-sitter worker and
 * its grammars.
 * A machine without Bun runs it.
 *
 *   bun scripts/build-binary.ts [target]      # default: this machine's
 *
 * Targets are Bun's: bun-linux-x64, bun-linux-arm64, bun-darwin-x64,
 * bun-darwin-arm64. The native library is the one of the machine it is
 * built on: a binary for another platform needs that platform's library in
 * node_modules (@opentui/core-<platform>) -- built there.
 *
 * Writes dist/delfin-tui-<os>-<arch>; the launcher (`delfin-agent tui`)
 * prefers it to running from source.
 */
import solidPlugin from "@opentui/solid/bun-plugin"
import { chmodSync, existsSync, mkdirSync, statSync } from "node:fs"
import { join } from "node:path"

const root = join(import.meta.dir, "..")
const here = `bun-${process.platform}-${process.arch}`
const target = process.argv[2] ?? here
const name = target.replace(/^bun-/, "")
const out = join(root, "dist", `delfin-tui-${name}`)
mkdirSync(join(root, "dist"), { recursive: true })

if (target !== here) {
  const lib = join(root, "node_modules", "@opentui", `core-${name}`)
  if (!existsSync(lib)) {
    console.error(`no native library for ${name} in node_modules (@opentui/core-${name}); build on that platform`)
    process.exit(2)
  }
}

// OpenTUI brings its tree-sitter worker, the grammars and the native
// library in with `import ... with { type: "file" }`: the compiler embeds
// them (overriding the worker's path breaks highlighting and markdown).
const started = Date.now()
const result = await Bun.build({
  entrypoints: [join(root, "src", "main.tsx")],
  target: "bun",
  plugins: [solidPlugin],
  minify: true,
  // Nothing of the folder it is started in: a project's bunfig.toml (its
  // preload would run code inside the window) or .env is not the window's.
  compile: { target: target as never, outfile: out, autoloadBunfig: false, autoloadDotenv: false },
})
if (!result.success) {
  for (const log of result.logs) console.error(log)
  process.exit(1)
}
chmodSync(out, 0o755)
const mb = (statSync(out).size / 1e6).toFixed(1)
console.log(`wrote ${out} (${mb} MB) in ${((Date.now() - started) / 1000).toFixed(1)} s`)
