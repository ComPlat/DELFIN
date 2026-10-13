/**
 * Replay a recorded spike run into xterm.js (headless) and check the
 * scrollback it ends up with.
 *
 * xterm.js is the terminal inside VS Code and most web terminals, and it is
 * one of the emulators where OpenTUI issue #1630 reports rows being deleted
 * instead of pushed into scrollback when the footer grows. Recording the raw
 * bytes in a real pseudo-terminal and feeding them to the real emulator tests
 * exactly that, with no window and no guessing.
 *
 * Usage (from tui/):
 *   script -q -e -f -O rec.raw -c "stty cols 110 rows 35; \
 *     SPIKE_REPORT=report.json bun spike/inline.tsx --auto" /dev/null
 *   bun spike/check_xterm.ts rec.raw report.json 110 35
 */
import { Terminal } from "@xterm/headless"

const [recording, reportPath, colsArg, rowsArg] = process.argv.slice(2)
if (!recording || !reportPath) {
  console.error("usage: bun spike/check_xterm.ts <recording> <report.json> <cols> <rows>")
  process.exit(2)
}
const cols = Number(colsArg ?? 110)
const rows = Number(rowsArg ?? 35)

let data = await Bun.file(recording).text()
// `script` frames the log with a "Script started/done" line on its own.
data = data.replace(/^Script started[^\n]*\n/, "").replace(/\nScript done[^\n]*\n?$/, "")

const term = new Terminal({ cols, rows, scrollback: 200_000, allowProposedApi: true })
await new Promise<void>((resolve) => term.write(data, resolve))

const buffer = term.buffer.normal
const lines: string[] = []
for (let i = 0; i < buffer.length; i++) lines.push(buffer.getLine(i)?.translateToString(true) ?? "")
const history = lines.join("\n")
await Bun.write(recording + ".xterm.txt", history)

const report = JSON.parse(await Bun.file(reportPath).text()) as { markers?: string[] }
const problems: string[] = []
let pos = 0
for (const marker of report.markers ?? []) {
  const at = history.indexOf(marker, pos)
  if (at < 0) problems.push(`marker missing or out of order: ${JSON.stringify(marker)}`)
  else pos = at
}
const numbers = [...history.matchAll(/replay line (\d+) —/g)].map((m) => Number(m[1]))
const inOrder = numbers.every((n, i) => n === i + 1)
if (numbers.length !== 1000 || !inOrder) {
  const seen = new Set(numbers)
  const missing = [...Array(1000).keys()].map((i) => i + 1).filter((n) => !seen.has(n))
  problems.push(
    `replay lines: ${numbers.length} seen, ${numbers.length - seen.size} duplicated, ` +
      `${missing.length} missing (first missing: ${missing.slice(0, 5).join(", ")})`,
  )
}

console.log(`xterm.js ${cols}x${rows}: ${buffer.length} buffer lines (scrollback + screen)`)
if (problems.length) {
  console.log("FAIL")
  for (const p of problems) console.log(`  - ${p}`)
  process.exit(1)
}
console.log("PASS: markers in order, 1000 replay lines exactly once")
