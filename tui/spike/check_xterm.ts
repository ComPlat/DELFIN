/**
 * Run the spike live against xterm.js (headless) and check the scrollback.
 *
 * xterm.js is the terminal inside VS Code and most web terminals, and the
 * emulator family where OpenTUI issue #1630 deletes rows instead of pushing
 * them into scrollback when the footer grows. Here it is a live terminal,
 * not a recording: the spike runs in a real pseudo-terminal (`script`), its
 * output goes into xterm.js, and what xterm.js answers to the spike's
 * queries (cursor position, device attributes) goes back -- after a delay,
 * to play a slow link (SSH, VS Code Remote) where the answer to the startup
 * cursor query arrives late.
 *
 *   --prefill N    N shell lines are on screen before the spike starts, so
 *                  it starts N rows down; they must survive above it
 *   --latency MS   delay of every reply from the terminal (default 0)
 *   --seed N       the spike's footer-growth fuzz seed
 *   --raw PATH     keep every byte the spike wrote, for a closer look
 *
 * Usage (from tui/): bun spike/check_xterm.ts [--cols 110] [--rows 35] ...
 */
import { Terminal } from "@xterm/headless"
import { mkdtempSync } from "node:fs"
import { tmpdir } from "node:os"
import { join, resolve } from "node:path"

const argv = process.argv.slice(2)
const opt = (name: string, fallback: number) => {
  const at = argv.indexOf(`--${name}`)
  return at >= 0 ? Number(argv[at + 1]) : fallback
}
const cols = opt("cols", 110)
const rows = opt("rows", 35)
const prefill = opt("prefill", 0)
const latency = opt("latency", 0)
const seed = opt("seed", 1)
const timeoutS = opt("timeout", 180)

const tui = resolve(import.meta.dir, "..")
const work = mkdtempSync(join(tmpdir(), "spike-xterm-"))
const reportPath = join(work, "report.json")
const stderrPath = join(work, "stderr.txt")
const bun = process.execPath

const term = new Terminal({ cols, rows, scrollback: 200_000, allowProposedApi: true })
const write = (data: string | Uint8Array) => new Promise<void>((done) => term.write(data, done))
for (let i = 1; i <= prefill; i++) await write(`shell line ${i}\r\n`)

const inner =
  `stty cols ${cols} rows ${rows}; cd '${tui}' && SPIKE_REPORT='${reportPath}' ` +
  `'${bun}' spike/inline.tsx --auto --seed ${seed} 2>'${stderrPath}'; echo SPIKE_EXIT=$?`
const child = Bun.spawn(["script", "-q", "-E", "never", "-e", "-c", inner, "/dev/null"], {
  stdin: "pipe",
  stdout: "pipe",
  stderr: "pipe",
  env: { ...process.env, TERM: "xterm-256color", COLORTERM: "truecolor" },
})

let replies = 0
term.onData((data) => {
  replies++
  setTimeout(() => {
    try {
      child.stdin.write(data)
      child.stdin.flush()
    } catch {
      // The spike has exited; a late reply has nowhere to go.
    }
  }, latency)
})

const rawAt = argv.indexOf("--raw")
const raw: Uint8Array[] = []
const killer = setTimeout(() => child.kill(), timeoutS * 1000)
for await (const chunk of child.stdout) {
  if (rawAt >= 0) raw.push(chunk)
  await write(chunk)
}
if (rawAt >= 0) await Bun.write(argv[rawAt + 1]!, new Blob(raw))
const exitCode = await child.exited
clearTimeout(killer)

const buffer = term.buffer.normal
const lines: string[] = []
for (let i = 0; i < buffer.length; i++) lines.push(buffer.getLine(i)?.translateToString(true) ?? "")
const history = lines.join("\n")
await Bun.write(join(work, "history.txt"), history)

// The rules live in one place: scrollback_check.py.
const exitSeen = /SPIKE_EXIT=(\d+)/.exec(history)
const problems: string[] = []
if (!exitSeen) problems.push("the spike never finished (no SPIKE_EXIT)")
else if (exitSeen[1] !== "0") problems.push(`spike exited with ${exitSeen[1]}`)
if (exitCode !== 0) problems.push(`script exited with ${exitCode}`)
const stderr = await Bun.file(stderrPath).exists() ? (await Bun.file(stderrPath).text()).trim() : ""
if (stderr) problems.push(`stderr: ${stderr.slice(0, 300)}`)
if (!(await Bun.file(reportPath).exists())) problems.push("no report written")

const label = `xterm.js ${cols}x${rows}, latency ${latency} ms, from row ${prefill + 1}, seed ${seed}`
let checkOut = ""
if (await Bun.file(reportPath).exists()) {
  const checker = Bun.spawnSync(
    ["python3", "-I", join(tui, "spike", "scrollback_check.py"), join(work, "history.txt"), reportPath,
     "--prefill", String(prefill), "--replays", "0", "--width", String(cols)],
    { stdout: "pipe", stderr: "pipe" },
  )
  checkOut = checker.stdout.toString().trim()
  if (checker.exitCode !== 0) problems.push(...checkOut.split("\n").slice(1).filter((l) => l.trim() && l.trim() !== "FAIL"))
}
console.log(`${label}: ${checkOut.split("\n")[0] ?? ""}; ${replies} replies sent; work dir ${work}`)
if (problems.length) {
  console.log("FAIL")
  for (const p of problems) console.log(`  ${p.startsWith("  -") ? p.trim() : `- ${p}`}`)
  process.exit(1)
}
console.log("PASS")
