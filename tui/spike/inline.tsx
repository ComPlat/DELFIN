/**
 * Phase 0 spike: OpenTUI's split-footer (inline) mode, on the real
 * scrollback core (src/scrollback/).
 *
 * Throwaway UI, real mechanics. It plays fake agent turns through the same
 * commit queue, transcript, markdown block splitting and resize replay the
 * TUI will use, and exercises every inline risk in docs/TUI_PLAN.md:
 *
 *   - answers committed block by block, the unfinished block live (as
 *     markdown) in the footer;
 *   - tool rows, a highlighted diff;
 *   - a footer that grows and shrinks at random (approval panels of random
 *     height, the full-height agent view), from any start row;
 *   - divider + replay, a realistic 30-turn replay, a 1,000-line replay;
 *   - resize: clear and reprint at the new width, also mid-turn;
 *   - a clean exit.
 *
 *   bun spike/inline.tsx            interactive
 *   bun spike/inline.tsx --auto     scripted run, then exit (test mode: tags)
 *   bun spike/inline.tsx --keys     interactive, logging what each key did
 *   --seed N                        footer-growth fuzz seed (default 1)
 *
 * SPIKE_REPORT=path: written at exit -- timings, errors, replays, the key
 * log and ``expected``: every tag the scrollback must hold, in order.
 */
import { SyntaxStyle, TextAttributes, createCliRenderer, parseColor, type CliRenderer, type KeyEvent } from "@opentui/core"
import { render, useKeyboard, writeSolidToScrollback } from "@opentui/solid"
import { For, Show, createEffect, createRoot, createSignal, type JSX } from "solid-js"
import { rename } from "node:fs/promises"
import { AnswerBlocks, splitMarkdown } from "../src/scrollback/blocks"
import { CommitQueue } from "../src/scrollback/queue"
import { replayOnWidthChange } from "../src/scrollback/resize"
import { Transcript, type Item } from "../src/scrollback/transcript"
import { commitSurface, diffBlock, isDrawableDiff, markdownBlock } from "../src/scrollback/writers"
import { Tagger, answer as fakeAnswer, diff as fakeDiff, rng, tagsIn } from "./content"

// The kitty-graphics probe retitles tmux panes; this TUI draws no images.
process.env.OPENTUI_GRAPHICS ??= "0"

const AUTO = process.argv.includes("--auto")
const SEED = Number(process.argv[process.argv.indexOf("--seed") + 1]) || 1
const REPORT_PATH = process.env.SPIKE_REPORT ?? ""
const REPLAY_LIMIT = 400

// --------------------------------------------------------------------------
// Look: "DELFIN blue, restrained"
// --------------------------------------------------------------------------

const C = {
  accent: "#6290C0", accentSoft: "#8FB2DA", text: "#D5DAE1", bright: "#F1F4F8",
  muted: "#8B94A3", faint: "#555D6B", ok: "#7FB685", warn: "#D9A65B", err: "#D46F6A",
} as const
const color = (hex: string) => parseColor(hex)
const heading = (hex: string) => ({ fg: color(hex), bold: true })
const MD_STYLES = {
  default: { fg: color(C.text) },
  "markup.heading": heading(C.accentSoft),
  "markup.heading.1": heading(C.accent), "markup.heading.2": heading(C.accentSoft),
  "markup.heading.3": heading(C.accentSoft), "markup.heading.4": heading(C.text),
  "markup.heading.5": heading(C.text), "markup.heading.6": heading(C.text),
  "markup.strong": { fg: color(C.bright), bold: true },
  "markup.italic": { fg: color(C.text), italic: true },
  "markup.list": { fg: color(C.accent) },
  "markup.raw": { fg: color("#C9B98F") },
  "markup.link": { fg: color(C.accentSoft), underline: true },
  "markup.link.label": { fg: color(C.accentSoft), underline: true },
  "markup.link.url": { fg: color(C.muted), underline: true },
  conceal: { fg: color(C.faint) },
  keyword: { fg: color(C.accentSoft), bold: true },
  string: { fg: color("#A9C28E") },
  comment: { fg: color(C.faint), italic: true },
  number: { fg: color(C.warn) },
  function: { fg: color(C.bright) },
  type: { fg: color("#C9B98F") },
  variable: { fg: color(C.text) },
  property: { fg: color(C.text) },
}
const DIFF_COLORS = { text: C.text, lineNumber: C.faint, added: "#1D2B22", removed: "#2E1F22",
                      addedSign: C.ok, removedSign: C.err }
const SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"

// --------------------------------------------------------------------------
// Items and how each one is written
// --------------------------------------------------------------------------

type Fields = Record<string, any>

/** The strings an item shows, in display order (where its tags are). */
function shown(item: Item): string {
  const f = item as Fields
  switch (item.kind) {
    case "banner": return `${f.line1}\n${f.line2}`
    case "tool": return [f.arg, f.summary, ...(f.tail ?? [])].join("\n")
    case "lines": return (f.lines as string[]).join("\n")
    default: return String(f.text ?? "")
  }
}

function Row(props: { item: Item; width: number }): JSX.Element {
  const f = props.item as Fields
  switch (props.item.kind) {
    case "banner":
      return (
        <box flexDirection="column" paddingTop={1} paddingBottom={1}>
          <text>
            <span style={{ fg: C.accent, attributes: TextAttributes.BOLD }}>{"◆ DELFIN"}</span>
            <span style={{ fg: C.muted }}>{"  " + f.line1}</span>
          </text>
          <text fg={C.faint}>{"  " + f.line2}</text>
        </box>
      )
    case "user":
      return (
        <box paddingTop={1}>
          <text wrapMode="word">
            <span style={{ fg: C.accent, attributes: TextAttributes.BOLD }}>{"❯ "}</span>
            <span style={{ fg: C.bright }}>{f.text}</span>
          </text>
        </box>
      )
    case "tool":
      return (
        <box flexDirection="column">
          <text wrapMode="char">
            <span style={{ fg: f.ok === false ? C.err : C.ok }}>{"⏺ "}</span>
            <span style={{ fg: C.text, attributes: TextAttributes.BOLD }}>{f.name}</span>
            <span style={{ fg: C.muted }}>{"  " + f.arg}</span>
          </text>
          <text fg={C.muted}>{"  ⎿  " + f.summary}</text>
          <For each={(f.tail ?? []) as string[]}>{(line) => <text fg={C.faint}>{"     " + line}</text>}</For>
        </box>
      )
    case "divider": {
      const label = ` ${f.text} `
      const side = Math.max(2, Math.floor((props.width - label.length) / 2))
      return (
        <box paddingTop={1}>
          <text>
            <span style={{ fg: C.accent }}>{"━".repeat(side)}</span>
            <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>{label}</span>
            <span style={{ fg: C.accent }}>{"━".repeat(Math.max(0, props.width - side - label.length))}</span>
          </text>
        </box>
      )
    }
    case "notice":
      return <text fg={C.warn}>{"! " + f.text}</text>
    case "turnEnd":
      return <text fg={C.faint}>{"  " + f.text}</text>
    case "lines":
      return (
        <box flexDirection="column">
          <For each={f.lines as string[]}>{(line) => <text fg={C.muted}>{line}</text>}</For>
        </box>
      )
    default:
      return <text>{String(f.text ?? "")}</text>
  }
}

let renderer!: CliRenderer
let queue!: CommitQueue
let transcript!: Transcript
const report: Record<string, unknown> = { errors: [] as string[], keyLog: [] as string[] }
let expected: string[] = []

async function writeItem(item: Item): Promise<void> {
  const f = item as Fields
  if (item.kind === "markdown") {
    await commitSurface(renderer, markdownBlock(f.text, { styles: MD_STYLES, marginTop: 1 }))
  } else if (item.kind === "diff") {
    writeSolidToScrollback(renderer, () => (
      <text>
        <span style={{ fg: C.ok }}>{"⏺ "}</span>
        <span style={{ fg: C.text, attributes: TextAttributes.BOLD }}>{"edit_file"}</span>
        <span style={{ fg: C.muted }}>{"  " + f.path}</span>
      </text>
    ))
    if (isDrawableDiff(f.text)) {
      await commitSurface(renderer, diffBlock(f.text, { styles: MD_STYLES, filetype: "python", colors: DIFF_COLORS }))
    } else {
      // Never "Error parsing diff" in the transcript: plain lines instead.
      const lines = String(f.text).split("\n")
      writeSolidToScrollback(renderer, () => <Row item={{ id: item.id, kind: "lines", lines }} width={0} />)
    }
  } else {
    writeSolidToScrollback(renderer, (ctx) => <Row item={item} width={ctx.width} />)
  }
  expected.push(...tagsIn(shown(item)))
}

let nextId = 0
const add = (kind: string, fields: Fields) => transcript.append({ id: `i${++nextId}`, kind, ...fields })

// --------------------------------------------------------------------------
// Footer state
// --------------------------------------------------------------------------

type Posture = "plan" | "default" | "acceptEdits"
interface Approval { tool: string; target: string; preview: string[]; resolve: (a: "y" | "n") => void }

const [phase, setPhase] = createSignal<{ label: string; t0: number } | null>(null)
const [tail, setTail] = createSignal("")
const [approval, setApproval] = createSignal<Approval | null>(null)
const [agentsOpen, setAgentsOpen] = createSignal(false)
const [agentIndex, setAgentIndex] = createSignal(0)
const [posture, setPosture] = createSignal<Posture>("plan")
const [tick, setTick] = createSignal(0)
const [termRows, setTermRows] = createSignal(24)
const [termCols, setTermCols] = createSignal(80)
const [draft, setDraft] = createSignal("")

const SESSIONS = [
  { name: "tddft-benzene", state: "input", detail: "approve bash: orca input.inp", age: "2m" },
  { name: "fix-xtb-parser", state: "working", detail: "edit_file delfin/xtb.py · $0.02", age: "5m" },
  { name: "redox-table", state: "done", detail: "wrote reports/redox.md", age: "12m" },
  { name: "dashboard", state: "elsewhere", detail: "kit.deepseek-v4-flash · idle", age: "20m" },
]
const GROUPS = [
  { state: "input", title: "NEEDS INPUT", glyph: "●", tint: C.warn },
  { state: "working", title: "WORKING", glyph: "◐", tint: C.accent },
  { state: "done", title: "DONE", glyph: "✓", tint: C.ok },
  { state: "elsewhere", title: "ELSEWHERE", glyph: "○", tint: C.faint },
]

const TAIL_MAX = 6
const tailRows = () => Math.min(TAIL_MAX, tail().trimEnd().split("\n").filter((l) => l.trim()).length)

function Rule(): JSX.Element {
  return <text fg={C.faint}>{"─".repeat(Math.max(1, termCols() - 1))}</text>
}

let tailStyle: SyntaxStyle | undefined

function Footer(props: { renderer: CliRenderer }): JSX.Element {
  tailStyle ??= SyntaxStyle.fromStyles(MD_STYLES, props.renderer.nativeScene!)
  createEffect(() => {
    let rows: number
    if (agentsOpen()) {
      rows = 2 + SESSIONS.length + GROUPS.length + 3
    } else {
      rows = 4 + (phase() ? 1 : 0) + tailRows()
      const req = approval()
      if (req) rows += req.preview.length + 3
    }
    props.renderer.footerHeight = Math.max(4, Math.min(rows, termRows() - 1))
  })
  useKeyboard((key: KeyEvent) => onKey(key))
  return (
    <box flexDirection="column" width="100%" height="100%" justifyContent="flex-end">
      <Show when={agentsOpen()} fallback={<SessionFooter />}>
        <AgentView />
      </Show>
    </box>
  )
}

function SessionFooter(): JSX.Element {
  const elapsed = () => {
    tick()
    const p = phase()
    return p ? ((performance.now() - p.t0) / 1000).toFixed(1) : "0.0"
  }
  return (
    <box flexDirection="column" width="100%">
      <Show when={tailRows() > 0}>
        <scrollbox height={tailRows()} stickyScroll stickyStart="bottom">
          <markdown content={tail()} streaming syntaxStyle={tailStyle!} conceal width="100%" />
        </scrollbox>
      </Show>
      <Show when={phase()}>
        <text>
          <span style={{ fg: C.accent }}>{SPINNER[tick() % SPINNER.length] + " "}</span>
          <span style={{ fg: C.text }}>{phase()?.label ?? ""}</span>
          <span style={{ fg: C.faint }}>{`  ${elapsed()}s · kit.deepseek-v4-flash · esc to interrupt`}</span>
        </text>
      </Show>
      <Show when={approval()}>
        {(req: () => Approval) => (
          <box flexDirection="column" border borderStyle="rounded" borderColor={C.accent}
               title={` ${req().tool} · ${req().target} `} titleAlignment="left" paddingLeft={1} paddingRight={1}>
            <For each={req().preview}>
              {(line) => <text fg={line.startsWith("+") ? C.ok : line.startsWith("-") ? C.err : C.text}>{line}</text>}
            </For>
            <text>
              <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>{"[y]"}</span>
              <span style={{ fg: C.text }}>{" yes   "}</span>
              <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>{"[n]"}</span>
              <span style={{ fg: C.text }}>{" no   "}</span>
              <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>{"[e]"}</span>
              <span style={{ fg: C.text }}>{" stop asking about writes   "}</span>
              <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>{"[a]"}</span>
              <span style={{ fg: C.text }}>{" abort"}</span>
            </text>
          </box>
        )}
      </Show>
      <Rule />
      <box flexDirection="row">
        <text fg={C.accent} attributes={TextAttributes.BOLD}>{"❯ "}</text>
        <input
          ref={(el) => (inputRef = el)}
          flexGrow={1}
          focused={!approval() && !agentsOpen()}
          placeholder="ask DELFIN…  (/approve /replay 1000 /redraw /quit, ← agents)"
          placeholderColor={C.faint}
          textColor={C.bright}
          backgroundColor="transparent"
          focusedBackgroundColor="transparent"
          cursorColor={C.accent}
          onInput={(value: string) => setDraft(value)}
          // The binding types onSubmit as an intersection of the event and the
          // value signatures; at run time the input passes its value.
          onSubmit={((value: string) => void submit(value)) as never}
        />
      </box>
      <Rule />
      <text fg={C.faint}>{`  ${posture()} · esc interrupt · shift+tab approval mode · ← agents · /help`}</text>
    </box>
  )
}

function AgentView(): JSX.Element {
  return (
    <box flexDirection="column" width="100%">
      <text>
        <span style={{ fg: C.accent, attributes: TextAttributes.BOLD }}>{"◆ DELFIN agents"}</span>
        <span style={{ fg: C.faint }}>{`   ${SESSIONS.length} sessions · `}</span>
        <span style={{ fg: C.warn }}>{"1 needs input"}</span>
      </text>
      <Rule />
      <For each={GROUPS}>
        {(group) => (
          <>
            <text fg={C.faint} attributes={TextAttributes.BOLD}>{group.title}</text>
            <For each={SESSIONS.filter((s) => s.state === group.state)}>
              {(s) => {
                const selected = () => SESSIONS[agentIndex()] === s
                return (
                  <text>
                    <span style={{ fg: selected() ? C.accent : C.faint }}>{selected() ? "▸ " : "  "}</span>
                    <span style={{ fg: group.tint }}>{group.glyph + " "}</span>
                    <span style={{ fg: selected() ? C.bright : C.text, attributes: selected() ? TextAttributes.BOLD : 0 }}>
                      {s.name.padEnd(18)}
                    </span>
                    <span style={{ fg: C.muted }}>{s.detail.padEnd(Math.max(10, termCols() - 34))}</span>
                    <span style={{ fg: C.faint }}>{s.age.padStart(4)}</span>
                  </text>
                )
              }}
            </For>
          </>
        )}
      </For>
      <Rule />
      <text fg={C.faint}><span style={{ fg: C.accent }}>{"❯ "}</span>{"type to start a new session…"}</text>
      <text fg={C.faint}>{"  ↑↓ select · → open · space peek · ctrl+x stop · esc back"}</text>
    </box>
  )
}

// --------------------------------------------------------------------------
// Behaviour
// --------------------------------------------------------------------------

let inputRef: { value: string } | undefined
let turnActive = false
let shuttingDown = false
const tagger = new Tagger(AUTO)
const sleep = (ms: number) => new Promise<void>((r) => setTimeout(r, ms))
const keyLog = report.keyLog as string[]

function onKey(key: KeyEvent): void {
  if (key.ctrl && (key.name === "c" || key.name === "d")) {
    key.preventDefault()
    void shutdown()
    return
  }
  const req = approval()
  if (req) {
    if (key.name === "y" || key.name === "n") {
      keyLog.push(`approval:${key.name}`)
      setApproval(null)
      req.resolve(key.name)
    }
    key.preventDefault()
    return
  }
  if (agentsOpen()) {
    if (key.name === "up") setAgentIndex((i) => Math.max(0, i - 1))
    else if (key.name === "down") setAgentIndex((i) => Math.min(SESSIONS.length - 1, i + 1))
    else if (key.name === "right" || key.name === "return" || key.name === "enter") {
      keyLog.push(`agents:open-session:${SESSIONS[agentIndex()]!.name}`)
      void openSession(SESSIONS[agentIndex()]!.name)
    } else if (key.name === "escape" || key.name === "left") {
      keyLog.push("agents:close")
      setAgentsOpen(false)
    }
    key.preventDefault()
    return
  }
  if (key.name === "left" && draft() === "") {
    key.preventDefault()
    keyLog.push("agents:open")
    setAgentsOpen(true)
    return
  }
  if (key.name === "tab" && key.shift) {
    key.preventDefault()
    const order: Posture[] = ["plan", "default", "acceptEdits"]
    const next = order[(order.indexOf(posture()) + 1) % order.length]!
    keyLog.push(`posture:${next}`)
    setPosture(next)
    return
  }
  if (key.name === "escape" && turnActive) keyLog.push("turn:interrupt")
}

function ask(tool: string, target: string, preview: string[], autoMs: number): Promise<"y" | "n"> {
  return new Promise((resolve) => {
    const req: Approval = { tool, target, preview, resolve }
    setApproval(req)
    if (AUTO) setTimeout(() => { if (approval() === req) { setApproval(null); resolve("y") } }, autoMs)
  })
}

async function fakeTurn(prompt: string, turn: number, r: () => number): Promise<void> {
  turnActive = true
  const t0 = performance.now()
  void add("user", { text: tagger.line(prompt) })
  setPhase({ label: "thinking…", t0 })
  await sleep(300)
  const tools = [["read_file", "delfin/xtb.py", "212 lines · 7.9 kB"],
                 ["grep_file", "FINAL SINGLE POINT ENERGY", "3 matches in 2 files"],
                 ["bash", "pytest -q tests/test_xtb.py", "exit 0 · 0.42 s"]] as const
  for (const [name, arg, summary] of tools) {
    setPhase({ label: `${name}  ${arg}`, t0: performance.now() })
    await sleep(120)
    void add("tool", { name, arg: tagger.line(arg), summary: tagger.line(summary),
                       tail: name === "bash" ? [tagger.line("3 passed in 0.42s")] : [] })
    // The footer grows and shrinks at random while rows are committed: the
    // case opentui#1630 is about, from wherever the footer happens to sit.
    const size = 3 + Math.floor(r() * 12)
    const preview = Array.from({ length: size }, (_, i) => `${i % 2 ? "+" : " "} preview line ${i + 1} of ${size}`)
    await ask("bash", arg, preview, 150 + Math.floor(r() * 150))
  }
  if (turn === 1) void add("diff", { path: "delfin/xtb.py", text: fakeDiff(tagger) })
  setPhase({ label: "writing…", t0: performance.now() })
  const blocks = new AnswerBlocks((block) => void add("markdown", { text: block }), setTail)
  const text = fakeAnswer(tagger, turn)
  for (let i = 0; i < text.length;) {
    const n = 3 + Math.floor(r() * 18)
    blocks.append(text.slice(i, i + n))
    i += n
    await sleep(8)
  }
  blocks.finish()
  const secs = ((performance.now() - t0) / 1000).toFixed(1)
  await add("turnEnd", { text: tagger.line(`✓ ${secs}s · ↑12.4k ↓1.1k tokens · kit.deepseek-v4-flash`) })
  setPhase(null)
  setTail("")
  turnActive = false
}

async function openSession(name: string): Promise<void> {
  setAgentsOpen(false)
  void add("divider", { text: tagger.line(`${name} · last 2 turns`) })
  void add("user", { text: tagger.line("the xtb parser fails on ORCA 6.1 outputs, fix it") })
  void add("tool", { name: "grep_file", arg: tagger.line("TOTAL ENERGY"), summary: tagger.line("3 matches") })
  await add("markdown", { text: tagger.line("Found it: the energy column moved in ORCA 6.1.") })
}

/** `turns` realistic turns through markdown, diff and tree-sitter: what a real replay costs. */
async function realisticReplay(turns: number): Promise<number> {
  const t0 = performance.now()
  void add("divider", { text: tagger.line(`long-session · ${turns} turns`) })
  for (let k = 1; k <= turns; k++) {
    void add("user", { text: tagger.line(`turn ${k}: check the parser again`) })
    void add("tool", { name: "read_file", arg: tagger.line("delfin/xtb.py"), summary: tagger.line("212 lines") })
    if (k % 3 === 0) void add("diff", { path: "delfin/xtb.py", text: fakeDiff(tagger) })
    for (const block of splitMarkdown(fakeAnswer(tagger, k), true).blocks) void add("markdown", { text: block })
  }
  await queue.drain()
  await renderer.idle()
  return performance.now() - t0
}

async function replayLines(n: number): Promise<number> {
  const t0 = performance.now()
  for (let i = 0; i < n; i += 10) {
    const lines: string[] = []
    for (let j = i; j < Math.min(n, i + 10); j++) {
      lines.push(tagger.line(`replay line ${String(j + 1).padStart(4, "0")} — tool output and answer text`))
    }
    void add("lines", { lines })
  }
  await queue.drain()
  await renderer.idle()
  return performance.now() - t0
}

async function submit(value: string): Promise<void> {
  const text = value.trim()
  if (inputRef) inputRef.value = ""
  setDraft("")
  if (!text) return
  if (text === "/quit" || text === "/exit") return shutdown()
  if (text === "/approve") {
    await ask("bash", "orca input.inp", ["orca input.inp > orca.out", "(starts a 4-core ORCA job)"], 900)
    return
  }
  if (text.startsWith("/replay")) {
    const n = Number(text.split(/\s+/)[1] ?? 1000) || 1000
    const ms = await replayLines(n)
    await add("notice", { text: `replayed ${n} lines in ${ms.toFixed(0)} ms` })
    return
  }
  if (text === "/redraw") return transcript.replay(REPLAY_LIMIT)
  await fakeTurn(text, 1, rng(SEED))
}

async function shutdown(): Promise<void> {
  if (shuttingDown) return
  shuttingDown = true
  await queue.drain()
  if (REPORT_PATH) {
    report.expected = expected
    report.replays = transcript.replayCount
    report.items = transcript.items.length
    await Bun.write(REPORT_PATH + ".tmp", JSON.stringify(report, null, 2))
    await rename(REPORT_PATH + ".tmp", REPORT_PATH)
  }
  renderer.destroy()
  // OpenTUI flushes queued scrollback after destroy; the process ends on its
  // own once that is written. The timer is only a backstop.
  setTimeout(() => process.exit(0), 5000).unref()
}

async function autoScript(): Promise<void> {
  const r = rng(SEED)
  await sleep(300)
  for (let turn = 1; turn <= 3; turn++) await fakeTurn(`fix the xtb parser, round ${turn}`, turn, r)
  setAgentsOpen(true)
  await sleep(500)
  setAgentIndex(1)
  await sleep(200)
  await openSession("fix-xtb-parser")
  report.realReplay30Ms = Math.round(await realisticReplay(30))
  report.replay1000Ms = Math.round(await replayLines(1000))
  await fakeTurn("summarise what changed for the changelog", 4, r)
  await sleep(300)
  await shutdown()
}

async function main(): Promise<void> {
  renderer = await createCliRenderer({
    screenMode: "split-footer",
    footerHeight: 4,
    externalOutputMode: "capture-stdout",
    useMouse: false,
    exitOnCtrlC: false,
    targetFps: 30,
    maxFps: 60,
    consoleMode: "disabled",
    useKittyKeyboard: {},
  })
  const target = {
    requestLive: () => renderer.requestLive(),
    dropLive: () => renderer.dropLive(),
    resetSplitFooterForReplay: (options?: { clearSavedLines?: boolean }) => {
      expected = []
      renderer.resetSplitFooterForReplay(options)
    },
  }
  queue = new CommitQueue(target, (error) => (report.errors as string[]).push(String(error)))
  transcript = new Transcript(queue, target, writeItem)
  setTermRows(renderer.terminalHeight)
  setTermCols(renderer.terminalWidth)
  renderer.on("resize", () => {
    setTermRows(renderer.terminalHeight)
    setTermCols(renderer.terminalWidth)
  })
  replayOnWidthChange(renderer, transcript, { limit: REPLAY_LIMIT })
  let spinner: ReturnType<typeof setInterval> | undefined
  createRoot(() =>
    createEffect(() => {
      if (phase() && !spinner) spinner = setInterval(() => setTick((t) => t + 1), 80)
      if (!phase() && spinner) {
        clearInterval(spinner)
        spinner = undefined
      }
    }),
  )
  await render(() => <Footer renderer={renderer} />, renderer)
  void add("banner", { line1: tagger.line("agent · kit/kit.deepseek-v4-flash · solo"),
                       line2: tagger.line("~/project/delfin · git tui-rewrite · approval plan · isolation bwrap") })
  if (AUTO) void autoScript()
}

await main()
