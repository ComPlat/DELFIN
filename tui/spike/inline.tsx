/**
 * Phase 0 spike: can OpenTUI's split-footer (inline) mode carry the
 * delfin-agent TUI?
 *
 * Throwaway code with one job: exercise every inline risk named in
 * docs/TUI_PLAN.md (phase 0) on fake data, so the go/no-go rests on what real
 * terminals did rather than on what the documentation promises:
 *
 *   - a markdown answer streamed into native scrollback, block by block, with
 *     the unfinished tail live in the footer;
 *   - tool rows and a highlighted diff;
 *   - a footer that grows and shrinks (approval panel, full-height agent view);
 *   - divider + replay when another session is opened;
 *   - a 1,000-line replay, with its time measured;
 *   - resize: clear screen and scrollback, reprint the session at the new width;
 *   - a clean exit.
 *
 * Run from tui/:  bun spike/inline.tsx          interactive
 *                 bun spike/inline.tsx --auto   scripted run, then exit
 * SPIKE_REPORT=path writes timings and the expected scrollback markers as JSON.
 */
import {
  DiffRenderable,
  MarkdownRenderable,
  SyntaxStyle,
  TextAttributes,
  createCliRenderer,
  getTreeSitterClient,
  parseColor,
  type CliRenderer,
  type KeyEvent,
  type ScrollbackSurface,
} from "@opentui/core"
import { render, useKeyboard, writeSolidToScrollback } from "@opentui/solid"
import { For, Show, createEffect, createRoot, createSignal, type JSX } from "solid-js"

// --------------------------------------------------------------------------
// Look: "DELFIN blue, restrained" — one accent from the logo, neutral text,
// colour only for state.
// --------------------------------------------------------------------------

const C = {
  accent: "#6290C0",
  accentSoft: "#8FB2DA",
  text: "#D5DAE1",
  bright: "#F1F4F8",
  muted: "#8B94A3",
  faint: "#555D6B",
  ok: "#7FB685",
  warn: "#D9A65B",
  err: "#D46F6A",
} as const

const color = (hex: string) => parseColor(hex)

const MD_STYLES = {
  default: { fg: color(C.text) },
  "markup.heading": { fg: color(C.accentSoft), bold: true },
  "markup.heading.1": { fg: color(C.accent), bold: true },
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

const SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
const AUTO = process.argv.includes("--auto")
const REPORT_PATH = process.env.SPIKE_REPORT ?? ""
const REPLAY_CAP = 400

// --------------------------------------------------------------------------
// Transcript items. Everything written to scrollback is recorded here so a
// resize can reprint the session at the new width.
// --------------------------------------------------------------------------

type Item =
  | { kind: "banner" }
  | { kind: "user"; text: string }
  | { kind: "tool"; name: string; arg: string; summary: string; tail?: string[]; ok?: boolean }
  | { kind: "markdown"; text: string }
  | { kind: "diff"; path: string; diff: string; filetype: string }
  | { kind: "notice"; text: string }
  | { kind: "divider"; label: string }
  | { kind: "turnEnd"; text: string }
  | { kind: "lines"; lines: string[] }

function Row(props: { item: Item; width: number }): JSX.Element {
  const item = props.item
  switch (item.kind) {
    case "banner":
      return (
        <box flexDirection="column" paddingTop={1} paddingBottom={1}>
          <text>
            <span style={{ fg: C.accent, attributes: TextAttributes.BOLD }}>◆ DELFIN</span>
            <span style={{ fg: C.muted }}>{"  agent · kit/kit.deepseek-v4-flash · solo"}</span>
          </text>
          <text fg={C.faint}>{"  ~/project/delfin · git tui-rewrite · approval plan · isolation bwrap"}</text>
        </box>
      )
    case "user":
      return (
        <box paddingTop={1}>
          <text wrapMode="word">
            <span style={{ fg: C.accent, attributes: TextAttributes.BOLD }}>{"❯ "}</span>
            <span style={{ fg: C.bright }}>{item.text}</span>
          </text>
        </box>
      )
    case "tool":
      return (
        <box flexDirection="column">
          <text wrapMode="char">
            <span style={{ fg: item.ok === false ? C.err : C.ok }}>{"⏺ "}</span>
            <span style={{ fg: C.text, attributes: TextAttributes.BOLD }}>{item.name}</span>
            <span style={{ fg: C.muted }}>{"  " + item.arg}</span>
          </text>
          <text fg={C.muted}>{"  ⎿  " + item.summary}</text>
          <For each={item.tail ?? []}>{(line) => <text fg={C.faint}>{"     " + line}</text>}</For>
        </box>
      )
    case "notice":
      return <text fg={C.warn}>! {item.text}</text>
    case "divider": {
      const label = ` ${item.label} `
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
    case "turnEnd":
      return <text fg={C.faint}>{"  " + item.text}</text>
    case "lines":
      return (
        <box flexDirection="column">
          <For each={item.lines}>{(line) => <text fg={C.muted}>{line}</text>}</For>
        </box>
      )
    default:
      return <text />
  }
}

// --------------------------------------------------------------------------
// Scrollback: one ordered queue. A streaming answer holds the queue until it
// finishes, so a tool row can never land in the middle of a paragraph.
// --------------------------------------------------------------------------

class MarkdownStream {
  private surface?: ScrollbackSurface
  private md?: MarkdownRenderable
  private style?: SyntaxStyle
  private content = ""
  private committedBlocks = 0
  private committedRow = 0
  private dirty = false
  private pumping: Promise<void> = Promise.resolve()
  private opened = false
  private finishing = false
  private resolveDone!: (text: string) => void
  readonly done = new Promise<string>((resolve) => (this.resolveDone = resolve))

  constructor(
    private renderer: CliRenderer,
    private onTail: (tail: string) => void,
  ) {}

  /** Called when the queue reaches this stream. */
  open(): void {
    this.surface = this.renderer.createScrollbackSurface({ startOnNewLine: true })
    this.style = SyntaxStyle.fromStyles(MD_STYLES, this.surface.renderContext.nativeScene!)
    this.md = new MarkdownRenderable(this.surface.renderContext, {
      id: `answer-${Date.now()}`,
      content: "",
      syntaxStyle: this.style,
      width: "100%",
      marginTop: 1,
      streaming: true,
      internalBlockMode: "top-level",
      conceal: true,
      tableOptions: { widthMode: "content" },
      treeSitterClient: getTreeSitterClient(),
    })
    this.surface.root.add(this.md)
    this.opened = true
    this.schedule()
  }

  append(delta: string): void {
    this.content += delta
    this.schedule()
  }

  finish(): Promise<string> {
    this.finishing = true
    this.schedule()
    return this.done
  }

  private schedule(): void {
    if (!this.opened) return
    if (this.dirty) return
    this.dirty = true
    this.pumping = this.pumping.then(() => this.flush())
  }

  private async flush(): Promise<void> {
    this.dirty = false
    const surface = this.surface!
    const md = this.md!
    const final = this.finishing
    md.content = this.content
    md.streaming = !final
    await surface.settle()
    if (final) {
      if (surface.height > this.committedRow) surface.commitRows(this.committedRow, surface.height)
      this.onTail("")
      surface.destroy()
      this.style!.destroy()
      this.resolveDone(this.content)
      return
    }
    const states = md._blockStates
    const stable = md._stableBlockCount
    if (stable > this.committedBlocks) {
      const next = states[stable]
      const last = states[stable - 1]!
      const end = next ? next.renderable.y : last.renderable.y + last.renderable.height
      if (end > this.committedRow) {
        surface.commitRows(this.committedRow, end)
        this.committedRow = end
      }
      this.committedBlocks = stable
    }
    this.onTail(states.slice(this.committedBlocks).map((s) => s.tokenRaw).join(""))
    if (this.dirty) return // a newer flush is already queued behind us
  }
}

class Scrollback {
  readonly history: Item[] = []
  private chain: Promise<void> = Promise.resolve()
  private pending = 0

  constructor(
    private renderer: CliRenderer,
    private onError: (error: unknown) => void,
  ) {}

  push(item: Item, record = true): Promise<void> {
    if (record) this.history.push(item)
    return this.enqueue(() => this.write(item))
  }

  stream(onTail: (tail: string) => void): MarkdownStream {
    const stream = new MarkdownStream(this.renderer, onTail)
    this.enqueue(async () => {
      stream.open()
      const text = await stream.done
      this.history.push({ kind: "markdown", text })
    })
    return stream
  }

  /** Resolves once everything queued so far has been written. */
  drain(): Promise<void> {
    return this.enqueue(() => undefined)
  }

  private enqueue(task: () => Promise<void> | void): Promise<void> {
    // Demand-driven rendering stalls on a large burst of scrollback commits,
    // so the renderer runs live while anything is queued and drops back to
    // demand-driven (zero idle CPU) once the queue is empty.
    if (this.pending++ === 0) this.renderer.requestLive()
    this.chain = this.chain
      .then(task)
      .catch(this.onError)
      .finally(() => {
        if (--this.pending === 0) this.renderer.dropLive()
      })
    return this.chain
  }

  private async write(item: Item): Promise<void> {
    if (item.kind === "markdown") return this.writeMarkdown(item.text)
    if (item.kind === "diff") return this.writeDiff(item)
    writeSolidToScrollback(this.renderer, (ctx) => <Row item={item} width={ctx.width} />)
  }

  private async writeMarkdown(text: string): Promise<void> {
    const stream = new MarkdownStream(this.renderer, () => undefined)
    stream.open()
    stream.append(text)
    await stream.finish()
  }

  private async writeDiff(item: Extract<Item, { kind: "diff" }>): Promise<void> {
    writeSolidToScrollback(this.renderer, () => (
      <text>
        <span style={{ fg: C.ok }}>{"⏺ "}</span>
        <span style={{ fg: C.text, attributes: TextAttributes.BOLD }}>edit_file</span>
        <span style={{ fg: C.muted }}>{"  " + item.path}</span>
      </text>
    ))
    const surface = this.renderer.createScrollbackSurface({ startOnNewLine: true })
    const style = SyntaxStyle.fromStyles(MD_STYLES, surface.renderContext.nativeScene!)
    const diff = new DiffRenderable(surface.renderContext, {
      id: `diff-${Date.now()}`,
      diff: item.diff,
      view: "unified",
      filetype: item.filetype,
      syntaxStyle: style,
      showLineNumbers: true,
      wrapMode: "char",
      width: "100%",
      fg: C.text,
      lineNumberFg: C.faint,
      addedBg: "#1D2B22",
      removedBg: "#2E1F22",
      addedContentBg: "#1D2B22",
      removedContentBg: "#2E1F22",
      addedSignColor: C.ok,
      removedSignColor: C.err,
      treeSitterClient: getTreeSitterClient(),
    })
    surface.root.add(diff)
    surface.render()
    await surface.settle()
    surface.commitRows(0, surface.height)
    surface.destroy()
    style.destroy()
  }
}

// --------------------------------------------------------------------------
// Fake data
// --------------------------------------------------------------------------

const ANSWER = `## What was wrong

ORCA 6.1 prints the final energy as \`FINAL SINGLE POINT ENERGY\` with an extra
column, so the parser picked up the **wrong field** and every redox potential
came out shifted by the same constant.

- \`delfin/xtb.py\` now reads the value with a regular expression
- the old ORCA 5 layout still parses, so archived runs are unaffected
- three new tests cover both layouts and a truncated file

| Output | Before | After |
| --- | --- | --- |
| ORCA 5.0.4 | ok | ok |
| ORCA 6.1.1 | wrong field | ok |
| truncated | crash | clear error |

\`\`\`typescript
const energy = Number(line.match(/ENERGY\\s+(-?\\d+\\.\\d+)/)?.[1])
\`\`\`

Run \`pytest tests/test_xtb.py\` to check it yourself — all three pass.
`

const EDIT_DIFF = `--- a/delfin/xtb.py
+++ b/delfin/xtb.py
@@ -41,5 +41,8 @@ def read_final_energy(path):
     for line in lines:
         if "FINAL SINGLE POINT ENERGY" in line:
-            energy = float(line.split()[3])
+            match = _ENERGY_RE.search(line)
+            if match is None:
+                raise ValueError(f"no energy on line: {line!r}")
+            energy = float(match.group(1))
             break
     return energy
`

interface FakeSession {
  name: string
  state: "input" | "working" | "done" | "elsewhere"
  detail: string
  age: string
}

const SESSIONS: FakeSession[] = [
  { name: "tddft-benzene", state: "input", detail: "approve bash: orca input.inp", age: "2m" },
  { name: "fix-xtb-parser", state: "working", detail: "edit_file delfin/xtb.py · $0.02", age: "5m" },
  { name: "redox-table", state: "done", detail: "wrote reports/redox.md", age: "12m" },
  { name: "dashboard", state: "elsewhere", detail: "kit.deepseek-v4-flash · idle", age: "20m" },
]

const GROUPS: { state: FakeSession["state"]; title: string; glyph: string; tint: string }[] = [
  { state: "input", title: "NEEDS INPUT", glyph: "●", tint: C.warn },
  { state: "working", title: "WORKING", glyph: "◐", tint: C.accent },
  { state: "done", title: "DONE", glyph: "✓", tint: C.ok },
  { state: "elsewhere", title: "ELSEWHERE", glyph: "○", tint: C.faint },
]

const sleep = (ms: number) => new Promise<void>((resolve) => setTimeout(resolve, ms))

function chunk(text: string, size: number): string[] {
  const out: string[] = []
  for (let i = 0; i < text.length; i += size) out.push(text.slice(i, i + size))
  return out
}

// --------------------------------------------------------------------------
// State shared by the footer and the scripted turns
// --------------------------------------------------------------------------

type Posture = "plan" | "default" | "acceptEdits"

interface Approval {
  tool: string
  target: string
  preview: string[]
  resolve: (answer: "y" | "n" | "a") => void
}

const [phase, setPhase] = createSignal<{ label: string; t0: number } | null>(null)
const [tail, setTail] = createSignal("")
const [approval, setApproval] = createSignal<Approval | null>(null)
const [agentsOpen, setAgentsOpen] = createSignal(false)
const [agentIndex, setAgentIndex] = createSignal(0)
const [posture, setPosture] = createSignal<Posture>("plan")
const [tick, setTick] = createSignal(0)
const [termRows, setTermRows] = createSignal(24)
const [termCols, setTermCols] = createSignal(80)
const [hint, setHint] = createSignal("")

const report: Record<string, unknown> = { markers: [] as string[] }
const marker = (text: string) => (report.markers as string[]).push(text)

// --------------------------------------------------------------------------
// Footer
// --------------------------------------------------------------------------

function Rule(): JSX.Element {
  return <text fg={C.faint}>{"─".repeat(Math.max(1, termCols() - 1))}</text>
}

function ApprovalPanel(props: { req: Approval }): JSX.Element {
  return (
    <box
      flexDirection="column"
      border
      borderStyle="rounded"
      borderColor={C.accent}
      title={` ${props.req.tool} · ${props.req.target} `}
      titleAlignment="left"
      paddingLeft={1}
      paddingRight={1}
    >
      <For each={props.req.preview}>
        {(line) => (
          <text fg={line.startsWith("+") ? C.ok : line.startsWith("-") ? C.err : C.text}>{line}</text>
        )}
      </For>
      <text>
        <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>[y]</span>
        <span style={{ fg: C.text }}>{" yes   "}</span>
        <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>[n]</span>
        <span style={{ fg: C.text }}>{" no   "}</span>
        <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>[e]</span>
        <span style={{ fg: C.text }}>{" stop asking about writes   "}</span>
        <span style={{ fg: C.accentSoft, attributes: TextAttributes.BOLD }}>[a]</span>
        <span style={{ fg: C.text }}>{" abort"}</span>
      </text>
    </box>
  )
}

function LiveRow(): JSX.Element {
  const elapsed = () => {
    tick()
    const p = phase()
    return p ? ((performance.now() - p.t0) / 1000).toFixed(1) : "0.0"
  }
  const tailLines = () => {
    const width = Math.max(10, termCols() - 6)
    return tail()
      .trimEnd()
      .split("\n")
      .filter((line) => line.trim())
      .slice(-3)
      .map((line) => (line.length > width ? "…" + line.slice(line.length - width + 1) : line))
  }
  return (
    <box flexDirection="column">
      <text>
        <span style={{ fg: C.accent }}>{SPINNER[tick() % SPINNER.length] + " "}</span>
        <span style={{ fg: C.text }}>{phase()?.label}</span>
        <span style={{ fg: C.faint }}>{`  ${elapsed()}s · kit.deepseek-v4-flash · esc to interrupt`}</span>
      </text>
      <For each={tailLines()}>{(line) => <text fg={C.muted}>{"  " + line}</text>}</For>
    </box>
  )
}

function SessionFooter(): JSX.Element {
  return (
    <box flexDirection="column" width="100%">
      <Show when={phase()}>
        <LiveRow />
      </Show>
      <Show when={approval()}>{(req: () => Approval) => <ApprovalPanel req={req()} />}</Show>
      <Rule />
      <box flexDirection="row">
        <text fg={C.accent} attributes={TextAttributes.BOLD}>{"❯ "}</text>
        <input
          ref={(el) => (inputRef = el)}
          flexGrow={1}
          focused={!approval() && !agentsOpen()}
          placeholder="ask DELFIN…  (/approve /diff /replay 1000 /quit, ← agents)"
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
      <text fg={C.faint}>
        {hint() || `  ${posture()} · esc interrupt · shift+tab approval mode · ← agents · /help`}
      </text>
    </box>
  )
}

function AgentView(): JSX.Element {
  const needs = SESSIONS.filter((s) => s.state === "input").length
  const flat = () => SESSIONS
  return (
    <box flexDirection="column" width="100%">
      <text>
        <span style={{ fg: C.accent, attributes: TextAttributes.BOLD }}>◆ DELFIN agents</span>
        <span style={{ fg: C.faint }}>{`   ${SESSIONS.length} sessions · `}</span>
        <span style={{ fg: C.warn }}>{`${needs} needs input`}</span>
      </text>
      <Rule />
      <For each={GROUPS}>
        {(group) => (
          <Show when={SESSIONS.some((s) => s.state === group.state)}>
            <text fg={C.faint} attributes={TextAttributes.BOLD}>
              {group.title}
            </text>
            <For each={SESSIONS.filter((s) => s.state === group.state)}>
              {(s) => {
                const selected = () => flat()[agentIndex()] === s
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
          </Show>
        )}
      </For>
      <Rule />
      <text fg={C.faint}>
        <span style={{ fg: C.accent }}>{"❯ "}</span>{"type to start a new session…"}
      </text>
      <text fg={C.faint}>{"  ↑↓ select · → open · space peek · ctrl+x stop · esc back"}</text>
    </box>
  )
}

function App(props: { renderer: CliRenderer }): JSX.Element {
  createEffect(() => {
    // The footer is exactly as tall as what it shows: agent view, or live
    // row + tail + approval panel + composer.
    let rows: number
    if (agentsOpen()) {
      rows = 2 + SESSIONS.length + GROUPS.filter((g) => SESSIONS.some((s) => s.state === g.state)).length + 3
    } else {
      rows = 4
      if (phase()) rows += 1 + Math.min(3, tail().trimEnd().split("\n").filter((l) => l.trim()).length)
      const req = approval()
      if (req) rows += req.preview.length + 3
    }
    props.renderer.footerHeight = Math.max(4, Math.min(rows, termRows() - 1))
  })

  useKeyboard((key: KeyEvent) => onKey(props.renderer, key))

  return (
    <box flexDirection="column" width="100%" height="100%" justifyContent="flex-end">
      <Show when={agentsOpen()} fallback={<SessionFooter />}>
        <AgentView />
      </Show>
    </box>
  )
}

// --------------------------------------------------------------------------
// Behaviour
// --------------------------------------------------------------------------

let renderer!: CliRenderer
let scrollback!: Scrollback
let inputRef: { value: string } | undefined
const [draft, setDraft] = createSignal("")
let turnActive = false
let interrupted = false
let shuttingDown = false
let replayPending = false

function onKey(r: CliRenderer, key: KeyEvent): void {
  if (key.ctrl && (key.name === "c" || key.name === "d")) {
    key.preventDefault()
    void shutdown()
    return
  }
  const req = approval()
  if (req) {
    if (key.name === "y" || key.name === "e") answer(req, "y")
    else if (key.name === "n" || key.name === "escape") answer(req, "n")
    else if (key.name === "a") answer(req, "a")
    key.preventDefault()
    return
  }
  if (agentsOpen()) {
    if (key.name === "up") setAgentIndex((i) => Math.max(0, i - 1))
    else if (key.name === "down") setAgentIndex((i) => Math.min(SESSIONS.length - 1, i + 1))
    else if (key.name === "right" || key.name === "return" || key.name === "enter") void openSession(SESSIONS[agentIndex()]!)
    else if (key.name === "escape" || key.name === "left") setAgentsOpen(false)
    key.preventDefault()
    return
  }
  if (key.name === "left" && draft() === "") {
    key.preventDefault()
    setAgentsOpen(true)
    return
  }
  if (key.name === "tab" && key.shift) {
    key.preventDefault()
    const order: Posture[] = ["plan", "default", "acceptEdits"]
    setPosture(order[(order.indexOf(posture()) + 1) % order.length]!)
    return
  }
  if (key.name === "escape" && turnActive) {
    interrupted = true
    setHint("  ! esc — ending this turn after the running step")
  }
  void r
}

function answer(req: Approval, value: "y" | "n" | "a"): void {
  setApproval(null)
  req.resolve(value)
}

function ask(tool: string, target: string, preview: string[], autoMs: number): Promise<"y" | "n" | "a"> {
  return new Promise((resolve) => {
    const req: Approval = { tool, target, preview, resolve }
    setApproval(req)
    if (AUTO) setTimeout(() => approval() === req && answer(req, "y"), autoMs)
  })
}

async function runTool(name: string, arg: string, ms: number, summary: string, tail?: string[]): Promise<void> {
  setPhase({ label: `${name}  ${arg}`, t0: performance.now() })
  await sleep(ms)
  marker(arg)
  await scrollback.push({ kind: "tool", name, arg, summary, tail })
}

async function fakeTurn(prompt: string): Promise<void> {
  if (turnActive) return
  turnActive = true
  interrupted = false
  const t0 = performance.now()
  marker(prompt)
  await scrollback.push({ kind: "user", text: prompt })
  setPhase({ label: "thinking…", t0 })
  await sleep(700)
  if (!interrupted) await runTool("read_file", "delfin/xtb.py", 400, "212 lines · 7.9 kB")
  if (!interrupted) await runTool("grep_file", "FINAL SINGLE POINT ENERGY", 300, "3 matches in 2 files")
  if (!interrupted) {
    setPhase({ label: "waiting for approval", t0: performance.now() })
    const verdict = await ask("edit_file", "delfin/xtb.py", EDIT_DIFF.split("\n").slice(5, 12), 900)
    if (verdict === "y") {
      marker("raise ValueError(f\"no energy on line")
      await scrollback.push({ kind: "diff", path: "delfin/xtb.py", diff: EDIT_DIFF, filetype: "python" })
    } else {
      await scrollback.push({ kind: "notice", text: "edit_file refused" })
      interrupted = verdict === "a"
    }
  }
  if (!interrupted) {
    await runTool("bash", "pytest -q tests/test_xtb.py", 1200, "exit 0 · 0.42 s", [
      "tests/test_xtb.py ...                                   [100%]",
      "3 passed in 0.42s",
    ])
  }
  if (!interrupted) {
    setPhase({ label: "writing…", t0: performance.now() })
    const stream = scrollback.stream(setTail)
    for (const piece of chunk(ANSWER, 7)) {
      if (interrupted) break
      stream.append(piece)
      await sleep(12)
    }
    await stream.finish()
    marker("What was wrong")
  }
  const secs = ((performance.now() - t0) / 1000).toFixed(1)
  const end = interrupted ? `✗ interrupted after ${secs}s` : `✓ ${secs}s · ↑12.4k ↓1.1k tokens · kit.deepseek-v4-flash`
  await scrollback.push({ kind: "turnEnd", text: end })
  setPhase(null)
  setTail("")
  setHint("")
  turnActive = false
  if (replayPending) void replayAll()
}

async function openSession(s: FakeSession): Promise<void> {
  setAgentsOpen(false)
  marker(`${s.name} · last 2 turns`)
  await scrollback.push({ kind: "divider", label: `${s.name} · last 2 turns` })
  await scrollback.push({ kind: "user", text: "the xtb parser fails on ORCA 6.1 outputs, fix it" })
  await scrollback.push({ kind: "tool", name: "grep_file", arg: "TOTAL ENERGY", summary: "3 matches" })
  await scrollback.push({ kind: "markdown", text: "Found it: the energy column moved in ORCA 6.1. Fixing the parser next." })
}

async function replayLines(n: number): Promise<number> {
  const t0 = performance.now()
  const width = String(n).length
  for (let i = 0; i < n; i += 10) {
    const lines: string[] = []
    for (let j = i; j < Math.min(n, i + 10); j++) {
      lines.push(`replay line ${String(j + 1).padStart(width, "0")} — tool output and answer text of a long session`)
    }
    void scrollback.push({ kind: "lines", lines })
  }
  await scrollback.drain()
  await renderer.idle()
  const ms = performance.now() - t0
  marker(`replay line ${String(1).padStart(width, "0")}`)
  marker(`replay line ${String(n).padStart(width, "0")}`)
  await scrollback.push({ kind: "notice", text: `replayed ${n} lines in ${ms.toFixed(0)} ms` })
  return ms
}

async function replayAll(): Promise<void> {
  if (turnActive) {
    replayPending = true // never wipe what a running turn is still writing
    return
  }
  replayPending = false
  await scrollback.drain()
  renderer.resetSplitFooterForReplay({ clearSavedLines: true })
  for (const item of scrollback.history.slice(-REPLAY_CAP)) void scrollback.push(item, false)
  await scrollback.drain()
}

async function submit(value: string): Promise<void> {
  const text = value.trim()
  if (inputRef) inputRef.value = ""
  setDraft("")
  if (!text) return
  if (text === "/quit" || text === "/exit") return shutdown()
  if (text === "/approve") {
    const verdict = await ask("bash", "orca input.inp", ["orca input.inp > orca.out", "(starts a 4-core ORCA job in this folder)"], 900)
    await scrollback.push({ kind: "notice", text: `approval answered: ${verdict}` })
    return
  }
  if (text === "/diff") {
    await scrollback.push({ kind: "diff", path: "delfin/xtb.py", diff: EDIT_DIFF, filetype: "python" })
    return
  }
  if (text.startsWith("/replay")) {
    const n = Number(text.split(/\s+/)[1] ?? 1000) || 1000
    await replayLines(n)
    return
  }
  if (text === "/redraw") return replayAll()
  await fakeTurn(text)
}

async function shutdown(): Promise<void> {
  if (shuttingDown) return
  shuttingDown = true
  await scrollback.drain()
  if (REPORT_PATH) {
    report.historyItems = scrollback.history.length
    await Bun.write(REPORT_PATH, JSON.stringify(report, null, 2))
  }
  renderer.destroy()
  // OpenTUI flushes queued scrollback in a microtask after destroy; the
  // process ends on its own once that output is written. The timer is only
  // a backstop and does not keep the process alive.
  setTimeout(() => process.exit(0), 5000).unref()
}

async function autoScript(): Promise<void> {
  await sleep(400)
  await fakeTurn("the xtb parser fails on ORCA 6.1 outputs — fix it and run the tests")
  // The tallest footer this TUI will ever have: the agent view.
  setAgentsOpen(true)
  await sleep(800)
  setAgentIndex(1)
  await sleep(300)
  await openSession(SESSIONS[1]!)
  report.replay1000Ms = await replayLines(1000)
  await sleep(300)
  await fakeTurn("summarise what changed for the changelog")
  await sleep(600)
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
  scrollback = new Scrollback(renderer, (error) => {
    report.errors = [...((report.errors as string[]) ?? []), String(error)]
  })
  setTermRows(renderer.terminalHeight)
  setTermCols(renderer.terminalWidth)

  let resizeTimer: ReturnType<typeof setTimeout> | undefined
  // "resize" also fires when the FOOTER changes height (the render region
  // is the footer in split-footer mode). Only a new terminal WIDTH changes
  // how committed rows wrap, so only that triggers the clear-and-replay;
  // anything else would wipe the user's scrollback for nothing.
  let replayWidth = renderer.terminalWidth
  renderer.on("resize", () => {
    setTermRows(renderer.terminalHeight)
    setTermCols(renderer.terminalWidth)
    if (renderer.terminalWidth === replayWidth) return
    clearTimeout(resizeTimer)
    resizeTimer = setTimeout(() => {
      replayWidth = renderer.terminalWidth
      report.replays = ((report.replays as number) ?? 0) + 1
      void replayAll()
    }, 250)
  })
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

  await render(() => <App renderer={renderer} />, renderer)
  await scrollback.push({ kind: "banner" })
  if (AUTO) void autoScript()
}

await main()
