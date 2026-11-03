/**
 * The live area under the scrollback: what the session in view is doing
 * now, the question it waits on, and the line the user types into.
 *
 *   message   a short note of the window's own (one row)
 *   tail      the answer's unfinished block, as markdown (≤ 6 rows), after
 *             the blocks that finished but are not in the rows yet
 *   status    spinner · what it does · elapsed · tokens · cost · esc
 *   tasks     the session's task list while a task is open (tasks.tsx)
 *   background  what the session still has out: shells, watched jobs,
 *             background agents, wake-ups, with their ages (background.tsx)
 *   approval  the open question, with exactly the keys the session offers
 *   composer  ❯ the prompt (grows to 8 rows; completion and search above);
 *             "!" in front: a shell line, said on the box's top border.
 *             It keeps the keys while a question is open (app.tsx).
 *   hints     posture · model · keys -- with the other sessions of the
 *             service ("2 more sessions (1 waits for you) ← agents"); at
 *             its right, quietly, what needs the user (◉ 3 · ✗ 1), or
 *             for a moment what just arrived
 *
 * The agent view (←) takes the whole terminal but its last row: its list
 * in place of the tail, status and question, the composer under it.
 *
 * Every line is cut to the terminal's width, so the footer knows exactly
 * how many rows it takes (`footerRows`), at any width: a footer taller
 * than its content wastes rows, a shorter one draws over itself.
 */
import { TextAttributes, infoStringToFiletype, type CliRenderer, type Renderable, type ScrollBoxRenderable,
         type SyntaxStyle } from "@opentui/core"
import { For, Show, createEffect, createMemo, createSignal, on, onCleanup, type Accessor, type JSX } from "solid-js"
import { AgentViewPanel, agentHelp, takesAnswer, type AgentView } from "../agents/view"
import { SwarmPanel, type SwarmView } from "../swarm/panel"
import { ComposerView, type Composer } from "../composer/composer"
import type { Keymap } from "../keymap/keymap"
import type { Asked, SessionState, TurnLive } from "../store/state"
import { GUTTER, SPINNER, postureLook, type Theme } from "../theme/theme"
import { fit, width } from "../text"
import { clickable } from "../mouse"
import { plainIfLate } from "../scrollback/highlighter"
import { SegmentedMarkdown, holdsTable } from "../scrollback/segmented"
import { TranscriptView, themedStyle, type ScreenItems } from "../scrollback/screenview"
import { fenceOf, liveBlocks, liveRows } from "./live"
import { ApprovalPanel, approvalRows, type ApprovalUi } from "./approval"
import { PickerView, type Picker } from "./picker"
import { TasksPanel, taskLines } from "./tasks"
import { BackgroundPanel, backgroundLines } from "./background"

export { TAIL_MAX } from "./live"

export interface FooterProps {
  renderer: CliRenderer
  theme: () => Theme
  keymap: Keymap
  session: () => SessionState | null
  connection: () => string
  /** Starting the first session: said in the hints line. */
  starting?: () => boolean
  message: () => string
  /** The first step onto bypass waiting for y: shown in the message's
   *  row (in red), in its place, so the footer's height is the same. */
  offer?: () => string
  /** The offer's clickable part, and what a click on it does. */
  offerYes?: string
  onOfferClick?: () => void
  composer: Composer
  approval: ApprovalUi
  picker?: Picker
  /** The task list folded to one line (ctrl+t). */
  tasksFolded?: () => boolean
  /** The agent view; while open it has the footer. */
  agents?: AgentView
  /** The swarm view (SW2b); while open it has the footer. */
  swarm?: SwarmView
  /** The view is drawn on the alternate screen: all of it is the view's,
   *  and the footer keeps its height for the way back. */
  agentsOnOwnScreen?: boolean
  /** What needs the user, for the right end of the hints line: what waits
   *  for an answer (other sessions here, questions parked anywhere), what
   *  failed, and a new one's note, shown for a few seconds instead. */
  attention?: () => { asks: number; failed: number; note: string }
  /** The other sessions running in the agent service (not the one in
   *  view, none asleep or ended), and how many of them wait for the user:
   *  said in the hints line beside the agent view's key. */
  others?: () => { running: number; waiting: number }
  /** A click on that part of the hints line (full screen): the agent view. */
  onOthersClick?: () => void
  focused: () => boolean
  /** Full screen: the footer is the window -- the session's rows (or the
   *  agent view) above, the composer at the bottom -- and the answer's
   *  tail streams at the rows' end instead of above the composer. */
  screen?: ScreenItems
  /** The rows' scroll box, for the keys that scroll it. */
  scrollRef?: (box: ScrollBoxRenderable) => void
  /** The wheel turned over the rows. */
  onUserScroll?: () => void
}

/** What the prompt box's top border says in shell mode, the longest that
 *  fits `room` columns (the words; the border adds a blank each side). */
export const SHELL_HINTS = [
  "shell · runs through the agent's checks · esc or backspace on an empty line to leave",
  "shell · runs through the agent's checks · esc to leave",
  "shell · esc to leave",
  "shell",
]

export function shellHint(room: number): string {
  return SHELL_HINTS.find((h) => width(h) <= room) ?? ""
}

/** The hints line's last part: the agent view's key, with how many other
 *  sessions run and wait -- the longest wording that fits `room` columns,
 *  the shortest when none does. */
export function othersHint(others: { running: number; waiting: number } | undefined, key: string,
                           room: number): string {
  const tail = `${key} agents`
  const n = others?.running ?? 0
  if (n <= 0) return tail
  const w = Math.min(n, others?.waiting ?? 0)
  const forms = [
    `${n} more session${n === 1 ? "" : "s"}${w ? ` (${w} wait${w === 1 ? "s" : ""} for you)` : ""} ${tail}`,
    `${n} more${w ? ` (${w} waiting)` : ""} ${tail}`,
    `${n} more ${tail}`,
  ]
  return forms.find((f) => width(f) <= room) ?? forms[forms.length - 1]!
}

/** Rows the footer takes with the agent view open: all but one. */
export function agentFooterRows(term: { height: number }): number {
  return Math.max(4, term.height - 1)
}


/** The rows everything but the approval panel takes. */
function fixedRows(s: SessionState | null, message: string, composerRows: number, pickerRows: number,
                   taskRows: number): number {
  return 3 + composerRows + pickerRows + taskRows + (s?.turn || s?.wait ? 1 : 0) + (message ? 1 : 0)
}

/** What the approval panel may take: the rest of the terminal (less one
 *  row: the footer never fills it), so its keys are always on screen. */
export function approvalBudget(s: SessionState | null, term: { width: number; height: number }, message: string,
                               composerRows = 1, pickerRows = 0, taskRows = 0): number {
  // Less a block on its way into the rows (there for a moment): the panel
  // must not push it out of sight.
  const live = s ? liveRows(s, term.width - GUTTER) : 0
  return Math.max(3, term.height - 1 - fixedRows(s, message, composerRows, pickerRows, taskRows) - live)
}

/** Rows the footer needs; never more than the terminal has (less one). */
export function footerRows(s: SessionState | null, term: { width: number; height: number }, message: string,
                           composerRows = 1, approval: ApprovalUi | null = null, pickerRows = 0,
                           taskRows = 0): number {
  let rows = fixedRows(s, message, composerRows, pickerRows, taskRows)
  const asked = s?.pending[0]
  // With a question open the turn waits: its tail is empty (a question
  // closes the text part) -- what is left is a block on its way into the
  // rows, there for a moment.
  if (s) rows += liveRows(s, term.width - GUTTER)
  rows += approvalRows(asked, approval?.state ?? null, term.width - GUTTER,
                       approvalBudget(s, term, message, composerRows, pickerRows, taskRows))
  return Math.max(4, Math.min(rows, Math.max(4, term.height - 1)))
}

function elapsed(startedAt: number, now: number): string {
  const secs = Math.max(0, (now - startedAt) / 1000)
  return secs < 60 ? `${secs.toFixed(1)}s` : `${Math.floor(secs / 60)}m${String(Math.floor(secs % 60)).padStart(2, "0")}s`
}

/** This turn's tokens and cost, as the REPL's turn line shows them: the
 *  difference from where the turn began; output estimated from the
 *  streamed characters (marked ~) until the provider counts it. */
export function turnUsage(turn: TurnLive, usage: { inputTokens: number; outputTokens: number; costUsd: number }): string[] {
  const tin = Math.max(0, usage.inputTokens - turn.base.inputTokens)
  const tout = Math.max(0, usage.outputTokens - turn.base.outputTokens)
  const out = tout > 0 ? `↓${compact(tout)}` : turn.streamed > 0 ? `↓~${compact(Math.max(1, Math.floor(turn.streamed / 4)))}` : "↓…"
  const parts = [`↑${compact(tin)} ${out}`]
  const cost = usage.costUsd - turn.base.costUsd
  if (cost > 0) parts.push(`$${cost.toFixed(4)}`)
  return parts
}

export function compact(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`
  if (n >= 1_000) return `${(n / 1_000).toFixed(1)}k`
  return String(n)
}

/**
 * The terminal's size, as a signal. Not useTerminalDimensions: in
 * split-footer mode that is the footer's own area. Set only when it
 * changes -- the resize event also fires when the footer changes height,
 * which this signal is used to compute.
 */
export function terminalSize(renderer: CliRenderer): Accessor<{ width: number; height: number }> {
  const [size, setSize] = createSignal({ width: renderer.terminalWidth, height: renderer.terminalHeight })
  const onResize = () => {
    const now = size()
    if (now.width !== renderer.terminalWidth || now.height !== renderer.terminalHeight) {
      setSize({ width: renderer.terminalWidth, height: renderer.terminalHeight })
    }
  }
  renderer.on("resize", onResize)
  onCleanup(() => renderer.off("resize", onResize))
  return size
}

export function Footer(props: FooterProps): JSX.Element {
  const t = props.theme
  /** The message row: the bypass offer while it waits, else the note. */
  const note = () => props.offer?.() || props.message()
  // Only the offer's "[y] switch" takes a click, never the rest of the row.
  const offerHit = clickable((col) => {
    if (!props.offer?.() || !props.offerYes) return null
    const shown = fit(note(), cols() - 1)
    const at = shown.indexOf(props.offerYes)
    return at >= 0 && col >= at && col < at + props.offerYes.length ? "yes" : null
  }, () => props.onOfferClick?.())
  const term = terminalSize(props.renderer)
  const [tick, setTick] = createSignal(0)
  const [now, setNow] = createSignal(Date.now())
  // The agent view's ages: once a second is enough (its rows are drawn
  // again on each change of it).
  const [clock, setClock] = createSignal(Date.now())
  const viewOpen = () => Boolean(props.agents?.state.open || props.swarm?.state.open)
  let timer: ReturnType<typeof setInterval> | undefined
  createEffect(() => {
    const running = Boolean(props.session()?.turn) || viewOpen()
    if (running && !timer) {
      timer = setInterval(() => {
        setTick((n) => n + 1)
        setNow(Date.now())
        if (Date.now() - clock() >= 1000) setClock(Date.now())
      }, 80)
      ;(timer as { unref?: () => void }).unref?.()        // never what keeps a window alive
    } else if (!running && timer) {
      clearInterval(timer)
      timer = undefined
    }
  })
  onCleanup(() => timer && clearInterval(timer))
  // The background list's ages: once a second while anything is out.
  const [bgNow, setBgNow] = createSignal(Date.now())
  let bgTimer: ReturnType<typeof setInterval> | undefined
  createEffect(() => {
    const out = (props.session()?.background?.length ?? 0) > 0
    if (out && !bgTimer) {
      setBgNow(Date.now())
      bgTimer = setInterval(() => setBgNow(Date.now()), 1000)
      ;(bgTimer as { unref?: () => void }).unref?.()
    } else if (!out && bgTimer) {
      clearInterval(bgTimer)
      bgTimer = undefined
    }
  })
  onCleanup(() => bgTimer && clearInterval(bgTimer))

  const style = themedStyle(props.renderer, t)
  const asked = () => props.session()?.pending[0]
  // What the footer draws in: the terminal's width less the gutter at the
  // left; every line is cut a column short of the right edge as well.
  const cols = () => Math.max(10, term().width - GUTTER)
  const rule = () => "─".repeat(Math.max(1, cols() - 1))
  // The live blocks by key (strings: <For> keeps an element per key, so a
  // block keeps its element from its first delta until it leaves).
  const live = createMemo(() => liveBlocks(props.session(), cols()))
  const liveKeys = createMemo(() => live().map((b) => b.key), [],
                              { equals: (a, b) => a.length === b.length && a.every((k, i) => k === b[i]) })
  const liveText = (key: string) => live().find((b) => b.key === key)?.text ?? ""
  const liveHeight = () => liveRows(props.session(), cols())
  /** The prompt box's top border, with the session's name on it at the
   *  right (/rename, -n); not in the agent view,
   *  whose box is for another session. */
  const topRule = () => {
    const s = props.session()
    const name = viewOpen() ? "" : (s?.title || s?.summary?.title || "")
    const total = Math.max(1, cols() - 1)
    const chip = name && total >= 12 ? ` ${fit(name, Math.max(4, Math.floor(total / 2)))} ` : ""
    // Shell mode: at the left, what the line is and how to leave it -- as
    // much as fits beside the name ("──" before it, one "─" after).
    const said = props.composer.shell() ? shellHint(total - width(chip) - (chip ? 1 : 0) - 5) : ""
    const head = said ? "── " : ""
    const words = said ? `${said} ` : ""
    const used = width(head) + width(words) + width(chip) + (chip ? 1 : 0)
    return { head, mode: said ? "shell" : "", words: words.slice(said ? 5 : 0),
             left: "─".repeat(Math.max(said ? 0 : 1, total - used)), chip, right: chip ? "─" : "" }
  }
  /** The box's rules: in the shell's colour while the line is a command. */
  const ruleTone = () => props.composer.shell() ? t().shell : t().faint

  /** The task list's lines: none in the agent view, which is another
   *  session's place. */
  const tasks = createMemo(() => viewOpen() ? []
    : taskLines(props.session()?.tasks ?? [], Boolean(props.tasksFolded?.()), cols(), term().height,
                props.keymap.label("tasks") || "ctrl+t"))
  /** What is still out: not in the agent view either. */
  const background = createMemo(() => viewOpen() ? []
    : backgroundLines(props.session()?.background ?? [], cols(), bgNow()))
  /** The rows the two lists above the prompt take. */
  const listRows = () => tasks().length + background().length
  const full = () => Boolean(props.screen)
  createEffect(() => {
    // Full screen has no footer of a height: the window is all of it.
    if (full()) return
    props.renderer.footerHeight = viewOpen() && !props.agentsOnOwnScreen
      ? agentFooterRows(term())
      : footerRows(props.session(), term(), note(), props.composer.height(),
                   props.approval, props.picker?.height() ?? 0, listRows())
  })
  /** The agent view's list: its rows less the message, the composer and
   *  its rules, and the hints. */
  const panelRows = () => Math.max(2, (props.agentsOnOwnScreen || full() ? term().height : agentFooterRows(term()))
    - (note() ? 1 : 0) - 3 - props.composer.height())

  const status = (turn: TurnLive) => {
    const s = props.session()
    const parts = [elapsed(turn.startedAt, now())]
    const parts2 = s ? turnUsage(turn, s.usage) : []
    parts.push(...parts2, `${props.keymap.label("interrupt")} to interrupt`)
    const doing = s?.wait || turn.label
    const label = fit(doing, Math.max(8, cols() - 4 - parts.join(" · ").length - 4))
    return { label, rest: fit(`  ${parts.join(" · ")}`, Math.max(0, cols() - 3 - label.length)) }
  }

  // A wrapped draft takes another number of rows at another width -- known
  // once the textarea is laid out at it, a frame later.
  createEffect(() => {
    term()
    props.composer.remeasure()
    for (const ms of [40, 150]) (setTimeout(() => props.composer.remeasure(), ms) as { unref?: () => void }).unref?.()
  })

  const budget = () => approvalBudget(props.session(), term(), note(), props.composer.height(),
                                      props.picker?.height() ?? 0, listRows())

  // Nothing in the footer shrinks: what does not fit is cut at the top
  // (the composer and the keys are at the bottom), never drawn over.
  const placeholder = () => {
    const v = props.agents
    if (props.swarm?.state.open) {
      if (props.swarm.state.pendingAnswer) return "the answer…"
      if (props.swarm.state.pendingStart) return "y to start · n lets it go"
      return "a new swarm's goal… (or /swarm again to close)"
    }
    if (v?.state.open) {
      if (v.state.renaming) return "the new name…"
      const peeked = v.state.peek ? v.rows().find((r) => r.id === v.state.peek!.id) : undefined
      if (peeked?.kind === "attention") {
        return takesAnswer(peeked.item) ? `the answer for ${peeked.name}…`
          : `${props.keymap.label("agentsOpen")} opens it in a session`
      }
      if (peeked) return peeked.kind === "elsewhere" ? `a message for ${peeked.name}…` : `reply to ${peeked.name}…`
      return "a new session: what should it do?"
    }
    if (props.composer.shell()) {
      return props.session()?.turn ? "a shell command… (runs when the turn ends)" : "a shell command…"
    }
    // A question open: the prompt keeps the keys (app.tsx), and says what
    // a line does here.
    const a = asked()
    if (a?.kind === "ask") return "type your own answer… (or pick above)"
    if (a && props.approval.state.mode === "reason" && props.approval.state.requestId === a.requestId) {
      return "the reason… (enter refuses with it)"
    }
    if (a) return "a message… (the keys above answer on an empty prompt)"
    return props.starting?.() ? "starting the session…"
      : props.session()?.turn
      ? `queue a message… (${props.keymap.label("steer")} sends it into the turn)` : "ask DELFIN…"
  }

  return (
    <box flexDirection="column" width="100%" height="100%" justifyContent="flex-end" overflow="hidden"
         paddingLeft={GUTTER}>
      <Show when={viewOpen()}>
        <Show when={props.swarm?.state.open} fallback={
          <AgentViewPanel view={props.agents!} theme={t} width={cols} height={panelRows} tick={tick} now={clock}
                          keymap={props.keymap} />
        }>
          <SwarmPanel view={props.swarm!} theme={t} width={cols} height={panelRows} />
        </Show>
      </Show>
      {/* Mounted for as long as the window is full screen -- hidden while
          the agent view has the screen, so closing the view keeps its
          scroll position and does not lay out every row again. */}
      <Show when={full()}>
        <TranscriptView renderer={props.renderer} screen={props.screen!} theme={t} width={cols} height={() => term().height}
                        visible={() => !viewOpen()} ref={props.scrollRef} onUserScroll={props.onUserScroll} />
      </Show>
      <Show when={note()}>
        <text flexShrink={0} fg={props.offer?.() ? t().err : t().warn} wrapMode="none"
              attributes={props.offer?.() ? TextAttributes.BOLD : 0} ref={offerHit.ref}
              onMouseDown={offerHit.onMouseDown} onMouseUp={offerHit.onMouseUp}
              onMouseDrag={offerHit.onMouseDrag}>{fit(note(), cols() - 1)}</text>
      </Show>
      <Show when={!viewOpen()}>
      {/* The answer's unfinished block, under the rows in both screens: a
          delta lays out this box only, never the whole conversation. A
          block that finishes stays until its rows are drawn, with the
          element that drew it while it streamed -- a new one shows
          nothing until highlighted -- and goes in the frame they come. */}
      <Show when={liveHeight() > 0}>
        {/* Its last rows, at the bottom: what does not fit is cut at the top. */}
        <box flexShrink={0} height={liveHeight()} flexDirection="column" justifyContent="flex-end" overflow="hidden">
          <For each={liveKeys()}>
            {(key, i) => <LiveBlockView text={() => liveText(key)} first={() => i() === 0} syntax={style} width={cols} />}
          </For>
        </box>
      </Show>
      <Show when={!props.session()?.turn && props.session()?.wait}>
        <text flexShrink={0} wrapMode="none" fg={t().warn}>{fit(`… ${props.session()?.wait ?? ""}`, cols() - 1)}</text>
      </Show>
      <Show when={props.session()?.turn}>
        {(turn: () => TurnLive) => (
          <text flexShrink={0} wrapMode="none">
            <span style={{ fg: t().accent }}>{SPINNER[tick() % SPINNER.length] + " "}</span>
            <span style={{ fg: t().text }}>{status(turn()).label}</span>
            <span style={{ fg: t().faint }}>{status(turn()).rest}</span>
          </text>
        )}
      </Show>
      <TasksPanel lines={tasks} theme={t} />
      <BackgroundPanel lines={background} theme={t} />
      <Show when={asked()} keyed>
        {(a: Asked) => <ApprovalPanel asked={a} ui={props.approval} theme={t} width={cols} maxRows={budget} />}
      </Show>
      <Show when={props.picker?.state.open && !asked()}>
        <PickerView picker={props.picker!} theme={t} width={cols} />
      </Show>
      </Show>
      <text flexShrink={0} wrapMode="none">
        <span style={{ fg: ruleTone() }}>{topRule().head}</span>
        <span style={{ fg: t().shell, bold: true }}>{topRule().mode}</span>
        <span style={{ fg: t().muted }}>{topRule().words}</span>
        <span style={{ fg: ruleTone() }}>{topRule().left}</span>
        <span style={{ fg: t().accentSoft }}>{topRule().chip}</span>
        <span style={{ fg: ruleTone() }}>{topRule().right}</span>
      </text>
      <ComposerView
        composer={props.composer}
        theme={t}
        width={cols}
        mark={() => {
          // The prompt marker in the posture's colour (plan, accept edits).
          const look = postureLook(viewOpen() ? "" : props.session()?.approval ?? "", t())
          return look.marked ? look.color : t().accent
        }}
        focused={props.focused}
        pasteActive={props.focused}
        placeholder={placeholder}
      />
      <text flexShrink={0} fg={ruleTone()} wrapMode="none">{rule()}</text>
      <Hints {...props} cols={cols} />
    </box>
  )
}

/**
 * One live block. Markdown, not drawn `streaming` -- OpenTUI's preview of
 * a streaming block is its raw text ("## " before a heading) until each
 * delta's highlighting lands; this keeps the last highlighted text. A
 * fenced code block is drawn as code that streams: as part of markdown it
 * would lose its colours for a frame on every delta.
 */
function LiveBlockView(props: { text: () => string; first: () => boolean; syntax: () => SyntaxStyle;
                         width: () => number }): JSX.Element {
  const fence = () => fenceOf(props.text())
  // Drawn without a preview, a block shows nothing until its colours come:
  // if they are late (a busy highlighter), plain text meanwhile.
  let element: Renderable | undefined
  const take = (r: Renderable) => { element = r }
  createEffect(on(props.text, () => queueMicrotask(() => {
    if (element && !element.isDestroyed) plainIfLate(element)
  })))
  return (
    // Content last: an element builds its blocks when it is set, with the
    // style it has by then (props go in in order).
    <Show when={fence()} fallback={
      // A block that holds a table: in segments, as the rows draw it.
      <Show when={holdsTable(props.text())} fallback={
        <markdown ref={take} syntaxStyle={props.syntax()} conceal width="100%" flexShrink={0}
                  marginTop={props.first() ? 0 : 1} tableOptions={{ widthMode: "content" }} content={props.text()} />
      }>
        <SegmentedMarkdown ref={take} text={props.text} width={props.width} syntax={props.syntax}
                           marginTop={() => props.first() ? 0 : 1} />
      </Show>
    }>
      <code ref={take} syntaxStyle={props.syntax()} streaming drawUnstyledText={false} conceal={false} width="100%"
            flexShrink={0} marginTop={props.first() ? 0 : 1} filetype={infoStringToFiletype(fence()?.lang ?? "")}
            content={fence()?.code ?? ""} />
    </Show>
  )
}

/** The right end of the hints line: the marks in their colours, the rest
 *  faint -- `◉ 3 · ✗ 1`, or a new item's `◉ relax Fe complex asks`. */
export function attentionCorner(a: { asks: number; failed: number; note: string } | undefined,
                                room: number): Corner {
  const none = { mark: "", text: "", sep: "", mark2: "", text2: "", note: false }
  if (!a || room < 6) return none
  if (a.note) {
    const mark = a.note.slice(0, 1)
    return { ...none, mark, text: fit(a.note.slice(1), Math.max(0, room - 1)), note: true }
  }
  const parts = [a.asks ? { mark: "◉", text: ` ${a.asks}` } : null, a.failed ? { mark: "✗", text: ` ${a.failed}` } : null]
    .filter(Boolean) as { mark: string; text: string }[]
  if (!parts.length) return none
  const [first, second] = parts
  return { mark: first!.mark, text: first!.text, sep: second ? " · " : "", mark2: second?.mark ?? "", text2: second?.text ?? "",
           note: false }
}

export interface Corner {
  mark: string
  text: string
  sep: string
  mark2: string
  text2: string
  /** A new item's note (its words muted), not the count (faint). */
  note: boolean
}

function Hints(props: FooterProps & { cols: () => number }): JSX.Element {
  const t = props.theme
  /** What needs the user, at the right end (not while the agent view is
   *  open: its own head says it). At most half the line. */
  const corner = createMemo(() => props.agents?.state.open ? attentionCorner(undefined, 0)
    : attentionCorner(props.attention?.(), Math.floor((props.cols() - 1) / 2)))
  const cornerWidth = () => {
    const c = corner()
    return width(c.mark + c.text + c.sep + c.mark2 + c.text2)
  }
  /** The line in three parts -- what is wrong first, then the posture in
   *  its colour, then the model and the keys -- cut from the right, then
   *  spaces up to the corner. */
  const line = createMemo(() => {
    const room = Math.max(0, props.cols() - 1 - (cornerWidth() ? cornerWidth() + 2 : 0))
    if (props.agents?.state.open) {
      return { lead: fit("  " + agentHelp(props.agents, props.keymap, props.cols() - 3), room), posture: null, rest: "",
               agents: "" }
    }
    const s = props.session()
    const lead: string[] = []
    const conn = props.connection()
    if (conn !== "up") lead.push(conn === "reconnecting" ? "reconnecting to the agent service…" : conn)
    if (props.starting?.()) lead.push("starting the session…")
    if (s?.state === "ended") lead.push("session ended")
    const rest: string[] = []
    if (s?.auto) rest.push(s.auto)
    if (s?.model) rest.push(s.model)
    if (s?.usage.costUsd && !s.turn) rest.push(`$${s.usage.costUsd.toFixed(4)}`)
    rest.push(`${props.keymap.label("cycleMode")} mode`)
    // Last: the agent view's key, with the other sessions -- in as many
    // words as the room left for it takes.
    const agentsIn = (used: number) => othersHint(props.others?.(), props.keymap.label("agents"),
                                                  Math.max(0, room - used - 3))
    const look = s?.approval ? postureLook(s.approval, t()) : null
    const head = "  " + lead.map((part) => `${part} · `).join("")
    if (!look) {
      const before = head + rest.join(" · ")
      const agents = agentsIn(width(before))
      return { lead: fit(`${before} · ${agents}`, room), posture: null, rest: "", agents }
    }
    // The posture keeps its room: a note before it is cut first.
    const postureText = fit(look.text, Math.max(0, room - 2))
    const leadText = fit(head, Math.max(0, room - width(postureText)))
    const before = ` · ${rest.join(" · ")}`
    const agents = agentsIn(width(leadText) + width(postureText) + width(before))
    const restText = fit(`${before} · ${agents}`, Math.max(0, room - width(leadText) - width(postureText)))
    return { lead: leadText, posture: postureText ? { ...look, text: postureText } : null, rest: restText, agents }
  })
  /** A click on the agents part (full screen): the agent view, as its key. */
  const agentsHit = clickable((col) => {
    const l = line()
    if (!l.agents || !props.onOthersClick) return null
    const shown = l.lead + (l.posture?.text ?? "") + l.rest
    const at = shown.lastIndexOf(l.agents)
    if (at < 0) return null
    const from = width(shown.slice(0, at))
    return col >= from && col < from + width(l.agents) ? "agents" : null
  }, () => props.onOthersClick?.())
  const gap = () => {
    if (!cornerWidth()) return ""
    const used = width(line().lead) + width(line().posture?.text ?? "") + width(line().rest)
    return " ".repeat(Math.max(2, props.cols() - 1 - used - cornerWidth()))
  }
  const markTone = (mark: string) => mark === "✗" ? t().err : t().warn
  return (
    <text flexShrink={0} wrapMode="none" ref={agentsHit.ref} onMouseDown={agentsHit.onMouseDown}
          onMouseUp={agentsHit.onMouseUp} onMouseDrag={agentsHit.onMouseDrag}>
      <span style={{ fg: t().faint }}>{line().lead}</span>
      <span style={{ fg: line().posture?.color ?? t().faint, bold: Boolean(line().posture?.marked) }}>
        {line().posture?.text ?? ""}
      </span>
      <span style={{ fg: t().faint }}>{line().rest}</span>
      <span style={{ fg: t().faint }}>{gap()}</span>
      <span style={{ fg: markTone(corner().mark) }}>{corner().mark}</span>
      <span style={{ fg: corner().note ? t().muted : t().faint }}>{corner().text}</span>
      <span style={{ fg: t().faint }}>{corner().sep}</span>
      <span style={{ fg: markTone(corner().mark2) }}>{corner().mark2}</span>
      <span style={{ fg: t().faint }}>{corner().text2}</span>
    </text>
  )
}
