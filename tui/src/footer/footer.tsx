/**
 * The live area under the scrollback: what the session in view is doing
 * now, the question it waits on, and the line the user types into.
 *
 *   message   a short note of the window's own (one row)
 *   tail      the answer's unfinished block, as markdown (≤ 6 rows)
 *   status    spinner · what it does · elapsed · tokens · cost · esc
 *   approval  the open question, with exactly the keys the session offers
 *   composer  ❯ the prompt (grows to 8 rows; completion and search above)
 *   hints     posture · model · keys (◉ when another session waits)
 *
 * The agent view (←) takes the whole terminal but its last row: its list
 * in place of the tail, status and question, the composer under it.
 *
 * Every line is cut to the terminal's width, so the footer knows exactly
 * how many rows it takes (`footerRows`), at any width: a footer taller
 * than its content wastes rows, a shorter one draws over itself.
 */
import { SyntaxStyle, type CliRenderer } from "@opentui/core"
import { Show, createEffect, createMemo, createSignal, onCleanup, type Accessor, type JSX } from "solid-js"
import { AgentViewPanel, agentHelp, type AgentView } from "../agents/view"
import { ComposerView, type Composer } from "../composer/composer"
import type { Keymap } from "../keymap/keymap"
import type { Asked, SessionState, TurnLive } from "../store/state"
import { GUTTER, SPINNER, markdownStyles, postureLook, type Theme } from "../theme/theme"
import { fit, width } from "../text"
import { ApprovalPanel, approvalRows, type ApprovalUi } from "./approval"
import { PickerView, type Picker } from "./picker"

export const TAIL_MAX = 6

export interface FooterProps {
  renderer: CliRenderer
  theme: () => Theme
  keymap: Keymap
  session: () => SessionState | null
  connection: () => string
  /** Starting the first session: said in the hints line. */
  starting?: () => boolean
  message: () => string
  composer: Composer
  approval: ApprovalUi
  picker?: Picker
  /** The agent view; while open it has the footer. */
  agents?: AgentView
  /** The view is drawn on the alternate screen: all of it is the view's,
   *  and the footer keeps its height for the way back. */
  agentsOnOwnScreen?: boolean
  /** Sessions other than the one in view that wait for the user. */
  othersWaiting?: () => number
  focused: () => boolean
}

/** Rows the footer takes with the agent view open: all but one. */
export function agentFooterRows(term: { height: number }): number {
  return Math.max(4, term.height - 1)
}

function tailRows(tail: string): number {
  return Math.min(TAIL_MAX, tail.trimEnd().split("\n").filter((l) => l.trim()).length)
}

/** The rows everything but the approval panel takes. */
function fixedRows(s: SessionState | null, message: string, composerRows: number, pickerRows: number): number {
  return 3 + composerRows + pickerRows + (s?.turn || s?.wait ? 1 : 0) + (message ? 1 : 0)
}

/** What the approval panel may take: the rest of the terminal (less one
 *  row: the footer never fills it), so its keys are always on screen. */
export function approvalBudget(s: SessionState | null, term: { height: number }, message: string,
                               composerRows = 1, pickerRows = 0): number {
  return Math.max(3, term.height - 1 - fixedRows(s, message, composerRows, pickerRows))
}

/** Rows the footer needs; never more than the terminal has (less one). */
export function footerRows(s: SessionState | null, term: { width: number; height: number }, message: string,
                           composerRows = 1, approval: ApprovalUi | null = null, pickerRows = 0): number {
  let rows = fixedRows(s, message, composerRows, pickerRows)
  const asked = s?.pending[0]
  // With a question open the turn waits: its tail is empty (a question
  // closes the text part), and the panel has the room.
  if (s && !asked) rows += tailRows(s.tail)
  rows += approvalRows(asked, approval?.state ?? null, term.width - GUTTER,
                       approvalBudget(s, term, message, composerRows, pickerRows))
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
  const term = terminalSize(props.renderer)
  const [tick, setTick] = createSignal(0)
  const [now, setNow] = createSignal(Date.now())
  // The agent view's ages: once a second is enough (its rows are drawn
  // again on each change of it).
  const [clock, setClock] = createSignal(Date.now())
  const viewOpen = () => Boolean(props.agents?.state.open)
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

  const style = createMemo(() => SyntaxStyle.fromStyles(markdownStyles(t()), props.renderer.nativeScene!))
  const asked = () => props.session()?.pending[0]
  // What the footer draws in: the terminal's width less the gutter at the
  // left; every line is cut a column short of the right edge as well.
  const cols = () => Math.max(10, term().width - GUTTER)
  const rule = () => "─".repeat(Math.max(1, cols() - 1))
  /** The prompt box's top border, with the session's name on it at the
   *  right (/rename, -n), as Claude Code shows it; not in the agent view,
   *  whose box is for another session. */
  const topRule = () => {
    const s = props.session()
    const name = viewOpen() ? "" : (s?.title || s?.summary?.title || "")
    const total = Math.max(1, cols() - 1)
    if (!name || total < 12) return { left: rule(), chip: "", right: "" }
    const chip = ` ${fit(name, Math.max(4, Math.floor(total / 2)))} `
    return { left: "─".repeat(Math.max(1, total - width(chip) - 1)), chip, right: "─" }
  }

  createEffect(() => {
    props.renderer.footerHeight = viewOpen() && !props.agentsOnOwnScreen
      ? agentFooterRows(term())
      : footerRows(props.session(), term(), props.message(), props.composer.height(),
                   props.approval, props.picker?.height() ?? 0)
  })
  /** The agent view's list: its rows less the message, the composer and
   *  its rules, and the hints. */
  const panelRows = () => Math.max(2, (props.agentsOnOwnScreen ? term().height : agentFooterRows(term()))
    - (props.message() ? 1 : 0) - 3 - props.composer.height())

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

  const budget = () => approvalBudget(props.session(), term(), props.message(), props.composer.height(),
                                      props.picker?.height() ?? 0)

  // Nothing in the footer shrinks: what does not fit is cut at the top
  // (the composer and the keys are at the bottom), never drawn over.
  const placeholder = () => {
    const v = props.agents
    if (v?.state.open) {
      if (v.state.renaming) return "the new name…"
      const peeked = v.state.peek ? v.rows().find((r) => r.id === v.state.peek!.id) : undefined
      if (peeked) return peeked.kind === "elsewhere" ? `a message for ${peeked.name}…` : `reply to ${peeked.name}…`
      return "a new session: what should it do?"
    }
    return props.starting?.() ? "starting the session…"
      : props.session()?.turn
      ? `queue a message… (${props.keymap.label("steer")} sends it into the turn)` : "ask DELFIN…"
  }

  return (
    <box flexDirection="column" width="100%" height="100%" justifyContent="flex-end" overflow="hidden"
         paddingLeft={GUTTER}>
      <Show when={viewOpen()}>
        <AgentViewPanel view={props.agents!} theme={t} width={cols} height={panelRows} tick={tick} now={clock}
                        keymap={props.keymap} />
      </Show>
      <Show when={props.message()}>
        <text flexShrink={0} fg={t().warn} wrapMode="none">{fit(props.message(), cols() - 1)}</text>
      </Show>
      <Show when={!viewOpen()}>
      <Show when={tailRows(props.session()?.tail ?? "") > 0 && !asked()}>
        <scrollbox flexShrink={0} height={tailRows(props.session()?.tail ?? "")} stickyScroll stickyStart="bottom">
          <markdown content={props.session()?.tail ?? ""} streaming syntaxStyle={style()} conceal width="100%" />
        </scrollbox>
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
      <Show when={asked()} keyed>
        {(a: Asked) => <ApprovalPanel asked={a} ui={props.approval} theme={t} width={cols} maxRows={budget} />}
      </Show>
      <Show when={props.picker?.state.open && !asked()}>
        <PickerView picker={props.picker!} theme={t} width={cols} />
      </Show>
      </Show>
      <text flexShrink={0} wrapMode="none">
        <span style={{ fg: t().faint }}>{topRule().left}</span>
        <span style={{ fg: t().accentSoft }}>{topRule().chip}</span>
        <span style={{ fg: t().faint }}>{topRule().right}</span>
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
        focused={() => props.focused() && (viewOpen() || !asked())}
        pasteActive={() => props.focused() && (viewOpen() || !asked())}
        placeholder={placeholder}
      />
      <text flexShrink={0} fg={t().faint} wrapMode="none">{rule()}</text>
      <Hints {...props} cols={cols} />
    </box>
  )
}

function Hints(props: FooterProps & { cols: () => number }): JSX.Element {
  const t = props.theme
  // Another session waits for the user: said first, where the eye is.
  const badge = () => {
    if (props.agents?.state.open) return ""
    const n = props.othersWaiting?.() ?? 0
    return n ? `  ◉ ${n} waiting` : ""
  }
  /** The line in three parts -- what is wrong first, then the posture in
   *  its colour, then the model and the keys -- cut from the right. */
  const line = createMemo(() => {
    const room = Math.max(0, props.cols() - 1 - width(badge()))
    if (props.agents?.state.open) {
      return { lead: fit("  " + agentHelp(props.agents, props.keymap, props.cols() - 3), room), posture: null, rest: "" }
    }
    const s = props.session()
    const lead: string[] = []
    const conn = props.connection()
    if (conn !== "up") lead.push(conn === "reconnecting" ? "reconnecting to the agent service…" : conn)
    if (props.starting?.()) lead.push("starting the session…")
    if (s?.state === "ended") lead.push("session ended")
    const rest: string[] = []
    if (s?.model) rest.push(s.model)
    if (s?.usage.costUsd && !s.turn) rest.push(`$${s.usage.costUsd.toFixed(4)}`)
    rest.push(`${props.keymap.label("cycleMode")} mode`, `${props.keymap.label("agents")} agents`)
    const look = s?.approval ? postureLook(s.approval, t()) : null
    const head = (badge() ? "  · " : "  ") + lead.map((part) => `${part} · `).join("")
    if (!look) return { lead: fit(head + rest.join(" · "), room), posture: null, rest: "" }
    // The posture keeps its room: a note before it is cut first.
    const postureText = fit(look.text, Math.max(0, room - 2))
    const leadText = fit(head, Math.max(0, room - width(postureText)))
    const restText = fit(rest.length ? ` · ${rest.join(" · ")}` : "", Math.max(0, room - width(leadText) - width(postureText)))
    return { lead: leadText, posture: postureText ? { ...look, text: postureText } : null, rest: restText }
  })
  return (
    <text flexShrink={0} wrapMode="none">
      <span style={{ fg: t().warn }}>{badge()}</span>
      <span style={{ fg: t().faint }}>{line().lead}</span>
      <span style={{ fg: line().posture?.color ?? t().faint, bold: Boolean(line().posture?.marked) }}>
        {line().posture?.text ?? ""}
      </span>
      <span style={{ fg: t().faint }}>{line().rest}</span>
    </text>
  )
}
