/**
 * The live area under the scrollback: what the session in view is doing
 * now, the question it waits on, and the line the user types into.
 *
 *   message   a short note of the window's own (one row)
 *   tail      the answer's unfinished block, as markdown (≤ 6 rows)
 *   status    spinner · what it does · elapsed · tokens · cost · esc
 *   approval  the open question, with exactly the keys the session offers
 *   composer  ❯ the prompt (grows to 8 rows; completion and search above)
 *   hints     posture · model · keys
 *
 * Every line is cut to the terminal's width, so the footer knows exactly
 * how many rows it takes (`footerRows`), at any width: a footer taller
 * than its content wastes rows, a shorter one draws over itself.
 */
import { SyntaxStyle, type CliRenderer } from "@opentui/core"
import { Show, createEffect, createMemo, createSignal, onCleanup, type Accessor, type JSX } from "solid-js"
import { ComposerView, type Composer } from "../composer/composer"
import type { Keymap } from "../keymap/keymap"
import type { Asked, SessionState, TurnLive } from "../store/state"
import { SPINNER, markdownStyles, type Theme } from "../theme/theme"
import { fit } from "../text"
import { ApprovalPanel, approvalRows, type ApprovalUi } from "./approval"

export const TAIL_MAX = 6

export interface FooterProps {
  renderer: CliRenderer
  theme: () => Theme
  keymap: Keymap
  session: () => SessionState | null
  connection: () => string
  message: () => string
  composer: Composer
  approval: ApprovalUi
  focused: () => boolean
}

function tailRows(tail: string): number {
  return Math.min(TAIL_MAX, tail.trimEnd().split("\n").filter((l) => l.trim()).length)
}

/** Rows the footer needs; never more than the terminal has (less one). */
export function footerRows(s: SessionState | null, term: { width: number; height: number }, message: string,
                           composerRows = 1, approval: ApprovalUi | null = null): number {
  let rows = 3 + composerRows                    // rule, composer, rule, hints
  if (s?.turn || s?.wait) rows += 1
  if (s) rows += tailRows(s.tail)
  rows += approvalRows(s?.pending[0], approval?.state ?? null, term.width)
  if (message) rows += 1
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
  let timer: ReturnType<typeof setInterval> | undefined
  createEffect(() => {
    const running = Boolean(props.session()?.turn)
    if (running && !timer) {
      timer = setInterval(() => {
        setTick((n) => n + 1)
        setNow(Date.now())
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
  const cols = () => Math.max(10, term().width)
  const rule = () => "─".repeat(Math.max(1, cols() - 1))

  createEffect(() => {
    props.renderer.footerHeight = footerRows(props.session(), term(), props.message(), props.composer.height(),
                                             props.approval)
  })

  const status = (turn: TurnLive) => {
    const s = props.session()
    const parts = [elapsed(turn.startedAt, now())]
    const parts2 = s ? turnUsage(turn, s.usage) : []
    parts.push(...parts2, `${props.keymap.label("interrupt")} to interrupt`)
    const doing = s?.wait || turn.label
    const label = fit(doing, Math.max(8, cols() - 4 - parts.join(" · ").length - 4))
    return { label, rest: fit(`  ${parts.join(" · ")}`, Math.max(0, cols() - 3 - label.length)) }
  }

  return (
    <box flexDirection="column" width="100%" height="100%" justifyContent="flex-end">
      <Show when={props.message()}>
        <text fg={t().warn} wrapMode="none">{fit(props.message(), cols() - 1)}</text>
      </Show>
      <Show when={tailRows(props.session()?.tail ?? "") > 0}>
        <scrollbox height={tailRows(props.session()?.tail ?? "")} stickyScroll stickyStart="bottom">
          <markdown content={props.session()?.tail ?? ""} streaming syntaxStyle={style()} conceal width="100%" />
        </scrollbox>
      </Show>
      <Show when={!props.session()?.turn && props.session()?.wait}>
        <text wrapMode="none" fg={t().warn}>{fit(`… ${props.session()?.wait ?? ""}`, cols() - 1)}</text>
      </Show>
      <Show when={props.session()?.turn}>
        {(turn: () => TurnLive) => (
          <text wrapMode="none">
            <span style={{ fg: t().accent }}>{SPINNER[tick() % SPINNER.length] + " "}</span>
            <span style={{ fg: t().text }}>{status(turn()).label}</span>
            <span style={{ fg: t().faint }}>{status(turn()).rest}</span>
          </text>
        )}
      </Show>
      <Show when={asked()} keyed>
        {(a: Asked) => <ApprovalPanel asked={a} ui={props.approval} theme={t} width={cols} />}
      </Show>
      <text fg={t().faint} wrapMode="none">{rule()}</text>
      <ComposerView
        composer={props.composer}
        theme={t}
        focused={() => props.focused() && !asked()}
        pasteActive={() => props.focused() && !asked()}
        placeholder={() => props.session()?.turn
          ? `queue a message… (${props.keymap.label("steer")} sends it into the turn)` : "ask DELFIN…"}
      />
      <text fg={t().faint} wrapMode="none">{rule()}</text>
      <Hints {...props} cols={cols} />
    </box>
  )
}

function Hints(props: FooterProps & { cols: () => number }): JSX.Element {
  const t = props.theme
  const text = () => {
    const s = props.session()
    const parts: string[] = []
    const conn = props.connection()
    if (conn !== "up") parts.push(conn === "reconnecting" ? "reconnecting to the agent service…" : conn)
    if (s?.state === "ended") parts.push("session ended")
    if (s?.approval) parts.push(s.approval)
    if (s?.model) parts.push(s.model)
    if (s?.usage.costUsd && !s.turn) parts.push(`$${s.usage.costUsd.toFixed(4)}`)
    parts.push(`${props.keymap.label("cycleMode")} mode`, `${props.keymap.label("agents")} agents`)
    return fit("  " + parts.join(" · "), props.cols() - 1)
  }
  return <text fg={t().faint} wrapMode="none">{text()}</text>
}
