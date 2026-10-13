/**
 * The live area under the scrollback: what the session in view is doing
 * now, the question it waits on, and the line the user types into.
 *
 *   tail      the answer's unfinished block, as markdown (≤ 6 rows)
 *   status    spinner · what it does · elapsed · model · esc to interrupt
 *   approval  the open question, with exactly the keys the session offers
 *   composer  ❯ the prompt
 *   hints     posture · keys
 *
 * The footer sizes itself to its content (FooterSize); everything it shows
 * comes from the store, so it redraws only when the store changes.
 */
import { SyntaxStyle, TextAttributes, type CliRenderer, type InputRenderable } from "@opentui/core"
import { For, Show, createEffect, createMemo, createSignal, onCleanup, type JSX } from "solid-js"
import type { Keymap } from "../keymap/keymap"
import type { Asked, SessionState, TurnLive } from "../store/state"
import { SPINNER, markdownStyles, type Theme } from "../theme/theme"

export const TAIL_MAX = 6
export const PREVIEW_MAX = 12

export interface FooterProps {
  renderer: CliRenderer
  theme: () => Theme
  keymap: Keymap
  session: () => SessionState | null
  connection: () => string
  message: () => string
  onSubmit: (text: string) => void
  inputRef?: (input: InputRenderable) => void
  focused: () => boolean
}

function tailRows(tail: string): number {
  return Math.min(TAIL_MAX, tail.trimEnd().split("\n").filter((l) => l.trim()).length)
}

function previewLines(asked: Asked): string[] {
  const lines = String(asked.preview ?? "").split("\n")
  if (lines.length <= PREVIEW_MAX) return lines
  return [...lines.slice(0, PREVIEW_MAX - 1), `… ${lines.length - PREVIEW_MAX + 1} more lines`]
}

/** How many rows the approval panel takes. */
export function approvalRows(asked: Asked | undefined): number {
  if (!asked) return 0
  if (asked.kind === "ask") return 3 + (asked.choices?.length ?? 0) + (asked.question ? 1 : 0)
  return 3 + previewLines(asked).length
}

export function footerRows(s: SessionState | null, termRows: number, message: string): number {
  let rows = 4                                   // rule, composer, rule, hints
  if (s?.turn) rows += 1
  if (s) rows += tailRows(s.tail)
  rows += approvalRows(s?.pending[0])
  if (message) rows += 1
  return Math.max(4, Math.min(rows, Math.max(4, termRows - 1)))
}

function elapsed(startedAt: number, now: number): string {
  const secs = Math.max(0, (now - startedAt) / 1000)
  return secs < 60 ? `${secs.toFixed(1)}s` : `${Math.floor(secs / 60)}m${String(Math.floor(secs % 60)).padStart(2, "0")}s`
}

export function Footer(props: FooterProps): JSX.Element {
  const t = props.theme
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
    } else if (!running && timer) {
      clearInterval(timer)
      timer = undefined
    }
  })
  onCleanup(() => timer && clearInterval(timer))

  const style = createMemo(() => SyntaxStyle.fromStyles(markdownStyles(t()), props.renderer.nativeScene!))
  const asked = () => props.session()?.pending[0]
  const rule = () => "─".repeat(Math.max(1, props.renderer.terminalWidth - 1))

  createEffect(() => {
    props.renderer.footerHeight = footerRows(props.session(), props.renderer.terminalHeight, props.message())
  })

  return (
    <box flexDirection="column" width="100%" height="100%" justifyContent="flex-end">
      <Show when={props.message()}>
        <text fg={t().warn}>{props.message()}</text>
      </Show>
      <Show when={tailRows(props.session()?.tail ?? "") > 0}>
        <scrollbox height={tailRows(props.session()?.tail ?? "")} stickyScroll stickyStart="bottom">
          <markdown content={props.session()?.tail ?? ""} streaming syntaxStyle={style()} conceal width="100%" />
        </scrollbox>
      </Show>
      <Show when={props.session()?.turn}>
        {(turn: () => TurnLive) => (
          <text>
            <span style={{ fg: t().accent }}>{SPINNER[tick() % SPINNER.length] + " "}</span>
            <span style={{ fg: t().text }}>{turn().label}</span>
            <span style={{ fg: t().faint }}>
              {`  ${elapsed(turn().startedAt, now())} · ${props.session()?.model ?? ""} · ${props.keymap.label("interrupt")} to interrupt`}
            </span>
          </text>
        )}
      </Show>
      <Show when={asked()}>{(a: () => Asked) => <Approval asked={a()} theme={t} />}</Show>
      <text fg={t().faint}>{rule()}</text>
      <box flexDirection="row">
        <text fg={t().accent} attributes={TextAttributes.BOLD}>{"❯ "}</text>
        <input
          ref={(el: InputRenderable) => props.inputRef?.(el)}
          flexGrow={1}
          focused={props.focused() && !asked()}
          placeholder={props.session()?.turn ? "queue a message… (ctrl+g sends it into the turn)" : "ask DELFIN…"}
          placeholderColor={t().faint}
          textColor={t().bright}
          backgroundColor="transparent"
          focusedBackgroundColor="transparent"
          cursorColor={t().accent}
          onSubmit={((value: string) => props.onSubmit(value)) as never}
        />
      </box>
      <text fg={t().faint}>{rule()}</text>
      <Hints {...props} />
    </box>
  )
}

function Approval(props: { asked: Asked; theme: () => Theme }): JSX.Element {
  const t = props.theme
  const a = props.asked
  const title = a.kind === "ask" ? ` ${a.header || "question"} ` : a.kind === "plan" ? " plan " : ` ${a.tool} `
  return (
    <box flexDirection="column" border borderStyle="rounded" borderColor={t().warn}
         title={title} titleAlignment="left" paddingLeft={1} paddingRight={1}>
      <Show when={a.kind === "ask"} fallback={
        <For each={previewLines(a)}>
          {(line) => <text fg={line.startsWith("+") ? t().ok : line.startsWith("-") ? t().err : t().text}>{line}</text>}
        </For>
      }>
        <text fg={t().bright}>{a.question ?? ""}</text>
        <For each={a.choices ?? []}>
          {(choice, i) => (
            <text>
              <span style={{ fg: t().accentSoft, attributes: TextAttributes.BOLD }}>{`  ${i() + 1}. `}</span>
              <span style={{ fg: t().text }}>{choice}</span>
              <span style={{ fg: t().faint }}>{a.choiceDetails?.[i()] ? `  ${a.choiceDetails[i()]}` : ""}</span>
            </text>
          )}
        </For>
      </Show>
      <text>
        <For each={a.kind === "ask" ? [{ key: "esc", label: "dismiss" }] : a.options}>
          {(o) => (
            <>
              <span style={{ fg: t().accentSoft, attributes: TextAttributes.BOLD }}>{`[${o.key}]`}</span>
              <span style={{ fg: t().text }}>{` ${o.label}   `}</span>
            </>
          )}
        </For>
      </text>
    </box>
  )
}

function Hints(props: FooterProps): JSX.Element {
  const t = props.theme
  const text = () => {
    const s = props.session()
    const parts: string[] = []
    if (props.connection() !== "up") parts.push(props.connection() === "reconnecting" ? "reconnecting to the agent service…" : props.connection())
    if (s?.approval) parts.push(s.approval)
    if (s?.usage.costUsd) parts.push(`$${s.usage.costUsd.toFixed(4)}`)
    parts.push(`${props.keymap.label("interrupt")} interrupt`, `${props.keymap.label("cycleMode")} mode`,
               `${props.keymap.label("agents")} agents`)
    return "  " + parts.join(" · ")
  }
  return <text fg={t().faint}>{text()}</text>
}
