/**
 * The open question, in the footer, and the keys that answer it.
 *
 *   confirm   the REPL's keys, exactly the ones the session offers:
 *             y yes · n no · A always this command · k always commands like
 *             this · R always read here · e stop asking about writes ·
 *             a abort · ? explain;   esc = n, as in the REPL;
 *             N = no, and say why (the reason goes to the model)
 *   plan      d approve → default · e approve → acceptEdits · n reject
 *             (N: reject and say why)
 *   ask       1-9 pick (several: toggle, enter sends) · esc dismisses
 *
 * A and k save a rule that is global and permanent: a second step shows
 * the rule and its reach, and only y saves it (n or esc approve once).
 *
 * Every piece of state belongs to one request (by its id): picks, the
 * second step and a half-typed reason never carry over to the next
 * question. Keys typed in a burst when a question opens (the end of a
 * message being typed) are not taken as answers: an answer is a key after
 * a pause.
 *
 * The panel is built as a list of lines that each fit the width, so the
 * footer knows exactly how many rows it takes (`approvalLines`).
 */
import { TextAttributes, type KeyEvent } from "@opentui/core"
import { For, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import type { Asked } from "../store/state"
import type { Theme } from "../theme/theme"
import { fit, width, wrap } from "../text"

/** Lines of a preview shown at most (a diff, a command). */
export const PREVIEW_MAX = 12
/** A key this soon after the question opened, or after the key before it,
 *  is typeahead, not an answer. */
export const ANSWER_PAUSE_MS = 350

type Tone = "text" | "bright" | "muted" | "faint" | "ok" | "err" | "warn" | "key"

export interface Span {
  text: string
  tone: Tone
  bold?: boolean
}
export type Line = Span[]

export interface ApprovalState {
  requestId: string
  mode: "options" | "help" | "persist" | "reason"
  choice: string           // the key the second step is for (A, k, n)
  picked: number[]         // a multi-select question's picks (1-based)
  reason: string
  openedAt: number
  lastKeyAt: number
  swallowed: boolean
  /** Answered from here; waiting for the session to settle it. */
  sent: boolean
}

function fresh(requestId: string, now: number): ApprovalState {
  return { requestId, mode: "options", choice: "", picked: [], reason: "", openedAt: now,
           lastKeyAt: 0, swallowed: false, sent: false }
}

/** Spans laid out in rows of `max` columns, broken between spans. */
export function flow(spans: Span[], max: number): Line[] {
  const rows: Line[] = []
  let row: Line = []
  let used = 0
  for (const span of spans) {
    const w = width(span.text)
    if (used > 0 && used + w > max) {
      rows.push(row)
      row = []
      used = 0
    }
    row.push(used === 0 ? { ...span, text: fit(span.text.trimStart(), max) } : span)
    used += w
  }
  if (row.length) rows.push(row)
  return rows
}

// -- what the panel shows --------------------------------------------------------

function previewLines(asked: Asked, max: number): Line[] {
  const lines = String(asked.preview ?? "").split("\n")
  const shown = lines.length <= PREVIEW_MAX ? lines
    : [...lines.slice(0, PREVIEW_MAX - 1), `… ${lines.length - PREVIEW_MAX + 1} more lines`]
  return shown.map((line) => [{
    text: fit(line, max),
    tone: line.startsWith("+") && !line.startsWith("+++") ? "ok"
      : line.startsWith("-") && !line.startsWith("---") ? "err" : "text",
  }])
}

/** The panel's lines for `asked` in state `ui`, inside `max` columns. */
export function approvalLines(asked: Asked, ui: ApprovalState | null, max: number): Line[] {
  const mine = ui?.requestId === asked.requestId ? ui : null
  const mode = mine?.mode ?? "options"
  const lines: Line[] = []
  if (mine?.sent) {
    // The answer is on its way: the panel says so until the session settles it.
    const body = asked.kind === "ask" ? wrap(asked.question || "(no question)", max).map((row) =>
      [{ text: row, tone: "bright" as Tone }]) : previewLines(asked, max)
    return [...body, [{ text: fit("answered — waiting for the session…", max), tone: "faint" }]]
  }
  if (asked.kind === "ask") {
    for (const row of wrap(asked.question || "(no question given)", max)) lines.push([{ text: row, tone: "bright" }])
    const choices = asked.choices ?? []
    choices.slice(0, 9).forEach((choice, i) => {
      const n = i + 1
      const on = ui?.requestId === asked.requestId && ui.picked.includes(n)
      const mark = asked.multiSelect ? (on ? "[x] " : "[ ] ") : ""
      const detail = asked.choiceDetails?.[i] ? `  ${asked.choiceDetails[i]}` : ""
      const head = `  ${n}. ${mark}`
      lines.push([{ text: head, tone: "key", bold: true },
                  { text: fit(choice + detail, Math.max(1, max - width(head))), tone: on ? "bright" : "text" }])
    })
    const keys: Span[] = asked.multiSelect
      ? [{ text: "[1-9] pick  ", tone: "key", bold: true }, { text: "[enter] send  ", tone: "key", bold: true },
         { text: "[esc] dismiss", tone: "key", bold: true }]
      : choices.length
        ? [{ text: "[1-9] answer  ", tone: "key", bold: true }, { text: "[esc] dismiss", tone: "key", bold: true }]
        : [{ text: "[enter] or [esc] dismiss", tone: "key", bold: true }]
    lines.push(...flow(keys, max))
    return lines
  }
  if (mode === "persist") {
    // k: the pattern it saves; A: the command itself, exactly.
    const command = String(asked.preview ?? "").split("\n")[0]!.replace(/^\$\s*/, "")
    const pattern = ui!.choice === "k" ? String(asked.suggestion ?? "") : `${command}  (exactly this command)`
    lines.push([{ text: "persist  ", tone: "warn", bold: true }, { text: fit(pattern, max - 9), tone: "bright" }])
    for (const row of wrap("Every future command matching it runs without asking — in this project and "
                           + "every other, in this session and every future one.", max)) {
      lines.push([{ text: row, tone: "muted" }])
    }
    lines.push(...flow([{ text: "[y] persist   ", tone: "key", bold: true },
                        { text: "[n] just this once", tone: "key", bold: true }], max))
    return lines
  }
  if (mode === "help") {
    for (const o of asked.options) {
      const head = `[${o.key}] `
      const rest = wrap(`${o.label}${o.detail ? ` — ${o.detail}` : ""}`, Math.max(1, max - width(head)))
      rest.forEach((row, i) => lines.push([{ text: i === 0 ? head : " ".repeat(width(head)), tone: "key", bold: true },
                                            { text: row, tone: "text" }]))
    }
    if (offers(asked, "n")) {
      lines.push([{ text: "[N] ", tone: "key", bold: true }, { text: "no, and say why: the reason goes to the model", tone: "text" }])
    }
    if (asked.flags?.outsideRead && !asked.readDir) {
      for (const row of wrap("no 'always' here: a home, system or key directory is approved one read at a time.", max)) {
        lines.push([{ text: row, tone: "muted" }])
      }
    }
    lines.push([{ text: "any key: back", tone: "faint" }])
    return lines
  }
  lines.push(...previewLines(asked, max))
  if (mode === "reason") {
    lines.push([{ text: asked.kind === "plan" ? "why not? " : "why? ", tone: "warn", bold: true },
                { text: fit(ui!.reason, max - 10) + "▏", tone: "bright" }])
    lines.push([{ text: fit("enter refuses with this reason · esc back", max), tone: "faint" }])
    return lines
  }
  const keys: Span[] = asked.options.map((o) => ({ text: `[${o.key}] ${o.label}   `, tone: "key" as Tone, bold: true }))
  if (offers(asked, "n")) keys.splice(keys.findIndex((k) => k.text.startsWith("[n]")) + 1, 0,
                                      { text: "[N] no + why   ", tone: "key", bold: true })
  lines.push(...flow(keys, max))
  return lines
}

function offers(asked: Asked, key: string): boolean {
  return asked.options.some((o) => o.key === key)
}

/** Rows the panel takes: its lines and its border. */
export function approvalRows(asked: Asked | undefined, ui: ApprovalState | null, termWidth: number): number {
  if (!asked) return 0
  return 2 + approvalLines(asked, ui, innerWidth(termWidth)).length
}

export function innerWidth(termWidth: number): number {
  return Math.max(10, termWidth - 4)             // the border and one column of padding each side
}

// -- the keys ----------------------------------------------------------------------

export interface ApprovalDeps {
  /** Send the answer; false (or a promise of it) when it did not go out. */
  reply(asked: Asked, params: Record<string, unknown>): Promise<boolean> | boolean | void
  say(text: string): void
  now?: () => number
  pauseMs?: number
}

export interface ApprovalUi {
  state: ApprovalState
  /** A key while `asked` is open; true when it was used (or swallowed). */
  handleKey(asked: Asked, key: KeyEvent): boolean
  /** Make sure the state is `asked`'s (the question changed). */
  sync(asked: Asked | undefined): void
}

function keyChar(key: KeyEvent): string {
  return key.sequence && key.sequence.length === 1 && key.sequence >= " " ? key.sequence : key.name
}

export function createApprovalUi(deps: ApprovalDeps): ApprovalUi {
  const now = deps.now ?? (() => Date.now())
  const pause = deps.pauseMs ?? ANSWER_PAUSE_MS
  const [state, setState] = createStore<ApprovalState>(fresh("", 0))

  const sync = (asked: Asked | undefined) => {
    if (asked && state.requestId !== asked.requestId) setState(fresh(asked.requestId, now()))
    if (!asked && state.requestId) setState(fresh("", 0))
  }

  const answer = (asked: Asked, params: Record<string, unknown>) => {
    setState({ mode: "options", choice: "", reason: "", sent: true })
    void Promise.resolve(deps.reply(asked, params)).then((ok) => {
      // Not delivered: the question is the user's to answer again.
      if (ok === false && state.requestId === asked.requestId) setState({ sent: false })
    })
  }

  const handleKey = (asked: Asked, key: KeyEvent): boolean => {
    sync(asked)
    if (state.sent) return true                    // answered: nothing more to say
    const t = now()
    const burst = t - state.openedAt < pause || t - state.lastKeyAt < pause
    setState({ lastKeyAt: t })
    // Typing a reason is typing: no pause needed between its keys.
    if (burst && state.mode !== "reason") {
      if (!state.swallowed) {
        setState({ swallowed: true })
        deps.say("a question just opened: keys typed with it were not used as answers")
      }
      return true
    }
    const char = keyChar(key)
    if (asked.kind === "ask") return askKey(asked, key, char)
    switch (state.mode) {
      case "help":
        setState({ mode: "options" })
        return true
      case "persist":
        if (char === "y") answer(asked, { choice: state.choice, persist: true })
        else if (char === "n" || key.name === "escape") answer(asked, { choice: state.choice, persist: false })
        return true
      case "reason":
        if (key.name === "escape") setState({ mode: "options", reason: "" })
        else if (key.name === "return") answer(asked, { choice: "n", reason: state.reason.trim() })
        else if (key.name === "backspace") setState({ reason: state.reason.slice(0, -1) })
        else if (key.name === "space") setState({ reason: state.reason + " " })
        else if (!key.ctrl && !key.meta && key.sequence && key.sequence.length >= 1 && !/[\x00-\x1f\x7f]/.test(key.sequence)) {
          setState({ reason: state.reason + key.sequence })
        }
        return true
    }
    if (key.name === "escape") {
      // Esc refuses, as in the REPL.
      if (offers(asked, "n")) answer(asked, { choice: "n" })
      return true
    }
    if (char === "N" && offers(asked, "n")) {
      setState({ mode: "reason", reason: "", choice: "n" })
      return true
    }
    if (char === "?" && offers(asked, "?")) {
      setState({ mode: "help" })
      return true
    }
    if (!offers(asked, char)) return true         // a key the question did not offer means nothing
    if ((char === "A" || char === "k") && asked.kind === "confirm") {
      setState({ mode: "persist", choice: char })
      return true
    }
    answer(asked, { choice: char })
    return true
  }

  const askKey = (asked: Asked, key: KeyEvent, char: string): boolean => {
    const choices = (asked.choices ?? []).slice(0, 9)
    if (key.name === "escape" || (key.name === "return" && choices.length === 0)) {
      answer(asked, { answers: [] })
      return true
    }
    const n = Number(char)
    if (Number.isInteger(n) && n >= 1 && n <= choices.length) {
      if (!asked.multiSelect) {
        answer(asked, { answers: [choices[n - 1]] })
        return true
      }
      const picked = state.picked.includes(n) ? state.picked.filter((x) => x !== n) : [...state.picked, n].sort()
      setState({ picked })
      return true
    }
    if (asked.multiSelect && key.name === "return") {
      answer(asked, { answers: state.picked.map((i) => choices[i - 1]) })
      return true
    }
    return true
  }

  return { state, handleKey, sync }
}

// -- drawing -------------------------------------------------------------------------

export function ApprovalPanel(props: { asked: Asked; ui: ApprovalUi; theme: () => Theme; width: () => number }): JSX.Element {
  const t = props.theme
  const tone = (s: Span) => {
    const theme = t()
    return s.tone === "key" ? theme.accentSoft : theme[s.tone]
  }
  const title = () => {
    const a = props.asked
    return a.kind === "ask" ? ` ${a.header || "question"} ` : a.kind === "plan" ? " plan " : ` ${a.tool} `
  }
  return (
    <box flexDirection="column" border borderStyle="rounded" borderColor={t().warn}
         title={title()} titleAlignment="left" paddingLeft={1} paddingRight={1}>
      <For each={approvalLines(props.asked, props.ui.state, innerWidth(props.width()))}>
        {(line) => (
          <text wrapMode="none">
            <For each={line}>
              {(span) => <span style={{ fg: tone(span), attributes: span.bold ? TextAttributes.BOLD : 0 }}>{span.text}</span>}
            </For>
          </text>
        )}
      </For>
    </box>
  )
}
