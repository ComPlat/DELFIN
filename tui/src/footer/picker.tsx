/**
 * A list to pick from, in the footer: the command palette (ctrl+p) and the
 * native pickers (model, permissions, effort, mode, resume).
 *
 *   type      narrows the list (every word must appear, in any order)
 *   ↑ ↓       select;  enter picks;  esc closes
 *
 * The list comes from the session (session.choices, session.commands,
 * sessions.past) and may take a moment: the picker opens at once and says
 * "loading…" until it is there.
 */
import type { KeyEvent } from "@opentui/core"
import { For, Show, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import { fit, width } from "../text"
import type { Theme } from "../theme/theme"

/** Rows of choices shown at once. */
export const PICKER_ROWS = 8

export interface Choice {
  label: string
  detail?: string
  /** What picking it does is the opener's: the choice is handed back. */
  value: string
  /** Marked as what is in force now (the current model, the posture). */
  current?: boolean
}

export interface PickerState {
  open: boolean
  title: string
  loading: boolean
  error: string
  choices: Choice[]
  query: string
  selected: number
}

export interface Picker {
  state: PickerState
  /** Open with a title, and the choices (or a promise of them); `more`
   *  are added when they arrive, the choices are pickable at once. */
  open(title: string, choices: Choice[] | Promise<Choice[]>, onPick: (choice: Choice) => void,
       more?: Promise<Choice[]>): void
  close(): void
  /** The choices that match what is typed, in their order. */
  shown(): Choice[]
  handleKey(key: KeyEvent): boolean
  /** Rows it takes in the footer. */
  height(): number
}

/** Every word of `query` in `text`, case-insensitively. */
export function matches(text: string, query: string): boolean {
  const hay = text.toLowerCase()
  return query.toLowerCase().split(/\s+/).filter(Boolean).every((word) => hay.includes(word))
}

export function createPicker(): Picker {
  const [state, setState] = createStore<PickerState>({ open: false, title: "", loading: false, error: "",
                                                       choices: [], query: "", selected: 0 })
  let pick: ((choice: Choice) => void) | null = null
  let generation = 0

  const shown = () => state.choices.filter((c) => matches(`${c.label} ${c.detail ?? ""}`, state.query))

  const open = (title: string, choices: Choice[] | Promise<Choice[]>, onPick: (choice: Choice) => void,
                more?: Promise<Choice[]>) => {
    const mine = ++generation
    pick = onPick
    const now = Array.isArray(choices) ? choices : []
    setState({ open: true, title, loading: !Array.isArray(choices), error: "", choices: now, query: "",
               selected: Math.max(0, now.findIndex((c) => c.current)) })
    if (!Array.isArray(choices)) {
      void choices.then((list) => {
        if (mine !== generation) return
        const current = Math.max(0, list.findIndex((c) => c.current))
        setState({ loading: false, choices: list, selected: current })
      }).catch((error) => {
        if (mine !== generation) return
        setState({ loading: false, error: String((error as Error)?.message ?? error) })
      })
    }
    void more?.then((list) => {
      if (mine === generation) setState({ choices: [...state.choices, ...list] })
    }).catch(() => {})
  }

  const close = () => {
    generation++
    pick = null
    setState({ open: false, choices: [], query: "", selected: 0, error: "" })
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    const list = shown()
    if (key.name === "escape" || (key.ctrl && (key.name === "c" || key.name === "p"))) {
      close()
      return true
    }
    if (key.name === "up" || key.name === "down") {
      if (list.length) {
        const n = list.length
        setState({ selected: (state.selected + (key.name === "up" ? n - 1 : 1)) % n })
      }
      return true
    }
    if (key.name === "return" || key.name === "tab") {
      const chosen = list[Math.min(state.selected, list.length - 1)]
      if (!chosen) return true
      const then = pick
      close()
      then?.(chosen)
      return true
    }
    let query = state.query
    if (key.name === "backspace") query = query.slice(0, -1)
    else if (key.name === "space") query += " "
    else if (!key.ctrl && !key.meta && key.sequence && key.sequence.length >= 1 && !/[\x00-\x1f\x7f]/.test(key.sequence)) {
      query += key.sequence
    } else return true                     // nothing else leaves the picker
    setState({ query, selected: 0 })
    return true
  }

  const height = () => {
    if (!state.open) return 0
    const n = state.loading || state.error ? 1 : Math.max(1, Math.min(PICKER_ROWS, shown().length))
    return n + 3                                     // the border, and the filter line
  }

  return { state, open, close, shown, handleKey, height }
}

export function PickerView(props: { picker: Picker; theme: () => Theme; width: () => number }): JSX.Element {
  const t = props.theme
  const p = props.picker
  const inner = () => Math.max(10, props.width() - 4)
  const rows = () => {
    const list = p.shown()
    const sel = Math.min(p.state.selected, Math.max(0, list.length - 1))
    const first = Math.max(0, Math.min(sel - PICKER_ROWS + 1, list.length - PICKER_ROWS))
    return list.slice(first, first + PICKER_ROWS).map((c, i) => ({ c, on: first + i === sel }))
  }
  const line = (c: Choice, on: boolean) => {
    const mark = on ? "› " : "  "
    const label = fit(c.label, Math.max(4, Math.floor(inner() * 0.5)))
    const detail = c.detail ? fit(c.detail, Math.max(0, inner() - width(mark) - width(label) - 2 - (c.current ? 9 : 0))) : ""
    return { mark, label, detail, now: c.current ? "  (now)" : "" }
  }
  return (
    <box flexShrink={0} flexDirection="column" border borderStyle="rounded" borderColor={t().accent}
         title={` ${p.state.title} `} titleAlignment="left" paddingLeft={1} paddingRight={1}>
      <text wrapMode="none">
        <span style={{ fg: t().accentSoft }}>{"filter: "}</span>
        <span style={{ fg: t().bright }}>{fit(p.state.query, inner() - 30) + "▏"}</span>
        <span style={{ fg: t().faint }}>{"  ↑↓ enter · esc"}</span>
      </text>
      <Show when={!p.state.loading && !p.state.error} fallback={
        <text fg={p.state.error ? t().err : t().faint}>{p.state.error ? fit(p.state.error, inner()) : "loading…"}</text>
      }>
        <Show when={p.shown().length} fallback={<text fg={t().faint}>{"nothing matches"}</text>}>
          <For each={rows()}>
            {(row) => {
              const l = line(row.c, row.on)
              return (
                <text wrapMode="none">
                  <span style={{ fg: row.on ? t().accent : t().faint }}>{l.mark}</span>
                  <span style={{ fg: row.on ? t().bright : t().text, bold: Boolean(row.on) }}>{l.label}</span>
                  <span style={{ fg: t().ok }}>{l.now}</span>
                  <span style={{ fg: t().faint }}>{l.detail ? `  ${l.detail}` : ""}</span>
                </text>
              )
            }}
          </For>
        </Show>
      </Show>
    </box>
  )
}
