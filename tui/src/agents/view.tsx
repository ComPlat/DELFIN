/**
 * The agent view: ← on an empty prompt. Every session the agent service
 * runs, grouped the way Claude Code groups them, plus the other open
 * DELFIN sessions (dashboards, `delfin-agent chat`) as rows to watch and
 * message:
 *
 *   PINNED · NEEDS INPUT · WORKING · DONE · ELSEWHERE
 *
 *   ↑ ↓          select              → / enter   open it in this window
 *   space        peek at it (and reply to it, without opening it)
 *   ctrl+x       stop it (twice); an ended one leaves the list
 *   ctrl+r       rename it           ctrl+t      pin it
 *   alt+1…9      open the n-th       esc / ←     back
 *   typing + enter                   a new session, in this window's folder
 *
 * A row says what the session is doing in a few words, from its summary
 * (no model writes them): the question it waits on, the tool it runs, the
 * last line of its answer. Rows keep their place within a group (newest
 * session first) and the selection follows its session, so the list does
 * not move under the cursor while sessions work.
 */
import { TextAttributes, type KeyEvent } from "@opentui/core"
import { For, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import type { Keymap } from "../keymap/keymap"
import { askedWhat } from "../scrollback/present"
import { toolHeadline } from "../scrollback/tools"
import type { AppState, SessionState } from "../store/state"
import { fit, width } from "../text"
import { SPINNER, type Theme } from "../theme/theme"

export type Group = "pinned" | "needs" | "working" | "done" | "elsewhere"

export const GROUP_TITLES: Record<Group, string> = {
  pinned: "PINNED", needs: "NEEDS INPUT", working: "WORKING", done: "DONE", elsewhere: "ELSEWHERE · watch and message",
}

const ORDER: Group[] = ["pinned", "needs", "working", "done", "elsewhere"]

export interface AgentRow {
  kind: "session" | "elsewhere"
  /** The session id, or an elsewhere session's presence key. */
  id: string
  group: Group
  /** What the session needs or does, whatever group it is shown in.
   *  "sleeping": no worker (it slept, or a service before stopped); it
   *  wakes when opened or sent a message. */
  status: Exclude<Group, "pinned" | "elsewhere"> | "ended" | "sleeping" | "elsewhere"
  name: string
  activity: string
  age: string
  cost: string
  current: boolean
  workspace: string
}

export function ageOf(seconds: number, now: number): string {
  const d = Math.max(0, now / 1000 - seconds)
  if (d < 60) return `${Math.floor(d)}s`
  if (d < 3600) return `${Math.floor(d / 60)}m`
  if (d < 86400) return `${Math.floor(d / 3600)}h`
  return `${Math.floor(d / 86400)}d`
}

function statusOf(s: SessionState): AgentRow["status"] {
  if (s.state === "ended") return "ended"
  if (s.state === "sleeping") return "sleeping"
  const pending = s.pending.length || Number(s.summary?.pending ?? 0)
  if (pending > 0 || s.state === "needs_input") return "needs"
  if (s.state === "running" || s.state === "starting" || s.turn) return "working"
  return "done"
}

/** Whether a session waits for the user: a question is open. */
export function needsInput(s: SessionState): boolean {
  return statusOf(s) === "needs"
}

function activityOf(s: SessionState, status: AgentRow["status"]): string {
  const summary = (s.summary ?? {}) as Record<string, any>
  if (status === "ended") return `ended${s.ended?.reason ? ` (${s.ended.reason})` : ""}`
  if (status === "sleeping") {
    const last = String(summary.answer || summary.lastPrompt || "").split("\n")[0]
    return `asleep${last ? ` · ${last}` : ""}`
  }
  if (status === "needs") return `waiting for you: ${s.pending[0] ? askedWhat(s.pending[0]) : summary.waiting || "a question"}`
  if (status === "working") return String(summary.activity || s.turn?.label || "working")
  return String(summary.answer || summary.lastPrompt || "idle").split("\n")[0]!
}

/** The rows, in their groups; within a group the newest session first. */
export function agentRows(state: AppState, pins: ReadonlySet<string>, now: number = Date.now()): AgentRow[] {
  const rows: (AgentRow & { created: number; at: number })[] = []
  state.order.forEach((id, at) => {
    const s = state.sessions[id]
    if (!s) return
    const status = statusOf(s)
    const summary = (s.summary ?? {}) as Record<string, any>
    const cost = Math.max(Number(summary.costUsd ?? 0), s.usage.costUsd)
    rows.push({
      kind: "session", id, status,
      group: pins.has(id) ? "pinned" : status === "ended" || status === "sleeping" ? "done" : status,
      // An untitled session is named by what it was first asked: a name
      // that does not change with every prompt.
      name: s.title || String(summary.title || "") || String(summary.firstPrompt || summary.lastPrompt || "").split("\n")[0]
        || id.slice(0, 8),
      activity: activityOf(s, status),
      age: summary.lastActivity ? ageOf(Number(summary.lastActivity), now) : "",
      cost: cost > 0 ? `$${cost.toFixed(2)}` : "",
      current: state.current === id,
      workspace: s.workspace || String(summary.workspace || ""),
      created: Number(summary.createdAt ?? 0), at,
    })
  })
  state.elsewhere.forEach((e, at) => {
    rows.push({
      kind: "elsewhere", id: e.key, group: "elsewhere", status: "elsewhere",
      name: e.title || e.key, activity: [e.workspace, e.host].filter(Boolean).join(" · "),
      age: "", cost: "", current: false, workspace: e.workspace, created: 0, at,
    })
  })
  rows.sort((a, b) => ORDER.indexOf(a.group) - ORDER.indexOf(b.group)
    || b.created - a.created || b.at - a.at)
  return rows.map(({ created: _c, at: _a, ...row }) => row)
}

export interface Peek {
  id: string
  lines: string[]
  loading: boolean
}

export interface AgentViewState {
  open: boolean
  /** The selected row's id: the selection follows its session. */
  selectedId: string | null
  /** Where the selection was, for when its row leaves the list. */
  selectedAt: number
  peek: Peek | null
  /** ctrl+x pressed once on this row: a second one within 3 s stops it. */
  armed: { id: string; at: number } | null
  /** ctrl+r on this row: what is typed next is its new name. */
  renaming: string | null
  /** Rows taken out of the list (ended sessions removed with ctrl+x). */
  hidden: string[]
}

export interface AgentViewDeps {
  state: () => AppState
  pins: () => ReadonlySet<string>
  keymap: Keymap
  togglePin(id: string): void
  open(row: AgentRow): void
  peek(id: string): Promise<string[]>
  stop(id: string): void
  rename(id: string, title: string): void
  message(row: AgentRow, text: string): void
  newSession(text: string): void
  draft(): string
  setDraft(text: string): void
  say(text: string): void
}

export interface AgentView {
  state: AgentViewState
  rows(now?: number): AgentRow[]
  /** The selected row's index in rows(). */
  selected(list?: AgentRow[]): number
  open(): void
  close(): void
  /** A key while the view is open; true when it was the view's. */
  handleKey(key: KeyEvent): boolean
  /** Enter with text typed: a new session, or a reply to the peeked one. */
  submit(text: string): boolean
  /** Fetch the peek again (its session moved on). */
  refreshPeek(): void
}

export function createAgentView(deps: AgentViewDeps): AgentView {
  const [state, setState] = createStore<AgentViewState>({
    open: false, selectedId: null, selectedAt: 0, peek: null, armed: null, renaming: null, hidden: [],
  })
  const km = deps.keymap
  const rows = (now?: number) => agentRows(deps.state(), deps.pins(), now).filter((r) => !state.hidden.includes(r.id))
  const selected = (list = rows()) => {
    const at = list.findIndex((r) => r.id === state.selectedId)
    return at >= 0 ? at : Math.min(state.selectedAt, Math.max(0, list.length - 1))
  }
  const selectedRow = (list = rows()) => list[selected(list)]
  const select = (list: AgentRow[], at: number) => {
    const row = list[at]
    setState({ selectedId: row?.id ?? null, selectedAt: at, armed: null })
  }

  const open = () => {
    const list = rows()
    const current = Math.max(0, list.findIndex((r) => r.current))
    setState({ open: true, peek: null, armed: null, renaming: null })
    select(list, current)
  }
  const close = () => {
    if (state.renaming) deps.setDraft("")
    setState({ open: false, peek: null, armed: null, renaming: null })
  }

  let asked = 0
  const peekAt = (row: AgentRow, quiet = false) => {
    if (row.kind === "elsewhere") {
      setState({ peek: { id: row.id, loading: false, lines: [row.name, row.activity,
        "a session elsewhere: type and press enter to message it"] } })
      return
    }
    if (row.status === "ended") {
      setState({ peek: { id: row.id, loading: false, lines: [row.activity, "→ resumes it in this window"] } })
      return
    }
    const mine = ++asked
    if (!quiet) setState({ peek: { id: row.id, loading: true, lines: [] } })
    const done = (lines: string[]) => {
      if (mine === asked && state.peek?.id === row.id) setState({ peek: { id: row.id, loading: false, lines } })
    }
    void deps.peek(row.id).then(done, (error) => done([String(error?.message ?? error)]))
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    const list = rows()
    const row = selectedRow(list)
    // Keys that edit text are the composer's while something is typed.
    const typing = deps.draft().length > 0
    if (key.name === "escape" || (km.matches("agentsBack", key) && !typing)) {
      if (state.renaming) {
        setState({ renaming: null })
        deps.setDraft("")
      } else if (state.peek) setState({ peek: null })
      else close()
      return true
    }
    if (km.matches("agentsUp", key) || km.matches("agentsDown", key)) {
      if (list.length) {
        const n = list.length
        select(list, (selected(list) + (km.matches("agentsUp", key) ? n - 1 : 1)) % n)
        if (state.peek) peekAt(selectedRow()!)
      }
      return true
    }
    if ((key.meta || key.option) && /^[1-9]$/.test(key.name)) {
      const target = list.filter((r) => r.kind === "session")[Number(key.name) - 1]
      if (target) {
        close()
        if (!target.current) deps.open(target)
      }
      return true
    }
    if (!row) return false
    if (state.renaming && km.matches("submit", key)) {
      const name = deps.draft().trim()
      if (name) deps.rename(state.renaming, name)
      setState({ renaming: null })
      deps.setDraft("")
      return true
    }
    if (km.matches("agentsOpen", key) && !typing) {
      if (row.kind === "elsewhere") {
        deps.say("a session elsewhere is watched, not opened: type and press enter to message it")
        return true
      }
      close()
      if (!row.current) deps.open(row)
      return true
    }
    if (km.matches("agentsPeek", key) && !typing) {
      if (state.peek?.id === row.id) setState({ peek: null })
      else peekAt(row)
      return true
    }
    if (km.matches("agentsStop", key)) {
      if (row.kind === "elsewhere") {
        deps.say("a session elsewhere is stopped where it runs")
        return true
      }
      if (row.status === "ended") {
        setState({ hidden: [...state.hidden, row.id], armed: null, peek: null })
        return true
      }
      const armed = state.armed
      if (armed && armed.id === row.id && Date.now() - armed.at < 3000) {
        setState({ armed: null })
        deps.stop(row.id)
        deps.say(`stopping ${row.name}: it ends and is saved`)
      } else {
        setState({ armed: { id: row.id, at: Date.now() } })
        deps.say(`${km.label("agentsStop")} again to stop ${row.name} (it ends and is saved)`)
      }
      return true
    }
    if (km.matches("agentsRename", key) && row.kind === "session") {
      setState({ renaming: row.id, peek: null })
      deps.setDraft(row.name)
      return true
    }
    if (km.matches("agentsPin", key) && row.kind === "session") {
      deps.togglePin(row.id)
      return true
    }
    return false                        // everything else is typing
  }

  const submit = (text: string): boolean => {
    if (!state.open) return false
    const body = text.trim()
    if (!body) return true
    if (state.renaming) {
      deps.rename(state.renaming, body)
      setState({ renaming: null })
      return true
    }
    const row = selectedRow()
    if (state.peek && row && state.peek.id === row.id && row.status !== "ended") {
      deps.message(row, body)
      return true
    }
    close()
    deps.newSession(body)
    return true
  }

  const refreshPeek = () => {
    const peek = state.peek
    if (!state.open || !peek || peek.loading) return
    const row = rows().find((r) => r.id === peek.id)
    if (row) peekAt(row, true)
  }

  return { state, rows, selected, open, close, handleKey, submit, refreshPeek }
}

/** The hints line while the view is open: what the keys do now, as many
 *  as fit in `room` columns -- the most needed ones first. */
export function agentHelp(view: AgentView, km: Keymap, room = Infinity): string {
  if (view.state.renaming) return "type the new name · enter renames · esc back"
  if (view.state.peek) return "type to reply to it · enter sends · space closes · esc back"
  const parts = [
    { at: 0, text: "↑↓ select" },
    { at: 1, text: `${km.label("agentsOpen")} open` },
    { at: 2, text: `${km.label("agentsPeek")} peek` },
    { at: 3, text: `${km.label("agentsStop")} stop` },
    { at: 4, text: `${km.label("agentsRename")} rename` },
    { at: 5, text: `${km.label("agentsPin")} pin` },
    { at: 6, text: "esc back" },
    { at: 7, text: "type to start a new session" },
  ]
  const need = [7, 1, 2, 6, 3, 0, 4, 5]
  const kept: typeof parts = []
  for (const at of need) {
    const next = [...kept, parts[at]!].sort((a, b) => a.at - b.at)
    if (width(next.map((p) => p.text).join(" · ")) > room) continue
    kept.splice(0, kept.length, ...next)
  }
  return kept.map((p) => p.text).join(" · ")
}

/** How many lines of a peek are shown under its row. */
export const PEEK_LINES = 12

export function AgentViewPanel(props: { view: AgentView; theme: () => Theme; width: () => number;
                                        height: () => number; tick: () => number; now: () => number;
                                        keymap?: Keymap }): JSX.Element {
  const t = props.theme
  const v = props.view
  const inner = () => Math.max(20, props.width() - 2)
  /** The list as lines: group titles, rows, and the peek under its row. */
  const lines = () => {
    const out: { kind: "title" | "row" | "peek" | "empty"; text?: string; row?: AgentRow; index?: number }[] = []
    const list = v.rows(props.now())
    const at = v.selected(list)
    let group: Group | null = null
    list.forEach((row, index) => {
      if (row.group !== group) {
        group = row.group
        out.push({ kind: "title", text: GROUP_TITLES[row.group] })
      }
      out.push({ kind: "row", row, index })
      if (v.state.peek?.id === row.id) {
        const peek = v.state.peek
        const room = Math.max(3, Math.min(PEEK_LINES, Math.floor(props.height() / 2)))
        for (const line of peek.loading ? ["loading…"] : peek.lines.slice(-room)) out.push({ kind: "peek", text: line })
      }
    })
    if (!list.length) out.push({ kind: "empty", text: "no sessions yet: type below to start one" })
    // Keep the selection (and its peek) in view: a window of the list around it.
    const room = Math.max(1, props.height() - 1)
    const sel = out.findIndex((l) => l.kind === "row" && l.index === at)
    const peekEnd = sel >= 0 && v.state.peek ? sel + 1 + out.slice(sel + 1).findIndex((l) => l.kind !== "peek") : sel
    const want = Math.max(sel, peekEnd < sel ? out.length - 1 : peekEnd)
    let first = Math.max(0, Math.min(sel - Math.floor(room / 3), out.length - room))
    if (want - first >= room) first = want - room + 1
    return out.slice(Math.max(0, first), Math.max(0, first) + room).map((line) => ({ ...line, selected: line.index === at }))
  }
  const glyph = (row: AgentRow) => row.status === "elsewhere" ? "◇"
    : row.status === "ended" ? "◼"
    : row.status === "sleeping" ? "☾"
    : row.status === "needs" ? "◉"
    : row.status === "working" ? SPINNER[props.tick() % SPINNER.length]!
    : "✓"
  const tone = (row: AgentRow) => row.status === "elsewhere" ? t().muted
    : row.status === "ended" || row.status === "sleeping" ? t().faint
    : row.status === "needs" ? t().warn
    : row.status === "done" ? t().ok : t().accent
  const summary = () => {
    const list = v.rows().filter((r) => r.kind === "session")
    const needs = list.filter((r) => r.status === "needs").length
    const working = list.filter((r) => r.status === "working").length
    return [`${list.length} session${list.length === 1 ? "" : "s"}`,
            needs ? `${needs} need${needs === 1 ? "s" : ""} input` : "",
            working ? `${working} working` : ""].filter(Boolean).join(" · ")
  }
  return (
    <box flexDirection="column" height={props.height()} flexShrink={0}>
      <text wrapMode="none">
        <span style={{ fg: t().accent, attributes: TextAttributes.BOLD }}>{"◆ AGENTS  "}</span>
        <span style={{ fg: t().muted }}>{summary()}</span>
      </text>
      <For each={lines()}>
        {(line) => {
          if (line.kind === "title") return <text fg={t().faint} wrapMode="none">{fit(` ${line.text}`, inner())}</text>
          if (line.kind === "empty") return <text fg={t().faint} wrapMode="none">{fit(`   ${line.text}`, inner())}</text>
          if (line.kind === "peek") return <text fg={t().muted} wrapMode="none">{`      │ ${fit(line.text ?? "", inner() - 8)}`}</text>
          const row = line.row!
          const on = line.selected
          // The marks change without the row (it is drawn again only when
          // the list does): read where they are drawn.
          const renaming = () => v.state.renaming === row.id
          const armed = () => v.state.armed?.id === row.id
          const right = [row.age, row.cost].filter(Boolean).join("  ")
          const name = () => fit(row.name, Math.max(8, Math.floor(inner() * 0.32)))
          const activity = () => renaming() ? "renaming: type the new name below"
            : armed() ? `${props.keymap?.label("agentsStop") ?? "ctrl+x"} again stops it` : row.activity
          const room = () => Math.max(0, inner() - 5 - width(name()) - 2 - (right ? width(right) + 2 : 0))
          return (
            <text wrapMode="none">
              <span style={{ fg: on ? t().accent : t().faint }}>{on ? " › " : "   "}</span>
              <span style={{ fg: tone(row) }}>{glyph(row) + " "}</span>
              <span style={{ fg: on ? t().bright : t().text, attributes: on || row.current ? TextAttributes.BOLD : 0 }}>{name()}</span>
              <span style={{ fg: armed() ? t().warn : renaming() ? t().accent : t().muted }}>{"  " + fit(activity(), room()).padEnd(room())}</span>
              <span style={{ fg: t().faint }}>{right ? "  " + right : ""}</span>
            </text>
          )
        }}
      </For>
    </box>
  )
}

/** A peek at a session: its last turn as a few lines (session.replay). */
export function peekLines(events: Record<string, any>[]): string[] {
  const out: string[] = []
  const started = new Map<string, Record<string, any>>()
  const open = new Map<string, Record<string, any>>()
  for (const e of events) {
    switch (e.type) {
      case "user.message":
        out.push(`❯ ${String(e.text ?? "").split("\n")[0]}`)
        break
      case "tool.started":
        started.set(String(e.callId), e)
        break
      case "tool.ended": {
        const s = started.get(String(e.callId)) ?? {}
        const head = toolHeadline(String(e.name ?? s.name ?? "tool"), s.input ?? {})
        out.push(`⏺ ${head.name}${head.arg ? `  ${head.arg}` : ""}`)
        break
      }
      case "text.ended":
        for (const line of String(e.text ?? "").split("\n").map((l) => l.trim()).filter(Boolean).slice(0, 6)) out.push(line)
        break
      case "notice":
        for (const line of String(e.text ?? "").split("\n").map((l) => l.trim()).filter(Boolean).slice(0, 3)) out.push(line)
        break
      case "permission.asked":
        open.set(String(e.requestId), e)
        break
      case "permission.resolved":
        open.delete(String(e.requestId))
        break
      case "turn.ended":
        out.push(e.stopped ? "◼ stopped" : e.error ? `✗ ${String(e.error).split("\n")[0]}` : `✓ ${(Number(e.durationMs ?? 0) / 1000).toFixed(1)}s`)
        break
    }
  }
  for (const asked of open.values()) out.push(`◉ waiting for you: ${askedWhat(asked as never)}`)
  return out.length ? out : ["(nothing yet)"]
}
