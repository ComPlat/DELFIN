/**
 * The agent view: ← on an empty prompt. Every session the agent service
 * runs, grouped by state, what needs the user from
 * every other DELFIN process (attention.ts: parked questions, failed
 * jobs, alerts, notices), and the other open DELFIN sessions (dashboards,
 * `delfin-agent chat`) as rows to watch and message:
 *
 *   PINNED · NEEDS YOU · WORKING · DONE · NOTICES · ELSEWHERE
 *
 *   ↑ ↓          select              → / enter   open it in this window
 *   space        peek at it (and reply to it, without opening it)
 *   ctrl+x       stop it (twice); an ended one leaves the list; an item
 *                is dismissed (what waits or failed: twice)
 *   ctrl+r       rename it           ctrl+t      pin it
 *   alt+1…9      open the n-th       esc / ←     back
 *   typing + enter                   a new session, in this window's folder
 *                                    (with an item's peek open: its answer)
 *
 * A row says what the session is doing in a few words, from its summary
 * (no model writes them): the question it waits on, the tool it runs, the
 * last line of its answer. Rows keep their place within a group (newest
 * session first) and the selection follows its session, so the list does
 * not move under the cursor while sessions work.
 */
import type { KeyEvent } from "@opentui/core"
import { For, createEffect, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import type { Keymap } from "../keymap/keymap"
import type { AttentionItem } from "../protocol/generated"
import { askedWhat } from "../scrollback/present"
import { toolHeadline } from "../scrollback/tools"
import type { AppState, SessionState } from "../store/state"
import { fit, middle, width, wrap } from "../text"
import { SPINNER, type Theme } from "../theme/theme"
import { MARK_GLYPH, ageOf, itemPeek, markOf, shortPath, whereOf } from "./attention"

export { ageOf } from "./attention"

export type Group = "pinned" | "needs" | "working" | "done" | "notices" | "elsewhere"

export const GROUP_TITLES: Record<Group, string> = {
  pinned: "PINNED", needs: "NEEDS YOU", working: "WORKING", done: "DONE", notices: "NOTICES",
  elsewhere: "ELSEWHERE · watch and message",
}

const ORDER: Group[] = ["pinned", "needs", "working", "done", "notices", "elsewhere"]

export interface AgentRow {
  /** attention: an item of what needs the user (attention.ts), from any
   *  DELFIN process; never one of this service's sessions. */
  kind: "session" | "elsewhere" | "attention"
  /** The session id, an elsewhere session's presence key, or
   *  "attention:" and the item's id. */
  id: string
  group: Group
  /** What the session needs or does, whatever group it is shown in.
   *  "sleeping": no worker (it slept, or a service before stopped); it
   *  wakes when opened or sent a message. "item": an attention row. */
  status: "needs" | "working" | "done" | "ended" | "sleeping" | "elsewhere" | "item"
  name: string
  activity: string
  age: string
  cost: string
  current: boolean
  workspace: string
  /** Where an item is (open on another machine, its folder), at the right. */
  where?: string
  item?: AttentionItem
}

/** What a session of this service needs or does. */
type SessionStatus = "needs" | "working" | "done" | "ended" | "sleeping"

function statusOf(s: SessionState): SessionStatus {
  if (s.state === "ended") return "ended"
  if (s.state === "sleeping") return "sleeping"
  const pending = s.pending.length || Number(s.summary?.pending ?? 0)
  if (pending > 0 || s.state === "needs_input") return "needs"
  if (s.state === "running" || s.state === "starting" || s.turn) return "working"
  return "done"
}

/** What a session is called: its name -- an untitled one by what it was
 *  first asked (a name that does not change with every prompt). */
export function sessionName(s: SessionState): string {
  const summary = (s.summary ?? {}) as Record<string, any>
  return s.title || String(summary.title || "") || String(summary.firstPrompt || summary.lastPrompt || "").split("\n")[0]
    || s.id.slice(0, 8)
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

/** Within a group: this service's sessions, then items, then sessions elsewhere. */
const RANK: Record<AgentRow["kind"], number> = { session: 0, attention: 1, elsewhere: 2 }

/** The rows, in their groups; within a group the newest session first,
 *  then the items in the service's order (what blocks, what failed, the
 *  newest first). */
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
      name: sessionName(s),
      activity: activityOf(s, status),
      age: summary.lastActivity ? ageOf(Number(summary.lastActivity), now) : "",
      cost: cost > 0 ? `$${cost.toFixed(2)}` : "",
      current: state.current === id,
      workspace: s.workspace || String(summary.workspace || ""),
      created: Number(summary.createdAt ?? 0), at,
    })
  })
  ;(state.attention ?? []).forEach((item, at) => {
    rows.push({
      kind: "attention", id: `attention:${item.id}`, group: item.group === "needs" ? "needs" : "notices",
      status: "item", name: item.name || item.kind, activity: item.text,
      age: item.createdAt ? ageOf(item.createdAt, now) : "", cost: "", current: false,
      workspace: item.workspace, where: whereOf(item, state.host), item, created: 0, at,
    })
  })
  state.elsewhere.forEach((e, at) => {
    rows.push({
      kind: "elsewhere", id: e.key, group: "elsewhere", status: "elsewhere",
      name: e.title || e.key, activity: [e.workspace, e.host].filter(Boolean).join(" · "),
      age: "", cost: "", current: false, workspace: e.workspace, created: 0, at,
    })
  })
  rows.sort((a, b) => ORDER.indexOf(a.group) - ORDER.indexOf(b.group) || RANK[a.kind] - RANK[b.kind]
    || (a.kind === "attention" ? a.at - b.at : b.created - a.created || b.at - a.at))
  return rows.map(({ created: _c, at: _a, ...row }) => row)
}

/** Whether an item takes an answer typed in the view: someone waits on
 *  it (a parked question, a confirmation, a plan). */
export function takesAnswer(item: AttentionItem | undefined): boolean {
  return Boolean(item?.blocking)
}

export interface Peek {
  id: string
  lines: string[]
  loading: boolean
  /** An item's peek: its lines wrapped to the width, from the top (a
   *  session's shows its last lines, cut). */
  wrap?: boolean
}

export interface AgentViewState {
  open: boolean
  /** The selected row's id: the selection follows its session. */
  selectedId: string | null
  /** Where the selection was, for when its row leaves the list. */
  selectedAt: number
  peek: Peek | null
  /** ctrl+x pressed once on this row: a second one within 3 s stops it
   *  (an item: dismisses it). */
  armed: { id: string; at: number } | null
  /** ctrl+r on this row: what is typed next is its new name. */
  renaming: string | null
  /** Rows taken out of the list (ended sessions removed with ctrl+x;
   *  items dismissed, until the service's list says so too). */
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
  /** → on an item: its session, a copy, the diagnosis, or a new session
   *  about it (attention.ts). */
  openItem(item: AttentionItem): void
  /** ctrl+x on an item: cleared, never answered; false when it was not. */
  dismiss(item: AttentionItem): Promise<boolean>
  /** Typing + enter on a peeked item that waits: its answer, for the next
   *  session in its folder; false when it was not taken. */
  answer(item: AttentionItem, text: string): Promise<boolean>
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
  /** ctrl+c: end a rename (the draft from before it comes back); false
   *  when there was none. */
  cancelRename(): boolean
}

export function createAgentView(deps: AgentViewDeps): AgentView {
  const [state, setState] = createStore<AgentViewState>({
    open: false, selectedId: null, selectedAt: 0, peek: null, armed: null, renaming: null, hidden: [],
  })
  const km = deps.keymap
  // A row taken out stays out while its session has ended; one that runs
  // again (resumed, woken) is listed again. A dismissed item stays out.
  const rows = (now?: number) => agentRows(deps.state(), deps.pins(), now)
    .filter((r) => !(state.hidden.includes(r.id) && (r.status === "ended" || r.kind === "attention")))
  /** What was typed when a rename began: it comes back when the rename ends. */
  let stash = ""
  const endRename = () => {
    setState({ renaming: null })
    deps.setDraft(stash)
    stash = ""
  }
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
    if (state.renaming) endRename()
    setState({ open: false, peek: null, armed: null, renaming: null })
  }

  let asked = 0
  const peekAt = (row: AgentRow, quiet = false) => {
    if (row.kind === "attention") {
      setState({ peek: { id: row.id, loading: false, wrap: true,
                         lines: itemPeek(row.item!, deps.state().host) } })
      return
    }
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
      if (state.renaming) endRename()
      else if (state.peek) setState({ peek: null })
      else close()
      return true
    }
    if (km.matches("agentsUp", key) || km.matches("agentsDown", key)) {
      justDismissed = false
      const peeked = state.peek ? list.some((r) => r.id === state.peek!.id) : false
      // A reply being typed stays with the session it is for -- while that
      // is in the list (an item answered elsewhere leaves it).
      if (typing && peeked) return true
      if (state.peek && !peeked) setState({ peek: null })
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
    if (km.matches("agentsRename", key) || km.matches("agentsPin", key)) {
      // The view's keys on every row (never the prompt's history search).
      if (!row || row.kind !== "session") {
        deps.say(!row ? "no session to rename or pin"
          : row.kind === "attention" ? `this is no session: ${km.label("agentsOpen")} opens it, ${km.label("agentsStop")} dismisses it`
          : "a session elsewhere is renamed and pinned where it runs")
        return true
      }
    }
    if (!row) return false
    if (state.renaming && km.matches("submit", key)) {
      const id = state.renaming
      const name = deps.draft().trim()
      if (name) deps.rename(id, name)
      endRename()
      return true
    }
    if (km.matches("agentsOpen", key) && !typing) {
      if (row.kind === "attention") {
        close()
        deps.openItem(row.item!)
        return true
      }
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
      if (row.kind === "attention") {
        dismissKey(row)
        return true
      }
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
    if (km.matches("agentsRename", key)) {
      if (!state.renaming) stash = deps.draft()
      setState({ renaming: row.id, peek: null })
      deps.setDraft(row.name)
      return true
    }
    if (km.matches("agentsPin", key)) {
      deps.togglePin(row.id)
      return true
    }
    return false                        // everything else is typing
  }

  /** ctrl+x on an item: a notice goes at once; what waits or failed is
   *  asked twice -- whoever waits on it gets no answer. Hidden at once,
   *  back in the list when the dismiss did not go through. */
  const dismissKey = (row: AgentRow) => {
    const item = row.item!
    const go = () => {
      // The row below moves up under the selection: a ctrl+x pressed twice
      // by habit must not dismiss it too.
      justDismissed = true
      setState({ hidden: [...state.hidden, row.id], armed: null, peek: state.peek?.id === row.id ? null : state.peek })
      void deps.dismiss(item).then((ok) => {
        if (ok) deps.say(`dismissed: ${row.name}`)
        else setState({ hidden: state.hidden.filter((h) => h !== row.id) })
      })
    }
    if (markOf(item) === "notice" && !justDismissed) {
      go()
      return
    }
    const armed = state.armed
    if (armed && armed.id === row.id && Date.now() - armed.at < 3000) {
      go()
      return
    }
    justDismissed = false
    setState({ armed: { id: row.id, at: Date.now() } })
    deps.say(item.blocking ? `${km.label("agentsStop")} again dismisses it: nobody waiting on it gets an answer`
      : `${km.label("agentsStop")} again dismisses it`)
  }
  /** A dismissal moved the next row under the selection (cleared when the
   *  selection is moved). */
  let justDismissed = false

  /** Typing + enter with an item's peek open: the answer to what waits;
   *  anything else takes none (the text stays). */
  const answerItem = (row: AgentRow | undefined, body: string): boolean => {
    if (!row) {
      deps.say("that item is not in the list any more (answered or dismissed elsewhere)")
      return false
    }
    if (!takesAnswer(row.item)) {
      deps.say(`this item takes no answer: ${km.label("agentsOpen")} opens it in a session`)
      return false
    }
    const item = row.item!
    setState({ hidden: [...state.hidden, row.id], peek: null })
    void deps.answer(item, body).then((ok) => {
      if (ok) {
        const where = item.workspace ? ` in ${shortPath(item.workspace)}` : ""
        deps.say(`answered: the next session started${where} gets it`)
        return
      }
      // Not taken (the note says why): back in the list, peeked again and
      // the text in the box -- enter tries the answer again, never a new
      // session with it.
      setState({ hidden: state.hidden.filter((h) => h !== row.id) })
      if (state.open && !state.peek) peekAt(row)
      if (!deps.draft()) deps.setDraft(body)
    })
    return true
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
    if (state.peek?.id.startsWith("attention:")) {
      return answerItem(rows().find((r) => r.id === state.peek!.id), body)
    }
    if (state.peek) {
      // A reply: to the session peeked at, whatever is selected now.
      const target = rows().find((r) => r.id === state.peek!.id)
      if (!target || target.status === "ended") {
        deps.say(`${target?.name ?? "that session"} has ended: → resumes it, then send this again`)
        return false                      // not sent: the text comes back
      }
      deps.message(target, body)
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
  // An item answered or dismissed elsewhere leaves the list: its peek goes
  // with it -- unless something is typed for it (enter then says so, and
  // keeps the text).
  createEffect(() => {
    const peek = state.peek
    if (!state.open || !peek) return
    if (!rows().some((r) => r.id === peek.id) && !deps.draft()) setState({ peek: null })
  })

  const cancelRename = () => {
    if (!state.renaming) return false
    endRename()
    return true
  }

  return { state, rows, selected, open, close, handleKey, submit, refreshPeek, cancelRename }
}

/** The hints line while the view is open: what the keys do now, as many
 *  as fit in `room` columns -- the most needed ones first. */
export function agentHelp(view: AgentView, km: Keymap, room = Infinity): string {
  if (view.state.renaming) return "type the new name · enter renames · esc back"
  const list = view.rows()
  const row = list[view.selected(list)]
  const peeked = view.state.peek ? list.find((r) => r.id === view.state.peek!.id) : undefined
  if (view.state.peek && !peeked) {
    return fitParts(["that item is gone: answered or dismissed elsewhere", "esc closes it", "↑↓ select"], room)
  }
  if (peeked?.kind === "attention") {
    const pick = [takesAnswer(peeked.item) ? "type the answer · enter sends it" : "",
                  `${km.label("agentsOpen")} open`, `${km.label("agentsStop")} dismiss`, "space closes", "esc back"]
      .filter(Boolean)
    return fitParts(pick, room)
  }
  if (view.state.peek) return "type to reply to it · enter sends · space closes · esc back"
  const item = row?.kind === "attention"
  const parts = [
    { at: 0, text: "↑↓ select" },
    { at: 1, text: `${km.label("agentsOpen")} open` },
    { at: 2, text: `${km.label("agentsPeek")} peek` },
    { at: 3, text: `${km.label("agentsStop")} ${item ? "dismiss" : "stop"}` },
    { at: 4, text: item ? "" : `${km.label("agentsRename")} rename` },
    { at: 5, text: item ? "" : `${km.label("agentsPin")} pin` },
    { at: 6, text: "esc back" },
    { at: 7, text: "type to start a new session" },
  ].filter((p) => p.text)
  const need = [7, 1, 2, 6, 3, 0, 4, 5]
  const kept: typeof parts = []
  for (const at of need) {
    const part = parts.find((p) => p.at === at)
    if (!part) continue
    const next = [...kept, part].sort((a, b) => a.at - b.at)
    if (width(next.map((p) => p.text).join(" · ")) > room) continue
    kept.splice(0, kept.length, ...next)
  }
  return kept.map((p) => p.text).join(" · ")
}

/** As many of `parts` as fit in `room` columns, in their order. */
function fitParts(parts: string[], room: number): string {
  const kept: string[] = []
  for (const part of parts) {
    if (width([...kept, part].join(" · ")) > room) break
    kept.push(part)
  }
  return kept.join(" · ")
}

/** An item's peek, wrapped: from its top, cut where there is no more room
 *  -- its last line (what → does) always shown. */
export function wrappedPeek(lines: string[], cols: number, room: number): string[] {
  const body = lines.slice(0, -1).flatMap((l) => wrap(l, Math.max(10, cols)))
  const last = lines.length ? wrap(lines[lines.length - 1]!, Math.max(10, cols)) : []
  const keep = Math.max(1, room - last.length)
  const head = body.length > keep ? [...body.slice(0, keep - 1), "…"] : body
  return [...head, ...last]
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
        const shown = peek.loading ? ["loading…"]
          : peek.wrap ? wrappedPeek(peek.lines, inner() - 8, room) : peek.lines.slice(-room)
        for (const line of shown) out.push({ kind: "peek", text: line })
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
  const glyph = (row: AgentRow) => row.kind === "attention" ? MARK_GLYPH[markOf(row.item!)]
    : row.status === "elsewhere" ? "◇"
    : row.status === "ended" ? "◼"
    : row.status === "sleeping" ? "☾"
    : row.status === "needs" ? "◉"
    : row.status === "working" ? SPINNER[props.tick() % SPINNER.length]!
    : "✓"
  const tone = (row: AgentRow) => row.kind === "attention"
      ? { asks: t().warn, failed: t().err, notice: t().faint }[markOf(row.item!)]
    : row.status === "elsewhere" ? t().muted
    : row.status === "ended" || row.status === "sleeping" ? t().faint
    : row.status === "needs" ? t().warn
    : row.status === "done" ? t().ok : t().accent
  const summary = () => {
    const all = v.rows()
    const list = all.filter((r) => r.kind === "session")
    // What needs the user: sessions here that wait, and every item of NEEDS YOU.
    const needs = list.filter((r) => r.status === "needs").length
      + all.filter((r) => r.kind === "attention" && r.group === "needs").length
    const working = list.filter((r) => r.status === "working").length
    return [`${list.length} session${list.length === 1 ? "" : "s"}`,
            needs ? `${needs} need${needs === 1 ? "s" : ""} you` : "",
            working ? `${working} working` : ""].filter(Boolean).join(" · ")
  }
  return (
    <box flexDirection="column" height={props.height()} flexShrink={0}>
      <text wrapMode="none">
        <span style={{ fg: t().accent, bold: true }}>{"◆ AGENTS  "}</span>
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
          // A notice is quiet: what it says in the faint colour.
          const quiet = row.kind === "attention" && markOf(row.item!) === "notice"
          const where = () => row.where ? middle(row.where, Math.max(8, Math.floor(inner() * 0.24))) : ""
          const right = () => [where(), row.age, row.cost].filter(Boolean).join("  ")
          const name = () => fit(row.name, Math.max(8, Math.floor(inner() * 0.32)))
          const stopLabel = props.keymap?.label("agentsStop") ?? "ctrl+x"
          const activity = () => renaming() ? "renaming: type the new name below"
            : armed() ? `${stopLabel} again ${row.kind === "attention" ? "dismisses" : "stops"} it` : row.activity
          const room = () => Math.max(0, inner() - 5 - width(name()) - 2 - (right() ? width(right()) + 2 : 0))
          return (
            <text wrapMode="none">
              <span style={{ fg: on ? t().accent : t().faint }}>{on ? " › " : "   "}</span>
              <span style={{ fg: tone(row) }}>{glyph(row) + " "}</span>
              <span style={{ fg: on ? t().bright : quiet ? t().muted : t().text, bold: Boolean(on || row.current) }}>{name()}</span>
              <span style={{ fg: armed() ? t().warn : renaming() ? t().accent : quiet ? t().faint : t().muted }}>{"  " + fit(activity(), room()).padEnd(room())}</span>
              <span style={{ fg: t().faint }}>{right() ? "  " + right() : ""}</span>
            </text>
          )
        }}
      </For>
    </box>
  )
}

/** A line of an answer as a peek shows it: without the markdown that
 *  marks headings, emphasis and code (a peek is plain text). */
export function plainLine(line: string): string {
  return line.trim()
    .replace(/^#{1,6}\s+/, "")
    .replace(/^[-*+]\s+/, "• ")
    .replace(/\*\*|__|`/g, "")
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
        for (const line of String(e.text ?? "").split("\n").map(plainLine).filter(Boolean).slice(0, 6)) out.push(line)
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
      case "context":
        out.push(`⏺ context  ${String(e.title ?? "")}`)
        break
    }
  }
  for (const asked of open.values()) out.push(`◉ waiting for you: ${askedWhat(asked as never)}`)
  return out.length ? out : ["(nothing yet)"]
}
