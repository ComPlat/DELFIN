/**
 * What the TUI knows about its sessions, and how each event changes it.
 *
 * `applyEvent` is a plain function over a plain object, so it is tested
 * without Solid; the app runs it inside a Solid store's `produce` on a batch
 * of events (store.ts). What becomes scrollback is not state -- committed
 * rows cannot change -- and is the presenter's job (present.ts).
 */
import type { AnyEvent, AttentionItem, SessionSummary } from "../protocol/generated"
import { toolHeadline } from "../scrollback/tools"

export interface Asked {
  requestId: string
  kind: string
  tool: string
  preview: string
  options: { key: string; label: string; detail?: string }[]
  flags: Record<string, boolean>
  question?: string
  header?: string
  choices?: string[]
  choiceDetails?: string[]
  multiSelect?: boolean
  suggestion?: string
  readDir?: string
}

export interface TurnLive {
  id: string
  source: string
  startedAt: number          // ms, the window's clock
  label: string              // what the status line says it is doing
  /** The session's totals when the turn began: the status line shows this
   *  turn's tokens and cost, not the whole conversation's (as the REPL). */
  base: { inputTokens: number; outputTokens: number; costUsd: number }
  /** Answer characters streamed: an estimate of the output tokens until
   *  the provider reports them (OpenAI-compatible streams report at the end). */
  streamed: number
}

/** One row of the session's task list (session.info / session.updated
 *  `tasks`), in the order the session sends them. */
export interface TaskRow {
  num: number
  status: "pending" | "in_progress" | "blocked" | "completed"
  /** The model's text, one plain line: drawn as it is, never as markup. */
  label: string
}

const TASK_STATUSES = new Set(["pending", "in_progress", "blocked", "completed"])

/** The rows a session sent, each field as its type says: an unknown status
 *  is "pending", control and invisible characters become spaces (the
 *  service cleans them too; a window does not rely on it), at most 160
 *  characters. Anything that is not a row is left out. */
export function taskRows(raw: unknown): TaskRow[] {
  if (!Array.isArray(raw)) return []
  const out: TaskRow[] = []
  for (const r of raw as Record<string, unknown>[]) {
    if (!r || typeof r !== "object") continue
    const label = str(r.label)
      .replace(/[\u0000-\u001f\u007f-\u009f\u200b-\u200f\u2028-\u202e\u2060-\u206f\ufeff]/g, " ")
      .replace(/\s+/g, " ").trim().slice(0, 160)
    if (!label) continue
    const status = TASK_STATUSES.has(str(r.status)) ? str(r.status) as TaskRow["status"] : "pending"
    out.push({ num: Number(r.num) || 0, status, label })
  }
  return out
}

/** One piece of work the session still has out (session.info /
 *  session.updated `background`): a shell, a watched job, a background
 *  agent, a wake-up. No age: `since` (started) or `at` (fires), in epoch
 *  seconds, and the window counts the time itself. */
export interface BackgroundRow {
  group: string              // shells | watches | agents | wakeups | errors
  id: string
  kind: string               // "Shell", "Server", "Watch · CI", "Agent", "Wake-up"
  /** What it is: a description or the command, one plain line. */
  label: string
  url: string
  state: string
  since: number | null
  at: number | null
}

const one = (value: unknown, max: number) => str(value)
  .replace(/[\u0000-\u001f\u007f-\u009f\u200b-\u200f\u2028-\u202e\u2060-\u206f\ufeff]/g, " ")
  .replace(/\s+/g, " ").trim().slice(0, max)
const seconds = (value: unknown) => typeof value === "number" && Number.isFinite(value) ? value : null

/** The rows a session sent, each field as its type says (control and
 *  invisible characters made spaces); a row without an id or a label is
 *  left out. */
export function backgroundRows(raw: unknown): BackgroundRow[] {
  if (!Array.isArray(raw)) return []
  const out: BackgroundRow[] = []
  for (const r of raw as Record<string, unknown>[]) {
    if (!r || typeof r !== "object") continue
    const id = one(r.id, 64)
    const label = one(r.label, 160) || id
    if (!id) continue
    out.push({ group: one(r.group, 20), id, kind: one(r.kind, 40) || "Job", label, url: one(r.url, 200),
               state: one(r.state, 40), since: seconds(r.since), at: seconds(r.at) })
  }
  return out
}

export interface SessionState {
  id: string
  title: string
  workspace: string
  provider: string
  model: string
  mode: string
  approval: string
  effort: string
  state: string              // idle | running | needs_input | ended | starting
  pending: Asked[]
  usage: { inputTokens: number; outputTokens: number; costUsd: number }
  turn: TurnLive | null
  tail: string               // the answer's unfinished block (footer)
  /** Which block of the answer the tail is (its part and number): the
   *  footer draws a block with one element from its first delta until its
   *  rows are drawn. */
  tailKey: string
  /** Blocks of the answer that are finished but whose rows are not drawn
   *  yet (their highlighting, a commit): still shown with the tail, so a
   *  block never leaves the screen on its way into the rows. */
  held: { key: string; text: string }[]
  thinking: string           // the live reasoning (footer, when shown)
  lastSeq: number
  summary: SessionSummary | null
  ended: { reason: string; exitCode?: number } | null
  /** A progress line the next one replaces (waiting out an endpoint outage). */
  wait: string
  /** The session's task list, as it sent it last. */
  tasks: TaskRow[]
  /** Why shift+tab steps over bypassPermissions ("" when it may stop
   *  there); null until the session has said it -- a service that never
   *  says is treated as stepping over it. */
  bypassBlocked: string | null
  /** Bypass was confirmed once in this session: the cycle passes through
   *  it without asking again. */
  bypassConfirmed: boolean
  /** Auto mode's one-line state ("" while off), as the session sent it. */
  auto: string
  /** What the session still has out, as it sent it last. */
  background: BackgroundRow[]
}

export interface Elsewhere {
  key: string
  sessionId?: string | null
  title?: string | null
  workspace: string
  host: string
}

export interface AppState {
  connection: "connecting" | "up" | "reconnecting" | "closed"
  mayAct: boolean
  host: string
  sessions: Record<string, SessionState>
  order: string[]
  current: string | null
  elsewhere: Elsewhere[]
  /** The first session is being started (the window is up already). */
  starting: boolean
  /** What needs the user from every DELFIN process -- parked questions,
   *  failed jobs, alerts, notices (the service's attention feed), in the
   *  service's order. Empty with a service before serviceVersion 3. */
  attention: AttentionItem[]
  /** The version of that list (the service counts its changes): a list
   *  and an event that crossed on the socket -- the newer one stays. -1
   *  until the service of this connection has said one. */
  attentionVersion: number
}

export function emptyState(): AppState {
  return { connection: "connecting", mayAct: true, host: "", sessions: {}, order: [],
           current: null, elsewhere: [], starting: false, attention: [], attentionVersion: -1 }
}

/** The attention items a service sent (attention.list, attention.updated),
 *  each field as its type says -- an item that has no id is left out. */
export function attentionItems(raw: unknown): AttentionItem[] {
  if (!Array.isArray(raw)) return []
  const out: AttentionItem[] = []
  for (const r of raw as Record<string, unknown>[]) {
    if (!r || typeof r !== "object" || !r.id) continue
    out.push({
      id: str(r.id), kind: str(r.kind) || "notice", group: r.group === "needs" ? "needs" : "notices",
      blocking: Boolean(r.blocking), name: str(r.name), text: str(r.text), detail: str(r.detail),
      workspace: str(r.workspace), host: str(r.host), sessionId: str(r.sessionId),
      sessionTitle: str(r.sessionTitle), sessionOpen: Boolean(r.sessionOpen), sessionSaved: Boolean(r.sessionSaved),
      createdAt: Number(r.createdAt ?? 0) || 0,
      options: Array.isArray(r.options) ? r.options.map(str) : [],
      job: r.job && typeof r.job === "object" ? r.job as Record<string, unknown> : null,
    })
  }
  return out
}

export function newSession(id: string, info: Record<string, unknown> = {}): SessionState {
  const s: SessionState = {
    id, title: "", workspace: "", provider: "", model: "", mode: "", approval: "",
    effort: "", state: "starting", pending: [], usage: { inputTokens: 0, outputTokens: 0, costUsd: 0 },
    turn: null, tail: "", tailKey: "", held: [], thinking: "", lastSeq: 0, summary: null, ended: null, wait: "",
    tasks: [], bypassBlocked: null, bypassConfirmed: false, auto: "", background: [],
  }
  applyInfo(s, info)
  return s
}

function str(value: unknown): string {
  return typeof value === "string" ? value : value == null ? "" : String(value)
}

export function applyInfo(s: SessionState, info: Record<string, unknown>): void {
  if ("title" in info && info.title) s.title = str(info.title)
  if ("workspace" in info) s.workspace = str(info.workspace)
  if ("provider" in info) s.provider = str(info.provider)
  if ("model" in info) s.model = str(info.model)
  if ("mode" in info) s.mode = str(info.mode)
  if ("approval" in info) s.approval = str(info.approval)
  if ("effort" in info) s.effort = str(info.effort)
  if ("tasks" in info && info.tasks != null) s.tasks = taskRows(info.tasks)
  if ("bypassBlocked" in info && info.bypassBlocked != null) s.bypassBlocked = str(info.bypassBlocked)
  if ("bypassConfirmed" in info && info.bypassConfirmed != null) s.bypassConfirmed = Boolean(info.bypassConfirmed)
  if ("auto" in info) s.auto = oneLine(str(info.auto))
  if ("background" in info && info.background != null) s.background = backgroundRows(info.background)
}

/** The session's text as one plain line: control and invisible characters
 *  made spaces, runs of space one, at most 160 characters. */
function oneLine(text: string): string {
  return text.replace(/[\u0000-\u001f\u007f-\u009f\u200b-\u200f\u2028-\u202e\u2060-\u206f\ufeff]/g, " ")
    .replace(/\s+/g, " ").trim().slice(0, 160)
}

/** Register a session (create, attach, a summary for one not seen yet). */
export function ensureSession(state: AppState, id: string,
                              info: Record<string, unknown> = {}): SessionState {
  let s = state.sessions[id]
  if (!s) {
    s = newSession(id, info)
    state.sessions[id] = s
    state.order.push(id)
  } else {
    applyInfo(s, info)
  }
  return s
}

/** A session that ended -- on its own, or because its worker died: nothing
 *  runs any more and nothing is waiting for an answer. */
function endSession(s: SessionState, e: Record<string, any>): void {
  s.state = "ended"
  s.ended = { reason: str(e.reason), exitCode: e.exitCode }
  s.turn = null
  s.tail = ""
  s.thinking = ""
  s.pending = []
  s.wait = ""
}

/** What the status line says a running tool does: its headline. */
function toolLabel(event: Record<string, any>): string {
  const head = toolHeadline(String(event.name ?? "tool"), (event.input ?? {}) as Record<string, unknown>)
  return head.arg ? `${head.name}  ${head.arg.slice(0, 80)}` : head.name
}

/** Apply one event. Unknown events and events for unknown sessions are
 *  ignored (a summary registers a session; nothing else does). */
export function applyEvent(state: AppState, event: AnyEvent, now: number = Date.now()): void {
  const e = event as Record<string, any>
  const type = String(e.type ?? "")
  const id = String(e.session ?? "")
  if (type === "session.summary") {
    const summary = e.summary as SessionSummary
    const s = ensureSession(state, summary.sessionId, summary as unknown as Record<string, unknown>)
    s.summary = summary
    s.state = summary.state
    // Listed as running (resumed, woken): what ended it is history.
    if (summary.state !== "ended") s.ended = null
    return
  }
  if (type === "session.removed") {
    const s = state.sessions[id]
    if (!s) return
    endSession(s, e)
    return
  }
  if (type === "attention.updated") {
    const version = typeof e.version === "number" ? e.version : null
    if (version !== null && version < state.attentionVersion) return     // older than what is shown
    state.attention = attentionItems(e.items)
    if (version !== null) state.attentionVersion = version
    return
  }
  const s = state.sessions[id]
  if (!s) return
  if (typeof e.seq === "number" && e.seq > s.lastSeq) s.lastSeq = e.seq
  switch (type) {
    case "session.info":
    case "session.updated":
      applyInfo(s, e)
      break
    case "session.state":
      s.state = str(e.state)
      break
    case "session.ended":
      endSession(s, e)
      break
    case "turn.started":
      s.turn = { id: str(e.turnId), source: str(e.source), startedAt: now, label: "thinking…",
                 base: { ...s.usage }, streamed: 0 }
      s.tail = ""
      s.thinking = ""
      s.wait = ""
      break
    case "text.delta":
      if (s.turn) s.turn.streamed += str(e.delta).length
      s.wait = ""
      break
    case "wait":
      s.wait = str(e.text)
      break
    case "reasoning.delta":
      s.thinking += str(e.delta)
      if (s.turn) s.turn.label = "thinking…"
      break
    case "reasoning.ended":
      s.thinking = ""
      break
    case "text.started":
      if (s.turn) s.turn.label = "writing…"
      break
    case "tool.started":
      if (s.turn) s.turn.label = toolLabel(e)
      break
    case "tool.ended":
    case "tool.denied":
      if (s.turn) s.turn.label = "thinking…"
      break
    case "turn.ended":
      s.turn = null
      s.tail = ""
      s.thinking = ""
      s.wait = ""
      break
    case "usage":
      // Session totals, sampled during a turn.
      s.usage.inputTokens = Number(e.inputTokens ?? s.usage.inputTokens)
      s.usage.outputTokens = Number(e.outputTokens ?? s.usage.outputTokens)
      s.usage.costUsd = Number(e.costUsd ?? s.usage.costUsd)
      break
    case "permission.asked":
      if (!s.pending.some((p) => p.requestId === e.requestId)) s.pending.push(e as unknown as Asked)
      break
    case "permission.resolved":
      s.pending = s.pending.filter((p) => p.requestId !== e.requestId)
      break
  }
}
