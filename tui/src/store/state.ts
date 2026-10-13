/**
 * What the TUI knows about its sessions, and how each event changes it.
 *
 * `applyEvent` is a plain function over a plain object, so it is tested
 * without Solid; the app runs it inside a Solid store's `produce` on a batch
 * of events (store.ts). What becomes scrollback is not state -- committed
 * rows cannot change -- and is the presenter's job (present.ts).
 */
import type { AnyEvent, SessionSummary } from "../protocol/generated"

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
  thinking: string           // the live reasoning (footer, when shown)
  lastSeq: number
  summary: SessionSummary | null
  ended: { reason: string; exitCode?: number } | null
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
  messages: string[]         // the window's own notes (config, connection)
}

export function emptyState(): AppState {
  return { connection: "connecting", mayAct: true, host: "", sessions: {}, order: [],
           current: null, elsewhere: [], messages: [] }
}

export function newSession(id: string, info: Record<string, unknown> = {}): SessionState {
  const s: SessionState = {
    id, title: "", workspace: "", provider: "", model: "", mode: "", approval: "",
    effort: "", state: "starting", pending: [], usage: { inputTokens: 0, outputTokens: 0, costUsd: 0 },
    turn: null, tail: "", thinking: "", lastSeq: 0, summary: null, ended: null,
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

function toolLabel(event: Record<string, any>): string {
  const input = (event.input ?? {}) as Record<string, unknown>
  const arg = ["path", "command", "pattern", "query", "url", "file"].map((k) => input[k])
    .find((v) => typeof v === "string" && v)
  return arg ? `${event.name}  ${String(arg).split("\n")[0]!.slice(0, 80)}` : String(event.name)
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
    return
  }
  if (type === "session.removed") {
    const s = state.sessions[id]
    if (!s) return
    s.state = "ended"
    s.ended = { reason: str(e.reason), exitCode: e.exitCode }
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
      s.state = "ended"
      s.ended = { reason: str(e.reason), exitCode: e.exitCode }
      break
    case "turn.started":
      s.turn = { id: str(e.turnId), source: str(e.source), startedAt: now, label: "thinking…" }
      s.tail = ""
      s.thinking = ""
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
