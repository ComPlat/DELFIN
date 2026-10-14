/**
 * What the TUI knows about its sessions, and how each event changes it.
 *
 * `applyEvent` is a plain function over a plain object, so it is tested
 * without Solid; the app runs it inside a Solid store's `produce` on a batch
 * of events (store.ts). What becomes scrollback is not state -- committed
 * rows cannot change -- and is the presenter's job (present.ts).
 */
import type { AnyEvent, SessionSummary } from "../protocol/generated"
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
}

export function emptyState(): AppState {
  return { connection: "connecting", mayAct: true, host: "", sessions: {}, order: [],
           current: null, elsewhere: [], starting: false }
}

export function newSession(id: string, info: Record<string, unknown> = {}): SessionState {
  const s: SessionState = {
    id, title: "", workspace: "", provider: "", model: "", mode: "", approval: "",
    effort: "", state: "starting", pending: [], usage: { inputTokens: 0, outputTokens: 0, costUsd: 0 },
    turn: null, tail: "", tailKey: "", held: [], thinking: "", lastSeq: 0, summary: null, ended: null, wait: "",
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
