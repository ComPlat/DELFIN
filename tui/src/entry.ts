/**
 * How the window enters its first session, as `claude` does:
 *
 *   delfin-agent tui                a new session (the default, unchanged)
 *   delfin-agent tui --continue     the latest session of this folder
 *   delfin-agent tui --resume       the resume picker first
 *   delfin-agent tui --resume <id>  that session
 *
 * The launcher hands the flags over in DELFIN_TUI_START (cli_service.py
 * _start_settings): `resume` (an id), `continue`, and `pick` for a bare
 * `--resume`. DELFIN_TUI_ENTRY=continue|resume says the same from the
 * environment (for a launcher that does not pass the flag on).
 *
 * And the rows of the resume picker (the bare `--resume`, `/resume`, the
 * palette's "resume…"): every session this service runs, then the saved
 * ones of the folder -- each with its title, folder, age and state.
 */
import { shortPath } from "./agents/attention"

/** What a bare `--resume` hands over as its id (cli_service.RESUME_PICK). */
export const RESUME_PICK = ":pick"

export type Entry =
  | { kind: "new" }
  | { kind: "resume"; id: string }
  | { kind: "continue" }
  | { kind: "pick" }

/** The one line the first screen of a new session shows. */
export const START_HINT = "/resume opens another session · delfin-agent tui --continue opens the latest one here"

/** How the window starts, from its start settings and the environment. An
 *  id to resume wins over --continue, which wins over the picker. */
export function entryOf(start: Record<string, unknown>,
                        env: Record<string, string | undefined> = process.env): Entry {
  const resume = String(start.resume ?? "").trim()
  if (resume && resume !== RESUME_PICK) return { kind: "resume", id: resume }
  if (start.continue === true) return { kind: "continue" }
  if (start.pick === true || resume === RESUME_PICK) return { kind: "pick" }
  const said = String(env.DELFIN_TUI_ENTRY ?? "").trim().toLowerCase()
  if (said === "continue") return { kind: "continue" }
  if (said === "resume" || said === "pick") return { kind: "pick" }
  return { kind: "new" }
}

/** The start settings the service is given for `entry`: the picker's own
 *  markers left out (session.create knows no `pick`), `continue` set when
 *  only the environment asked for it. */
export function settingsFor(start: Record<string, unknown>, entry: Entry): Record<string, unknown> {
  const { pick: _pick, ...rest } = start
  if (rest.resume === RESUME_PICK || entry.kind !== "resume") delete rest.resume
  if (entry.kind === "continue") rest.continue = true
  else delete rest.continue
  return rest
}

// -- the resume picker's rows ---------------------------------------------------

/** running: a worker works or waits in the service; waiting: it waits for
 *  an answer; sleeping: kept by the service without a worker; saved: in
 *  the session store only. */
export type ResumeState = "running" | "waiting" | "sleeping" | "saved"

export interface ResumeRow {
  id: string
  workspace: string
  title: string
  state: ResumeState
  /** When it was last active, in seconds since the epoch (0: unknown). */
  at: number
  /** The picker's second column: state · folder · age · messages · short id. */
  detail: string
  /** What the filter reads besides title and detail: the whole id and
   *  the whole folder. */
  search: string
  /** The session in view. */
  current: boolean
}

/** "now", "5m ago", "3h ago", "2d ago" (from seconds since the epoch). */
export function ageOf(at: number, now: number = Date.now() / 1000): string {
  if (!at) return ""
  const d = Math.max(0, now - at)
  if (d < 60) return "now"
  if (d < 3600) return `${Math.floor(d / 60)}m ago`
  if (d < 86400) return `${Math.floor(d / 3600)}h ago`
  return `${Math.floor(d / 86400)}d ago`
}

/** A folder with the home as ~, cut in the middle to `width` characters
 *  (the end, the folder's own name, stays). */
export function shortFolder(path: string, width = 32, home: string = process.env.HOME ?? ""): string {
  const text = shortPath(String(path ?? ""), home)
  const chars = [...text]
  if (width <= 0 || chars.length <= width) return text
  const head = Math.floor((width - 1) / 3)
  return chars.slice(0, head).join("") + "…" + chars.slice(chars.length - (width - 1 - head)).join("")
}

function stateOf(summary: Record<string, unknown>): ResumeState {
  const state = String(summary.state ?? "")
  if (state === "sleeping" || state === "ended") return "sleeping"
  if (Number(summary.pending ?? 0) > 0 || state === "needs_input") return "waiting"
  return "running"
}

const firstLine = (text: unknown) => String(text ?? "").split("\n")[0]!.trim()

/**
 * The picker's rows: the sessions the service runs (sessions.list) on top,
 * then the sleeping ones and the saved ones of the folder (sessions.past)
 * -- each group newest first. A saved session the service runs or keeps is
 * listed once, as the service's.
 */
export function resumeRowsOf(listed: unknown[], past: unknown[], currentId = "",
                             now: number = Date.now() / 1000,
                             home: string = process.env.HOME ?? ""): ResumeRow[] {
  const rows: ResumeRow[] = []
  const seen = new Set<string>()
  const row = (id: string, workspace: string, title: string, state: ResumeState, at: number,
               extra: string[] = []): ResumeRow => ({
    id, workspace, title: title || "(untitled)", state, at, current: id === currentId,
    detail: [state === "waiting" ? "waits for you" : state, workspace ? shortFolder(workspace, 32, home) : "",
             ageOf(at, now), ...extra, id.slice(0, 8)].filter(Boolean).join(" · "),
    search: `${id} ${workspace}`,
  })
  for (const raw of listed) {
    const s = (raw ?? {}) as Record<string, unknown>
    const id = String(s.sessionId ?? "")
    if (!id || seen.has(id)) continue
    seen.add(id)
    const at = Number(s.lastActivity ?? s.createdAt ?? 0) || 0
    rows.push(row(id, String(s.workspace ?? ""), firstLine(s.title || s.firstPrompt), stateOf(s), at))
  }
  for (const raw of past) {
    const s = (raw ?? {}) as Record<string, unknown>
    const id = String(s.sessionId ?? "")
    if (!id || seen.has(id) || s.running) continue
    seen.add(id)
    const messages = Number(s.messages ?? 0)
    rows.push(row(id, String(s.workspace ?? ""), firstLine(s.title), "saved", Number(s.updatedAt ?? 0) || 0,
                  messages ? [`${messages} message${messages === 1 ? "" : "s"}`] : []))
  }
  const live = (r: ResumeRow) => (r.state === "running" || r.state === "waiting" ? 0 : 1)
  return rows.sort((a, b) => live(a) - live(b) || b.at - a.at)
}
