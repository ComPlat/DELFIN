/**
 * What needs the user from every DELFIN process, as the agent view shows
 * it: the service's attention items -- parked questions and
 * confirmations, failed and finished jobs, security alerts, budget
 * notices, written by dashboards, the REPL, the job monitor and the
 * scheduler -- as rows of NEEDS YOU and NOTICES, what a peek at one says,
 * and what → does with it.
 *
 * → opens every item in this window (main.tsx `openItem`):
 *   - its session, when the store keeps it: resumed -- or, while another
 *     process still has it open, a copy (two programs writing one
 *     conversation would lose one side's saves);
 *   - a failed job the job monitor diagnosed: the diagnosis session;
 *   - anything else: a new session in its folder with the item as its
 *     brief (one row, ctrl+o shows it all), and a prompt suggested in the
 *     prompt box. Nothing is sent until enter.
 */
import type { AttentionItem } from "../protocol/generated"
import type { AppState } from "../store/state"

/** How old something is, in its largest unit: 12s, 4m, 2h, 3d. */
export function ageOf(seconds: number, now: number): string {
  const d = Math.max(0, now / 1000 - seconds)
  if (d < 60) return `${Math.floor(d)}s`
  if (d < 3600) return `${Math.floor(d / 60)}m`
  if (d < 86400) return `${Math.floor(d / 3600)}h`
  return `${Math.floor(d / 86400)}d`
}

/** asks: someone waits for an answer; failed: a job, a run, an alert;
 *  notice: what happened and needs nothing. */
export type ItemMark = "asks" | "failed" | "notice"

export function markOf(item: AttentionItem): ItemMark {
  if (item.blocking) return "asks"
  if (item.kind === "failed" || item.kind === "alert") return "failed"
  return "notice"
}

export const MARK_GLYPH: Record<ItemMark, string> = { asks: "◉", failed: "✗", notice: "·" }

/** A path with the home folder as ~. */
export function shortPath(path: string, home: string = process.env.HOME ?? ""): string {
  if (home && (path === home || path.startsWith(home.endsWith("/") ? home : `${home}/`))) {
    return `~${path.slice(home.replace(/\/$/, "").length)}`
  }
  return path
}

/** Where an item is: open in another process (and on which machine), else
 *  its folder -- and the machine, when it is not this one. */
export function whereOf(item: AttentionItem, thisHost: string): string {
  const host = item.host && item.host !== thisHost ? item.host : ""
  if (item.sessionOpen) return host ? `open on ${host}` : "open elsewhere"
  return [item.workspace ? shortPath(item.workspace) : "", host].filter(Boolean).join(" · ")
}

/** A job's fields (the job monitor's finding, or what the inbox title said). */
export interface JobOf {
  id: string
  state: string
  folder: string
  signatures: string[]
  diagnosisSession: string
  summary: string
  recommendation: string
  /** The prepared fix in a line ("maxcore 3000 -> 4500"), and its diff. */
  fix: string
  diff: string
}

export function jobOf(item: AttentionItem): JobOf | null {
  const j = item.job as Record<string, any> | null | undefined
  if (!j || !j.id) return null
  const assessed = (j.fix && typeof j.fix === "object" ? j.fix : {}) as Record<string, any>
  const prepared = (assessed.fix && typeof assessed.fix === "object" ? assessed.fix : {}) as Record<string, any>
  return {
    id: String(j.id), state: String(j.state ?? ""), folder: String(j.folder ?? "") || item.workspace,
    signatures: Array.isArray(j.signatures) ? j.signatures.map(String) : [],
    diagnosisSession: String(j.diagnosisSession ?? ""), summary: String(j.summary ?? ""),
    recommendation: String(assessed.recommendation ?? ""),
    fix: String(prepared.summary_line ?? prepared.kind ?? ""), diff: String(prepared.diff ?? ""),
  }
}

/** What → does with an item. */
export type Opening =
  /** Its session, resumed here -- a copy while another process has it. */
  | { how: "session"; id: string; copy: boolean; title: string }
  /** The job monitor's diagnosis of a failed job (a saved session). */
  | { how: "diagnosis"; id: string; title: string }
  /** A new session about it, which waits for the user. */
  | { how: "new" }

export function openingOf(item: AttentionItem): Opening {
  const job = jobOf(item)
  if (item.kind === "failed" && job?.diagnosisSession) {
    return { how: "diagnosis", id: job.diagnosisSession, title: `the diagnosis of job ${job.id}` }
  }
  if (item.sessionId && item.sessionSaved) {
    return { how: "session", id: item.sessionId, copy: Boolean(item.sessionOpen),
             title: item.sessionTitle || `session ${item.sessionId.slice(0, 8)}` }
  }
  return { how: "new" }
}

/** What → does, in words (the last line of a peek). */
export function openingText(item: AttentionItem, thisHost: string): string {
  const o = openingOf(item)
  if (o.how === "diagnosis") return "→ opens the job monitor's diagnosis of it"
  if (o.how === "session") {
    if (!o.copy) return `→ resumes ${o.title} here`
    const host = item.host && item.host !== thisHost ? ` on ${item.host}` : ""
    return `→ opens a copy of ${o.title} (it is still open${host})`
  }
  const folder = item.workspace ? ` in ${shortPath(item.workspace)}` : ""
  return `→ starts a session about it${folder}, which waits for you`
}

function when(seconds: number, now: number): string {
  if (!seconds) return ""
  return `${ageOf(seconds, now)} ago`
}

/** A peek at an item: everything it says, a line each (wrapped where it
 *  is drawn), then what the keys do with it. */
export function itemPeek(item: AttentionItem, thisHost: string, now: number = Date.now()): string[] {
  const out: string[] = []
  // "in ~/calc on login02", "in ~/calc", "on login02"
  const where = [item.workspace ? `in ${shortPath(item.workspace)}` : "", item.host ? `on ${item.host}` : ""]
    .filter(Boolean).join(" ")
  const job = jobOf(item)
  if (item.kind === "question") {
    out.push(item.detail || item.text.replace(/^question:\s*/, ""))
    item.options.forEach((o, i) => out.push(`  ${i + 1}. ${o}`))
    const who = item.sessionTitle ? `asked by ${item.sessionTitle}` : "asked by a session"
    out.push([who, where, when(item.createdAt, now)].filter(Boolean).join(" · "))
    out.push(`type the answer and press enter: it goes to the next session started in ${
      item.workspace ? shortPath(item.workspace) : "its folder"}`)
  } else if (job) {
    out.push([`job ${job.id}`, job.state, job.folder ? shortPath(job.folder) : ""].filter(Boolean).join(" · "))
    if (job.signatures.length) out.push(`errors: ${job.signatures.join(", ")}`)
    if (job.summary) out.push(job.summary)
    if (job.fix) {
      out.push(`prepared fix: ${job.fix} (not applied: the dashboard's /fix applies it)`)
    } else if (job.recommendation) {
      out.push(`recommendation: ${job.recommendation}`)
    }
    if (item.detail && !job.summary) out.push(...item.detail.split("\n").filter((l) => l.trim()).slice(0, 6))
    out.push([item.host ? `on ${item.host}` : "", when(item.createdAt, now)].filter(Boolean).join(" · "))
  } else {
    // Its row says its name: the peek, what it says in full.
    for (const line of (item.detail || item.text).split("\n").filter((l) => l.trim()).slice(0, 8)) out.push(line)
    const from = item.sessionTitle ? `from ${item.sessionTitle}` : ""
    out.push([from, where, when(item.createdAt, now)].filter(Boolean).join(" · "))
  }
  out.push(openingText(item, thisHost))
  return out.filter((l) => l.trim())
}

/** What a new session about an item is given to read (session.create's
 *  brief): the item in full, a field a line. */
export function briefOf(item: AttentionItem, now: number = Date.now()): string {
  const job = jobOf(item)
  const lines: string[] = []
  const kind: Record<string, string> = {
    question: "A question a DELFIN session asked, which nobody has answered yet.",
    confirm: "A DELFIN session asked for a confirmation, which nobody has given yet.",
    plan: "A DELFIN session put up a plan for approval, which nobody has answered yet.",
    alert: "A security alert from DELFIN's sandbox.",
    budget: "A budget warning from a DELFIN run.",
  }
  if (job) {
    lines.push(`Job ${job.id} ${item.kind === "finished" ? "finished" : "failed"}${job.state ? ` (${job.state})` : ""}.`)
    if (job.folder) lines.push(`Folder: ${job.folder}`)
  } else {
    lines.push(kind[item.kind] ?? `${item.name}.`)
    if (kind[item.kind]) lines.push(`About: ${item.name}`)
    if (item.workspace) lines.push(`Folder: ${item.workspace}`)
  }
  if (item.host) lines.push(`Machine: ${item.host}`)
  if (item.createdAt) {
    const at = new Date(item.createdAt * 1000)
    lines.push(`When: ${at.toISOString().replace("T", " ").slice(0, 16)} UTC (${when(item.createdAt, now)})`)
  }
  if (item.sessionTitle || item.sessionId) {
    lines.push(`From the session: ${item.sessionTitle || item.sessionId}${item.sessionId ? ` (${item.sessionId})` : ""}`)
  }
  if (item.kind === "question") {
    lines.push(`The question: ${item.detail || item.text.replace(/^question:\s*/, "")}`)
    if (item.options.length) lines.push(`Its choices: ${item.options.join(" | ")}`)
  }
  if (job) {
    if (job.signatures.length) lines.push(`Errors found in its output: ${job.signatures.join(", ")}`)
    if (job.summary) lines.push(`The job monitor's diagnosis: ${job.summary}`)
    if (job.recommendation) lines.push(`Recommendation: ${job.recommendation}`)
    if (job.fix) {
      lines.push(`A fix is prepared (not applied): ${job.fix}`)
      if (job.diff) lines.push(job.diff)
    }
  }
  if (item.detail && item.kind !== "question" && item.detail !== job?.summary) {
    lines.push("", "Details:", item.detail)
  }
  return lines.join("\n")
}

/** The prompt a new session about an item suggests (in the prompt box;
 *  sent only when the user presses enter). */
export function suggestedPrompt(item: AttentionItem): string {
  const job = jobOf(item)
  if (job && item.kind === "failed") {
    return `Job ${job.id} failed${job.state ? ` with ${job.state}` : ""}. What went wrong, and what should we change?`
  }
  if (job && item.kind === "finished") return `Job ${job.id} finished. Summarise its results.`
  switch (item.kind) {
    case "alert": return "What does this alert mean, and does anything need to be done?"
    case "budget": return "What used up the budget, and how should we go on?"
    case "question": return "Help me answer this question."
    default: return "What happened here, and is anything to be done?"
  }
}

/** The quiet count at the bottom right: what waits for an answer (here
 *  and anywhere else) and what failed. Notices are not counted. */
export function attentionCounts(state: AppState, waitingHere: number): { asks: number; failed: number } {
  let asks = waitingHere
  let failed = 0
  for (const item of state.attention) {
    const mark = markOf(item)
    if (mark === "asks") asks++
    else if (mark === "failed") failed++
  }
  return { asks, failed }
}

/** How a new item names itself at the bottom right for a few seconds. */
export function arrivalNote(item: AttentionItem): string {
  const mark = markOf(item)
  if (mark === "asks") {
    const who = item.sessionTitle || item.name
    return `◉ ${who} ${item.kind === "question" ? "asks" : "waits for you"}`
  }
  return `✗ ${item.name}`
}
