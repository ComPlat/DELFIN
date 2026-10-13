/**
 * The session in view, as scrollback: each event becomes the rows a person
 * reads, in the order it happened.
 *
 * Answers are committed one finished markdown block at a time while the
 * unfinished block stays live in the footer (`tail`); a replay, which has
 * no deltas, commits each answer's blocks from its final text. Tool calls
 * become one row when they end (the footer shows them while they run).
 * An approval leaves one line saying how it was answered.
 */
import type { AnyEvent } from "../protocol/generated"
import type { AppState } from "../store/state"
import { AnswerBlocks, splitMarkdown } from "./blocks"
import type { Transcript } from "./transcript"

type Fields = Record<string, any>

export interface PresentOptions {
  showThinking: boolean
}

const DIFF_TOOLS = new Set(["edit_file", "write_file", "apply_patch", "multi_edit", "str_replace"])

function compact(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`
  if (n >= 1_000) return `${(n / 1_000).toFixed(1)}k`
  return String(n)
}

function firstLine(text: string, max = 100): string {
  const line = text.split("\n").map((l) => l.trim()).find((l) => l) ?? ""
  return line.length > max ? line.slice(0, max - 1) + "…" : line
}

export function toolArg(input: Fields): string {
  for (const key of ["path", "command", "pattern", "query", "url", "file", "name"]) {
    const value = input?.[key]
    if (typeof value === "string" && value.trim()) return firstLine(value, 90)
  }
  return ""
}

/** The one-line result of a tool call, and a short tail when it ran a command. */
export function toolSummary(output: string, meta: Fields): { summary: string; tail: string[]; ok: boolean } {
  let ok = meta?.ok !== false && !meta?.error
  let text = output ?? ""
  let tail: string[] = []
  const trimmed = text.trim()
  if (trimmed.startsWith("{") && trimmed.includes("exit_code")) {
    try {
      const run = JSON.parse(trimmed)
      const code = Number(run.exit_code ?? 0)
      ok = ok && code === 0
      const out = String(run.stdout ?? "") + (run.stderr ? `\n${run.stderr}` : "")
      const lines = out.split("\n").filter((l) => l.trim())
      tail = lines.slice(-3).map((l) => firstLine(l, 110))
      return { summary: `exit ${code}${lines.length ? ` · ${lines.length} line${lines.length === 1 ? "" : "s"}` : ""}`, tail, ok }
    } catch {
      // not the shell's envelope after all
    }
  }
  const lines = text.split("\n").filter((l) => l.trim())
  const size = Number(meta?.chars ?? text.length)
  let summary = lines.length <= 1 ? firstLine(text) || "(no output)"
    : `${lines.length} lines · ${firstLine(text, 70)}`
  if (meta?.truncated) summary += ` · ${compact(size)} chars`
  if (meta?.error) summary = String(meta.error)
  return { summary, tail, ok }
}

export class Presenter {
  private blocks: AnswerBlocks | null = null
  private tools = new Map<string, Fields>()
  private asked = new Map<string, Fields>()
  private bannered = false
  private n = 0

  constructor(readonly sessionId: string, private readonly transcript: Transcript,
              private readonly options: PresentOptions) {}

  private add(kind: string, fields: Fields): void {
    void this.transcript.append({ id: `${this.sessionId}:${++this.n}`, kind, ...fields })
  }

  /** Every event of a batch; only this session's are shown. */
  handle(events: AnyEvent[], draft: AppState): void {
    for (const event of events) {
      const e = event as Fields
      if (e.session === this.sessionId) this.one(e, draft)
    }
  }

  private one(e: Fields, draft: AppState): void {
    const s = draft.sessions[this.sessionId]
    const setTail = (tail: string) => {
      if (s) s.tail = tail
    }
    switch (e.type) {
      case "session.info":
        if (!this.bannered) {
          this.bannered = true
          const lines = (e.banner as string[] | undefined)?.filter((l) => l && !l.startsWith("delfin-agent")) ?? []
          const head = [e.provider && e.model ? `${e.provider}/${e.model}` : e.model, e.title]
            .filter(Boolean).join(" · ")
          this.add("banner", { head, lines })
          for (const notice of (e.notices as string[] | undefined) ?? []) this.add("notice", { text: notice })
        }
        break
      case "user.message":
        this.add("user", { text: e.text, source: e.source })
        break
      case "text.started":
        this.blocks = new AnswerBlocks((block) => this.add("markdown", { text: block }), setTail)
        break
      case "text.delta":
        if (!this.blocks) this.blocks = new AnswerBlocks((block) => this.add("markdown", { text: block }), setTail)
        this.blocks.append(String(e.delta ?? ""))
        break
      case "text.ended":
        if (this.blocks) {
          this.blocks.finish()
          this.blocks = null
        } else {
          for (const block of splitMarkdown(String(e.text ?? ""), true).blocks) {
            this.add("markdown", { text: block })
          }
        }
        setTail("")
        break
      case "reasoning.ended":
        if (this.options.showThinking && String(e.text ?? "").trim()) {
          this.add("reasoning", { text: firstLine(String(e.text), 400) })
        }
        break
      case "tool.started":
        this.tools.set(String(e.callId), e)
        break
      case "tool.ended": {
        const started = this.tools.get(String(e.callId)) ?? {}
        this.tools.delete(String(e.callId))
        const input = (started.input ?? {}) as Fields
        const output = String(e.output ?? "")
        if (DIFF_TOOLS.has(String(e.name)) && /^(---|diff --git|\+\+\+|@@)/m.test(output)) {
          this.add("diff", { name: e.name, path: toolArg(input), text: output, filetype: filetypeOf(toolArg(input)) })
        } else {
          const { summary, tail, ok } = toolSummary(output, (e.meta ?? {}) as Fields)
          this.add("tool", { name: e.name, arg: toolArg(input), summary, tail, ok })
        }
        break
      }
      case "tool.denied":
        this.add("denied", { name: e.name })
        break
      case "notice":
        this.add("notice", { text: e.text })
        break
      case "error":
        this.add("error", { text: e.text })
        break
      case "permission.asked":
        this.asked.set(String(e.requestId), e)
        break
      case "permission.resolved": {
        const asked = this.asked.get(String(e.requestId))
        this.asked.delete(String(e.requestId))
        const what = asked ? askedWhat(asked) : "question"
        const ok = e.decision === "approved" || e.decision === "answered"
        const how = e.by === "abort" ? "refused (the turn ended)" : e.by === "expired" ? "expired" : String(e.decision)
        this.add("approval", { text: `${how}: ${what.trim()}`, ok })
        break
      }
      case "turn.ended": {
        const secs = (Number(e.durationMs ?? 0) / 1000).toFixed(1)
        if (e.error) {
          this.add("turnEnd", { text: `✗ ${secs}s · ${firstLine(String(e.error), 120)}`, ok: false })
        } else {
          const cost = Number(e.costUsd ?? 0)
          const parts = [`✓ ${secs}s`, `↑${compact(Number(e.inputTokens ?? 0))} ↓${compact(Number(e.outputTokens ?? 0))}`]
          if (cost > 0) parts.push(`$${cost.toFixed(4)}`)
          this.add("turnEnd", { text: parts.join(" · "), ok: true })
        }
        this.blocks = null
        break
      }
      case "session.ended":
        this.add("notice", { text: `session ended (${e.reason ?? "closed"})` })
        break
    }
  }
}

/** What a question was about, in a few words: the command, the file, the
 *  question itself. */
export function askedWhat(asked: Fields): string {
  if (asked.kind === "ask") return String(asked.question ?? "question")
  if (asked.kind === "plan") return "the plan"
  const preview = String(asked.preview ?? "")
  if (asked.flags?.shell) return `${asked.tool} ${firstLine(preview.replace(/^\$\s*/, ""), 80)}`
  const target = /^\+\+\+ (?:b\/)?(.+)$/m.exec(preview)?.[1] ?? /^--- (?:a\/)?(.+)$/m.exec(preview)?.[1]
  if (target) return `${asked.tool} ${target.trim()}`
  return `${asked.tool} ${firstLine(preview, 70)}`.trim()
}

export function filetypeOf(path: string): string {
  const ext = path.split(".").pop()?.toLowerCase() ?? ""
  return ({ py: "python", ts: "typescript", tsx: "tsx", js: "javascript", md: "markdown",
            sh: "bash", json: "json", zig: "zig" } as Record<string, string>)[ext] ?? "text"
}
