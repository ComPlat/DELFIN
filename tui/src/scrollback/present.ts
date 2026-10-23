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
import type { AppState, SessionState } from "../store/state"
import { AnswerBlocks, splitMarkdown } from "./blocks"
import { resultBody, shortToolName, stripControl, toolHeadline, toolResult } from "./tools"
import type { Transcript } from "./transcript"

type Fields = Record<string, any>

export interface PresentOptions {
  showThinking: boolean
  /** A message as the user typed it, when this window sent it with its
   *  pastes put back (the placeholder reads better than the paste). */
  shownAs?: (text: string) => string | undefined
  /** A finished block's rows are drawn: it leaves the footer (the
   *  session's `held`). */
  release?: (session: string, key: string) => void
}

function compact(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`
  if (n >= 1_000) return `${(n / 1_000).toFixed(1)}k`
  return String(n)
}

function firstLine(text: string, max = 100): string {
  const line = text.split("\n").map((l) => l.trim()).find((l) => l) ?? ""
  return line.length > max ? line.slice(0, max - 1) + "…" : line
}

export class Presenter {
  /** The answer parts streaming now, by partId. */
  private parts = new Map<string, AnswerBlocks>()
  /** Parts seen to end: a delta of one (after a reconnect, the live stream
   *  goes on where the replay ended) is not shown again. */
  private endedParts = new Set<string>()
  private tools = new Map<string, Fields>()
  private asked = new Map<string, Fields>()
  private bannered = false
  private ended = false
  private lastError = ""
  private lastOutput: { name: string; output: string } | null = null
  private n = 0
  /** The session's state in the batch being handled: an answer's blocks
   *  finish over many batches. */
  private s: SessionState | undefined
  private static made = 0
  private readonly no = ++Presenter.made

  /** `groups` numbers the turns (shared by the window's presenters): a
   *  replay capped by turns keeps the newest. */
  constructor(readonly sessionId: string, private readonly transcript: Transcript,
              private readonly options: PresentOptions, private readonly groups: { n: number } = { n: 0 }) {}

  /** Resolves once the item's rows are drawn. */
  private add(kind: string, fields: Fields): Promise<void> {
    // A message that starts a turn starts its group; a queued or steered
    // one stays in the turn it was sent into.
    if (kind === "user" && (fields.source === "user" || fields.source === "harness" || !fields.source)) this.groups.n++
    return this.transcript.append({ id: `${this.sessionId}:${++this.n}`, kind, group: this.groups.n, ...fields })
  }

  /** A block of a streamed answer is final: into the rows -- and on in the
   *  footer, under the key it streamed with, until its rows are drawn. */
  private hold(text: string, key: string): void {
    this.s?.held?.push({ key, text })
    // The writer lets go of it as its rows are drawn (in one step with
    // them); this is for a write that failed or drew nothing.
    void this.add("markdown", { text, held: key, session: this.sessionId })
      .finally(() => this.options.release?.(this.sessionId, key))
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
    this.s = s
    const setTail = (tail: string) => {
      if (s) s.tail = tail
    }
    // Said to have ended, and it goes on (woken, resumed: what its new run
    // logs -- the store lets only events not seen yet through): its next
    // end is the window's to say again.
    if (this.ended && typeof e.seq === "number" && e.seq > 0 && e.type !== "session.ended") this.ended = false
    switch (e.type) {
      case "session.info":
        if (!this.bannered) {
          this.bannered = true
          this.add("banner", bannerOf(e))
          for (const notice of (e.notices as string[] | undefined) ?? []) this.add("notice", { text: notice })
        }
        break
      case "user.message":
        this.add("user", { text: this.options.shownAs?.(String(e.text ?? "")) ?? e.text, source: e.source })
        break
      case "text.started":
        this.part(String(e.partId ?? ""))
        break
      case "text.delta": {
        const id = String(e.partId ?? "")
        if (this.endedParts.has(id)) break
        this.part(id)?.append(String(e.delta ?? ""))
        break
      }
      case "text.ended": {
        const id = String(e.partId ?? "")
        if (this.endedParts.has(id)) break            // shown already
        this.endedParts.add(id)
        const streamed = this.parts.get(id)
        this.parts.delete(id)
        if (streamed) {
          // The whole text is the truth: deltas missed while the window
          // reconnected are committed from it.
          streamed.finishWith(String(e.text ?? ""), () => this.add("notice", {
            text: "(part of this answer arrived while the window was reconnecting; here it is in full)",
          }))
        } else {
          for (const block of splitMarkdown(String(e.text ?? ""), true).blocks) {
            this.add("markdown", { text: block })
          }
        }
        setTail("")
        break
      }
      case "reasoning.ended":
        // As the REPL shows it (--verbose): one dim line, whitespace run
        // together, cut in the middle to the width when it is drawn.
        if (this.options.showThinking && String(e.text ?? "").trim()) {
          this.add("reasoning", { text: stripControl(String(e.text)).replace(/\s+/g, " ").trim() })
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
        const head = toolHeadline(String(e.name ?? started.name ?? "tool"), input)
        this.lastOutput = { name: head.name, output }
        this.add("tool", { name: head.name, arg: head.arg,
                           result: toolResult(String(e.name ?? ""), output, (e.meta ?? {}) as Fields, input) })
        break
      }
      case "tool.denied":
        this.add("denied", { name: shortToolName(String(e.name ?? "")) })
        break
      case "context": {
        // What the session was started to look into (an item of the agent
        // view): one row like a tool call's, all of it with the details.
        const text = stripControl(String(e.text ?? "")).replace(/\n+$/, "")
        const lines = text.split("\n")
        const n = lines.filter((l) => l.trim()).length
        this.lastOutput = { name: "context", output: text }
        this.add("tool", { name: "context", arg: String(e.title ?? ""),
                           result: { kind: "text", ok: true, wrap: true, lines,
                                     summary: `${n} line${n === 1 ? "" : "s"}, given to the model with the first prompt` } })
        break
      }
      case "notice":
        this.add("notice", { text: e.text })
        break
      case "error":
        this.lastError = String(e.text ?? "")
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
        // A part the turn left open (it stopped mid-answer) keeps what it
        // streamed.
        for (const [id, blocks] of this.parts) {
          this.endedParts.add(id)
          blocks.finish()
        }
        this.parts.clear()
        setTail("")
        const secs = (Number(e.durationMs ?? 0) / 1000).toFixed(1)
        if (e.stopped) {
          // Ended by a cancel: neither done nor failed.
          this.add("turnEnd", { text: `◼ stopped · ${secs}s`, ok: true, stopped: true })
        } else if (e.error) {
          // The error row already said it: the line ends the turn only.
          const said = this.lastError && String(e.error).includes(firstLine(this.lastError, 400))
          this.add("turnEnd", { text: said ? `✗ ${secs}s` : `✗ ${secs}s · ${firstLine(String(e.error), 120)}`, ok: false })
        } else {
          const cost = Number(e.costUsd ?? 0)
          const parts = [`✓ ${secs}s`, `↑${compact(Number(e.inputTokens ?? 0))} ↓${compact(Number(e.outputTokens ?? 0))}`]
          if (cost > 0) parts.push(`$${cost.toFixed(4)}`)
          this.add("turnEnd", { text: parts.join(" · "), ok: true })
        }
        this.lastError = ""
        break
      }
      case "session.ended":
      case "session.removed":
        // The worker's own last word, or the service's when the worker
        // died without one: said once.
        if (!this.ended) {
          this.ended = true
          for (const blocks of this.parts.values()) blocks.finish()
          this.parts.clear()
          setTail("")
          const reason = String(e.reason ?? "closed")
          // A reason that says how the worker ended says the code already.
          const code = e.exitCode != null && e.exitCode !== 0 && !/worker/.test(reason) ? `, exit ${e.exitCode}` : ""
          this.add("ended", { text: `session ended (${reason}${code})` })
        }
        break
    }
  }

  /** The last tool result in full (the palette; ctrl+o under GNU screen),
   *  as its output (not its packaging) -- what the session forwarded. */
  expandLast(): void {
    if (!this.lastOutput) {
      this.add("notice", { text: "nothing to expand yet" })
      return
    }
    const body = stripControl(resultBody(this.lastOutput.output)).replace(/\n+$/, "")
    this.add("lines", { title: `${this.lastOutput.name}, in full`, lines: body ? body.split("\n") : ["(no output)"] })
  }

  private part(id: string): AnswerBlocks | null {
    let blocks = this.parts.get(id)
    if (!blocks) {
      // This presenter's own: a session shown again gets a new presenter,
      // whose parts start their count anew while the last one's blocks
      // may still be held.
      const key = (n: number) => `${this.no}.${id}#${n}`
      blocks = new AnswerBlocks((block, n) => this.hold(block, key(n)), (tail, n) => {
        if (!this.s) return
        this.s.tail = tail
        this.s.tailKey = key(n)
      })
      this.parts.set(id, blocks)
    }
    return blocks
  }
}

/**
 * The banner the REPL prints (cli._startup_banner), as the window draws it:
 * its wordmark box is the window's own head ("◆ DELFIN"), its first line
 * ("delfin-agent · kit/model · solo") the head's text, and every other line
 * -- workspace and git, approval, isolation, watcher, grants, the session
 * -- kept as it is: they are the lines that decide safety.
 */
export function bannerOf(info: Fields): { head: string; lines: string[] } {
  const raw = ((info.banner as string[] | undefined) ?? []).map((l) => stripControl(String(l ?? "")).trimEnd())
  // Chat's key hint ("esc interrupt · … · ctrl+d exit") is chat's: the
  // window has its own, and ctrl+d does not end the session here.
  const lines = raw.filter((l) => l.trim() && !/^[╭│╰]/.test(l) && !/^esc interrupt\b/.test(l.trim()))
  const first = lines.findIndex((l) => l.startsWith("delfin-agent"))
  let head = [info.provider && info.model ? `${info.provider}/${info.model}` : info.model, info.mode]
    .filter(Boolean).join(" · ")
  if (first >= 0) {
    head = lines[first]!.replace(/^delfin-agent\s*·\s*/, "")
    lines.splice(first, 1)
  }
  if (info.title) head += ` · ${info.title}`
  return { head, lines }
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
