/**
 * Which sessions of this window were never asked anything: a new session
 * (not resumed, not --continue) with no message and no turn yet. Such a
 * session is ended when the window lets it go -- leaving the window,
 * opening another session in it -- instead of running on in the agent
 * service with a worker of its own, and its worker does not save it
 * (worker.py): Claude Code keeps no empty conversation either.
 *
 * Input: the events the window hears (user.message, turn.started) and the
 * session's summary (firstPrompt, lastInput -- what another window sent
 * counts too). Output: isEmpty(id).
 */

/** What of a summary says the session was asked something. */
export interface AskedSummary {
  firstPrompt?: string | null
  lastInput?: number | null
}

export class EmptySessions {
  private fresh = new Set<string>()
  private used = new Set<string>()
  private gone = new Set<string>()

  /** This window started `id` as a new session. */
  started(id: string): void {
    if (id && !this.used.has(id)) this.fresh.add(id)
  }

  /** An event from the service: a message or a turn marks its session. */
  saw(event: { type?: unknown; session?: unknown } | null | undefined): void {
    const type = String(event?.type ?? "")
    if (type !== "user.message" && type !== "turn.started") return
    const id = String(event?.session ?? "")
    if (!id) return
    this.used.add(id)
    this.fresh.delete(id)
  }

  /** Started here as new, and nothing asked of it: no message, no turn,
   *  nothing its summary names. */
  isEmpty(id: string, summary?: AskedSummary | null): boolean {
    if (!id || !this.fresh.has(id) || this.used.has(id)) return false
    if (summary && (String(summary.firstPrompt ?? "").trim() || Number(summary.lastInput ?? 0) > 0)) return false
    return true
  }

  /** `id` was ended for being empty. */
  ended(id: string): void {
    this.fresh.delete(id)
    this.gone.add(id)
  }

  /** Whether `id` was ended for being empty (the leaving note then names
   *  no way back to it). */
  wasEnded(id: string): boolean {
    return this.gone.has(id)
  }
}
