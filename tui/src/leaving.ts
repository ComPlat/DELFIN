/**
 * What the window prints when it has given the terminal back: what became
 * of the session in view, and the command that opens exactly that session
 * again. Written after the renderer restored the terminal (main.tsx), so it
 * stays in the shell's scrollback in both screen modes.
 *
 * Input: the session in view (its id, title, state, why it ended) and the
 * connection's problem, if the service went away. Output: "" or whole
 * lines, each ending in "\n", free of control characters (a title is the
 * model's or the user's text and must not move the cursor).
 */

export interface LeavingSession {
  id: string
  title?: string | null
  state?: string
  ended?: { reason?: string | null } | null
  /** Ended on leaving because nothing was asked of it (empty.ts): not
   *  kept, so no way back to it is printed. */
  discarded?: boolean
  /** Keep session (the service's summary `keep`): false, the service ends
   *  the session once this window, which holds it, has gone. */
  kept?: boolean
}

/** Whether a session runs on after its windows leave: the service's Keep
 *  session flag (`session.setKeep`, `/keep`). A service before the switch
 *  kept every session and says nothing. */
export function keptOf(summary?: { keep?: boolean | null } | null): boolean {
  return summary?.keep !== false
}

/** The session.create params of a session this window starts: Keep
 *  session off, held by this window -- leaving it (ctrl+c twice, ctrl+d, a
 *  closed terminal) ends the session, saved; `/keep on` keeps it running in
 *  the agent service. A service before the switch ignores the field and
 *  keeps every session, and keptOf reads that from its summaries. */
export function startedHere<T extends Record<string, unknown>>(params: T): T & { keep: false } {
  return { ...params, keep: false }
}

/** What leaving the window does to the session in view, as ctrl+c and the
 *  menu say it. */
export function leaveWords(kept: boolean): string {
  return kept ? "the session keeps running in the agent service"
    : "the session ends -- /keep on keeps it"
}

/** Characters a session id is made of; anything else is not printed as a
 *  command to paste into a shell. */
const ID = /^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$/

/** The shell command that comes back to session `id` (`tui --resume`, the
 *  launcher's flag; cli.py `tui --session/--resume`). */
export function resumeCommand(id: string): string {
  return `delfin-agent tui --resume ${id}`
}

const plain = (text: string) => text.replace(/[\u0000-\u001f\u007f-\u009f]/g, " ")

export function leavingNote(input: { session?: LeavingSession; problem?: string }): string {
  const lines: string[] = []
  if (input.problem) lines.push(`DELFIN: ${plain(input.problem)}`)
  const s = input.session
  if (s?.id && s.discarded) {
    // Nothing was asked of it: like Claude Code, nothing is kept.
    if (!input.problem) lines.push("DELFIN: nothing was asked in this session, so it was not kept.")
  } else if (s?.id) {
    const name = plain(s.title || `session ${s.id.slice(0, 8)}`)
    if (!input.problem) {
      lines.push(s.state === "ended"
        ? `DELFIN: ${name} has ended${s.ended?.reason ? ` (${plain(s.ended.reason)})` : ""}.`
        : s.kept === false
          ? `DELFIN: ${name} ends with this window (Keep session is off). It is saved.`
          : `DELFIN: ${name} keeps running in the agent service (or ← in any window).`)
    }
    // The same words chat's session list uses ("resume with: ...").
    if (ID.test(s.id)) lines.push(`resume with: ${resumeCommand(s.id)}`)
  }
  return lines.map((l) => `${l}\n`).join("")
}
