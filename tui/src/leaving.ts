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
  if (s?.id) {
    const name = plain(s.title || `session ${s.id.slice(0, 8)}`)
    if (!input.problem) {
      lines.push(s.state === "ended"
        ? `DELFIN: ${name} has ended${s.ended?.reason ? ` (${plain(s.ended.reason)})` : ""}.`
        : `DELFIN: ${name} keeps running in the agent service (or ← in any window).`)
    }
    // The same words chat's session list uses ("resume with: ...").
    if (ID.test(s.id)) lines.push(`resume with: ${resumeCommand(s.id)}`)
  }
  return lines.map((l) => `${l}\n`).join("")
}
