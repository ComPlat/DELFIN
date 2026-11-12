/**
 * The window's own `/commands`: what the session cannot do for it.
 *
 * `/rename [name]` names the session -- the prompt box, the agent view and
 * the terminal's title show the name, and the session store saves it
 * (the session keeps the address it started with: docs/TUI_PROTOCOL.md).
 * Without a name it is named after what it was first asked. `/name [name]`
 * does the same, and without a name puts the default back (folder name and
 * short id, as the service's summary names it).
 *
 * `/resume [id]` opens another session in this window: bare, the picker;
 * with an id (the argument menu offers the sessions), that session.
 *
 * `/clusters` and a bare `/jobs` open a screen of the window's own
 * (clusters/view.tsx): DELFIN's clusters, where one is logged into, and
 * what runs in the background. `/jobs` with an argument (`/jobs watch on`)
 * is the session's, as typed.
 */

/** What makes a line a shell command: repl_commands.SHELL_PREFIX (a test
 *  holds the two equal). */
export const SHELL_PREFIX = "!"

/** Where the line after a leading "!" begins (blanks before the "!" are
 *  stripped from a sent line, as the session strips them), else -1. */
export function shellLead(text: string): number {
  const m = /^\s*!/.exec(text)
  return m ? m[0].length : -1
}

/** As session.commands lists the session's: what completion offers. */
export const WINDOW_COMMANDS = [
  // /name is read too (nameOf) but takes no row of its own: the menu shows
  // seven rows, and the session's commands come after these.
  { name: "/rename", summary: "name this session (no name: after its first prompt; /name: back to the default)",
    category: "session", takesArgs: true, source: "window" },
  { name: "/resume", summary: "open another session in this window (running here or saved)",
    category: "session", takesArgs: true, source: "window" },
  { name: "/clusters", summary: "DELFIN's clusters: logged in, logged out or refused; log in to one",
    category: "workspace", takesArgs: false, source: "window" },
  // Its arguments are the session's /jobs's (the argument menu asks it).
  { name: "/jobs", summary: "what runs in the background: the agents', the job monitor's, each cluster's",
    category: "workspace", takesArgs: true, source: "window" },
  { name: "/merge", summary: "merge a session worktree back into this checkout (conflicts shown, nothing overwritten)",
    category: "session", takesArgs: true, source: "window" },
  { name: "/swarm", summary: "manage swarms: packages, members, what waits on you (SW2b)",
    category: "session", takesArgs: false, source: "window" },
]

/** The screen a bare `/clusters` or `/jobs` line opens, or null. */
export function screenOf(text: string): "clusters" | "jobs" | null {
  const m = /^\/(clusters|jobs)$/.exec(text.trim())
  return m ? m[1] as "clusters" | "jobs" : null
}

/** The session a `/resume <id>` line asks for, or null for another line
 *  (a bare /resume opens the picker). */
export function resumeOf(text: string): string | null {
  const m = /^\/resume\s+(\S+)\s*$/.exec(text.trim())
  return m ? m[1]! : null
}

/** The name a `/rename` line asks for ("" for "after its first prompt"),
 *  or null when the line is not a `/rename`. */
export function renameOf(text: string): string | null {
  const m = /^\/rename(?:\s+([\s\S]*))?$/.exec(text.trim())
  return m ? (m[1] ?? "").replace(/\s+/g, " ").trim() : null
}

/** The name a `/name` line gives ("" for the default: the service's
 *  folder name and short id), or null when the line is not a `/name`.
 *  The same name every window shows: the service keeps it. */
export function nameOf(text: string): string | null {
  const m = /^\/name(?:\s+([\s\S]*))?$/.exec(text.trim())
  return m ? (m[1] ?? "").replace(/\s+/g, " ").trim() : null
}

/** A `/merge [path]` line: the path it asks to merge back, "" for the
 *  picker, null when the line is not a `/merge`. */
export function mergeOf(text: string): string | null {
  const m = /^\/merge(?:\s+(\S[\s\S]*))?$/.exec(text.trim())
  return m ? (m[1] ?? "").trim() : null
}

/** NEEDS_LOCAL_PROJECT (-32006): session.create refused a cwd without a
 *  repository -- the user may consent to making it a local project. */
export const NEEDS_LOCAL_PROJECT = -32006

/** Whether a refused session.create is that plain-folder case: the retry
 *  params the consent asks for, else null (another refusal, no error). */
export function worktreeRetryOf(error:
        { code?: number; message?: string } | null | undefined):
        { initProject: true } | null {
  return error && error.code === NEEDS_LOCAL_PROJECT ? { initProject: true } : null
}

/** A name from what a session was first asked: its opening words -- at
 *  most five, about 40 characters -- without a closing stop or colon. */
export function nameFromPrompt(prompt: string): string {
  const words = String(prompt ?? "").replace(/\s+/g, " ").trim().split(" ").filter(Boolean)
  let name = ""
  for (const word of words.slice(0, 5)) {
    const next = name ? `${name} ${word}` : word
    if (next.length > 40 && name) break
    name = next
  }
  return [...name].slice(0, 40).join("").replace(/[\s.,;:!?…—–-]+$/u, "")
}
