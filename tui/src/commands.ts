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
  { name: "/merge", summary: "merge a session worktree back into this checkout (conflicts shown, nothing overwritten)",
    category: "session", takesArgs: true, source: "window" },
  { name: "/swarm", summary: "manage swarms: packages, members, what waits on you (SW2b)",
    category: "session", takesArgs: false, source: "window" },
]

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

// -- unknown /commands ----------------------------------------------------------

/** Commands the window reads itself that take no row of their own in
 *  WINDOW_COMMANDS: /name (nameOf), /keys and /editor, and the ones whose
 *  picker it opens (/model, /effort, /mode; the session knows them too). */
export const WINDOW_ONLY_NAMES = ["/name", "/keys", "/editor", "/model", "/effort", "/mode"]

/** What a command name may look like: a letter, then letters, digits,
 *  "_", "-", "." or ":" (a plugin's "/plugin:skill"). */
const COMMAND_SHAPE = /^\/[A-Za-z][\w.:-]*$/

/** Edit distance (insert, delete, substitute), case-insensitive. */
function distance(a: string, b: string): number {
  const s = a.toLowerCase()
  const t = b.toLowerCase()
  let prev = Array.from({ length: t.length + 1 }, (_, j) => j)
  for (let i = 1; i <= s.length; i++) {
    const row = [i]
    for (let j = 1; j <= t.length; j++) {
      row[j] = Math.min(prev[j]! + 1, row[j - 1]! + 1, prev[j - 1]! + (s[i - 1] === t[j - 1] ? 0 : 1))
    }
    prev = row
  }
  return prev[t.length]!
}

/** The known command `head` was probably meant to be, or "": the closest
 *  within two edits (one for names of up to three letters), or the only
 *  one it begins (a name cut short, "/hel"). */
export function nearestCommand(head: string, known: readonly string[]): string {
  const name = head.toLowerCase()
  let best = ""
  let bestD = Infinity
  for (const k of known) {
    const d = distance(name, k)
    if (d < bestD) {
      best = k
      bestD = d
    }
  }
  const limit = name.length <= 4 ? 1 : 2
  if (best && bestD <= limit) return best
  const started = known.filter((k) => name.length >= 3 && k.toLowerCase().startsWith(name))
  return started.length === 1 ? started[0]! : ""
}

/** A line refused as an unknown command: its first word, and the command
 *  it was probably meant to be ("" when none is close). */
export interface UnknownCommand {
  head: string
  near: string
}

/**
 * Whether a line is an unknown `/command`, which is answered in the window
 * and never sent (Claude Code, the dashboard's "Unknown command", the REPL's
 * "did you mean"). `known` are the names the window and the session know
 * (WINDOW_COMMANDS, WINDOW_ONLY_NAMES, session.commands); `exists` says
 * whether a path is there.
 *
 * Not a command, so sent as typed:
 *   - a line that does not start with "/", or only "/";
 *   - a path: "/home/me/x.py is broken" -- a first word with another "/",
 *     unless what comes before that "/" is a known command ("/help/keys",
 *     two commands run together), and a first word that names something on
 *     disk ("/tmp", "/scratch");
 *   - a first word no command could be called ("/.hidden", "/* ... *\/",
 *     "/123").
 */
export function unknownCommandOf(text: string, known: readonly string[],
                                 exists: (path: string) => boolean = () => false): UnknownCommand | null {
  const line = text.trim()
  if (!line.startsWith("/") || line.length < 2) return null
  const first = line.split(/\s+/)[0]!
  const names = new Set(known.map((k) => k.toLowerCase()))
  const second = first.indexOf("/", 1)
  if (second > 0) {
    const lead = first.slice(0, second)
    if (COMMAND_SHAPE.test(lead) && names.has(lead.toLowerCase())) {
      return { head: first, near: known.find((k) => k.toLowerCase() === lead.toLowerCase()) ?? lead }
    }
    return null                                  // a path
  }
  if (!COMMAND_SHAPE.test(first)) return null
  if (names.has(first.toLowerCase())) return null
  let there = false
  try {
    there = exists(first)
  } catch {
    there = false
  }
  if (there) return null                         // "/tmp": a path after all
  return { head: first, near: nearestCommand(first, known) }
}

/** What the window says about an unknown command. */
export function unknownCommandNote(unknown: UnknownCommand): string {
  const near = unknown.near ? ` — did you mean ${unknown.near}?` : ""
  return `Unknown command: ${unknown.head}${near} /help lists them. Nothing was sent.`
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
