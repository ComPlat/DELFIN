/**
 * Completing what is being typed: a `/command` at the start of the prompt,
 * an `@path` anywhere in it.
 *
 * The lists are the session's (session.commands, session.completePath):
 * the commands a typed line would run, the paths under its workspace.
 * This module only finds what to complete and orders the candidates; the
 * composer asks, shows and inserts.
 */

export interface Command {
  name: string
  summary: string
  category: string
  takesArgs: boolean
  source: string
}

export interface Token {
  kind: "command" | "path"
  /** The token as typed, with its `/` or `@`. */
  text: string
  /** Where it starts and ends in the prompt (offsets). */
  start: number
  end: number
}

export interface Candidate {
  /** What replaces the token. */
  insert: string
  /** What the list shows. */
  label: string
  detail: string
  /** A command that takes no arguments runs when Enter picks it. */
  runs: boolean
  /** A directory: the completion continues inside it. */
  descends: boolean
}

/** The token the cursor is at the end of, if it is one to complete. */
export function tokenAt(text: string, cursor: number): Token | null {
  const before = text.slice(0, cursor)
  const start = Math.max(before.lastIndexOf(" "), before.lastIndexOf("\n"), before.lastIndexOf("\t")) + 1
  const word = before.slice(start)
  // The token ends where the word under the cursor ends.
  let end = cursor
  while (end < text.length && !/\s/.test(text[end]!)) end++
  if (word.startsWith("@")) return { kind: "path", text: text.slice(start, end), start, end }
  // A command only as the prompt's first word, and not a path (/home/x).
  if (start === 0 && /^\/[\w.-]*$/.test(word) && !text.slice(0, end).slice(1).includes("/")) {
    return { kind: "command", text: text.slice(0, end), start: 0, end }
  }
  return null
}

/** Commands for `typed` (with its slash): names that start with it first,
 *  then -- from two letters on, where it means something -- names that
 *  contain it, each in the session's order. */
export function matchCommands(commands: Command[], typed: string, limit = 50): Candidate[] {
  const q = typed.replace(/^\//, "").toLowerCase()
  const starts: Command[] = []
  const contains: Command[] = []
  for (const c of commands) {
    const name = c.name.slice(1).toLowerCase()
    if (name.startsWith(q)) starts.push(c)
    else if (q.length >= 2 && name.includes(q)) contains.push(c)
  }
  return [...starts, ...contains].slice(0, limit).map((c) => ({
    insert: c.takesArgs ? `${c.name} ` : c.name,
    label: c.name,
    detail: c.summary,
    runs: !c.takesArgs,
    descends: false,
  }))
}

/** Path candidates as the session listed them (relative to its workspace). */
export function pathCandidates(paths: string[]): Candidate[] {
  return paths.map((p) => ({
    insert: `@${p}`,
    label: `@${p}`,
    detail: p.endsWith("/") ? "folder" : "",
    runs: false,
    descends: p.endsWith("/"),
  }))
}

/** The longest prefix every candidate's insert shares (Tab's first step). */
export function commonPrefix(candidates: Candidate[]): string {
  if (candidates.length === 0) return ""
  let prefix = candidates[0]!.insert
  for (const c of candidates.slice(1)) {
    let i = 0
    while (i < prefix.length && i < c.insert.length && prefix[i] === c.insert[i]) i++
    prefix = prefix.slice(0, i)
  }
  return prefix
}

/** `text` with `token` replaced by `insert`; and where the cursor goes. */
export function applyCandidate(text: string, token: Token, insert: string): { text: string; cursor: number } {
  const after = text.slice(token.end)
  // A file or a command with arguments: a space after it, unless one is there.
  const spaced = !insert.endsWith("/") && !insert.endsWith(" ") && token.kind === "path" && !after.startsWith(" ")
    ? insert + " " : insert
  const out = text.slice(0, token.start) + spaced + after
  return { text: out, cursor: token.start + spaced.length }
}
