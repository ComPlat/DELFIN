/**
 * The prompts typed before, for ↑/↓ and ctrl+r.
 *
 * The window keeps its own file, `~/.delfin/tui_history.jsonl` (one JSON
 * object per line, so a prompt may span lines and keep its paste
 * placeholders), and reads the REPL's readline history
 * (`~/.delfin/agent_repl_history`) in front of it, read-only: what was
 * typed into `delfin-agent chat` comes back here too. Appends only, one
 * line per prompt, so two windows never overwrite each other's prompts.
 */
import { appendFileSync, chmodSync, existsSync, mkdirSync, readFileSync, renameSync, writeFileSync } from "node:fs"
import { homedir } from "node:os"
import { dirname, join } from "node:path"

export interface Entry {
  /** The prompt as typed, paste placeholders included. */
  text: string
  /** What each placeholder stands for (paste.ts). */
  pastes?: Record<string, string>
}

/** Entries kept in memory and, after a trim, in the file. */
export const HISTORY_MAX = 1000

export function historyPaths(env: Record<string, string | undefined> = process.env) {
  const home = env.HOME || homedir()
  // The window's own state follows DELFIN_SCRATCH_STATE, as every other
  // ~/.delfin sink does; the REPL's history is always the real one.
  const state = env.DELFIN_SCRATCH_STATE || home
  return {
    own: env.DELFIN_TUI_HISTORY || join(state, ".delfin", "tui_history.jsonl"),
    repl: join(home, ".delfin", "agent_repl_history"),
  }
}

/** readline's file: one entry per line. libedit (macOS) writes a header
 *  and octal escapes; GNU readline may write `#<time>` lines. */
export function parseReadline(text: string): string[] {
  const lines = text.split("\n")
  const libedit = lines[0] === "_HiStOrY_V2_"
  const out: string[] = []
  for (const raw of libedit ? lines.slice(1) : lines) {
    if (!raw.trim() || /^#\d+$/.test(raw)) continue
    out.push(libedit ? raw.replace(/\\(\d{3})/g, (_, oct) => String.fromCharCode(parseInt(oct, 8))) : raw)
  }
  return out
}

export function parseOwn(text: string): Entry[] {
  const out: Entry[] = []
  for (const line of text.split("\n")) {
    if (!line.trim()) continue
    try {
      const e = JSON.parse(line)
      if (e && typeof e.text === "string" && e.text) {
        const entry: Entry = { text: e.text }
        if (e.pastes && typeof e.pastes === "object") entry.pastes = e.pastes
        out.push(entry)
      }
    } catch {
      // A torn last line (a window killed mid-write) costs that line only.
    }
  }
  return out
}

function samePastes(a: Record<string, string> | undefined, b: Record<string, string> | undefined): boolean {
  const ka = Object.keys(a ?? {})
  const kb = Object.keys(b ?? {})
  return ka.length === kb.length && ka.every((k) => (a ?? {})[k] === (b ?? {})[k])
}

export class History {
  readonly entries: Entry[] = []

  constructor(private readonly ownPath: string | null = null, entries: Entry[] = []) {
    for (const e of entries) this.push(e)
  }

  static load(paths = historyPaths()): History {
    const entries: Entry[] = []
    try {
      if (existsSync(paths.repl)) for (const text of parseReadline(readFileSync(paths.repl, "utf8"))) entries.push({ text })
    } catch {
      // Unreadable: the window's own history still works.
    }
    let own: Entry[] = []
    try {
      if (existsSync(paths.own)) own = parseOwn(readFileSync(paths.own, "utf8"))
    } catch {
      own = []
    }
    const history = new History(paths.own, [...entries, ...own])
    if (own.length > 2 * HISTORY_MAX) history.trimFile(own.slice(-HISTORY_MAX))
    return history
  }

  private push(entry: Entry): boolean {
    const last = this.entries[this.entries.length - 1]
    // The same prompt twice in a row is kept once -- unless what its
    // placeholders stand for differs (two different logs pasted).
    if (last && last.text === entry.text && samePastes(last.pastes, entry.pastes)) return false
    this.entries.push(entry)
    if (this.entries.length > HISTORY_MAX) this.entries.splice(0, this.entries.length - HISTORY_MAX)
    return true
  }

  /** Remember a prompt (in memory and in the window's file). */
  add(entry: Entry): void {
    if (!entry.text.trim() || !this.push(entry)) return
    if (!this.ownPath) return
    try {
      mkdirSync(dirname(this.ownPath), { recursive: true, mode: 0o700 })
      const line = JSON.stringify(entry.pastes && Object.keys(entry.pastes).length
        ? { text: entry.text, pastes: entry.pastes, ts: Date.now() / 1000 }
        : { text: entry.text, ts: Date.now() / 1000 }) + "\n"
      // Prompts describe the user's project: the file is theirs alone.
      appendFileSync(this.ownPath, line, { mode: 0o600 })
    } catch {
      // A full disk or a read-only home must not cost the prompt.
    }
  }

  private trimFile(keep: Entry[]): void {
    if (!this.ownPath) return
    try {
      const tmp = `${this.ownPath}.${process.pid}.tmp`
      writeFileSync(tmp, keep.map((e) => JSON.stringify(e)).join("\n") + "\n", { mode: 0o600 })
      chmodSync(tmp, 0o600)
      renameSync(tmp, this.ownPath)
    } catch {
      // Left long; the next start tries again.
    }
  }

  /** The newest entry before `before` (exclusive) whose text contains
   *  `query`, case-insensitively; -1 when there is none. */
  search(query: string, before: number = this.entries.length): number {
    const needle = query.toLowerCase()
    for (let i = Math.min(before, this.entries.length) - 1; i >= 0; i--) {
      if (this.entries[i]!.text.toLowerCase().includes(needle)) return i
    }
    return -1
  }
}

/**
 * Walking the history with ↑/↓, the way readline does: the draft being
 * typed is kept when the walk starts and comes back when it walks past
 * the newest entry.
 */
export class HistoryWalk {
  private index: number
  private draft: Entry | null = null

  constructor(private readonly history: History) {
    this.index = history.entries.length
  }

  get walking(): boolean {
    return this.draft !== null
  }

  /** The older entry, or null at the oldest. */
  older(current: Entry): Entry | null {
    if (this.index <= 0) return null
    if (this.draft === null) {
      this.draft = current
      this.index = this.history.entries.length
    }
    this.index--
    return this.history.entries[this.index] ?? null
  }

  /** The newer entry, the draft after the newest, or null when not walking. */
  newer(): Entry | null {
    if (this.draft === null) return null
    this.index++
    if (this.index >= this.history.entries.length) {
      const draft = this.draft
      this.reset()
      return draft
    }
    return this.history.entries[this.index] ?? null
  }

  /** Stop walking (the user typed, or sent the prompt). */
  reset(): void {
    this.draft = null
    this.index = this.history.entries.length
  }
}
