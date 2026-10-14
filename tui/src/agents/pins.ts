/**
 * Pinned sessions (the agent view's ctrl+t), kept across windows in the
 * window's own state file, `~/.delfin/tui_state.json` (DELFIN_SCRATCH_STATE
 * moves it, as it moves every ~/.delfin sink). A file that cannot be read
 * costs the pins, never the start.
 */
import { existsSync, mkdirSync, readFileSync, renameSync, writeFileSync } from "node:fs"
import { homedir } from "node:os"
import { dirname, join } from "node:path"

export function statePath(env: Record<string, string | undefined> = process.env): string {
  const home = env.DELFIN_SCRATCH_STATE || env.HOME || homedir()
  return env.DELFIN_TUI_STATE || join(home, ".delfin", "tui_state.json")
}

export class Pins {
  private readonly set: Set<string>

  constructor(private readonly path: string | null = null, initial: string[] = []) {
    this.set = new Set(initial)
  }

  static load(path = statePath()): Pins {
    try {
      if (existsSync(path)) {
        const data = JSON.parse(readFileSync(path, "utf8"))
        if (Array.isArray(data?.pinned)) return new Pins(path, data.pinned.filter((p: unknown) => typeof p === "string"))
      }
    } catch {
      // unreadable: no pins
    }
    return new Pins(path)
  }

  has(id: string): boolean {
    return this.set.has(id)
  }

  get all(): ReadonlySet<string> {
    return this.set
  }

  /** Pin or unpin `id` -- as this window shows it; true when pinned now. */
  toggle(id: string): boolean {
    const on = !this.set.has(id)
    if (this.path) {
      // What another window pinned or unpinned since is kept: the file is
      // read again, and only this one pin changes in it.
      const fresh = Pins.load(this.path)
      this.set.clear()
      for (const pin of fresh.set) this.set.add(pin)
    }
    this.set.delete(id)
    if (on) this.set.add(id)
    // Pins of sessions long gone are let go, oldest first.
    for (const pin of [...this.set].slice(0, Math.max(0, this.set.size - PINS_MAX))) this.set.delete(pin)
    this.save()
    return on
  }

  private save(): void {
    if (!this.path) return
    try {
      mkdirSync(dirname(this.path), { recursive: true, mode: 0o700 })
      const tmp = `${this.path}.${process.pid}.tmp`
      writeFileSync(tmp, JSON.stringify({ pinned: [...this.set] }, null, 2) + "\n", { mode: 0o600 })
      renameSync(tmp, this.path)
    } catch {
      // a read-only home keeps the pins for this window only
    }
  }
}

/** Pins kept at most. */
export const PINS_MAX = 100
