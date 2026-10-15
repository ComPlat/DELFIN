/**
 * Pinned sessions (the agent view's ctrl+t), kept across windows in the
 * window's own state file (config/state.ts). A file that cannot be read
 * costs the pins, never the start.
 */
import { readState, statePath, updateState } from "../config/state"

export { statePath }

export class Pins {
  private readonly set: Set<string>

  constructor(private readonly path: string | null = null, initial: string[] = []) {
    this.set = new Set(initial)
  }

  static load(path = statePath()): Pins {
    const pinned = readState(path).pinned
    return new Pins(path, Array.isArray(pinned) ? pinned.filter((p): p is string => typeof p === "string") : [])
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
    if (this.path) updateState({ pinned: [...this.set] }, this.path)
  }
}

/** Pins kept at most. */
export const PINS_MAX = 100
