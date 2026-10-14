/**
 * The session as the scrollback shows it: an ordered list of items, each
 * written once through the commit queue and kept so it can be written
 * again -- after a resize, at the new width.
 *
 * An item is recorded when it is WRITTEN, inside its queue task, not when it
 * is queued: the history is then in the order the terminal received it, and
 * a replay queued meanwhile can neither miss an item nor write one twice.
 * The replay itself is one queue task -- wipe, then every item again -- so
 * nothing can be written between the wipe and the reprint.
 */
import { CommitQueue } from "./queue"

export interface Item {
  /** Unique within the session. */
  id: string
  kind: string
  /** The turn the item belongs to (0: before the first), for a replay
   *  that is capped by turns. */
  group?: number
  [field: string]: unknown
}

export interface ReplayTarget {
  resetSplitFooterForReplay(options?: { clearSavedLines?: boolean }): void
}

/** `first`: the item's first write (not a replay's). */
export type ItemWriter = (item: Item, how: { first: boolean }) => void | Promise<void>

export class Transcript {
  readonly items: Item[] = []
  private replays = 0

  /** `keep`: what was written is kept for a replay (the scrollback's; the
   *  full screen's rows live in its own list and are never replayed). */
  constructor(
    private readonly queue: CommitQueue,
    private readonly target: ReplayTarget,
    private readonly write: ItemWriter,
    private readonly keep: boolean = true,
  ) {}

  /** Write `item` after everything queued before it, and keep it. */
  append(item: Item): Promise<void> {
    return this.queue.run(async () => {
      await this.write(item, { first: true })
      if (this.keep) this.items.push(item)
    })
  }

  /**
   * Clear the screen and the terminal's scrollback, then write the last
   * `limit` items again -- or, with `unit: "turns"`, the items of the last
   * `limit` turns (all of them when `limit` <= 0). When turns were left
   * out, `lead` gives what to write first (the banner, a note).
   */
  replay(limit = 0, unit: "items" | "turns" = "items",
         lead?: (cutTurns: number, first: Item | undefined) => Item[]): Promise<void> {
    return this.queue.run(async () => {
      this.replays++
      this.target.resetSplitFooterForReplay({ clearSavedLines: true })
      for (const item of this.pick(limit, unit, lead)) await this.write(item, { first: false })
    })
  }

  private pick(limit: number, unit: "items" | "turns",
               lead?: (cutTurns: number, first: Item | undefined) => Item[]): Item[] {
    if (limit <= 0) return this.items.slice()
    if (unit === "items") return this.items.slice(-limit)
    const last = this.items.reduce((n, item) => Math.max(n, item.group ?? 0), 0)
    const from = last - limit + 1
    if (from <= 1) return this.items.slice()
    const kept = this.items.filter((item) => (item.group ?? 0) >= from)
    return [...(lead?.(from - 1, this.items[0]) ?? []), ...kept]
  }

  get replayCount(): number {
    return this.replays
  }
}
