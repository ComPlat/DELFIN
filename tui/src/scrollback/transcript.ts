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
  [field: string]: unknown
}

export interface ReplayTarget {
  resetSplitFooterForReplay(options?: { clearSavedLines?: boolean }): void
}

export type ItemWriter = (item: Item) => void | Promise<void>

export class Transcript {
  readonly items: Item[] = []
  private replays = 0

  constructor(
    private readonly queue: CommitQueue,
    private readonly target: ReplayTarget,
    private readonly write: ItemWriter,
  ) {}

  /** Write `item` after everything queued before it, and keep it. */
  append(item: Item): Promise<void> {
    return this.queue.run(async () => {
      await this.write(item)
      this.items.push(item)
    })
  }

  /**
   * Clear the screen and the terminal's scrollback, then write the last
   * `limit` items again (all of them when `limit` <= 0).
   */
  replay(limit = 0): Promise<void> {
    return this.queue.run(async () => {
      this.replays++
      this.target.resetSplitFooterForReplay({ clearSavedLines: true })
      const items = limit > 0 ? this.items.slice(-limit) : this.items.slice()
      for (const item of items) await this.write(item)
    })
  }

  get replayCount(): number {
    return this.replays
  }
}
