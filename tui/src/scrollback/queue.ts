/**
 * The scrollback commit queue. Every write to the terminal's scrollback goes
 * through here, one at a time, in the order it was queued.
 *
 * Three properties, each one a lesson from phase 0 (docs/TUI_REWRITE.md):
 *
 * - Order. A task starts only when the previous one has finished, so a tool
 *   row can never land in the middle of an answer that is still being
 *   committed, and a replay is never interleaved with anything.
 * - Liveness. OpenTUI flushes queued scrollback only when frames run, and
 *   demand-driven rendering stalls on a burst of commits. The renderer is
 *   held live while anything is queued and let go when the queue empties,
 *   so an idle TUI costs no CPU.
 * - It never stalls. A task that throws is reported and the next one runs;
 *   waits inside tasks are bounded with `bounded` (a highlighter that never
 *   settles must not hold the transcript forever).
 */

export interface LiveTarget {
  requestLive(): void
  dropLive(): void
}

export class CommitQueue {
  private tail: Promise<void> = Promise.resolve()
  private pending = 0

  constructor(
    private readonly target: LiveTarget,
    private readonly onError: (error: unknown) => void = () => {},
  ) {}

  /** Tasks queued or running. */
  get size(): number {
    return this.pending
  }

  /** Run `task` after everything queued before it. Never rejects. */
  run(task: () => void | Promise<void>): Promise<void> {
    if (this.pending++ === 0) this.target.requestLive()
    const done = this.tail
      .then(task)
      .catch((error) => {
        try {
          this.onError(error)
        } catch {
          // An error reporter that throws must not stall the queue either.
        }
      })
      .finally(() => {
        if (--this.pending === 0) this.target.dropLive()
      })
    this.tail = done
    return done
  }

  /** Resolves once everything queued so far has been written. */
  drain(): Promise<void> {
    return this.run(() => {})
  }
}

/** `promise`, or `undefined` once `ms` have passed -- whichever comes first. */
export function bounded<T>(promise: Promise<T>, ms: number): Promise<T | undefined> {
  let timer: ReturnType<typeof setTimeout> | undefined
  const timeout = new Promise<undefined>((resolve) => {
    timer = setTimeout(() => resolve(undefined), ms)
  })
  return Promise.race([promise, timeout]).finally(() => clearTimeout(timer))
}
