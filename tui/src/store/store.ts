/**
 * The app's Solid store: events in, applied in batches.
 *
 * A turn streams hundreds of deltas a second; applying each on its own
 * would re-render the footer for every one. Events are queued and applied
 * together every ~16 ms (one frame), inside one `produce`, followed by the
 * presenter for the session in view, so scrollback and footer change in
 * the same step and in the order the events came.
 *
 * A durable event the window has already seen (same session, a `seq` not
 * above the highest one seen) is dropped on arrival: a window that attaches
 * again after a reconnect gets the session's replay, and must not show it
 * twice.
 */
import { batch } from "solid-js"
import { createStore, produce, type SetStoreFunction } from "solid-js/store"
import type { AnyEvent } from "../protocol/generated"
import { applyEvent, emptyState, type AppState } from "./state"

export type BatchListener = (events: AnyEvent[], draft: AppState) => void

export interface Clock {
  setTimeout(fn: () => void, ms: number): unknown
  clearTimeout(handle: unknown): void
  now(): number
}

const realClock: Clock = {
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (handle) => clearTimeout(handle as ReturnType<typeof setTimeout>),
  now: () => Date.now(),
}

export interface AppStore {
  state: AppState
  /** Queue one event (applied on the next flush). */
  push(event: AnyEvent): void
  /** Apply everything queued now. */
  flush(): void
  /** A change that is not an event (connection state, the current session). */
  update(fn: (draft: AppState) => void): void
  /** Run for every batch, inside the same produce, after the reducer. */
  onBatch(listener: BatchListener): () => void
}

export function createAppStore(clock: Clock = realClock, frameMs = 16): AppStore {
  const [state, setState] = createStore<AppState>(emptyState())
  return wrap(state, setState, clock, frameMs)
}

function wrap(state: AppState, setState: SetStoreFunction<AppState>,
              clock: Clock, frameMs: number): AppStore {
  let queue: AnyEvent[] = []
  let timer: unknown = null
  const seen = new Map<string, number>()
  const listeners = new Set<BatchListener>()

  const flush = () => {
    if (timer !== null) {
      clock.clearTimeout(timer)
      timer = null
    }
    if (queue.length === 0) return
    const events = queue
    queue = []
    const now = clock.now()
    batch(() => {
      setState(produce((draft) => {
        for (const event of events) applyEvent(draft, event, now)
        for (const listener of listeners) listener(events, draft)
      }))
    })
  }

  return {
    state,
    push(event) {
      const e = event as { session?: string; seq?: number }
      if (typeof e.seq === "number" && e.seq > 0 && e.session) {
        const last = seen.get(e.session) ?? 0
        if (e.seq <= last) return
        seen.set(e.session, e.seq)
      }
      queue.push(event)
      if (timer === null) timer = clock.setTimeout(flush, frameMs)
    },
    flush,
    update(fn) {
      setState(produce(fn))
    },
    onBatch(listener) {
      listeners.add(listener)
      return () => listeners.delete(listener)
    },
  }
}
