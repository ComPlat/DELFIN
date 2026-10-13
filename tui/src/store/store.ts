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
  /**
   * Hold every event pushed from now on, until `release`: while a request
   * that brings a replay (create, attach) is in flight. Its answer and the
   * live events behind it can arrive in one read; the live ones would
   * otherwise be applied first -- a durable one making the whole replay
   * look seen, a delta landing before the answer it belongs to.
   */
  hold(): void
  /** The replay first, then what was held (each hold has its release). */
  release(replay?: AnyEvent[]): void
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
  let holds = 0
  let held: AnyEvent[] = []

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

  const store: AppStore = {
    state,
    push(event) {
      if (holds > 0) {
        held.push(event)
        return
      }
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
    hold() {
      holds++
    },
    release(replay = []) {
      // In front of everything held: a replay must precede its session's
      // live events; other sessions keep their own order.
      held = [...replay, ...held]
      holds = Math.max(0, holds - 1)
      if (holds > 0) return
      const all = held
      held = []
      for (const event of all) store.push(event)
    },
  }
  return store
}

/**
 * A request whose answer brings a replay (attach, create): live events
 * that arrive with or behind the answer wait until the replay is in.
 * `register` puts the session in place and returns the replay.
 */
export async function withReplay<T>(store: AppStore, ask: () => Promise<T>,
                                    register: (answer: T) => AnyEvent[]): Promise<T> {
  store.hold()
  let replay: AnyEvent[] = []
  try {
    const answer = await ask()
    replay = register(answer)
    return answer
  } finally {
    store.release(replay)
  }
}

