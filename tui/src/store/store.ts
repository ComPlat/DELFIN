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
   * Hold the events `test` picks (every one by default) from now on, until
   * this hold's `release`: while a request that brings a replay (create,
   * attach) is in flight. Its answer and the live events behind it can
   * arrive in one read; the live ones would otherwise be applied first --
   * a durable one making the whole replay look seen, a delta landing
   * before the answer it belongs to. Only the session the replay is about
   * is held: a create may take minutes, and the other sessions go on.
   */
  hold(test?: (event: AnyEvent) => boolean): Hold
  /** Its replay first, then what it held (other holds keep theirs). */
  release(hold: Hold, replay?: AnyEvent[]): void
  /** Forget which of a session's events were seen: its replay is shown
   *  again (switching back to a session shown before). */
  unsee(session: string): void
}

/** One hold: what it holds. */
export interface Hold {
  test: (event: AnyEvent) => boolean
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
  let holds: Hold[] = []
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
      if (holds.some((h) => h.test(event))) {
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
    hold(test = () => true) {
      const hold: Hold = { test }
      holds.push(hold)
      return hold
    },
    unsee(session) {
      seen.delete(session)
    },
    release(hold, replay = []) {
      holds = holds.filter((h) => h !== hold)
      // In front of what it held: a replay must precede its session's live
      // events. What another hold still wants stays held, in order.
      const keep: AnyEvent[] = []
      const go: AnyEvent[] = []
      for (const event of held) (holds.some((h) => h.test(event)) ? keep : go).push(event)
      held = keep
      for (const event of [...replay, ...go]) store.push(event)
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
                                    register: (answer: T) => AnyEvent[],
                                    test?: (event: AnyEvent) => boolean): Promise<T> {
  const hold = store.hold(test)
  let replay: AnyEvent[] = []
  try {
    const answer = await ask()
    replay = register(answer)
    return answer
  } finally {
    store.release(hold, replay)
  }
}

/** What a hold for one session picks: its own events. */
export function eventsOf(session: string): (event: AnyEvent) => boolean {
  return (event) => (event as { session?: string }).session === session
}

/** What a hold for a session still to be created picks: the worker
 *  events of a session this window does not know yet -- the new one's,
 *  which can come in the same read as the answer. The service's own
 *  (summaries, removals) pass: they register a session, and another
 *  window's new session must not wait for this one's start. */
export function eventsOfUnknown(store: AppStore): (event: AnyEvent) => boolean {
  return (event) => {
    const e = event as { session?: string; type?: string }
    return Boolean(e.session) && !store.state.sessions[e.session!]
      && e.type !== "session.summary" && e.type !== "session.removed"
  }
}

/**
 * What the service lists (a hello, sessions.list) is what this window
 * knows: each listed session's summary applied, and a session it knew
 * that is no longer listed ended -- after a service restart its rows
 * would otherwise stay, waiting and working for ever.
 */
export function listed(store: AppStore, summaries: { sessionId?: string }[], now: number = Date.now()): void {
  const ids = new Set<string>()
  for (const summary of summaries) {
    if (!summary?.sessionId) continue
    ids.add(summary.sessionId)
    store.push({ type: "session.summary", session: summary.sessionId, ts: now / 1000, summary } as never)
  }
  for (const [id, s] of Object.entries(store.state.sessions)) {
    if (!ids.has(id) && s.state !== "ended") {
      store.push({ type: "session.removed", session: id, ts: now / 1000,
                   reason: "not in the agent service any more" } as never)
    }
  }
}

