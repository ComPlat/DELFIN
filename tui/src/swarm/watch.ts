/**
 * The read-only /swarm watch view (package SWV, phase 2b): a table of one
 * row per member of a swarm (or operator wave), Enter on a member opens its
 * live feed, Esc back. Read-only like attaching with tmux -r -- there is NO
 * input into a member session; the view only asks the observe.* requests
 * (swarm_rpc.py > observe.py) and renders what comes back.
 *
 * It mirrors the SW2b command view's conventions (deps injection, keymap,
 * note/draft) so the two can share the /swarm slot in app.tsx: the watch
 * opens the table on request, keeps its own selected index and, on
 * agentsOpen, swaps in the member's feed (member_feed, newest last) with
 * agentsBack returning to the table.
 */
import type { KeyEvent } from "@opentui/core"
import { createStore } from "solid-js/store"
import type { Keymap } from "../keymap/keymap"

/** One row of swarm.observe.overview (observe.py swarm_overview): package,
 *  role, name, model, state, commits on its branch beyond the base, last
 *  tool call, last commit subject, age of the last activity, pending
 *  approvals count. All strings/numbers, JSON-able. */
export interface MemberRow extends Record<string, unknown> {
  package: string
  role: string
  name: string
  model: string
  state: string
  commits_on_branch: number
  last_tool_call: Record<string, unknown>
  last_commit_subject: string
  last_activity_age_s: number
  pending_approvals_count: number
}

/** One entry of swarm.observe.feed (observe.py member_feed): a tool call
 *  (with its result status), a refusal or a commit. Newest last. */
export interface FeedEntry extends Record<string, unknown> {
  at: string
  kind: string // "tool" | "refusal" | "commit"
  text: string
}

export interface SwarmWatchDeps {
  keymap: Keymap
  now: () => number
  /** swarm.observe.overview: one row per member of the source swarm. */
  askOverview(): Promise<MemberRow[]>
  /** swarm.observe.feed: a member's recent steps, newest last. */
  askFeed(name: string): Promise<FeedEntry[]>
  note(text: string): void
  draft(): string
}

export interface SwarmWatchState {
  open: boolean
  loading: boolean
  rows: MemberRow[]
  selected: number
  /** The feed of the member currently stepped into; null = table. */
  feed: { name: string; entries: FeedEntry[]; loading: boolean } | null
  feedSelected: number
  error: string | null
  /** Every member's recent steps, cached by name so Enter can open a feed
   *  synchronously (the watch is read-only; the cache is re-prefetched on
   *  every open/refresh). */
  feeds: Record<string, FeedEntry[]>
}

/** The member row in one table line: what the footer and the selected row
 *  show for one member, read off the raw row. */
export function memberLine(row: MemberRow): string {
  const tool =
    row.last_tool_call && typeof row.last_tool_call.tool === "string"
      ? String(row.last_tool_call.tool)
      : ""
  const extra = [
    row.state,
    row.commits_on_branch ? `${row.commits_on_branch} commits` : "",
    tool ? `↳ ${tool}` : "",
  ].filter(Boolean).join(" · ")
  return `${row.name}  ${extra}`
}

/** The liveness of one member, read-only: never a process grep (agreed
 *  practice 6) -- the age comes straight from the row's
 *  ``last_activity_age_s``, already computed by observe.py from the
 *  store; the states come from observe.py's projection too. */
export function ageText(seconds: number | undefined): string {
  if (typeof seconds !== "number" || !Number.isFinite(seconds)) return ""
  const age = Math.max(0, Math.round(seconds))
  if (age < 60) return `${age}s`
  if (age < 3600) return `${Math.round(age / 60)}m`
  return `${Math.round(age / 3600)}h`
}

export function createSwarmWatch(deps: SwarmWatchDeps) {
  const [state, setState] = createStore<SwarmWatchState>({
    open: false,
    loading: false,
    rows: [],
    selected: 0,
    feed: null,
    feedSelected: 0,
    error: null,
    feeds: {},
  })

  const select = (list: unknown[], at: number, n = 1): number =>
    list.length ? ((at + n) % list.length + list.length) % list.length : 0

  const open = async (): Promise<void> => {
    setState({ open: true, loading: true, error: null })
    try {
      const rows = await deps.askOverview()
      setState({ rows: [...(Array.isArray(rows) ? rows : [])], feed: null,
                 selected: 0, loading: false })
      // Prefetch every member's feed so Enter can open a member's recent
      // steps synchronously -- the watch is read-only, so the reads happen
      // as the table opens. A member whose feed cannot be read keeps an
      // empty entry and the failure is surfaced through the note channel,
      // never by wedging the table.
      const feeds: Record<string, FeedEntry[]> = {}
      await Promise.all(
        rows.map(async (row) => {
          try {
            feeds[row.name] = await deps.askFeed(row.name)
          } catch (error) {
            feeds[row.name] = []
            deps.note(`feed of ${row.name} failed: ${String((error as Error)?.message ?? error)}`)
          }
        }),
      )
      setState({ feeds })
    } catch (error) {
      setState({
        loading: false,
        error: `swarm watch: ${String((error as Error)?.message ?? error)}`,
      })
      deps.note(`swarm watch failed: ${String((error as Error)?.message ?? error)}`)
    }
  }

  const close = (): void => {
    setState({ open: false, feed: null })
  }

  /** Enter on a member: open its live feed, read-only, straight from the
   *  cache prefetched when the table opened -- so the feed shows the
   *  member's recent steps the moment Enter is pressed. Nothing is ever
   *  sent into the session. */
  const openFeed = (): void => {
    const row = state.rows[state.selected]
    if (!row || state.feed) return
    const entries = state.feeds[row.name] ?? []
    setState({ feed: { name: row.name, entries, loading: false }, feedSelected: 0 })
    deps.note(`${row.name}: no input is accepted in the watch view`)
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    const km = deps.keymap
    const typing = deps.draft().length > 0
    if (km.matches("agentsBack", key) && !typing) {
      if (state.feed) setState({ feed: null }) // back to the table
      else close() // leaving the view
      return true
    }
    if (km.matches("agentsUp", key) || km.matches("agentsDown", key)) {
      if (typing) return true
      const up = km.matches("agentsUp", key)
      if (state.feed) {
        setState({ feedSelected: select(state.feed.entries, state.feedSelected, up ? -1 : 1) })
      } else {
        setState({ selected: select(state.rows, state.selected, up ? -1 : 1) })
      }
      return true
    }
    if (km.matches("agentsOpen", key) && !typing) {
      if (state.feed) {
        // In the feed, Enter only peeks the selected entry; it never
        // edits or answers -- the watch is read-only.
        const e = state.feed.entries[state.feedSelected]
        if (e) deps.note(`[${e.kind}] ${e.text}`)
      } else {
        openFeed()
      }
      return true
    }
    return false
  }

  const refresh = (): Promise<void> => open()

  return { state, open, close, refresh, handleKey, memberLine }
}
