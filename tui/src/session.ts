/**
 * The scrollback machinery for the session in view, in one place: the
 * commit queue, the transcript, the presenter that turns the session's
 * events into rows, and the re-wrap after a width change. The window
 * (main.tsx) and the headless tests build it the same way.
 *
 * Resizing. After a width change settles, the screen and the terminal's
 * scrollback are cleared and the last `replayTurns` turns are written
 * again at the new width (the user's choice, docs/TUI_REWRITE.md). Under
 * GNU screen, which cannot clear its scrollback (no ED3), a replay would
 * duplicate every row: a divider says the width changed instead, and the
 * rows above keep their old width.
 */
import type { CliRenderer } from "@opentui/core"
import { CommitQueue } from "./scrollback/queue"
import { Presenter } from "./scrollback/present"
import { replayOnWidthChange } from "./scrollback/resize"
import { itemWriter } from "./scrollback/rows"
import { Transcript, type Item } from "./scrollback/transcript"
import type { AppStore } from "./store/store"
import type { Theme } from "./theme/theme"

export interface SessionViewOptions {
  showThinking: boolean
  /** Turns written again after a resize or a redraw (0: all). */
  replayTurns: number
  /** Print a divider instead of replaying after a resize (GNU screen).
   *  Defaults to whether $STY says GNU screen. */
  dividerOnResize?: boolean
}

export interface SessionView {
  queue: CommitQueue
  transcript: Transcript
  /** Prompts sent with pastes put back -> as typed (see App's `sentAs`). */
  sentAs: Map<string, string>
  /** Put session `id` in view: its events become rows from now on. With
   *  `divider`, a divider row says which session follows. */
  show(id: string, options?: { divider?: string }): void
  /** A note of the window's own (a config problem, a lost session), as a
   *  row in the scrollback, after what is already queued. */
  note(text: string): void
  /** ctrl+o: the last tool result of the session in view, in full. */
  expandLast(): void
  /** Thinking shown or hidden from now on. */
  setShowThinking(show: boolean): void
  readonly showThinking: boolean
  redraw(): Promise<void>
  stop(): void
}

export function startSessionView(renderer: CliRenderer, store: AppStore, theme: () => Theme,
                                 options: SessionViewOptions): SessionView {
  const target = {
    requestLive: () => renderer.requestLive(),
    dropLive: () => renderer.dropLive(),
    resetSplitFooterForReplay: (o?: { clearSavedLines?: boolean }) => renderer.resetSplitFooterForReplay(o),
  }
  // Shared by every presenter of this window: turns are numbered across
  // the sessions shown in it, so a replay by turns keeps the newest.
  const groups = { n: 0 }
  let notes = 0
  const own = (kind: string, fields: Record<string, unknown>): Item =>
    ({ id: `window:${++notes}`, kind, group: groups.n, ...fields })
  const note = (text: string) => {
    void transcript.append(own("notice", { text }))
  }
  // A row that cannot be drawn is said once: if saying it fails too (the
  // same fault), that is not said again -- or the queue would never empty.
  let reporting = false
  const queue = new CommitQueue(target, (error) => {
    if (reporting) return
    reporting = true
    void transcript.append(own("notice", { text: `could not draw a row: ${String(error)}` }))
      .finally(() => { reporting = false })
  })
  const transcript = new Transcript(queue, target, itemWriter(renderer, theme))
  const presentOptions = { showThinking: options.showThinking, shownAs: (text: string) => {
    const shown = sentAs.get(text)
    if (shown !== undefined) sentAs.delete(text)
    return shown
  } }

  /** When turns were left out of a replay: the banner of the session they
   *  belonged to, and a note saying how many. */
  const lead = (cut: number): Item[] => {
    const kept = transcript.items.filter((i) => (i.group ?? 0) >= cut + 1)
    const firstKept = transcript.items.indexOf(kept[0]!)
    const banner = [...transcript.items.slice(0, Math.max(0, firstKept))].reverse().find((i) => i.kind === "banner")
    const say = own("notice", { text: `(${cut} earlier turn${cut === 1 ? "" : "s"} not shown again after the resize)` })
    return banner ? [banner, say] : [say]
  }
  const replay = () => transcript.replay(options.replayTurns, "turns", lead)
  const divider = options.dividerOnResize ?? Boolean(process.env.STY)
  const stopResize = replayOnWidthChange(renderer, transcript, {
    onSettled: (width) => {
      if (divider) {
        void transcript.append(own("divider", { text: `${width} columns now · the rows above keep their width` }))
      } else void replay()
    },
  })
  let off: (() => void) | null = null
  let presenter: Presenter | null = null
  const sentAs = new Map<string, string>()
  return {
    queue,
    transcript,
    sentAs,
    show(id, how = {}) {
      off?.()
      if (how.divider) {
        groups.n++
        void transcript.append(own("divider", { text: how.divider }))
      }
      const shown = new Presenter(id, transcript, presentOptions, groups)
      presenter = shown
      off = store.onBatch((events, draft) => shown.handle(events, draft))
      store.update((d) => { d.current = id })
    },
    note,
    expandLast: () => presenter?.expandLast(),
    setShowThinking(show) {
      presentOptions.showThinking = show
    },
    get showThinking() {
      return presentOptions.showThinking
    },
    redraw: replay,
    stop() {
      off?.()
      stopResize()
    },
  }
}
