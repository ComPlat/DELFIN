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
 *
 * Full screen (`screen` given): the rows go to the window's own list
 * (scrollback/screenview.tsx) instead of the terminal's scrollback; they
 * re-wrap by themselves, so nothing is replayed, and a session put in view
 * starts the list again -- its own rows, not a divider under another's.
 */
import type { CliRenderer } from "@opentui/core"
import { createSignal } from "solid-js"
import { CommitQueue } from "./scrollback/queue"
import { Presenter } from "./scrollback/present"
import { replayOnWidthChange } from "./scrollback/resize"
import { itemWriter } from "./scrollback/rows"
import type { ScreenItems } from "./scrollback/screenview"
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
  /** Full screen: the rows' list in the window (none: the scrollback). */
  screen?: ScreenItems
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
  /** The last tool result of the session in view, in full (the palette). */
  expandLast(): void
  /** Every tool call's details shown. */
  details(): boolean
  /** Whether they can be shown (not under GNU screen, inline: the rows
   *  above cannot be written again). */
  readonly canShowDetails: boolean
  /** ctrl+o: every tool call's details on or off -- in place on the full
   *  screen; inline the session is written again with them (under GNU
   *  screen, which cannot, the last result is shown in full instead). */
  toggleDetails(): void
  /** Lines of the window's own (/keys), under a title. */
  lines(title: string, lines: string[]): void
  /** Thinking shown or hidden from now on. */
  setShowThinking(show: boolean): void
  readonly showThinking: boolean
  redraw(): Promise<void>
  stop(): void
}

/** How long ctrl+o waits, inline, before it writes the session again: a
 *  second press meanwhile takes the first back. */
export const DETAILS_SETTLE_MS = 150

/** Full screen's ctrl+l: every cell written again, so what a terminal left
 *  of an earlier frame goes (xterm leaves pixels of a block glyph beside
 *  it). OpenTUI keeps its own flag for that -- it repaints so after a
 *  resize or a resume -- and no call to set it. */
function repaint(renderer: CliRenderer): void {
  if (renderer.isDestroyed) return             // a timer after the window closed
  ;(renderer as unknown as { forceFullRepaintRequested: boolean }).forceFullRepaintRequested = true
  renderer.requestRender()
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
  // Inline: whether rows are written with their details (full screen: the
  // list's own, see screenview.tsx), and with which the scrollback was
  // last written as a whole.
  const [details, setDetails] = createSignal(false)
  let written = false
  let detailsTimer: ReturnType<typeof setTimeout> | undefined
  /** A finished block's rows are drawn: the footer lets go of it. */
  const release = (session: string, key: string) => {
    if (!store.state.sessions[session]?.held.some((h) => h.key === key)) return
    store.update((d) => {
      const held = d.sessions[session]?.held
      const at = held?.findIndex((h) => h.key === key) ?? -1
      if (at >= 0) held!.splice(at, 1)
    })
  }
  const screen = options.screen
  const transcript = screen
    ? new Transcript(queue, { resetSplitFooterForReplay: () => {} }, (item) => screen.push(item), false)
    : new Transcript(queue, target, itemWriter(renderer, theme, {
        release: (item) => release(String(item.session ?? ""), String(item.held ?? "")),
        details,
      }))
  const presentOptions = { showThinking: options.showThinking, shownAs: (text: string) => {
    const shown = sentAs.get(text)
    if (shown !== undefined) sentAs.delete(text)
    return shown
  }, release }

  /** When turns were left out of a replay: the banner of the session they
   *  belonged to, and a note saying how many. */
  const lead = (cut: number): Item[] => {
    const kept = transcript.items.filter((i) => (i.group ?? 0) >= cut + 1)
    const firstKept = transcript.items.indexOf(kept[0]!)
    const banner = [...transcript.items.slice(0, Math.max(0, firstKept))].reverse().find((i) => i.kind === "banner")
    const say = own("notice", { text: `(${cut} earlier turn${cut === 1 ? "" : "s"} not written again)` })
    return banner ? [banner, say] : [say]
  }
  const replay = () => transcript.replay(options.replayTurns, "turns", lead)
  const divider = options.dividerOnResize ?? Boolean(process.env.STY)
  // Full screen re-wraps by itself: no replay after a resize.
  const stopResize = screen ? () => {} : replayOnWidthChange(renderer, transcript, {
    onSettled: (width, widthChanged) => {
      if (divider) {
        // GNU screen cannot clear its scrollback: a replay would print
        // every row twice. A width change is said; a height change is
        // left as it is.
        if (widthChanged) {
          void transcript.append(own("divider", { text: `${width} columns now · the rows above keep their width` }))
        }
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
      if (screen && how.divider) {
        // Its own rows from the top: the list is the session's, not a log.
        groups.n++
        void queue.run(() => screen.clear())
      } else if (how.divider) {
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
    details: () => (screen ? screen.details : details()),
    get canShowDetails() {
      return Boolean(screen) || !divider
    },
    toggleDetails() {
      if (screen) screen.setDetails(!screen.details)
      else if (divider) presenter?.expandLast()
      else {
        setDetails(!details())
        // Pressed again and again: written once, when the presses stop --
        // a replay is the whole session (resizes settle the same way).
        clearTimeout(detailsTimer)
        detailsTimer = setTimeout(() => {
          if (details() === written) return
          written = details()
          void replay()
        }, DETAILS_SETTLE_MS)
        ;(detailsTimer as { unref?: () => void }).unref?.()
      }
    },
    lines(title, lines) {
      void transcript.append(own("lines", { title, lines }))
    },
    setShowThinking(show) {
      presentOptions.showThinking = show
    },
    get showThinking() {
      return presentOptions.showThinking
    },
    // ctrl+l: the session written again -- or, where the scrollback cannot
    // be cleared (GNU screen), the footer drawn again.
    redraw: () => (screen ? Promise.resolve(repaint(renderer))
      : divider ? Promise.resolve(renderer.requestRender()) : replay()),
    stop() {
      off?.()
      stopResize()
    },
  }
}
