/**
 * The scrollback machinery for the session in view, in one place: the
 * commit queue, the transcript, the presenter that turns the session's
 * events into rows, and the re-wrap after a width change. The window
 * (main.tsx) and the headless tests build it the same way.
 */
import type { CliRenderer } from "@opentui/core"
import { CommitQueue } from "./scrollback/queue"
import { Presenter } from "./scrollback/present"
import { replayOnWidthChange } from "./scrollback/resize"
import { itemWriter } from "./scrollback/rows"
import { Transcript } from "./scrollback/transcript"
import type { AppStore } from "./store/store"
import type { Theme } from "./theme/theme"

export interface SessionViewOptions {
  showThinking: boolean
  /** Items written again after a resize (0: all). */
  replayLimit: number
}

export interface SessionView {
  queue: CommitQueue
  transcript: Transcript
  /** Prompts sent with pastes put back -> as typed (see App's `sentAs`). */
  sentAs: Map<string, string>
  /** Put session `id` in view: its events become rows from now on. */
  show(id: string): void
  /** A note of the window's own (a config problem, a lost session), as a
   *  row in the scrollback, after what is already queued. */
  note(text: string): void
  /** ctrl+o: the last tool result of the session in view, in full. */
  expandLast(): void
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
  let notes = 0
  const note = (text: string) => {
    void transcript.append({ id: `note:${++notes}`, kind: "notice", text })
  }
  const queue = new CommitQueue(target, (error) => note(`could not draw a row: ${String(error)}`))
  const transcript = new Transcript(queue, target, itemWriter(renderer, theme))
  const stopResize = replayOnWidthChange(renderer, transcript, { limit: options.replayLimit })
  let off: (() => void) | null = null
  let presenter: Presenter | null = null
  const sentAs = new Map<string, string>()
  const shownAs = (text: string) => {
    const shown = sentAs.get(text)
    if (shown !== undefined) sentAs.delete(text)
    return shown
  }
  return {
    queue,
    transcript,
    sentAs,
    show(id) {
      off?.()
      const shown = new Presenter(id, transcript, { showThinking: options.showThinking, shownAs })
      presenter = shown
      off = store.onBatch((events, draft) => shown.handle(events, draft))
      store.update((d) => { d.current = id })
    },
    note,
    expandLast: () => presenter?.expandLast(),
    redraw: () => transcript.replay(options.replayLimit),
    stop() {
      off?.()
      stopResize()
    },
  }
}
