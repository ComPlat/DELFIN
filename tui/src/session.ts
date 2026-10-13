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
  /** Put session `id` in view: its events become rows from now on. */
  show(id: string): void
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
  const queue = new CommitQueue(target, (error) => store.update((d) => { d.messages.push(String(error)) }))
  const transcript = new Transcript(queue, target, itemWriter(renderer, theme))
  const stopResize = replayOnWidthChange(renderer, transcript, { limit: options.replayLimit })
  let off: (() => void) | null = null
  return {
    queue,
    transcript,
    show(id) {
      off?.()
      const presenter = new Presenter(id, transcript, { showThinking: options.showThinking })
      off = store.onBatch((events, draft) => presenter.handle(events, draft))
      store.update((d) => { d.current = id })
    },
    redraw: () => transcript.replay(options.replayLimit),
    stop() {
      off?.()
      stopResize()
    },
  }
}
