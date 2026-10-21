/**
 * The agent view on the terminal's alternate screen: opening it leaves
 * the session's scrollback exactly as it is, and closing it brings the
 * screen back as it was -- the way `less` or an editor leave a shell.
 *
 * Drawn in the footer instead (it grows to the terminal's height), the
 * view would push the session's last rows into the terminal's history,
 * and on closing the prompt would sit at the top of an empty screen.
 *
 * While the alternate screen is up nothing is written to the scrollback
 * (OpenTUI writes it only in split-footer mode): the commit queue is held
 * and what came meanwhile is written on the way back. Back in split-footer
 * mode, OpenTUI seeds its scrollback from the cursor, which on the
 * alternate screen stood somewhere else: it is told where the footer was.
 *
 * Three things the way there and back waits for, each a lost or a stale
 * row when it did not:
 * - before the seed is taken, OpenTUI has *presented* what was committed
 *   (commitRows only queues rows; the footer's offset moves when a frame
 *   writes them);
 * - the way back restores the footer at the height it had when the view
 *   opened, and applies the height it wants now only after that -- so
 *   OpenTUI's own resize path scrolls rows into history (it grew) or
 *   clears the ones no longer covered (it shrank);
 * - the scrollback is let go only once split-footer mode is applied
 *   (OpenTUI defers a switch asked for while a frame is drawn), and is let
 *   go even when a switch fails.
 */
import type { CliRenderer } from "@opentui/core"
import type { CommitQueue } from "../scrollback/queue"

export interface AgentScreen {
  /** On the alternate screen now. */
  readonly active: boolean
  /** Go to the alternate screen (once what is queued is written). */
  enter(): Promise<void>
  /** Come back, and write what was held meanwhile. */
  leave(): void
}

interface SplitInternals {
  renderOffset: number
  rendering: boolean
  _screenMode: string
  _splitHeight: number
  _footerHeight: number
  externalOutputQueue: { size: number }
  getSplitCursorSeedRows(): number
  /** The screen mode, without the public setter's rule that captured
   *  stdout needs split-footer mode: the capture stays on throughout (no
   *  scrollback is written meanwhile), so the way back is one reset. */
  applyScreenMode(mode: "alternate-screen" | "split-footer", emitResize?: boolean, requestRender?: boolean,
                  footerHeight?: number): void
}

/** How long the seed stays told after the way back is applied. */
const SEED_HOLD_MS = 500
/** The longest wait for OpenTUI to present, or to apply a mode. */
const WAIT_STEPS = 400
const STEP_MS = 5

export function createAgentScreen(renderer: CliRenderer, queue: CommitQueue): AgentScreen {
  const internals = renderer as unknown as SplitInternals
  let want = false
  let shown = false
  let release: (() => void) | null = null
  let seedTimer: ReturnType<typeof setTimeout> | undefined
  let busy: Promise<void> = Promise.resolve()
  let heightAtShow = 0

  const tellSeed = (rows: number) => {
    clearTimeout(seedTimer)
    internals.getSplitCursorSeedRows = () => Math.max(0, Math.min(rows, renderer.terminalHeight))
  }
  const forgetSeed = () => {
    // The renderer's own method again (it lives on the prototype).
    delete (internals as unknown as Record<string, unknown>).getSplitCursorSeedRows
  }
  const step = () => new Promise<void>((resolve) => setTimeout(resolve, STEP_MS))
  /** Until OpenTUI has written what is committed and no frame is drawn. */
  const presented = async () => {
    for (let i = 0; i < WAIT_STEPS && (internals.externalOutputQueue.size > 0 || internals.rendering); i++) {
      renderer.requestRender()
      await step()
    }
  }
  /** Until OpenTUI has applied `mode` (it defers one asked for mid-frame). */
  const applied = async (mode: string) => {
    for (let i = 0; i < WAIT_STEPS && internals._screenMode !== mode; i++) {
      renderer.requestRender()
      await step()
    }
  }
  const letGo = () => {
    const done = release
    release = null
    done?.()
  }

  const show = async () => {
    let started!: () => void
    const ready = new Promise<void>((resolve) => { started = resolve })
    const held = new Promise<void>((resolve) => { release = resolve })
    // Everything queued before is written first; everything after waits.
    void queue.run(() => {
      started()
      return held
    })
    await ready
    try {
      await presented()
      heightAtShow = internals._splitHeight
      // The rows above the footer now: where it goes back to. (The reset
      // takes the rows in use, not the footer's 1-based top line.)
      tellSeed(internals.renderOffset)
      internals.applyScreenMode("alternate-screen")
      await applied("alternate-screen")
      shown = true
    } catch (error) {
      forgetSeed()
      letGo()
      throw error
    }
  }
  const hide = async () => {
    // What the footer wants now (the window kept setting it meanwhile; on
    // the alternate screen that is only stored).
    const wanted = internals._footerHeight
    try {
      internals.applyScreenMode("split-footer", true, true, heightAtShow || wanted)
      await applied("split-footer")
      if (renderer.footerHeight !== wanted) renderer.footerHeight = wanted
      renderer.requestRender()
    } finally {
      shown = internals._screenMode !== "split-footer"
      seedTimer = setTimeout(forgetSeed, SEED_HOLD_MS)
      ;(seedTimer as { unref?: () => void }).unref?.()
      letGo()
    }
  }
  // One change at a time, until the screen is where it is wanted (a quick
  // ← esc ← ends on the alternate screen).
  const settle = () => {
    busy = busy.then(async () => {
      for (let round = 0; want !== shown && round < 8; round++) {
        if (want) await show()
        else await hide()
      }
    }).catch(() => {
      // A renderer gone meanwhile (the window left): nothing to switch.
    })
    return busy
  }

  return {
    get active() {
      return shown
    },
    enter() {
      want = true
      return settle()
    },
    leave() {
      want = false
      void settle()
    },
  }
}
