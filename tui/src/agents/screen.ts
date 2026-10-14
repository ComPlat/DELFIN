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
  getSplitCursorSeedRows(): number
  /** The screen mode, without the public setter's rule that captured
   *  stdout needs split-footer mode: the capture stays on throughout (no
   *  scrollback is written meanwhile), so the way back is one reset. */
  applyScreenMode(mode: "alternate-screen" | "split-footer"): void
}

/** How long the seed stays told after the way back: the modes are applied
 *  at the next frame, which may wait for one being drawn. */
const SEED_HOLD_MS = 500

export function createAgentScreen(renderer: CliRenderer, queue: CommitQueue): AgentScreen {
  const internals = renderer as unknown as SplitInternals
  let want = false
  let shown = false
  let release: (() => void) | null = null
  let seedTimer: ReturnType<typeof setTimeout> | undefined
  let busy: Promise<void> = Promise.resolve()

  const tellSeed = (rows: number) => {
    clearTimeout(seedTimer)
    internals.getSplitCursorSeedRows = () => Math.max(0, Math.min(rows, renderer.terminalHeight))
  }
  const forgetSeed = () => {
    // The renderer's own method again (it lives on the prototype).
    delete (internals as unknown as Record<string, unknown>).getSplitCursorSeedRows
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
    // The rows above the footer now: where it goes back to. (The reset
    // takes the rows in use, not the footer's 1-based top line.)
    tellSeed(internals.renderOffset)
    internals.applyScreenMode("alternate-screen")
    shown = true
  }
  const hide = () => {
    internals.applyScreenMode("split-footer")
    shown = false
    renderer.requestRender()
    seedTimer = setTimeout(forgetSeed, SEED_HOLD_MS)
    ;(seedTimer as { unref?: () => void }).unref?.()
    const done = release
    release = null
    done?.()
  }
  // One change at a time, until the screen is where it is wanted (a quick
  // ← esc ← ends on the alternate screen).
  const settle = () => {
    busy = busy.then(async () => {
      while (want !== shown) {
        if (want) await show()
        else hide()
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
