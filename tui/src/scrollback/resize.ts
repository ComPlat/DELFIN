/**
 * Re-wrapping the session after the terminal changes width.
 *
 * Committed rows are hard-wrapped at the width they were written at; no
 * terminal re-wraps them. After a WIDTH change settles, the screen and the
 * terminal's scrollback are cleared and the transcript is written again at
 * the new width (the user's choice, docs/TUI_REWRITE.md).
 *
 * Only the width counts. In split-footer mode OpenTUI's ``resize`` event
 * also fires when the footer changes height, and a height change needs no
 * re-wrap; treating it as one wiped the scrollback on every footer change
 * (measured in phase 0).
 */
import type { Transcript } from "./transcript"

export interface WidthSource {
  readonly terminalWidth: number
  /** The terminal's own height (not the footer's): a change of it moves
   *  the footer, and OpenTUI's split footer leaves rows behind (measured
   *  in tmux), so it is replayed too. */
  readonly terminalHeight?: number
  on(event: "resize", listener: () => void): unknown
  off(event: "resize", listener: () => void): unknown
}

export interface ReplayPolicy {
  /** Quiet time after the last width change before replaying. */
  settleMs?: number
  /** Items to write again; older ones are left out (0 = all). */
  limit?: number
  /** Instead of the default replay: what to do once the width settled
   *  (the session view replays by turns, or prints a divider where the
   *  terminal cannot clear its scrollback). */
  onSettled?: (width: number, widthChanged: boolean) => void
}

/** Watch ``source``; replay ``transcript`` after each settled width change. */
export function replayOnWidthChange(
  source: WidthSource,
  transcript: Transcript,
  policy: ReplayPolicy = {},
): () => void {
  const settleMs = policy.settleMs ?? 250
  let width = source.terminalWidth
  let height = source.terminalHeight
  let timer: ReturnType<typeof setTimeout> | undefined
  const same = () => source.terminalWidth === width && source.terminalHeight === height
  const onResize = () => {
    if (same()) {
      // Back to the size the scrollback was written at (or only the
      // footer changed height): nothing to redo.
      clearTimeout(timer)
      timer = undefined
      return
    }
    clearTimeout(timer)
    timer = setTimeout(() => {
      timer = undefined
      if (same()) return
      const widthChanged = source.terminalWidth !== width
      width = source.terminalWidth
      height = source.terminalHeight
      if (policy.onSettled) policy.onSettled(width, widthChanged)
      else void transcript.replay(policy.limit ?? 0)
    }, settleMs)
  }
  source.on("resize", onResize)
  return () => {
    clearTimeout(timer)
    source.off("resize", onResize)
  }
}
