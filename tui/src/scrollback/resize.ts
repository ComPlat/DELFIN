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
  onSettled?: (width: number) => void
}

/** Watch ``source``; replay ``transcript`` after each settled width change. */
export function replayOnWidthChange(
  source: WidthSource,
  transcript: Transcript,
  policy: ReplayPolicy = {},
): () => void {
  const settleMs = policy.settleMs ?? 250
  let width = source.terminalWidth
  let timer: ReturnType<typeof setTimeout> | undefined
  const onResize = () => {
    if (source.terminalWidth === width) {
      // Back to the width the scrollback was written at: nothing to redo.
      clearTimeout(timer)
      timer = undefined
      return
    }
    clearTimeout(timer)
    timer = setTimeout(() => {
      timer = undefined
      if (source.terminalWidth === width) return
      width = source.terminalWidth
      if (policy.onSettled) policy.onSettled(width)
      else void transcript.replay(policy.limit ?? 0)
    }, settleMs)
  }
  source.on("resize", onResize)
  return () => {
    clearTimeout(timer)
    source.off("resize", onResize)
  }
}
