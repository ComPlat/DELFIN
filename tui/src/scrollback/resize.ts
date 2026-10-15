/**
 * Re-wrapping the session after the terminal changes width.
 *
 * Committed rows are hard-wrapped at the width they were written at; no
 * terminal re-wraps them. After a WIDTH change settles, the screen and the
 * terminal's scrollback are cleared and the transcript is written again at
 * the new width (the user's choice, docs/TUI_REWRITE.md).
 *
 * Only the terminal's size counts. In split-footer mode OpenTUI's
 * ``resize`` event also fires when the footer changes height, and a footer
 * change needs no re-wrap; treating it as one wiped the scrollback on every
 * footer change (measured in phase 0).
 *
 * A size that comes back to the one the scrollback was written at is
 * replayed too, once it settles. Between the two the window was drawn at
 * another size: OpenTUI clears the rows above the footer when the width
 * shrinks, a reflowing terminal (tmux, VTE, kitty) rewraps every row, and
 * a replay still being written when the terminal changed width was written
 * at a width the terminal no longer had -- OpenTUI holds a new size back
 * until the queued rows are out. tmux passes a pane's new size on at most
 * every 250 ms, so a resize storm there often settles once during the storm;
 * one that ended at the width of such a replay kept that replay's misplaced
 * footer rows in the scrollback (test/e2e/check_resize_tmux.py). The same
 * holds for a size the terminal reported that OpenTUI never applied (see
 * ``ReplayPolicy.reports``).
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
  /** The terminal's own size reports (SIGWINCH), as a subscription that
   *  returns its unsubscribe. OpenTUI holds a size reported while queued
   *  rows are being written and applies only the last one, with no
   *  ``resize`` event when that is the size it already has: a pane that
   *  went to another width and back during a replay is seen only here. */
  reports?: (listener: (width: number, height: number) => void) => () => void
}

/** Watch ``source``; replay ``transcript`` after each settled width change. */
export function replayOnWidthChange(
  source: WidthSource,
  transcript: Transcript,
  policy: ReplayPolicy = {},
): () => void {
  const settleMs = policy.settleMs ?? 250
  // The size the scrollback was last written at, and the size of the last
  // event: an event at the size of the one before is the footer's own.
  let width = source.terminalWidth
  let height = source.terminalHeight
  let seenWidth = width
  let seenHeight = height
  let reportedWidth = width
  let reportedHeight = height
  let timer: ReturnType<typeof setTimeout> | undefined
  const onResize = () => {
    if (source.terminalWidth === seenWidth && source.terminalHeight === seenHeight) return
    seenWidth = source.terminalWidth
    seenHeight = source.terminalHeight
    arm()
  }
  const onReport = (w: number, h: number) => {
    if (w === reportedWidth && h === reportedHeight) return
    reportedWidth = w
    reportedHeight = h
    arm()
  }
  const arm = () => {
    clearTimeout(timer)
    timer = setTimeout(() => {
      timer = undefined
      const widthChanged = source.terminalWidth !== width
      width = source.terminalWidth
      height = source.terminalHeight
      if (policy.onSettled) policy.onSettled(width, widthChanged)
      else void transcript.replay(policy.limit ?? 0)
    }, settleMs)
  }
  source.on("resize", onResize)
  const stopReports = policy.reports?.(onReport)
  return () => {
    clearTimeout(timer)
    source.off("resize", onResize)
    stopReports?.()
  }
}
