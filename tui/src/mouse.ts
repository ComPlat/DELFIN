/**
 * A click, as the window takes one: the left button pressed and let go on
 * the same cell of one row, within CLICK_MAX_MS. Not a release that ends a
 * drag begun elsewhere, not another button, not the wheel, never a move
 * (the renderer reports no moves: enableMouseMovement is off, so nothing
 * hovers). OpenTUI starts a selection when the button goes down on text;
 * a press and release on one cell selects nothing.
 *
 * What the click means is the row's: `actAt(col)` names what is under the
 * column (an option's key, a list row), and the press and the release must
 * name the same thing.
 *
 * Mouse reports are on only in the full-screen window (main.tsx): inline,
 * the terminal's own scrollback and selection stay the terminal's. With
 * them on, shift+drag selects with the terminal's own selection.
 */
import type { MouseEvent, Renderable } from "@opentui/core"

/** Longer than this between press and release is a hold, not a click. */
export const CLICK_MAX_MS = 1500
/** A click this soon after what it would answer appeared (a question, a
 *  list, a second step) is not taken: the hand was already on its way. */
export const CLICK_PAUSE_MS = 300

export interface ClickHandlers {
  ref: (r: Renderable) => void
  onMouseDown: (event: MouseEvent) => void
  onMouseUp: (event: MouseEvent) => void
  onMouseDrag: (event: MouseEvent) => void
}

export function clickable(actAt: (col: number) => string | null, onClick: (act: string) => void,
                          now: () => number = Date.now): ClickHandlers {
  let row: Renderable | undefined
  let press: { x: number; y: number; at: number; act: string } | null = null
  const at = (event: MouseEvent) => (row ? actAt(event.x - row.x) : null)
  return {
    ref: (r) => { row = r },
    onMouseDown: (event) => {
      const act = event.button === 0 ? at(event) : null
      press = act ? { x: event.x, y: event.y, at: now(), act } : null
    },
    onMouseDrag: (event) => {
      // Moved while held: a drag (a selection), never a click.
      if (press && (event.x !== press.x || event.y !== press.y)) press = null
    },
    onMouseUp: (event) => {
      const was = press
      press = null
      if (!was || event.button !== 0 || event.x !== was.x || event.y !== was.y) return
      if (now() - was.at > CLICK_MAX_MS || at(event) !== was.act) return
      onClick(was.act)
    },
  }
}

/** What a row of spans has under column `col`: the `act` of the span
 *  there (its text, not the spaces after it), or null. */
export function spanAt<T extends { text: string; act?: string }>(spans: readonly T[], col: number,
                                                                  widthOf: (s: string) => number): string | null {
  let x = 0
  for (const span of spans) {
    const w = widthOf(span.text)
    if (col >= x && col < x + widthOf(span.text.trimEnd())) return span.act ?? null
    x += w
  }
  return null
}
