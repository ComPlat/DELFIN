/**
 * The answer's live blocks in the footer: the blocks that finished but are
 * not in the rows yet (held), then the unfinished one -- which of them are
 * drawn, and how many rows they take.
 *
 * The rows are counted, not measured: inline, the footer's height must be
 * known before the frame that draws it, and a box sized a frame late would
 * show the text a row off for that frame. A block's rows are its lines as
 * drawn -- a heading without its #s, emphasis without its marks, a link as
 * its text, a fence's lines not at all -- each wrapped at a word, a column
 * short of the width (a count a row too high leaves a blank row above the
 * text; one too low would cut the top of a block being handed over).
 */
import type { SessionState } from "../store/state"
import { width } from "../text"

/** Rows the live blocks take at most: what is above is cut off. */
export const TAIL_MAX = 6

export interface LiveBlock {
  key: string
  text: string
}

/** A fenced code block, as the markdown renderer draws it: its lines,
 *  without the fences, highlighted for its language. */
export function fenceOf(text: string): { lang: string; code: string } | null {
  const open = /^ {0,3}(`{3,}|~{3,})([^\n]*)(?:\n|$)/.exec(text)
  if (!open) return null
  const mark = open[1]!
  if (mark[0] === "`" && open[2]!.includes("`")) return null        // inline code, not a fence
  let body = text.slice(open[0].length).replace(/\n$/, "")
  const close = new RegExp(`(?:^|\\n) {0,3}\\${mark[0]}{${mark.length},}[ \\t]*$`).exec(body)
  if (close) body = body.slice(0, close.index)
  return { lang: open[2]!.trim(), code: body }
}

/** A line of markdown as drawn: no heading #s, quote marks, emphasis
 *  marks or link targets. */
function drawn(line: string): string {
  return line
    .replace(/^ {0,3}#{1,6}(?:[ \t]+|$)/, "")
    .replace(/[ \t]+#+[ \t]*$/, "")
    .replace(/^\s*>\s?/, "")
    .replace(/!?\[([^\]]*)\]\([^)]*\)/g, "$1")
    .replace(/(\*\*|__|~~)/g, "")
    .replace(/`+/g, "")
}

/** Rows `line` takes wrapped at a word to `cols` columns (its indent
 *  counted, a tab as four). */
export function wrapRows(line: string, cols: number): number {
  if (cols < 2) return 1
  const text = line.replace(/\t/g, "    ")
  const indent = text.length - text.trimStart().length
  let rows = 1
  let at = Math.min(indent, cols - 1)
  for (const word of text.trimStart().split(/ +/)) {
    const w = width(word)
    if (at > 0 && at + 1 + w <= cols) {
      at += 1 + w
      continue
    }
    if (at > indent || (at > 0 && at + w > cols)) rows++
    // A word longer than a row is cut into rows of its own.
    rows += Math.max(0, Math.ceil(w / cols) - 1)
    at = w > cols ? w % cols || cols : w
  }
  return rows
}

/** The fence that opens a code block on `line`, if one does (backticks
 *  in a ``` fence's info string make it inline code instead). */
function fenceMark(line: string): string | null {
  const open = /^ {0,3}(`{3,}|~{3,})(.*)$/.exec(line)
  if (!open) return null
  return open[1]![0] === "`" && open[2]!.includes("`") ? null : open[1]!
}

/** Rows a block takes at `cols` columns. */
export function blockRows(text: string, cols: number): number {
  const room = Math.max(2, cols - 1)
  const lines = text.replace(/\s+$/, "").split("\n")
  while (lines.length > 0 && !lines[0]!.trim()) lines.shift()
  // A table: a border above and below, and a rule between every two rows.
  const table = lines.filter((l) => l.trim().startsWith("|")).length
  if (table >= 2 && table === lines.filter((l) => l.trim()).length) return 2 * table - 1
  let rows = 0
  let fence: string | null = null
  for (const line of lines) {
    const mark = fenceMark(line)
    if (fence === null && mark) {
      fence = mark
      continue
    }
    if (fence !== null && mark && mark[0] === fence[0] && mark.length >= fence.length && !line.trim().slice(mark.length)) {
      fence = null
      continue
    }
    if (fence !== null) rows += wrapRows(line, room)
    else if (!line.trim()) rows += 1
    // A quote's bar and its space take two columns.
    else rows += wrapRows(drawn(line), /^\s*>/.test(line) ? Math.max(2, room - 2) : room)
  }
  return rows
}

/** The live blocks to draw, by key: held ones, then the unfinished one
 *  (less a heading's "##" that has no words yet -- it would be drawn as it
 *  is until they come). With `cols`, only the trailing ones that fill the
 *  footer's TAIL_MAX rows: those above would be cut off anyway, and each
 *  drawn costs a highlight (held blocks pile up while the rows cannot be
 *  written: the agent view has the terminal). */
export function liveBlocks(s: SessionState | null, cols = 0): LiveBlock[] {
  if (!s) return []
  const blocks: LiveBlock[] = [...(s.held ?? [])]
  const tail = s.tail.replace(/(^|\n) {0,3}#{1,6}[ \t]*$/, "$1")
  if (tail.trim() && !blocks.some((b) => b.key === s.tailKey)) blocks.push({ key: s.tailKey, text: tail })
  if (cols <= 0) return blocks
  let rows = 0
  let from = blocks.length
  while (from > 0 && rows < TAIL_MAX) {
    from--
    rows += blockRows(blocks[from]!.text, cols) + (from < blocks.length - 1 ? 1 : 0)
  }
  return blocks.slice(from)
}

/** Rows the live blocks take (a row between two), at most TAIL_MAX. */
export function liveRows(s: SessionState | null, cols: number): number {
  const blocks = liveBlocks(s, cols)
  const rows = blocks.reduce((n, b, i) => n + blockRows(b.text, cols) + (i > 0 ? 1 : 0), 0)
  return Math.min(TAIL_MAX, rows)
}
