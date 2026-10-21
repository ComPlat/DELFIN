/**
 * A pipe table in an answer, laid out for the window by DELFIN rather than
 * by OpenTUI's markdown renderer.
 *
 * Input: the lines of a markdown block (blocks.ts) and the columns there
 * are. Output: the block cut into markdown and table segments
 * (`tableSegments`), and a table as rows of styled spans (`drawTable`), or
 * null when it cannot fit -- the caller then shows the source lines.
 *
 * Recognition is delfin/agent/markdown_table.py's, ported line by line, so
 * the dashboard, the plain terminal and the window take the same lines for
 * a table: the alignment row must have the header's cell count, an escaped
 * \| or a pipe inside a code span stays in its cell, short rows are padded,
 * a column holding only numbers is right-aligned unless the alignment row
 * says otherwise, and a table ends at its first line without a pipe.
 *
 * Layout: a box (┌┬┐ │ ├┼┤ └┴┘), one space of padding in each cell, a rule
 * under the header and -- when a body row takes more than one line --
 * between body rows. A token that holds a number (a numeric cell whole; a
 * word that starts with a digit or a sign and a digit) is never split. When
 * the table is wider than the columns, text columns give way first, at
 * word boundaries (widest first); then a number column's header, at its
 * words; then words are cut (after a "/", "-", "_", "." where one fits),
 * down to MIN_COLUMN. Still too wide: null.
 *
 * Rejected: patching OpenTUI's table (TextTableRenderable and the list
 * coalescing in MarkdownRenderable). It has no column alignment and wraps
 * every cell alike, and a table in a list item never reaches it -- the
 * list is drawn as highlighted markdown source. The patch would be edits
 * to two generated bundles (chunk-bun, chunk-node) redone at every OpenTUI
 * release.
 */
import { width as cols } from "../text"

export type Align = "l" | "c" | "r" | ""

const DELIM_CELL = /^:?-+:?$/
const NUMERIC = /^[~≈±+\-−]?\s?\d[\d.,'  ]*(?:[eE][+\-−]?\d+)?\s?%?$/u

/** The cells of a pipe-table row, inline markup left in place. */
export function splitRow(line: string): string[] {
  let s = line.trim()
  if (s.startsWith("|")) s = s.slice(1)
  if (s.endsWith("|") && !s.endsWith("\\|")) s = s.slice(0, -1)
  const cells: string[] = []
  let buf = ""
  let codeRun = 0
  let i = 0
  while (i < s.length) {
    const ch = s[i]!
    if (ch === "\\" && s[i + 1] === "|") {
      buf += codeRun ? "|" : "\\|"
      i += 2
      continue
    }
    if (ch === "`") {
      let j = i
      while (j < s.length && s[j] === "`") j++
      const run = j - i
      if (codeRun === 0) {
        if (s.indexOf("`".repeat(run), j) >= 0) codeRun = run
      } else if (run === codeRun) {
        codeRun = 0
      }
      buf += s.slice(i, j)
      i = j
      continue
    }
    if (ch === "|" && codeRun === 0) {
      cells.push(buf)
      buf = ""
      i++
      continue
    }
    buf += ch
    i++
  }
  cells.push(buf)
  return cells.map((c) => c.trim())
}

/** Whether `line` holds an unescaped pipe. */
export function hasPipe(line: string): boolean {
  return /(^|[^\\])\|/.test(line)
}

/** Column alignments when `header` over `delim` opens a table, else null. */
export function alignment(header: string, delim: string): Align[] | null {
  if (!hasPipe(header)) return null
  if (!delim.includes("-") || !/^[\s|:\-]+$/.test(delim)) return null
  const dcells = splitRow(delim)
  if (dcells.length === 0 || !dcells.every((c) => DELIM_CELL.test(c))) return null
  if (dcells.length !== splitRow(header).length) return null
  return dcells.map((c) => c.startsWith(":") && c.endsWith(":") ? "c" : c.endsWith(":") ? "r" : c.startsWith(":") ? "l" : "")
}

/** A cell holding one number (sign, digits, separators, exponent, %). */
export function isNumeric(cell: string): boolean {
  return NUMERIC.test(cell.replace(/`/g, "").trim())
}

/** Every row as wide as the widest; a column of numbers aligned right. */
export function normalise(header: string[], rows: string[][], aligns: Align[]) {
  const n = Math.max(header.length, ...rows.map((r) => r.length))
  const pad = (r: string[]) => [...r, ...Array(n - r.length).fill("")]
  const body = rows.map(pad)
  const al = [...aligns, ...Array(Math.max(0, n - aligns.length)).fill("")] as Align[]
  for (let c = 0; c < n; c++) {
    if (al[c]) continue
    const values = body.map((r) => r[c]!).filter((v) => v)
    if (values.length > 0 && values.every(isNumeric)) al[c] = "r"
  }
  return { header: pad(header), rows: body, aligns: al }
}

// -- segments ----------------------------------------------------------------

export interface Segment {
  kind: "markdown" | "table"
  /** The segment's source, as written (a table's lines with their indent). */
  text: string
  /** Columns a table stands in from the block's left edge. */
  indent: number
  /** A blank line stood before this segment in the source. */
  gapBefore: boolean
}

const LIST_ITEM = /^ {0,3}(?:[-*+]|\d{1,9}[.)])(?:[ \t]|$)/
const FENCE = /^\s*(`{3,}|~{3,})/

const indentOf = (line: string) => line.length - line.trimStart().length

/**
 * The block cut at every table it holds. A table counts at an indent of
 * up to three columns, and deeper inside a list (the block is a list: it
 * starts with an item) -- deeper outside one it is an indented code block.
 * Lines inside a fence and quoted lines are left to the markdown renderer.
 * A block with no table is one markdown segment, its text unchanged.
 */
export function tableSegments(block: string): Segment[] {
  const lines = block.split("\n")
  const inList = LIST_ITEM.test(lines[0] ?? "")
  const out: Segment[] = []
  let from = 0
  let fence: string | null = null
  const flush = (to: number) => {
    let a = from
    let b = to
    while (a < b && !lines[a]!.trim()) a++
    while (b > a && !lines[b - 1]!.trim()) b--
    if (b > a) out.push({ kind: "markdown", text: lines.slice(a, b).join("\n"), indent: 0, gapBefore: a > from })
  }
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]!
    const mark = FENCE.exec(line)
    if (fence !== null) {
      if (mark && mark[1]![0] === fence[0] && mark[1]!.length >= fence.length && !line.trim().slice(mark[1]!.length).trim()) fence = null
      continue
    }
    if (mark) {
      fence = mark[1]!
      continue
    }
    if (i + 1 >= lines.length || line.trimStart().startsWith(">")) continue
    const indent = indentOf(line)
    if (indent > 3 && !inList) continue
    if (!alignment(line, lines[i + 1]!)) continue
    let end = i + 2
    while (end < lines.length && lines[end]!.trim() && hasPipe(lines[end]!) && !FENCE.test(lines[end]!)) end++
    flush(i)
    out.push({ kind: "table", text: lines.slice(i, end).join("\n"), indent,
               gapBefore: i > 0 && !lines[i - 1]!.trim() })
    from = end
    i = end - 1
  }
  if (out.length === 0) return [{ kind: "markdown", text: block, indent: 0, gapBefore: false }]
  flush(lines.length)
  return out
}

// -- cells -------------------------------------------------------------------

export type SpanStyle = "border" | "header" | "text" | "strong" | "code"

export interface Span {
  text: string
  style: SpanStyle
}

/** One character of a cell, with its style and its columns. */
interface Glyph {
  ch: string
  style: SpanStyle
  w: number
}

const GRAPHEMES = new Intl.Segmenter(undefined, { granularity: "grapheme" })

const MARKUP = /`([^`]+)`|\*\*(.+?)\*\*|__(.+?)__|\[([^\]\n]+)\]\([^)\s]+\)|\\([\\|*_`[\]])/g

/** A cell's visible characters: code spans, strong text and links without
 *  their markers, escapes as the character they escape. */
function glyphs(cell: string, base: SpanStyle): Glyph[] {
  const out: Glyph[] = []
  const push = (text: string, style: SpanStyle) => {
    for (const { segment } of GRAPHEMES.segment(text)) out.push({ ch: segment, style, w: cols(segment) })
  }
  let at = 0
  for (const m of cell.matchAll(MARKUP)) {
    push(cell.slice(at, m.index), base)
    if (m[1] !== undefined) push(m[1], base === "header" ? "header" : "code")
    else if (m[2] !== undefined || m[3] !== undefined) push((m[2] ?? m[3])!, base === "header" ? "header" : "strong")
    else if (m[4] !== undefined) push(m[4], base)
    else push(m[5]!, base)
    at = m.index! + m[0].length
  }
  push(cell.slice(at), base)
  return out
}

interface Word {
  glyphs: Glyph[]
  w: number
  /** A number, or the whole of a number cell: never split. */
  rigid: boolean
}

const NUMBER_WORD = /^[(~≈±+\-−]?\d/u

function words(g: Glyph[], wholeCell: boolean): Word[] {
  if (wholeCell) {
    const trimmed = trimGlyphs(g)
    return trimmed.length ? [{ glyphs: trimmed, w: sum(trimmed), rigid: true }] : []
  }
  const out: Word[] = []
  let cur: Glyph[] = []
  const end = () => {
    if (cur.length) {
      const text = cur.map((x) => x.ch).join("")
      out.push({ glyphs: cur, w: sum(cur), rigid: NUMBER_WORD.test(text) })
    }
    cur = []
  }
  for (const x of g) {
    if (x.ch === " " || x.ch === "\t") end()
    else cur.push(x)
  }
  end()
  return out
}

function trimGlyphs(g: Glyph[]): Glyph[] {
  let a = 0
  let b = g.length
  while (a < b && g[a]!.ch.trim() === "") a++
  while (b > a && g[b - 1]!.ch.trim() === "") b--
  return g.slice(a, b)
}

const sum = (g: Glyph[]) => g.reduce((n, x) => n + x.w, 0)

const CUT_AFTER = new Set(["/", "-", "_", ".", ","])

/** The narrowest a column is drawn: narrower, a cell of prose turns into
 *  a column of syllables, and the source reads better. */
const MIN_COLUMN = 6

/** The lines of a cell at `w` columns: words packed greedily, a word
 *  wider than `w` cut into pieces (never a rigid one: `w` is at least the
 *  widest of those). */
function wrapCell(ws: Word[], w: number): Glyph[][] {
  const lines: Glyph[][] = []
  let line: Glyph[] = []
  let lw = 0
  const space = (style: SpanStyle): Glyph => ({ ch: " ", style, w: 1 })
  for (const word of ws) {
    if (lw > 0 && lw + 1 + word.w <= w) {
      line.push(space(word.glyphs[0]!.style), ...word.glyphs)
      lw += 1 + word.w
      continue
    }
    if (lw > 0) {
      lines.push(line)
      line = []
      lw = 0
    }
    if (word.w <= w) {
      line = [...word.glyphs]
      lw = word.w
      continue
    }
    // Cut after the last "/", "-", "_" or "." that fits, else where the
    // row ends ("PBE0/" "def2-SVP" rather than "PBE0/d" "ef2-SVP").
    let rest = word.glyphs
    while (sum(rest) > w) {
      let fit = 0
      let used = 0
      let cut = 0
      for (const x of rest) {
        if (used + x.w > w) break
        used += x.w
        fit++
        if (CUT_AFTER.has(x.ch) && fit < rest.length) cut = fit
      }
      const at = cut > 0 ? cut : Math.max(1, fit)
      lines.push(rest.slice(0, at))
      rest = rest.slice(at)
    }
    line = [...rest]
    lw = sum(rest)
  }
  if (line.length || lines.length === 0) lines.push(line)
  return lines
}

// -- layout ------------------------------------------------------------------

/**
 * The table in `lines` (header, alignment row, body) drawn in at most
 * `width` columns, one array of spans per screen row; null when the lines
 * are not a table or it cannot fit without splitting a number.
 */
export function drawTable(lines: string[], width: number): Span[][] | null {
  const src = lines.map((l) => l.trim()).filter((l) => l)
  if (src.length < 2) return null
  const aligns0 = alignment(src[0]!, src[1]!)
  if (!aligns0) return null
  const { header, rows, aligns } = normalise(splitRow(src[0]!), src.slice(2).map(splitRow), aligns0)
  const n = header.length
  const numeric = aligns.map((_, c) => {
    const values = rows.map((r) => r[c]!).filter((v) => v)
    return values.length > 0 && values.every(isNumeric)
  })
  const head = header.map((h) => words(glyphs(h, "header"), false))
  const body = rows.map((r) => r.map((v, c) => words(glyphs(v, "text"), numeric[c]! && isNumeric(v))))
  const all = (c: number) => [head[c]!, ...body.map((r) => r[c]!)]
  const natural = (ws: Word[]) => ws.reduce((m, x, i) => m + x.w + (i > 0 ? 1 : 0), 0)
  const want = Array.from({ length: n }, (_, c) => Math.max(1, ...all(c).map(natural)))
  // The narrowest a column may get: at word boundaries (soft), and at all (hard).
  const soft = Array.from({ length: n }, (_, c) => Math.max(1, ...all(c).flat().map((x) => x.w)))
  const hard = Array.from({ length: n }, (_, c) => Math.max(1, Math.min(MIN_COLUMN, want[c]!),
                                                           ...all(c).flat().filter((x) => x.rigid).map((x) => x.w)))
  const overhead = 3 * n + 1
  const widths = [...want]
  const total = () => widths.reduce((a, b) => a + b, 0) + overhead
  // Text columns first, down to their longest word; then a number
  // column's header down to its longest word; then words are cut, text
  // columns first. Widest first each time.
  const stages: Array<[number[], (c: number) => boolean]> = [
    [soft, (c) => !numeric[c]],
    [soft, (c) => numeric[c]!],
    [hard, (c) => !numeric[c]],
    [hard, (c) => numeric[c]!],
  ]
  for (const [floor, may] of stages) {
    while (total() > width) {
      let pick = -1
      for (let c = 0; c < n; c++) {
        if (!may(c) || widths[c]! <= floor[c]!) continue
        if (pick < 0 || widths[c]! > widths[pick]!) pick = c
      }
      if (pick < 0) break
      widths[pick]!--
    }
  }
  if (total() > width) return null

  const border = (l: string, m: string, r: string): Span[] =>
    [{ text: l + widths.map((w) => "─".repeat(w + 2)).join(m) + r, style: "border" }]
  const row = (cells: Word[][]): Span[][] => {
    const wrapped = cells.map((ws, c) => wrapCell(ws, widths[c]!))
    const height = Math.max(...wrapped.map((w) => w.length))
    const out: Span[][] = []
    for (let k = 0; k < height; k++) {
      const spans: Span[] = [{ text: "│", style: "border" }]
      for (let c = 0; c < n; c++) {
        const g = wrapped[c]![k] ?? []
        const gap = widths[c]! - sum(g)
        const left = aligns[c] === "r" ? gap : aligns[c] === "c" ? Math.floor(gap / 2) : 0
        spans.push({ text: " ".repeat(1 + left), style: "text" })
        for (const x of g) {
          const last = spans[spans.length - 1]!
          if (last.style === x.style) last.text += x.ch
          else spans.push({ text: x.ch, style: x.style })
        }
        spans.push({ text: " ".repeat(gap - left + 1), style: "text" }, { text: "│", style: "border" })
      }
      out.push(spans)
    }
    return out
  }
  const bodyRows = body.map(row)
  const ruled = bodyRows.some((r) => r.length > 1)
  const out: Span[][] = [border("┌", "┬", "┐"), ...row(head), border("├", "┼", "┤")]
  bodyRows.forEach((r, i) => {
    if (i > 0 && ruled) out.push(border("├", "┼", "┤"))
    out.push(...r)
  })
  out.push(border("└", "┴", "┘"))
  return out
}

/**
 * A table that could not be drawn, as the model wrote it: its source
 * lines, each wrapped at spaces only, so a number stays on one row (a
 * terminal's word wrap also breaks after "." and "-"). A run of text
 * wider than `width` is cut where the row ends.
 */
export function sourceRows(text: string, width: number): string[] {
  const room = Math.max(1, width)
  const out: string[] = []
  for (const line of text.split("\n")) {
    let row = ""
    for (const word of line.split(/(?<= )/)) {
      if (cols(row + word.trimEnd()) <= room) {
        row += word
        continue
      }
      if (row.trim()) out.push(row.trimEnd())
      row = ""
      let rest = word
      while (cols(rest.trimEnd()) > room) {
        let cut = ""
        for (const { segment } of GRAPHEMES.segment(rest)) {
          if (cols(cut + segment) > room) break
          cut += segment
        }
        out.push(cut)
        rest = rest.slice(cut.length)
      }
      row = rest
    }
    out.push(row.trimEnd())
  }
  return out
}

