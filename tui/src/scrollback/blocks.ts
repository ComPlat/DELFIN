/**
 * Splitting a streaming markdown answer into blocks that are final.
 *
 * Committed scrollback cannot be changed, so an answer is committed one
 * top-level block at a time -- a heading, a paragraph, a list, a table, a
 * code block -- each as its own immutable snapshot, once it is complete.
 * Every block but the last is complete: the lexer only starts a new block
 * after the previous one ended. The last block stays in the footer, live,
 * until the next one begins or the answer ends.
 *
 * Committing whole blocks rather than row ranges of one growing surface is
 * what lets a resize re-wrap a session mid-turn: every committed piece is
 * an item that can be rendered again at the new width.
 *
 * Uses `marked`, the lexer OpenTUI's own markdown renderer uses (the version
 * it pins), so the blocks are the blocks it would draw.
 */
import { marked, type Token } from "marked"

export interface MarkdownSplit {
  /** Complete top-level blocks, in order, without trailing blank lines. */
  blocks: string[]
  /** The unfinished rest, exactly as it arrived. Empty once final. */
  tail: string
  /** How much of the text the blocks cover (where the tail starts). */
  consumed: number
}

/** A last line that may still become a list item ("2", "2.", "-"): the
 *  list before it is not complete yet -- one more character can make
 *  them one list. */
const MAYBE_ITEM = /^\s{0,3}(\d{1,9}[.)]?|[-*+])\s*$/

function lex(text: string): Token[] {
  try {
    return marked.lexer(text)
  } catch {
    // The lexer has never been seen to throw; if it does, keep the whole
    // text in the tail rather than commit something wrong.
    return []
  }
}

/** Where a table token's raw text stops being the table: the start of its
 *  first body line with no unescaped pipe, or -1 when every line has one.
 *  The table extension takes such a line as one more row; a model means
 *  the sentence after the table. The dashboard and the plain terminal end
 *  a table at the same line (delfin/agent/markdown_table.py). */
export function tableEnd(raw: string): number {
  const lines = raw.split("\n")
  let at = 0
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]!
    if (i >= 2 && line.trim() && !/(^|[^\\])\|/.test(line)) return at
    at += line.length + 1
  }
  return -1
}

/** Split `text`; `final` commits the last block too (the answer ended). */
export function splitMarkdown(text: string, final: boolean): MarkdownSplit {
  const normal = text.replace(/\r\n?/g, "\n")
  const tokens = lex(normal)
  if (tokens.length === 0) {
    return final ? { blocks: normal.trim() ? [normal.replace(/\n+$/, "")] : [], tail: "", consumed: normal.length }
                 : { blocks: [], tail: normal, consumed: 0 }
  }
  let lastBlock = tokens.length - 1
  while (lastBlock > 0 && tokens[lastBlock]!.type === "space") lastBlock--
  let stableEnd = final ? tokens.length : lastBlock
  if (!final && MAYBE_ITEM.test(tokens[lastBlock]!.raw)) {
    // Hold the list (and the space after it) back with the would-be item.
    let prev = lastBlock - 1
    while (prev >= 0 && tokens[prev]!.type === "space") prev--
    if (prev >= 0 && tokens[prev]!.type === "list") stableEnd = prev
  }
  const blocks: string[] = []
  let offset = 0
  for (let i = 0; i < stableEnd; i++) {
    const token = tokens[i]!
    const cut = token.type === "table" ? tableEnd(token.raw) : -1
    if (cut > 0) {
      // The table, then the rest split again from its first pipe-less line.
      blocks.push(token.raw.slice(0, cut).replace(/\n+$/, ""))
      const start = offset + cut
      const rest = splitMarkdown(normal.slice(start), final)
      return { blocks: blocks.concat(rest.blocks), tail: rest.tail, consumed: start + rest.consumed }
    }
    offset += token.raw.length
    if (token.type === "space") continue
    const raw = token.raw.replace(/\n+$/, "")
    if (raw.trim()) blocks.push(raw)
  }
  return { blocks, tail: final ? "" : normal.slice(offset), consumed: final ? normal.length : offset }
}

/**
 * Feeds an answer delta by delta and calls `commit` once for each block as
 * it becomes final, and `onTail` with whatever is still unfinished (and
 * the number that block will have when it is committed).
 *
 * What is committed is a prefix of the text, and only the text after it
 * is split again: a later delta can change how earlier text lexes (a
 * loose list going on after a blank line merges into one list), and a
 * block committed is never taken back -- so what follows it must never be
 * counted as already committed either. At worst one list arrives as two.
 */
export class AnswerBlocks {
  private text = ""
  private offset = 0                  // characters of `text` committed
  private n = 0
  private cr = ""                     // a CR whose LF may still come

  constructor(
    private readonly commit: (block: string, index: number) => void,
    private readonly onTail: (tail: string, index: number) => void = () => {},
  ) {}

  get content(): string {
    return this.text
  }

  /** The part of the text already in the scrollback. */
  get committedText(): string {
    return this.text.slice(0, this.offset)
  }

  append(delta: string): void {
    if (!delta) return
    // A CR at the end of a delta may be the first half of a CRLF the next
    // delta completes: held back until it is known which.
    let piece = this.cr + delta
    this.cr = ""
    if (piece.endsWith("\r")) {
      this.cr = "\r"
      piece = piece.slice(0, -1)
    }
    this.text += piece.replace(/\r\n?/g, "\n")
    this.flush(false)
  }

  finish(): void {
    if (this.cr) this.text += "\n"
    this.cr = ""
    this.flush(true)
  }

  /**
   * The part ended with `full` as its whole text (text.ended). When what
   * streamed is the start of it, the rest is committed; when deltas went
   * missing (a reconnect mid-answer), what was committed stays and the
   * rest of `full` after it is committed -- or, if even the committed
   * part does not match, `full` is committed whole after a note.
   */
  finishWith(full: string, onMismatch: () => void = () => {}): void {
    this.cr = ""
    const text = full.replace(/\r\n?/g, "\n")
    if (!text.startsWith(this.committedText)) {
      onMismatch()
      this.text = text
      this.offset = 0
    } else {
      this.text = text
    }
    this.flush(true)
  }

  private flush(final: boolean): void {
    const rest = this.text.slice(this.offset)
    const { blocks, tail, consumed } = splitMarkdown(rest, final)
    for (const block of blocks) this.commit(block, this.n++)
    this.offset += consumed
    this.onTail(tail, this.n)
  }
}
