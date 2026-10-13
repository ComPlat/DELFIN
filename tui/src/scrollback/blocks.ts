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
}

function lex(text: string): Token[] {
  try {
    return marked.lexer(text)
  } catch {
    // The lexer has never been seen to throw; if it does, keep the whole
    // text in the tail rather than commit something wrong.
    return []
  }
}

/** Split `text`; `final` commits the last block too (the answer ended). */
export function splitMarkdown(text: string, final: boolean): MarkdownSplit {
  const normal = text.replace(/\r\n?/g, "\n")
  const tokens = lex(normal)
  if (tokens.length === 0) {
    return final ? { blocks: normal.trim() ? [normal.replace(/\n+$/, "")] : [], tail: "" }
                 : { blocks: [], tail: normal }
  }
  let lastBlock = tokens.length - 1
  while (lastBlock > 0 && tokens[lastBlock]!.type === "space") lastBlock--
  const stableEnd = final ? tokens.length : lastBlock
  const blocks: string[] = []
  let offset = 0
  for (let i = 0; i < stableEnd; i++) {
    const token = tokens[i]!
    offset += token.raw.length
    if (token.type === "space") continue
    const raw = token.raw.replace(/\n+$/, "")
    if (raw.trim()) blocks.push(raw)
  }
  return { blocks, tail: final ? "" : normal.slice(offset) }
}

/**
 * Feeds an answer delta by delta and calls `commit` once for each block as
 * it becomes final, and `onTail` with whatever is still unfinished.
 */
export class AnswerBlocks {
  private text = ""
  private committed = 0

  constructor(
    private readonly commit: (block: string, index: number) => void,
    private readonly onTail: (tail: string) => void = () => {},
  ) {}

  get content(): string {
    return this.text
  }

  append(delta: string): void {
    if (!delta) return
    this.text += delta
    this.flush(false)
  }

  finish(): void {
    this.flush(true)
  }

  private flush(final: boolean): void {
    const { blocks, tail } = splitMarkdown(this.text, final)
    // A later delta can merge what looked like two blocks (a loose list
    // continuing after a blank line). What is committed stays committed;
    // only blocks past that point are new.
    for (let i = this.committed; i < blocks.length; i++) this.commit(blocks[i]!, i)
    this.committed = Math.max(this.committed, blocks.length)
    this.onTail(tail)
  }
}
