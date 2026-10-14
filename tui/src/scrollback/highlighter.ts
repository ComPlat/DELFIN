/**
 * Waiting for the syntax highlighting of a markdown block -- OpenTUI's
 * tree-sitter client, in a worker, answering a few ms after a block is
 * made (and much later on a busy machine) -- and drawing it plain when
 * the answer does not come.
 */
import { CodeRenderable, MarkdownRenderable, type Renderable } from "@opentui/core"
import { bounded } from "./queue"

/** How long a block waits for its highlighting at most: a highlighter that
 *  does not answer costs colour, not the text. */
export const HIGHLIGHT_WAIT_MS = 1000

function codeBlocks(root: Renderable): CodeRenderable[] {
  const found: CodeRenderable[] = []
  const visit = (r: Renderable) => {
    if (r.isDestroyed) return
    if (r instanceof CodeRenderable) found.push(r)
    for (const child of r.getChildren() as Renderable[]) visit(child)
  }
  visit(root)
  return found
}

/** Starts what a frame's update pass would (`onUpdate` is that hook; it
 *  only starts what is due). */
function start(code: CodeRenderable): void {
  ;(code as unknown as { onUpdate(dt: number): void }).onUpdate(0)
}

/**
 * Resolves once every code block under `root` is highlighted, or after
 * `ms` -- and then what is not is drawn as plain text (`drawPlain`). The
 * highlighting is started here, at once: OpenTUI starts it in a frame's
 * update pass, which leaves out what a scroll box culls (rows below the
 * view while the reader is scrolled up) and what is hidden (the agent
 * view has the screen) -- such a block, shown later, would be blank for a
 * frame. A markdown element's style reaches its blocks in that pass too:
 * it is passed on first (`refreshStyles`), or they would be highlighted
 * without one -- as plain text, for good.
 */
export async function whenHighlighted(root: Renderable, ms = HIGHLIGHT_WAIT_MS): Promise<void> {
  const until = Date.now() + ms
  const styled = new Set<Renderable>()
  for (;;) {
    const pending: Promise<unknown>[] = []
    // A markdown element builds its blocks when its content is set, which
    // may come after the element is mounted: not built yet, not started.
    const visit = (r: Renderable) => {
      if (r.isDestroyed) return
      if (r instanceof MarkdownRenderable && !styled.has(r)) {
        if (r.content && r.getChildren().length === 0) pending.push(Bun.sleep(1))
        else {
          r.refreshStyles()
          styled.add(r)
        }
      }
      for (const child of r.getChildren() as Renderable[]) visit(child)
    }
    visit(root)
    for (const code of codeBlocks(root)) {
      if (!code.syntaxStyle) pending.push(Bun.sleep(16))         // not styled yet: not started
      else {
        start(code)
        if (code.isHighlighting) pending.push(code.highlightingDone)
      }
    }
    const left = until - Date.now()
    if (pending.length === 0) return
    if (left <= 0) {
      drawPlain(root)
      return
    }
    await bounded(Promise.all(pending), left)
  }
}

/** Every code block under `root` that is still waiting for its colours
 *  draws its text plain meanwhile (a block drawn before its highlighting
 *  shows nothing). */
export function drawPlain(root: Renderable): void {
  for (const code of codeBlocks(root)) {
    if (!code.isHighlighting || code.drawUnstyledText) continue
    code.drawUnstyledText = true
    start(code)
  }
}
