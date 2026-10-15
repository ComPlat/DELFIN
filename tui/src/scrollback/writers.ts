/**
 * Writing one item into the scrollback: render it off-screen at the current
 * width, wait for syntax highlighting (bounded), commit every row at once.
 *
 * Highlighting runs in a tree-sitter worker and arrives asynchronously; a
 * snapshot taken before it lands is committed unhighlighted and can never
 * be repainted. ``commitSurface`` therefore renders on a ScrollbackSurface,
 * awaits ``settle()`` -- but never longer than the bound, because a worker
 * that does not answer must cost colour, not the rest of the transcript.
 */
import {
  BoxRenderable,
  DiffRenderable,
  MarkdownRenderable,
  SyntaxStyle,
  getTreeSitterClient,
  type RenderContext,
  type Renderable,
  type ScrollbackSurface,
} from "@opentui/core"
import { parsePatch } from "diff"
import { drawPlain } from "./highlighter"
import { bounded } from "./queue"
import { tableSegments } from "./table"
import { blockRenderable } from "./tableview"

export interface SurfaceTarget {
  createScrollbackSurface(options?: { startOnNewLine?: boolean }): ScrollbackSurface
}

/** What tells whether the renderer is between frames with nothing waiting
 *  to be written (OpenTUI keeps both to itself). */
export interface QuietTarget {
  readonly rendering?: boolean
  readonly externalOutputQueue?: { readonly size: number }
  /** "explicit_suspended" while the terminal is handed over (/editor):
   *  no frame comes then, and nothing is waited for. */
  readonly controlState?: string
}

/** How long a hand-off waits for the renderer to be quiet before it goes
 *  anyway (no frame may come at all: the window is suspended or closing). */
const QUIET_WAIT_MS = 60

/**
 * Resolves when no frame is being drawn and nothing waits to be written.
 * A change of the footer's height made then is applied at once; made
 * while a frame is drawn or rows wait to go out, OpenTUI applies it a
 * frame later -- the footer would have let go of a block's text and kept
 * its rows for a frame, blank above the block's rows.
 */
export async function whenQuiet(target: QuietTarget, ms = QUIET_WAIT_MS): Promise<void> {
  const until = Date.now() + ms
  while ((target.rendering || (target.externalOutputQueue?.size ?? 0) > 0) && Date.now() < until
         && target.controlState !== "explicit_suspended") {
    await new Promise<void>((resolve) => setImmediate(resolve))
  }
}

/** A footer that shows the text of a block being committed: it lets go of
 *  it right before the block's rows are queued, in the same step. */
export interface HandOff {
  renderer: QuietTarget
  release: () => void
}

export interface Built {
  root: Renderable
  dispose?: () => void
}

/** How long a commit waits for highlighting before it goes without. */
export const SETTLE_MS = 1500

export async function commitSurface(
  target: SurfaceTarget,
  build: (ctx: RenderContext) => Built,
  settleMs: number = SETTLE_MS,
  handOff?: HandOff,
): Promise<void> {
  const surface = target.createScrollbackSurface({ startOnNewLine: true })
  let built: Built | undefined
  try {
    built = build(surface.renderContext)
    surface.root.add(built.root)
    surface.render()
    try {
      await bounded(surface.settle(settleMs), settleMs + 250)
    } catch {
      // The highlighter did not answer in time (a busy node, a large
      // code block): the rows go without colour -- never without text,
      // which is what a block drawn before its highlighting shows.
      drawPlain(built.root)
      surface.render()
    }
    if (handOff) {
      await whenQuiet(handOff.renderer)
      handOff.release()
    }
    if (surface.isDestroyed) return
    try {
      if (surface.height > 0) surface.commitRows(0, surface.height)
    } catch {
      // The terminal changed size meanwhile: drawn again at the new one.
      surface.render()
      if (surface.height > 0) surface.commitRows(0, surface.height)
    }
  } finally {
    surface.destroy()
    built?.dispose?.()
  }
}

/** True when ``text`` is a unified diff OpenTUI's diff view can draw. */
export function isDrawableDiff(text: string): boolean {
  try {
    const patches = parsePatch(text)
    return patches.length > 0 && patches.every((p) => p.hunks.length > 0)
  } catch {
    return false
  }
}

export type StyleMap = Parameters<typeof SyntaxStyle.fromStyles>[0]

export interface MarkdownOptions {
  styles: StyleMap
  marginTop?: number
  /** Columns kept free at the left (the theme's GUTTER). */
  gutter?: number
}

/** One complete markdown block (see blocks.ts), as a committable surface.
 *  A block that holds a table is drawn in segments: its tables by
 *  table.ts, the rest by OpenTUI's markdown element (tableview.ts). */
export function markdownBlock(raw: string, options: MarkdownOptions) {
  return (ctx: RenderContext): Built => {
    const style = SyntaxStyle.fromStyles(options.styles, ctx.nativeScene!)
    const segments = tableSegments(raw)
    if (segments.some((s) => s.kind === "table")) {
      const gutter = options.gutter ?? 0
      const root = blockRenderable(ctx, segments, { syntax: style, width: ctx.width - gutter,
                                                    marginTop: options.marginTop ?? 0, paddingLeft: gutter })
      return { root, dispose: () => style.destroy() }
    }
    const markdown = new MarkdownRenderable(ctx, {
      content: raw,
      syntaxStyle: style,
      width: "100%",
      marginTop: options.marginTop ?? 0,
      streaming: false,
      conceal: true,
      tableOptions: { widthMode: "content" },
      treeSitterClient: getTreeSitterClient(),
    })
    if (!options.gutter) return { root: markdown, dispose: () => style.destroy() }
    // In a box that keeps the gutter free: the markdown wraps at what is
    // left, as the live tail in the footer does.
    const root = new BoxRenderable(ctx, { width: "100%", flexDirection: "column", paddingLeft: options.gutter })
    root.add(markdown)
    return { root, dispose: () => style.destroy() }
  }
}

export interface DiffOptions {
  styles: StyleMap
  filetype: string
  colors: {
    text: string
    lineNumber: string
    added: string
    removed: string
    addedSign: string
    removedSign: string
  }
}

/** A validated unified diff, with line numbers and highlighting. */
export function diffBlock(diff: string, options: DiffOptions) {
  return (ctx: RenderContext): Built => {
    const style = SyntaxStyle.fromStyles(options.styles, ctx.nativeScene!)
    const c = options.colors
    const root = new DiffRenderable(ctx, {
      diff,
      view: "unified",
      filetype: options.filetype,
      syntaxStyle: style,
      showLineNumbers: true,
      wrapMode: "char",
      width: "100%",
      fg: c.text,
      lineNumberFg: c.lineNumber,
      addedBg: c.added,
      removedBg: c.removed,
      addedContentBg: c.added,
      removedContentBg: c.removed,
      addedSignColor: c.addedSign,
      removedSignColor: c.removedSign,
      treeSitterClient: getTreeSitterClient(),
    })
    return { root, dispose: () => style.destroy() }
  }
}
