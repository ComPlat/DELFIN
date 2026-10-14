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
import { bounded } from "./queue"

export interface SurfaceTarget {
  createScrollbackSurface(options?: { startOnNewLine?: boolean }): ScrollbackSurface
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
): Promise<void> {
  const surface = target.createScrollbackSurface({ startOnNewLine: true })
  let built: Built | undefined
  try {
    built = build(surface.renderContext)
    surface.root.add(built.root)
    surface.render()
    await bounded(surface.settle(settleMs), settleMs + 250)
    if (surface.height > 0) surface.commitRows(0, surface.height)
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

/** One complete markdown block (see blocks.ts), as a committable surface. */
export function markdownBlock(raw: string, options: MarkdownOptions) {
  return (ctx: RenderContext): Built => {
    const style = SyntaxStyle.fromStyles(options.styles, ctx.nativeScene!)
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
