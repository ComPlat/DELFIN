/**
 * A table segment (table.ts) as OpenTUI text: the drawn rows in the
 * markdown style's colours, or -- when the table cannot fit -- its source
 * lines as the model wrote them.
 *
 * Input: the segment, the columns the block has (the right-hand column
 * kept free, as for the rest of the window), the markdown style. Output:
 * the rows a text element shows, already wrapped (it draws them as they
 * are); `blockRenderable` builds a
 * whole block that holds tables (writers.ts), the Solid view does the same
 * for the footer and the full screen (segmented.tsx).
 */
import {
  BoxRenderable,
  MarkdownRenderable,
  StyledText,
  TextAttributes,
  TextRenderable,
  getTreeSitterClient,
  type RenderContext,
  type Renderable,
  type SyntaxStyle,
  type TextChunk,
} from "@opentui/core"
import { drawTable, sourceRows, type Segment, type SpanStyle } from "./table"

/** The style each kind of span takes from the markdown style. */
const STYLE_OF: Record<SpanStyle, string> = {
  border: "conceal", header: "markup.strong", text: "default", strong: "markup.strong", code: "markup.raw",
}

export interface TableView {
  /** The rows, already wrapped: drawn with wrapMode "none". */
  content: StyledText | string
  /** Columns in from the block's left edge (a table under a list item). */
  indent: number
}

/** `segment` at `cols` columns (the block's width less the free column). */
export function tableView(segment: Segment, cols: number, syntax: SyntaxStyle): TableView {
  const drawn = drawTable(segment.text.split("\n"), cols - segment.indent)
  if (!drawn) return { content: sourceRows(segment.text, cols).join("\n"), indent: 0 }
  const chunks: TextChunk[] = []
  drawn.forEach((line, i) => {
    if (i > 0) chunks.push({ __isChunk: true, text: "\n" })
    for (const span of line) {
      const style = syntax.getStyle(STYLE_OF[span.style]) ?? syntax.getStyle("default")
      chunks.push({ __isChunk: true, text: span.text, fg: style?.fg,
                    attributes: style?.bold ? TextAttributes.BOLD : 0 })
    }
  })
  return { content: new StyledText(chunks), indent: segment.indent }
}

export interface BlockOptions {
  syntax: SyntaxStyle
  /** The block's width in columns (its box's, the gutter left out). */
  width: number
  marginTop: number
  paddingLeft: number
}

/** A block cut into segments (at least one a table), as one column of
 *  markdown elements and table texts, spaced as the source was. */
export function blockRenderable(ctx: RenderContext, segments: Segment[], options: BlockOptions): Renderable {
  const root = new BoxRenderable(ctx, { width: "100%", flexDirection: "column", flexShrink: 0,
                                        marginTop: options.marginTop, paddingLeft: options.paddingLeft })
  segments.forEach((segment, i) => {
    const marginTop = i > 0 && segment.gapBefore ? 1 : 0
    if (segment.kind === "markdown") {
      root.add(new MarkdownRenderable(ctx, {
        content: segment.text, syntaxStyle: options.syntax, width: "100%", marginTop, streaming: false,
        conceal: true, tableOptions: { widthMode: "content" }, treeSitterClient: getTreeSitterClient(),
      }))
      return
    }
    const view = tableView(segment, options.width - 1, options.syntax)
    const box = new BoxRenderable(ctx, { width: "100%", flexShrink: 0, marginTop, paddingLeft: view.indent, paddingRight: 1 })
    box.add(new TextRenderable(ctx, { content: view.content, wrapMode: "none", width: "100%", flexShrink: 0 }))
    root.add(box)
  })
  return root
}
