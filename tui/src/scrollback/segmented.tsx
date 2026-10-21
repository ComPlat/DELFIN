/**
 * A markdown block that holds a table, in the footer and the full screen:
 * its markdown segments as markdown elements, its tables as text laid out
 * by table.ts at the width there is (tableview.ts). A block without a
 * table never comes here -- the sites draw it as one markdown element.
 */
import type { Renderable, SyntaxStyle, TextRenderable } from "@opentui/core"
import { For, createEffect, createMemo, type JSX } from "solid-js"
import { tableSegments, type Segment } from "./table"
import { tableView } from "./tableview"

/** Whether `text` holds a table this view draws. */
export function holdsTable(text: string): boolean {
  return tableSegments(text).some((s) => s.kind === "table")
}

export function SegmentedMarkdown(props: {
  text: () => string
  /** The columns the block has. */
  width: () => number
  syntax: () => SyntaxStyle
  marginTop?: () => number
  ref?: (r: Renderable) => void
}): JSX.Element {
  const segments = createMemo(() => tableSegments(props.text()))
  const top = (segment: Segment, i: number) => i > 0 && segment.gapBefore ? 1 : 0
  return (
    <box ref={props.ref} width="100%" flexDirection="column" flexShrink={0} marginTop={props.marginTop?.() ?? 0}>
      <For each={segments()}>
        {(segment, i) => {
          if (segment.kind === "markdown") {
            return <markdown syntaxStyle={props.syntax()} conceal width="100%" flexShrink={0} marginTop={top(segment, i())}
                             tableOptions={{ widthMode: "content" }} content={segment.text} />
          }
          const view = createMemo(() => tableView(segment, props.width() - 1, props.syntax()))
          // Styled text goes in through the element: as a prop it would be
          // turned into a string.
          let node: TextRenderable | undefined
          createEffect(() => {
            const content = view().content
            if (node && !node.isDestroyed) node.content = content
          })
          return (
            <box width="100%" flexShrink={0} marginTop={top(segment, i())} paddingLeft={view().indent} paddingRight={1}>
              <text ref={(r: TextRenderable) => { node = r }} width="100%" flexShrink={0} wrapMode="none" />
            </box>
          )
        }}
      </For>
    </box>
  )
}
