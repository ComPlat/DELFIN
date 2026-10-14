/**
 * The session in view when the window has the whole terminal (the default;
 * `--inline` keeps the terminal's own scrollback instead): its rows in a
 * scroll box that sticks to the bottom, the footer under it (footer.tsx) --
 * the answer's unfinished block there, as inline: a delta then lays out
 * that block, not every row of the conversation.
 *
 * The rows are the scrollback's (rows.tsx), so both screens look alike; a
 * markdown block is a live markdown element here and re-wraps by itself
 * when the width changes (no replay). The wheel scrolls it, PgUp/PgDn and
 * ctrl+home/end too (app.tsx). Scrolled up, it stays where it is while the
 * session goes on; scrolled back to the end, it follows again.
 */
import { SyntaxStyle, type CliRenderer, type Renderable, type ScrollBoxRenderable } from "@opentui/core"
import { For, Show, createEffect, createMemo, on, onCleanup, onMount, type Accessor, type JSX } from "solid-js"
import { createStore, produce } from "solid-js/store"
import { markdownStyles, type Theme } from "../theme/theme"
import { HIGHLIGHT_WAIT_MS, whenHighlighted } from "./highlighter"
import { Row } from "./rows"
import type { Item } from "./transcript"

/** Rows kept in view: older ones leave, so a long session costs what a
 *  few hundred rows do -- laid out again when one comes, held in memory --
 *  not all of it (a note says how many left; `delfin-agent report` has
 *  the whole session). */
export const MAX_ITEMS = 600

/** How long a held block is kept out of sight at most (its highlighting
 *  gives up sooner, and the block is shown plain). */
export const HOLD_MS = HIGHLIGHT_WAIT_MS + 500

export interface ScreenItems {
  readonly items: readonly Item[]
  /** Rows let go from the top to stay within the limit. */
  readonly dropped: number
  /** Counts the clears (a session put in view): the view follows again. */
  readonly cleared: number
  /** A block the footer still shows (`held` set: it streamed there), not
   *  shown here yet: it takes no rows until its highlighting is in (a
   *  block drawn before that shows nothing for a frame). */
  held(id: string): boolean
  /** Resolves once the item is shown: a held block once `reveal`d, or
   *  after HOLD_MS; any other row at once (a block that did not stream --
   *  a replay -- has no copy to hand over: it is drawn as it comes). */
  push(item: Item): Promise<void>
  reveal(id: string): void
  clear(): void
}

export function createScreenItems(max = MAX_ITEMS, holdMs = HOLD_MS): ScreenItems {
  const [state, setState] = createStore<{ items: Item[]; dropped: number; cleared: number; held: Record<string, boolean> }>(
    { items: [], dropped: 0, cleared: 0, held: {} })
  const waiting = new Map<string, { done: () => void; timer: ReturnType<typeof setTimeout> }>()
  const reveal = (id: string) => {
    const wait = waiting.get(id)
    if (!wait) return
    waiting.delete(id)
    clearTimeout(wait.timer)
    setState("held", id, undefined!)
    wait.done()
  }
  return {
    get items() { return state.items },
    get dropped() { return state.dropped },
    get cleared() { return state.cleared },
    held: (id) => Boolean(state.held[id]),
    push(item) {
      const hold = item.kind === "markdown" && Boolean(item.held)
      const shown = hold
        ? new Promise<void>((done) => {
            const timer = setTimeout(() => reveal(item.id), holdMs)
            ;(timer as { unref?: () => void }).unref?.()
            waiting.set(item.id, { done, timer })
          })
        : Promise.resolve()
      setState(produce((s) => {
        if (hold) s.held[item.id] = true
        s.items.push(item)
        if (s.items.length > max) {
          const over = s.items.length - max
          for (const gone of s.items.splice(0, over)) delete s.held[gone.id]
          s.dropped += over
        }
      }))
      return shown
    },
    reveal,
    clear() {
      for (const id of [...waiting.keys()]) reveal(id)
      setState((s) => ({ items: [], dropped: 0, cleared: s.cleared + 1, held: {} }))
    },
  }
}


export interface TranscriptViewProps {
  renderer: CliRenderer
  screen: ScreenItems
  theme: () => Theme
  /** The width it draws in (the window's, less the gutter): a row is
   *  cut to it, less the scroll bar. */
  width: Accessor<number>
  /** Shown (hidden while the agent view has the screen). */
  visible: Accessor<boolean>
  ref?: (box: ScrollBoxRenderable) => void
  /** The wheel turned over it (the app repaints once it rests). */
  onUserScroll?: () => void
}

export function TranscriptView(props: TranscriptViewProps): JSX.Element {
  const t = props.theme
  const style = themedStyle(props.renderer, t)
  const inner = () => Math.max(10, props.width() - 1)
  // A quiet scroll bar: its thumb in the faint tone, no track (made once:
  // a prop object in the JSX would be made again on every read).
  const bar = createMemo(() => ({ trackOptions: { foregroundColor: t().faint, backgroundColor: "transparent" } }))
  let box: ScrollBoxRenderable | undefined
  // Another session in view: it is followed from its end, also when the
  // reader had scrolled up in the one before (OpenTUI keeps a manual
  // scroll until something scrolls -- a list that fits never does).
  createEffect(on(() => props.screen.cleared, () => {
    if (!box || box.isDestroyed) return
    box.stickyScroll = false
    box.stickyScroll = true
    box.scrollTo(box.scrollHeight)
  }, { defer: true }))
  return (
    <scrollbox ref={(b: ScrollBoxRenderable) => { box = b; props.ref?.(b) }} flexGrow={1} flexShrink={1} stickyScroll stickyStart="bottom" viewportCulling
               verticalScrollbarOptions={bar()} visible={props.visible()} focusable={false}
               onMouseScroll={() => props.onUserScroll?.()}>
      <Show when={props.screen.dropped > 0}>
        <text fg={t().faint}>
          {`(${props.screen.dropped} earlier rows are not kept in the window; \`delfin-agent report\` has them all)`}
        </text>
      </Show>
      <For each={props.screen.items}>
        {(item) => item.kind === "markdown"
          ? <Block item={item} screen={props.screen} syntax={style} />
          // A row takes its colours when it is made: made again when the
          // theme changes (the terminal said light after the first frame).
          : <Show when={t()} keyed>{(theme: Theme) => <Row item={item} theme={theme} width={inner()} />}</Show>}
      </For>
    </scrollbox>
  )
}

/** A finished markdown block. Held (the footer still shows it), it is in
 *  a box of no rows until it is highlighted, then shown -- the footer
 *  lets go of it in the same frame. Its highlighting is started at once
 *  in any case, so a block out of sight is ready when it comes into view. */
function Block(props: { item: Item; screen: ScreenItems; syntax: Accessor<SyntaxStyle> }): JSX.Element {
  const held = () => props.screen.held(props.item.id)
  let markdown: Renderable | undefined
  onMount(() => {
    void whenHighlighted(markdown!).then(() => props.screen.reveal(props.item.id))
  })
  return (
    <box flexShrink={0} height={held() ? 0 : "auto"} overflow="hidden" marginTop={held() ? 0 : 1}>
      {/* Content last: the element builds its blocks when it is set, with
          the style it has by then (props go in in order). */}
      <markdown ref={(m: Renderable) => { markdown = m }} syntaxStyle={props.syntax()} conceal width="100%"
                tableOptions={{ widthMode: "content" }} content={String(props.item.text ?? "")} />
    </box>
  )
}

/** The markdown style for the theme, made again when the theme changes.
 *  Every one made is freed with the view, not before: a block out of view
 *  takes the new one only when it is next drawn, and may draw with the
 *  old one until then (a theme changes once or twice in a window's life). */
export function themedStyle(renderer: CliRenderer, theme: () => Theme): Accessor<SyntaxStyle> {
  const made: SyntaxStyle[] = []
  const style = createMemo(() => {
    const next = SyntaxStyle.fromStyles(markdownStyles(theme()), renderer.nativeScene!)
    made.push(next)
    return next
  })
  onCleanup(() => { for (const one of made) one.destroy() })
  return style
}
