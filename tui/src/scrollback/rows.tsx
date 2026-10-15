/**
 * How each transcript item looks in the scrollback.
 *
 * An item is written once (and again after a resize, see transcript.ts), so
 * a renderer is a pure function of the item, the theme and the width.
 * Markdown and diffs are rendered on a surface first, so highlighting can
 * settle before the rows are committed (writers.ts); everything else is a
 * few styled lines.
 */
import { RGBA, TextAttributes, type CliRenderer, type MouseEvent, type Renderable } from "@opentui/core"
import { writeSolidToScrollback } from "@opentui/solid"
import { For, Show, type JSX } from "solid-js"
import { fit, middle, width } from "../text"
import { DOLPHIN, DOLPHIN_WIDTH, LOGO_WIDTH, logoDepth, logoRuns } from "../theme/logo"
import { GUTTER, markdownStyles, type Theme } from "../theme/theme"
import type { DiffRow, ToolResult } from "./tools"
import type { Item } from "./transcript"
import { SETTLE_MS, commitSurface, markdownBlock, type QuietTarget } from "./writers"

type Fields = Record<string, any>

/** A message longer than this shows its first lines and how many more. */
export const USER_MAX_LINES = 12

export function cappedLines(text: string, max: number): string {
  const lines = text.split("\n")
  if (lines.length <= max) return text
  return [...lines.slice(0, max - 1), `… ${lines.length - max + 1} more lines`].join("\n")
}

/** The banner: the dolphin, this much space, then the words (which need
 *  about this many columns beside it, or the dolphin is left out). */
const LOGO_GAP = 3
const BANNER_TEXT_MIN = 32

/** The marker in front of a message, by who sent it. */
const SOURCE_MARK: Record<string, string> = {
  user: "❯ ", queued: "❯ ", steer: "↳ ", harness: "◇ ",
}

export interface RowProps {
  item: Item
  theme: Theme
  width: number
  /** A tool call's details: what it printed, an edit's diff (ctrl+o). */
  expanded?: boolean
  /** A click on a tool call (full screen): its details on or off; the
   *  row clicked, to keep it where it is. */
  onToggle?: (row: Renderable) => void
}

export function Row(props: RowProps): JSX.Element {
  const f = props.item as Fields
  const t = props.theme
  switch (props.item.kind) {
    case "banner": {
      // The dolphin beside what the session is, where both have room;
      // the words alone (under the mark ◆) where they do not.
      // A function: a full-screen row follows the width (rows.tsx is drawn
      // once into the scrollback, but lives on in the window).
      // In colour the logo itself (the disc, half blocks); in mono braille.
      const depth = t.name === "mono" ? "mono" : logoDepth()
      const logoWidth = depth === "mono" ? DOLPHIN_WIDTH : LOGO_WIDTH
      const logo = () => props.width >= logoWidth + LOGO_GAP + BANNER_TEXT_MIN
      const paint = (c: string | number | null) => c === null ? undefined : typeof c === "number" ? RGBA.fromIndex(c) : c
      return (
        <box flexDirection="row" paddingTop={1} paddingBottom={1}>
          <Show when={logo()}>
            <box flexDirection="column" flexShrink={0} width={logoWidth + LOGO_GAP}>
              <Show when={depth !== "mono"} fallback={
                <For each={DOLPHIN}>{(row) => <text fg={t.accent} wrapMode="none">{row}</text>}</For>
              }>
                <For each={logoRuns(depth)}>
                  {(runs) => (
                    <text wrapMode="none">
                      <For each={runs}>{(run) => <span style={{ fg: paint(run.fg), bg: paint(run.bg) }}>{run.text}</span>}</For>
                    </text>
                  )}
                </For>
              </Show>
            </box>
          </Show>
          {/* A basis of 0: the words wrap at what is left, and every wrapped
              row is counted (an auto basis measured them unwrapped). */}
          <box flexDirection="column" flexGrow={1} flexShrink={1} flexBasis={0} paddingTop={logo() ? 1 : 0}>
            <text>
              <span style={{ fg: t.accent, bold: true }}>{logo() ? "DELFIN" : "◆ DELFIN"}</span>
              <span style={{ fg: t.muted }}>{"  " + String(f.head ?? "")}</span>
            </text>
            <For each={(f.lines ?? []) as string[]}>
              {(line) => <text fg={t.faint}>{(logo() ? "" : "  ") + line}</text>}
            </For>
          </box>
        </box>
      )
    }
    case "user": {
      const source = String(f.source ?? "user")
      const tint = source === "harness" ? t.muted : t.bright
      const note = source === "queued" ? "  (queued)" : source === "steer" ? "  (into the turn)" : ""
      // The marker in a column of its own: every further line, and every
      // wrapped one, hangs under the first.
      return (
        <box paddingTop={1} flexDirection="row">
          <text width={2} fg={t.accent} attributes={TextAttributes.BOLD}>{SOURCE_MARK[source] ?? "❯ "}</text>
          <text flexGrow={1} wrapMode="word">
            <span style={{ fg: tint }}>{cappedLines(String(f.text ?? ""), USER_MAX_LINES)}</span>
            <span style={{ fg: t.faint }}>{note}</span>
          </text>
        </box>
      )
    }
    case "reasoning":
      return <text fg={t.faint} wrapMode="none">{"· " + middle(String(f.text ?? ""), props.width - 2)}</text>
    case "tool": {
      const r = (f.result ?? { kind: "text", summary: "", lines: [], ok: true }) as ToolResult
      const failed = r.kind === "blocked" || (r.kind === "text" && !r.ok)
      const name = String(f.name ?? "")
      // A click toggles its details: the left button pressed and let go on
      // the same cell of this row -- not a release that ends a drag begun
      // elsewhere, nor another button. (OpenTUI starts a selection when
      // the button goes down on text, so a release always looks like a
      // drag's end; the view also checks that nothing was selected.)
      let row: Renderable | undefined
      let press: { x: number; y: number; at: number } | null = null
      const down = (event: MouseEvent) => {
        press = event.button === 0 ? { x: event.x, y: event.y, at: Date.now() } : null
      }
      const up = (event: MouseEvent) => {
        const was = press
        press = null
        if (!was || !row || event.button !== 0 || event.x !== was.x || event.y !== was.y) return
        if (Date.now() - was.at < 1500) props.onToggle?.(row)
      }
      return (
        <box ref={(r: Renderable) => { row = r }} flexDirection="column"
             onMouseDown={props.onToggle ? down : undefined} onMouseUp={props.onToggle ? up : undefined}>
          <text wrapMode="none">
            <span style={{ fg: failed ? t.err : t.ok }}>{"⏺ "}</span>
            <span style={{ fg: t.text, bold: true }}>{name}</span>
            <span style={{ fg: t.muted }}>{f.arg ? "  " + middle(String(f.arg), props.width - name.length - 5) : ""}</span>
          </text>
          <ToolResultRows result={r} theme={t} width={props.width} expanded={Boolean(props.expanded)} />
        </box>
      )
    }
    case "denied":
      return (
        <text>
          <span style={{ fg: t.err }}>{"⏺ "}</span>
          <span style={{ fg: t.text, bold: true }}>{String(f.name ?? "")}</span>
          <span style={{ fg: t.err }}>{"  refused"}</span>
        </text>
      )
    case "approval":
      return (
        <text fg={f.ok ? t.muted : t.warn}>
          {(f.ok ? "  ✓ " : "  ✗ ") + String(f.text ?? "")}
        </text>
      )
    case "notice":
      return (
        <box flexDirection="column">
          <For each={String(f.text ?? "").split("\n")}>
            {(line) => <text fg={t.muted}>{line}</text>}
          </For>
        </box>
      )
    case "error":
      return <text fg={t.err} wrapMode="word">{"✗ " + String(f.text ?? "")}</text>
    case "ended":
      return (
        <box paddingTop={1}>
          <text fg={t.muted}>{"◼ " + String(f.text ?? "")}</text>
        </box>
      )
    case "turnEnd":
      return <text fg={f.ok === false ? t.err : f.stopped ? t.warn : t.faint}>{"  " + String(f.text ?? "")}</text>
    case "divider": {
      // A column short of the right edge, as the gutter is at the left;
      // the label in columns (a name in CJK takes two each), cut to fit.
      const span = () => props.width - GUTTER
      const label = () => ` ${fit(String(f.text ?? ""), Math.max(4, span() - 6))} `
      const side = () => Math.max(2, Math.floor((span() - width(label())) / 2))
      return (
        <box paddingTop={1}>
          <text>
            <span style={{ fg: t.accent }}>{"━".repeat(side())}</span>
            <span style={{ fg: t.accentSoft, bold: true }}>{label()}</span>
            <span style={{ fg: t.accent }}>{"━".repeat(Math.max(0, span() - side() - width(label())))}</span>
          </text>
        </box>
      )
    }
    case "lines":
      return (
        <box flexDirection="column" paddingTop={f.title ? 1 : 0}>
          <Show when={f.title}>
            <text fg={t.accentSoft} attributes={TextAttributes.BOLD}>{"  " + String(f.title)}</text>
          </Show>
          <For each={(f.lines ?? []) as string[]}>{(line) => <text fg={t.muted} wrapMode="char">{"  " + line}</text>}</For>
        </box>
      )
    default:
      return <text fg={t.text}>{String(f.text ?? "")}</text>
  }
}

/** Under a tool's headline: what it did -- then, with its details, the
 *  lines that say most (what it printed, the diff). */
function ToolResultRows(props: { result: ToolResult; theme: Theme; width: number; expanded: boolean }): JSX.Element {
  const t = props.theme
  const r = props.result
  const room = () => Math.max(10, props.width - 6)
  if (r.kind === "blocked") {
    return <text fg={t.err} wrapMode="none">{"  ⎿ " + fit(`blocked: ${r.reason}`, props.width - 4)}</text>
  }
  if (r.kind === "edit") {
    const tone = (row: DiffRow) => row.kind === "add" ? t.ok : row.kind === "delete" ? t.err
      : row.kind === "hunk" || row.kind === "file" ? t.accentSoft : t.faint
    const bg = (row: DiffRow) => row.kind === "add" ? t.added : row.kind === "delete" ? t.removed : undefined
    return (
      <box flexDirection="column">
        <text wrapMode="none">
          <span style={{ fg: t.muted }}>{"  ⎿ "}</span>
          <span style={{ fg: t.text, bold: true }}>{`${r.action} ${middle(r.label, Math.max(8, props.width - 18 - r.action.length))}`}</span>
          <span style={{ fg: t.muted }}>{" ("}</span>
          <span style={{ fg: t.ok }}>{`+${r.adds}`}</span>
          <span style={{ fg: t.err }}>{` -${r.dels}`}</span>
          <span style={{ fg: t.muted }}>{")"}</span>
        </text>
        <Show when={props.expanded}>
          <For each={r.rows}>
            {(row) => (
              <text wrapMode="none">
                <span>{"      "}</span>
                <span style={{ fg: tone(row), bg: bg(row) }}>{fit(row.text, room())}</span>
              </text>
            )}
          </For>
          <Show when={r.omitted > 0}>
            <text fg={t.faint}>{`      … ${r.omitted} more diff lines`}</text>
          </Show>
        </Show>
      </box>
    )
  }
  return (
    <box flexDirection="column">
      <text fg={r.ok ? t.muted : t.err} wrapMode="none">{"  ⎿ " + fit(r.summary, props.width - 4)}</text>
      <Show when={props.expanded}>
        <Show when={r.wrap} fallback={
          <For each={r.lines}>{(line) => <text fg={t.faint} wrapMode="none">{"      " + fit(line, room())}</text>}</For>
        }>
          <box flexDirection="column" paddingLeft={6}>
            <For each={r.lines}>{(line) => <text fg={t.faint} wrapMode="word">{line}</text>}</For>
          </box>
        </Show>
      </Show>
    </box>
  )
}

/** The strings an item shows, in display order (what a test reads back). */
export function shownText(item: Item): string {
  const f = item as Fields
  switch (item.kind) {
    case "banner": return [f.head, ...(f.lines ?? [])].join("\n")
    case "tool": {
      const r = f.result as ToolResult | undefined
      const rest = !r ? [] : r.kind === "blocked" ? [`blocked: ${r.reason}`]
        : r.kind === "edit" ? [`${r.action} ${r.label}`, ...r.rows.map((row) => row.text)]
        : [r.summary, ...r.lines]
      return [f.name, f.arg, ...rest].filter(Boolean).join("\n")
    }
    case "lines": return (f.lines as string[]).join("\n")
    default: return String(f.text ?? "")
  }
}

export function itemWriter(renderer: CliRenderer, theme: () => Theme,
                           options: { release?: (item: Item) => void; details?: () => boolean } = {}) {
  return async (item: Item, how?: { first: boolean }): Promise<void> => {
    const f = item as Fields
    const t = theme()
    if (item.kind === "markdown") {
      // A block the footer still shows (held, written the first time):
      // the footer lets go of it in the same step its rows are queued.
      const release = options.release
      const handOff = how?.first && f.held && release
        ? { renderer: renderer as unknown as QuietTarget, release: () => release(item) } : undefined
      await commitSurface(renderer, markdownBlock(String(f.text ?? ""),
                                                  { styles: markdownStyles(t), marginTop: 1, gutter: GUTTER }),
                          SETTLE_MS, handOff)
      return
    }
    // Every row starts a column in (the gutter); its width is what is left.
    const expanded = Boolean(options.details?.())
    writeSolidToScrollback(renderer, (ctx) => (
      <box paddingLeft={GUTTER} flexDirection="column">
        <Row item={item} theme={t} width={ctx.width - GUTTER} expanded={expanded} />
      </box>
    ))
  }
}
