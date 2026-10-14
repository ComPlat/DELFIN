/**
 * How each transcript item looks in the scrollback.
 *
 * An item is written once (and again after a resize, see transcript.ts), so
 * a renderer is a pure function of the item, the theme and the width.
 * Markdown and diffs are rendered on a surface first, so highlighting can
 * settle before the rows are committed (writers.ts); everything else is a
 * few styled lines.
 */
import { TextAttributes, type CliRenderer } from "@opentui/core"
import { writeSolidToScrollback } from "@opentui/solid"
import { For, Show, type JSX } from "solid-js"
import { fit, middle } from "../text"
import { markdownStyles, type Theme } from "../theme/theme"
import type { DiffRow, ToolResult } from "./tools"
import type { Item } from "./transcript"
import { commitSurface, markdownBlock } from "./writers"

type Fields = Record<string, any>

/** A message longer than this shows its first lines and how many more. */
export const USER_MAX_LINES = 12

export function cappedLines(text: string, max: number): string {
  const lines = text.split("\n")
  if (lines.length <= max) return text
  return [...lines.slice(0, max - 1), `… ${lines.length - max + 1} more lines`].join("\n")
}

/** The marker in front of a message, by who sent it. */
const SOURCE_MARK: Record<string, string> = {
  user: "❯ ", queued: "❯ ", steer: "↳ ", harness: "◇ ",
}

function Row(props: { item: Item; theme: Theme; width: number }): JSX.Element {
  const f = props.item as Fields
  const t = props.theme
  switch (props.item.kind) {
    case "banner":
      return (
        <box flexDirection="column" paddingTop={1} paddingBottom={1}>
          <text>
            <span style={{ fg: t.accent, attributes: TextAttributes.BOLD }}>{"◆ DELFIN"}</span>
            <span style={{ fg: t.muted }}>{"  " + String(f.head ?? "")}</span>
          </text>
          <For each={(f.lines ?? []) as string[]}>
            {(line) => <text fg={t.faint}>{"  " + line}</text>}
          </For>
        </box>
      )
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
      return (
        <box flexDirection="column">
          <text wrapMode="none">
            <span style={{ fg: failed ? t.err : t.ok }}>{"⏺ "}</span>
            <span style={{ fg: t.text, attributes: TextAttributes.BOLD }}>{name}</span>
            <span style={{ fg: t.muted }}>{f.arg ? "  " + middle(String(f.arg), props.width - name.length - 5) : ""}</span>
          </text>
          <ToolResultRows result={r} theme={t} width={props.width} />
        </box>
      )
    }
    case "denied":
      return (
        <text>
          <span style={{ fg: t.err }}>{"⏺ "}</span>
          <span style={{ fg: t.text, attributes: TextAttributes.BOLD }}>{String(f.name ?? "")}</span>
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
      const label = ` ${String(f.text ?? "")} `
      const side = Math.max(2, Math.floor((props.width - label.length) / 2))
      return (
        <box paddingTop={1}>
          <text>
            <span style={{ fg: t.accent }}>{"━".repeat(side)}</span>
            <span style={{ fg: t.accentSoft, attributes: TextAttributes.BOLD }}>{label}</span>
            <span style={{ fg: t.accent }}>{"━".repeat(Math.max(0, props.width - side - label.length))}</span>
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

/** Under a tool's headline: what it did, then the lines that say most. */
function ToolResultRows(props: { result: ToolResult; theme: Theme; width: number }): JSX.Element {
  const t = props.theme
  const r = props.result
  const room = Math.max(10, props.width - 6)
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
          <span style={{ fg: t.text, attributes: TextAttributes.BOLD }}>{`${r.action} ${middle(r.label, Math.max(8, props.width - 18 - r.action.length))}`}</span>
          <span style={{ fg: t.muted }}>{" ("}</span>
          <span style={{ fg: t.ok }}>{`+${r.adds}`}</span>
          <span style={{ fg: t.err }}>{` -${r.dels}`}</span>
          <span style={{ fg: t.muted }}>{")"}</span>
        </text>
        <For each={r.rows}>
          {(row) => (
            <text wrapMode="none">
              <span>{"      "}</span>
              <span style={{ fg: tone(row), bg: bg(row) }}>{fit(row.text, room)}</span>
            </text>
          )}
        </For>
        <Show when={r.omitted > 0}>
          <text fg={t.faint}>{`      … ${r.omitted} more diff lines (ctrl+o shows the result)`}</text>
        </Show>
      </box>
    )
  }
  return (
    <box flexDirection="column">
      <text fg={r.ok ? t.muted : t.err} wrapMode="none">{"  ⎿ " + fit(r.summary, props.width - 4)}</text>
      <For each={r.lines}>{(line) => <text fg={t.faint} wrapMode="none">{"      " + fit(line, room)}</text>}</For>
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

export function itemWriter(renderer: CliRenderer, theme: () => Theme) {
  return async (item: Item): Promise<void> => {
    const f = item as Fields
    const t = theme()
    if (item.kind === "markdown") {
      await commitSurface(renderer, markdownBlock(String(f.text ?? ""), { styles: markdownStyles(t), marginTop: 1 }))
      return
    }
    writeSolidToScrollback(renderer, (ctx) => <Row item={item} theme={t} width={ctx.width} />)
  }
}
