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
import { For, type JSX } from "solid-js"
import { diffColors, markdownStyles, type Theme } from "../theme/theme"
import type { Item } from "./transcript"
import { commitSurface, diffBlock, isDrawableDiff, markdownBlock } from "./writers"

type Fields = Record<string, any>

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
      return (
        <box paddingTop={1}>
          <text wrapMode="word">
            <span style={{ fg: t.accent, attributes: TextAttributes.BOLD }}>{SOURCE_MARK[source] ?? "❯ "}</span>
            <span style={{ fg: tint }}>{String(f.text ?? "")}</span>
            <span style={{ fg: t.faint }}>{note}</span>
          </text>
        </box>
      )
    }
    case "reasoning":
      return <text fg={t.faint} wrapMode="word">{"  ∴ " + String(f.text ?? "")}</text>
    case "tool":
      return (
        <box flexDirection="column">
          <text wrapMode="char">
            <span style={{ fg: f.ok === false ? t.err : t.ok }}>{"⏺ "}</span>
            <span style={{ fg: t.text, attributes: TextAttributes.BOLD }}>{String(f.name ?? "")}</span>
            <span style={{ fg: t.muted }}>{f.arg ? "  " + String(f.arg) : ""}</span>
          </text>
          <text fg={t.muted} wrapMode="char">{"  ⎿  " + String(f.summary ?? "")}</text>
          <For each={(f.tail ?? []) as string[]}>{(line) => <text fg={t.faint}>{"     " + line}</text>}</For>
        </box>
      )
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
    case "turnEnd":
      return <text fg={f.ok === false ? t.err : t.faint}>{"  " + String(f.text ?? "")}</text>
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
        <box flexDirection="column">
          <For each={(f.lines ?? []) as string[]}>{(line) => <text fg={t.muted}>{line}</text>}</For>
        </box>
      )
    default:
      return <text fg={t.text}>{String(f.text ?? "")}</text>
  }
}

/** The strings an item shows, in display order (what a test reads back). */
export function shownText(item: Item): string {
  const f = item as Fields
  switch (item.kind) {
    case "banner": return [f.head, ...(f.lines ?? [])].join("\n")
    case "tool": return [f.name, f.arg, f.summary, ...(f.tail ?? [])].filter(Boolean).join("\n")
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
    if (item.kind === "diff") {
      writeSolidToScrollback(renderer, () => (
        <text>
          <span style={{ fg: t.ok }}>{"⏺ "}</span>
          <span style={{ fg: t.text, attributes: TextAttributes.BOLD }}>{String(f.name ?? "edit")}</span>
          <span style={{ fg: t.muted }}>{"  " + String(f.path ?? "")}</span>
        </text>
      ))
      const text = String(f.text ?? "")
      if (isDrawableDiff(text)) {
        await commitSurface(renderer, diffBlock(text, {
          styles: markdownStyles(t), filetype: String(f.filetype ?? "text"), colors: diffColors(t),
        }))
      } else {
        // Never "Error parsing diff" in a transcript: the lines as they are.
        const lines = text.split("\n")
        writeSolidToScrollback(renderer, (ctx) => <Row item={{ id: item.id, kind: "lines", lines }} theme={t} width={ctx.width} />)
      }
      return
    }
    writeSolidToScrollback(renderer, (ctx) => <Row item={item} theme={t} width={ctx.width} />)
  }
}
