/**
 * The /swarm watch drawn in the window (package SWV, phase 2b): one row
 * per member -- name, state, commits, last tool call -- the selection
 * marked, a member's feed opened with Enter in place of the table. The
 * SwarmPanel is the pattern (tui/src/swarm/panel.tsx): plain text rows,
 * nothing that does not fit cut at the top. Read-only by construction:
 * the panel renders the watch state; it takes no input of its own.
 */
import { For, Show, type JSX } from "solid-js"
import { fit, width } from "../text"
import type { Theme } from "../theme/theme"
import { memberLine, createSwarmWatch } from "./watch"

export type SwarmWatchView = ReturnType<typeof createSwarmWatch>

export function SwarmWatchPanel(props: { view: SwarmWatchView; theme: () => Theme;
                                         width: () => number; height: () => number }): JSX.Element {
  const t = props.theme
  const v = props.view
  const inner = () => Math.max(20, props.width() - 2)
  const feed = () => v.state.feed
  const head = (): string => feed()
    ? ` member ${feed()!.name} — recent steps (read-only)`
    : ` ${v.state.rows.length} member${v.state.rows.length === 1 ? "" : "s"} · read-only · enter opens a member's steps`
  const lines = (): { text: string; kind: "row" | "title"; on?: boolean }[] => {
    if (feed()) {
      const out: { text: string; kind: "row" | "title"; on?: boolean }[] = []
      const entries = feed()!.entries
      if (!entries.length) out.push({ kind: "row", text: "  no recent steps" })
      entries.forEach((e, i) => {
        out.push({ kind: "row", on: i === v.state.feedSelected,
                   text: `  ${e.kind === "refusal" ? "✗" : e.kind === "commit" ? "+" : "·"} ${e.text}` })
      })
      return out
    }
    const out: { text: string; kind: "row" | "title"; on?: boolean }[] = []
    if (v.state.loading) out.push({ kind: "row", text: "  reading the swarm…" })
    if (v.state.error) out.push({ kind: "row", text: `  ${v.state.error}` })
    v.state.rows.forEach((row, i) => {
      out.push({ kind: "row", on: i === v.state.selected, text: `  ${memberLine(row)}` })
    })
    if (!v.state.loading && !v.state.rows.length && !v.state.error)
      out.push({ kind: "row", text: "  no members" })
    return out
  }
  const room = () => Math.max(1, props.height() - 1)
  const shown = () => lines().slice(Math.max(0, lines().length - room()))
  return (
    <box flexDirection="column" height={props.height()} flexShrink={0}>
      <text wrapMode="none">
        <span style={{ fg: t().accent, bold: true }}>{"◆ SWARM WATCH"}</span>
        <span style={{ fg: t().muted }}>{fit(head(), Math.max(8, inner() - width("◆ SWARM WATCH") - 2))}</span>
      </text>
      <Show when={!v.state.loading || feed()}>
        <For each={shown()}>{(line) => {
          const tone = line.kind === "title" ? t().faint
            : line.on ? t().accent : t().text
          return <text wrapMode="none" fg={tone}>{fit(line.text, inner())}</text>
        }}</For>
      </Show>
    </box>
  )
}
