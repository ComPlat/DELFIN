/**
 * The /swarm view drawn in the window (SW2b): the head with the goal and
 * the counts, the packages as cards (phase, members with health), the
 * "waits on you" queue, and the held plan / held answer. The agents
 * view's panel is the pattern (tui/src/agents/view.tsx AgentViewPanel):
 * plain text rows, the selection marked, nothing that does not fit cut
 * at the top.
 */
import { For, type JSX } from "solid-js"
import { fit, width } from "../text"
import type { Theme } from "../theme/theme"
import { createSwarmView } from "./view"

export type SwarmView = ReturnType<typeof createSwarmView>

export function SwarmPanel(props: { view: SwarmView; theme: () => Theme; width: () => number;
                                    height: () => number }): JSX.Element {
  const t = props.theme
  const v = props.view
  const inner = () => Math.max(20, props.width() - 2)
  const goal = () => v.state.status?.goal ?? ""
  const counts = () => {
    const rows = v.rows()
    const waits = v.waits()
    return [`${rows.length} package${rows.length === 1 ? "" : "s"}`,
            waits.length ? `${waits.length} wait${waits.length === 1 ? "s" : ""} on you` : ""]
      .filter(Boolean).join(" · ")
  }
  /** The lines: the head, the cards, the queue, the held plan or answer. */
  const lines = (): { text: string; kind: "head" | "title" | "row" | "held" | "queue";
                      on?: boolean }[] => {
    const out: { text: string; kind: "head" | "title" | "row" | "held" | "queue"; on?: boolean }[] = []
    const rows = v.rows()
    rows.forEach((row, i) => {
      out.push({ kind: "row", on: i === v.state.selected && !v.state.queueOpen,
                 text: `  ${row.name}  ·  ${row.phase}  ·  ${row.health}` })
    })
    if (!rows.length) out.push({ kind: "row", text: "  no swarm: type a goal in the prompt, then y" })
    if (v.waits().length) {
      out.push({ kind: "title", text: " WAITS ON YOU" })
      v.waits().forEach((item, i) => {
        out.push({ kind: "queue", on: v.state.queueOpen && i === v.state.queueSelected,
                   text: `  ${item.question ?? item.item_id}${item.package ? `  (${item.package})` : ""}` })
      })
    }
    if (v.state.pendingStart) {
      out.push({ kind: "held", text: `  held plan: y to start · n to let it go` })
    } else if (v.state.pendingAnswer) {
      out.push({ kind: "held", text: `  held answer: type it, enter sends` })
    }
    return out
  }
  const room = () => Math.max(1, props.height() - 1)
  const shown = () => lines().slice(Math.max(0, lines().length - room()))
  return (
    <box flexDirection="column" height={props.height()} flexShrink={0}>
      <text wrapMode="none">
        <span style={{ fg: t().accent, bold: true }}>{"◆ SWARM  "}</span>
        <span style={{ fg: t().muted }}>{fit(goal(), Math.max(8, inner() - width("◆ SWARM  ") - width(counts()) - 4))}</span>
        <span style={{ fg: t().faint }}>{`  ${counts()}`}</span>
      </text>
      <For each={shown()}>{(line) => {
        const tone = line.kind === "title" ? t().faint
          : line.kind === "held" ? t().warn
          : line.on ? t().accent : t().text
        return <text wrapMode="none" fg={tone}>{fit(line.text, inner())}</text>
      }}</For>
    </box>
  )
}
