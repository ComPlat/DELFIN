/**
 * What the session still has out, above the prompt, while anything is:
 *
 *   ⚙ Shell   pytest tests/ -q                       4m12s  [235d02ab]
 *   ⚙ Server  http://127.0.0.1:8000                    40s  [9c1e77f0]
 *   ⚙ Agent   explore · where is the parser            1m05s [a1b2c3d4]
 *     … 2 more (/bash lists them)
 *
 * The rows come from the session (`background` on session.info and
 * session.updated: background_view.snapshot, the dashboard panel's
 * rows) without ages; the age is counted here, so the session sends the
 * list when it changes and not every second. The id in brackets is what
 * `/bash <id>` takes. A label is the model's or a command's text: drawn
 * as a plain string.
 */
import { For, Show, type JSX } from "solid-js"
import type { BackgroundRow } from "../store/state"
import type { Theme } from "../theme/theme"
import { fit, width } from "../text"

/** Rows the list takes at most, its "… more" line included. */
export const BACKGROUND_MAX = 4

export interface BackgroundLine {
  kind: string
  text: string
  age: string
  more: boolean
}

export function age(seconds: number): string {
  const s = Math.max(0, Math.floor(seconds))
  if (s < 60) return `${s}s`
  if (s < 3600) return `${Math.floor(s / 60)}m${String(s % 60).padStart(2, "0")}s`
  return `${Math.floor(s / 3600)}h${String(Math.floor((s % 3600) / 60)).padStart(2, "0")}m`
}

function when(row: BackgroundRow, nowS: number): string {
  if (row.at != null) return row.at > nowS ? `in ${age(row.at - nowS)}` : "due"
  const ran = row.since != null ? age(nowS - row.since) : ""
  return row.state && row.group === "watches" ? [row.state, ran].filter(Boolean).join(" ") : ran
}

/** The lines the list draws in `cols` columns (each a column short of
 *  them) at `nowMs`: none while nothing is out. */
export function backgroundLines(rows: readonly BackgroundRow[], cols: number, nowMs: number): BackgroundLine[] {
  if (!rows.length) return []
  const room = Math.max(10, cols - 1)
  const nowS = nowMs / 1000
  const shown = rows.length > BACKGROUND_MAX ? rows.slice(0, BACKGROUND_MAX - 1) : rows
  const kindWidth = Math.min(12, Math.max(...shown.map((r) => width(r.kind))))
  const lines: BackgroundLine[] = shown.map((r) => {
    const k = fit(r.kind, kindWidth)
    const kind = `⚙ ${k}${" ".repeat(Math.max(0, kindWidth - width(k)))}  `
    // A narrow terminal keeps the name and the age; the id goes first.
    const right = room >= 50 ? `  ${when(r, nowS)}  [${r.id.slice(0, 8)}]` : `  ${when(r, nowS)}`
    const name = r.url || r.label
    const text = fit(name, Math.max(4, room - width(kind) - width(right)))
    const pad = Math.max(0, room - width(kind) - width(right) - width(text))
    return { kind, text: text + " ".repeat(pad), age: right, more: false }
  })
  if (rows.length > shown.length) {
    lines.push({ kind: "  ", text: fit(`… ${rows.length - shown.length} more (/bash lists them)`, room - 2),
                 age: "", more: true })
  }
  return lines
}

export function BackgroundPanel(props: { lines: () => BackgroundLine[]; theme: () => Theme }): JSX.Element {
  const t = props.theme
  return (
    <Show when={props.lines().length}>
      <box flexShrink={0} flexDirection="column">
        <For each={props.lines()}>
          {(l) => (
            <text wrapMode="none">
              <span style={{ fg: t().muted }}>{l.kind}</span>
              <span style={{ fg: l.more ? t().faint : t().text }}>{l.text}</span>
              <span style={{ fg: t().faint }}>{l.age}</span>
            </text>
          )}
        </For>
      </box>
    </Show>
  )
}
