/**
 * The session's task list, above the prompt, while any task is open:
 *
 *   ⎿ ◼ writing the parser          the one in hand: bold, in the accent
 *     ☐ add tests for parse()       to do
 *     ⏸ deploy — waiting on review  waiting on someone
 *     ☒ read the input files        done: faint, struck through
 *     … 3 more (ctrl+t folds)
 *
 * ctrl+t folds it to one line (`▸ 5 tasks, 1 done · ◼ writing the parser`),
 * and a terminal under FOLD_BELOW rows gets that line only. The rows come
 * in the session's order (in hand, waiting, to do, done), so what a long
 * list leaves out is the done end. Every mark is one column wide in every
 * terminal (no emoji), and carries the state without colour (NO_COLOR).
 *
 * A label is the model's text: drawn as a plain string, never as markdown.
 */
import { TextAttributes } from "@opentui/core"
import { For, Show, type JSX } from "solid-js"
import type { TaskRow } from "../store/state"
import type { Theme } from "../theme/theme"
import { fit } from "../text"

/** Rows the list takes at most, its "… more" line included. */
export const TASKS_MAX = 6
/** Below this many terminal rows the list is the folded line only. */
export const FOLD_BELOW = 20

const MARK: Record<TaskRow["status"], string> = { in_progress: "◼", pending: "☐", blocked: "⏸", completed: "☒" }

export interface TaskLine {
  lead: string
  mark: string
  text: string
  status: TaskRow["status"] | "more" | "folded"
}

const open = (rows: readonly TaskRow[]) => rows.some((r) => r.status !== "completed")

/** The lines the list draws in `cols` columns (each a column short of
 *  them): none while no task is open. */
export function taskLines(rows: readonly TaskRow[], folded: boolean, cols: number, termHeight: number,
                          foldKey = "ctrl+t"): TaskLine[] {
  if (!open(rows)) return []
  const room = Math.max(4, cols - 1)
  if (folded || termHeight < FOLD_BELOW) {
    const done = rows.filter((r) => r.status === "completed").length
    const current = rows.find((r) => r.status === "in_progress")
    const head = `${rows.length} task${rows.length === 1 ? "" : "s"}, ${done} done`
    const text = fit(current ? `${head} · ${MARK.in_progress} ${current.label}` : head, room - 2)
    return [{ lead: "▸ ", mark: "", text, status: "folded" }]
  }
  const shown = rows.length > TASKS_MAX ? rows.slice(0, TASKS_MAX - 1) : rows
  const lines: TaskLine[] = shown.map((r, i) => ({
    lead: i === 0 ? "⎿ " : "  ", mark: `${MARK[r.status]} `,
    text: fit(r.label, Math.max(1, room - 4)), status: r.status,
  }))
  if (rows.length > shown.length) {
    lines.push({ lead: "  ", mark: "", text: fit(`… ${rows.length - shown.length} more (${foldKey} folds)`, room - 2),
                 status: "more" })
  }
  return lines
}

export function TasksPanel(props: { lines: () => TaskLine[]; theme: () => Theme }): JSX.Element {
  const t = props.theme
  const tone = (l: TaskLine) => l.status === "in_progress" ? t().accent : l.status === "blocked" ? t().warn
    : l.status === "completed" || l.status === "more" ? t().faint : l.status === "folded" ? t().muted : t().text
  const markTone = (l: TaskLine) => l.status === "in_progress" ? t().accent : l.status === "blocked" ? t().warn : t().faint
  const attrs = (l: TaskLine) => l.status === "in_progress" ? TextAttributes.BOLD
    : l.status === "completed" ? TextAttributes.STRIKETHROUGH : 0
  return (
    <Show when={props.lines().length}>
      <box flexShrink={0} flexDirection="column">
        <For each={props.lines()}>
          {(l) => (
            <text wrapMode="none">
              <span style={{ fg: t().faint }}>{l.lead}</span>
              <span style={{ fg: markTone(l) }}>{l.mark}</span>
              <span style={{ fg: tone(l), attributes: attrs(l) }}>{l.text}</span>
            </text>
          )}
        </For>
      </box>
    </Show>
  )
}
