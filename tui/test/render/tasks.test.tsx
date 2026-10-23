/**
 * The session's task list above the prompt, as the session sends it
 * (`tasks` on session.info and session.updated): one row per task with its
 * mark -- ◼ the one in hand (bold, the accent), ☐ to do, ⏸ waiting, ☒ done
 * (struck through) -- while any task is open; ctrl+t folds it to one line.
 * A label is the model's text: drawn as it is, never as markdown, with
 * control characters made spaces.
 */
import { afterEach, expect, test } from "bun:test"
import { emptyState, ensureSession, applyEvent, taskRows } from "../../src/store/state"
import { taskLines, TASKS_MAX } from "../../src/footer/tasks"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const ev = (type: string, fields: Record<string, unknown>) => ({ type, session: "s1", ts: 0, ...fields }) as never
const ROWS = [
  { num: 2, status: "in_progress", label: "writing the parser" },
  { num: 3, status: "pending", label: "add **tests** for `parse()`" },
  { num: 1, status: "completed", label: "read the input files" },
]

test("the rows arrive with the session's facts, each field as its type says", () => {
  const state = emptyState()
  ensureSession(state, "s1")
  applyEvent(state, ev("session.updated", { seq: 1, tasks: ROWS }))
  expect(state.sessions.s1!.tasks.map((t) => t.status)).toEqual(["in_progress", "pending", "completed"])
  // Another fact changing leaves the list as it is.
  applyEvent(state, ev("session.updated", { seq: 2, model: "m2" }))
  expect(state.sessions.s1!.tasks.length).toBe(3)
  expect(taskRows([{ num: "4", status: "odd", label: "a\x1b[31mb\nc\u0007d" }, "junk", null])).toEqual([
    { num: 4, status: "pending", label: "a [31mb c d" },
  ])
})

test("the list: the task in hand first and marked, the rest under it", async () => {
  const w = open = await mountWindow({ width: 100, height: 30 })
  w.store.push(ev("session.updated", { seq: 1, tasks: ROWS }))
  await w.settle()
  const footer = w.footer().join("\n")
  expect(footer).toContain("⎿ ◼ writing the parser")
  expect(footer).toContain("  ☐ add **tests** for `parse()`")     // as written: no markdown
  expect(footer).toContain("  ☒ read the input files")
  const lines = w.footer()
  const at = lines.findIndex((l) => l.includes("◼ writing the parser"))
  expect(lines.findIndex((l) => l.includes("☒ read the input files"))).toBe(at + 2)
  expect(lines.slice(at).some((l) => l.startsWith("❯"))).toBe(true)      // above the prompt
})

test("ctrl+t folds the list to one line and opens it again", async () => {
  const w = open = await mountWindow({ width: 100, height: 30 })
  w.store.push(ev("session.updated", { seq: 1, tasks: ROWS }))
  await w.settle()
  w.keys.pressKey("t", { ctrl: true })
  await w.settle()
  let footer = w.footer().join("\n")
  expect(footer).not.toContain("☐ add")
  expect(footer).toContain("▸ 3 tasks, 1 done · ◼ writing the parser")
  w.keys.pressKey("t", { ctrl: true })
  await w.settle()
  footer = w.footer().join("\n")
  expect(footer).toContain("☐ add **tests**")
  expect(w.calls.some((c) => c.method === "session.toggleTasks")).toBe(false)
})

test("nothing open, nothing shown", async () => {
  const w = open = await mountWindow({ width: 100, height: 30 })
  w.store.push(ev("session.updated", { seq: 1, tasks: [{ num: 1, status: "completed", label: "done and dusted" }] }))
  await w.settle()
  expect(w.footer().join("\n")).not.toContain("done and dusted")
})

test("a long list keeps to its rows; a short terminal gets the folded line", () => {
  const many = taskRows(Array.from({ length: 12 }, (_, i) => ({ num: i + 1, status: i ? "pending" : "in_progress",
                                                               label: `task ${i + 1}` })))
  const lines = taskLines(many, false, 80, 40)
  expect(lines.length).toBe(TASKS_MAX)
  expect(lines.at(-1)!.text).toBe(`… ${12 - TASKS_MAX + 1} more (ctrl+t folds)`)
  expect(taskLines(many, false, 80, 16).length).toBe(1)
  expect(taskLines(many, true, 80, 40).length).toBe(1)
  // Every line fits the width.
  for (const l of taskLines(taskRows([{ num: 1, status: "pending", label: "x".repeat(300) }]), false, 40, 40)) {
    expect(Bun.stringWidth(l.lead + l.mark + l.text)).toBeLessThanOrEqual(39)
  }
})
