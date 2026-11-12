/**
 * What the session still has out, above the prompt (`background` on
 * session.info and session.updated): one row per shell, watched job,
 * background agent or wake-up, with its kind, its name, an age the window
 * counts itself, and the id `/bash <id>` takes. Nothing out, no rows.
 */
import { afterEach, expect, test } from "bun:test"
import { emptyState, ensureSession, applyEvent, backgroundRows } from "../../src/store/state"
import { backgroundLines, BACKGROUND_MAX, age } from "../../src/footer/background"
import { width } from "../../src/text"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const ev = (type: string, fields: Record<string, unknown>) => ({ type, session: "s1", ts: 0, ...fields }) as never
const NOW = Date.now() / 1000
const SHELL = { group: "shells", id: "235d02ab99", kind: "Shell", label: "pytest tests/ -q", url: "", state: "",
                since: NOW - 252, at: null }

test("the rows arrive with the session's facts, each field as its type says", () => {
  const state = emptyState()
  ensureSession(state, "s1")
  applyEvent(state, ev("session.updated", { seq: 1, background: [SHELL] }))
  expect(state.sessions.s1!.background.map((r) => r.id)).toEqual(["235d02ab99"])
  // Another fact changing leaves the list as it is; an empty list clears it.
  applyEvent(state, ev("session.updated", { seq: 2, model: "m2" }))
  expect(state.sessions.s1!.background.length).toBe(1)
  applyEvent(state, ev("session.updated", { seq: 3, background: [] }))
  expect(state.sessions.s1!.background.length).toBe(0)
  expect(backgroundRows([{ id: "x1", label: "a\x1b[31mb\nc", since: "soon" }, { label: "no id" }, "junk", null]))
    .toEqual([{ group: "", id: "x1", kind: "Job", label: "a [31mb c", url: "", state: "", since: null, at: null }])
})

test("a running shell is on screen above the prompt, with its age and its id", async () => {
  const w = open = await mountWindow({ width: 100, height: 30 })
  w.store.push(ev("session.updated", { seq: 1, background: [SHELL] }))
  await w.settle()
  const lines = w.footer()
  const at = lines.findIndex((l) => l.includes("⚙ Shell") && l.includes("pytest tests/ -q"))
  expect(at).toBeGreaterThanOrEqual(0)
  expect(lines[at]).toMatch(/4m1\ds {2}\[235d02ab\]/)
  expect(lines.slice(at).some((l) => l.startsWith("❯"))).toBe(true)      // above the prompt
  w.store.push(ev("session.updated", { seq: 2, background: [] }))
  await w.settle()
  expect(w.footer().join("\n")).not.toContain("⚙ Shell")
})

test("a long list keeps to its rows, and every line fits the terminal", () => {
  const many = backgroundRows(Array.from({ length: 9 }, (_, i) => ({ ...SHELL, id: `job${i}`,
                                                                     label: `run number ${i} `.repeat(20) })))
  for (const cols of [30, 60, 120]) {
    const lines = backgroundLines(many, cols, Date.now())
    expect(lines.length).toBe(BACKGROUND_MAX)
    expect(lines[lines.length - 1]!.text).toContain("6 more")
    for (const l of lines) expect(width(l.kind + l.text + l.age)).toBeLessThanOrEqual(cols - 1)
  }
  expect(backgroundLines([], 80, Date.now())).toEqual([])
})

test("ages and wake-ups read as times", () => {
  expect([age(5), age(252), age(7300)]).toEqual(["5s", "4m12s", "2h01m"])
  const wake = backgroundRows([{ group: "wakeups", id: "w1", kind: "Wake-up", label: "check CI", at: NOW + 600 }])
  expect(backgroundLines(wake, 80, NOW * 1000)[0]!.age).toContain("in 10m00s")
})
