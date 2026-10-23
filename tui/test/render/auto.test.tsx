/**
 * Auto mode in the mode line: the session sends its one-line state as
 * `auto` on session.info / session.updated ("" while off), and the hints
 * line draws it after the posture, before the model. It is the session's
 * text: drawn as one plain line, control characters made spaces.
 */
import { afterEach, expect, test } from "bun:test"
import { emptyState, ensureSession, applyEvent } from "../../src/store/state"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const ev = (type: string, fields: Record<string, unknown>) => ({ type, session: "s1", ts: 0, ...fields }) as never

test("the state arrives with the session's facts and goes when it is off", () => {
  const state = emptyState()
  ensureSession(state, "s1")
  expect(state.sessions.s1!.auto).toBe("")
  applyEvent(state, ev("session.updated", { seq: 1, auto: "loop table 3/12 · kept · 1 need you" }))
  expect(state.sessions.s1!.auto).toBe("loop table 3/12 · kept · 1 need you")
  applyEvent(state, ev("session.updated", { seq: 2, model: "m2" }))
  expect(state.sessions.s1!.auto).toBe("loop table 3/12 · kept · 1 need you")
  applyEvent(state, ev("session.updated", { seq: 3, auto: "a\x1b[31mb\nc" }))
  expect(state.sessions.s1!.auto).toBe("a [31mb c")
  applyEvent(state, ev("session.updated", { seq: 4, auto: "" }))
  expect(state.sessions.s1!.auto).toBe("")
})

test("the hints line shows it after the posture", async () => {
  const w = open = await mountWindow({ width: 140, height: 30 })
  w.store.push(ev("session.updated", { seq: 1, approval: "acceptEdits", model: "m1",
                                       auto: "loop table 3/12 · kept" }))
  await w.settle()
  const hints = w.footer().find((l) => l.includes("m1")) ?? ""
  expect(hints).toContain("loop table 3/12 · kept")
  expect(hints.indexOf("loop table")).toBeLessThan(hints.indexOf("m1"))
  w.store.push(ev("session.updated", { seq: 2, auto: "" }))
  await w.settle()
  expect(w.footer().join("\n")).not.toContain("loop table")
})
