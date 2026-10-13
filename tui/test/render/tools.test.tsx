/**
 * Tool calls in the window: each becomes its headline and result row in
 * the scrollback, ctrl+o writes the last result out in full, and ctrl+t
 * asks the session to show or hide its task list.
 */
import { afterEach, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const ev = (type: string, fields: Record<string, unknown>) => ({ type, session: "s1", ts: 0, ...fields }) as never

test("each call is its headline and its result, in order", async () => {
  const w = open = await mountWindow({ width: 100 })
  let seq = 0
  const call = (callId: string, name: string, input: object, output: string, meta: object = { ok: true }) => {
    w.store.push(ev("tool.started", { seq: ++seq, turnId: "t1", callId, name, input }))
    w.store.push(ev("tool.ended", { seq: ++seq, turnId: "t1", callId, name, output, meta }))
  }
  w.store.push(ev("turn.started", { seq: ++seq, turnId: "t1", source: "user" }))
  call("c1", "bash", { command: "pytest -q" }, JSON.stringify({ exit_code: 1, stdout: "1 failed, 2 passed\n" }))
  call("c2", "edit_file", { path: "a.py" }, "--- a/a.py\n+++ b/a.py\n@@ -1,1 +1,1 @@\n-x = 1\n+x = 2")
  call("c3", "read_file", { path: "/etc/shadow" }, "", { ok: false, error: "outside the workspace" })
  call("c4", "mcp__delfin-docs__search_docs", { query: "SMD" }, "3 results\nSMD\nCPCM\nnone")
  w.store.push(ev("tool.denied", { seq: ++seq, turnId: "t1", name: "bash" }))
  await w.settle()
  const text = w.t.externalOutput.takeText()
  if (process.env.SHOW) console.log(text)
  const order = ["⏺ bash  pytest -q", "⎿ exit 1 · 1 line", "1 failed, 2 passed",
                 "⏺ edit_file  a.py", "⎿ Edited a.py (+1 -1)", "1 -x = 1", "1 +x = 2",
                 "⏺ read_file  /etc/shadow", "⎿ blocked: outside the workspace",
                 "⏺ delfin-docs:search_docs  SMD", "⎿ 4+ lines, 23 B", "⏺ bash  refused"]
  let at = 0
  for (const part of order) {
    const found = text.indexOf(part, at)
    expect(found).toBeGreaterThanOrEqual(at)
    at = found
  }
})

test("ctrl+o writes the last result out in full; ctrl+t asks for the task list", async () => {
  const w = open = await mountWindow({ width: 100 })
  w.t.mockInput.pressKey("o", { ctrl: true })
  await w.settle()
  expect(w.t.externalOutput.takeText()).toContain("nothing to expand yet")
  const lines = Array.from({ length: 30 }, (_, i) => `result line ${i}`).join("\n")
  w.store.push(ev("tool.started", { seq: 1, turnId: "t1", callId: "c1", name: "bash", input: { command: "ls" } }))
  w.store.push(ev("tool.ended", { seq: 2, turnId: "t1", callId: "c1", name: "bash",
                                  output: JSON.stringify({ exit_code: 0, stdout: lines }), meta: { ok: true } }))
  await w.settle()
  expect(w.t.externalOutput.takeText()).not.toContain("result line 5")      // three lines only
  w.t.mockInput.pressKey("o", { ctrl: true })
  await w.settle()
  const expanded = w.t.externalOutput.takeText()
  expect(expanded).toContain("bash, in full")
  for (const i of [0, 5, 29]) expect(expanded).toContain(`result line ${i}`)
  w.t.mockInput.pressKey("t", { ctrl: true })
  await w.settle()
  expect(w.calls.at(-1)).toEqual({ method: "session.toggleTasks", params: { sessionId: "s1" } })
})
