/**
 * Tool calls in the window: each becomes its headline and result row in
 * the scrollback, and ctrl+o writes the last result out in full. (ctrl+t
 * folds the task list in the window: render/tasks.test.tsx.)
 */
import { afterEach, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const ev = (type: string, fields: Record<string, unknown>) => ({ type, session: "s1", ts: 0, ...fields }) as never

/** The calls every kind of tool row is drawn for, in one turn. */
function calls(w: Awaited<ReturnType<typeof mountWindow>>) {
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
}
const inOrder = (text: string, order: string[]) => {
  let at = 0
  for (const part of order) {
    const found = text.indexOf(part, at)
    expect(found).toBeGreaterThanOrEqual(at)
    at = found
  }
}
const HEADLINES = ["⏺ bash  pytest -q", "⎿ exit 1 · 1 line", "⏺ edit_file  a.py", "⎿ Edited a.py (+1 -1)",
                   "⏺ read_file  /etc/shadow", "⎿ blocked: outside the workspace",
                   "⏺ delfin-docs:search_docs  SMD", "⎿ 4+ lines, 23 B", "⏺ bash  refused"]

test("each call is its headline and what came back, in order -- nothing under it", async () => {
  const w = open = await mountWindow({ width: 100 })
  calls(w)
  await w.settle()
  const text = w.t.externalOutput.takeText()
  if (process.env.SHOW) console.log(text)
  inOrder(text, HEADLINES)
  for (const detail of ["1 failed, 2 passed", "1 -x = 1", "1 +x = 2"]) expect(text).not.toContain(detail)
})

test("with its details, a call shows what it printed and an edit its diff", async () => {
  const w = open = await mountWindow({ width: 100, height: 40, full: true })
  w.screen!.setDetails(true)
  calls(w)
  await w.settle()
  const text = w.screenLines().join("\n")
  inOrder(text, ["⏺ bash  pytest -q", "⎿ exit 1 · 1 line", "1 failed, 2 passed",
                 "⏺ edit_file  a.py", "⎿ Edited a.py (+1 -1)", "1 -x = 1", "1 +x = 2",
                 "⏺ read_file  /etc/shadow", "⎿ blocked: outside the workspace",
                 "⏺ delfin-docs:search_docs  SMD", "⎿ 4+ lines, 23 B", "⏺ bash  refused"])
})

test("the palette writes the last result out in full", async () => {
  const w = open = await mountWindow({ width: 100, height: 40 })
  const expandLast = async () => {
    w.t.mockInput.pressKey("p", { ctrl: true })
    await w.settle()
    await w.t.mockInput.typeText("expand the last")
    w.t.mockInput.pressEnter()
    await w.settle()
  }
  await expandLast()
  expect(w.t.externalOutput.takeText()).toContain("nothing to expand yet")
  const lines = Array.from({ length: 30 }, (_, i) => `result line ${i}`).join("\n")
  w.store.push(ev("tool.started", { seq: 1, turnId: "t1", callId: "c1", name: "bash", input: { command: "ls" } }))
  w.store.push(ev("tool.ended", { seq: 2, turnId: "t1", callId: "c1", name: "bash",
                                  output: JSON.stringify({ exit_code: 0, stdout: lines }), meta: { ok: true } }))
  await w.settle()
  expect(w.t.externalOutput.takeText()).not.toContain("result line 5")      // its headline only
  await expandLast()
  const expanded = w.t.externalOutput.takeText()
  expect(expanded).toContain("bash, in full")
  for (const i of [0, 5, 29]) expect(expanded).toContain(`result line ${i}`)
})
