/**
 * B05: tool rows against repl_render.py.
 *  (a) the REPL cuts a long argument in its MIDDLE (truncate_middle: "the
 *      tail ... says which file is about to be written"); the TUI cuts the
 *      end, so a long path loses its file name -- in the headline and in
 *      the "Edited <path>" line.
 *  (b) the REPL masks a credential in a command headline
 *      (_headline_value -> cli_approve.redact); the TUI prints it whole.
 *      The worker's redactor does not catch an assignment such as
 *      `export API_TOKEN=hunter2xyz` (checked in b05_redact.py).
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const deep = "/pfs/work7/workspace/scratch/ab1234-delfin-project/calculations/benzene_dimer/opt_freq/run_042/orca_input.inp"
const diff = [`--- a/${deep}`, `+++ b/${deep}`, "@@ -1,2 +1,2 @@", " ! B3LYP def2-SVP", "-%pal nprocs 4 end", "+%pal nprocs 8 end"].join("\n")

test("(a) a long path keeps its file name", async () => {
  w = await mount({ width: 80 })
  w.store.push({ type: "tool.started", session: "s1", ts: 0, seq: 1, turnId: "t1", callId: "c1", name: "edit_file",
                 input: { path: deep } } as never)
  w.store.push({ type: "tool.ended", session: "s1", ts: 0, seq: 2, turnId: "t1", callId: "c1", name: "edit_file",
                 output: diff, meta: { ok: true } } as never)
  await w.settle()
  const rows = w.t.externalOutput.takeText().split("\n").filter((l) => /⏺|⎿/.test(l))
  console.log(rows.join("\n"))
  for (const row of rows) expect(row).toContain("orca_input.inp")
})

test("(b) a credential in a command headline", async () => {
  w = await mount({ width: 100 })
  const command = "export API_TOKEN=hunter2xyz && curl -s https://api.example.org/v1/jobs"
  w.store.push({ type: "tool.started", session: "s1", ts: 0, seq: 1, turnId: "t1", callId: "c1", name: "bash",
                 input: { command } } as never)
  w.store.push({ type: "tool.ended", session: "s1", ts: 0, seq: 2, turnId: "t1", callId: "c1", name: "bash",
                 output: '{"exit_code": 0, "stdout": "[]\\n"}', meta: { ok: true } } as never)
  await w.settle()
  const rows = w.t.externalOutput.takeText().split("\n").filter((l) => /⏺/.test(l))
  console.log(rows.join("\n"))
  expect(rows.join("\n")).not.toContain("hunter2xyz")
})
