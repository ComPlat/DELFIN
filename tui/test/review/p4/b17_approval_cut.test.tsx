/**
 * B17: the approval preview cuts a long line at its END (fit()), where the
 * REPL's render_request cuts it in the middle (truncate_middle: "cutting
 * the tail hides which file is about to be written"). A long command's
 * target, or a deep path's file name, is not on screen when the user
 * answers. And the preview stops at 12 lines (REPL: 24) with no way to see
 * the rest -- a plan of 30 lines is approved after reading 11.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const opts = [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "A", label: "always this exact command" },
              { key: "a", label: "abort" }, { key: "?", label: "explain" }]

test("a long command at 80 columns", async () => {
  w = await mount({ width: 80, height: 30 })
  const command = "cd /pfs/work7/workspace/scratch/ab1234-delfin/calculations/benzene_dimer && rm -rf ~/.ssh/known_hosts"
  w.store.push({ type: "permission.asked", session: "s1", ts: 0, seq: 100, requestId: "c1", kind: "confirm", tool: "bash",
                 preview: `$ ${command}`, flags: { shell: true }, options: opts } as never)
  await w.settle()
  const panel = w.footer().filter((l) => l.startsWith("│") || l.startsWith("╭"))
  console.log(panel.join("\n"))
  expect(panel.join("\n")).toContain("known_hosts")
})

test("a plan of 30 lines", async () => {
  w = await mount({ width: 100, height: 45 })
  const plan = Array.from({ length: 30 }, (_, i) => `${i + 1}. step ${i + 1}${i === 25 ? ": delete the old calculations under ~/runs" : ""}`).join("\n")
  w.store.push({ type: "permission.asked", session: "s1", ts: 0, seq: 100, requestId: "pl1", kind: "plan", tool: "exit_plan_mode",
                 preview: plan, flags: {}, options: [{ key: "d", label: "approve → default" },
                 { key: "e", label: "approve → acceptEdits" }, { key: "n", label: "reject" }] } as never)
  await w.settle()
  const panel = w.footer().filter((l) => l.startsWith("│"))
  // The REPL's 24 lines, the rest announced -- and one key away.
  expect(panel.join("\n")).toContain("more lines (ctrl+o shows all)")
  w.keys.pressKey("o", { ctrl: true })
  await w.settle()
  expect(w.t.externalOutput.takeText()).toContain("26. step 26: delete the old calculations under ~/runs")
  expect(w.replies()).toEqual([])
})
