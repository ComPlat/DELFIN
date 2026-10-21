/**
 * B06: keys while the palette is open.
 *  (a) ctrl+c is meant to close the picker (picker.tsx handleKey), but the
 *      app's ladder takes it first: it cancels the running turn, or clears
 *      the draft typed behind the palette; the palette stays open.
 *  (b) a command "with arguments" picked from the palette replaces the
 *      draft (composer.setText): what was typed is gone, not in history.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const COMMANDS = [{ name: "/help", summary: "List these commands", category: "session", takesArgs: true, source: "builtin" }]
function client(calls: any[]) {
  return { request: async (method: string, params: any = {}) => {
    calls.push({ method, params })
    if (method === "session.commands") return { commands: COMMANDS }
    return { accepted: true, stopping: true }
  } }
}

test("(a1) ctrl+c in the palette, a draft behind it", async () => {
  const calls: any[] = []
  w = await mount({ client: client(calls) })
  await w.keys.typeText("a long draft I am still writing")
  w.keys.pressKey("p", { ctrl: true })
  await w.settle()
  expect(w.footer().join("\n")).toContain("commands")
  w.keys.pressKey("c", { ctrl: true })                    // close the palette
  await w.settle()
  const shown = w.footer().join("\n")
  console.log(shown.split("\n").filter((l) => l.trim()).join("\n"))
  expect(shown).toContain("a long draft I am still writing")
  expect(shown).not.toContain("filter:")
})

test("(a2) ctrl+c in the palette during a turn", async () => {
  const calls: any[] = []
  w = await mount({ client: client(calls), turn: true })
  w.keys.pressKey("p", { ctrl: true })
  await w.settle()
  w.keys.pressKey("c", { ctrl: true })                    // close the palette
  await w.settle()
  console.log("calls:", JSON.stringify(calls.map((c) => c.method)), "palette open:", w.footer().join("\n").includes("filter:"))
  expect(calls.some((c) => c.method === "session.cancel")).toBe(false)
})

test("(b) a command with arguments from the palette over a draft", async () => {
  const calls: any[] = []
  w = await mount({ client: client(calls) })
  await w.keys.typeText("a long draft I am still writing")
  w.keys.pressKey("p", { ctrl: true })
  await w.settle()
  await w.keys.typeText("/help")
  await w.settle()
  w.keys.pressEnter()
  await w.settle()
  const line = w.footer().find((l) => l.startsWith("❯"))!
  console.log("composer:", JSON.stringify(line.trimEnd()), "history:", JSON.stringify(w.history.entries))
  w.keys.pressArrow("up")
  await w.settle()
  console.log("after ↑:", JSON.stringify(w.footer().find((l) => l.startsWith("❯"))!.trimEnd()))
  expect(w.frame()).toContain("a long draft I am still writing")
})
