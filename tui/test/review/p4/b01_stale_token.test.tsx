/**
 * B01: a path completion list is applied with the token it was opened for.
 * The user types "@del" (list opens), types "f" (a new request goes out,
 * the old list stays on screen), and presses Tab before the answer: the
 * candidate replaces only the OLD token, and the "f" is left behind it.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

test("tab while the answer for the edited token is still on its way", async () => {
  let slow: ((v: unknown) => void) | null = null
  const client = {
    request: async (method: string, params: any = {}) => {
      if (method === "session.completePath") {
        const all = ["delfin/", "delta.txt", "docs/"]
        const frag = String(params.fragment)
        const answer = { paths: all.filter((p) => p.startsWith(frag) && p !== frag) }
        if (frag === "delf") return new Promise((resolve) => { slow = () => resolve(answer) })  // in flight
        return answer
      }
      if (method === "session.commands") return { commands: [] }
      return { accepted: true }
    },
  }
  w = await mount({ client })
  await w.keys.typeText("look at @del")
  await w.settle()
  expect(w.frame()).toContain("@delta.txt")             // the list for "@del"
  await w.keys.typeText("f")                            // "@delf": asked again, answer pending
  await w.settle()
  w.keys.pressArrow("down")                             // the second candidate, @delta.txt
  w.keys.pressTab()
  await w.settle()
  const line = w.frame().split("\n").find((l) => l.startsWith("❯ "))!
  console.log("composer line:", JSON.stringify(line.trimEnd()))
  // Expected: nothing applied for a token that changed, or "look at @delf…" kept intact.
  expect(line).not.toContain("@delta.txt f")
  ;(slow as any)?.()
})

test("enter on the stale list for a path does the same", async () => {
  let release: (() => void) | null = null
  const client = {
    request: async (method: string, params: any = {}) => {
      if (method === "session.completePath") {
        const frag = String(params.fragment)
        const all = ["delfin/", "delfin/xtb.py", "docs/"]
        const answer = { paths: all.filter((p) => p.startsWith(frag) && p !== frag) }
        if (frag !== "del") return new Promise((resolve) => { release = () => resolve(answer) })
        return answer
      }
      return { accepted: true }
    },
  }
  w = await mount({ client })
  await w.keys.typeText("see @del")
  await w.settle()
  await w.keys.typeText("fin/x")                        // typed on; the old list for "@del" still shown
  await w.settle()
  w.keys.pressEnter()
  await w.settle()
  const line = w.frame().split("\n").find((l) => l.startsWith("❯ "))!
  console.log("composer line:", JSON.stringify(line.trimEnd()), "prompts:", JSON.stringify(w.prompts()))
  expect(line).not.toContain("@delfin/fin/x")
  ;(release as any)?.()
})
