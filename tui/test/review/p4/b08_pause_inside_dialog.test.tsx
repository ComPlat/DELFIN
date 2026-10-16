/**
 * B08: the typeahead pause also applies between the keys of ONE dialog
 * (`t - state.lastKeyAt < pause`), long after the question opened. A
 * multi-select answered at a normal typing pace loses picks silently (the
 * note is said once per question, and says "a question just opened").
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const multi = {
  type: "permission.asked", session: "s1", ts: 0, seq: 100, requestId: "q1", kind: "ask", tool: "ask_user_question",
  preview: "", options: [], flags: {}, question: "Which spectra should I compute?", header: "Pick",
  choices: ["IR", "NMR", "UV-Vis"], choiceDetails: ["", "", ""], multiSelect: true,
}

test("1 and 3, 200 ms apart, two seconds after the question opened", async () => {
  w = await mount({ answerPauseMs: 350 })
  w.store.push(multi as never)
  await w.settle()
  await Bun.sleep(2000)                         // the user reads the question
  w.keys.pressKey("1")
  await Bun.sleep(200)
  w.keys.pressKey("3")
  await Bun.sleep(600)                          // looks, then sends
  w.keys.pressEnter()
  await w.settle()
  console.log("replies:", JSON.stringify(w.replies()))
  console.log(w.footer().filter((l) => l.trim()).join("\n"))
  expect(w.replies()).toEqual([{ sessionId: "s1", requestId: "q1", answers: ["IR", "UV-Vis"] }])
})

test("A then y, 250 ms apart: typing (\"Ay\"), not a saved rule -- y after the hold saves it", async () => {
  w = await mount({ answerPauseMs: 350 })
  w.store.push({ type: "permission.asked", session: "s1", ts: 0, seq: 100, requestId: "c1", kind: "confirm", tool: "bash",
                 preview: "$ pytest -q", flags: { shell: true }, options: [{ key: "y", label: "yes" }, { key: "n", label: "no" },
                 { key: "A", label: "always this exact command" }, { key: "a", label: "abort" }] } as never)
  await w.settle()
  await Bun.sleep(1500)
  w.keys.pressKey("A", { shift: true })
  await Bun.sleep(250)
  w.keys.pressKey("y")
  await w.settle()
  await Bun.sleep(650)
  await w.settle()
  // A key with another after it within HOLD_MS is the start of a message
  // ("An..." would be A, then n "just this once"): both go into the prompt.
  expect(w.replies()).toEqual([])
  expect(w.footer().join("\n")).toContain("❯ Ay")
  w.keys.pressBackspace()
  w.keys.pressBackspace()
  await w.settle()
  await Bun.sleep(1100)
  w.keys.pressKey("A", { shift: true })
  await Bun.sleep(650)                          // alone: the second step shows
  await w.settle()
  expect(w.footer().join("\n")).toContain("persist")
  w.keys.pressKey("y")
  await Bun.sleep(650)
  await w.settle()
  console.log("replies:", JSON.stringify(w.replies()))
  console.log(w.footer().filter((l) => l.trim()).slice(0, 3).join("\n"))
  expect(w.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "A", persist: true }])
}, 10000)
