/**
 * B02: the history walk and its draft.
 *  (a) readline (the REPL) keeps the line being typed while you look at
 *      history and edit an entry; the TUI forgets it on the first edit.
 *  (b) two prompts that read the same but carry different pastes: the
 *      history keeps the FIRST paste, so ↑ brings back the wrong content.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })
const composerLine = (frame: string) => frame.split("\n").find((l) => l.startsWith("❯ "))?.trimEnd()

test("(a) the draft after an edit of a recalled entry", async () => {
  w = await mount({ history: [{ text: "older prompt" }] })
  await w.keys.typeText("my long careful draft")
  w.keys.pressArrow("up")
  await w.settle()
  expect(composerLine(w.frame())).toBe("❯ older prompt")
  await w.keys.typeText("!")                         // a small edit of the recalled entry
  await w.settle()
  w.keys.pressArrow("down")                          // readline: back to the draft
  await w.settle()
  w.keys.pressArrow("down")
  await w.settle()
  console.log("after ↓↓:", composerLine(w.frame()))
  w.keys.pressArrow("up")
  w.keys.pressArrow("up")
  await w.settle()
  console.log("after ↑↑:", composerLine(w.frame()))
  w.keys.pressArrow("down")
  w.keys.pressArrow("down")
  w.keys.pressArrow("down")
  await w.settle()
  console.log("after ↓↓↓:", composerLine(w.frame()))
  expect(w.frame()).toContain("my long careful draft")
})

test("(b) same text, other paste", async () => {
  w = await mount()
  const logA = Array.from({ length: 30 }, (_, i) => `A ${i} first failure`).join("\n")
  const logB = Array.from({ length: 30 }, (_, i) => `B ${i} second failure`).join("\n")
  await w.keys.typeText("why? ")
  await w.keys.pasteBracketedText(logA)
  w.keys.pressEnter()
  await w.settle()
  await w.keys.typeText("why? ")
  await w.keys.pasteBracketedText(logB)
  w.keys.pressEnter()
  await w.settle()
  expect(w.prompts().map((p) => p.text.split("\n")[0])).toEqual(["why? A 0 first failure", "why? B 0 second failure"])
  w.keys.pressArrow("up")                            // the last prompt, to send it again
  await w.settle()
  w.keys.pressEnter()
  await w.settle()
  const again = w.prompts().at(-1)!.text as string
  console.log("re-sent:", JSON.stringify(again.split("\n")[0]), "entries:", w.history.entries.length)
  expect(again).toBe(`why? ${logB}`)
})
