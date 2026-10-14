/**
 * B14: a question that opens while a message is being typed. Keys in a
 * burst are swallowed (and said once), but the first option key after any
 * ordinary pause in the typing (a breath between words, > 350 ms) answers
 * the question: here "and" -> [a] abort (refuse and end the turn). The
 * rest of the message goes nowhere.
 */
import { afterEach, expect, test } from "bun:test"
import { FIXTURE, mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const edit = FIXTURE.events.find((e: any) => e.type === "permission.asked" && e.kind === "confirm")
const before = FIXTURE.events.slice(0, FIXTURE.events.indexOf(edit))

async function typeAt(text: string, gapMs: number) {
  for (const ch of text) {
    if (ch === " ") w!.keys.pressKey(" ")
    else w!.keys.pressKey(ch, ch !== ch.toLowerCase() ? { shift: true } : undefined)
    await Bun.sleep(gapMs)
  }
}

test("typing a message, a pause between words, then 'and'", async () => {
  w = await mount({ answerPauseMs: 350 })
  for (const e of before) w.store.push(e)                 // a turn runs
  await w.settle()
  await typeAt("please re", 90)                             // queueing a message...
  w.store.push(edit)                                        // ...the edit question opens
  await w.settle()
  await typeAt("run the parser tests ", 90)                 // swallowed: a burst
  await Bun.sleep(600)                                      // a breath
  await typeAt("and check the energy", 90)               // after a breath: "and"
  await w.settle()
  console.log("replies:", JSON.stringify(w.replies()))
  console.log(w.footer().filter((l) => l.trim()).slice(0, 2).join("\n"))
  expect(w.replies()).toEqual([])
}, 20000)

test("typing a message into a question that has been open a while answers nothing", async () => {
  w = await mount({ answerPauseMs: 350 })
  w.store.push({ type: "permission.asked", session: "s1", ts: 0, seq: 100, requestId: "c1", kind: "confirm", tool: "bash",
                 preview: "$ pytest -q", flags: { shell: true },
                 options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "a", label: "abort" }] } as never)
  await w.settle()
  await Bun.sleep(1500)                                  // the question has been on screen
  await w.keys.typeText("bash again")                    // ... and the user types a message anyway
  await w.settle()
  expect(w.replies()).toEqual([])
  expect(w.footer().join("\n")).toContain("were not taken as answers")
  await Bun.sleep(1100)
  w.keys.pressKey("y")                                   // after a pause: a deliberate answer
  await w.settle()
  expect(w.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "y" }])
})
