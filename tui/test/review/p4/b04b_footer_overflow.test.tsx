/**
 * B04: when the footer needs more rows than the terminal has, footerRows
 * caps it at height-1 and Yoga squeezes the children: rows are drawn over
 * each other. The approval's diff then reads as lines that do not exist.
 * Here at the classic 80x24 with a 3-line draft, and a fourth line typed
 * on after the question opened (the prompt keeps the keys).
 */
import { afterEach, expect, test } from "bun:test"
import { FIXTURE, mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const edit = FIXTURE.events.find((e: any) => e.type === "permission.asked" && e.kind === "confirm")
const upToEdit = FIXTURE.events.slice(0, FIXTURE.events.indexOf(edit) + 1)

test("80x24: a draft of three lines, then the question", async () => {
  w = await mount({ width: 80, height: 24, answerPauseMs: 350 })
  for (const e of upToEdit.slice(0, -1)) w.store.push(e)
  await w.settle()
  await w.keys.typeText("first line of a queued message")
  w.keys.pressEnter({ meta: true })
  await w.keys.typeText("second line")
  w.keys.pressEnter({ meta: true })
  await w.keys.typeText("third")
  w.store.push(edit)
  await w.settle()
  w.keys.pressEnter({ meta: true })            // typed on: the draft grows under the panel
  await w.keys.typeText("fourth")
  await w.settle()
  const frame = w.frame()
  console.log(`--- footerHeight=${w.t.renderer.footerHeight}\n${frame}`)
  // Nothing is drawn over anything: every preview line on screen is a
  // line of the preview, as it is, and what did not fit is announced.
  const real = new Set(String(edit.preview).split("\n").map((l: string) => l.trimEnd()))
  const panel = frame.split("\n").filter((l) => l.startsWith("│ ")).map((l) => l.slice(2).replace(/\s*│\s*$/, ""))
  const drawn = panel.filter((l) => !/^(\[|… \d+ more lines|keys typed)/.test(l))
  expect(drawn.length).toBeGreaterThan(5)
  for (const line of drawn) expect(real.has(line.trimEnd())).toBe(true)
  expect(frame).toMatch(/… \d+ more lines \(ctrl\+o shows all\)/)
  expect(frame).toContain("[y] yes")
})
