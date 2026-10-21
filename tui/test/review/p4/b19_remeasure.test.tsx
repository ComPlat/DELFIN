/**
 * B19: composer.remeasure() exists for "a width change" but nothing calls
 * it: after a resize the prompt keeps the row count of the old width until
 * the next keystroke -- narrower, a wrapped draft shows fewer rows than it
 * has (its start scrolls out of view); wider, blank rows stay.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

test("a three-row draft, then the terminal gets narrower", async () => {
  w = await mount({ width: 100, height: 30 })
  const draft = "please look at the conformer search output and tell me which of the twelve structures is the lowest in energy, then plot the spectra"
  await w.keys.typeText(draft)
  await w.settle()
  const before = w.footer().filter((l) => l.trim()).length
  w.t.resize(50, 30)
  await Bun.sleep(400)
  await w.settle()
  const shown = w.footer()
  console.log(`footer rows before ${before}; after the resize:\n` + shown.join("\n"))
  expect(shown.join(" ")).toContain("please look at")
})
