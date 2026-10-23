/**
 * B10: completion rows are <text> without wrapMode="none": a long command
 * summary (a skill's description is a sentence or two) wraps, while
 * composer.height() counts one row per candidate -- the footer is too
 * short for what it draws, and rows are drawn over each other. Same for
 * the ctrl+r line with a long query.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const COMMANDS = [
  { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" },
  { name: "/spectra", summary: "Compute IR, Raman and UV-Vis spectra for the selected conformers with ORCA, then plot them against the experimental reference and write a short report", category: "skills", takesArgs: true, source: "skill" },
  { name: "/session", summary: "List or search past sessions", category: "history", takesArgs: true, source: "builtin" },
]

test("a long summary in the /command list at 80 columns", async () => {
  const client = { request: async (method: string) => method === "session.commands" ? { commands: COMMANDS } : { accepted: true } }
  w = await mount({ width: 80, height: 24, client })
  await w.keys.typeText("/s")
  await w.settle()
  const frame = w.frame().split("\n")
  const footer = frame.slice(frame.length - w.t.renderer.footerHeight - 1)
  console.log(`footerHeight=${w.t.renderer.footerHeight}\n` + footer.join("\n"))
  expect(footer.some((l) => l.startsWith("❯ /s"))).toBe(true)
  expect(footer.filter((l) => /^─{20,}/.test(l)).length).toBe(2)
  for (const name of ["/status", "/spectra", "/session"]) expect(footer.join("\n")).toContain(name)
})

test("a long ctrl+r query at 40 columns", async () => {
  w = await mount({ width: 40, height: 24, history: [{ text: "run the conformer search on the benzene dimer again" }] })
  w.keys.pressKey("r", { ctrl: true })
  await w.keys.typeText("conformer search on the benzene")
  await w.settle()
  const frame = w.frame().split("\n")
  const footer = frame.slice(frame.length - w.t.renderer.footerHeight - 1)
  console.log(`footerHeight=${w.t.renderer.footerHeight}\n` + footer.join("\n"))
  expect(footer.some((l) => l.startsWith("❯ "))).toBe(true)
  expect(footer.filter((l) => /^─{20,}/.test(l)).length).toBe(2)
})
