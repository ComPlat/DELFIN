/**
 * The gutter: everything the window draws starts a column in (theme
 * GUTTER) -- the scrollback's rows, its markdown and its dividers, and the
 * footer -- and what spans the width ends a column short of the right
 * edge, so the two margins match.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { GUTTER } from "../../src/theme/theme"
import { FIXTURE, mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const PAD = " ".repeat(GUTTER)

async function conversation(width: number) {
  w = await mountWindow({ width, height: 30 })
  // The session's header first, as the service sends it (the banner).
  w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, ...FIXTURE.session } as never)
  for (const e of FIXTURE.events) w.store.push(e)
  await w.settle()
  w.view.show("s1", { divider: "now: the parser" })
  await w.settle()
  return { scrollback: w.t.externalOutput.takeText().split("\n"), frame: w.t.captureCharFrame().split("\n") }
}

describe("the gutter", () => {
  for (const width of [100, 40]) {
    test(`keeps column 0 free in the scrollback at ${width} columns`, async () => {
      const { scrollback } = await conversation(width)
      const drawn = scrollback.filter((l) => l.trim())
      expect(drawn.length).toBeGreaterThan(20)
      for (const line of drawn) expect(line.startsWith(PAD)).toBe(true)
      // The markdown answer is in it too, a column in.
      expect(drawn.some((l) => l.startsWith(PAD + "The workspace holds one README"))).toBe(true)
    })

    test(`keeps column 0 free in the footer at ${width} columns`, async () => {
      const { frame } = await conversation(width)
      const footer = frame.slice(frame.length - w!.t.renderer.footerHeight - 1).filter((l) => l.trim())
      expect(footer.some((l) => l.startsWith(PAD + "❯ "))).toBe(true)
      for (const line of footer) expect(line.startsWith(PAD)).toBe(true)
    })
  }

  test("a divider and the composer's rules end a column short of the right edge", async () => {
    const { scrollback, frame } = await conversation(80)
    const divider = scrollback.find((l) => l.includes("now: the parser"))!
    expect(divider.startsWith(PAD + "━")).toBe(true)
    expect(divider.trimEnd().length).toBe(80 - 1)
    const rules = frame.filter((l) => /^\s─{10,}/.test(l))
    expect(rules.length).toBe(2)
    for (const rule of rules) expect(rule.trimEnd().length).toBe(80 - 1)
  })
})

describe("the banner", () => {
  test("has the dolphin beside the session's words, where there is room", async () => {
    const { scrollback } = await conversation(100)
    const { DOLPHIN } = await import("../../src/theme/logo")
    const at = scrollback.findIndex((l) => l.includes(DOLPHIN[1]!.trimEnd()))
    expect(at).toBeGreaterThanOrEqual(0)
    expect(scrollback[at]).toContain("DELFIN  fake/fake-model")
    expect(scrollback[at - 1]!.trimEnd()).toBe(PAD + DOLPHIN[0]!.trimEnd())
    for (let i = 2; i < DOLPHIN.length; i++) expect(scrollback[at + i - 1]!.startsWith(PAD + DOLPHIN[i]!.trimEnd())).toBe(true)
  })

  test("is the words alone, under the mark, where there is not", async () => {
    const { scrollback } = await conversation(40)
    expect(scrollback.some((l) => l.includes("▟"))).toBe(false)
    expect(scrollback.some((l) => l.startsWith(PAD + "◆ DELFIN"))).toBe(true)
  })
})

describe("the banner's long lines", () => {
  // What a real session's banner says: long lines that wrap at a narrow width.
  const info = {
    provider: "kit", model: "kit.deepseek-v4-flash", mode: "solo", title: "solvent scan for benzene",
    banner: ["delfin-agent · kit/kit.deepseek-v4-flash · solo",
             "workspace  ~/project/delfin  (git tui-rewrite · 3 uncommitted)",
             "approval   default  [from ~/.delfin/settings.json]",
             "isolation  on (bwrap) — a command the agent runs writes only inside the workspace (--no-isolate)"],
  }
  for (const width of [40, 52, 60, 100]) {
    test(`are all there, wrapped, at ${width} columns`, async () => {
      w = await mountWindow({ width, height: 30 })
      w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, workspace: "/w", ...info } as never)
      await w.settle()
      // Read without the wrapping and without the dolphin beside it.
      const flat = (t: string) => t.replace(/[\u2580-\u259f]/g, "").replace(/\s+/g, "")
      const text = flat(w.t.externalOutput.takeText())
      for (const words of ["3 uncommitted)", "[from ~/.delfin/settings.json]", "writes only inside the workspace (--no-isolate)",
                           "solvent scan for benzene"]) {
        expect(text).toContain(flat(words))
      }
    })
  }
})
