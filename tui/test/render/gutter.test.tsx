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

  test("draws the dolphin in braille cells, an empty cell a space", async () => {
    const { DOLPHIN, DOLPHIN_WIDTH } = await import("../../src/theme/logo")
    expect(DOLPHIN.length).toBeGreaterThanOrEqual(6)
    expect(DOLPHIN.length).toBeLessThanOrEqual(8)
    expect(DOLPHIN_WIDTH).toBeGreaterThanOrEqual(18)
    expect(DOLPHIN_WIDTH).toBeLessThanOrEqual(22)
    for (const row of DOLPHIN) {
      expect(row.length).toBe(DOLPHIN_WIDTH)
      expect(row).toMatch(/^[ \u2801-\u28ff]+$/)
    }
  })

  const MARK_GLYPHS: Record<string, RegExp> = {
    octant: /^[ \u{1cd00}-\u{1cde5}\u2580-\u2590\u2596-\u259f\u{1cea0}\u{1cea3}\u{1cea8}\u{1ceab}\u{1fb82}\u{1fb85}\u{1fbe6}\u{1fbe7}]+$/u,
    sextant: /^[ \u{1fb00}-\u{1fb3b}\u258c\u2590\u2588]+$/u,
    halfblock: /^[ \u2580\u2584\u2588]+$/u,
  }
  const cols = (text: string) => [...text].length

  for (const cells of ["octant", "sextant", "halfblock"] as const) {
    test(`is the dolphin and its waves in ${cells} cells where the terminal draws colours`, async () => {
      const logo = await import("../../src/theme/logo")
      logo.setLogoDepth("rgb")
      logo.setLogoCells(cells)
      try {
        const mark = logo.MARK[cells]
        expect(logo.LOGO_WIDTH).toBe(40)
        expect(mark.length).toBe(cells === "halfblock" ? 17 : 18)
        for (const line of mark) {
          expect(cols(line)).toBe(logo.LOGO_WIDTH)
          // One column a glyph, as the window and the terminal count them.
          expect(Bun.stringWidth(line)).toBe(logo.LOGO_WIDTH)
          expect(line).toMatch(MARK_GLYPHS[cells]!)
        }
        const { scrollback } = await conversation(100)
        const first = mark.findIndex((l) => l.trim())
        const at = scrollback.findIndex((l) => l.startsWith(PAD + mark[first]!.trimEnd()))
        expect(at).toBeGreaterThanOrEqual(0)
        for (let i = first; i < mark.length; i++) expect(scrollback[at + i - first]!.startsWith(PAD + mark[i]!.trimEnd())).toBe(true)
        // The session's words beside it, where one column a glyph puts them.
        const words = scrollback.find((l) => l.includes("DELFIN  fake/fake-model"))!
        expect(cols(words.slice(0, words.indexOf("DELFIN  fake")))).toBe(PAD.length + logo.LOGO_WIDTH + 3)
        expect(scrollback.some((l) => /[\u2801-\u28ff]/.test(l))).toBe(false)
        // One run a line, in the logo's blue on the terminal's background.
        for (const d of ["rgb", "256"] as const) {
          const runs = logo.logoRuns(d)
          expect(runs.map((r) => r.map((x) => x.text).join(""))).toEqual([...mark])
          for (const r of runs.flat()) expect([r.fg, r.bg]).toEqual([d === "rgb" ? logo.LOGO_RGB : logo.LOGO_256, null])
        }
      } finally {
        logo.setLogoCells("octant")
        logo.setLogoDepth("mono")
      }
    })
  }

  test("maps the 256 octant patterns to 256 glyphs, the 230 of the block in order", async () => {
    const { octant, sextant, halfBlock } = await import("../../src/theme/logo")
    const all = Array.from({ length: 256 }, (_, i) => octant(i))
    expect(new Set(all).size).toBe(256)
    const block = all.filter((g) => g.codePointAt(0)! >= 0x1cd00 && g.codePointAt(0)! <= 0x1cde5)
    expect(block.map((g) => g.codePointAt(0))).toEqual(Array.from({ length: 230 }, (_, i) => 0x1cd00 + i))
    expect(octant(0b100)).toBe("\u{1cd00}")       // BLOCK OCTANT-3
    expect(octant(0b11111101)).toBe("\u{1cde4}")  // BLOCK OCTANT-1345678
    for (const g of all) expect(Bun.stringWidth(g)).toBe(1)
    const six = Array.from({ length: 64 }, (_, i) => sextant(i))
    expect(new Set(six).size).toBe(64)
    for (const g of six) expect(Bun.stringWidth(g)).toBe(1)
    expect([0, 1, 2, 3].map(halfBlock)).toEqual([" ", "\u2580", "\u2584", "\u2588"])
  })

  test("is beside the words in a window of 76 columns or more, the words alone in 75", async () => {
    // The row is the window less the gutter: 40 for the mark, 3 between,
    // 32 for the words (the 30-column disc needed 66).
    const logo = await import("../../src/theme/logo")
    logo.setLogoDepth("rgb")
    try {
      for (const [width, beside] of [[75, false], [76, true]] as const) {
        const { scrollback } = await conversation(width)
        expect(scrollback.some((l) => /[\u{1cd00}-\u{1cde5}]/u.test(l))).toBe(beside)
        expect(scrollback.some((l) => l.startsWith(PAD + "\u25c6 DELFIN"))).toBe(!beside)
        await w?.close()
        w = null
      }
    } finally {
      logo.setLogoDepth("mono")
    }
  })

  test("picks the cells by DELFIN_TUI_LOGO: octants, sextants or half blocks", async () => {
    const { logoCellsFor, logoDepthFor } = await import("../../src/theme/logo")
    expect(logoCellsFor({})).toBe("octant")
    expect(logoCellsFor({ TMUX: "/tmp/tmux-1/default,1,0" })).toBe("octant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "sextant" })).toBe("sextant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "halfblock" })).toBe("halfblock")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "rgb" })).toBe("octant")
    // A kind of cell leaves the depth to the terminal.
    expect(logoDepthFor({ COLORTERM: "truecolor", DELFIN_TUI_LOGO: "sextant" }, "dark")).toBe("rgb")
    expect(logoDepthFor({ TERM: "xterm-256color", DELFIN_TUI_LOGO: "halfblock" }, "dark")).toBe("256")
    expect(logoDepthFor({ TERM: "xterm", DELFIN_TUI_LOGO: "sextant" }, "dark")).toBe("mono")
  })

  test("picks the picture by COLORTERM and TERM, braille under NO_COLOR and the mono theme", async () => {
    const { logoDepthFor } = await import("../../src/theme/logo")
    expect(logoDepthFor({ COLORTERM: "truecolor" }, "dark")).toBe("rgb")
    expect(logoDepthFor({ COLORTERM: "24bit", TERM: "xterm" }, "light")).toBe("rgb")
    expect(logoDepthFor({ TERM: "xterm-256color" }, "dark")).toBe("256")
    expect(logoDepthFor({ TERM: "tmux-256color" }, "dark")).toBe("256")
    expect(logoDepthFor({ TERM: "xterm" }, "dark")).toBe("mono")
    expect(logoDepthFor({}, "dark")).toBe("mono")
    expect(logoDepthFor({ COLORTERM: "truecolor", NO_COLOR: "1" }, "dark")).toBe("mono")
    expect(logoDepthFor({ COLORTERM: "truecolor" }, "mono")).toBe("mono")
    expect(logoDepthFor({ COLORTERM: "truecolor", DELFIN_TUI_LOGO: "mono" }, "dark")).toBe("mono")
    expect(logoDepthFor({ TERM: "xterm", DELFIN_TUI_LOGO: "256" }, "dark")).toBe("256")
    expect(logoDepthFor({ TERM: "xterm", DELFIN_TUI_LOGO: "rgb", NO_COLOR: "1" }, "dark")).toBe("mono")
  })

  test("is the logo's blue: 67 (#5f87af) in 256 colours", async () => {
    const { LOGO_RGB, LOGO_256 } = await import("../../src/theme/logo")
    expect(LOGO_RGB).toMatch(/^#[0-9a-f]{6}$/)
    expect(LOGO_256).toBe(67)
  })

  test("is the words alone, under the mark, where there is not", async () => {
    const { scrollback } = await conversation(40)
    expect(scrollback.some((l) => /[\u2801-\u28ff]/.test(l))).toBe(false)
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
      const flat = (t: string) => t.replace(/[\u2800-\u28ff]/g, "").replace(/\s+/g, "")
      const text = flat(w.t.externalOutput.takeText())
      for (const words of ["3 uncommitted)", "[from ~/.delfin/settings.json]", "writes only inside the workspace (--no-isolate)",
                           "solvent scan for benzene"]) {
        expect(text).toContain(flat(words))
      }
    })
  }
})
