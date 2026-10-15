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

async function conversation(width: number, height = 30) {
  w = await mountWindow({ width, height })
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
    octant: /^[ \u{1cd00}-\u{1cde5}▀-▐▖-▟\u{1cea0}\u{1cea3}\u{1cea8}\u{1ceab}\u{1fb82}\u{1fb85}\u{1fbe6}\u{1fbe7}]+$/u,
    sextant: /^[ \u{1fb00}-\u{1fb3b}▌▐█]+$/u,
    halfblock: /^[ ▀▄█]+$/u,
  }
  const cols = (text: string) => [...text].length
  // A terminal of these columns x rows draws this size (the window less the
  // gutter is the banner's width).
  const TERMINAL = { L: [100, 45], M: [100, 30], S: [100, 22] } as const
  const SHAPE = { L: [30, 13], M: [20, 9], S: [16, 7] } as const

  for (const cells of ["sextant", "octant", "halfblock"] as const) {
    for (const size of ["L", "M", "S"] as const) {
      test(`is the dolphin and its waves in ${cells} cells, size ${size}, where the terminal draws colours`, async () => {
        const logo = await import("../../src/theme/logo")
        logo.setLogoDepth("rgb")
        logo.setLogoCells(cells)
        try {
          const mark = logo.MARK[cells][size]
          const [width, rows] = SHAPE[size]
          expect(logo.pickedWidth(size)).toBe(width)
          expect(mark.length).toBe(rows)
          for (const line of mark) {
            expect(cols(line)).toBe(width)
            // One column a glyph, as the window and the terminal count them.
            expect(Bun.stringWidth(line)).toBe(width)
            expect(line).toMatch(MARK_GLYPHS[cells]!)
          }
          const [tw, th] = TERMINAL[size]
          expect(logo.pickLogo("rgb", tw - GUTTER, th)).toBe(size)
          const { scrollback } = await conversation(tw, th)
          const first = mark.findIndex((l) => l.trim())
          const at = scrollback.findIndex((l) => l.startsWith(PAD + mark[first]!.trimEnd()))
          expect(at).toBeGreaterThanOrEqual(0)
          for (let i = first; i < mark.length; i++) expect(scrollback[at + i - first]!.startsWith(PAD + mark[i]!.trimEnd())).toBe(true)
          // The session's words beside it, where one column a glyph puts them.
          const words = scrollback.find((l) => l.includes("DELFIN  fake/fake-model"))!
          expect(cols(words.slice(0, words.indexOf("DELFIN  fake")))).toBe(PAD.length + width + logo.LOGO_GAP)
          expect(scrollback.some((l) => /[⠁-⣿]/.test(l))).toBe(false)
          // One run a line, in the logo's blue on the terminal's background.
          for (const d of ["rgb", "256"] as const) {
            const runs = logo.logoRuns(d, size)
            expect(runs.map((r) => r.map((x) => x.text).join(""))).toEqual([...mark])
            for (const r of runs.flat()) expect([r.fg, r.bg]).toEqual([d === "rgb" ? logo.LOGO_RGB : logo.LOGO_256, null])
          }
        } finally {
          logo.setLogoCells("sextant")
          logo.setLogoDepth("mono")
        }
      })
    }
  }

  test("draws every glyph of every size and kind one column wide", async () => {
    const { MARK } = await import("../../src/theme/logo")
    const used = new Set<string>()
    for (const sizes of Object.values(MARK)) for (const lines of Object.values(sizes)) for (const line of lines) for (const g of line) used.add(g)
    expect(used.size).toBeGreaterThan(40)
    for (const g of used) expect([g, Bun.stringWidth(g)]).toEqual([g, 1])
  })

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
    expect([0, 1, 2, 3].map(halfBlock)).toEqual([" ", "▀", "▄", "█"])
  })

  test("takes the largest size the terminal's rows and width allow", async () => {
    // The banner's width is the window less the gutter: a size needs its
    // columns, 3 between, 32 for the words; and 40, 28 or 20 rows.
    const { pickLogo } = await import("../../src/theme/logo")
    const at = (columns: number, rows: number) => pickLogo("rgb", columns - GUTTER, rows)
    // Rows, with width to spare.
    expect([at(200, 39), at(200, 40)]).toEqual(["M", "L"])
    expect([at(200, 27), at(200, 28)]).toEqual(["S", "M"])
    expect([at(200, 19), at(200, 20)]).toEqual([null, "S"])
    // Width, one column below and at each threshold (66, 56, 52 columns).
    expect([at(65, 45), at(66, 45)]).toEqual(["M", "L"])
    expect([at(55, 45), at(56, 45)]).toEqual(["S", "M"])
    expect([at(51, 45), at(52, 45)]).toEqual([null, "S"])
    expect([at(55, 30), at(56, 30)]).toEqual(["S", "M"])
    expect([at(51, 22), at(52, 22)]).toEqual([null, "S"])
    // Never more than about a third of the rows.
    for (let rows = 1; rows < 120; rows++) {
      const p = at(300, rows)
      const tall = p === "L" ? 13 : p === "M" ? 9 : p === "S" ? 7 : 0
      expect(tall / rows).toBeLessThanOrEqual(0.35)
    }
    // Braille (mono): 8 rows in 24, 19 columns.
    expect([pickLogo("mono", 200, 23), pickLogo("mono", 200, 24)]).toEqual([null, "braille"])
    expect([pickLogo("mono", 53, 30), pickLogo("mono", 54, 30)]).toEqual([null, "braille"])
  })

  for (const [width, height, size] of [[66, 45, "L"], [65, 45, "M"], [100, 40, "L"], [100, 39, "M"],
                                       [100, 28, "M"], [100, 27, "S"], [52, 22, "S"], [51, 22, null], [100, 19, null]] as const) {
    test(`is ${size ?? "the words alone"} in a terminal of ${width} x ${height}`, async () => {
      const logo = await import("../../src/theme/logo")
      logo.setLogoDepth("rgb")
      try {
        const { scrollback } = await conversation(width, height)
        const sextants = scrollback.filter((l) => /[\u{1fb00}-\u{1fb3b}]/u.test(l))
        expect(scrollback.some((l) => l.startsWith(PAD + "◆ DELFIN"))).toBe(size === null)
        if (size === null) {
          expect(sextants.length).toBe(0)
        } else {
          const mark = logo.MARK.sextant[size]
          expect(scrollback.some((l) => l.startsWith(PAD + mark[0]!.trimEnd()))).toBe(true)
          const words = scrollback.find((l) => l.includes("DELFIN  fake/fake-model"))!
          expect(cols(words.slice(0, words.indexOf("DELFIN  fake")))).toBe(PAD.length + logo.pickedWidth(size) + logo.LOGO_GAP)
        }
      } finally {
        logo.setLogoDepth("mono")
      }
    })
  }

  test("full screen: picks the size again when the terminal is resized", async () => {
    const logo = await import("../../src/theme/logo")
    logo.setLogoDepth("rgb")
    try {
      w = await mountWindow({ full: true, width: 100, height: 45 })
      w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, ...FIXTURE.session } as never)
      await w.settle()
      const shows = (size: "L" | "M" | "S") => {
        const frame = w!.t.captureCharFrame()
        return logo.MARK.sextant[size].filter((l) => l.trim()).every((l) => frame.includes(l.trimEnd()))
      }
      expect([shows("L"), shows("M")]).toEqual([true, false])
      w.t.resize(100, 30)
      await w.settle()
      expect([shows("L"), shows("M")]).toEqual([false, true])
      w.t.resize(100, 19)
      await w.settle()
      expect(["L", "M", "S"].some((s) => shows(s as "L"))).toBe(false)
      expect(w.t.captureCharFrame()).toContain("\u25c6 DELFIN")
    } finally {
      logo.setLogoDepth("mono")
    }
  })

  test("is sextant cells unless DELFIN_TUI_LOGO asks for octants or half blocks", async () => {
    const { logoCellsFor, logoDepthFor } = await import("../../src/theme/logo")
    expect(logoCellsFor({})).toBe("sextant")
    // No guess from the terminal's name: it does not say what its font has.
    expect(logoCellsFor({ TERM_PROGRAM: "vscode" })).toBe("sextant")
    expect(logoCellsFor({ KITTY_WINDOW_ID: "1", TERM: "xterm-kitty" })).toBe("sextant")
    expect(logoCellsFor({ TERM_PROGRAM: "WezTerm" })).toBe("sextant")
    expect(logoCellsFor({ TMUX: "/tmp/tmux-1/default,1,0" })).toBe("sextant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "octant" })).toBe("octant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "octant", TERM_PROGRAM: "vscode" })).toBe("octant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "sextant" })).toBe("sextant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "halfblock" })).toBe("halfblock")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "rgb" })).toBe("sextant")
    expect(logoCellsFor({ DELFIN_TUI_LOGO: "quadrant" })).toBe("sextant")
    // A kind of cell leaves the depth to the terminal.
    expect(logoDepthFor({ COLORTERM: "truecolor", DELFIN_TUI_LOGO: "octant" }, "dark")).toBe("rgb")
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
