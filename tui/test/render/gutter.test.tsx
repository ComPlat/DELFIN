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

  test("is the logo's disc in half blocks where the terminal draws colours", async () => {
    const logo = await import("../../src/theme/logo")
    logo.setLogoDepth("rgb")
    try {
      const { scrollback } = await conversation(100)
      const rows = scrollback.filter((l) => /[\u2580\u2584]/.test(l))
      expect(rows.length).toBe(logo.LOGO_ROWS)
      expect(logo.LOGO_ROWS).toBe(15)
      expect(logo.LOGO_WIDTH).toBe(30)
      for (const row of rows) expect(row.startsWith(PAD)).toBe(true)
      // The session's words beside it; no braille.
      expect(rows.some((l) => l.includes("DELFIN  fake/fake-model"))).toBe(true)
      expect(scrollback.some((l) => /[\u2801-\u28ff]/.test(l))).toBe(false)
      // Every run is one glyph, its width the logo's; 24-bit and 256
      // colours cover the same cells.
      for (const d of ["rgb", "256"] as const) {
        const all = logo.logoRuns(d)
        expect(all.length).toBe(logo.LOGO_ROWS)
        all.forEach((runs, y) => {
          expect(runs.reduce((n, r) => n + r.text.length, 0)).toBe(logo.LOGO_WIDTH)
          for (const r of runs) expect(new Set(r.text).size).toBe(1)
          const shape = (rs: readonly { text: string }[]) => rs.map((r) => r.text).join("")
          expect(shape(runs)).toBe(shape(logo.LOGO_RUNS_RGB[y]!))
        })
      }
      // Outside the disc the terminal's own background: no run paints a
      // background behind a half block that has a transparent half.
      for (const runs of logo.LOGO_RUNS_RGB) {
        for (const r of runs) {
          if (r.text[0] === " ") expect([r.fg, r.bg]).toEqual([null, null])
          if (r.text[0] === "\u2584") expect(r.bg).toBeNull()
          if (r.text[0] === "\u2580" && r.bg === null) expect(r.fg).not.toBeNull()
        }
      }
      // The corners are outside the disc.
      expect(logo.LOGO_RUNS_RGB[0]![0]!.text[0]).toBe(" ")
      expect(logo.LOGO_RUNS_RGB[logo.LOGO_ROWS - 1]![0]!.text[0]).toBe(" ")
    } finally {
      logo.setLogoDepth("mono")
    }
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

  test("in 256 colours is the disc's blue, blues of its hue and greys, never lavender", async () => {
    const { LOGO_RUNS_256, LOGO_RUNS_RGB } = await import("../../src/theme/logo")
    const L = [0, 95, 135, 175, 215, 255]
    const count = new Map<number, number>()
    for (const runs of LOGO_RUNS_256) {
      for (const r of runs) {
        for (const c of [r.fg, r.bg]) {
          if (c === null) continue
          expect(typeof c).toBe("number")
          const i = c as number
          expect(i).toBeGreaterThanOrEqual(16)
          expect(i).toBeLessThanOrEqual(255)
          count.set(i, (count.get(i) ?? 0) + r.text.length)
          if (i >= 232 || i === 231 || i === 16) continue
          const [red, green, blue] = [L[Math.floor((i - 16) / 36)]!, L[Math.floor((i - 16) / 6) % 6]!, L[(i - 16) % 6]!]
          // A cube colour is a blue: red below green below blue (lavender
          // has red = green, a green-grey green = blue).
          if (!(red === green && green === blue)) expect(red < green && green < blue).toBe(true)
        }
      }
    }
    // The disc is mostly its blue, 67 (#5f87af); the dolphin white, 231.
    const top = [...count.entries()].sort((a, b) => b[1] - a[1])[0]![0]
    expect(top).toBe(67)
    expect(count.has(231)).toBe(true)
    // In 24-bit the pixels are hex colours.
    for (const runs of LOGO_RUNS_RGB) for (const r of runs) for (const c of [r.fg, r.bg]) {
      if (c !== null) expect(c).toMatch(/^#[0-9a-f]{6}$/)
    }
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
