/**
 * A streamed answer on its way from the live tail into the rows, frame by
 * frame, in both screens: a heading never shows its "## " while it
 * streams, a code block keeps its colours from delta to delta, and a
 * finished block is on the screen exactly once in every frame -- held in
 * the footer (with the tail beside it, whole) while its rows are being
 * drawn, never gone, never there twice, and the footer never left with
 * blank rows. A slow highlighter costs colour, not text. Full screen also:
 * a block that comes while it cannot be seen (scrolled up, the agent view
 * open) is highlighted when it is.
 *
 * A highlighter that holds back its answers for chosen text (`gate`) is
 * what a busy machine does; it keeps a hand-off in its middle for as many
 * frames as a test needs.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { CodeRenderable, getTreeSitterClient, rgbToHex, type Renderable } from "@opentui/core"
import { DARK } from "../../src/theme/theme"
import { FIXTURE, mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  gate?.open()
  await w?.close()
  w = null
  gate?.restore()
  gate = null
})

/** The highlighter every element uses (OpenTUI's one client), but an
 *  answer for text `hold` picks waits until `open`. */
function gated() {
  const real = getTreeSitterClient() as unknown as {
    highlightOnce(content: string, filetype: string): Promise<unknown>
  }
  const own = Object.getOwnPropertyDescriptor(real, "highlightOnce")
  const original = real.highlightOnce.bind(real)
  let holds: ((content: string) => boolean) | null = null
  let waiting: (() => void)[] = []
  real.highlightOnce = async (content: string, filetype: string) => {
    const result = await original(content, filetype)
    while (holds?.(content)) await new Promise<void>((resolve) => waiting.push(resolve))
    return result
  }
  return {
    hold(test: (content: string) => boolean) { holds = test },
    open() {
      holds = null
      const all = waiting
      waiting = []
      for (const go of all) go()
    },
    restore() {
      if (own) Object.defineProperty(real, "highlightOnce", own)
      else delete (real as { highlightOnce?: unknown }).highlightOnce
    },
  }
}
let gate: ReturnType<typeof gated> | null = null

// Above the fixture's: the store drops an event whose seq it has seen.
let seq = 10_000
const ev = (type: string, fields: Record<string, unknown> = {}) =>
  ({ type, session: "s1", ts: 0, ...fields }) as never

async function streaming(full: boolean, size = { width: 80, height: 30 }) {
  committed = ""
  gate = gated()
  w = await mountWindow({ full, ...size })
  w.store.push(ev("session.info", { seq: 0, ...FIXTURE.session }))
  w.store.push(ev("user.message", { seq: ++seq, text: "summarise the run", source: "user" }))
  w.store.push(ev("turn.started", { seq: ++seq, turnId: "t1", source: "user" }))
  w.store.push(ev("text.started", { partId: "p1" }))
  await w.settle()
  return w
}
const delta = (text: string) => w!.store.push(ev("text.delta", { partId: "p1", delta: text }))
/** Inline: the rows committed to the terminal's scrollback so far. */
let committed = ""
/** The next frame as the window draws it: what came in applied, whatever
 *  was waiting (a highlight, a row's commit) given a moment, one render --
 *  and, inline, the scrollback above it. */
async function frame(): Promise<string> {
  w!.store.flush()
  await Bun.sleep(1)
  await w!.t.renderOnce()
  const external = (w!.t as unknown as { externalOutput?: { takeText(): string } }).externalOutput
  if (external) committed += external.takeText()
  return committed + "\n" + w!.t.captureCharFrame()
}
const count = (text: string, part: string) => text.split(part).length - 1
/** Frames until no code block is being highlighted: the highlighter kept
 *  up (it takes a few ms; the next delta comes a frame or more later). */
async function highlighted(): Promise<void> {
  const busy = (r: Renderable): boolean =>
    (r instanceof CodeRenderable && r.isHighlighting) || (r.getChildren() as Renderable[]).some(busy)
  for (let i = 0; i < 200; i++) {
    await w!.t.renderOnce()
    if (!busy(w!.t.renderer.root)) return
    await Bun.sleep(5)
  }
}
/** Inline: the footer's rows as the last frame painted them, top first
 *  (the frame holds the footer alone, at the height it was painted with). */
const footerRows = () => w!.t.captureCharFrame().split("\n").slice(0, -1)

for (const full of [true, false]) {
  describe(full ? "full screen" : "inline", () => {
    test("a heading streaming in never shows its ## (each frame keeps the last highlighted text)", async () => {
      await streaming(full)
      const text = "## Summary of the run\n\nUsing SMD **bold** here, and more words after it."
      const raw: string[] = []
      for (let i = 0; i < text.length; i += 3) {
        delta(text.slice(i, i + 3))
        const f = await frame()                         // the frame right after it
        if (/##/.test(f)) raw.push(`after ${JSON.stringify(text.slice(0, i + 3))}`)
        await highlighted()
      }
      expect(raw).toEqual([])
      await w!.settle()
      const end = await frame()
      expect(end).toContain("Summary of the run")
      expect(end).toContain("Using SMD bold here")
    })

    test("a code block streaming in keeps its colours from delta to delta", async () => {
      await streaming(full)
      // TypeScript: a grammar OpenTUI ships (the window adds none).
      delta("Here:\n\n```typescript\nfunction area(r: number) {\n  return 3.14 * r * r\n")
      await w!.settle()
      await highlighted()
      // "return": a keyword no other row on the screen has.
      const colourOf = (word: string) => {
        const span = w!.t.captureSpans().lines.flatMap((l) => l.spans).find((s) => s.text.includes(word))
        return span ? rgbToHex(span.fg).toLowerCase() : "(not drawn)"
      }
      const highlightedColour = colourOf("return")
      expect(highlightedColour).not.toBe("(not drawn)")
      expect(highlightedColour).not.toBe(String(DARK.text).toLowerCase())   // highlighted at all
      const seen: string[] = []
      for (const piece of ["}\n", "console.log(", "area(2))\n", "const x = area", "(3)\n"]) {
        delta(piece)
        await frame()                                            // the frame right after it
        seen.push(colourOf("return"))
        await highlighted()
      }
      expect(seen).toEqual(seen.map(() => highlightedColour))
    })

    test("a finished block is on the screen once in every frame of its way into the rows", async () => {
      await streaming(full)
      delta("## Summary of the run")
      await w!.settle()
      await highlighted()
      // Its rows cannot be drawn yet (the highlighter is busy with them):
      // the footer keeps it, the next block streaming under it.
      gate!.hold((content) => content.includes("Summary of the run"))
      delta("\n\nThe first words")
      const seen: number[] = []
      // Until the next block shows (its own first highlighting) -- the
      // finished one still held: both whole in the footer.
      let both = false
      for (let i = 0; i < 300 && !both; i++) {
        const f = await frame()
        seen.push(count(f, "Summary of the run") + (f.includes("## Summary") ? 100 : 0))
        both = f.includes("The first words")
        if (!both) await Bun.sleep(5)
      }
      expect(both).toBe(true)
      expect(w!.store.state.sessions.s1!.held.map((h) => h.text)).toEqual(["## Summary of the run"])
      if (full) expect(w!.screen!.held(w!.screen!.items.at(-1)!.id)).toBe(true)
      else expect(committed).not.toContain("Summary of the run")
      for (let i = 0; i < 5; i++) seen.push(count(await frame(), "Summary of the run"))
      gate!.open()
      const blank: number[] = []
      for (let i = 0; i < 20; i++) {
        const f = await frame()
        seen.push(count(f, "Summary of the run") + (f.includes("## Summary") ? 100 : 0))
        if (!full && !footerRows()[0]!.trim()) blank.push(i)
      }
      expect(seen).toEqual(seen.map(() => 1))
      expect(blank).toEqual([])                                  // the footer never left with blank rows
      expect(w!.store.state.sessions.s1!.held).toEqual([])
      if (full) expect(w!.screen!.items.some((i) => i.text === "## Summary of the run")).toBe(true)
      else expect(committed).toContain("Summary of the run")
      // The answer ends: its last block goes the same way -- on the screen
      // once in every frame (its new words a frame or two after the rest:
      // shown when highlighted, as every delta's are).
      w!.store.push(ev("text.ended", { seq: ++seq, partId: "p1",
        text: "## Summary of the run\n\nThe first words, and a second one starts here." }))
      w!.store.push(ev("turn.ended", { seq: ++seq, turnId: "t1", durationMs: 10 }))
      const last: number[] = []
      let end = ""
      for (let i = 0; i < 40; i++) {
        end = await frame()
        last.push(count(end, "The first words"))
        await Bun.sleep(3)
      }
      expect(last).toEqual(last.map(() => 1))
      expect(count(end, "a second one starts here")).toBe(1)
      expect(w!.store.state.sessions.s1!.held).toEqual([])
    })

    test("a block on its way stays on the screen while a question opens", async () => {
      await streaming(full)
      delta("Before the question.")
      await w!.settle()
      await highlighted()
      gate!.hold((content) => content.includes("Before the question"))
      w!.store.push(ev("text.ended", { seq: ++seq, partId: "p1", text: "Before the question." }))
      w!.store.push(ev("permission.asked", { seq: ++seq, requestId: "q1", kind: "confirm", tool: "bash",
        preview: "$ ls", options: [{ key: "y", label: "yes", detail: "this one action" },
                                   { key: "n", label: "no", detail: "refused" }] }))
      const seen: number[] = []
      for (let i = 0; i < 8; i++) seen.push(count(await frame(), "Before the question"))
      gate!.open()
      for (let i = 0; i < 20; i++) seen.push(count(await frame(), "Before the question"))
      expect(seen).toEqual(seen.map(() => 1))
    })

    test("a live block whose colours are late shows its words plain meanwhile, then keeps its colours", async () => {
      await streaming(full)
      gate!.hold((content) => content.includes("Streaming words arrive"))
      delta("Streaming words arrive here")
      const t0 = Date.now()
      let shownAt = -1
      for (let i = 0; i < 80 && shownAt < 0; i++) {
        if ((await frame()).includes("Streaming words arrive here")) shownAt = Date.now() - t0
        else await Bun.sleep(10)
      }
      expect(shownAt).toBeGreaterThanOrEqual(0)                  // never blank for long
      expect(shownAt).toBeLessThan(800)
      gate!.open()
      await highlighted()
      // From its first colours on, a delta keeps the last highlighted text.
      delta(" and **more** after")
      expect(await frame()).not.toContain("**")
      await highlighted()
      expect(await frame()).toContain("Streaming words arrive here and more after")
    })

    test("a highlighter that does not answer costs the block its colours, not its words", async () => {
      await streaming(full)
      delta("Words before the stall.")
      await w!.settle()
      await highlighted()
      gate!.hold((content) => content.includes("Words before the stall"))
      w!.store.push(ev("text.ended", { seq: ++seq, partId: "p1", text: "Words before the stall." }))
      w!.store.push(ev("turn.ended", { seq: ++seq, turnId: "t1", durationMs: 10 }))
      const seen: number[] = []
      const until = Date.now() + 2600
      while (Date.now() < until) {
        seen.push(count(await frame(), "Words before the stall"))
        await Bun.sleep(40)
      }
      expect(seen).toEqual(seen.map(() => 1))
      expect(w!.store.state.sessions.s1!.held).toEqual([])        // in the rows, plain
    }, 15_000)
  })
}

/** The code blocks drawing `text`'s rows, and whether each is highlighted
 *  already -- so the first frame that shows them has their colours. */
function blockState(text: string): boolean[] {
  const found: boolean[] = []
  const visit = (r: Renderable) => {
    if (r instanceof CodeRenderable && r.content.includes(text)) {
      const dirty = (r as unknown as { _highlightsDirty: boolean })._highlightsDirty
      found.push(!dirty && !r.isHighlighting)
    }
    for (const child of r.getChildren() as Renderable[]) visit(child)
  }
  visit(w!.t.renderer.root)
  return found
}

describe("full screen, out of sight", () => {
  const answer = (id: string, text: string) => {
    w!.store.push(ev("turn.started", { seq: ++seq, turnId: `t${id}`, source: "user" }))
    w!.store.push(ev("text.ended", { seq: ++seq, partId: `p${id}`, text }))
    w!.store.push(ev("turn.ended", { seq: ++seq, turnId: `t${id}`, durationMs: 10 }))
  }
  const firstFrame = async () => {
    await Bun.sleep(40)                                   // a lone Esc is held 20 ms
    w!.store.flush()
    await w!.t.renderOnce()
    return w!.t.captureCharFrame()
  }

  test("a block that came while the reader was scrolled up is highlighted when they come back", async () => {
    await streaming(true, { width: 80, height: 16 })
    for (const e of FIXTURE.events) w!.store.push(e)
    // Far more rows than the view holds: the end is well out of it.
    for (let i = 0; i < 60; i++) w!.store.push(ev("notice", { seq: ++seq, text: `filler row ${i}` }))
    await w!.settle()
    w!.keys.pressKey("\x1b[1;5H")                         // ctrl+home
    await w!.settle()
    answer("9", "## Late heading\n\nA paragraph under it, with `code`.")
    await w!.settle()
    await Bun.sleep(30)
    expect(w!.t.captureCharFrame()).not.toContain("Late heading")   // still scrolled up
    // Below the view (culled: no frame updates it) and highlighted all the same.
    expect(blockState("Late heading")).toEqual([true])
    w!.keys.pressKey("\x1b[1;5F")                         // ctrl+end
    const f = await firstFrame()
    expect(f).toContain("Late heading")
    expect(f).not.toContain("## Late")
    expect(f).not.toContain("`code`")
  })

  test("a block that came while the agent view was open is highlighted when it closes", async () => {
    await streaming(true)
    for (const e of FIXTURE.events) w!.store.push(e)
    await w!.settle()
    w!.keys.pressArrow("left")
    await w!.settle()
    expect(w!.t.captureCharFrame()).toContain("◆ AGENTS")
    answer("8", "## Heading while away\n\nWritten while the agent view was open.")
    await w!.settle()
    await Bun.sleep(30)
    // Hidden (no frame updates it) and highlighted all the same.
    expect(blockState("Heading while away")).toEqual([true])
    w!.t.mockInput.pressEscape()
    const f = await firstFrame()
    expect(f).not.toContain("◆ AGENTS")
    expect(f).toContain("Heading while away")
    expect(f).not.toContain("## Heading")
  })

  test("a replay asks for the newest blocks' colours first: the end is what is in view", async () => {
    await streaming(true)
    const asked: string[] = []
    const real = getTreeSitterClient() as unknown as { highlightOnce(c: string, f: string): Promise<unknown> }
    const wrapped = real.highlightOnce
    real.highlightOnce = (content: string, filetype: string) => {
      const n = /answer number (\d+)/.exec(content)?.[1]
      if (n) asked.push(n)
      return wrapped.call(real, content, filetype)
    }
    for (let n = 1; n <= 6; n++) answer(`r${n}`, `The answer number ${n}, replayed.`)
    await w!.settle()
    await Bun.sleep(30)
    expect(asked.slice(0, 6)).toEqual(["6", "5", "4", "3", "2", "1"])
  })

  test("a replayed answer does not wait for its highlighting: the rows after it come at once", async () => {
    await streaming(true)
    gate!.hold((content) => content.includes("Replayed heading"))
    answer("7", "## Replayed heading\n\nIts paragraph.")
    w!.store.push(ev("notice", { seq: ++seq, text: "a row after the answer" }))
    await w!.settle()
    expect(w!.t.captureCharFrame()).toContain("a row after the answer")
  })
})
