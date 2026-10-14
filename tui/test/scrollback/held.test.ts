/**
 * A streamed answer's finished block between the footer and the rows:
 * held in the session's state until its rows are drawn, then let go; the
 * full screen's list keeps a held block out of sight until it is
 * highlighted (bounded); an inline hand-off waits for the renderer to be
 * between frames (bounded too). What the footer draws of the live blocks
 * and how many rows it gives them.
 */
import { describe, expect, test } from "bun:test"
import { Presenter } from "../../src/scrollback/present"
import { createScreenItems } from "../../src/scrollback/screenview"
import { Transcript, type Item } from "../../src/scrollback/transcript"
import { CommitQueue } from "../../src/scrollback/queue"
import { whenQuiet } from "../../src/scrollback/writers"
import { TAIL_MAX, blockRows, fenceOf, liveBlocks, liveRows, wrapRows } from "../../src/footer/live"
import { newSession } from "../../src/store/state"

const ev = (type: string, fields: Record<string, unknown> = {}) => ({ type, session: "s1", ts: 0, ...fields }) as never

describe("held blocks", () => {
  test("a finished block stays in the footer, by the key it streamed with, until its rows are drawn", async () => {
    const written: Item[] = []
    let draw: () => void = () => {}
    const queue = new CommitQueue({ requestLive() {}, dropLive() {} })
    const transcript = new Transcript(queue, { resetSplitFooterForReplay() {} },
      (item) => item.kind === "markdown" ? new Promise<void>((done) => { draw = () => { written.push(item); done() } }) : undefined)
    const released: string[] = []
    const s = newSession("s1")
    const draft = { sessions: { s1: s } } as never
    const p = new Presenter("s1", transcript, { showThinking: false,
      release: (_sid, key) => { released.push(key); s.held = s.held.filter((h) => h.key !== key) } })
    p.handle([ev("text.started", { partId: "p1" }), ev("text.delta", { partId: "p1", delta: "## Head" })], draft)
    expect(s.tail).toBe("## Head")
    const key0 = s.tailKey
    expect(key0).toMatch(/^\d+\.p1#0$/)
    expect(liveBlocks(s).map((b) => b.key)).toEqual([key0])
    p.handle([ev("text.delta", { partId: "p1", delta: "ing\n\nThe next" })], draft)
    // The heading is final: held under the tail's old key, the tail is the next block.
    expect(s.held).toEqual([{ key: key0, text: "## Heading" }])
    expect(s.tail).toBe("The next")
    const key1 = s.tailKey
    expect(key1).toBe(key0.replace(/#0$/, "#1"))
    expect(liveBlocks(s).map((b) => [b.key, b.text])).toEqual([[key0, "## Heading"], [key1, "The next"]])
    await Bun.sleep(5)
    expect(released).toEqual([])                                  // not drawn yet
    draw()
    await queue.drain()
    expect(written.map((i) => i.text)).toEqual(["## Heading"])
    expect(written[0]!.held).toBe(key0)
    expect(written[0]!.session).toBe("s1")
    expect(released).toEqual([key0])
    expect(liveBlocks(s).map((b) => b.key)).toEqual([key1])
  })

  test("a session shown again streams under keys of its own: never one still held", () => {
    const transcript = { append: () => new Promise<void>(() => {}) }       // never drawn
    const s = newSession("s1")
    const draft = { sessions: { s1: s } } as never
    const first = new Presenter("s1", transcript as never, { showThinking: false })
    first.handle([ev("text.started", { partId: "p1" }), ev("text.delta", { partId: "p1", delta: "One.\n\nTwo" })], draft)
    const again = new Presenter("s1", transcript as never, { showThinking: false })
    again.handle([ev("text.started", { partId: "p1" }), ev("text.delta", { partId: "p1", delta: "Three" })], draft)
    const keys = liveBlocks(s).map((b) => b.key)
    expect(new Set(keys).size).toBe(keys.length)
    expect(liveBlocks(s).map((b) => b.text)).toContain("Three")
  })

  test("a heading's ## without words is not drawn yet", () => {
    const s = newSession("s1")
    s.tail = "Some words.\n\n## "
    s.tailKey = "p1#0"
    expect(liveBlocks(s)).toEqual([{ key: "p1#0", text: "Some words.\n\n" }])
    s.tail = "##"
    expect(liveBlocks(s)).toEqual([])
    s.tail = "## A"
    expect(liveBlocks(s)).toEqual([{ key: "p1#0", text: "## A" }])
  })

  test("a block that was never streamed (a replay) is not held", () => {
    const appended: Item[] = []
    const transcript = { append: async (item: Item) => { appended.push(item) } }
    const s = newSession("s1")
    const p = new Presenter("s1", transcript as never, { showThinking: false })
    p.handle([ev("text.ended", { seq: 3, partId: "p9", text: "## Old\n\nFrom the replay." })], { sessions: { s1: s } } as never)
    expect(appended.map((i) => i.text)).toEqual(["## Old", "From the replay."])
    expect(appended.every((i) => i.held === undefined)).toBe(true)
    expect(s.held).toEqual([])
  })
})

describe("the full screen's list", () => {
  test("a held block is held until revealed; any other row is shown at once", async () => {
    const screen = createScreenItems(600, 5_000)
    let shown = false
    void screen.push({ id: "m1", kind: "markdown", text: "## A", held: "k1" }).then(() => { shown = true })
    await screen.push({ id: "r1", kind: "notice", text: "a row" })
    // A block that did not stream (a replay) has nothing to hand over.
    await screen.push({ id: "m2", kind: "markdown", text: "## Replayed" })
    expect(screen.held("m2")).toBe(false)
    expect(screen.held("m1")).toBe(true)
    expect(screen.held("r1")).toBe(false)
    await Bun.sleep(5)
    expect(shown).toBe(false)
    screen.reveal("m1")
    await Bun.sleep(0)
    expect(shown).toBe(true)
    expect(screen.held("m1")).toBe(false)
  })

  test("a block whose highlighting never comes is shown after the bound", async () => {
    const screen = createScreenItems(600, 40)
    const t0 = Date.now()
    await screen.push({ id: "m1", kind: "markdown", text: "## A", held: "k1" })
    expect(Date.now() - t0).toBeGreaterThanOrEqual(35)
    expect(screen.held("m1")).toBe(false)
  })

  test("a clear lets every held block go", async () => {
    const screen = createScreenItems(600, 5_000)
    const pending = [screen.push({ id: "m1", kind: "markdown", text: "a", held: "k1" }),
                     screen.push({ id: "m2", kind: "markdown", text: "b", held: "k2" })]
    screen.clear()
    await Promise.all(pending)
    expect(screen.items).toEqual([])
    expect(screen.held("m1")).toBe(false)
  })

  test("a block dropped at the cap leaves no hold behind", async () => {
    const screen = createScreenItems(2, 5_000)
    void screen.push({ id: "m1", kind: "markdown", text: "a", held: "k1" })
    await screen.push({ id: "r1", kind: "notice", text: "b" })
    await screen.push({ id: "r2", kind: "notice", text: "c" })
    expect(screen.items.map((i) => i.id)).toEqual(["r1", "r2"])
    expect(screen.held("m1")).toBe(false)
  })
})

describe("the live blocks' rows", () => {
  test("a held block and the tail get their rows and the one between", () => {
    const s = newSession("s1")
    s.held = [{ key: "a", text: "## Summary of the run" }]
    s.tail = "The first words"
    s.tailKey = "b"
    expect(liveRows(s, 80)).toBe(3)                      // heading, the row between, the tail
  })

  test("a long paragraph counts its wrapped rows; a fence its code lines", () => {
    expect(wrapRows("word ".repeat(40).trim(), 50)).toBe(4)           // 199 columns at 50
    expect(wrapRows("x".repeat(120), 50)).toBe(3)
    expect(blockRows("**bold** and `code` and [a link](https://example.org/a/long/target)", 80)).toBe(1)
    expect(blockRows("```python\ndef f(x):\n    return x\n```", 80)).toBe(2)
    expect(blockRows("```python\ndef f(x):", 80)).toBe(1)
    expect(blockRows("- one\n- two\n- three", 80)).toBe(3)
  })

  test("tables, indents, tabs, quotes and inline backticks get the rows they are drawn in", () => {
    expect(blockRows("| Name | Value |\n|---|---|\n| alpha | 1 |", 80)).toBe(5)       // borders and a rule
    expect(blockRows("def f():\n" + " ".repeat(8) + "x".repeat(72), 80)).toBe(3)
    expect(wrapRows("\t" + "y".repeat(76), 79)).toBe(2)
    expect(blockRows("> " + "w".repeat(78), 80)).toBe(2)
    expect(blockRows("```inline``` is how you write it", 80)).toBe(1)
  })

  test("only the trailing blocks that fill the footer are drawn", () => {
    const s = newSession("s1")
    s.held = Array.from({ length: 20 }, (_, i) => ({ key: `k${i}`, text: `block number ${i}` }))
    const shown = liveBlocks(s, 80)
    expect(shown.at(-1)!.key).toBe("k19")
    expect(shown.length).toBeLessThanOrEqual(Math.ceil(TAIL_MAX / 2) + 1)
    expect(liveRows(s, 80)).toBe(TAIL_MAX)
  })

  test("a fenced block is drawn as code: its lines, without the fences", () => {
    expect(fenceOf("```python\ndef f(x):\n    return x\n```\n")).toEqual({ lang: "python", code: "def f(x):\n    return x" })
    expect(fenceOf("~~~\nplain\n")).toEqual({ lang: "", code: "plain" })
    expect(fenceOf("```py")).toEqual({ lang: "py", code: "" })
    expect(fenceOf("Some ```inline``` text")).toBeNull()
    expect(fenceOf("```a`b\nx")).toBeNull()
  })
})

describe("the room an open question gets", () => {
  test("is less the rows of a block on its way into the rows, so the panel does not push it out", async () => {
    const { approvalBudget } = await import("../../src/footer/footer")
    const s = newSession("s1")
    const without = approvalBudget(s, { width: 80, height: 30 }, "")
    s.held = [{ key: "k", text: "Before the question, two lines\nof text." }]
    expect(approvalBudget(s, { width: 80, height: 30 }, "")).toBe(without - 2)
  })
})

describe("whenQuiet", () => {
  test("waits while a frame is drawn or rows wait to be written", async () => {
    const target = { rendering: true, externalOutputQueue: { size: 1 } }
    let quiet = false
    const done = whenQuiet(target, 500).then(() => { quiet = true })
    await Bun.sleep(10)
    expect(quiet).toBe(false)
    target.rendering = false
    await Bun.sleep(10)
    expect(quiet).toBe(false)                             // rows still wait
    target.externalOutputQueue.size = 0
    await done
    expect(quiet).toBe(true)
  })

  test("waits for nothing while the terminal is handed over (no frame comes)", async () => {
    const t0 = Date.now()
    await whenQuiet({ rendering: false, externalOutputQueue: { size: 3 }, controlState: "explicit_suspended" }, 500)
    expect(Date.now() - t0).toBeLessThan(50)
  })

  test("goes anyway when no frame comes", async () => {
    const t0 = Date.now()
    await whenQuiet({ rendering: true }, 40)
    expect(Date.now() - t0).toBeGreaterThanOrEqual(35)
  })
})
