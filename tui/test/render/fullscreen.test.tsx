/**
 * The window full screen (the default): the conversation in the window,
 * scrolling inside it -- the wheel, PgUp/PgDn, ctrl+home/end -- with the
 * composer fixed at the bottom; scrolled up it stays put while the session
 * goes on, at the end it follows; a drag selects and is copied; the agent
 * view takes the same screen; a session put in view starts the list
 * again; a resize re-wraps without a replay.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { FIXTURE, mountWindow } from "../harness"

const PAGE_UP = "\x1b[5~"
const PAGE_DOWN = "\x1b[6~"
const CTRL_HOME = "\x1b[1;5H"
const CTRL_END = "\x1b[1;5F"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

async function conversation(size: { width?: number; height?: number } = {}) {
  w = await mountWindow({ full: true, width: size.width ?? 100, height: size.height ?? 30 })
  w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, ...FIXTURE.session } as never)
  for (const e of FIXTURE.events) w.store.push(e)
  await w.settle()
  return w
}
const rows = () => w!.screenLines().slice(0, -1)          // the frame ends with a newline
const shown = () => rows().join("\n")
const key = async (sequence: string) => {
  w!.keys.pressKey(sequence)
  await Bun.sleep(20)
  await w!.settle()
}
let seq = 1000
const turn = (text: string) => {
  const s = ++seq
  w!.store.push({ type: "user.message", session: "s1", ts: 0, seq: s, text, source: "user" } as never)
  w!.store.push({ type: "turn.started", session: "s1", ts: 0, seq: s + 0.1, turnId: `t${s}`, source: "user" } as never)
  w!.store.push({ type: "text.ended", session: "s1", ts: 0, seq: s + 0.2, partId: `p${s}`, text: `answer to ${text}` } as never)
  w!.store.push({ type: "turn.ended", session: "s1", ts: 0, seq: s + 0.3, turnId: `t${s}`, durationMs: 100 } as never)
}

describe("full screen", () => {
  test("a finished answer block shows its text at once (no blank while it is parsed)", async () => {
    await conversation()
    expect(shown()).toContain("Using SMD.")
  })

  test("the conversation fills the window, the composer at the bottom", async () => {
    await conversation()
    const lines = rows()
    expect(lines.at(-1)).toContain("shift+tab mode")
    expect(lines.slice(-4).some((l) => l.startsWith("❯ "))).toBe(true)
    // The conversation's end just above them, and its rows from the top.
    expect(shown()).toContain("Using SMD.")
    expect(lines[0]!.trim()).not.toBe("")
  })

  test("PgUp scrolls back, PgDn and ctrl+end come back, ctrl+home goes to the start", async () => {
    await conversation({ height: 20 })
    expect(shown()).not.toContain("DELFIN  fake/fake-model")       // the banner is far up
    await key(CTRL_HOME)
    expect(shown()).toContain("DELFIN  fake/fake-model")
    expect(shown()).not.toContain("Using SMD.")
    await key(CTRL_END)
    expect(shown()).toContain("Using SMD.")
    await key(PAGE_UP)
    const up = shown()
    expect(up).not.toContain("Using SMD.")
    await key(PAGE_DOWN)
    await key(PAGE_DOWN)
    expect(shown()).toContain("Using SMD.")
  })

  test("scrolled up it stays put while the session goes on; at the end it follows", async () => {
    await conversation({ height: 20 })
    await key(CTRL_HOME)
    turn("while scrolled up")
    await w!.settle()
    expect(shown()).toContain("DELFIN  fake/fake-model")
    expect(shown()).not.toContain("answer to while scrolled up")
    await key(CTRL_END)
    expect(shown()).toContain("answer to while scrolled up")
    turn("at the end")
    await w!.settle()
    expect(shown()).toContain("answer to at the end")
  })

  test("the wheel scrolls it", async () => {
    await conversation({ height: 20 })
    for (let i = 0; i < 40; i++) await w!.t.mockMouse.scroll(10, 5, "up")
    await w!.settle()
    expect(shown()).toContain("DELFIN  fake/fake-model")
    for (let i = 0; i < 60; i++) await w!.t.mockMouse.scroll(10, 5, "down")
    await w!.settle()
    expect(shown()).toContain("Using SMD.")
  })

  test("a drag selects, and what it selected is copied when it ends -- the prompt keeps the keys", async () => {
    await conversation()
    const copied: string[] = []
    const renderer = w!.t.renderer as unknown as { copyToClipboardOSC52(text: string): boolean }
    renderer.copyToClipboardOSC52 = (text: string) => { copied.push(text); return true }
    const lines = w!.t.captureCharFrame().split("\n")
    const row = lines.findIndex((l) => l.includes("Using SMD."))
    expect(row).toBeGreaterThan(0)
    const col = lines[row]!.indexOf("Using SMD.")
    await w!.t.mockMouse.drag(col, row, col + "Using SMD".length, row)
    await w!.settle()
    expect(copied.length).toBe(1)
    expect(copied[0]).toContain("Using SMD")
    expect(w!.screenLines().join("\n")).toContain("copied ")
    // Typing goes on into the prompt: the drag took no focus.
    await w!.keys.typeText("after the copy")
    await w!.settle()
    expect(rows().some((l) => l.startsWith("❯ after the copy"))).toBe(true)
  })

  test("a click in the conversation takes no keys, and the next ctrl+c does its work", async () => {
    await conversation()
    w!.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 3000, turnId: "t30", source: "user" } as never)
    await w!.settle()
    await w!.t.mockMouse.click(10, 5)
    await Bun.sleep(20)
    await w!.settle()
    await w!.keys.typeText("typed")
    await w!.settle()
    expect(rows().some((l) => l.startsWith("❯ typed"))).toBe(true)
    w!.keys.pressKey("c", { ctrl: true })
    await Bun.sleep(20)
    await w!.settle()
    expect(w!.calls.some((c) => c.method === "session.cancel")).toBe(true)
  })

  test("the agent view takes the screen above the composer, and gives it back", async () => {
    await conversation()
    w!.keys.pressArrow("left")
    await w!.settle()
    expect(shown()).toContain("◆ AGENTS")
    expect(shown()).not.toContain("Using SMD.")
    await w!.escape()
    expect(shown()).not.toContain("◆ AGENTS")
    expect(shown()).toContain("Using SMD.")
  })

  test("a session put in view starts the list again: its own rows, no divider", async () => {
    await conversation()
    w!.store.update((d) => {
      d.sessions.s2 = { ...d.sessions.s1!, id: "s2", title: "the other one" }
    })
    w!.view.show("s2", { divider: "now: the other one" })
    w!.store.push({ type: "user.message", session: "s2", ts: 0, seq: 1, text: "hello from two", source: "user" } as never)
    await w!.settle()
    expect(shown()).toContain("hello from two")
    expect(shown()).not.toContain("Using SMD.")
    expect(shown()).not.toContain("━")
  })

  test("a resize re-wraps the rows without a replay: a tool's lines follow the width both ways", async () => {
    w = await mountWindow({ full: true, width: 60, height: 30 })
    w.screen!.setDetails(true)                                      // what the tool printed, shown
    const long = "x".repeat(40) + " the end of a long result line " + "y".repeat(30)
    w.store.push({ type: "tool.started", session: "s1", ts: 0, seq: 1, turnId: "t1", callId: "c1", name: "bash",
                   input: { command: "cat log" } } as never)
    w.store.push({ type: "tool.ended", session: "s1", ts: 0, seq: 2, turnId: "t1", callId: "c1", name: "bash",
                   ok: true, output: `${long}\nsecond line` } as never)
    await w.settle()
    expect(shown()).not.toContain("the end of a long result")      // cut at 60
    w.t.resize(120, 30)
    await Bun.sleep(300)
    await w.settle()
    expect(shown()).toContain("the end of a long result")           // room at 120
    w.t.resize(60, 30)
    await Bun.sleep(300)
    await w.settle()
    expect(shown()).not.toContain("the end of a long result")
    expect(shown()).toContain("…")                                   // cut again, and says so
    expect(w.view.transcript.replayCount).toBe(0)
  })

  test("the answer streams in its own box under the rows; done, it joins them", async () => {
    await conversation()
    w!.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 2000, turnId: "t9", source: "user" } as never)
    w!.store.push({ type: "text.started", session: "s1", ts: 0, partId: "p9" } as never)
    w!.store.push({ type: "text.delta", session: "s1", ts: 0, partId: "p9", delta: "Streaming words arrive" } as never)
    await w!.settle()
    const lines = rows()
    const at = lines.findIndex((l) => l.includes("Streaming words arrive"))
    const composer = lines.findLastIndex((l) => l.startsWith("❯ "))
    expect(at).toBeGreaterThan(0)
    expect(at).toBeLessThan(composer)
    // Not one of the rows: the rows hold no such item yet.
    expect(w!.screen!.items.some((i) => String(i.text ?? "").includes("Streaming words arrive"))).toBe(false)
    w!.store.push({ type: "text.ended", session: "s1", ts: 0, seq: 2001, partId: "p9", text: "Streaming words arrive." } as never)
    w!.store.push({ type: "turn.ended", session: "s1", ts: 0, seq: 2002, turnId: "t9", durationMs: 10 } as never)
    await w!.settle()
    expect(w!.screen!.items.some((i) => String(i.text ?? "").includes("Streaming words arrive"))).toBe(true)
    expect(shown()).toContain("Streaming words arrive.")
  })

  test("a session put in view while scrolled up is followed from its end", async () => {
    await conversation({ height: 20 })
    await key(CTRL_HOME)
    w!.store.update((d) => { d.sessions.s2 = { ...d.sessions.s1!, id: "s2", title: "the other one" } })
    w!.view.show("s2", { divider: "now: the other one" })
    for (let i = 0; i < 30; i++) {
      w!.store.push({ type: "notice", session: "s2", ts: 0, seq: 10 + i, text: `row number ${i} of the other` } as never)
    }
    await w!.settle()
    expect(shown()).toContain("row number 29 of the other")
  })

  test("a theme that changes after the first frame reaches every row", async () => {
    const { createSignal } = await import("solid-js")
    const { DARK, LIGHT } = await import("../../src/theme/theme")
    const { rgbToHex } = await import("@opentui/core")
    const [theme, setTheme] = createSignal(DARK)
    w = await mountWindow({ full: true, width: 100, height: 30, theme })
    w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, ...FIXTURE.session } as never)
    for (const e of FIXTURE.events) w.store.push(e)
    await w.settle()
    const colourOf = (text: string) => {
      const span = w!.t.captureSpans().lines.flatMap((l) => l.spans).find((s) => s.text.includes(text))
      return span ? rgbToHex(span.fg).toLowerCase() : "(not drawn)"
    }
    const before = colourOf("which solvent")                        // a user row, in view
    setTheme(LIGHT)
    await w.settle()
    const after = colourOf("which solvent")
    expect(after).not.toBe("(not drawn)")
    expect(after).not.toBe(before)
  })

  test("closing the agent view gives the conversation back where it was", async () => {
    await conversation({ height: 20 })
    await key(CTRL_HOME)
    // The conversation and the composer -- not the hints line under them,
    // nor the scroll bar at the right, which say what changed meanwhile
    // (the fixture's session ends: a row more at the end).
    const conversationRows = () => rows().slice(0, -1).map((l) => l.slice(0, -1)).join("\n")
    const before = conversationRows()
    w!.keys.pressArrow("left")
    await w!.settle()
    expect(shown()).toContain("◆ AGENTS")
    await w!.escape()
    expect(conversationRows()).toBe(before)
  })
})
