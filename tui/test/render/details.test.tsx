/**
 * A tool call is its headline and one line saying what came back; what it
 * printed, or an edit's diff, only when asked for. ctrl+o shows the
 * details of every tool call -- on the full screen in place, the reader's
 * place kept; inline by writing the session again -- and hides them
 * again; a click on one (full screen) shows its own. Under GNU screen,
 * which cannot write the session again, ctrl+o shows the last result in
 * full instead.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { DETAILS_SETTLE_MS } from "../../src/session"
import { FIXTURE, mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const CTRL_HOME = "\x1b[1;5H"
let committed = ""

async function conversation(options: { full: boolean; height?: number; divider?: boolean }) {
  committed = ""
  w = await mountWindow({ full: options.full, width: 100, height: options.height ?? 40, divider: options.divider })
  w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, ...FIXTURE.session } as never)
  for (const e of FIXTURE.events) w.store.push(e)
  await w.settle()
  return w
}
/** What the window shows: inline, the scrollback written so far (a
 *  replay wipes it: what it wrote again is all there is) and the footer. */
function shown(): string {
  const external = (w!.t as unknown as { externalOutput?: { takeText(): string } }).externalOutput
  if (external && !w!.screen) committed += external.takeText()
  return committed + "\n" + w!.t.captureCharFrame()
}
const ctrlO = async () => {
  w!.keys.pressKey("o", { ctrl: true })
  await Bun.sleep(20)
  await w!.settle()
}

const DETAIL = "# tui_demo"                                 // read_file README.md printed it
const DIFF = "_ENERGY_RE.search(line)"                     // edit_file's diff

for (const full of [true, false]) {
  describe(full ? "full screen" : "inline", () => {
    test("a tool call shows its headline and what came back, nothing under it", async () => {
      await conversation({ full })
      const text = shown()
      expect(text).toContain("⏺ read_file  README.md")
      expect(text).toContain("⎿ 1 line, 11 B")
      expect(text).toContain("⎿ Edited delfin/xtb.py (+4 -1)")
      expect(text).not.toContain(DETAIL)
      expect(text).not.toContain(DIFF)
    })

    if (full) {
      test("ctrl+o shows every tool call's details, and again hides them", async () => {
        await conversation({ full })
        await ctrlO()
        let text = shown()
        expect(text).toContain(DETAIL)
        expect(text).toContain(DIFF)
        expect(text).toContain("⎿ 1 line, 11 B")
        await ctrlO()
        text = shown()
        expect(text).not.toContain(DETAIL)
        expect(text).not.toContain(DIFF)
        expect(text).toContain("⏺ edit_file  delfin/xtb.py")
      })
    } else {
      // The session written again is the terminal's (a replay: the test
      // terminal has none; check_tui_tmux.py does it live).
      test("ctrl+o writes the session again with the details, and rows after it have theirs", async () => {
        await conversation({ full })
        const replays = w!.view.transcript.replayCount
        await ctrlO()
        expect(w!.view.details()).toBe(true)
        await Bun.sleep(DETAILS_SETTLE_MS + 50)                // written once the presses stop
        await w!.settle()
        expect(w!.view.transcript.replayCount).toBe(replays + 1)
        committed = ""
        w!.store.push({ type: "tool.started", session: "s1", ts: 0, seq: 9001, turnId: "t9", callId: "c9", name: "bash",
                        input: { command: "ls" } } as never)
        w!.store.push({ type: "tool.ended", session: "s1", ts: 0, seq: 9002, turnId: "t9", callId: "c9", name: "bash",
                        ok: true, output: "molecule.xyz\nCONTROL.txt" } as never)
        await w!.settle()
        expect(shown()).toContain("CONTROL.txt")
        await ctrlO()
        expect(w!.view.details()).toBe(false)
      })

      test("ctrl+o pressed again before the session is written takes it back: nothing written", async () => {
        await conversation({ full })
        const replays = w!.view.transcript.replayCount
        await ctrlO()
        await ctrlO()
        await Bun.sleep(DETAILS_SETTLE_MS + 50)
        await w!.settle()
        expect(w!.view.transcript.replayCount).toBe(replays)
        expect(w!.view.details()).toBe(false)
      })
    }
  })
}

describe("full screen", () => {
  test("a click on a tool call shows its own details, and again hides them", async () => {
    await conversation({ full: true })
    const lines = () => w!.t.captureCharFrame().split("\n")
    const row = lines().findIndex((l) => l.includes("⏺ edit_file"))
    expect(row).toBeGreaterThan(0)
    await w!.t.mockMouse.click(6, row)
    await Bun.sleep(20)
    await w!.settle()
    expect(shown()).toContain(DIFF)
    expect(shown()).not.toContain(DETAIL)                   // the other one stays as it was
    const again = lines().findIndex((l) => l.includes("⏺ edit_file"))
    await w!.t.mockMouse.click(6, again)
    await Bun.sleep(20)
    await w!.settle()
    expect(shown()).not.toContain(DIFF)
  })

  test("a drag over a tool call selects and copies; it does not open it", async () => {
    await conversation({ full: true })
    const copied: string[] = []
    const renderer = w!.t.renderer as unknown as { copyToClipboardOSC52(text: string): boolean }
    renderer.copyToClipboardOSC52 = (text: string) => { copied.push(text); return true }
    const lines = w!.t.captureCharFrame().split("\n")
    const row = lines.findIndex((l) => l.includes("⏺ edit_file"))
    const col = lines[row]!.indexOf("edit_file")
    await w!.t.mockMouse.drag(col, row, col + "edit_file".length, row)
    await Bun.sleep(20)
    await w!.settle()
    expect(copied.join("")).toContain("edit_file")
    expect(shown()).not.toContain(DIFF)
  })

  const topRow = () => w!.t.captureCharFrame().split("\n")[0]!.trim()
  /** Down from the start until `text` is the top row. */
  async function scrollTo(text: string) {
    w!.keys.pressKey(CTRL_HOME)
    await Bun.sleep(20)
    await w!.settle()
    for (let i = 0; i < 80 && !topRow().includes(text); i++) {
      await w!.t.mockMouse.scroll(10, 5, "down")
      await w!.settle()
    }
    expect(topRow()).toContain(text)
  }

  test("scrolled up, ctrl+o keeps the row at the top where it was -- with a tool call above it that grows", async () => {
    await conversation({ full: true, height: 16 })
    // Below the first read_file, whose "# tui_demo" opens above it.
    await scrollTo("Summary")
    await ctrlO()
    expect(topRow()).toContain("Summary")
    await ctrlO()
    expect(topRow()).toContain("Summary")
  })

  test("scrolled into an opened diff, ctrl+o keeps the row under it in place, never jumps past it", async () => {
    await conversation({ full: true, height: 16 })
    await ctrlO()
    await scrollTo("42          if")                        // inside the edit's diff
    const rows = () => w!.t.captureCharFrame().split("\n").map((l) => l.trim())
    const under = rows().findIndex((l) => l.includes("The parser now reads"))
    expect(under).toBeGreaterThan(0)
    await ctrlO()                                           // the diff closes above it
    expect(rows().findIndex((l) => l.includes("The parser now reads"))).toBe(under)
  })

  test("a click on a tool call at the end keeps its headline in view; its details open under it", async () => {
    await conversation({ full: true, height: 16 })
    const rows = () => w!.t.captureCharFrame().split("\n")
    const row = rows().findIndex((l) => l.includes("⏺ edit_file"))
    expect(row).toBeGreaterThanOrEqual(0)
    await w!.t.mockMouse.click(6, row)
    await Bun.sleep(20)
    await w!.settle()
    expect(rows().findIndex((l) => l.includes("⏺ edit_file"))).toBe(row)
    expect(shown()).toContain("@@ -41,5 +41,8 @@")
  })

  test("only a left click on the row opens it: not another button, nor a drag that ends on it", async () => {
    await conversation({ full: true })
    const lines = w!.t.captureCharFrame().split("\n")
    const row = lines.findIndex((l) => l.includes("⏺ edit_file"))
    await w!.t.mockMouse.click(6, row, 2)                        // the right button
    await Bun.sleep(20)
    await w!.settle()
    expect(shown()).not.toContain(DIFF)
    // From the gutter (no text: no selection starts) onto the row.
    await w!.t.mockMouse.drag(0, row - 3, 6, row)
    await Bun.sleep(20)
    await w!.settle()
    expect(shown()).not.toContain(DIFF)
  })

  test("at the end, ctrl+o keeps following the end", async () => {
    await conversation({ full: true, height: 16 })
    await ctrlO()
    expect(shown()).toContain("Using SMD.")
  })
})

describe("inline under GNU screen", () => {
  test("ctrl+o shows the last result in full (the session cannot be written again)", async () => {
    await conversation({ full: false, divider: true })
    const replays = w!.view.transcript.replayCount
    await ctrlO()
    await Bun.sleep(DETAILS_SETTLE_MS + 50)
    const text = shown()
    expect(w!.view.transcript.replayCount).toBe(replays)
    expect(text).toContain("edit_file, in full")
  })

  test("the palette offers no details to show there", async () => {
    await conversation({ full: false, divider: true })
    w!.keys.pressKey("p", { ctrl: true })
    await w!.settle()
    await w!.keys.typeText("tool details")
    await w!.settle()
    expect(w!.footer().join("\n")).not.toContain("show tool details")
  })
})
