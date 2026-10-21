/**
 * The other sessions, seen from the session in view: the hints line says
 * how many more run in the agent service and how many of them wait for
 * the user ("2 more sessions (1 waits for you) ← agents"), from the
 * summaries every window hears; nothing when there are none. Sleeping and
 * ended sessions do not run. Where the line is short it says less. A click
 * on it (full screen) opens the agent view, as ← does. Another session
 * that starts to wait says so once in the message row.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { othersHint } from "../../src/footer/footer"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
})

const now = Date.now() / 1000

async function mount(options: Parameters<typeof mountWindow>[0] = {}) {
  w = await mountWindow({ width: 110, height: 30, arrival: { shownMs: 5000, quietMs: 0 }, ...options })
  ;(w.t.renderer as unknown as { copyToClipboardOSC52(text: string): boolean }).copyToClipboardOSC52 = () => true
  w.store.push({ type: "session.summary", session: "s1", ts: 0,
                 summary: { sessionId: "s1", title: "this window", workspace: "/w", state: "idle", pending: 0,
                            createdAt: 900, lastActivity: now } } as never)
  await w.settle()
  return w
}
const other = async (id: string, title: string, state: string, pending = 0) => {
  w!.store.push({ type: "session.summary", session: id, ts: 0,
                  summary: { sessionId: id, title, workspace: "/w", state, pending, createdAt: 901,
                             lastActivity: now } } as never)
  await w!.settle()
}
const hintsLine = () => w!.footer().filter((l) => l.trim()).at(-1) ?? ""
const screen = () => w!.t.captureCharFrame()

describe("the other sessions in the hints line", () => {
  test("none: the line says only where the agents are", async () => {
    await mount()
    expect(hintsLine()).toContain("← agents")
    expect(hintsLine()).not.toContain("more session")
  })

  test("counts what runs elsewhere in the service and what waits, live", async () => {
    await mount()
    await other("o1", "scan benzene", "running")
    expect(hintsLine()).toContain("· 1 more session ← agents")
    await other("o2", "relax Fe complex", "needs_input", 1)
    expect(hintsLine()).toContain("· 2 more sessions (1 waits for you) ← agents")
    await other("o3", "asleep one", "sleeping")
    await other("o4", "ended one", "ended")
    expect(hintsLine()).toContain("· 2 more sessions (1 waits for you) ← agents")
    // Answered elsewhere: it works again.
    await other("o2", "relax Fe complex", "running", 0)
    expect(hintsLine()).toContain("· 2 more sessions ← agents")
    w!.store.push({ type: "session.removed", session: "o1", ts: 0, reason: "shutdown" } as never)
    w!.store.push({ type: "session.removed", session: "o2", ts: 0, reason: "shutdown" } as never)
    await w!.settle()
    expect(hintsLine()).not.toContain("more session")
    expect(hintsLine()).toContain("← agents")
  })

  test("where the line is short it says less, never more than fits", () => {
    expect(othersHint({ running: 2, waiting: 1 }, "←", 200)).toBe("2 more sessions (1 waits for you) ← agents")
    expect(othersHint({ running: 3, waiting: 2 }, "←", 200)).toBe("3 more sessions (2 wait for you) ← agents")
    expect(othersHint({ running: 2, waiting: 1 }, "←", 30)).toBe("2 more (1 waiting) ← agents")
    expect(othersHint({ running: 2, waiting: 1 }, "←", 16)).toBe("2 more ← agents")
    expect(othersHint({ running: 0, waiting: 0 }, "←", 200)).toBe("← agents")
  })

  test("a click on it opens the agent view (full screen)", async () => {
    await mount({ full: true })
    await other("o1", "scan benzene", "running")
    const rows = screen().split("\n")
    const y = rows.findIndex((r) => r.includes("1 more session"))
    expect(y).toBeGreaterThanOrEqual(0)
    const x = rows[y]!.indexOf("1 more session")
    await w!.t.mockMouse.click(x + 2, y)
    await Bun.sleep(20)
    await w!.settle()
    expect(screen()).toContain("◆ AGENTS")
  })

  test("a click elsewhere on the line does not", async () => {
    await mount({ full: true })
    await other("o1", "scan benzene", "running")
    const rows = screen().split("\n")
    const y = rows.findIndex((r) => r.includes("1 more session"))
    await w!.t.mockMouse.click(3, y)
    await Bun.sleep(20)
    await w!.settle()
    expect(screen()).not.toContain("◆ AGENTS")
  })
})

describe("another session that starts to wait", () => {
  test("says so once in the message row, for a moment; a new wait says it again", async () => {
    await mount({ arrival: { shownMs: 200, quietMs: 0 } })
    const said = () => w!.footer().filter((l) => l.includes("waits for you —"))
    await other("o2", "relax Fe complex", "needs_input", 1)
    expect(said()).toHaveLength(1)
    expect(said()[0]).toContain("◉ session 'relax Fe complex' waits for you — ← to open")
    await Bun.sleep(300)
    await w!.settle()
    expect(said()).toHaveLength(0)
    // The same question, summarised again: not said again.
    await other("o2", "relax Fe complex", "needs_input", 1)
    expect(said()).toHaveLength(0)
    // Answered, then a new question: said again.
    await other("o2", "relax Fe complex", "running", 0)
    await other("o2", "relax Fe complex", "needs_input", 1)
    expect(said()).toHaveLength(1)
  })

  test("the session in view is not announced", async () => {
    await mount()
    await other("s1", "this window", "needs_input", 1)
    expect(w!.footer().join("\n")).not.toContain("waits for you —")
  })
})
