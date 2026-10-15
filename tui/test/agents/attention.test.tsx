/**
 * What needs the user from every DELFIN process, in the window: the
 * service's attention items become rows of NEEDS YOU (what waits, what
 * failed) and NOTICES (quiet, above ELSEWHERE); space peeks at one; →
 * opens it (the window decides how: its session, a copy, the diagnosis, a
 * new session with the item as its brief and a suggested prompt in the
 * box); typing + enter answers what waits; ctrl+x dismisses (what waits or
 * failed: twice). Seen from anywhere: a quiet count at the right end of
 * the hints line, a new item's note there for a moment, the count in the
 * terminal's title. A session started about an item shows it as one row.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { parseColor, rgbToHex, type CapturedSpan } from "@opentui/core"
import { attentionCounts, briefOf, itemPeek, openingOf, suggestedPrompt, whereOf } from "../../src/agents/attention"
import { agentRows, peekLines } from "../../src/agents/view"
import { attentionCorner } from "../../src/footer/footer"
import type { AttentionItem } from "../../src/protocol/generated"
import { Presenter } from "../../src/scrollback/present"
import { shownText } from "../../src/scrollback/rows"
import type { Item } from "../../src/scrollback/transcript"
import { applyEvent, attentionItems, emptyState, ensureSession } from "../../src/store/state"
import { DARK } from "../../src/theme/theme"
import { mountWindow, type Call } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
})

const HOME = process.env.HOME ?? "/home/u"
const now = Date.now() / 1000

function item(extra: Partial<AttentionItem> = {}): AttentionItem {
  return {
    id: "att-1", kind: "notice", group: "notices", blocking: false, name: "something", text: "it happened",
    detail: "", workspace: "", host: "", sessionId: "", sessionTitle: "", sessionOpen: false, sessionSaved: false,
    createdAt: now - 120, options: [], job: null, ...extra,
  }
}

const question = item({
  id: "att-q", kind: "question", group: "needs", blocking: true, name: "scan of benzene",
  text: "question: Which functional for the excited states?", detail: "Which functional for the excited states?",
  options: ["B3LYP", "CAM-B3LYP"], workspace: `${HOME}/calc/benzene`, host: "node2", sessionId: "dash-7",
  sessionTitle: "scan of benzene", sessionSaved: true, sessionOpen: true, createdAt: now - 7200,
})
const failed = item({
  id: "att-f", kind: "failed", group: "needs", name: "job 4711 failed",
  text: "OUT_OF_MEMORY · diagnosis ready · a fix is prepared", detail: "out of memory",
  workspace: `${HOME}/calc/benzene`, host: "node2", createdAt: now - 720,
  job: { id: "4711", state: "OUT_OF_MEMORY", folder: `${HOME}/calc/benzene`, signatures: ["oom-kill"],
         diagnosisSession: "diag-4711", summary: "the SCF ran out of memory",
         fix: { fix_class: "resource", recommendation: "raise maxcore",
                fix: { kind: "memory", summary_line: "maxcore 3000 -> 4500", diff: "-%maxcore 3000\n+%maxcore 4500" } } },
})
const finished = item({
  id: "att-n", kind: "finished", group: "notices", name: "job 4702 finished", text: "COMPLETED",
  workspace: `${HOME}/calc/water`, createdAt: now - 3600,
  job: { id: "4702", state: "COMPLETED", folder: `${HOME}/calc/water` },
})

describe("an attention item", () => {
  test("arrives as the service sends it: every field typed, one without an id left out", () => {
    const got = attentionItems([{ id: "a", kind: "question", group: "needs", blocking: 1, options: ["x", 2],
                                  createdAt: "12" }, { kind: "notice" }, null, "nonsense"])
    expect(got).toHaveLength(1)
    expect(got[0]).toMatchObject({ id: "a", kind: "question", group: "needs", blocking: true, options: ["x", "2"],
                                   createdAt: 12, name: "", sessionOpen: false })
    expect(attentionItems({ not: "a list" })).toEqual([])
    const state = emptyState()
    applyEvent(state, { type: "attention.updated", session: "", ts: 0, items: [question, failed] } as never)
    expect(state.attention.map((i) => i.id)).toEqual(["att-q", "att-f"])
    applyEvent(state, { type: "attention.updated", session: "", ts: 0, items: [] } as never)
    expect(state.attention).toEqual([])
  })

  test("is a row of NEEDS YOU or NOTICES: after this service's sessions, in the service's order", () => {
    const state = emptyState()
    const add = (id: string, extra: Record<string, unknown>) => {
      const s = ensureSession(state, id, {})
      s.summary = { sessionId: id, createdAt: 1, lastActivity: now, ...extra } as never
      s.state = String(extra.state ?? "idle")
    }
    add("here", {})
    add("asks", { state: "needs_input", pending: 1 })
    state.attention = [question, failed, finished]
    state.elsewhere = [{ key: "dash-1", title: "dashboard", workspace: "/w", host: "login2" }]
    state.host = "node1"
    const rows = agentRows(state, new Set())
    expect(rows.map((r) => [r.group, r.id])).toEqual([
      ["needs", "asks"], ["needs", "attention:att-q"], ["needs", "attention:att-f"],
      ["done", "here"], ["notices", "attention:att-n"], ["elsewhere", "dash-1"],
    ])
    const q = rows.find((r) => r.id === "attention:att-q")!
    expect(q).toMatchObject({ kind: "attention", status: "item", name: "scan of benzene", age: "2h",
                              where: "open on node2" })
    expect(q.item).toBe(question)
    expect(rows.find((r) => r.id === "attention:att-f")!.where).toBe("~/calc/benzene · node2")
  })

  test("says where it is: open elsewhere, or its folder and the machine when it is not this one", () => {
    expect(whereOf(question, "node2")).toBe("open elsewhere")
    expect(whereOf({ ...question, sessionOpen: false }, "node2")).toBe("~/calc/benzene")
    expect(whereOf({ ...question, sessionOpen: false }, "node1")).toBe("~/calc/benzene · node2")
  })

  test("is opened as its session, a copy of it, the job's diagnosis, or a new session", () => {
    expect(openingOf(question)).toEqual({ how: "session", id: "dash-7", copy: true, title: "scan of benzene" })
    expect(openingOf({ ...question, sessionOpen: false })).toMatchObject({ how: "session", copy: false })
    expect(openingOf(failed)).toEqual({ how: "diagnosis", id: "diag-4711", title: "the diagnosis of job 4711" })
    expect(openingOf(finished)).toEqual({ how: "new" })
    // A session the store does not have (never saved) cannot be resumed.
    expect(openingOf({ ...question, sessionSaved: false })).toEqual({ how: "new" })
  })

  test("a peek at a question: the question, its choices, who asked, where an answer goes, what → does", () => {
    const lines = itemPeek(question, "node1", now * 1000)
    expect(lines).toEqual([
      "Which functional for the excited states?",
      "  1. B3LYP",
      "  2. CAM-B3LYP",
      "asked by scan of benzene · in ~/calc/benzene on node2 · 2h ago",
      "type the answer and press enter: it goes to the next session started in ~/calc/benzene",
      "→ opens a copy of scan of benzene (it is still open on node2)",
    ])
  })

  test("a peek at a failed job: its state, errors, the diagnosis and the prepared fix", () => {
    const lines = itemPeek(failed, "node2", now * 1000)
    expect(lines).toEqual([
      "job 4711 · OUT_OF_MEMORY · ~/calc/benzene",
      "errors: oom-kill",
      "the SCF ran out of memory",
      "prepared fix: maxcore 3000 -> 4500 (not applied: the dashboard's /fix applies it)",
      "on node2 · 12m ago",
      "→ opens the job monitor's diagnosis of it",
    ])
    expect(itemPeek(finished, "node1", now * 1000).at(-1)).toBe("→ starts a session about it in ~/calc/water, which waits for you")
    // A notice: what it says (its row has its name), from where and when.
    const budget = item({ kind: "budget", name: "budget nearly spent", text: "long refactor · 80% of $2.00",
                          detail: "long refactor · 80% of $2.00", host: "node2", createdAt: now - 10800 })
    expect(itemPeek(budget, "node1", now * 1000)).toEqual([
      "long refactor · 80% of $2.00", "on node2 · 3h ago", "→ starts a session about it, which waits for you"])
  })

  test("a new session about it reads all of it first, and suggests a prompt", () => {
    const brief = briefOf(failed, now * 1000)
    for (const part of ["Job 4711 failed (OUT_OF_MEMORY).", `Folder: ${HOME}/calc/benzene`, "Machine: node2",
                        "Errors found in its output: oom-kill", "The job monitor's diagnosis: the SCF ran out of memory",
                        "Recommendation: raise maxcore", "A fix is prepared (not applied): maxcore 3000 -> 4500",
                        "+%maxcore 4500", "Details:\nout of memory"]) {
      expect(brief).toContain(part)
    }
    expect(suggestedPrompt(failed)).toBe("Job 4711 failed with OUT_OF_MEMORY. What went wrong, and what should we change?")
    expect(suggestedPrompt(finished)).toBe("Job 4702 finished. Summarise its results.")
    const asked = briefOf(question)
    expect(asked).toContain("The question: Which functional for the excited states?")
    expect(asked).toContain("Its choices: B3LYP | CAM-B3LYP")
    expect(asked).toContain("From the session: scan of benzene (dash-7)")
  })

  test("counts at the bottom right: what waits (here and elsewhere) and what failed, no notices", () => {
    const state = emptyState()
    state.attention = [question, failed, finished]
    expect(attentionCounts(state, 2)).toEqual({ asks: 3, failed: 1 })
    expect(attentionCorner({ asks: 3, failed: 1, note: "" }, 40)).toMatchObject({ mark: "◉", text: " 3", sep: " · ",
                                                                                   mark2: "✗", text2: " 1" })
    expect(attentionCorner({ asks: 0, failed: 2, note: "" }, 40)).toMatchObject({ mark: "✗", text: " 2", sep: "" })
    expect(attentionCorner({ asks: 0, failed: 0, note: "" }, 40).mark).toBe("")
    expect(attentionCorner({ asks: 1, failed: 0, note: "◉ relax Fe complex asks" }, 40))
      .toMatchObject({ mark: "◉", text: " relax Fe complex asks", note: true })
  })
})

// -- in the window ------------------------------------------------------------

async function mount(options: { client?: (call: Call) => unknown; onOpenItem?: (item: AttentionItem) => Promise<string | null>;
                                quietMs?: number; shownMs?: number } = {}) {
  const calls: Call[] = []
  const client = {
    request: async (method: string, params: Record<string, any> = {}) => {
      calls.push({ method, params })
      const own = options.client?.({ method, params })
      if (own !== undefined) return own
      if (method === "sessions.list") return { sessions: [], elsewhere: [] }
      if (method === "session.replay") return { events: [] }
      if (method === "attention.answer" || method === "attention.dismiss") return { ok: true }
      return { accepted: true, applied: true, message: "" }
    },
  }
  w = await mountWindow({ width: 110, height: 30, client, onOpenItem: options.onOpenItem,
                          arrival: { shownMs: options.shownMs ?? 5000, quietMs: options.quietMs ?? 600_000 } })
  w.store.push({ type: "session.summary", session: "s1", ts: 0,
                 summary: { sessionId: "s1", title: "this window", workspace: "/w", state: "idle", pending: 0,
                            createdAt: 900, lastActivity: now } } as never)
  await w.settle()
  return { w, calls }
}
const items = async (list: AttentionItem[]) => {
  w!.store.push({ type: "attention.updated", session: "", ts: 0, items: list } as never)
  await w!.settle()
}
const panel = () => w!.footer().join("\n")
/** The hints line: the footer's last line with something on it. */
const hintsLine = () => w!.footer().filter((l) => l.trim()).at(-1) ?? ""
const called = (calls: Call[], method: string) => calls.filter((c) => c.method === method).map((c) => c.params)
const left = async () => {
  w!.keys.pressArrow("left")
  await w!.settle()
}
/** Select the row whose line holds `text` (from the top of the list). */
const select = async (text: string) => {
  for (let i = 0; i < 12; i++) {
    if (panel().split("\n").some((l) => l.startsWith(" › ") && l.includes(text))) return
    w!.keys.pressArrow("down")
    await w!.settle()
  }
  throw new Error(`no row with ${text}:\n${panel()}`)
}
const hex = (c: Parameters<typeof rgbToHex>[0]) => rgbToHex(c).toLowerCase()
const spans = (): CapturedSpan[] => w!.t.captureSpans().lines.flatMap((l) => l.spans)

describe("the agent view's NEEDS YOU and NOTICES", () => {
  test("list every item with its mark, what it needs, where and how old; notices are quiet", async () => {
    await mount()
    await items([question, failed, finished])
    await left()
    const text = panel()
    expect(text).toContain("◆ AGENTS  1 session · 2 need you")
    expect(text).toMatch(/NEEDS YOU\s*\n.*◉ scan of benzene\s+question: Which functional for the excited states\?\s+open on node2\s+2h/)
    expect(text).toMatch(/✗ job 4711 failed\s+OUT_OF_MEMORY · diagnosis ready · a fix is prepared\s+~\/calc\/benzene · node2\s+12m/)
    expect(text).toMatch(/NOTICES\s*\n.*· job 4702 finished\s+COMPLETED\s+~\/calc\/water\s+1h/)
    expect(text.indexOf("NEEDS YOU")).toBeLessThan(text.indexOf("DONE"))
    expect(text.indexOf("DONE")).toBeLessThan(text.indexOf("NOTICES"))
    // The marks in their colours; a notice's words faint.
    const all = spans()
    expect(hex(all.find((s) => s.text.startsWith("◉ "))!.fg)).toBe(hex(parseColor(DARK.warn)))
    expect(hex(all.find((s) => s.text.startsWith("✗ "))!.fg)).toBe(hex(parseColor(DARK.err)))
    expect(hex(all.find((s) => s.text.includes("COMPLETED"))!.fg)).toBe(hex(parseColor(DARK.faint)))
  })

  test("space peeks at an item without asking the service; the hints say what its keys do", async () => {
    const { calls } = await mount()
    await items([failed])
    await left()
    await select("job 4711 failed")
    expect(panel()).toContain("↑↓ select · → open · space peek · ctrl+x dismiss · esc back")
    w!.keys.pressKey(" ")
    await w!.settle()
    const text = panel()
    expect(text).toContain("│ job 4711 · OUT_OF_MEMORY · ~/calc/benzene")
    expect(text).toContain("│ prepared fix: maxcore 3000 -> 4500")
    expect(text).toContain("│ → opens the job monitor's diagnosis of it")
    expect(called(calls, "session.replay")).toEqual([])
    expect(text).toContain("→ open · ctrl+x dismiss · space closes · esc back")
    expect(text).toContain("❯ → opens it in a session")
  })

  test("→ on an item hands it to the window to open, and its suggested prompt waits in the box", async () => {
    const opened: AttentionItem[] = []
    await mount({ onOpenItem: async (it) => { opened.push(it); return suggestedPrompt(it) } })
    await items([failed, finished])
    await left()
    await select("job 4702 finished")
    w!.keys.pressArrow("right")
    await w!.settle()
    expect(opened.map((i) => i.id)).toEqual(["att-n"])
    const text = panel()
    expect(text).not.toContain("◆ AGENTS")                       // the view closed
    expect(text).toContain("❯ Job 4702 finished. Summarise its results.")
  })

  test("typing + enter on a peeked question answers it; the row goes and a note says where the answer goes", async () => {
    const { calls } = await mount()
    await items([question, finished])
    await left()
    await select("scan of benzene")
    w!.keys.pressKey(" ")
    await w!.settle()
    expect(panel()).toContain("type the answer · enter sends it")
    expect(panel()).toContain("❯ the answer for scan of benzene…")
    await w!.keys.typeText("CAM-B3LYP, with D3")
    w!.keys.pressEnter()
    await w!.settle()
    expect(called(calls, "attention.answer")).toEqual([{ id: "att-q", text: "CAM-B3LYP, with D3" }])
    const text = panel()
    expect(text).toContain("answered: the next session started in ~/calc/benzene gets it")
    expect(text).not.toContain("scan of benzene  question")
    expect(called(calls, "session.prompt")).toEqual([])            // never a session's prompt
  })

  test("an answer the inbox did not take puts the row and the text back", async () => {
    await mount({ client: (c) => c.method === "attention.answer" ? { ok: false, error: "item att-q is no longer pending" } : undefined })
    await items([question])
    await left()
    await select("scan of benzene")
    w!.keys.pressKey(" ")
    await w!.settle()
    await w!.keys.typeText("B3LYP")
    w!.keys.pressEnter()
    await w!.settle()
    await w!.settle()
    const text = panel()
    expect(text).toContain("not answered: item att-q is no longer pending")
    expect(text).toContain("◉ scan of benzene")
    expect(text).toContain("❯ B3LYP")
  })

  test("typing + enter on an item that waits on nobody sends nothing and keeps the text", async () => {
    const { calls } = await mount()
    await items([finished])
    await left()
    await select("job 4702 finished")
    w!.keys.pressKey(" ")
    await w!.settle()
    await w!.keys.typeText("thanks")
    w!.keys.pressEnter()
    await w!.settle()
    expect(called(calls, "attention.answer")).toEqual([])
    expect(called(calls, "session.create")).toEqual([])
    expect(panel()).toContain("this item takes no answer: → opens it in a session")
    expect(panel()).toContain("❯ thanks")
  })

  test("ctrl+x dismisses a notice at once, and what waits only when pressed twice", async () => {
    const { calls } = await mount()
    await items([question, finished])
    await left()
    await select("job 4702 finished")
    w!.keys.pressKey("x", { ctrl: true })
    await w!.settle()
    expect(called(calls, "attention.dismiss")).toEqual([{ id: "att-n" }])
    expect(panel()).not.toContain("job 4702 finished  COMPLETED")
    expect(panel()).toContain("dismissed: job 4702 finished")
    await select("scan of benzene")
    w!.keys.pressKey("x", { ctrl: true })
    await w!.settle()
    expect(called(calls, "attention.dismiss")).toHaveLength(1)      // not yet
    expect(panel()).toContain("ctrl+x again dismisses it: nobody waiting on it gets an answer")
    w!.keys.pressKey("x", { ctrl: true })
    await w!.settle()
    expect(called(calls, "attention.dismiss")).toEqual([{ id: "att-n" }, { id: "att-q" }])
    expect(panel()).not.toContain("◉ scan of benzene")
  })

  test("a dismiss that did not go through brings the row back", async () => {
    await mount({ client: (c) => c.method === "attention.dismiss" ? { ok: false, error: "item att-n is no longer pending" } : undefined })
    await items([finished])
    await left()
    await select("job 4702 finished")
    w!.keys.pressKey("x", { ctrl: true })
    await w!.settle()
    await w!.settle()
    expect(panel()).toContain("not dismissed: item att-n is no longer pending")
    expect(panel()).toContain("· job 4702 finished")
  })

  test("an item is not renamed or pinned", async () => {
    const { calls } = await mount()
    await items([finished])
    await left()
    await select("job 4702 finished")
    w!.keys.pressKey("r", { ctrl: true })
    await w!.settle()
    expect(panel()).toContain("this is no session: → opens it, ctrl+x dismisses it")
    w!.keys.pressKey("t", { ctrl: true })
    await w!.settle()
    expect(called(calls, "session.setTitle")).toEqual([])
    expect(panel()).not.toContain("PINNED")
  })
})

describe("what needs the user, seen from anywhere in the window", () => {
  test("a quiet count at the right end of the hints line, the marks in colour; the title counts too", async () => {
    await mount()
    const titles: string[] = []
    const renderer = w!.t.renderer as unknown as { setTerminalTitle(t: string): void }
    const original = renderer.setTerminalTitle.bind(renderer)
    renderer.setTerminalTitle = (t: string) => { titles.push(t); original(t) }
    await items([question, failed, finished])
    const hints = hintsLine()!
    expect(hints).toMatch(/← agents +◉ 1 · ✗ 1 *$/)
    const all = spans()
    const mark = all.find((s) => s.text === "◉")!
    expect(hex(mark.fg)).toBe(hex(parseColor(DARK.warn)))
    expect(hex(all.find((s) => s.text === "✗")!.fg)).toBe(hex(parseColor(DARK.err)))
    expect(hex(all.find((s) => s.text === " 1")!.fg)).toBe(hex(parseColor(DARK.faint)))
    expect(titles.at(-1)).toBe("(2) DELFIN · this window")
    await items([])
    expect(hintsLine()).not.toMatch(/[◉✗] \d/)
    expect(titles.at(-1)).toBe("DELFIN · this window")
  })

  test("a new item that waits or failed names itself there for a moment, then the count is back", async () => {
    await mount({ quietMs: 0, shownMs: 300 })
    await items([{ ...question, createdAt: now }])
    expect(hintsLine()).toMatch(/◉ scan of benzene asks *$/)
    await Bun.sleep(400)
    await w!.settle()
    expect(hintsLine()).toMatch(/← agents +◉ 1 *$/)
    // A notice is quiet; an old item the window had not heard of is no news.
    await items([{ ...question, createdAt: now }, { ...finished, createdAt: now }, { ...failed, createdAt: now - 3600 }])
    expect(hintsLine()).toMatch(/← agents +◉ 1 · ✗ 1 *$/)
    await items([{ ...question, createdAt: now }, { ...failed, id: "att-f2", name: "job 4712 failed", createdAt: now }])
    expect(hintsLine()).toMatch(/✗ job 4712 failed *$/)
  })

  test("what the window finds when it opens is no news", async () => {
    await mount({ quietMs: 600_000 })
    await items([{ ...question, createdAt: now }])
    expect(hintsLine()).toMatch(/← agents +◉ 1 *$/)
  })

  test("another session here that starts to wait names itself in the message row; the corner counts it", async () => {
    await mount({ quietMs: 0 })
    w!.store.push({ type: "session.summary", session: "q9", ts: 0,
                    summary: { sessionId: "q9", title: "relax Fe complex", workspace: "/w", state: "needs_input",
                               pending: 1, createdAt: 901, lastActivity: now } } as never)
    await w!.settle()
    expect(panel()).toContain("◉ session 'relax Fe complex' waits for you — ← to open")
    expect(hintsLine()).toMatch(/1 more session \(1 waits for you\) ← agents +◉ 1 *$/)
  })
})

describe("a session started about an item", () => {
  test("shows its brief as one row, like a tool call's; ctrl+o shows all of it", () => {
    const added: Item[] = []
    const transcript = { append: async (i: Item) => { added.push(i) } } as never
    const presenter = new Presenter("s1", transcript, { showThinking: false })
    const state = emptyState()
    ensureSession(state, "s1", {})
    presenter.handle([{ type: "context", session: "s1", ts: 0, seq: 1, title: "job 4711 failed",
                        text: "Job 4711 failed (OUT_OF_MEMORY).\nFolder: /calc/benzene\n" } as never], state)
    expect(added).toHaveLength(1)
    const row = added[0] as Record<string, any>
    expect(row.kind).toBe("tool")
    expect(row.name).toBe("context")
    expect(row.arg).toBe("job 4711 failed")
    expect(row.result).toMatchObject({ kind: "text", ok: true, wrap: true,
                                       summary: "2 lines, given to the model with the first prompt",
                                       lines: ["Job 4711 failed (OUT_OF_MEMORY).", "Folder: /calc/benzene"] })
    expect(shownText(added[0]!)).toContain("Folder: /calc/benzene")
    expect(peekLines([{ type: "context", title: "job 4711 failed", text: "…" }])).toEqual(["⏺ context  job 4711 failed"])
  })

  test("its row in the window: the headline and the summary; the details wrap, not cut", async () => {
    w = await mountWindow({ width: 60, height: 30, full: true })
    const long = "The job monitor's diagnosis: the SCF ran out of memory in the second cycle of the excited-state run"
    w.store.push({ type: "context", session: "s1", ts: 0, seq: 1, title: "job 4711 failed",
                   text: `Job 4711 failed (OUT_OF_MEMORY).\n${long}` } as never)
    await w.settle()
    let screen = w.screenLines().join("\n")
    expect(screen).toContain("⏺ context  job 4711 failed")
    expect(screen).toContain("⎿ 2 lines, given to the model with the first prompt")
    expect(screen).not.toContain("second cycle")
    w.keys.pressKey("o", { ctrl: true })
    await w.settle()
    screen = w.screenLines().join("\n")
    expect(screen).toContain("Job 4711 failed (OUT_OF_MEMORY).")
    expect(screen).toContain("excited-state run")              // wrapped onto the next row, all of it
  })
})

describe("what the review found", () => {
  test("an answer the inbox did not take is tried again by enter -- never a new session with it", async () => {
    let refused = 1
    const news: string[] = []
    const { calls } = await mount({ client: (c) => c.method === "attention.answer" && refused-- > 0
      ? { ok: false, error: "the inbox is locked" } : undefined })
    w!.store.update((d) => { d.starting = false })
    await items([question])
    await left()
    await select("scan of benzene")
    w!.keys.pressKey(" ")
    await w!.settle()
    await w!.keys.typeText("B3LYP")
    w!.keys.pressEnter()
    await w!.settle()
    await w!.settle()
    expect(panel()).toContain("type the answer · enter sends it")    // its peek is back
    w!.keys.pressEnter()                                             // try again
    await w!.settle()
    await w!.settle()
    expect(called(calls, "attention.answer")).toEqual([{ id: "att-q", text: "B3LYP" }, { id: "att-q", text: "B3LYP" }])
    expect(called(calls, "session.create")).toEqual([])
    expect(news).toEqual([])
    expect(panel()).toContain("answered: the next session started in ~/calc/benzene gets it")
  })

  test("ctrl+x pressed twice on a notice dismisses it, not the one that moves up under it", async () => {
    const { calls } = await mount()
    const n2 = { ...finished, id: "att-n2", name: "job 4703 finished", createdAt: finished.createdAt - 10 }
    await items([finished, n2])
    await left()
    await select("job 4702 finished")
    w!.keys.pressKey("x", { ctrl: true })
    await w!.settle()
    w!.keys.pressKey("x", { ctrl: true })                            // by habit
    await w!.settle()
    expect(called(calls, "attention.dismiss")).toEqual([{ id: "att-n" }])
    expect(panel()).toContain("ctrl+x again dismisses it")
    expect(panel()).toContain("· job 4703 finished")
    w!.keys.pressKey("x", { ctrl: true })                            // meant
    await w!.settle()
    expect(called(calls, "attention.dismiss")).toEqual([{ id: "att-n" }, { id: "att-n2" }])
    // Moved on purpose, a notice goes at once again.
    const n3 = { ...finished, id: "att-n3", name: "job 4704 finished" }
    const n4 = { ...finished, id: "att-n4", name: "job 4705 finished", createdAt: finished.createdAt - 20 }
    await items([n3, n4])
    await select("job 4705 finished")
    w!.keys.pressKey("x", { ctrl: true })
    await w!.settle()
    expect(called(calls, "attention.dismiss").at(-1)).toEqual({ id: "att-n4" })
  })

  test("/rename on an item renames nothing", async () => {
    const { calls } = await mount()
    await items([question])
    await left()
    await select("scan of benzene")
    await w!.keys.typeText("/rename foo")
    w!.keys.pressEnter()
    await w!.settle()
    expect(called(calls, "session.setTitle")).toEqual([])
    expect(panel()).toContain("this is no session: → opens it, ctrl+x dismisses it")
  })

  test("an item answered elsewhere while peeked: its peek goes, or says it is gone while something is typed", async () => {
    await mount()
    await items([question, failed])
    await left()
    await select("scan of benzene")
    w!.keys.pressKey(" ")
    await w!.settle()
    await items([failed])                                            // answered in the dashboard
    expect(panel()).not.toContain("│ Which functional")
    expect(panel()).toContain("a new session: what should it do?")
    // With an answer half typed: kept, and the hints say what happened.
    await items([question, failed])
    await select("scan of benzene")
    w!.keys.pressKey(" ")
    await w!.settle()
    await w!.keys.typeText("CAM")
    await items([failed])
    expect(panel()).toContain("that item is gone: answered or dismissed elsewhere")
    w!.keys.pressEnter()
    await w!.settle()
    expect(panel()).toContain("❯ CAM")                                // not sent anywhere
    w!.keys.pressArrow("down")                                       // moving on closes it
    await w!.settle()
    expect(panel()).not.toContain("that item is gone")
  })

  test("of two items that arrive together, the one that waits names itself", async () => {
    await mount({ quietMs: 0 })
    await items([{ ...question, createdAt: now }, { ...failed, createdAt: now }])
    expect(hintsLine()).toMatch(/◉ scan of benzene asks \(\+1 more\) *$/)
  })

  test("an update older than the list shown is not applied", () => {
    const state = emptyState()
    applyEvent(state, { type: "attention.updated", session: "", ts: 0, items: [question, failed], version: 7 } as never)
    applyEvent(state, { type: "attention.updated", session: "", ts: 0, items: [question], version: 6 } as never)
    expect(state.attention.map((i) => i.id)).toEqual(["att-q", "att-f"])
    applyEvent(state, { type: "attention.updated", session: "", ts: 0, items: [failed], version: 8 } as never)
    expect(state.attention.map((i) => i.id)).toEqual(["att-f"])
    expect(state.attentionVersion).toBe(8)
  })

  test("a blank line of a context keeps its row", async () => {
    w = await mountWindow({ width: 60, height: 30, full: true })
    w.store.push({ type: "context", session: "s1", ts: 0, seq: 1, title: "job 4711 failed",
                   text: "Job 4711 failed.\n\nDetails:\nout of memory" } as never)
    await w.settle()
    w.keys.pressKey("o", { ctrl: true })
    await w.settle()
    const lines = w.screenLines().map((l) => l.trimEnd())
    const at = lines.findIndex((l) => l.includes("Job 4711 failed."))
    expect(lines[at + 1]).toBe("")
    expect(lines[at + 2]).toContain("Details:")
  })
})
