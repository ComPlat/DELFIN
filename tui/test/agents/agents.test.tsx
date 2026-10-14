/**
 * The agent view in the window, driven with keys: ← opens it on an empty
 * prompt; its rows are grouped by what each session needs and keep their
 * place; → opens a session, space peeks at it (a reply goes to it, not to
 * the session in view), typing + enter starts a new one, ctrl+x stops,
 * ctrl+r renames, ctrl+t pins. Letters are always typing. Rendered
 * headless through the real App against a fake agent service.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"
import { Pins } from "../../src/agents/pins"
import { agentRows, peekLines } from "../../src/agents/view"
import { Keymap } from "../../src/keymap/keymap"
import { emptyState, ensureSession } from "../../src/store/state"
import { mountWindow, type Call } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
const dirs: string[] = []
afterEach(async () => {
  await w?.close()
  w = null
  for (const d of dirs.splice(0)) rmSync(d, { recursive: true, force: true })
})

const summary = (id: string, extra: Record<string, unknown> = {}) => ({
  sessionId: id, title: null, workspace: `/work/${id}`, provider: "kit", model: "kit.test", approval: "default",
  state: "idle", pending: 0, costUsd: 0, createdAt: 1000, lastActivity: Date.now() / 1000, pid: 1,
  lastPrompt: null, activity: null, waiting: null, answer: null, ...extra,
})

/** The window with three more sessions in the service, as main.tsx
 *  registers them from the hello. */
async function mount(options: { client?: (call: Call) => unknown; pins?: Pins; onNewSession?: () => Promise<string | null>;
                                height?: number; elsewhere?: Record<string, unknown>[] } = {}) {
  const opened: [string, string][] = []
  const calls: Call[] = []
  // What the fake service runs: what the test announced (a real one lists
  // every session it has; one not listed has ended).
  const running = new Map<string, Record<string, unknown>>()
  const client = {
    request: async (method: string, params: Record<string, any> = {}) => {
      const call = { method, params }
      calls.push(call)
      const own = options.client?.(call)
      if (own !== undefined) return own
      if (method === "sessions.list") return { sessions: [...running.values()], elsewhere: options.elsewhere ?? [] }
      if (method === "session.replay") return { events: [] }
      if (method === "session.setTitle") return { title: params.title }
      return { accepted: true, applied: true, message: "" }
    },
  }
  w = await mountWindow({ width: 100, height: options.height ?? 30, client, pins: options.pins,
                          onResume: (id, ws) => { opened.push([id, ws]) }, onNewSession: options.onNewSession })
  const push = (s: Record<string, unknown>) => {
    running.set(String(s.sessionId), s)
    w!.store.push({ type: "session.summary", session: s.sessionId, ts: 0, summary: s } as never)
  }
  push(summary("s1", { title: "this window", createdAt: 900 }))
  push(summary("w1", { title: "optimise benzene", state: "running", activity: "bash  orca opt.inp", createdAt: 1001 }))
  push(summary("q1", { title: "plan the scan", state: "needs_input", pending: 1, waiting: "bash: sbatch scan.sh", createdAt: 1002 }))
  push(summary("d1", { title: "write the report", answer: "The report is in report.md.", createdAt: 1003 }))
  await w.settle()
  return { w: w!, calls, opened, push }
}

const left = async () => {
  w!.keys.pressArrow("left")
  await w!.settle()
}
const panel = () => w!.footer().join("\n")
const called = (calls: Call[], method: string) => calls.filter((c) => c.method === method).map((c) => c.params)

describe("the agent view's rows", () => {
  test("are grouped by what each session needs, newest session first, pins on top", () => {
    const state = emptyState()
    const add = (id: string, extra: Record<string, unknown>) => {
      const s = ensureSession(state, id, {})
      s.summary = summary(id, extra) as never
      s.state = String(extra.state ?? "idle")
    }
    add("old", { createdAt: 1 })
    add("new", { createdAt: 5 })
    add("busy", { state: "running", createdAt: 2 })
    add("ask", { state: "needs_input", pending: 1, createdAt: 3 })
    add("pin", { createdAt: 0 })
    state.sessions.gone = { ...ensureSession(state, "gone", {}), state: "ended", ended: { reason: "shutdown" } }
    state.elsewhere = [{ key: "dash-1", title: "dashboard", workspace: "/w", host: "login2" }]
    state.current = "old"
    const rows = agentRows(state, new Set(["pin"]))
    expect(rows.map((r) => [r.group, r.id])).toEqual([
      ["pinned", "pin"], ["needs", "ask"], ["working", "busy"],
      ["done", "new"], ["done", "old"], ["done", "gone"], ["elsewhere", "dash-1"],
    ])
    expect(rows.find((r) => r.id === "old")!.current).toBe(true)
    expect(rows.find((r) => r.id === "gone")!.activity).toBe("ended (shutdown)")
    expect(rows.find((r) => r.id === "ask")!.activity).toContain("waiting for you")
  })

  test("a peek shows an answer as plain text", () => {
    expect(peekLines([{ type: "text.ended", text: "## Summary\n\nThe workspace holds **one** `README`.\n- first\n* second" }]))
      .toEqual(["Summary", "The workspace holds one README.", "• first", "• second"])
  })

  test("a peek reads a turn as a few lines", () => {
    expect(peekLines([
      { type: "user.message", text: "run it\nplease" },
      { type: "tool.started", callId: "c1", name: "bash", input: { command: "sbatch job.sh" } },
      { type: "tool.ended", callId: "c1", name: "bash" },
      { type: "text.ended", text: "Submitted.\n\nJob 4711 is queued." },
      { type: "permission.asked", requestId: "r1", kind: "confirm", tool: "bash", preview: "$ scancel 4711" },
      { type: "turn.ended", durationMs: 2300 },
    ])).toEqual(["❯ run it", "⏺ bash  sbatch job.sh", "Submitted.", "Job 4711 is queued.", "✓ 2.3s",
                 "◉ waiting for you: bash $ scancel 4711"])
  })
})

describe("the agent view", () => {
  test("← on an empty prompt opens it over the whole terminal, esc closes it", async () => {
    const { w } = await mount()
    const before = w.t.renderer.footerHeight
    await left()
    expect(w.t.renderer.footerHeight).toBe(29)
    const text = panel()
    expect(text).toContain("◆ AGENTS")
    expect(text).toContain("4 sessions · 1 needs input · 1 working")
    for (const title of ["NEEDS INPUT", "WORKING", "DONE", "plan the scan", "optimise benzene", "write the report"]) {
      expect(text).toContain(title)
    }
    expect(text.indexOf("NEEDS INPUT")).toBeLessThan(text.indexOf("WORKING"))
    expect(text).toContain("waiting for you: bash: sbatch scan.sh")
    expect(text).toContain("type to start a new session")              // the hints line is the view's
    await w.escape()
    expect(w.t.renderer.footerHeight).toBe(before)
    expect(panel()).not.toContain("◆ AGENTS")
  })

  test("← with something typed moves the cursor, it opens nothing", async () => {
    const { w } = await mount()
    await w.keys.typeText("abc")
    await left()
    expect(panel()).not.toContain("◆ AGENTS")
  })

  test("the selection starts on the session in view and follows its session", async () => {
    const { w, opened, push } = await mount()
    await left()
    expect(panel()).toMatch(/›\s+\S\s+this window/)
    w.keys.pressArrow("up")
    await w.settle()
    expect(panel()).toMatch(/›\s+\S\s+write the report/)
    // It starts to work: its row moves to WORKING, the selection with it.
    push(summary("d1", { title: "write the report", state: "running", activity: "writing…", createdAt: 1003 }))
    await w.settle()
    expect(panel()).toMatch(/›\s+\S\s+write the report\s+writing…/)
    w.keys.pressArrow("right")
    await w.settle()
    expect(opened).toEqual([["d1", "/work/d1"]])
    expect(panel()).not.toContain("◆ AGENTS")
  })

  test("enter on the session in view only closes the view", async () => {
    const { w, opened } = await mount()
    await left()
    w.keys.pressEnter()
    await w.settle()
    expect(opened).toEqual([])
    expect(panel()).not.toContain("◆ AGENTS")
  })

  test("letters are typing: j and k go into the prompt, typing + enter starts a new session", async () => {
    let started = 0
    const { w, calls } = await mount({
      onNewSession: async () => {
        started++
        w!.store.update((d) => { d.starting = true })
        await Bun.sleep(20)
        w!.store.update((d) => {
          ensureSession(d, "n1", { title: "" })
          d.current = "n1"
          d.starting = false
        })
        return "n1"
      },
    })
    await left()
    await w.keys.typeText("jk check the scan")
    await w.settle()
    expect(panel()).toContain("❯ jk check the scan")
    w.keys.pressEnter()
    await Bun.sleep(60)
    await w.settle()
    expect(started).toBe(1)
    expect(panel()).not.toContain("◆ AGENTS")
    expect(called(calls, "session.prompt")).toEqual([{ sessionId: "n1", text: "jk check the scan" }])
  })

  test("a new session that comes up while the view is open again gets its prompt once", async () => {
    let started = 0
    let up!: () => void
    const { w, calls } = await mount({
      onNewSession: async () => {
        started++
        w!.store.update((d) => { d.starting = true })
        await new Promise<void>((resolve) => { up = resolve })
        w!.store.update((d) => {
          ensureSession(d, "n1", {})
          d.current = "n1"
          d.starting = false
        })
        return "n1"
      },
    })
    await left()
    await w.keys.typeText("first job")
    w.keys.pressEnter()
    await w.settle()
    await left()                                        // the view again, while it starts
    expect(panel()).toContain("◆ AGENTS")
    up()
    await Bun.sleep(20)
    await w.settle()
    expect(started).toBe(1)
    expect(called(calls, "session.prompt")).toEqual([{ sessionId: "n1", text: "first job" }])
  })

  test("a new session that does not start gives the prompt back", async () => {
    const { w, calls } = await mount({ onNewSession: async () => null })
    await left()
    await w.keys.typeText("start the opt")
    w.keys.pressEnter()
    await Bun.sleep(20)
    await w.settle()
    expect(called(calls, "session.prompt")).toEqual([])
    expect(panel()).toContain("❯ start the opt")
  })

  test("space peeks; a reply goes to the peeked session, the view stays", async () => {
    const { w, calls } = await mount({
      client: ({ method }) => method === "session.replay" ? {
        events: [{ type: "user.message", text: "optimise it" },
                 { type: "text.ended", text: "Running the optimisation now." }],
      } : undefined,
    })
    await left()
    w.keys.pressArrow("down")                           // from the last row, round to the first
    w.keys.pressArrow("down")                           // working: optimise benzene
    await w.settle()
    w.keys.pressKey(" ")
    await Bun.sleep(10)
    await w.settle()
    expect(called(calls, "session.replay")).toEqual([{ sessionId: "w1", turns: 1 }])
    expect(panel()).toContain("│ ❯ optimise it")
    expect(panel()).toContain("│ Running the optimisation now.")
    expect(panel()).toContain("reply to optimise benzene…")
    await w.keys.typeText("use def2-TZVP")
    w.keys.pressEnter()
    await w.settle()
    expect(called(calls, "session.prompt")).toEqual([{ sessionId: "w1", text: "use def2-TZVP" }])
    expect(panel()).toContain("◆ AGENTS")
  })

  test("/rename as a reply names the peeked session, and sends it no prompt", async () => {
    const { w, calls } = await mount()
    await left()
    w.keys.pressArrow("down")
    w.keys.pressArrow("down")                           // working: optimise benzene
    await w.settle()
    w.keys.pressKey(" ")
    await Bun.sleep(10)
    await w.settle()
    await w.keys.typeText("/rename benzene opt")
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(called(calls, "session.setTitle")).toEqual([{ sessionId: "w1", title: "benzene opt" }])
    expect(called(calls, "session.prompt")).toEqual([])
  })

  test("/rename in the view with nothing peeked names the selected session, and starts none", async () => {
    let started = 0
    const { w, calls } = await mount({ onNewSession: async () => { started++; return null } })
    await left()
    w.keys.pressArrow("up")                             // write the report
    await w.settle()
    await w.keys.typeText("/rename the report")
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(called(calls, "session.setTitle")).toEqual([{ sessionId: "d1", title: "the report" }])
    expect(called(calls, "session.prompt")).toEqual([])
    expect(started).toBe(0)
  })

  test("/rename to a session elsewhere is not mailed to it", async () => {
    const { w, calls } = await mount({ elsewhere: [{ key: "dash1", title: "dashboard scan", workspace: "/w", host: "here" }] })
    await Bun.sleep(50)                                 // the view's poll of the list is in
    await w.settle()
    await left()
    for (let i = 0; i < 8 && !/›\s+\S\s+dashboard scan/.test(panel()); i++) {
      w.keys.pressArrow("down")
      await w.settle()
    }
    w.keys.pressKey(" ")
    await Bun.sleep(10)
    await w.settle()
    await w.keys.typeText("/rename x")
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(called(calls, "session.message")).toEqual([])
    expect(panel()).toContain("renamed where it runs")
  })

  test("ctrl+x twice stops a session; an ended one leaves the list", async () => {
    const { w, calls } = await mount()
    await left()
    w.keys.pressArrow("up")                             // write the report
    await w.settle()
    w.keys.pressKey("x", { ctrl: true })
    await w.settle()
    expect(panel()).toContain("ctrl+x again to stop write the report")
    expect(panel()).toMatch(/write the report\s+ctrl\+x again stops it/)
    expect(called(calls, "session.shutdown")).toEqual([])
    w.keys.pressKey("x", { ctrl: true })
    await w.settle()
    expect(called(calls, "session.shutdown")).toEqual([{ sessionId: "d1" }])
    w.store.push({ type: "session.removed", session: "d1", ts: 0, reason: "shutdown" } as never)
    await w.settle()
    expect(panel()).toContain("ended (shutdown)")
    w.keys.pressKey("x", { ctrl: true })
    await w.settle()
    expect(panel()).not.toContain("ended (shutdown)")
    expect(panel()).toContain("3 sessions")
  })

  test("ctrl+r renames with the name to edit; the name is not a prompt", async () => {
    const { w, calls } = await mount()
    await left()
    w.keys.pressArrow("up")
    await w.settle()
    w.keys.pressKey("r", { ctrl: true })
    await w.settle()
    expect(panel()).toContain("❯ write the report")
    expect(panel()).toContain("renaming: type the new name below")
    await w.keys.typeText(" v2")
    w.keys.pressEnter()
    await w.settle()
    expect(called(calls, "session.setTitle")).toEqual([{ sessionId: "d1", title: "write the report v2" }])
    expect(called(calls, "session.prompt")).toEqual([])
    expect(panel()).toContain("◆ AGENTS")
  })

  test("ctrl+t pins a session to the top, kept beside another window's pins", async () => {
    const dir = mkdtempSync(join(tmpdir(), "pins-"))
    dirs.push(dir)
    const file = join(dir, "tui_state.json")
    const pins = Pins.load(file)
    const { w } = await mount({ pins })
    // Another window pins a session meanwhile.
    writeFileSync(file, JSON.stringify({ pinned: ["elsewhere-pin"] }))
    await left()
    w.keys.pressArrow("up")
    await w.settle()
    w.keys.pressKey("t", { ctrl: true })
    await w.settle()
    const text = panel()
    expect(text.indexOf("PINNED")).toBeLessThan(text.indexOf("write the report"))
    expect(text.indexOf("write the report")).toBeLessThan(text.indexOf("NEEDS INPUT"))
    expect(JSON.parse(readFileSync(file, "utf8")).pinned.sort()).toEqual(["d1", "elsewhere-pin"])
  })

  test("other DELFIN sessions are listed to watch and message", async () => {
    const { w, calls } = await mount({
      client: ({ method }) => method === "sessions.list"
        ? { sessions: [], elsewhere: [{ key: "dash-7", title: "dashboard run", workspace: "/w/x", host: "login2" }] }
        : undefined,
    })
    await left()
    await Bun.sleep(10)
    await w.settle()
    expect(panel()).toContain("ELSEWHERE · watch and message")
    expect(panel()).toContain("dashboard run")
    w.keys.pressArrow("down")                           // the last row, after this window's
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(panel()).toContain("watched, not opened")
    w.keys.pressKey(" ")
    await w.settle()
    await w.keys.typeText("please pause")
    w.keys.pressEnter()
    await w.settle()
    expect(called(calls, "session.message")).toEqual([{ key: "dash-7", text: "please pause" }])
  })

  test("alt+n opens the n-th session; ctrl+c clears what is typed, then closes", async () => {
    const { w, opened } = await mount()
    await left()
    await w.keys.typeText("draft")
    w.keys.pressCtrlC()
    await w.settle()
    expect(panel()).toContain("◆ AGENTS")
    expect(panel()).not.toContain("❯ draft")
    w.keys.pressCtrlC()
    await w.settle()
    expect(panel()).not.toContain("◆ AGENTS")
    expect(w.quits()).toBe(0)
    await left()
    w.keys.pressKey("2", { meta: true })
    await w.settle()
    expect(opened).toEqual([["w1", "/work/w1"]])
  })

  test("with a question open, ← still opens the view, and letters do not answer it", async () => {
    const { w, calls } = await mount()
    w.store.push({ type: "permission.asked", session: "s1", ts: 0, seq: 5, requestId: "p1", kind: "confirm",
                   tool: "bash", preview: "$ rm -rf build", flags: {},
                   options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }] } as never)
    await w.settle()
    expect(panel()).toContain("rm -rf build")
    await left()
    await w.keys.typeText("yn")
    await w.settle()
    expect(called(calls, "permission.reply")).toEqual([])
    expect(panel()).toContain("❯ yn")
    w.keys.pressCtrlC()
    await w.escape()
    expect(panel()).toContain("rm -rf build")
  })

  test("another session waiting is said in the hints and the terminal title", async () => {
    const titles: string[] = []
    const { w, push } = await mount()
    const renderer = w.t.renderer as unknown as { setTerminalTitle(t: string): void }
    const original = renderer.setTerminalTitle.bind(renderer)
    renderer.setTerminalTitle = (t: string) => { titles.push(t); original(t) }
    push(summary("q2", { title: "second question", state: "needs_input", pending: 1, createdAt: 1004 }))
    await w.settle()
    expect(panel()).toContain("◉ 2 waiting  · ")
    expect(titles.at(-1)).toBe("◉ DELFIN · 2 need input")
    push(summary("q1", { title: "plan the scan", state: "idle", pending: 0, createdAt: 1002 }))
    push(summary("q2", { title: "second question", state: "idle", pending: 0, createdAt: 1004 }))
    await w.settle()
    expect(panel()).not.toContain("waiting ·")
    // Back to the session in view, by its name.
    expect(titles.at(-1)).toBe("DELFIN · this window")
  })

  test("a sleeping session is listed as asleep, and enter wakes it in this window", async () => {
    const { w, opened, push } = await mount()
    push(summary("z1", { title: "old scan", state: "sleeping", answer: "The scan finished.", createdAt: 1005 }))
    await w.settle()
    await left()
    expect(panel()).toMatch(/☾ old scan\s+asleep · The scan finished\./)
    expect(panel()).toContain("5 sessions · 1 needs input · 1 working")
    w.keys.pressArrow("up")                             // DONE: z1 is the newest there
    w.keys.pressArrow("up")
    await w.settle()
    expect(panel()).toMatch(/›\s+☾ old scan/)
    w.keys.pressEnter()
    await w.settle()
    expect(opened).toEqual([["z1", "/work/z1"]])
    expect(panel()).toContain("waking old scan…")
  })

  test("its keys can be rebound, and /keys lists them", () => {
    const km = new Keymap({ agentsPeek: ["ctrl+e"] })
    expect(km.matches("agentsPeek", { name: "e", ctrl: true, meta: false, shift: false })).toBe(true)
    expect(km.matches("agentsUp", { name: "k", ctrl: false, meta: false, shift: false })).toBe(false)
  })
})
