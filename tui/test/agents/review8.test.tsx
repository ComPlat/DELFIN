/**
 * The agent view's keys and state, as review 8 found them wanting --
 * each test was a failing repro before its fix. Headless: the real App on
 * OpenTUI's test renderer and a fake service, as test/agents does.
 */
import { afterEach, expect, test } from "bun:test"
import { ensureSession } from "../../src/store/state"
import { mountWindow, type Call } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
})

const summary = (id: string, extra: Record<string, unknown> = {}) => ({
  sessionId: id, title: null, workspace: `/work/${id}`, provider: "kit", model: "kit.test", approval: "default",
  state: "idle", pending: 0, costUsd: 0, createdAt: 1000, lastActivity: Date.now() / 1000, pid: 1,
  lastPrompt: null, activity: null, waiting: null, answer: null, ...extra,
})

async function mount(options: { client?: (call: Call) => unknown; onNewSession?: () => Promise<string | null> } = {}) {
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
      if (method === "sessions.list") return { sessions: [...running.values()], elsewhere: [] }
      if (method === "session.replay") return { events: [{ type: "user.message", text: `peeked ${params.sessionId}` }] }
      if (method === "session.setTitle") return { title: params.title }
      return { accepted: true, applied: true, message: "" }
    },
  }
  w = await mountWindow({ width: 100, height: 30, client, onResume: (id, ws) => { opened.push([id, ws]) },
                          onNewSession: options.onNewSession })
  const push = (s: Record<string, unknown>) => {
    running.set(String(s.sessionId), s)
    w!.store.push({ type: "session.summary", session: s.sessionId, ts: 0, summary: s } as never)
  }
  push(summary("s1", { title: "this window", createdAt: 900 }))
  push(summary("w1", { title: "optimise benzene", state: "running", activity: "bash  orca opt.inp", createdAt: 1001 }))
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

test("a session whose ended row was taken out and that runs again is listed again", async () => {
  const { w, push } = await mount()
  await left()
  w.keys.pressArrow("up")                       // write the report (d1)
  await w.settle()
  w.store.push({ type: "session.removed", session: "d1", ts: 0, reason: "shutdown" } as never)
  await w.settle()
  w.keys.pressKey("x", { ctrl: true })          // an ended row leaves the list
  await w.settle()
  expect(panel()).not.toContain("write the report")
  // The same session is resumed (by this window's /resume, or another window):
  // it runs again under its id and waits for an answer.
  push(summary("d1", { title: "write the report", state: "needs_input", pending: 1, waiting: "bash: make" }))
  await w.settle()
  const text = panel()
  console.log("R1 view:\n" + text)
  await w.escape()
  const hints = w.footer().join("\n")
  console.log("R1 hints after esc:\n" + hints)
  expect(hints).toMatch(/◉ 1 *$/m)              // the count at the right counts it ...
  expect(text).toContain("write the report")    // ... but the view does not show it
})

test("ctrl+c ends a rename: the next prompt renames nothing", async () => {
  const { w, calls } = await mount()
  await left()
  w.keys.pressArrow("up")                       // write the report
  await w.settle()
  w.keys.pressKey("r", { ctrl: true })
  await w.settle()
  w.keys.pressCtrlC()                           // "cancel": clears what is typed
  await w.settle()
  await w.keys.typeText("check the convergence")
  w.keys.pressEnter()
  await w.settle()
  console.log("R2 setTitle:", JSON.stringify(called(calls, "session.setTitle")))
  expect(called(calls, "session.setTitle")).toEqual([])
})

test("ctrl+r on an elsewhere row is the view's, never the history search", async () => {
  const { w } = await mount({
    client: ({ method }) => method === "sessions.list"
      ? { sessions: [], elsewhere: [{ key: "dash-7", title: "dashboard run", workspace: "/w/x", host: "login2" }] }
      : undefined,
  })
  // something in the prompt history
  await w.keys.typeText("rm the old runs")
  w.keys.pressEnter()
  await w.settle()
  await left()
  await Bun.sleep(10)
  await w.settle()
  w.keys.pressArrow("down")                     // round to ... the elsewhere row is last
  await w.settle()
  let tries = 0
  while (!/›\s+\S\s+dashboard run/.test(panel()) && tries++ < 6) {
    w.keys.pressArrow("down")
    await w.settle()
  }
  expect(panel()).toMatch(/›\s+\S\s+dashboard run/)
  w.keys.pressKey("r", { ctrl: true })
  await w.settle()
  const text = panel()
  console.log("R3 after ctrl+r on the elsewhere row:\n" + text.split("\n").slice(-6).join("\n"))
  expect(text).not.toContain("search history")
})

test("a reply typed for the peeked session stays with it when ↑ is pressed", async () => {
  const { w, calls } = await mount()
  await left()
  w.keys.pressArrow("up")                       // write the report (d1)
  await w.settle()
  w.keys.pressKey(" ")
  await Bun.sleep(10)
  await w.settle()
  await w.keys.typeText("add the figure")
  await w.settle()
  w.keys.pressArrow("up")                       // e.g. meant for the prompt's first line
  await Bun.sleep(10)
  await w.settle()
  w.keys.pressEnter()
  await w.settle()
  console.log("R4 prompts:", JSON.stringify(called(calls, "session.prompt")))
  expect(called(calls, "session.prompt")).toEqual([{ sessionId: "d1", text: "add the figure" }])
})

test("two new sessions started in a row: each prompt reaches its own session", async () => {
  // The same steps as main.tsx newSession, with a create that takes a moment.
  let n = 0
  const { w, calls } = await mount({
    onNewSession: async () => {
      const id = `n${++n}`
      await Bun.sleep(id === "n1" ? 40 : 80)
      w!.store.update((d) => { ensureSession(d, id, {}); d.current = id })
      return id
    },
  })
  await left()
  await w.keys.typeText("job one")
  w.keys.pressEnter()
  await w.settle()
  await left()
  await w.keys.typeText("job two")
  w.keys.pressEnter()
  await Bun.sleep(150)
  await w.settle()
  console.log("R5 prompts:", JSON.stringify(called(calls, "session.prompt")))
  expect(called(calls, "session.prompt")).toEqual([{ sessionId: "n1", text: "job one" },
                                                   { sessionId: "n2", text: "job two" }])
})

test("a new session that does not start, while something else is typed: its prompt is kept", async () => {
  const { w, calls } = await mount({
    onNewSession: async () => {
      await Bun.sleep(60)
      return null
    },
  })
  await left()
  await w.keys.typeText("start the geometry optimisation")
  w.keys.pressEnter()
  await w.settle()
  await w.keys.typeText("x")                    // typing goes on while it starts
  await Bun.sleep(100)
  await w.settle()
  const text = panel()
  console.log("R6 composer after the failure:", text.split("\n").filter((l) => l.startsWith("❯")))
  expect(text).toContain("start the geometry optimisation")
})

test("a new session that does not start: its prompt comes back in the view, never to the session in view", async () => {
  const { w, calls } = await mount({ onNewSession: async () => null })
  await left()
  await w.keys.typeText("start the opt")
  w.keys.pressEnter()
  await Bun.sleep(20)
  await w.settle()
  const open = panel().includes("◆ AGENTS")
  w.keys.pressEnter()                           // try again
  await w.settle()
  console.log("R7 view open after the failure:", open, "| prompts:", JSON.stringify(called(calls, "session.prompt")))
  expect(called(calls, "session.prompt")).toEqual([])
})

test("a session that ended with /exit elsewhere, resumed in this window, does not make it leave", async () => {
  const { w, push } = await mount()
  push(summary("b1", { title: "bee", createdAt: 1004 }))
  await w.settle()
  // Another window typed /exit in it: the service says it was removed, reason exit.
  w.store.push({ type: "session.removed", session: "b1", ts: 0, reason: "exit", exitCode: 0 } as never)
  await w.settle()
  // Resumed from this window's agent view (its ended row): main.tsx switchTo
  // registers it (ensureSession; state from the attach) and shows it.
  w.store.update((d) => {
    const s = ensureSession(d, "b1", { title: "bee" })
    s.state = "idle"
    s.pending = []
  })
  w.view.show("b1", { divider: "resumed: bee" })
  await w.settle()
  console.log("R8 quits:", w.quits(), "| ended kept:", JSON.stringify(w.store.state.sessions.b1!.ended))
  expect(w.quits()).toBe(0)
})

test("ctrl+r while a prompt is typed in the view keeps the draft for after the rename", async () => {
  const { w, calls } = await mount()
  await left()
  await w.keys.typeText("compare the two conformers")   // a new session's prompt, not sent yet
  await w.settle()
  w.keys.pressKey("r", { ctrl: true })                   // the history search in the normal prompt
  await w.settle()
  const during = panel().split("\n").filter((l) => l.startsWith("❯"))
  await w.escape()                                       // leave the rename
  const after = panel().split("\n").filter((l) => l.startsWith("❯"))
  console.log("R9 prompt row during the rename:", during, "| after esc:", after)
  expect(after.join("")).toContain("compare the two conformers")
})

test("while a new session starts, the other sessions' events are not held", async () => {
  const { eventsOfUnknown, withReplay } = await import("../../src/store/store")
  let releaseCreate!: () => void
  const { w, push } = await mount({
    onNewSession: async () => {
      // main.tsx newSession: withReplay(store, session.create, register, eventsOfUnknown)
      await withReplay(w!.store, () => new Promise<any>((resolve) => {
        releaseCreate = () => resolve({ sessionId: "n1", session: {}, state: "idle", pending: [], events: [] })
      }), (answer) => {
        w!.store.update((d) => { ensureSession(d, answer.sessionId, {}); d.current = answer.sessionId })
        return []
      }, eventsOfUnknown(w!.store))
      return "n1"
    },
  })
  await left()
  await w.keys.typeText("a slow-starting session")
  w.keys.pressEnter()
  await w.settle()
  // Meanwhile another session asks a question (its summary arrives).
  push(summary("q9", { title: "needs you", state: "needs_input", pending: 1, waiting: "bash: rm" }))
  await w.settle()
  const during = w.footer().join("\n")
  const seenDuring = Object.keys(w.store.state.sessions).includes("q9")
  releaseCreate()
  await Bun.sleep(30)
  await w.settle()
  const seenAfter = Object.keys(w.store.state.sessions).includes("q9")
  console.log(`R10 badge while the create is in flight: ${/◉ \d+ waiting/.test(during)}; q9 known: during ${seenDuring}, after ${seenAfter}`)
  expect(seenDuring).toBe(true)
})

test("a new session asked for while the first one starts follows it, with its own prompt", async () => {
  let started = 0
  const { w, calls } = await mount({
    onNewSession: async () => {
      started++
      w!.store.update((d) => { ensureSession(d, "n1", {}); d.current = "n1" })
      return "n1"
    },
  })
  w.store.update((d) => { d.starting = true })       // the first session is not up yet
  await left()
  await w.keys.typeText("job from the view")
  w.keys.pressEnter()
  await Bun.sleep(150)
  await w.settle()
  expect(started).toBe(0)                            // it waits for the first
  w.store.update((d) => { d.starting = false })
  await Bun.sleep(250)
  await w.settle()
  expect(started).toBe(1)
  expect(called(calls, "session.prompt")).toEqual([{ sessionId: "n1", text: "job from the view" }])
})
