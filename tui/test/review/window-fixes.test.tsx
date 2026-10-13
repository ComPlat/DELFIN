/**
 * What the foundation review found, kept fixed: the window shows what the
 * session did and answers only the question it shows.
 *
 * - a replay is applied before the live events that came with it;
 * - an answer whose deltas went missing is committed from its whole text;
 * - the approval panel follows the open question, its keys belong to it,
 *   Esc refuses, ctrl+c and ctrl+d still work, typeahead is not an answer;
 * - a prompt that could not be sent stays in the composer;
 * - the footer fits the terminal at any width, also after a resize;
 * - a session whose worker died is shown as ended, not running.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { ensureSession } from "../../src/store/state"
import { withReplay } from "../../src/store/store"
import { RpcClient } from "../../src/transport/rpc"
import { FIXTURE, mountWindow } from "../harness"

type Window = Awaited<ReturnType<typeof mountWindow>>
let open: Window | null = null
afterEach(async () => {
  await open?.close()
  open = null
})
async function mount(options: Parameters<typeof mountWindow>[0] = {}) {
  open = await mountWindow(options)
  return open
}

const edit = FIXTURE.events.find((e: any) => e.type === "permission.asked" && e.kind === "confirm")
const upToEdit = FIXTURE.events.slice(0, FIXTURE.events.indexOf(edit) + 1)
const bash = (id: string, seq: number, cmd: string, extra: Record<string, unknown> = {}) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "confirm", tool: "bash",
  preview: `$ ${cmd}`, flags: { shell: true, write: false, protected: false, outsideRead: false },
  options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "A", label: "always this exact command" },
            { key: "a", label: "abort" }, { key: "?", label: "explain" }],
  ...extra,
})
const resolved = (id: string, seq: number) =>
  ({ type: "permission.resolved", session: "s1", ts: 0, seq, requestId: id, decision: "approved", by: "outside" })
const multi = (id: string, seq: number, question: string, choices: string[]) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "ask", tool: "ask_user_question",
  preview: "", options: [], flags: {}, question, header: "Pick", choices, choiceDetails: choices.map(() => ""),
  multiSelect: true,
})

describe("a replay and the live events behind it", () => {
  const header = { type: "session.info", session: "s1", ts: 0, seq: 0, ...FIXTURE.session }
  const replay = [header, ...FIXTURE.events.filter((e: any) => typeof e.seq === "number" && e.seq <= 19)]
  const expected = ["❯ hello there", "The workspace holds", "❯ please edit the parser", "⏺ edit_file", "regular expression"]

  async function resumeWith(live: object[]) {
    const w = await mount({ sid: "other" })
    const sent: string[] = []
    const rpc = new RpcClient({ write: (l) => sent.push(l), close: () => {} })
    rpc.onEvent((event) => w.store.push(event))
    const attached = withReplay(w.store, () => rpc.request("session.attach", { sessionId: "s1", turns: 20 }),
      (answer: any) => {
        w.store.update((d) => { ensureSession(d, "s1", answer.session) })
        w.view.show("s1")
        return answer.events
      })
    const id = JSON.parse(sent[0]!).id
    const answer = { jsonrpc: "2.0", id, result: { session: FIXTURE.session, state: "running", pending: [],
                                                   events: replay, summary: {} } }
    // The answer and, right behind it, live events: one read.
    rpc.feed(JSON.stringify(answer) + "\n"
             + live.map((p) => JSON.stringify({ jsonrpc: "2.0", method: "event", params: p }) + "\n").join(""))
    await attached
    await w.settle()
    return w.t.externalOutput.takeText()
  }

  test("a live durable event in the same read does not hide the replay", async () => {
    const text = await resumeWith([{ type: "user.message", session: "s1", ts: 0, seq: 20,
                                     text: "ask me which solvent", source: "user" }])
    for (const part of [...expected, "❯ ask me which solvent"]) expect(text).toContain(part)
    expect(text.indexOf("◆ DELFIN")).toBeLessThan(text.indexOf("❯ hello there"))
  })

  test("a live delta in the same read does not take the replayed answer's place", async () => {
    const text = await resumeWith([{ type: "text.delta", session: "s1", ts: 0, partId: "t3.p9", delta: "LIVE " }])
    for (const part of expected) expect(text).toContain(part)
  })
})

describe("an answer whose deltas went missing", () => {
  test("is committed from its whole text", async () => {
    const w = await mount()
    w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
    w.store.push({ type: "text.started", session: "s1", ts: 0, turnId: "t1", partId: "p1" } as never)
    w.store.push({ type: "text.delta", session: "s1", ts: 0, partId: "p1", delta: "Working through a long calcul" } as never)
    await w.settle()
    // The window reconnects: the replay ends the part, then live deltas of it arrive.
    w.store.push({ type: "text.ended", session: "s1", ts: 0, seq: 2, turnId: "t1", partId: "p1",
                   text: "Working through a long calculation, step by step. Four sentences in all." } as never)
    w.store.push({ type: "text.delta", session: "s1", ts: 0, partId: "p1", delta: " DUPLICATE" } as never)
    await w.settle()
    const text = w.t.externalOutput.takeText()
    expect(text).toContain("Four sentences in all.")
    expect(text).not.toContain("DUPLICATE")
  })
})

describe("the approval panel", () => {
  for (const how of ["two questions open at once", "resolved and asked in one batch"]) {
    test(`follows the open question (${how}), and its keys answer that one`, async () => {
      const w = await mount({ height: 40 })
      for (const e of upToEdit) w.store.push(e)
      if (how === "two questions open at once") w.store.push(bash("p-bash", 100, "rm -rf build/ dist/") as never)
      await w.settle()
      expect(w.footer().join("\n")).toContain("FINAL SINGLE POINT ENERGY")
      w.store.push(resolved(edit.requestId, 101) as never)
      if (how !== "two questions open at once") w.store.push(bash("p-bash", 102, "rm -rf build/ dist/") as never)
      await w.settle()
      const shown = w.footer().join("\n")
      expect(shown).toContain("rm -rf build")
      expect(shown).not.toContain("FINAL SINGLE POINT ENERGY")
      w.t.mockInput.pressKey("y")
      await w.settle()
      expect(w.replies().at(-1)).toMatchObject({ requestId: "p-bash", choice: "y" })
    })
  }

  test("esc refuses, as in the REPL", async () => {
    const w = await mount()
    for (const e of upToEdit) w.store.push(e)
    await w.settle()
    w.t.mockInput.pressEscape()
    await Bun.sleep(40)
    await w.settle()
    expect(w.replies()).toEqual([{ sessionId: "s1", requestId: edit.requestId, choice: "n" }])
  })

  test("N refuses with a reason the model reads", async () => {
    const w = await mount()
    for (const e of upToEdit) w.store.push(e)
    await w.settle()
    w.t.mockInput.pressKey("N", { shift: true })
    await w.settle()
    expect(w.footer().join("\n")).toContain("why?")
    await w.t.mockInput.typeText("use the def2 parser")
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(w.replies()).toEqual([{ sessionId: "s1", requestId: edit.requestId, choice: "n",
                                   reason: "use the def2 parser" }])
  })

  test("? explains every option in the panel", async () => {
    const w = await mount({ height: 40 })
    for (const e of upToEdit) w.store.push(e)
    await w.settle()
    w.t.mockInput.pressKey("?")
    await w.settle()
    const shown = w.footer().join("\n")
    for (const option of edit.options) expect(shown).toContain(`[${option.key}] ${option.label}`)
    expect(w.replies()).toEqual([])
  })

  test("ctrl+c ends the turn and ctrl+d leaves, with a question open", async () => {
    const w = await mount()
    for (const e of upToEdit) w.store.push(e)
    await w.settle()
    w.t.mockInput.pressKey("c", { ctrl: true })
    await w.settle()
    expect(w.calls.some((c) => c.method === "session.cancel")).toBe(true)
    w.t.mockInput.pressKey("d", { ctrl: true })
    expect(w.quits()).toBe(1)
    expect(w.replies()).toEqual([])                     // never answered for the user
  })

  test("picks made for one question are not sent for the next", async () => {
    const w = await mount()
    w.store.push(multi("q1", 100, "Which properties?", ["IR", "NMR", "UV-Vis"]) as never)
    await w.settle()
    w.t.mockInput.pressKey("1")
    w.t.mockInput.pressKey("3")
    await w.settle()
    w.store.push({ type: "permission.resolved", session: "s1", ts: 0, seq: 101, requestId: "q1",
                   decision: "dismissed", by: "abort" } as never)
    w.store.push(multi("q2", 102, "Which files may I delete?", ["build/", "dist/", "src/"]) as never)
    await w.settle()
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(w.replies().at(-1)).toEqual({ sessionId: "s1", requestId: "q2", answers: [] })
  })

  test("the persist step of A belongs to its question, and only y persists", async () => {
    const w = await mount()
    w.store.push(bash("c1", 100, "pytest -q") as never)
    await w.settle()
    w.t.mockInput.pressKey("A", { shift: true })
    await w.settle()
    expect(w.footer().join("\n")).toContain("persist")
    w.store.push(resolved("c1", 101) as never)
    w.store.push(bash("c2", 102, "rm -rf ~/runs") as never)
    await w.settle()
    w.t.mockInput.pressKey("n")
    w.t.mockInput.pressKey("y")                         // answered already: not sent again
    await w.settle()
    expect(w.replies().at(-1)).toEqual({ sessionId: "s1", requestId: "c2", choice: "n" })
    expect(w.replies().filter((r) => r.requestId === "c2").length).toBe(1)
    w.store.push(resolved("c2", 103) as never)          // what the worker says next
    w.store.push(bash("c3", 104, "pytest -q") as never)
    await w.settle()
    w.t.mockInput.pressKey("A", { shift: true })
    w.t.mockInput.pressKey("y")
    await w.settle()
    expect(w.replies().at(-1)).toEqual({ sessionId: "s1", requestId: "c3", choice: "A", persist: true })
  })

  test("keys typed on when a question opens are not answers", async () => {
    const w = await mount({ answerPauseMs: 350 })
    for (const e of FIXTURE.events.slice(0, FIXTURE.events.indexOf(edit))) w.store.push(e)
    await w.settle()
    await w.t.mockInput.typeText("pl")
    w.store.push(edit)
    await w.settle()
    await w.t.mockInput.typeText("ease rerun")
    await w.settle()
    expect(w.replies()).toEqual([])
    expect(w.footer().join("\n")).toContain("were not used as answers")
    await Bun.sleep(400)                                  // a pause, then a deliberate key
    w.t.mockInput.pressKey("y")
    await w.settle()
    expect(w.replies()).toEqual([{ sessionId: "s1", requestId: edit.requestId, choice: "y" }])
  })
})

describe("the footer", () => {
  for (const width of [100, 80, 60]) {
    test(`fits the edit approval at ${width} columns`, async () => {
      const w = await mount({ width, height: 40 })
      for (const e of upToEdit) w.store.push(e)
      await w.settle()
      const footer = w.footer()
      expect(footer.some((l) => l.includes("[y] yes"))).toBe(true)
      expect(footer.some((l) => l.includes("❯"))).toBe(true)
      expect(footer.some((l) => l.includes("shift+tab"))).toBe(true)
      expect(footer.some((l) => l.includes("FINAL SINGLE POINT ENERGY") || l.includes("@@"))).toBe(true)
      for (const line of footer) expect(Bun.stringWidth(line.trimEnd())).toBeLessThanOrEqual(width)
    })
  }

  test("follows a resize: rules as wide as the terminal, nothing drawn over", async () => {
    const w = await mount({ width: 100, height: 30 })
    for (const e of FIXTURE.events.slice(0, 30)) w.store.push(e)
    await w.settle()
    w.t.resize(70, 30)
    await Bun.sleep(400)
    await w.settle()
    const footer = w.footer()
    const rules = footer.filter((l) => /^─{10,}/.test(l))
    expect(rules.length).toBe(2)
    for (const rule of rules) expect(rule.trimEnd().length).toBe(69)
    expect(footer.some((l) => l.startsWith("❯"))).toBe(true)
    expect(footer.some((l) => l.includes("shift+tab"))).toBe(true)
  })
})

describe("the composer when a prompt cannot be sent", () => {
  test("keeps the prompt", async () => {
    const rpc = new RpcClient(null)                       // the wire is gone (reconnecting)
    const w = await mount({ client: rpc })
    w.store.update((d) => { d.connection = "reconnecting" })
    await w.settle()
    await w.t.mockInput.typeText("please keep this long message")
    w.t.mockInput.pressEnter()
    await Bun.sleep(30)
    await w.settle()
    expect(w.footer().join("\n")).toContain("❯ please keep this long message")
  })
})

describe("a session whose worker died", () => {
  test("is shown as ended, and ctrl+c leaves instead of cancelling", async () => {
    const w = await mount()
    w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
    w.store.push({ type: "tool.started", session: "s1", ts: 0, seq: 2, turnId: "t1", callId: "c1", name: "bash",
                   input: { command: "sleep 20" } } as never)
    w.store.push(bash("p1", 3, "sleep 20") as never)
    await w.settle()
    w.store.push({ type: "session.removed", session: "s1", ts: 0, reason: "the worker was killed (SIGKILL)", exitCode: -9 } as never)
    await w.settle()
    const s = w.store.state.sessions.s1!
    expect(s.turn).toBe(null)
    expect(s.pending).toEqual([])
    expect(w.footer().join("\n")).not.toContain("to interrupt")
    expect(w.footer().join("\n")).toContain("session ended")
    expect(w.t.externalOutput.takeText()).toContain("session ended (the worker was killed (SIGKILL))")
    w.t.mockInput.pressKey("c", { ctrl: true })
    w.t.mockInput.pressKey("c", { ctrl: true })
    expect(w.quits()).toBe(1)
    expect(w.calls.some((c) => c.method === "session.cancel")).toBe(false)
  })

  test("a failed turn says its error once", async () => {
    const w = await mount()
    w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
    w.store.push({ type: "error", session: "s1", ts: 0, seq: 2, text: "scripted failure from the fake engine" } as never)
    w.store.push({ type: "turn.ended", session: "s1", ts: 0, seq: 3, turnId: "t1", durationMs: 400, inputTokens: 1,
                   outputTokens: 0, costUsd: 0, toolCalls: 0, error: "scripted failure from the fake engine" } as never)
    await w.settle()
    const text = w.t.externalOutput.takeText()
    expect(text.split("scripted failure").length - 1).toBe(1)
    expect(text).toContain("✗ 0.4s")
  })
})

describe("notes of the window's own", () => {
  test("are rows in the scrollback", async () => {
    const w = await mount()
    w.view.note("~/.delfin/tui.json: theme must be \"auto\", \"dark\" or \"light\"")
    await w.settle()
    expect(w.t.externalOutput.takeText()).toContain("theme must be")
  })
})
