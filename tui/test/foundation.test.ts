/**
 * The window's foundation, piece by piece: the keymap, the config file, the
 * store (reducer, batching, replay deduplication), the JSON-RPC framing and
 * the presenter's summaries.
 */
import { describe, expect, test } from "bun:test"
import { DEFAULTS, parseConfig } from "../src/config/config"
import { Keymap, describe as describeKey, parseBinding } from "../src/keymap/keymap"
import { askedWhat } from "../src/scrollback/present"
import { applyEvent, emptyState, ensureSession } from "../src/store/state"
import { createAppStore, type Clock } from "../src/store/store"
import { RpcClient, RpcError, type Wire } from "../src/transport/rpc"

const key = (name: string, mods: { ctrl?: boolean; shift?: boolean; meta?: boolean } = {}) =>
  ({ name, ctrl: Boolean(mods.ctrl), shift: Boolean(mods.shift), meta: Boolean(mods.meta), option: false })

describe("the keymap", () => {
  test("has today's keys by default", () => {
    const km = new Keymap()
    expect(km.matches("steer", key("g", { ctrl: true }))).toBe(true)
    expect(km.matches("cycleMode", key("tab", { shift: true }))).toBe(true)
    expect(km.matches("interrupt", key("escape"))).toBe(true)
    expect(km.matches("agents", key("left"))).toBe(true)
    expect(km.matches("steer", key("g"))).toBe(false)              // ctrl matters
  })

  test("takes the user's bindings and reports the ones it cannot read", () => {
    const km = new Keymap({ steer: ["alt+s"], nonsense: ["x"], quit: ["super+q"] })
    expect(km.matches("steer", key("s", { meta: true }))).toBe(true)
    expect(km.matches("steer", key("g", { ctrl: true }))).toBe(false)
    expect(km.matches("quit", key("d", { ctrl: true }))).toBe(true)  // the bad one is not taken
    expect(km.problems.join(" ")).toContain("nonsense")
    expect(km.problems.join(" ")).toContain("super+q")
  })

  test("reads and writes bindings the way a hint shows them", () => {
    expect(parseBinding("Ctrl+Shift+Enter")).toEqual({ name: "return", ctrl: true, alt: false, shift: true })
    expect(describeKey(parseBinding("escape")!)).toBe("esc")
    expect(describeKey(parseBinding("ctrl+g")!)).toBe("ctrl+g")
    expect(new Keymap().label("agents")).toBe("←")
  })
})

describe("the config file", () => {
  test("defaults when empty", () => {
    expect(parseConfig("{}").config).toEqual({ ...DEFAULTS, colors: {}, keys: {} })
  })

  test("keeps what is right and names what is wrong", () => {
    const { config, problems } = parseConfig(JSON.stringify({
      theme: "light", colors: { accent: "#112233", bad: "blue" }, keys: { steer: "alt+s" },
      replayTurns: 5, showThinking: "yes", unknown: 1,
    }))
    expect(config.theme).toBe("light")
    expect(config.colors).toEqual({ accent: "#112233" })
    expect(config.keys).toEqual({ steer: ["alt+s"] })
    expect(config.replayTurns).toBe(5)
    expect(config.showThinking).toBe(false)
    expect(problems.length).toBe(3)
  })

  test("a file that is not JSON costs its settings, not the start", () => {
    const { config, problems } = parseConfig("{ not json")
    expect(config.theme).toBe("auto")
    expect(problems[0]).toContain("not JSON")
  })
})

describe("the store", () => {
  test("follows a turn: state, label, tail, questions", () => {
    const state = emptyState()
    ensureSession(state, "s1", { model: "m", approval: "default" })
    const e = (type: string, fields: Record<string, unknown> = {}) =>
      applyEvent(state, { type, session: "s1", ts: 0, ...fields } as never, 1000)
    e("session.state", { state: "running" })
    e("turn.started", { turnId: "t1", source: "user" })
    expect(state.sessions.s1!.turn?.label).toBe("thinking…")
    e("tool.started", { name: "bash", input: { command: "pytest -q\nmore" } })
    // The REPL's headline: a command's lines run together, on one line.
    expect(state.sessions.s1!.turn?.label).toBe("bash  pytest -q more")
    e("permission.asked", { requestId: "p1", kind: "confirm", tool: "bash", preview: "$ ls", options: [], flags: {} })
    e("permission.asked", { requestId: "p1", kind: "confirm", tool: "bash", preview: "$ ls", options: [], flags: {} })
    expect(state.sessions.s1!.pending.length).toBe(1)
    e("permission.resolved", { requestId: "p1", decision: "approved", by: "client" })
    expect(state.sessions.s1!.pending.length).toBe(0)
    e("usage", { inputTokens: 10, outputTokens: 5, costUsd: 0.01 })
    e("turn.ended", { turnId: "t1", durationMs: 5, inputTokens: 1, outputTokens: 1, costUsd: 0.5, toolCalls: 1 })
    expect(state.sessions.s1!.turn).toBe(null)
    expect(state.sessions.s1!.usage.costUsd).toBe(0.01)        // totals, not summed twice
  })

  test("a summary registers a session; nothing else does", () => {
    const state = emptyState()
    applyEvent(state, { type: "turn.started", session: "nobody", ts: 0, turnId: "t1", source: "user" } as never)
    expect(state.order).toEqual([])
    applyEvent(state, { type: "session.summary", session: "s2", ts: 0,
                        summary: { sessionId: "s2", state: "idle", workspace: "/w", provider: "p", model: "m",
                                   approval: "plan", pending: 0, costUsd: 0, createdAt: 0, lastActivity: 0, pid: 1 } } as never)
    expect(state.order).toEqual(["s2"])
    expect(state.sessions.s2!.approval).toBe("plan")
  })

  test("applies events in one batch per frame, and drops a durable event seen before", () => {
    let pending: (() => void) | null = null
    const clock: Clock = { setTimeout: (fn) => { pending = fn; return 1 }, clearTimeout: () => { pending = null }, now: () => 0 }
    const store = createAppStore(clock, 16)
    store.update((d) => { ensureSession(d, "s1") })
    const batches: number[] = []
    store.onBatch((events) => batches.push(events.length))
    store.push({ type: "user.message", session: "s1", ts: 0, seq: 1, text: "a", source: "user" } as never)
    store.push({ type: "text.delta", session: "s1", ts: 0, partId: "p", delta: "x" } as never)
    store.push({ type: "user.message", session: "s1", ts: 0, seq: 1, text: "a", source: "user" } as never)
    expect(batches).toEqual([])                  // nothing until the frame
    ;(pending as unknown as () => void)()
    expect(batches).toEqual([2])                 // the repeat of seq 1 was dropped
    expect(store.state.sessions.s1!.lastSeq).toBe(1)
  })
})

describe("the JSON-RPC framing", () => {
  function pair() {
    const sent: string[] = []
    const wire: Wire = { write: (line) => sent.push(line), close: () => {} }
    return { rpc: new RpcClient(wire), sent }
  }

  test("matches answers to requests, however the bytes are split", async () => {
    const { rpc, sent } = pair()
    const a = rpc.request("sessions.list", {})
    const b = rpc.request("hello", { protocolVersion: 1 })
    const ids = sent.map((line) => JSON.parse(line).id)
    const bytes = new TextEncoder().encode(
      `{"jsonrpc":"2.0","id":${ids[1]},"result":{"ok":"Energie ≈ −76.4 Eₕ"}}\n` +
      `{"jsonrpc":"2.0","id":${ids[0]},"error":{"code":-32003,"message":"no"}}\n`)
    for (let i = 0; i < bytes.length; i += 7) rpc.feed(bytes.subarray(i, i + 7))   // split mid-character
    expect((await b) as unknown).toEqual({ ok: "Energie ≈ −76.4 Eₕ" })
    const err = await a.catch((e) => e)
    expect(err).toBeInstanceOf(RpcError)
    expect((err as RpcError).code).toBe(-32003)
  })

  test("hands events to every listener, and survives a bad line", () => {
    const { rpc } = pair()
    const got: string[] = []
    rpc.onEvent((e) => got.push(String((e as { type: string }).type)))
    rpc.feed('not json\n{"jsonrpc":"2.0","method":"event","params":{"type":"notice","session":"s"}}\n')
    expect(got).toEqual(["notice"])
  })

  test("fails whatever waits when the connection is lost", async () => {
    const { rpc } = pair()
    const waiting = rpc.request("sessions.list", {})
    rpc.lost("the service went away")
    const err = (await waiting.catch((e) => e)) as RpcError
    expect(err.message).toContain("went away")
    await expect(rpc.request("sessions.list", {})).rejects.toThrow("not connected")
  })
})

describe("what the presenter says about questions", () => {
  test("an approval is named by its command, its file or its question", () => {
    expect(askedWhat({ kind: "confirm", tool: "bash", preview: "$ make deploy", flags: { shell: true } })).toBe("bash make deploy")
    expect(askedWhat({ kind: "confirm", tool: "edit_file", preview: "--- a/x.py\n+++ b/x.py\n@@", flags: {} })).toBe("edit_file x.py")
    expect(askedWhat({ kind: "ask", question: "Which solvent?" })).toBe("Which solvent?")
    expect(askedWhat({ kind: "plan" })).toBe("the plan")
  })
})
