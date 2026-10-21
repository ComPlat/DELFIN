/**
 * The command palette and the native pickers: ctrl+p lists the window's
 * actions and the session's commands; the pickers list what the session
 * offers and apply a pick the way a typed command would. A bare /model
 * opens the model picker. After a resize the last turns are written again,
 * or, under GNU screen, a divider says the width changed.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { createPicker, matches } from "../../src/footer/picker"
import { Transcript, type Item } from "../../src/scrollback/transcript"
import { CommitQueue } from "../../src/scrollback/queue"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const key = (name: string, extra: Record<string, unknown> = {}) =>
  ({ name, ctrl: false, meta: false, shift: false, option: false, sequence: name.length === 1 ? name : "", ...extra }) as never

describe("a picker", () => {
  test("narrows by every word typed, selects, picks and closes", async () => {
    const picker = createPicker()
    const picked: string[] = []
    picker.open("model", Promise.resolve([
      { label: "kit.deepseek-v4-flash", detail: "KIT", value: "a", current: true },
      { label: "kit.glm-5.3-flash", detail: "KIT", value: "b" },
      { label: "azure.gpt-5.6-sol", detail: "Azure (EU) · $5.50 in", value: "c" },
    ]), (c) => picked.push(c.value))
    expect(picker.state.loading).toBe(true)
    await Bun.sleep(1)
    expect(picker.state.selected).toBe(0)                 // the current one first
    for (const ch of "glm") picker.handleKey(key(ch))
    expect(picker.shown().map((c) => c.value)).toEqual(["b"])
    picker.handleKey(key("return"))
    expect(picked).toEqual(["b"])
    expect(picker.state.open).toBe(false)
    expect(matches("azure.gpt-5.6-sol Azure (EU) · $5.50 in", "eu sol")).toBe(true)
  })

  test("esc closes it, and a list that fails says why", async () => {
    const picker = createPicker()
    picker.open("model", Promise.reject(new Error("the provider did not answer")), () => {})
    await Bun.sleep(1)
    expect(picker.state.error).toContain("did not answer")
    picker.handleKey(key("escape"))
    expect(picker.state.open).toBe(false)
  })
})

function client(calls: { method: string; params: any }[]) {
  return {
    request: async (method: string, params: any = {}) => {
      calls.push({ method, params })
      if (method === "session.commands") {
        return { commands: [{ name: "/help", summary: "List these commands", category: "session", takesArgs: true, source: "builtin" },
                            { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" }] }
      }
      if (method === "session.choices") {
        return { current: "kit.deepseek-v4-flash", choices: [{ value: "kit.deepseek-v4-flash", detail: "KIT" },
                                                            { value: "kit.glm-5.3-flash", detail: "KIT" }] }
      }
      if (method === "sessions.past") {
        return { sessions: [{ sessionId: "old1", title: "the parser fix", workspace: "/w", model: "kit.x",
                              updatedAt: Date.now() / 1000 - 7200, messages: 12, costUsd: 0, running: false }] }
      }
      return { accepted: true, approval: params.mode }
    },
  }
}

describe("the palette and the pickers in the window", () => {
  test("ctrl+p lists actions and commands; model… opens the model picker; a pick sends /model", async () => {
    const calls: { method: string; params: any }[] = []
    const w = open = await mountWindow({ client: client(calls), height: 40 })
    w.t.mockInput.pressKey("p", { ctrl: true })
    await w.settle()
    let footer = w.footer().join("\n")
    expect(footer).toContain("commands")
    expect(footer).toContain("model…")
    await w.t.mockInput.typeText("model")
    w.t.mockInput.pressEnter()
    await w.settle()
    footer = w.footer().join("\n")
    expect(footer).toContain("kit.glm-5.3-flash")
    expect(footer).toContain("(now)")
    w.t.mockInput.pressArrow("down")
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)).toEqual(["/model kit.glm-5.3-flash"])
  })

  test("a bare /model opens the picker instead of asking the session", async () => {
    const calls: { method: string; params: any }[] = []
    const w = open = await mountWindow({ client: client(calls), height: 40 })
    await w.t.mockInput.typeText("/model")
    w.t.mockInput.pressEscape()                           // close the completion list
    await Bun.sleep(40)
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(w.footer().join("\n")).toContain("kit.deepseek-v4-flash")
    expect(calls.some((c) => c.method === "session.prompt")).toBe(false)
  })

  test("a command without arguments runs from the palette; one with them goes into the prompt", async () => {
    const calls: { method: string; params: any }[] = []
    const w = open = await mountWindow({ client: client(calls), height: 40 })
    w.t.mockInput.pressKey("p", { ctrl: true })
    await w.settle()
    await w.t.mockInput.typeText("/status")
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)).toEqual(["/status"])
    w.t.mockInput.pressKey("p", { ctrl: true })
    await w.settle()
    await w.t.mockInput.typeText("/help")
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(w.footer().join("\n")).toContain("❯ /help")
  })

  test("permissions… sets the posture; resume… shows a saved session", async () => {
    const calls: { method: string; params: any }[] = []
    const w = open = await mountWindow({ client: client(calls), height: 40 })
    w.t.mockInput.pressKey("p", { ctrl: true })
    await w.t.mockInput.typeText("permissions")
    w.t.mockInput.pressEnter()
    await w.settle()
    await w.t.mockInput.typeText("accept")
    w.t.mockInput.pressEnter()
    await w.settle()
    expect(calls.at(-1)).toEqual({ method: "session.setMode", params: { sessionId: "s1", mode: "acceptEdits" } })
    w.t.mockInput.pressKey("p", { ctrl: true })
    await w.t.mockInput.typeText("resume")
    w.t.mockInput.pressEnter()
    await w.settle()
    const footer = w.footer().join("\n")
    expect(footer).toContain("the parser fix")
    expect(footer).toContain("12 messages")
  })
})

describe("after a resize", () => {
  function setup() {
    const written: string[] = []
    const target = { requestLive() {}, dropLive() {}, resetSplitFooterForReplay() { written.push("<wipe>") } }
    const queue = new CommitQueue(target)
    const transcript = new Transcript(queue, target, (item: Item) => { written.push(String(item.text ?? item.kind)) })
    return { written, transcript }
  }

  test("the last turns are written again, after their banner and a note of what was left out", async () => {
    const { written, transcript } = setup()
    transcript.append({ id: "b", kind: "banner", text: "banner", group: 0 })
    for (let g = 1; g <= 5; g++) {
      transcript.append({ id: `u${g}`, kind: "user", text: `turn ${g}`, group: g })
      transcript.append({ id: `a${g}`, kind: "markdown", text: `answer ${g}`, group: g })
    }
    await transcript.replay(2, "turns", (cut, first) => [first!, { id: "n", kind: "notice", text: `${cut} left out` }])
    expect(written.slice(written.indexOf("<wipe>"))).toEqual(
      ["<wipe>", "banner", "3 left out", "turn 4", "answer 4", "turn 5", "answer 5"])
  })

  test("under GNU screen a divider says so, and nothing is wiped", async () => {
    const w = open = await mountWindow({ width: 100 })
    w.view.stop()                                          // the harness's own view: replaced below
    const { startSessionView } = await import("../../src/session")
    const { DARK } = await import("../../src/theme/theme")
    const view = startSessionView(w.t.renderer, w.store, () => DARK, {
      showThinking: false, replayTurns: 20, dividerOnResize: true })
    const before = view.transcript.replayCount
    w.t.resize(80, 30)
    await Bun.sleep(400)
    await view.queue.drain()
    expect(view.transcript.replayCount).toBe(before)
    expect(view.transcript.items.at(-1)?.kind).toBe("divider")
    expect(String(view.transcript.items.at(-1)?.text)).toContain("80 columns")
    view.stop()
  })
})

describe("a row that cannot be drawn", () => {
  test("is said once, and the queue still empties", async () => {
    const w = open = await mountWindow({ width: 100 })
    w.view.stop()
    const { startSessionView } = await import("../../src/session")
    // A theme that breaks every row: drawing the note of the failure fails too.
    const view = startSessionView(w.t.renderer, w.store, () => (undefined as never), {
      showThinking: false, replayTurns: 20 })
    void view.transcript.append({ id: "x", kind: "notice", text: "boom" })
    const drained = await Promise.race([view.queue.drain().then(() => "drained"), Bun.sleep(2000).then(() => "stuck")])
    expect(drained).toBe("drained")
    view.stop()
  })
})
