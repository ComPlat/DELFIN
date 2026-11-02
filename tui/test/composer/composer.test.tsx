/**
 * The composer in the window, driven with keys the way a person types:
 * sending, new lines, the history, ctrl+r, completion, large pastes and
 * the ctrl+c / ctrl+d ladder. Rendered headless through the real App; the
 * agent service is a fake that records what it is asked.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { createTestRenderer } from "@opentui/core/testing"
import { render } from "@opentui/solid"
import { App } from "../../src/app"
import { History, type Entry } from "../../src/composer/history"
import { Keymap } from "../../src/keymap/keymap"
import { ensureSession } from "../../src/store/state"
import { createAppStore } from "../../src/store/store"
import { DARK } from "../../src/theme/theme"

type Call = { method: string; params: Record<string, any> }
type Setup = Awaited<ReturnType<typeof createTestRenderer>>

let current: Setup | null = null
afterEach(async () => {
  await current?.renderer.destroy()
  current = null
})

const COMMANDS = [
  { name: "/help", summary: "List these commands", category: "session", takesArgs: true, source: "builtin" },
  { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" },
  { name: "/session", summary: "List or search past sessions", category: "history", takesArgs: true, source: "builtin" },
]

async function mount(options: { history?: Entry[]; turn?: boolean } = {}) {
  const setup = await createTestRenderer({
    width: 80, height: 24, screenMode: "split-footer", footerHeight: 4,
    externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false,
  })
  current = setup
  const store = createAppStore()
  store.update((d) => {
    const s = ensureSession(d, "s1", { model: "kit.test", approval: "default" })
    s.state = "idle"
    if (options.turn) s.turn = { id: "t1", source: "user", startedAt: Date.now(), label: "thinking…",
                                 base: { inputTokens: 0, outputTokens: 0, costUsd: 0 }, streamed: 0 }
    d.current = "s1"
  })
  const calls: Call[] = []
  const client = {
    request: async (method: string, params: Record<string, any> = {}) => {
      calls.push({ method, params })
      if (method === "session.commands") return { commands: COMMANDS }
      if (method === "session.completePath") {
        const all = ["delfin/", "delfin/xtb.py", "docs/"]
        const frag = String(params.fragment)
        return { paths: all.filter((p) => p.startsWith(frag) && p !== frag) }
      }
      return { accepted: true }
    },
  }
  const history = new History(null, options.history ?? [])
  let quits = 0
  await render(() => (
    <App renderer={setup.renderer} store={store} keymap={new Keymap()} theme={() => DARK}
         client={client} history={history} onQuit={() => quits++} onRedraw={() => {}} />
  ), setup.renderer)
  await setup.renderOnce()
  const settle = async () => {
    for (let i = 0; i < 3; i++) {
      await Bun.sleep(5)
      await setup.renderOnce()
    }
  }
  // A lone Esc is held 20 ms by the input parser (it may start a sequence).
  const escape = async () => {
    setup.mockInput.pressEscape()
    await Bun.sleep(40)
    await settle()
  }
  const prompts = () => calls.filter((c) => c.method === "session.prompt").map((c) => c.params)
  const frame = () => setup.captureCharFrame()
  return { setup, keys: setup.mockInput, calls, prompts, history, frame, settle, escape, quits: () => quits }
}

describe("the composer", () => {
  test("sends what is typed with enter, and remembers it", async () => {
    const w = await mount()
    await w.keys.typeText("hello there")
    await w.settle()
    expect(w.frame()).toContain("❯ hello there")
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()).toEqual([{ sessionId: "s1", text: "hello there" }])
    expect(w.frame()).toContain("ask DELFIN")             // empty again
    expect(w.history.entries.map((e) => e.text)).toEqual(["hello there"])
  })

  test("alt+enter, ctrl+j and a trailing backslash start a new line", async () => {
    const w = await mount()
    await w.keys.typeText("one")
    w.keys.pressEnter({ meta: true })
    await w.keys.typeText("two")
    w.keys.pressKey("\n")                                 // ctrl+j without kitty keys
    await w.keys.typeText("three\\")
    w.keys.pressEnter()
    await w.keys.typeText("four")
    await w.settle()
    expect(w.prompts()).toEqual([])                       // nothing sent yet
    const shown = w.frame()
    for (const part of ["❯ one", "two", "three", "four"]) expect(shown).toContain(part)
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()[0]!.text).toBe("one\ntwo\nthree\nfour")
  })

  test("↑ and ↓ walk the history and come back to the draft", async () => {
    const w = await mount({ history: [{ text: "first" }, { text: "second" }] })
    await w.keys.typeText("draft")
    w.keys.pressArrow("up")
    await w.settle()
    expect(w.frame()).toContain("❯ second")
    w.keys.pressArrow("up")
    await w.settle()
    expect(w.frame()).toContain("❯ first")
    w.keys.pressArrow("down")
    w.keys.pressArrow("down")
    await w.settle()
    expect(w.frame()).toContain("❯ draft")
    expect(w.prompts()).toEqual([])
  })

  test("ctrl+r finds an older prompt; enter takes it into the prompt, not out", async () => {
    const w = await mount({ history: [{ text: "run orca opt" }, { text: "ls" }, { text: "orca freq job" }] })
    w.keys.pressKey("r", { ctrl: true })
    await w.keys.typeText("orca")
    await w.settle()
    expect(w.frame()).toContain("search history: orca")
    expect(w.frame()).toContain("❯ orca freq job")
    w.keys.pressKey("r", { ctrl: true })                 // the next older match
    await w.settle()
    expect(w.frame()).toContain("❯ run orca opt")
    w.keys.pressEnter()
    await w.settle()
    expect(w.frame()).not.toContain("search history")
    expect(w.frame()).toContain("❯ run orca opt")
    expect(w.prompts()).toEqual([])
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()[0]!.text).toBe("run orca opt")
  })

  test("esc leaves the search and puts the draft back", async () => {
    const w = await mount({ history: [{ text: "orca freq job" }] })
    await w.keys.typeText("my draft")
    w.keys.pressKey("r", { ctrl: true })
    await w.keys.typeText("orca")
    await w.escape()
    expect(w.frame()).toContain("❯ my draft")
  })

  test("a large paste stands in as a placeholder and goes out whole", async () => {
    const w = await mount()
    const log = Array.from({ length: 30 }, (_, i) => `line ${i} of the log`).join("\n")
    await w.keys.typeText("why does this fail? ")
    await w.keys.pasteBracketedText(log)
    await w.settle()
    expect(w.frame()).toContain("[Pasted text #1 +30 lines]")
    expect(w.frame()).not.toContain("line 29 of the log")
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()[0]!.text).toBe(`why does this fail? ${log}`)
    // The history keeps the placeholder, and what it stands for.
    expect(w.history.entries[0]!.text).toBe("why does this fail? [Pasted text #1 +30 lines]")
    expect(Object.values(w.history.entries[0]!.pastes ?? {})).toEqual([log])
  })

  test("a short paste goes in as it is", async () => {
    const w = await mount()
    await w.keys.pasteBracketedText("two\nlines")
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()[0]!.text).toBe("two\nlines")
  })
})

describe("completion in the composer", () => {
  test("a /command list opens as it is typed; tab takes one; enter runs one without arguments", async () => {
    const w = await mount()
    await w.keys.typeText("/s")
    await w.settle()
    expect(w.frame()).toContain("/status")
    expect(w.frame()).toContain("Model, mode, tokens, cost")
    expect(w.frame()).toContain("/session")
    // SW2b added /swarm, a window command that starts with "s": window
    // commands list before the session's (composer.tsx loadCommands),
    // so the first /s match is /swarm and the second down is /session.
    w.keys.pressArrow("down")
    w.keys.pressArrow("down")
    w.keys.pressTab()
    await w.settle()
    expect(w.frame()).toContain("❯ /session")
    w.keys.pressKey("u", { ctrl: true })                 // clear the line
    await w.keys.typeText("/sta")
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()).toEqual([{ sessionId: "s1", text: "/status" }])
  })

  test("an @path completes from the session's workspace, folder by folder", async () => {
    const w = await mount()
    await w.keys.typeText("look at @del")
    await w.settle()
    expect(w.calls.some((c) => c.method === "session.completePath" && c.params.fragment === "del")).toBe(true)
    expect(w.frame()).toContain("@delfin/")
    w.keys.pressTab()
    await w.settle()
    expect(w.frame()).toContain("❯ look at @delfin/")
    expect(w.frame()).toContain("@delfin/xtb.py")         // the list goes on inside
    w.keys.pressTab()
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()[0]!.text).toBe("look at @delfin/xtb.py")
  })

  test("esc closes the list; enter then sends what is typed", async () => {
    const w = await mount()
    await w.keys.typeText("/s")
    await w.settle()
    await w.escape()
    expect(w.frame()).not.toContain("Model, mode, tokens, cost")
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()[0]!.text).toBe("/s")
  })
})

describe("the keys around the composer", () => {
  test("ctrl+c clears what is typed, then twice leaves", async () => {
    const w = await mount()
    await w.keys.typeText("half a thought")
    w.keys.pressCtrlC()
    await w.settle()
    expect(w.frame()).toContain("ask DELFIN")
    expect(w.quits()).toBe(0)
    w.keys.pressCtrlC()
    await w.settle()
    expect(w.frame()).toContain("ctrl+c again to leave")
    w.keys.pressCtrlC()
    expect(w.quits()).toBe(1)
  })

  test("ctrl+d leaves on an empty prompt only", async () => {
    const w = await mount()
    await w.keys.typeText("ab")
    w.keys.pressArrow("left")
    w.keys.pressKey("d", { ctrl: true })                 // deletes the b
    await w.settle()
    expect(w.quits()).toBe(0)
    expect(w.frame()).toContain("❯ a")
    w.keys.pressKey("u", { ctrl: true })
    w.keys.pressKey("d", { ctrl: true })
    expect(w.quits()).toBe(1)
  })

  test("during a turn, ctrl+g steers what is typed into it; enter queues it", async () => {
    const w = await mount({ turn: true })
    await w.keys.typeText("use the def2 basis")
    w.keys.pressKey("g", { ctrl: true })
    await w.settle()
    await w.keys.typeText("then plot it")
    w.keys.pressEnter()
    await w.settle()
    expect(w.prompts()).toEqual([
      { sessionId: "s1", text: "use the def2 basis", how: "steer" },
      { sessionId: "s1", text: "then plot it" },
    ])
  })
})
