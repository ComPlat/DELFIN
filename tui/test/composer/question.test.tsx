/**
 * The prompt keeps the keys while a question is open: what is typed goes
 * into the prompt, and the question takes its keys only from an empty
 * prompt.
 *
 *   a question (ask)  typing + enter answers with the line (permission.reply
 *                     `text`); on an empty prompt 1-9 pick, ↑ ↓ move the
 *                     highlight, enter takes it, esc dismisses
 *   an approval       its keys (y n A k ...) on an empty prompt only, after
 *                     the pause; typing in progress, a paste, enter with a
 *                     line or a ctrl key never answer it -- a line is a
 *                     message (queued for the turn); N takes its reason
 *                     from the prompt
 *
 * After an answer the prompt still has the keys, and any draft it held.
 */
import { TextareaRenderable } from "@opentui/core"
import { afterEach, describe, expect, test } from "bun:test"
import type { Requester } from "../../src/app"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const ask = (id: string, seq: number, multi = false) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "ask", tool: "ask_user_question",
  preview: "", options: [], flags: {}, question: "Which solvent model?", header: "Solvent",
  choices: ["CPCM", "SMD", "none"], choiceDetails: ["conductor-like", "density-based", "gas phase"], multiSelect: multi,
})
const bash = (id: string, seq: number) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "confirm", tool: "bash",
  preview: "$ pytest -q", flags: { shell: true, write: false, protected: false, outsideRead: false },
  suggestion: "pytest *",
  options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "A", label: "always this exact command" },
            { key: "k", label: "always commands like this" }, { key: "a", label: "abort" }, { key: "?", label: "explain" }],
})

const edit = (id: string, seq: number) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "confirm", tool: "edit_file",
  preview: "--- a/x.py\n+++ b/x.py\n-a\n+b", flags: { shell: false, write: true, protected: false, outsideRead: false },
  options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "e", label: "stop asking about writes this session" },
            { key: "a", label: "abort" }, { key: "?", label: "explain" }],
})

async function mount(options: Parameters<typeof mountWindow>[0] = {}) {
  w = await mountWindow(options)
  ;(w.t.renderer as unknown as { copyToClipboardOSC52(text: string): boolean }).copyToClipboardOSC52 = () => true
  w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, sessionId: "s1", workspace: "/w",
                 provider: "fake", model: "m", mode: "solo", approval: "default", banner: [], notices: [] } as never)
  // A turn runs: the question belongs to it.
  w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
  await w.settle()
  return w
}
const asked = async (event: object) => {
  w!.store.push(event as never)
  await w!.settle()
}
async function typeAt(text: string, gapMs: number) {
  for (const ch of text) {
    w!.keys.pressKey(ch === " " ? " " : ch, ch !== ch.toLowerCase() ? { shift: true } : undefined)
    await Bun.sleep(gapMs)
  }
  await w!.settle()
}
const type = async (text: string) => {
  await w!.keys.typeText(text)
  await w!.settle()
}
/** The test renderer's names for the keys that are no character. */
const NAMED: Record<string, string> = { return: "RETURN", up: "ARROW_UP", down: "ARROW_DOWN", backspace: "BACKSPACE" }
const press = async (key: string, mods?: Record<string, boolean>) => {
  w!.keys.pressKey((NAMED[key] ?? key) as never, mods as never)
  await Bun.sleep(5)
  await w!.settle()
}
/** The prompt's line, as drawn. */
const promptLine = () => w!.footer().find((l) => l.startsWith("❯ ") || l.startsWith("! ")) ?? ""
const promptFocused = () => w!.t.renderer.currentFocusedRenderable instanceof TextareaRenderable
const prompts = () => w!.calls.filter((c) => c.method === "session.prompt").map((c) => c.params)

describe("a question keeps the prompt", () => {
  test("the prompt has the keys while it is open: what is typed is in the prompt", async () => {
    await mount()
    await asked(ask("q1", 2))
    expect(promptFocused()).toBe(true)
    await type("toluene")
    expect(promptLine()).toContain("❯ toluene")
    expect(w!.replies()).toEqual([])
  })

  test("enter with a line answers with the line; the prompt is empty and keeps the keys", async () => {
    await mount()
    await asked(ask("q1", 2))
    await type("toluene, implicit")
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", text: "toluene, implicit" }])
    expect(prompts()).toEqual([])
    expect(promptLine()).not.toContain("toluene")
    expect(promptFocused()).toBe(true)
  })

  test("digits after a letter are the line's: 'option 2' is a typed answer", async () => {
    await mount()
    await asked(ask("q1", 2))
    await type("option 2")
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", text: "option 2" }])
  })

  test("a number on an empty prompt picks that choice; the prompt keeps the keys", async () => {
    await mount()
    await asked(ask("q1", 2))
    await press("2")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", answers: ["SMD"] }])
    expect(promptLine()).not.toContain("2")
    expect(promptFocused()).toBe(true)
  })

  test("↓ moves the highlight and enter takes it", async () => {
    await mount()
    await asked(ask("q1", 2))
    expect(w!.footer().join("\n")).toContain("› 1. CPCM")
    await press("down")
    expect(w!.footer().join("\n")).toContain("› 2. SMD")
    await press("down")
    await press("up")
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", answers: ["SMD"] }])
  })

  test("several answers: space toggles the highlighted one, enter sends", async () => {
    await mount()
    await asked(ask("q1", 2, true))
    await press(" ")                                      // CPCM on
    await press("down")
    await press("down")
    await press(" ")                                      // none on
    await press("3")                                      // none off again
    expect(w!.footer().join("\n")).toContain("[x] CPCM")
    expect(w!.replies()).toEqual([])
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", answers: ["CPCM"] }])
  })

  test("an answer typed into an empty prompt that begins with a digit is the line's", async () => {
    await mount({ answerPauseMs: 350 })
    await asked(ask("q1", 2))
    await Bun.sleep(1200)
    await typeAt("2 of them", 70)
    await Bun.sleep(500)
    expect(w!.replies()).toEqual([])
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", text: "2 of them" }])
  })

  test("esc on an empty prompt dismisses it", async () => {
    await mount()
    await asked(ask("q1", 2))
    await w!.escape()
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", answers: [] }])
  })

  test("esc with a line typed neither dismisses nor ends the turn; the line stays", async () => {
    await mount()
    await asked(ask("q1", 2))
    await type("tol")
    await w!.escape()
    expect(w!.replies()).toEqual([])
    expect(w!.calls.filter((c) => c.method === "session.cancel")).toEqual([])
    expect(promptLine()).toContain("❯ tol")
  })

  test("a question that opens while the user types: the typing goes on, and answers with the line", async () => {
    await mount({ answerPauseMs: 350 })
    await typeAt("the ", 60)
    await asked(ask("q1", 2))
    await typeAt("gas phase please", 60)
    await Bun.sleep(400)
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", text: "the gas phase please" }])
  })

  test("enter at once as the question opens is not taken: the line stays in the prompt", async () => {
    await mount({ answerPauseMs: 350 })
    await type("run the tests")
    await asked(ask("q1", 2))
    await press("return")
    expect(w!.replies()).toEqual([])
    expect(prompts()).toEqual([])
    expect(promptLine()).toContain("❯ run the tests")
  })

  test("a click on a choice answers, and the draft stays in the prompt with the keys", async () => {
    await mount({ full: true })
    await type("keep this draft")
    await asked(ask("q1", 2))
    const rows = w!.t.captureCharFrame().split("\n")
    const y = rows.findIndex((r) => r.includes("density-based"))
    await w!.t.mockMouse.click(rows[y]!.indexOf("density-based") + 1, y)
    await Bun.sleep(20)
    await w!.settle()
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", answers: ["SMD"] }])
    expect(promptLine()).toContain("❯ keep this draft")
    expect(promptFocused()).toBe(true)
    await type(" and more")
    expect(promptLine()).toContain("❯ keep this draft and more")
  })

  test("a service that has no typed answer dismissed it: the line goes out as a message", async () => {
    const calls: { method: string; params: any }[] = []
    const client: Requester = {
      request: async (method: string, params: Record<string, any> = {}) => {
        calls.push({ method, params })
        return method === "permission.reply" ? { applied: true, message: "dismissed" } : { accepted: true }
      },
    }
    await mount({ client })
    await asked(ask("q1", 2))
    await type("toluene")
    await press("return")
    await Bun.sleep(20)
    await w!.settle()
    expect(calls.filter((c) => c.method === "permission.reply").map((c) => c.params))
      .toEqual([{ sessionId: "s1", requestId: "q1", text: "toluene" }])
    expect(calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)).toEqual(["toluene"])
  })
})

describe("an approval keeps the prompt, and typing never answers it", () => {
  test("typing in progress as it opens goes on into the prompt: no answer", async () => {
    await mount({ answerPauseMs: 350 })
    await typeAt("please re", 60)
    await asked(bash("c1", 2))
    expect(promptFocused()).toBe(true)
    await typeAt("yes and abort nothing", 60)
    await Bun.sleep(1200)
    await typeAt("y", 60)
    expect(w!.replies()).toEqual([])
    expect(promptLine()).toContain("please reyes and abort nothingy")
  })

  test("enter with a line is a message, never an answer", async () => {
    await mount()
    await asked(bash("c1", 2))
    await type("hold on, use -x")
    await press("return")
    expect(w!.replies()).toEqual([])
    expect(prompts().map((p) => p.text)).toEqual(["hold on, use -x"])
  })

  test("its key on an empty prompt answers, after the pause (unchanged)", async () => {
    await mount({ answerPauseMs: 350 })
    await asked(bash("c1", 2))
    await press("y")                                      // at once: a reflex
    expect(w!.replies()).toEqual([])
    await press("backspace")
    await Bun.sleep(1100)
    await press("y")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "y" }])
    expect(promptFocused()).toBe(true)
  })

  test("letters typed, then deleted, then y at once: typing, not an answer", async () => {
    await mount({ answerPauseMs: 350 })
    await asked(bash("c1", 2))
    await Bun.sleep(1200)                                 // open a while: no reflex left
    await type("hello")
    for (let i = 0; i < 5; i++) await press("backspace")
    await press("y")
    expect(w!.replies()).toEqual([])
    expect(promptLine()).toContain("❯ y")
  })

  test("a paste is never an answer", async () => {
    await mount()
    await asked(bash("c1", 2))
    await w!.keys.pasteBracketedText("y")
    await w!.settle()
    expect(w!.replies()).toEqual([])
    expect(promptLine()).toContain("❯ y")
  })

  test("ctrl and alt letters edit the prompt, never answer", async () => {
    await mount()
    await asked(bash("c1", 2))
    for (const [key, mods] of [["y", { ctrl: true }], ["a", { ctrl: true }], ["e", { ctrl: true }],
                               ["n", { ctrl: true }], ["y", { meta: true }]] as const) {
      await press(key, mods)
    }
    expect(w!.replies()).toEqual([])
  })

  test("N: the reason is typed in the prompt, enter refuses with it", async () => {
    await mount()
    await asked(bash("c1", 2))
    await press("N", { shift: true })
    expect(w!.replies()).toEqual([])
    expect(w!.footer().join("\n")).toContain("why?")
    await type("use pytest -x")
    expect(promptLine()).toContain("❯ use pytest -x")
    await press("return")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "n", reason: "use pytest -x" }])
    expect(prompts()).toEqual([])
  })

  test("A still takes its second step; typing in the prompt does not take it", async () => {
    await mount()
    await asked(bash("c1", 2))
    await press("A", { shift: true })
    expect(w!.replies()).toEqual([])
    expect(w!.footer().join("\n")).toContain("persist")
    await type("ok")
    expect(w!.replies()).toEqual([])
    for (let i = 0; i < 2; i++) await press("backspace")
    await press("y")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "A", persist: true }])
  })

  // A message typed into an empty prompt may begin with an option key
  // ("yes, but", "every", "Run"): the key answers only alone -- another
  // key soon after it makes both typing.
  for (const [words, question] of [["yes please", bash], ["every file is fine", edit], ["keep going", bash],
                                    ["abort? no", bash]] as const) {
    test(`a message typed into an empty prompt ("${words}") is typing, not an answer`, async () => {
      await mount({ answerPauseMs: 350 })
      await asked(question("c1", 2))
      await Bun.sleep(1200)                               // open a while: a deliberate first key
      await typeAt(words, 70)
      await Bun.sleep(500)
      expect(w!.replies()).toEqual([])
      expect(promptLine()).toContain(`❯ ${words}`)
    })
  }

  test("a lone key answers once nothing follows it; with enter at once, at once", async () => {
    await mount({ answerPauseMs: 350 })
    await asked(bash("c1", 2))
    await Bun.sleep(1200)
    await press("y")
    expect(w!.replies()).toEqual([])                      // held: a word may follow
    await Bun.sleep(450)
    await w!.settle()
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "y" }])
    expect(promptLine()).not.toContain("y")
    await asked(bash("c2", 3))
    await Bun.sleep(1200)
    await press("n")
    await press("return")
    expect(w!.replies().at(-1)).toEqual({ sessionId: "s1", requestId: "c2", choice: "n" })
  })

  test("the second step of A: 'yes' typed is typing, a lone y saves", async () => {
    await mount({ answerPauseMs: 350 })
    await asked(bash("c1", 2))
    await Bun.sleep(1200)
    await press("A", { shift: true })
    await Bun.sleep(1200)
    await typeAt("yes", 70)
    await Bun.sleep(500)
    expect(w!.replies()).toEqual([])
    for (let i = 0; i < 3; i++) await press("backspace")
    await Bun.sleep(1200)
    await press("y")
    await Bun.sleep(450)
    await w!.settle()
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "c1", choice: "A", persist: true }])
  })

  test("with a draft in the prompt, its keys are the draft's", async () => {
    await mount()
    await type("note: ")
    await asked(bash("c1", 2))
    await type("y")
    await w!.escape()
    expect(w!.replies()).toEqual([])
    expect(w!.calls.filter((c) => c.method === "session.cancel")).toEqual([])
    expect(promptLine()).toContain("❯ note: y")
  })
})
