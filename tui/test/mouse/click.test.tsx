/**
 * Clicks answer what keys answer -- the approval panel's options, a
 * question's choices, a picker's rows -- and nothing else does: not a
 * click that comes within the pause after the panel appeared or changed,
 * not the wheel, not a drag, not the right button, not a click beside an
 * option. Saving a rule still takes its second step.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const bash = (id: string, seq: number) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "confirm", tool: "bash",
  preview: "$ pytest -q", flags: { shell: true, write: false, protected: false, outsideRead: false },
  suggestion: "pytest *",
  options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "A", label: "always this exact command" },
            { key: "k", label: "always commands like this" }, { key: "a", label: "abort" }, { key: "?", label: "explain" }],
})
const ask = (id: string, seq: number, multi = false) => ({
  type: "permission.asked", session: "s1", ts: 0, seq, requestId: id, kind: "ask", tool: "ask_user_question",
  preview: "", options: [], flags: {}, question: "Which solvent model?", header: "Solvent",
  choices: ["CPCM", "SMD", "none"], choiceDetails: ["conductor-like", "density-based", "gas phase"], multiSelect: multi,
})

async function mount(options: Parameters<typeof mountWindow>[0] = {}) {
  w = await mountWindow({ full: true, ...options })
  // A drag selects and copies; the test terminal has no clipboard.
  ;(w.t.renderer as unknown as { copyToClipboardOSC52(text: string): boolean }).copyToClipboardOSC52 = () => true
  w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, sessionId: "s1", workspace: "/w",
                 provider: "fake", model: "m", mode: "solo", approval: "default", banner: [], notices: [] } as never)
  await w.settle()
  return w
}
const asked = async (event: object) => {
  w!.store.push(event as never)
  await w!.settle()
}
/** Where `text` is on the screen (its first cell). */
function where(text: string): { x: number; y: number } {
  const rows = w!.t.captureCharFrame().split("\n")
  const y = rows.findIndex((r) => r.includes(text))
  if (y < 0) throw new Error(`${text} is not on the screen:\n${rows.join("\n")}`)
  return { x: rows[y]!.indexOf(text), y }
}
const settle = async () => {
  await Bun.sleep(20)
  await w!.settle()
}
const click = async (text: string, dx = 1, button: 0 | 2 = 0) => {
  const at = where(text)
  await w!.t.mockMouse.click(at.x + dx, at.y, button)
  await settle()
}

describe("a click on an approval", () => {
  test("on an option answers it as its key does", async () => {
    await mount()
    await asked(bash("q1", 1))
    await click("[y] yes")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", choice: "y" }])
  })

  test("beside the options, on the preview or with another button: nothing", async () => {
    await mount()
    await asked(bash("q1", 1))
    await click("pytest -q")
    await click("[n] no", 1, 2)
    const at = where("[?] explain")
    await w!.t.mockMouse.click(at.x + "[?] explain".length + 1, at.y)   // the spaces after the last option
    await settle()
    expect(w!.replies()).toEqual([])
  })

  test("the wheel and a drag over an option answer nothing", async () => {
    await mount()
    await asked(bash("q1", 1))
    const at = where("[y] yes")
    await w!.t.mockMouse.scroll(at.x + 1, at.y, "down")
    await w!.t.mockMouse.drag(at.x, at.y, at.x + 4, at.y)
    await settle()
    expect(w!.replies()).toEqual([])
  })

  test("a click right after the question appeared is not an answer; one after the pause is", async () => {
    await mount({ clickPauseMs: 300 })
    await asked(bash("q1", 1))
    await click("[y] yes")
    expect(w!.replies()).toEqual([])
    await Bun.sleep(320)
    await click("[y] yes")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", choice: "y" }])
  })

  test("always-this-command still asks again before a rule is saved, and that step has its own pause", async () => {
    await mount({ clickPauseMs: 300 })
    await asked(bash("q1", 1))
    await Bun.sleep(320)
    await click("[k] always commands like this")
    expect(w!.replies()).toEqual([])
    expect(where("[y] persist").y).toBeGreaterThan(0)
    // The second click of a double-click lands on the new step: not taken.
    await click("[y] persist")
    expect(w!.replies()).toEqual([])
    await Bun.sleep(320)
    await click("[y] persist")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", choice: "k", persist: true }])
  })

  test("just this once, from the second step", async () => {
    await mount()
    await asked(bash("q1", 1))
    await click("[A] always this exact command")
    await click("[n] just this once")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q1", choice: "A", persist: false }])
  })

  test("[N] asks for the reason, which is typed", async () => {
    await mount()
    await asked(bash("q1", 1))
    await click("[N] no + why")
    expect(w!.replies()).toEqual([])
    expect(where("why?").y).toBeGreaterThan(0)
  })
})

describe("a click on a question", () => {
  test("on a choice's row answers it", async () => {
    await mount()
    await asked(ask("q2", 1))
    await click("density-based")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q2", answers: ["SMD"] }])
  })

  test("several: a click toggles, [enter] send sends", async () => {
    await mount()
    await asked(ask("q3", 1, true))
    await click("CPCM")
    await click("gas phase")
    expect(w!.replies()).toEqual([])
    expect(w!.t.captureCharFrame()).toContain("[x] CPCM")
    await click("[enter] send")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q3", answers: ["CPCM", "none"] }])
  })

  test("[esc] dismiss dismisses", async () => {
    await mount()
    await asked(ask("q2", 1))
    await click("[esc] dismiss")
    expect(w!.replies()).toEqual([{ sessionId: "s1", requestId: "q2", answers: [] }])
  })
})

describe("a click on a picker", () => {
  const client = (calls: { method: string; params: any }[]) => ({
    request: async (method: string, params: Record<string, any> = {}) => {
      calls.push({ method, params })
      if (method === "session.choices") {
        return { current: "small-model", choices: [{ value: "small-model", detail: "fast" },
                                                   { value: "large-model", detail: "careful" }] }
      }
      return { accepted: true }
    },
  })

  test("/model: a click on a row picks that model", async () => {
    const calls: { method: string; params: any }[] = []
    await mount({ client: client(calls) })
    w!.keys.typeText("/model")
    w!.keys.pressEnter()
    await settle()
    expect(where("large-model").y).toBeGreaterThan(0)
    await click("large-model")
    expect(calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)).toEqual(["/model large-model"])
  })

  test("a click as the list appears picks nothing", async () => {
    const calls: { method: string; params: any }[] = []
    await mount({ client: client(calls), clickPauseMs: 300 })
    w!.keys.typeText("/model")
    w!.keys.pressEnter()
    await settle()
    await click("large-model")
    expect(calls.filter((c) => c.method === "session.prompt")).toEqual([])
    await Bun.sleep(320)
    await click("large-model")
    expect(calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)).toEqual(["/model large-model"])
  })
})

describe("where the window takes the mouse", () => {
  test("full screen only, unless tui.json says no or the terminal reports none", async () => {
    const { mouseOn, parseConfig } = await import("../../src/config/config")
    expect(mouseOn(true, "full", { TERM: "xterm-256color" })).toBe(true)
    expect(mouseOn(true, "inline", { TERM: "xterm-256color" })).toBe(false)
    expect(mouseOn(false, "full", { TERM: "xterm-256color" })).toBe(false)
    expect(mouseOn(true, "full", { TERM: "dumb" })).toBe(false)
    expect(parseConfig('{"mouse": false}').config.mouse).toBe(false)
    expect(parseConfig('{"mouse": "no"}').problems).toEqual(["mouse must be true or false"])
  })
})
