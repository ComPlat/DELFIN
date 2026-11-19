/**
 * Entering a session the way `claude --resume` does: the resume picker
 * (a bare `delfin-agent tui --resume`, `/resume`, the palette) lists the
 * running sessions on top, then the saved ones, each with title, folder,
 * age and state; typing filters by title, folder or id; ↑ ↓ select, enter
 * opens, esc closes. At the start, esc leaves the window to a new session.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { createPicker } from "../../src/footer/picker"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const key = (name: string, extra: Record<string, unknown> = {}) =>
  ({ name, ctrl: false, meta: false, shift: false, option: false, sequence: name.length === 1 ? name : "", ...extra }) as never

describe("a picker's filter and its cancel", () => {
  test("the filter reads a choice's hidden search text too", () => {
    const picker = createPicker()
    picker.open("resume a session", [
      { label: "the parser fix", detail: "saved · ~/a", search: "0f3a9c21-77aa /home/u/a", value: "a" },
      { label: "benzene scan", detail: "running · ~/calc", search: "b81d0000-1111 /home/u/calc/benzene", value: "b" },
    ], () => {})
    for (const ch of "77aa") picker.handleKey(key(ch))
    expect(picker.shown().map((c) => c.value)).toEqual(["a"])
    for (let i = 0; i < 4; i++) picker.handleKey(key("backspace"))
    for (const ch of "calc/benz") picker.handleKey(key(ch))
    expect(picker.shown().map((c) => c.value)).toEqual(["b"])
  })

  test("esc cancels; a pick does not; a list opened in its place cancels it", () => {
    const picker = createPicker()
    const said: string[] = []
    picker.open("one", [{ label: "x", value: "x" }], (c) => said.push(`pick ${c.value}`), undefined,
                () => said.push("cancel one"))
    picker.handleKey(key("escape"))
    picker.open("two", [{ label: "y", value: "y" }], (c) => said.push(`pick ${c.value}`), undefined,
                () => said.push("cancel two"))
    picker.handleKey(key("return"))
    picker.open("three", [{ label: "z", value: "z" }], () => {}, undefined, () => said.push("cancel three"))
    picker.open("four", [{ label: "w", value: "w" }], () => {})
    picker.close()
    expect(said).toEqual(["cancel one", "pick y", "cancel three"])
  })
})

const NOW = Date.now() / 1000
function client(calls: { method: string; params: any }[]) {
  return {
    request: async (method: string, params: any = {}) => {
      calls.push({ method, params })
      if (method === "sessions.list") {
        return { sessions: [
          { sessionId: "run-1111-aaaa", title: "benzene scan", workspace: "/home/u/calc/benzene", state: "running",
            lastActivity: NOW - 20 },
        ] }
      }
      if (method === "sessions.past") {
        return { sessions: [
          { sessionId: "old-2222-bbbb", title: "the parser fix", workspace: "/w", model: "kit.x",
            updatedAt: NOW - 7200, messages: 12, running: false },
          { sessionId: "old-3333-cccc", title: "docs pass", workspace: "/w", model: "kit.x",
            updatedAt: NOW - 600, messages: 3, running: false },
        ] }
      }
      return { accepted: true }
    },
  }
}

describe("delfin-agent tui --resume (no id): the picker before any session", () => {
  test("it lists running then saved sessions, newest first, with state, folder and age", async () => {
    const calls: { method: string; params: any }[] = []
    const w = open = await mountWindow({ client: client(calls), height: 40, show: false,
                                        resumeAtStart: { workspace: "/w", onPick: () => {}, onCancel: () => {} } })
    await w.settle()
    const footer = w.footer().join("\n")
    expect(footer).toContain("resume a session")
    const order = ["benzene scan", "docs pass", "the parser fix"].map((t) => footer.indexOf(t))
    expect(order.every((at) => at >= 0)).toBe(true)
    expect([...order].sort((a, b) => a - b)).toEqual(order)
    expect(footer).toContain("running")
    expect(footer).toContain("saved")
    expect(footer).toContain("2h ago")
    expect(calls.find((c) => c.method === "sessions.past")?.params.workspace).toBe("/w")
  })

  test("typing a folder or an id narrows it; enter opens the one selected", async () => {
    const picked: string[] = []
    const w = open = await mountWindow({ client: client([]), height: 40, show: false,
                                        resumeAtStart: { workspace: "/w", onPick: (id, ws) => picked.push(`${id} ${ws}`),
                                                         onCancel: () => picked.push("cancel") } })
    await w.settle()
    await w.keys.typeText("3333-cc")
    await w.settle()
    let footer = w.footer().join("\n")
    expect(footer).toContain("docs pass")
    expect(footer).not.toContain("the parser fix")
    for (let i = 0; i < 7; i++) w.keys.pressBackspace()
    await w.keys.typeText("calc/benz")
    await w.settle()
    footer = w.footer().join("\n")
    expect(footer).toContain("benzene scan")
    expect(footer).not.toContain("docs pass")
    w.keys.pressEnter()
    await w.settle()
    expect(picked).toEqual(["run-1111-aaaa /home/u/calc/benzene"])
  })

  test("esc closes it and hands the start back (a new session)", async () => {
    const picked: string[] = []
    const w = open = await mountWindow({ client: client([]), height: 40, show: false,
                                        resumeAtStart: { workspace: "/w", onPick: (id) => picked.push(id),
                                                         onCancel: () => picked.push("cancel") } })
    await w.settle()
    await w.escape()
    expect(picked).toEqual(["cancel"])
    expect(w.footer().join("\n")).not.toContain("resume a session")
  })

  test("arrows then enter open the second row", async () => {
    const picked: string[] = []
    const w = open = await mountWindow({ client: client([]), height: 40, show: false,
                                        resumeAtStart: { workspace: "/w", onPick: (id) => picked.push(id),
                                                         onCancel: () => picked.push("cancel") } })
    await w.settle()
    w.keys.pressArrow("down")
    await w.settle()
    w.keys.pressEnter()
    await w.settle()
    expect(picked).toEqual(["old-3333-cccc"])
  })
})

describe("/resume in a session", () => {
  test("the same rows, filtered by the whole id", async () => {
    const resumed: string[] = []
    const w = open = await mountWindow({ client: client([]), height: 40, onResume: (id) => resumed.push(id) })
    await w.keys.typeText("/resume")
    w.keys.pressEscape()                                // the completion list
    await Bun.sleep(40)
    w.keys.pressEnter()
    await w.settle()
    expect(w.footer().join("\n")).toContain("benzene scan")
    await w.keys.typeText("old-2222-bbbb")
    await w.settle()
    const footer = w.footer().join("\n")
    expect(footer).toContain("the parser fix")
    expect(footer).not.toContain("benzene scan")
    w.keys.pressEnter()
    await w.settle()
    expect(resumed).toEqual(["old-2222-bbbb"])
  })
})
