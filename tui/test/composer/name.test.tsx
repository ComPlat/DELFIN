/**
 * One name per session, kept by the agent service and shown alike in every
 * window (TUI, dashboard, phone):
 *
 * - `/name NEW NAME` names the session in view (session.setTitle) and
 *   sends no prompt; `/name` alone puts the default back (the service's
 *   folder name and short id) -- unlike `/rename`, which names it after
 *   its first prompt;
 * - F2 renames the selected row in the agent view, as ctrl+r does (letters
 *   there are typing a new session's prompt, so no plain letter can);
 * - a row is called what the service's summary calls it (`name`), and a
 *   rename to the default seen as `title: ""` is no longer shown.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { sessionName } from "../../src/agents/view"
import { nameOf } from "../../src/commands"
import { Keymap, KEY_HELP } from "../../src/keymap/keymap"
import { applyInfo, newSession } from "../../src/store/state"
import { mountWindow, type Call } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await open?.close(); open = null })

function service() {
  const calls: Call[] = []
  return {
    calls,
    client: {
      request: async (method: string, params: Record<string, any> = {}) => {
        calls.push({ method, params })
        if (method === "session.setTitle") return { title: params.title }
        if (method === "sessions.list") return { sessions: [], elsewhere: [] }
        if (method === "session.replay") return { events: [] }
        if (method === "session.commands") return { commands: [] }
        return { applied: true, message: "", accepted: true }
      },
    },
  }
}

const sent = (calls: Call[], method: string) => calls.filter((c) => c.method === method).map((c) => c.params)

async function line(w: NonNullable<typeof open>, text: string) {
  await w.t.mockInput.typeText(text)
  w.t.mockInput.pressEscape()
  await Bun.sleep(40)
  w.t.mockInput.pressEnter()
  await w.settle()
}

describe("reading a /name line", () => {
  test("takes the name; alone it means the default; nothing else is one", () => {
    expect(nameOf("/name benzene opt")).toBe("benzene opt")
    expect(nameOf("  /name   two   words ")).toBe("two words")
    expect(nameOf("/name")).toBe("")
    expect(nameOf("/names x")).toBeNull()
    expect(nameOf("/rename x")).toBeNull()
    expect(nameOf("call it /name x")).toBeNull()
  })
})

describe("/name in the window", () => {
  test("names the session in view and sends no prompt", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    await line(w, "/name solvent scan")
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "s1", title: "solvent scan" }])
    expect(sent(s.calls, "session.prompt")).toEqual([])
    expect(w.footer().join("\n")).toContain("renamed: solvent scan")
  })

  test("alone, puts the default name back", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    await line(w, "/name")
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "s1", title: "" }])
    expect(sent(s.calls, "session.prompt")).toEqual([])
  })
})

describe("the name a row shows", () => {
  test("is the service's display name before the first prompt", () => {
    const s = newSession("abcdef0123", {})
    s.summary = { sessionId: "abcdef0123", title: null, name: "scan-a · abcdef01", firstPrompt: "hello there" } as never
    expect(sessionName(s)).toBe("scan-a · abcdef01")
  })

  test("a rename to the default clears the title the window held", () => {
    const s = newSession("abcdef0123", { title: "old name" })
    expect(s.title).toBe("old name")
    applyInfo(s, { title: "" })
    expect(s.title).toBe("")
    s.summary = { sessionId: "abcdef0123", title: null, name: "calc · abcdef01" } as never
    expect(sessionName(s)).toBe("calc · abcdef01")
  })
})

describe("F2 in the agent view", () => {
  test("is bound to the rename and listed with it", () => {
    const km = new Keymap()
    expect(km.matches("agentsRename", { name: "f2", ctrl: false, meta: false, shift: false })).toBe(true)
    expect(km.matches("agentsRename", { name: "r", ctrl: true, meta: false, shift: false })).toBe(true)
    // No plain letter: in the view, letters type a new session's prompt.
    expect(km.matches("agentsRename", { name: "r", ctrl: false, meta: false, shift: false })).toBe(false)
    expect(KEY_HELP.some((row) => row.action === "agentsRename")).toBe(true)
  })

  test("renames the selected session, the name ready to edit", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    w.store.push({ type: "session.summary", session: "d1", ts: 0,
                   summary: { sessionId: "d1", title: "write the report", name: "write the report",
                              workspace: "/work/d1", provider: "kit", model: "kit.test", approval: "default",
                              state: "idle", pending: 0, costUsd: 0, createdAt: 1003,
                              lastActivity: Date.now() / 1000, pid: 1 } } as never)
    await w.settle()
    w.keys.pressArrow("left")
    await w.settle()
    w.keys.pressArrow("up")
    await w.settle()
    w.keys.pressKey("f2")
    await w.settle()
    expect(w.footer().join("\n")).toContain("❯ write the report")
    await w.keys.typeText(" v2")
    w.keys.pressEnter()
    await w.settle()
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "d1", title: "write the report v2" }])
    expect(sent(s.calls, "session.prompt")).toEqual([])
  })
})
