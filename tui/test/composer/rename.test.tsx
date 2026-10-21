/**
 * /rename names the session in view: the window asks the service
 * (session.setTitle) at once -- also while a turn runs -- and never sends
 * the line as a prompt. Without a name it is named after its first
 * prompt's opening words. The name stands on the prompt box's top border,
 * at the right (not in the agent view, whose box is for another session),
 * and in the terminal's title; completion and the palette offer it.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { nameFromPrompt, renameOf } from "../../src/commands"
import { mountWindow, type Call } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await open?.close(); open = null })

/** A service that names a session as asked, and records what it is asked. */
function service() {
  const calls: Call[] = []
  return {
    calls,
    client: {
      request: async (method: string, params: Record<string, any> = {}) => {
        calls.push({ method, params })
        if (method === "session.setTitle") return { title: params.title }
        if (method === "session.commands") return { commands: [{ name: "/help", summary: "help", category: "x", takesArgs: false, source: "builtin" }] }
        return { applied: true, message: "" }
      },
    },
  }
}

const sent = (calls: Call[], method: string) => calls.filter((c) => c.method === method).map((c) => c.params)

async function line(w: NonNullable<typeof open>, text: string) {
  await w.t.mockInput.typeText(text)
  w.t.mockInput.pressEscape()                    // the completion list, if any
  await Bun.sleep(40)
  w.t.mockInput.pressEnter()
  await w.settle()
}

function asked(w: NonNullable<typeof open>, firstPrompt: string) {
  w.store.update((d) => { d.sessions.s1!.summary = { sessionId: "s1", state: "idle", firstPrompt } as never })
}

describe("reading a /rename line", () => {
  test("takes the name, or none, and nothing else", () => {
    expect(renameOf("/rename Parser fix")).toBe("Parser fix")
    expect(renameOf("  /rename   two   words ")).toBe("two words")
    expect(renameOf("/rename")).toBe("")
    expect(renameOf("/renamed x")).toBeNull()
    expect(renameOf("please /rename it")).toBeNull()
  })

  test("names after a first prompt by its opening words", () => {
    expect(nameFromPrompt("Read notes.txt and tell me in one sentence what it says.")).toBe("Read notes.txt and tell me")
    expect(nameFromPrompt("hello there.")).toBe("hello there")
    expect(nameFromPrompt("optimise:")).toBe("optimise")
    expect(nameFromPrompt("a".repeat(60))).toBe("a".repeat(40))
    expect(nameFromPrompt("")).toBe("")
  })
})

describe("/rename in the window", () => {
  test("names the session in view, and sends no prompt", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    await line(w, "/rename Parser fix")
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "s1", title: "Parser fix" }])
    expect(sent(s.calls, "session.prompt")).toEqual([])
    expect(w.footer().join("\n")).toContain("renamed: Parser fix")
  })

  test("without a name, after the session's first prompt", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    asked(w, "Read notes.txt and tell me in one sentence what it says.")
    await line(w, "/rename")
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "s1", title: "Read notes.txt and tell me" }])
  })

  test("without a name and nothing asked yet, it says so", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    await line(w, "/rename")
    expect(sent(s.calls, "session.setTitle")).toEqual([])
    expect(sent(s.calls, "session.prompt")).toEqual([])
    expect(w.footer().join("\n")).toContain("nothing to name it after yet")
  })

  test("while a turn runs: at once, not queued behind it", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
    await w.settle()
    await line(w, "/rename mid turn")
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "s1", title: "mid turn" }])
    expect(sent(s.calls, "session.prompt")).toEqual([])
  })

  test("the name stands on the prompt box's top border, at the right", async () => {
    const w = open = await mountWindow({ width: 80 })
    w.store.push({ type: "session.updated", session: "s1", ts: 0, seq: 5, title: "Parser fix" } as never)
    await w.settle()
    const top = w.footer().find((l) => l.includes(" Parser fix "))!
    expect(top).toBeDefined()
    // ───…─── Parser fix ─ : the rule, the name, one dash; a column short of the edge.
    expect(top.trimEnd()).toMatch(/^─{10,} Parser fix ─$/)
    expect(top.trimEnd().length).toBe(80 - 2)
  })

  test("a long name is cut on the border, which keeps its dashes", async () => {
    const w = open = await mountWindow({ width: 40 })
    w.store.push({ type: "session.updated", session: "s1", ts: 0, seq: 5, title: "a very long name for a session that goes on" } as never)
    await w.settle()
    const top = w.footer().find((l) => l.includes("─ a very"))!
    expect(top.trimEnd().length).toBe(40 - 2)
    expect(top).toContain("…")
    expect(top.trimEnd().endsWith(" ─")).toBe(true)
  })

  test("is offered by completion and by the palette", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    await w.t.mockInput.typeText("/ren")
    await Bun.sleep(30)
    await w.settle()
    expect(w.footer().join("\n")).toContain("/rename")
  })

  test("the terminal's title names the session in view", async () => {
    const titles: string[] = []
    const w = open = await mountWindow()
    const renderer = w.t.renderer as unknown as { setTerminalTitle(t: string): void }
    const original = renderer.setTerminalTitle.bind(renderer)
    renderer.setTerminalTitle = (t: string) => { titles.push(t); original(t) }
    w.store.push({ type: "session.updated", session: "s1", ts: 0, seq: 5, title: "Parser fix" } as never)
    await w.settle()
    expect(titles.at(-1)).toBe("DELFIN · Parser fix")
  })
})

describe("/rename, where it could go astray", () => {
  test("with the key that steers, mid-turn: renamed, never sent into the turn", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
    await w.settle()
    await w.t.mockInput.typeText("/rename steered name")
    await w.settle()
    w.t.mockInput.pressKey("g", { ctrl: true })
    await w.settle()
    expect(sent(s.calls, "session.setTitle")).toEqual([{ sessionId: "s1", title: "steered name" }])
    expect(sent(s.calls, "session.prompt")).toEqual([])
  })

  test("a session command of the same name: completion lists the window's, which enter runs", async () => {
    const calls: Call[] = []
    const client = {
      request: async (method: string, params: Record<string, any> = {}) => {
        calls.push({ method, params })
        if (method === "session.commands") {
          return { commands: [{ name: "/rename", summary: "a user command of its own", category: "user", takesArgs: true, source: "user" }] }
        }
        return { title: params.title }
      },
    }
    const w = open = await mountWindow({ client })
    await w.t.mockInput.typeText("/ren")
    await Bun.sleep(30)
    await w.settle()
    const shown = w.footer().join("\n")
    expect(shown).toContain("name this session")
    expect(shown).not.toContain("a user command of its own")
  })

  test("a service that cannot be asked: the line stays", async () => {
    const client = {
      request: async (method: string) => {
        if (method === "session.setTitle") throw new Error("the agent service is not reachable")
        return {}
      },
    }
    const w = open = await mountWindow({ client })
    await line(w, "/rename kept")
    expect(w.footer().join("\n")).toContain("❯ /rename kept")
  })

  test("the palette offers it", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client })
    w.t.mockInput.pressKey("p", { ctrl: true })
    await Bun.sleep(30)
    await w.settle()
    expect(w.footer().join("\n")).toContain("rename…")
  })

  test("a name in wide characters keeps a divider on one row", async () => {
    const w = open = await mountWindow({ width: 60 })
    w.view.show("s1", { divider: "now: 二酸化炭素の振動解析" })
    await w.settle()
    const rows = w.t.externalOutput.takeText().split("\n").filter((l) => l.includes("━"))
    expect(rows.length).toBe(1)
    expect(Bun.stringWidth(rows[0]!.trimEnd())).toBe(60 - 1)
  })
})

describe("the hints line, narrow", () => {
  test("keeps the posture when a connection note is long", async () => {
    const w = open = await mountWindow({ width: 40 })
    w.store.update((d) => {
      d.sessions.s1!.approval = "plan"
      d.connection = "reconnecting"
    })
    await w.settle()
    expect(w.footer().join("\n")).toContain("⏸ plan")
  })
})
