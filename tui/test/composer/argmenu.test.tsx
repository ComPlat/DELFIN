/**
 * The menu after `/command `, in the window: `/` lists the commands with
 * their one-line summaries; `/model ` the served models with the current
 * one marked; `/resume ` the sessions; a command whose values cannot be
 * listed its usage line. ↑↓ and enter pick, a click picks.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const COMMANDS = [
  { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" },
  { name: "/model", summary: "Show or switch the model", category: "setup", takesArgs: true, source: "builtin" },
  { name: "/compact", summary: "Summarise history to free space", category: "session", takesArgs: true, source: "builtin" },
  { name: "/git", summary: "status | diff | log | branch | role", category: "workspace", takesArgs: true, source: "builtin" },
]
const calls: { method: string; params: any }[] = []
const client = {
  request: async (method: string, params: Record<string, any> = {}) => {
    calls.push({ method, params })
    if (method === "session.commands") return { commands: COMMANDS }
    if (method === "session.choices" && params.what === "args") {
      if (params.command === "/model") {
        return { current: "kit.glm-5", choices: [{ value: "kit.glm-5", detail: "KIT" },
                                                  { value: "kit.deepseek-v4", detail: "KIT · cheap" }], usage: "" }
      }
      if (params.command === "/compact") return { current: "", choices: [], usage: "/compact [what the summary must keep]" }
      if (params.command === "/git") return { current: "", choices: ["status", "diff", "log"].map((value) => ({ value })), usage: "" }
      return { current: "", choices: [], usage: "" }
    }
    if (method === "sessions.list") return { sessions: [{ sessionId: "s1", title: "this one", workspace: "/w", state: "idle" },
                                                         { sessionId: "s2", title: "scan of benzene", workspace: "/w2", state: "idle", model: "m" }] }
    if (method === "sessions.past") return { sessions: [] }
    return { accepted: true }
  },
}
const prompts = () => calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)

async function mount(resumed: string[] = []) {
  calls.length = 0
  w = await mountWindow({ full: true, client, onResume: (id) => resumed.push(id) })
  w.store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, sessionId: "s1", workspace: "/w",
                 provider: "kit", model: "kit.glm-5", mode: "solo", approval: "default", banner: [], notices: [] } as never)
  await w.settle()
  return w
}
const type = async (text: string) => {
  w!.keys.typeText(text)
  await Bun.sleep(30)
  await w!.settle()
  await Bun.sleep(10)
  await w!.settle()
}
const screen = () => w!.t.captureCharFrame()

describe("the menu after a command", () => {
  test("/ lists every command with its summary", async () => {
    await mount()
    await type("/")
    for (const c of COMMANDS) expect(screen()).toContain(c.name)
    expect(screen()).toContain("Summarise history to free space")
    expect(screen()).toContain("/resume")                       // the window's own
  })

  test("/model : the served models, the current one marked; enter on another switches", async () => {
    await mount()
    await type("/model ")
    expect(screen()).toMatch(/kit\.glm-5\s+\(now\) KIT/)
    expect(screen()).toContain("kit.deepseek-v4")
    w!.keys.pressArrow("down")
    await w!.settle()
    w!.keys.pressEnter()
    await Bun.sleep(20)
    await w!.settle()
    expect(prompts()).toEqual(["/model kit.deepseek-v4"])
  })

  test("typing narrows it; a click picks", async () => {
    await mount()
    await type("/git d")
    expect(screen()).toContain("diff")
    expect(screen()).not.toMatch(/›\s+status/)
    const rows = screen().split("\n")
    const y = rows.findIndex((r) => /\s{2}diff\s*$/.test(r) || r.trimEnd().endsWith(" diff"))
    await w!.t.mockMouse.click(rows[y]!.indexOf("diff") + 1, y)
    await Bun.sleep(20)
    await w!.settle()
    expect(prompts()).toEqual(["/git diff"])
  })

  test("choosing a command that takes an argument opens its menu", async () => {
    await mount()
    await type("/mod")
    w!.keys.pressEnter()
    await Bun.sleep(30)
    await w!.settle()
    await Bun.sleep(10)
    await w!.settle()
    expect(screen()).toContain("❯ /model")
    expect(screen()).toContain("kit.deepseek-v4")
    expect(prompts()).toEqual([])
  })

  test("a command whose values cannot be listed: its usage, and enter sends what is typed", async () => {
    await mount()
    await type("/compact keep the energies")
    expect(screen()).toContain("usage: /compact [what the summary must keep]")
    w!.keys.pressEnter()
    await Bun.sleep(20)
    await w!.settle()
    expect(prompts()).toEqual(["/compact keep the energies"])
  })

  test("/resume : the sessions by name; picking one opens it in this window", async () => {
    const resumed: string[] = []
    await mount(resumed)
    await type("/resume benz")
    expect(screen()).toContain("scan of benzene")
    w!.keys.pressEnter()
    await Bun.sleep(20)
    await w!.settle()
    expect(resumed).toEqual(["s2"])
    expect(prompts()).toEqual([])
  })
})
