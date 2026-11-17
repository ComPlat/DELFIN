/**
 * An unknown /command is answered in the window and never sent -- not a
 * paid turn about a typo (Claude Code; the dashboard's "Unknown command";
 * the REPL's "did you mean"). The line stays in the prompt to be put right.
 * A path ("/home/…") is no command and goes out; so does every line when
 * the session cannot say which commands it has.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mountWindow, type Call } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await open?.close(); open = null })

const COMMANDS = [
  { name: "/help", summary: "List these commands", category: "session", takesArgs: true, source: "builtin" },
  { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" },
  { name: "/keys", summary: "Keys", category: "session", takesArgs: false, source: "builtin" },
]

function service(commands: unknown[] | null = COMMANDS) {
  const calls: Call[] = []
  return {
    calls,
    client: {
      request: async (method: string, params: Record<string, any> = {}) => {
        calls.push({ method, params })
        if (method === "session.commands") {
          if (commands === null) throw new Error("no method 'session.commands'")
          return { commands }
        }
        return { accepted: true, applied: true, message: "" }
      },
    },
  }
}

const prompts = (calls: Call[]) => calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text)

async function line(w: NonNullable<typeof open>, text: string) {
  await w.t.mockInput.typeText(text)
  w.t.mockInput.pressEscape()                    // the completion list, if any
  await Bun.sleep(40)
  w.t.mockInput.pressEnter()
  await Bun.sleep(20)
  await w.settle()
}

describe("an unknown /command in the window", () => {
  test("is said, not sent, and the line stays", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client, height: 30 })
    await line(w, "/nonexistent")
    expect(prompts(s.calls)).toEqual([])
    const footer = w.footer().join("\n")
    expect(footer).toContain("Unknown command: /nonexistent")
    expect(footer).toContain("❯ /nonexistent")
  })

  test("a typo names the command it was probably meant to be", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client, height: 30 })
    await line(w, "/stauts")
    expect(prompts(s.calls)).toEqual([])
    expect(w.footer().join("\n")).toContain("did you mean /status?")
  })

  test("two commands run together are not sent either", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client, height: 30 })
    await line(w, "/help/keys")
    expect(prompts(s.calls)).toEqual([])
    expect(w.footer().join("\n")).toContain("Unknown command: /help/keys")
  })

  test("a known command goes to the session", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client, height: 30 })
    await line(w, "/status")
    expect(prompts(s.calls)).toEqual(["/status"])
  })

  test("a path is no command: it goes out as typed", async () => {
    const s = service()
    const w = open = await mountWindow({ client: s.client, height: 30 })
    await line(w, "/home/someone/run.py fails, why?")
    expect(prompts(s.calls)).toEqual(["/home/someone/run.py fails, why?"])
  })

  test("a session that cannot list its commands: nothing is judged", async () => {
    const s = service(null)
    const w = open = await mountWindow({ client: s.client, height: 30 })
    await line(w, "/nonexistent")
    expect(prompts(s.calls)).toEqual(["/nonexistent"])
  })
})
