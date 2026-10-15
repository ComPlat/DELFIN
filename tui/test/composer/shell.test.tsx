/**
 * Shell mode: a "!" at the start of the prompt makes the line a shell
 * command (the session runs it through the agent's checks). The composer
 * shows it -- "!" as the marker in the shell's colour, the box's border
 * saying what the line is and how to leave -- and offers no completion;
 * backspace on the "!" or esc leave it; what is sent is the line exactly
 * as typed, "!" and all. Drawn headless through the real App.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { parseColor, rgbToHex, type CapturedSpan } from "@opentui/core"
import { SHELL_PREFIX, shellLead } from "../../src/commands"
import { SHELL_HINTS, shellHint } from "../../src/footer/footer"
import { shellCommand } from "../../src/scrollback/rows"
import { DARK, MONO, themeFor, type Theme } from "../../src/theme/theme"
import { width } from "../../src/text"
import { mountWindow, type Call } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const COMMANDS = [
  { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" },
]
const hex = (c: Parameters<typeof rgbToHex>[0]) => rgbToHex(c).toLowerCase()
const same = (a: Parameters<typeof rgbToHex>[0], b: string) => hex(a) === hex(parseColor(b))

async function mount(options: { width?: number; theme?: () => Theme; full?: boolean; history?: string[] } = {}) {
  const calls: Call[] = []
  w = await mountWindow({
    width: options.width ?? 100, theme: options.theme, full: options.full,
    client: {
      request: async (method: string, params: Record<string, any> = {}) => {
        calls.push({ method, params })
        if (method === "session.commands") return { commands: COMMANDS }
        if (method === "session.completePath") return { paths: ["delfin/", "docs/"] }
        return { accepted: true }
      },
    },
  })
  await w.settle()
  const settle = async () => {
    for (let i = 0; i < 3; i++) {
      await Bun.sleep(5)
      await w!.settle()
    }
  }
  const type = async (text: string) => {
    await w!.keys.typeText(text)
    await settle()
  }
  const footer = () => w!.footer().join("\n")
  const prompt = () => w!.footer().find((l) => l.startsWith("❯ ") || l.startsWith("! "))?.trimEnd() ?? ""
  const prompts = () => calls.filter((c) => c.method === "session.prompt").map((c) => c.params)
  const spans = (): CapturedSpan[] => w!.t.captureSpans().lines.flatMap((l) => l.spans)
  return { w, calls, settle, type, footer, prompt, prompts, spans }
}

describe("shell mode", () => {
  test("a ! at the start: the marker, the border and the placeholder say shell", async () => {
    const m = await mount()
    expect(m.prompt()).toStartWith("❯ ")
    await m.type("!")
    expect(m.prompt()).toBe("! a shell command…")
    expect(m.footer()).toContain("── shell · runs through the agent's checks · esc or backspace on an empty line to leave ─")
    // The marker and the word in the shell's colour, the rules too.
    const marker = m.spans().find((s) => s.text.startsWith("! "))!
    expect(same(marker.fg, DARK.shell as string)).toBe(true)
    const word = m.spans().find((s) => s.text === "shell")!
    expect(same(word.fg, DARK.shell as string)).toBe(true)
    await m.type("ls -la")
    expect(m.prompt()).toBe("! ls -la")             // the "!" is the marker, not the text
  })

  test("enter sends the line as typed, ! and all, and leaves shell mode", async () => {
    const m = await mount()
    await m.type("!ls -la")
    m.w.keys.pressEnter()
    await m.settle()
    expect(m.prompts()).toEqual([{ sessionId: "s1", text: "!ls -la" }])
    expect(m.prompt()).toStartWith("❯ ")
    expect(m.footer()).not.toContain("shell ·")
    // A blank after the "!" stays, as it was sent before the mode existed.
    await m.type("! git status")
    m.w.keys.pressEnter()
    await m.settle()
    expect(m.prompts().at(-1)).toEqual({ sessionId: "s1", text: "! git status" })
  })

  test("a ! with nothing after it sends nothing", async () => {
    const m = await mount()
    await m.type("!   ")
    m.w.keys.pressEnter()
    await m.settle()
    expect(m.prompts()).toEqual([])
    expect(m.prompt()).toBe("!")                    // still a shell line
    expect(m.footer()).toContain("── shell ·")
  })

  test("no completion in shell mode: no command list, no paths asked", async () => {
    const m = await mount()
    await m.type("/st")
    expect(m.footer()).toContain("/status")         // a prompt: the list
    m.w.keys.pressKey("\x15")                       // ctrl+u
    for (let i = 0; i < 3; i++) m.w.keys.pressBackspace()
    await m.settle()
    const before = m.calls.length
    await m.type("!/st")
    expect(m.footer()).not.toContain("Model, mode, tokens, cost")
    await m.type(" @del")
    m.w.keys.pressTab()
    await m.settle()
    expect(m.footer()).not.toContain("delfin/")
    expect(m.calls.slice(before).filter((c) => c.method === "session.completePath")).toEqual([])
    expect(m.prompt()).toBe("! /st @del")
  })

  test("backspace on the ! alone, or esc, goes back to a prompt", async () => {
    const m = await mount()
    await m.type("!")
    m.w.keys.pressBackspace()
    await m.settle()
    expect(m.prompt()).toStartWith("❯ ")
    expect(m.footer()).not.toContain("shell ·")
    expect(m.footer()).toContain("ask DELFIN")
    // esc with a command typed: the "!" goes, the words stay -- a prompt now.
    await m.type("!uname -a")
    await m.w.escape()
    expect(m.prompt()).toBe("❯ uname -a")
    expect(m.footer()).not.toContain("shell ·")
    m.w.keys.pressEnter()
    await m.settle()
    expect(m.prompts()).toEqual([{ sessionId: "s1", text: "uname -a" }])
  })

  test("backspace inside the command edits it; only at the start does it leave", async () => {
    const m = await mount()
    await m.type("!lss")
    m.w.keys.pressBackspace()
    await m.settle()
    expect(m.prompt()).toBe("! ls")
    m.w.keys.pressBackspace()
    m.w.keys.pressBackspace()
    await m.settle()
    expect(m.prompt()).toStartWith("! a shell command")
    m.w.keys.pressBackspace()
    await m.settle()
    expect(m.prompt()).toStartWith("❯ ")
  })

  test("a ! typed at the start of a draft makes it a shell line; the history brings one back as one", async () => {
    const m = await mount()
    await m.type("pwd")
    m.w.keys.pressKey("\x01")                       // ctrl+a: the start of the line
    await m.type("!")
    expect(m.prompt()).toBe("! pwd")
    m.w.keys.pressEnter()
    await m.settle()
    expect(m.prompts()).toEqual([{ sessionId: "s1", text: "!pwd" }])
    m.w.keys.pressArrow("up")
    await m.settle()
    expect(m.prompt()).toBe("! pwd")
    expect(m.footer()).toContain("shell ·")
  })

  test("a pasted !line is a shell line too", async () => {
    const m = await mount()
    await m.w.keys.pasteBracketedText("!echo hi")
    await m.settle()
    expect(m.prompt()).toBe("! echo hi")
    m.w.keys.pressEnter()
    await m.settle()
    expect(m.prompts()).toEqual([{ sessionId: "s1", text: "!echo hi" }])
  })

  test("the steer key does not steer a shell line into the turn: it is queued", async () => {
    const m = await mount()
    m.w.store.update((d) => {
      d.sessions.s1!.turn = { id: "t1", source: "user", startedAt: Date.now(), label: "thinking…",
                              base: { inputTokens: 0, outputTokens: 0, costUsd: 0 }, streamed: 0 }
    })
    await m.settle()
    await m.type("!")
    expect(m.prompt()).toBe("! a shell command… (runs when the turn ends)")
    await m.type("ls")
    m.w.keys.pressKey("\x07")                       // ctrl+g
    await m.settle()
    expect(m.prompts()).toEqual([{ sessionId: "s1", text: "!ls" }])
    expect(m.footer()).toContain("a shell line is not steered")
  })

  test("in the agent view a ! is text (a name, an answer, a message), not shell mode", async () => {
    const m = await mount()
    m.w.keys.pressArrow("left")                     // the agent view
    await m.settle()
    await m.type("!x")
    expect(m.footer()).toContain("❯ !x")
    expect(m.footer()).not.toContain("shell ·")
  })

  test("full screen: the same marker and border", async () => {
    const m = await mount({ full: true })
    await m.type("!df -h")
    const lines = m.w.screenLines()
    expect(lines.some((l) => l.startsWith("! df -h"))).toBe(true)
    expect(lines.some((l) => l.includes("── shell · runs through the agent's checks"))).toBe(true)
  })

  test("NO_COLOR (mono): the marker and the word are text, in the terminal's own colour", async () => {
    const m = await mount({ theme: () => MONO })
    await m.type("!ls")
    expect(m.prompt()).toBe("! ls")
    expect(m.footer()).toContain("── shell · runs through the agent's checks")
    expect(themeFor("mono").shell).toBe(MONO.shell)
    expect(MONO.shell).toBe(MONO.text)                // no colour of the window's own
  })

  test("a narrow terminal: the border's words shorten, the marker stays", async () => {
    const m = await mount({ width: 44 })
    await m.type("!ls")
    expect(m.prompt()).toBe("! ls")
    expect(m.footer()).toContain("── shell · esc to leave ─")
    expect(m.footer()).not.toContain("agent's checks")
  })
})

describe("shell mode, the parts", () => {
  test("the prefix is the session's (repl_commands.SHELL_PREFIX)", async () => {
    const py = await Bun.file(`${import.meta.dir}/../../../delfin/agent/repl_commands.py`).text()
    const m = /^SHELL_PREFIX = "(.*)"$/m.exec(py)
    expect(m?.[1]).toBe(SHELL_PREFIX)
  })

  test("shellLead: a leading ! (after blanks), nowhere else", () => {
    expect(shellLead("!ls")).toBe(1)
    expect(shellLead("  !ls")).toBe(3)
    expect(shellLead("ls !x")).toBe(-1)
    expect(shellLead("")).toBe(-1)
  })

  test("shellHint: the longest that fits, down to the word", () => {
    expect(shellHint(200)).toBe(SHELL_HINTS[0])
    for (const hint of SHELL_HINTS) expect(shellHint(width(hint))).toBe(hint)
    expect(shellHint(5)).toBe("shell")
    expect(shellHint(4)).toBe("")
  })

  test("shellCommand: a line the session runs as a command, as repl.py reads it", () => {
    expect(shellCommand("!ls -la")).toBe("ls -la")
    expect(shellCommand("  ! git status ")).toBe("git status")
    expect(shellCommand("!")).toBeNull()               // the model's, as before
    expect(shellCommand("why ! here")).toBeNull()
  })

  test("the scrollback heads a shell line with $ in the shell's colour", async () => {
    const m = await mount({ full: true })
    m.w.store.push({ type: "user.message", session: "s1", ts: 0, seq: 900, text: "!ls -la", source: "user" } as never)
    m.w.store.push({ type: "notice", session: "s1", ts: 0, seq: 901, text: "  README.md\n  delfin" } as never)
    await m.settle()
    const out = m.w.screenLines().join("\n")
    expect(out).toContain("$ ls -la")
    expect(out).not.toContain("❯ !ls -la")
    expect(out).toContain("README.md")
    const dollar = m.spans().find((s) => s.text.startsWith("$"))!
    expect(same(dollar.fg, DARK.shell as string)).toBe(true)
  })
})
