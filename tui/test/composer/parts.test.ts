/**
 * The composer's parts on their own: the history file and its walk, paste
 * placeholders, and what completion finds to complete.
 */
import { describe, expect, test } from "bun:test"
import { mkdtempSync, readFileSync, statSync, writeFileSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"
import { applyCandidate, commonPrefix, matchCommands, pathCandidates, tokenAt, type Command } from "../../src/composer/complete"
import { History, HISTORY_MAX, HistoryWalk, historyPaths, parseOwn, parseReadline } from "../../src/composer/history"
import { Pastes, shouldCollapse } from "../../src/composer/paste"

const tmp = () => mkdtempSync(join(tmpdir(), "tui-history-"))

describe("the prompt history", () => {
  test("reads the REPL's readline file in front of its own, and keeps prompts that span lines", () => {
    const dir = tmp()
    const paths = { repl: join(dir, "agent_repl_history"), own: join(dir, "tui_history.jsonl") }
    writeFileSync(paths.repl, "first in chat\n#1712345678\nsecond in chat\n")
    writeFileSync(paths.own, JSON.stringify({ text: "line one\nline two" }) + "\n{torn")
    const h = History.load(paths)
    expect(h.entries.map((e) => e.text)).toEqual(["first in chat", "second in chat", "line one\nline two"])
  })

  test("reads libedit's file (macOS) with its escapes", () => {
    expect(parseReadline("_HiStOrY_V2_\nls\\040-la\n")).toEqual(["ls -la"])
  })

  test("appends each prompt as one line, private to the user, never twice in a row", () => {
    const dir = tmp()
    const paths = { repl: join(dir, "none"), own: join(dir, "sub", "tui_history.jsonl") }
    const h = History.load(paths)
    h.add({ text: "hello" })
    h.add({ text: "hello" })
    h.add({ text: "see [Pasted text #1 +9 lines]", pastes: { "[Pasted text #1 +9 lines]": "a\nb" } })
    const lines = readFileSync(paths.own, "utf8").trim().split("\n")
    expect(lines.length).toBe(2)
    expect(parseOwn(lines.join("\n"))[1]!.pastes).toEqual({ "[Pasted text #1 +9 lines]": "a\nb" })
    expect(statSync(paths.own).mode & 0o777).toBe(0o600)
    expect(History.load(paths).entries.map((e) => e.text)).toEqual(["hello", "see [Pasted text #1 +9 lines]"])
  })

  test("keeps the newest entries when it grows past its limit", () => {
    const h = new History(null)
    for (let i = 0; i < HISTORY_MAX + 5; i++) h.add({ text: `p${i}` })
    expect(h.entries.length).toBe(HISTORY_MAX)
    expect(h.entries[0]!.text).toBe("p5")
  })

  test("follows DELFIN_SCRATCH_STATE for its own file, never for the REPL's", () => {
    const p = historyPaths({ HOME: "/h", DELFIN_SCRATCH_STATE: "/s" })
    expect(p).toEqual({ own: "/s/.delfin/tui_history.jsonl", repl: "/h/.delfin/agent_repl_history" })
  })

  test("walks back and forth, and brings the draft back past the newest", () => {
    const h = new History(null, [{ text: "one" }, { text: "two" }])
    const walk = new HistoryWalk(h)
    expect(walk.older({ text: "draft" })?.text).toBe("two")
    expect(walk.older({ text: "two" })?.text).toBe("one")
    expect(walk.older({ text: "one" })).toBe(null)           // the oldest
    expect(walk.newer()?.text).toBe("two")
    expect(walk.newer()?.text).toBe("draft")
    expect(walk.walking).toBe(false)
    expect(walk.newer()).toBe(null)
  })

  test("searches backwards, case-insensitively", () => {
    const h = new History(null, [{ text: "run ORCA opt" }, { text: "ls" }, { text: "orca freq" }])
    expect(h.search("orca")).toBe(2)
    expect(h.search("orca", 2)).toBe(0)
    expect(h.search("orca", 0)).toBe(-1)
  })
})

describe("paste placeholders", () => {
  test("a long paste stands in as one line and goes out whole", () => {
    const text = Array.from({ length: 40 }, (_, i) => `row ${i}`).join("\n")
    expect(shouldCollapse(text)).toBe(true)
    expect(shouldCollapse("one\ntwo")).toBe(false)
    const p = new Pastes()
    const label = p.add(text)
    expect(label).toBe("[Pasted text #1 +40 lines]")
    expect(p.expand(`look at ${label} please`)).toBe(`look at ${text} please`)
    expect(p.add("x\n".repeat(10))).toBe("[Pasted text #2 +11 lines]")
  })

  test("an edited placeholder stays as typed; a recalled one comes back", () => {
    const p = new Pastes()
    const label = p.add("a\n".repeat(20))
    expect(p.expand(label.replace("#1", "#9"))).toBe(label.replace("#1", "#9"))
    expect(p.used(`x ${label}`)).toEqual({ [label]: "a\n".repeat(20) })
    const q = new Pastes()
    q.restore({ [label]: "content" })
    expect(q.expand(label)).toBe("content")
    expect(q.add("b\n".repeat(20)).startsWith("[Pasted text #2")).toBe(true)
  })
})

describe("completion", () => {
  const commands: Command[] = [
    { name: "/help", summary: "List these commands", category: "session", takesArgs: true, source: "builtin" },
    { name: "/status", summary: "Model, mode, tokens, cost", category: "session", takesArgs: false, source: "builtin" },
    { name: "/session", summary: "List or search past sessions", category: "history", takesArgs: true, source: "builtin" },
    { name: "/memories", summary: "Stored memories", category: "history", takesArgs: true, source: "builtin" },
  ]

  test("a /command only as the first word, and never a path", () => {
    expect(tokenAt("/he", 3)).toEqual({ kind: "command", text: "/he", start: 0, end: 3 })
    expect(tokenAt("say /he", 7)).toBe(null)
    expect(tokenAt("/home/user/x", 12)).toBe(null)
    expect(tokenAt("/help me", 8)).toBe(null)
  })

  test("an @path anywhere, the whole word under the cursor", () => {
    expect(tokenAt("look at @delfin/x", 17)).toEqual({ kind: "path", text: "@delfin/x", start: 8, end: 17 })
    expect(tokenAt("look at @delfin/xtb.py now", 12)).toEqual({ kind: "path", text: "@delfin/xtb.py", start: 8, end: 22 })
    expect(tokenAt("mail@example", 12)).toBe(null)
  })

  test("names that start with what is typed come first, then names that contain it", () => {
    expect(matchCommands(commands, "/s").map((c) => c.label)).toEqual(["/status", "/session"])
    expect(matchCommands(commands, "/mem").map((c) => c.label)).toEqual(["/memories"])
    expect(matchCommands(commands, "/ess").map((c) => c.label)).toEqual(["/session"])
    const [status] = matchCommands(commands, "/sta")
    expect(status!.runs).toBe(true)                 // no arguments: Enter runs it
    expect(matchCommands(commands, "/he")[0]!.insert).toBe("/help ")
  })

  test("a path goes in with a space after it; a folder goes on", () => {
    const [file, dir] = pathCandidates(["delfin/xtb.py", "docs/"])
    const token = tokenAt("see @del", 8)!
    expect(applyCandidate("see @del", token, file!.insert)).toEqual({ text: "see @delfin/xtb.py ", cursor: 19 })
    expect(dir!.descends).toBe(true)
    expect(applyCandidate("@d", tokenAt("@d", 2)!, dir!.insert)).toEqual({ text: "@docs/", cursor: 6 })
    expect(commonPrefix(pathCandidates(["delfin/xtb.py", "delfin/orca.py"]))).toBe("@delfin/")
  })
})
