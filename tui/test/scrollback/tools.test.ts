/**
 * Tool rows read as in the REPL (delfin/agent/repl_render.py): the
 * argument that matters, an edit as what it did plus its diff, a call
 * that did not succeed as "blocked", a shell envelope as its exit code
 * and its last lines, a cut result as cut.
 */
import { describe, expect, test } from "bun:test"
import {
  EDIT_DIFF_ROWS, cutEnvelopeCarrying, envelopeCarrying, extractDiff, resultBody, shortToolName,
  toolHeadline, toolResult,
} from "../../src/scrollback/tools"

describe("a tool's headline", () => {
  test("names the argument that matters, on one line", () => {
    expect(toolHeadline("bash", { command: "pytest -q\n  --lf", description: "x" })).toEqual({ name: "bash", arg: "pytest -q --lf" })
    expect(toolHeadline("grep_file", { path: "delfin/", pattern: "ENERGY" }).arg).toBe("ENERGY")
    expect(toolHeadline("subagent", { subagent_type: "explore", description: "find the parser" }).arg).toBe("find the parser")
    expect(toolHeadline("web_fetch", JSON.stringify({ url: "https://example.org" })).arg).toBe("https://example.org")
    expect(toolHeadline("odd_tool", { a: [1], b: {} }).arg).toBe("(2 args)")
  })

  test("an MCP tool keeps its server in view", () => {
    expect(shortToolName("mcp__delfin-docs__read_file")).toBe("delfin-docs:read_file")
    expect(toolHeadline("mcp__delfin-docs__read_file", { path: "README.md" })).toEqual({ name: "delfin-docs:read_file", arg: "README.md" })
  })
})

describe("a tool's result", () => {
  test("a call that did not succeed is blocked, with why", () => {
    expect(toolResult("read_file", "", { ok: false, error: "outside\nthe workspace" }))
      .toEqual({ kind: "blocked", reason: "outside the workspace" })
  })

  test("a shell envelope is its exit code and its last lines", () => {
    const r = toolResult("bash", JSON.stringify({ exit_code: 1, stdout: "a\nb\nc\nd\n", stderr: "boom" }))
    expect(r).toEqual({ kind: "text", summary: "exit 1 · 5 lines", lines: ["c", "d", "boom"], ok: false, exitCode: 1 })
  })

  test("other output is its size, and its last three lines", () => {
    const r = toolResult("read_file", "one\ntwo\nthree\nfour\n", { chars: 19 })
    expect(r).toEqual({ kind: "text", summary: "4 lines, 19 B", lines: ["two", "three", "four"], ok: true })
    const cut = toolResult("read_file", "x\n".repeat(40), { chars: 9000, truncated: true })
    expect(cut.kind === "text" && cut.summary).toBe("40+ lines, 8.8 kB (truncated)")
    expect(toolResult("list_files", "", {})).toMatchObject({ summary: "(no output)" })
  })

  test("an envelope cut short says so, and shows its first lines", () => {
    const output = '{"exit_code": 0, "stdout": "' + "line\\n".repeat(400)
    const r = toolResult("bash", output)
    expect(r.kind === "text" && r.summary).toContain("cut for display")
    expect(cutEnvelopeCarrying('{"exit_code": 0, "stdout": "ok\\n", "cw')).toEqual(["ok\n", true])
    expect(envelopeCarrying('{"stdout": "136\\n", "cwd": "."}')).toBe("136\n")
    expect(resultBody('{"stdout": "136\\n"}')).toBe("136\n")
  })

  test("an edit is what it did, its counts and its numbered diff", () => {
    const diff = "--- a/delfin/xtb.py\n+++ b/delfin/xtb.py\n@@ -41,3 +41,4 @@\n def f():\n-    old\n+    new\n+    more\n     same\n\nTip: run the tests"
    const r = toolResult("edit_file", diff, { ok: true }, { path: "delfin/xtb.py" })
    expect(r.kind).toBe("edit")
    if (r.kind !== "edit") return
    expect([r.action, r.label, r.adds, r.dels]).toEqual(["Edited", "delfin/xtb.py", 2, 1])
    expect(r.rows.map((row) => row.text)).toEqual([
      "@@ -41,3 +41,4 @@", "   41  def f():", "   42 -    old", "   42 +    new", "   43 +    more", "   44      same"])
    expect(extractDiff(diff).at(-1)).toBe("     same")      // the advice after it is not the diff
  })

  test("a new file is Created; a long diff is cut and counted", () => {
    const lines = Array.from({ length: 60 }, (_, i) => `+line ${i}`)
    const diff = `--- /dev/null\n+++ b/new.py\n@@ -0,0 +1,60 @@\n${lines.join("\n")}`
    const r = toolResult("write_file", diff, {}, { path: "new.py" })
    expect(r.kind === "edit" && [r.action, r.rows.length, r.omitted]).toEqual(["Created", EDIT_DIFF_ROWS, 61 - EDIT_DIFF_ROWS])
  })

  test("apply_patch shows the patch it applied, from its input", () => {
    const patch = "*** Begin Patch\n*** Update File: a.py\n@@\n-x = 1\n+x = 2\n*** End Patch"
    const r = toolResult("apply_patch", JSON.stringify({ status: "applied", files_touched: ["a.py"] }), {}, { patch })
    expect(r.kind === "edit" && [r.label, r.adds, r.dels]).toEqual(["a.py", 1, 1])
    const staged = toolResult("apply_patch", JSON.stringify({ status: "staged" }), {}, { patch })
    expect(staged.kind).toBe("text")
  })
})
