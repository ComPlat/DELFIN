/**
 * After `/command ` the composer offers what the command takes: the
 * argument at the cursor is a token of its own, matched against the
 * command's values (every typed word, in any order), the one in force
 * marked; a command whose values cannot be listed shows its usage line,
 * which Enter never takes.
 */
import { describe, expect, test } from "bun:test"
import { applyCandidate, argumentAt, matchArguments } from "../../src/composer/complete"

describe("the argument token", () => {
  test("after a command and a space, to the end of the line", () => {
    expect(argumentAt("/model ", 7)).toEqual({ kind: "argument", command: "/model", text: "", start: 7, end: 7 })
    expect(argumentAt("/jobs watch o", 13)).toEqual({ kind: "argument", command: "/jobs", text: "watch o",
                                                       start: 6, end: 13 })
    expect(argumentAt("/model", 6)).toBe(null)                 // still the command
    expect(argumentAt("say /model x", 12)).toBe(null)          // a command only at the start
    expect(argumentAt("/model x\ny", 10)).toBe(null)           // one line only
    expect(argumentAt("/model abc", 3)).toBe(null)             // the cursor in the command
    expect(argumentAt("/home/user/x y", 14)).toBe(null)        // a path, not a command
  })
})

describe("the argument's values", () => {
  const values = [{ value: "kit.glm-5", detail: "KIT" }, { value: "kit.deepseek-v4", detail: "KIT · cheap" },
                  { value: "gpt-5" }]
  test("every one for nothing typed, the current one marked; narrowed by every typed word", () => {
    const all = matchArguments({ current: "gpt-5", choices: values, usage: "" }, "")
    expect(all.map((c) => c.label)).toEqual(["kit.glm-5", "kit.deepseek-v4", "gpt-5"])
    expect(all[2]!.detail).toBe("(now)")
    expect(all.every((c) => c.runs && !c.hint)).toBe(true)
    expect(matchArguments({ current: "", choices: values, usage: "" }, "kit cheap").map((c) => c.insert))
      .toEqual(["kit.deepseek-v4"])
  })

  test("a label of its own (a session's title) inserts the value", () => {
    const got = matchArguments({ current: "", choices: [{ value: "abc123", label: "relax Fe", detail: "running" }],
                                 usage: "" }, "relax")
    expect(got).toEqual([{ insert: "abc123", label: "relax Fe", detail: "running", runs: true, descends: false }])
  })

  test("no values: the usage line, as a hint", () => {
    const got = matchArguments({ current: "", choices: [], usage: "/compact [what the summary must keep]" }, "keep x")
    expect(got).toEqual([{ insert: "keep x", label: "/compact [what the summary must keep]", detail: "",
                           runs: false, descends: false, hint: true }])
    expect(matchArguments({ current: "", choices: [], usage: "" }, "")).toEqual([])
  })

  test("values and nothing matching: the usage, if any", () => {
    expect(matchArguments({ current: "", choices: values, usage: "/x <name>" }, "zzz").map((c) => c.label))
      .toEqual(["/x <name>"])
  })

  test("taking a value replaces the whole argument", () => {
    const token = argumentAt("/jobs watch o", 13)!
    expect(applyCandidate("/jobs watch o", token, "watch on")).toEqual({ text: "/jobs watch on", cursor: 14 })
  })
})
