import { describe, expect, test } from "bun:test"
import { WINDOW_COMMANDS, mergeOf, worktreeRetryOf } from "../src/commands"

describe("/merge", () => {
  test("is offered with its summary", () => {
    const entry = WINDOW_COMMANDS.find((c: { name: string }) => c.name === "/merge")
    expect(entry?.source).toEqual("window")
    expect(entry?.takesArgs).toBe(true)
    expect(entry?.summary).toMatch(/worktree/i)
  })

  test("a bare /merge asks for no path", () => {
    expect(mergeOf("/merge")).toEqual("")
  })

  test("the argument is the worktree path", () => {
    expect(mergeOf("/merge .delfin/worktrees/notes")).toEqual(".delfin/worktrees/notes")
  })

  test("not a /merge line", () => {
    expect(mergeOf("merge it")).toBeNull()
    expect(mergeOf("/resume x")).toBeNull()
  })
})

describe("worktreeRetryOf", () => {
  test("NEEDS_LOCAL_PROJECT asks to init the project", () => {
    expect(worktreeRetryOf({ code: -32006, message: "plain folder" })).toEqual({ initProject: true })
  })

  test("any other refusal is not retried", () => {
    expect(worktreeRetryOf({ code: -32001, message: "no" })).toBeNull()
    expect(worktreeRetryOf(null)).toBeNull()
  })
})
