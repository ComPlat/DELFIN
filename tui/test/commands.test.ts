import { describe, expect, test } from "bun:test"
import { WINDOW_COMMANDS, mergeOf, nearestCommand, unknownCommandNote, unknownCommandOf, worktreeRetryOf } from "../src/commands"

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

describe("an unknown /command", () => {
  const known = ["/help", "/keys", "/status", "/model", "/resume", "/rename", "/compact", "/plugin:skill"]
  const nothingThere = () => false

  test("is refused, with the command it was probably meant to be", () => {
    expect(unknownCommandOf("/nonexistent", known, nothingThere)).toEqual({ head: "/nonexistent", near: "" })
    expect(unknownCommandOf("/stauts", known, nothingThere)).toEqual({ head: "/stauts", near: "/status" })
    expect(unknownCommandOf("/hepl me", known, nothingThere)).toEqual({ head: "/hepl", near: "/help" })
    expect(unknownCommandOf("/compa", known, nothingThere)).toEqual({ head: "/compa", near: "/compact" })
  })

  test("two commands run together are one unknown command", () => {
    expect(unknownCommandOf("/help/keys", known, nothingThere)).toEqual({ head: "/help/keys", near: "/help" })
  })

  test("a known command, any case, with or without arguments, is no unknown one", () => {
    for (const line of ["/help", "/HELP", "/model kit.glm-5", "  /status  ", "/plugin:skill go"]) {
      expect(unknownCommandOf(line, known, nothingThere)).toBeNull()
    }
  })

  test("a path is no command: it goes out as typed", () => {
    for (const line of ["/home/me/run.py fails, why?", "/scratch/calc/out.log", "/.hidden", "//", "/*", "/123 apples"]) {
      expect(unknownCommandOf(line, known, nothingThere)).toBeNull()
    }
    // A single word that names something on disk is a path too.
    expect(unknownCommandOf("/tmp is full", known, (p) => p === "/tmp")).toBeNull()
    expect(unknownCommandOf("/tmpx", known, (p) => p === "/tmp")).not.toBeNull()
  })

  test("no slash, or only a slash, is not judged", () => {
    expect(unknownCommandOf("hello /nonexistent", known, nothingThere)).toBeNull()
    expect(unknownCommandOf("/", known, nothingThere)).toBeNull()
    expect(unknownCommandOf("", known, nothingThere)).toBeNull()
  })

  test("what the window says names the line, the guess and /help", () => {
    expect(unknownCommandNote({ head: "/stauts", near: "/status" }))
      .toBe("Unknown command: /stauts — did you mean /status? /help lists them. Nothing was sent.")
    expect(unknownCommandNote({ head: "/nonexistent", near: "" }))
      .toBe("Unknown command: /nonexistent /help lists them. Nothing was sent.")
  })

  test("nearestCommand: close names only", () => {
    expect(nearestCommand("/sttus", known)).toBe("/status")
    expect(nearestCommand("/xyzzy", known)).toBe("")
    expect(nearestCommand("/re", known)).toBe("")              // too short to guess between /resume and /rename
  })
})
