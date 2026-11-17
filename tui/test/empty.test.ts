/**
 * Which sessions of the window were never asked anything (empty.ts): a new
 * one with no message and no turn is ended when the window lets it go; a
 * resumed one, or one that was asked something, never.
 */
import { describe, expect, test } from "bun:test"
import { EmptySessions } from "../src/empty"

describe("empty sessions", () => {
  test("a new session is empty until a message or a turn reaches it", () => {
    const e = new EmptySessions()
    e.started("s1")
    expect(e.isEmpty("s1")).toBe(true)
    e.saw({ type: "user.message", session: "s1" })
    expect(e.isEmpty("s1")).toBe(false)
    e.started("s2")
    e.saw({ type: "turn.started", session: "s2" })
    expect(e.isEmpty("s2")).toBe(false)
  })

  test("other events leave it empty; another session's message does not count", () => {
    const e = new EmptySessions()
    e.started("s1")
    e.saw({ type: "session.info", session: "s1" })
    e.saw({ type: "notice", session: "s1" })
    e.saw({ type: "user.message", session: "s2" })
    expect(e.isEmpty("s1")).toBe(true)
  })

  test("what its summary says counts: another window may have asked it", () => {
    const e = new EmptySessions()
    e.started("s1")
    expect(e.isEmpty("s1", { firstPrompt: null, lastInput: null })).toBe(true)
    expect(e.isEmpty("s1", { firstPrompt: "relax Fe" })).toBe(false)
    expect(e.isEmpty("s1", { lastInput: 1700000000 })).toBe(false)
  })

  test("a session this window did not start new (resumed, attached) is never empty", () => {
    const e = new EmptySessions()
    expect(e.isEmpty("old")).toBe(false)
    expect(e.isEmpty("")).toBe(false)
  })

  test("asked before it was known as new: not empty", () => {
    const e = new EmptySessions()
    e.saw({ type: "user.message", session: "s1" })
    e.started("s1")
    expect(e.isEmpty("s1")).toBe(false)
  })

  test("ended for being empty: remembered for the leaving note, never ended twice", () => {
    const e = new EmptySessions()
    e.started("s1")
    e.ended("s1")
    expect(e.wasEnded("s1")).toBe(true)
    expect(e.isEmpty("s1")).toBe(false)
    expect(e.wasEnded("s2")).toBe(false)
  })
})
