import { describe, expect, test } from "bun:test"
import { leavingNote, resumeCommand } from "../../src/leaving"

const running = { id: "abc123def456", title: "relax Fe", state: "idle" as const, ended: null }

describe("what leaving the window prints", () => {
  test("the command that comes back to exactly this session", () => {
    expect(resumeCommand("abc123def456")).toBe("delfin-agent tui --resume abc123def456")
    const note = leavingNote({ session: running })
    expect(note).toContain("relax Fe keeps running in the agent service")
    expect(note).toContain("resume with: delfin-agent tui --resume abc123def456\n")
    expect(note.endsWith("\n")).toBe(true)
  })

  test("an ended session is saved: the command resumes it as well", () => {
    const note = leavingNote({ session: { ...running, state: "ended", ended: { reason: "exit" } } })
    expect(note).toContain("relax Fe has ended (exit)")
    expect(note).toContain("resume with: delfin-agent tui --resume abc123def456")
  })

  test("an untitled session is named by its id's head", () => {
    expect(leavingNote({ session: { ...running, title: "" } })).toContain("session abc123de keeps running")
  })

  test("the service gone: why, and still the way back", () => {
    const note = leavingNote({ session: running, problem: "the agent service stopped" })
    expect(note).toContain("DELFIN: the agent service stopped")
    expect(note).toContain("resume with: delfin-agent tui --resume abc123def456")
  })

  test("no session: nothing but the problem, if any", () => {
    expect(leavingNote({})).toBe("")
    expect(leavingNote({ problem: "x" })).toBe("DELFIN: x\n")
  })

  test("text from the session cannot move the cursor; an id that is no id is not printed", () => {
    const note = leavingNote({ session: { ...running, title: "a\x1b[2Jb\nc" } })
    expect(note).not.toContain("\x1b")
    expect(note.split("\n").length).toBe(3)
    expect(leavingNote({ session: { ...running, id: "x; rm -rf ~" } })).not.toContain("rm -rf")
  })
})
