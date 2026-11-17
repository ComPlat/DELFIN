import { describe, expect, test } from "bun:test"
import { readFileSync } from "node:fs"
import { join } from "node:path"
import { hostedOf, keptOf, leaveWords, leavingNote, resumeCommand, startedHere } from "../../src/leaving"

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

  test("a session nobody asked anything was not kept: no way back is printed", () => {
    const note = leavingNote({ session: { ...running, discarded: true } })
    expect(note).toBe("DELFIN: nothing was asked in this session, so it was not kept.\n")
    expect(note).not.toContain("resume with")
    expect(leavingNote({ session: { ...running, discarded: true }, problem: "gone" })).toBe("DELFIN: gone\n")
  })

  test("Keep session off: the session ends with the window, and says so", () => {
    const note = leavingNote({ session: { ...running, kept: false } })
    expect(note).toContain("relax Fe ends with this window (Keep session is off)")
    expect(note).not.toContain("keeps running")
    expect(note).toContain("resume with: delfin-agent tui --resume abc123def456")
    expect(leavingNote({ session: { ...running, kept: true } })).toContain("keeps running")
  })

  test("what ctrl+c and the menu say follows the service's keep flag", () => {
    expect(keptOf({ keep: false })).toBe(false)
    expect(keptOf({ keep: true })).toBe(true)
    // A service before the switch kept every session.
    expect(keptOf({})).toBe(true)
    expect(keptOf(undefined)).toBe(true)
    expect(leaveWords(true)).toBe("the session keeps running in the agent service")
    expect(leaveWords(false)).toBe("the session ends -- /keep on keeps it")
  })

  test("a session a dashboard runs goes on there, whatever this window does", () => {
    // Hosted by a dashboard's kernel (summary host "dashboard"): leaving
    // this window lets go of it; the dashboard still has it open.
    expect(hostedOf({ host: "dashboard" })).toBe(true)
    expect(hostedOf({})).toBe(false)
    expect(hostedOf(undefined)).toBe(false)
    expect(leaveWords(false, true)).toBe("the session goes on in the dashboard that runs it")
    expect(leaveWords(true, true)).toBe("the session goes on in the dashboard that runs it")
    const note = leavingNote({ session: { ...running, kept: false, hosted: true } })
    expect(note).toContain("relax Fe goes on in the dashboard that runs it")
    expect(note).not.toContain("ends with this window")
    expect(note).toContain("resume with: delfin-agent tui --resume abc123def456")
    // The leaving paths ask the summary whether a dashboard runs it.
    const app = readFileSync(join(import.meta.dir, "../../src/app.tsx"), "utf8")
    for (const call of app.split("leaveWords(").slice(1)) {
      expect(call.startsWith("keptOf(session()?.summary), hostedOf(session()?.summary))")).toBe(true)
    }
    const main = readFileSync(join(import.meta.dir, "../../src/main.tsx"), "utf8")
    expect(main).toContain("hosted: hostedOf(s.summary)")
  })

  test("a session this window starts is not kept: ctrl+c twice ends it", () => {
    expect(startedHere({ cwd: "/w", title: "t" })).toEqual({ cwd: "/w", title: "t", keep: false })
    // `delfin-agent tui --keep`: the window's sessions start kept.
    expect(startedHere({ cwd: "/w", keep: true })).toEqual({ cwd: "/w", keep: true })
    // Every session the window creates goes through it -- new, resumed from
    // the store, about an attention item.
    const src = readFileSync(join(import.meta.dir, "../../src/main.tsx"), "utf8")
    const creates = src.split('request("session.create",').slice(1)
    expect(creates.length).toBeGreaterThanOrEqual(3)
    for (const call of creates) expect(call.trimStart().startsWith("startedHere(")).toBe(true)
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
