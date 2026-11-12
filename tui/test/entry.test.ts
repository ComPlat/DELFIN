/**
 * How the window enters its first session (entry.ts): a new one, the
 * folder's latest (--continue), one by id (--resume <id>), or the resume
 * picker first (a bare --resume); and the picker's rows -- the running
 * sessions on top, then sleeping and saved ones, newest first, each with
 * title, folder, age and state.
 */
import { describe, expect, test } from "bun:test"
import { RESUME_PICK, START_HINT, ageOf, entryOf, neverAsked, resumeRowsOf, settingsFor, shortFolder } from "../src/entry"

describe("the entry", () => {
  test("no flag is a new session, as before", () => {
    expect(entryOf({ cwd: "/w" }, {})).toEqual({ kind: "new" })
    expect(settingsFor({ cwd: "/w", model: "m" }, { kind: "new" })).toEqual({ cwd: "/w", model: "m" })
  })

  test("--continue and --resume <id>", () => {
    expect(entryOf({ continue: true }, {})).toEqual({ kind: "continue" })
    expect(entryOf({ resume: "abc123" }, {})).toEqual({ kind: "resume", id: "abc123" })
    // An id wins over --continue (the service does the same).
    expect(entryOf({ resume: "abc123", continue: true }, {})).toEqual({ kind: "resume", id: "abc123" })
    expect(settingsFor({ resume: "abc123", cwd: "/w" }, { kind: "resume", id: "abc123" }))
      .toEqual({ resume: "abc123", cwd: "/w" })
  })

  test("a bare --resume opens the picker, and the service never sees its marker", () => {
    for (const start of [{ pick: true }, { resume: RESUME_PICK }]) {
      const entry = entryOf({ cwd: "/w", ...start }, {})
      expect(entry).toEqual({ kind: "pick" })
      expect(settingsFor({ cwd: "/w", ...start }, entry)).toEqual({ cwd: "/w" })
    }
  })

  test("DELFIN_TUI_ENTRY says the same; the flags win over it", () => {
    expect(entryOf({}, { DELFIN_TUI_ENTRY: "continue" })).toEqual({ kind: "continue" })
    expect(entryOf({}, { DELFIN_TUI_ENTRY: "resume" })).toEqual({ kind: "pick" })
    expect(entryOf({}, { DELFIN_TUI_ENTRY: "nonsense" })).toEqual({ kind: "new" })
    expect(entryOf({ resume: "x" }, { DELFIN_TUI_ENTRY: "continue" })).toEqual({ kind: "resume", id: "x" })
    // Asked for by the environment only: the service is told to continue.
    expect(settingsFor({ cwd: "/w" }, { kind: "continue" })).toEqual({ cwd: "/w", continue: true })
  })

  test("the first screen's hint names /resume and --continue in one line", () => {
    expect(START_HINT).toContain("/resume")
    expect(START_HINT).toContain("--continue")
    expect(START_HINT).not.toContain("\n")
  })
})

describe("the resume picker's rows", () => {
  const now = 1_000_000
  const home = "/home/u"
  const listed = [
    { sessionId: "run-old-0000", title: "old run", workspace: "/home/u/a", state: "idle", lastActivity: now - 7200 },
    { sessionId: "sleep-00000", title: "asleep", workspace: "/home/u/a", state: "sleeping", lastActivity: now - 60 },
    { sessionId: "run-new-0000", title: "", firstPrompt: "fix the parser\nplease", workspace: "/home/u/b",
      state: "running", lastActivity: now - 30 },
    { sessionId: "asks-000000", title: "asks", workspace: "/x", state: "idle", pending: 1, lastActivity: now - 600 },
  ]
  const past = [
    { sessionId: "saved-old-00", title: "older save", workspace: "/home/u/a", updatedAt: now - 3 * 86400, messages: 4 },
    { sessionId: "saved-new-00", title: "newer save", workspace: "/home/u/a", updatedAt: now - 300, messages: 1 },
    // Listed by the service already: once only.
    { sessionId: "sleep-00000", title: "asleep", workspace: "/home/u/a", updatedAt: now - 60, running: true },
  ]

  test("running ones on top, then sleeping and saved -- each group newest first", () => {
    const rows = resumeRowsOf(listed, past, "run-old-0000", now, home)
    expect(rows.map((r) => r.id)).toEqual(
      ["run-new-0000", "asks-000000", "run-old-0000", "sleep-00000", "saved-new-00", "saved-old-00"])
    expect(rows.map((r) => r.state)).toEqual(["running", "waiting", "running", "sleeping", "saved", "saved"])
    expect(rows.find((r) => r.current)?.id).toBe("run-old-0000")
  })

  test("each row: title, state, short folder, age, short id; whole id and folder for the filter", () => {
    const rows = resumeRowsOf(listed, past, "", now, home)
    const fresh = rows[0]!
    expect(fresh.title).toBe("fix the parser")                 // its first prompt's first line
    expect(fresh.detail).toBe("running · ~/b · now · run-new-")
    expect(fresh.search).toContain("run-new-0000")
    expect(fresh.search).toContain("/home/u/b")
    expect(rows[1]!.detail).toBe("waits for you · /x · 10m ago · asks-000")
    expect(rows.at(-1)!.detail).toBe("saved · ~/a · 3d ago · 4 messages · saved-ol")
    expect(rows.at(-2)!.detail).toBe("saved · ~/a · 5m ago · 1 message · saved-ne")
  })

  test("untitled rows and ages", () => {
    // The one in view is listed even when nothing was asked of it yet.
    const rows = resumeRowsOf([{ sessionId: "n1", workspace: "" }], [], "n1", now, home)
    expect(rows[0]!.title).toBe("(untitled)")
    expect(ageOf(now - 59, now)).toBe("now")
    expect(ageOf(now - 3 * 3600 - 5, now)).toBe("3h ago")
    expect(ageOf(0, now)).toBe("")
  })

  test("a session nobody asked anything, nor named, is not offered (unless in view)", () => {
    const empty = { sessionId: "empty-00000", workspace: "/home/u/a", state: "idle", lastActivity: now - 5 }
    const rows = resumeRowsOf([...listed, empty], past, "", now, home)
    expect(rows.map((r) => r.id)).not.toContain("empty-00000")
    expect(neverAsked(empty)).toBe(true)
    // Asked once (here or from another window), or named: offered.
    expect(neverAsked({ ...empty, lastInput: now - 60 })).toBe(false)
    expect(neverAsked({ ...empty, firstPrompt: "hi" })).toBe(false)
    expect(neverAsked({ ...empty, title: "benzene scan" })).toBe(false)
    expect(resumeRowsOf([empty], [], "empty-00000", now, home).map((r) => r.id)).toEqual(["empty-00000"])
  })

  test("a long folder is cut in the middle, its own name kept", () => {
    const cut = shortFolder("/home/u/projects/chemistry/very/deep/folder/benzene-scan", 24, home)
    expect([...cut].length).toBe(24)
    expect(cut.startsWith("~/")).toBe(true)
    expect(cut.endsWith("benzene-scan")).toBe(true)
    expect(shortFolder("/home/u", 24, home)).toBe("~")
  })
})
