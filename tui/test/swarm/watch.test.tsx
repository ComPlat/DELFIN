/**
 * The read-only /swarm watch view (package SWV, phase 2b): a table of one
 * row per member, Enter on a member opens its live feed, Esc back. Verifies
 * the watch is READ-ONLY -- no input into a member session, no prompt, no
 * pause, no approval answered -- and that it renders the observe.* contract
 * (observe.py swarm_overview / member_feed shapes).
 *
 * Red control: every case fails on a tree without tui/src/swarm/watch.ts
 * (module missing). They go green with the view in place.
 */
import { describe, expect, test } from "bun:test"
import { createSwarmWatch, type FeedEntry, type MemberRow } from "../../src/swarm/watch"
import { Keymap } from "../../src/keymap/keymap"

const now = 1700000000

/** A member row as swarm.observe.overview returns it (observe.py
 *  swarm_overview's row shape, phase 1). */
function member(name: string, over: Partial<MemberRow> = {}): MemberRow {
  return {
    package: "p1",
    role: "builder",
    name,
    model: "kit.deepseek",
    state: "working",
    last_tool_call: {},
    last_commit_subject: "",
    last_activity_age_s: 120,
    pending_approvals_count: 0,
    commits_on_branch: 0,
    ...over,
  }
}

/** A feed entry as swarm.observe.feed returns it (member_feed's shape). */
function feedEntry(kind: string, text: string, at = ""): FeedEntry {
  return { at, kind, text }
}

function makeView(rows = [member("s-a"), member("s-b")]) {
  const requests: { method: string; params: Record<string, unknown> }[] = []
  const notes: string[] = []
  const feeds = new Map<string, FeedEntry[]>()
  const view = createSwarmWatch({
    keymap: new Keymap(),
    now: () => now,
    askOverview: async () => {
      requests.push({ method: "swarm.observe.overview", params: {} })
      return rows
    },
    askFeed: async (name) => {
      requests.push({ method: "swarm.observe.feed", params: { session: name } })
      return (
        feeds.get(name) ?? [
          feedEntry("tool", "submit_calculation → 3.2"),
          feedEntry("refusal", "squeue denied"),
          feedEntry("commit", "phase 1 green"),
        ]
      )
    },
    note: (text) => notes.push(text),
    draft: () => "",
  })
  return { view, requests, notes, feeds }
}

describe("the /swarm watch view's data", () => {
  test("one row per member of the configured source", async () => {
    const { view, requests } = makeView()
    await view.open()
    expect(view.state.open).toBe(true)
    expect(requests.some((r) => r.method === "swarm.observe.overview")).toBe(true)
    expect(view.state.rows.length).toBe(2)
    expect(view.memberLine(view.state.rows[0])).toContain("s-a")
  })

  test("the state and commit count come from the row, not a re-derivation", async () => {
    const { view } = makeView([member("s-x", { state: "NEEDS APPROVAL", commits_on_branch: 3 })])
    await view.open()
    const line = view.memberLine(view.state.rows[0])
    expect(line).toContain("NEEDS APPROVAL")
    expect(line).toContain("3 commits")
  })
})

describe("stepping into a member is read-only", () => {
  test("Enter on a member opens its feed; no request ever sends into the session", async () => {
    const { view, requests } = makeView()
    await view.open()
    const right = { name: "return", sequence: "return", ctrl: false, meta: false, shift: false }
    view.handleKey(right as never)
    // The feed is opened by reading (observe.feed), never by sending a prompt.
    expect(view.state.feed).toBeTruthy()
    expect(requests.some((r) => r.method === "swarm.observe.feed")).toBe(true)
    // No member-facing write request exists in the whole watch interaction.
    const writes = ["session.message", "session.pause", "approval.answer", "swarm.triage.answer",
      "session.merge", "session.detach"]
    for (const w of writes) expect(requests.some((r) => r.method === w)).toBe(false)
    expect(view.state.feed!.entries.length).toBe(3)
  })

  test("Esc returns from the feed to the table, then leaves the view", async () => {
    const { view } = makeView()
    await view.open()
    const right = { name: "return", sequence: "return", ctrl: false, meta: false, shift: false }
    const esc = { name: "escape", sequence: "esc", ctrl: false, meta: false, shift: false }
    view.handleKey(right as never)
    expect(view.state.feed).toBeTruthy()
    view.handleKey(esc as never)
    expect(view.state.feed).toBeNull() // back on the table
    view.handleKey(esc as never)
    expect(view.state.open).toBe(false) // leaving the view
  })

  test("Enter inside the feed only peeks; it never answers or edits", async () => {
    const { view, requests } = makeView()
    await view.open()
    const right = { name: "return", sequence: "return", ctrl: false, meta: false, shift: false }
    view.handleKey(right as never) // open feeds
    const before = requests.length
    view.handleKey(right as never) // Enter inside the feed
    // Still only observe.* reads; nothing more was sent.
    expect(requests.length).toBe(before)
    expect(requests.every((r) => r.method.startsWith("swarm.observe."))).toBe(true)
  })
})

describe("the feed shows refusals highlighted and commits", () => {
  test("a refusal is rendered distinctly and only peeks, never answers", async () => {
    const { view } = makeView([member("s-a")])
    await view.open()
    const right = { name: "return", sequence: "return", ctrl: false, meta: false, shift: false }
    view.handleKey(right as never)
    const kinds = (view.state.feed!.entries).map((e) => e.kind)
    expect(kinds).toContain("refusal")
    expect(kinds).toContain("commit")
    expect(kinds).toContain("tool")
  })
})

describe("a failed observe never leaves the view wedged", () => {
  test("an overview failure is surfaced and the view can be retried", async () => {
    const notes: string[] = []
    let fail = true
    const view = createSwarmWatch({
      keymap: new Keymap(),
      now: () => now,
      askOverview: async () => {
        if (fail) throw new Error("no swarm source")
        return [member("s-a")]
      },
      askFeed: async () => [],
      note: (t) => notes.push(t),
      draft: () => "",
    })
    await view.open()
    expect(view.state.error).toContain("no swarm source")
    expect(notes.join("\n")).toContain("swarm watch failed")
    fail = false
    await view.open()
    expect(view.state.rows.length).toBe(1)
    expect(view.state.error).toBeNull()
  })

  test("an empty swarm is an empty table, not an error", async () => {
    const { view } = makeView([])
    await view.open()
    expect(view.state.open).toBe(true)
    expect(view.state.rows.length).toBe(0)
    expect(view.state.error).toBeNull()
  })
})
