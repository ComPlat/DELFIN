/**
 * The /swarm watch table on the SCREEN (package SWV, phase 2b): with the
 * watch open, the members' rows -- name, state, commits -- are on the
 * window's frame, and Esc leaves it again. Rendered headless through the
 * real App against a fake observe client (the same harness the agents
 * view's screen test uses).
 *
 * Red control: on the tree without a renderer for the watch the table is
 * never drawn -- the footer mounts only the agents and swarm panels
 * (app.tsx passes agents={agents} swarm={swarm}) -- so the frame shows
 * nothing of the members.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { Keymap } from "../../src/keymap/keymap"
import { mountWindow, type Call } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

/** The observe answers the real RPC rows of a scratch swarm: one running,
 *  one stopped, one with a PR. */
const members = [
  { package: "rr4", role: "run", name: "rr4-run", model: "agent rr4-run", state: "working",
    commits_on_branch: 0, last_tool_call: { tool: "edit_file", arg: "delfin/pipeline.py" },
    last_commit_subject: "", last_activity_age_s: 0, pending_approvals_count: 0 },
  { package: "rr4", role: "stop", name: "rr4-stop", model: "agent rr4-stop", state: "ended",
    commits_on_branch: 0, last_tool_call: { tool: "bash", arg: "git status" },
    last_commit_subject: "", last_activity_age_s: 3600, pending_approvals_count: 0 },
  { package: "rr4", role: "pr", name: "rr4-pr", model: "agent rr4-pr", state: "NEEDS APPROVAL",
    commits_on_branch: 1, last_tool_call: { tool: "git", arg: "commit -m 'rr4: fix'" },
    last_commit_subject: "rr4: fix", last_activity_age_s: 120, pending_approvals_count: 1 },
]

function client() {
  const handle = (call: Call) => {
    if (call.method === "swarm.observe.overview") return { rows: members }
    if (call.method === "swarm.observe.feed") return { entries: [
      { at: "2026-06-01T10:00:05Z", kind: "tool", text: "pytest -q" }] }
    return {}
  }
  return { request: async (method: string, params: Record<string, unknown> = {}) =>
             handle({ method, params } as Call) } as never
}

async function openWatch() {
  w = await mountWindow({ client: client(), width: 110, height: 34 })
  await w.keys.typeText("/swarm watch scratch/realrun/roster.json")
  await Bun.sleep(30)
  await w.settle()
  w.keys.pressEnter()
  await w.settle()
  await Bun.sleep(30)
  await w.settle()
  return w
}

describe("the /swarm watch on the screen", () => {
  test("the members' rows are on the frame with state and commits", async () => {
    const win = await openWatch()
    const frame = win.screenLines().join("\n")
    expect(frame).toContain("rr4-run")
    expect(frame).toContain("working")
    expect(frame).toContain("rr4-stop")
    expect(frame).toContain("rr4-pr")
    expect(frame).toContain("NEEDS APPROVAL")
    expect(frame).toContain("1 commits")
  })

  test("Esc leaves the watch: the frame shows no member rows again", async () => {
    const win = await openWatch()
    await win.escape()
    const frame = win.screenLines().join("\n")
    expect(frame).not.toContain("rr4-run")
  })
})
