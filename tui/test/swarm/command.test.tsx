/**
 * Typing /swarm opens the swarm view (SW2b): the window's own command, so
 * the line never reaches the session as a prompt. The routing and the
 * view are wired inside App itself (a protected-file change,
 * tui/src/app.tsx + tui/src/footer/footer.tsx, built by the operator from
 * .gate/px/orig_app.tsx + new_app.tsx and orig_footer.tsx + new_footer.tsx);
 * this test is red on the unpatched tree (the line goes out as
 * session.prompt and no swarm request is made) and green with the patched
 * copies (the panel shows the goal, swarm.status was requested).
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
  calls.length = 0
})

const SWARM_STATUS = {
  schema_version: 1, goal: "measure convergence",
  packages: [{ id: "p0", name: "package 0", phase: "builder",
               members: [{ id: "m1", name: "a", health: "working", activity: "editing" }] }],
  waiting_for_user: [{ item_id: "w0", package: "p0", question: "confirm package 0?" }],
}
const calls: { method: string; params: Record<string, any> }[] = []
const client = {
  request: async (method: string, params: Record<string, any> = {}) => {
    calls.push({ method, params })
    if (method === "swarm.status") return SWARM_STATUS
    if (method === "swarm.triage.list") return { items: SWARM_STATUS.waiting_for_user }
    return { applied: true, message: "", accepted: true }
  },
}

const type = async (text: string) => {
  w!.keys.typeText(text)
  await Bun.sleep(30)
  await w!.settle()
  await Bun.sleep(10)
  await w!.settle()
}

describe("typing /swarm", () => {
  test("opens the view: the panel shows the goal, no prompt goes out", async () => {
    w = await mountWindow({ client })
    await type("/swarm")
    w.keys.pressEnter()
    await Bun.sleep(50)
    await w.settle()
    const frame = w.t.captureCharFrame()
    expect(frame).toContain("◆ SWARM")
    expect(frame).toContain("measure convergence")
    expect(frame).toContain("confirm package 0?")
    const prompts = calls.filter((c) => c.method === "session.prompt")
    expect(prompts).toEqual([])
    expect(calls.some((c) => c.method === "swarm.status")).toBe(true)
  })

  test("with an argument the line goes to the session as typed", async () => {
    w = await mountWindow({ client })
    await type("/swarm status of the conformation swarm")
    w.keys.pressEnter()
    await Bun.sleep(30)
    await w.settle()
    expect(w.t.captureCharFrame()).not.toContain("◆ SWARM")
    expect(calls.some((c) => c.method === "session.prompt")).toBe(true)
  })

  test("/swarm is completed by the window's own commands", async () => {
    w = await mountWindow({ client })
    await type("/sw")
    const frame = w.t.captureCharFrame()
    expect(frame).toContain("/swarm")
  })
})
