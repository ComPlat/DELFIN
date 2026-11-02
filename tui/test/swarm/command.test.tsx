/**
 * Typing /swarm opens the swarm view (SW2b): the window's own command, so
 * the line never reaches the session as a prompt. The routing branch is a
 * protected-file change (tui/src/app.tsx) built by the operator from
 * .gate/px/orig_app.tsx + .gate/px/new_app.tsx (plus the harness's
 * onSwarm pass-through); this test is red on orig (the line goes out as
 * session.prompt) and green on new (onSwarm is called, no prompt request).
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
})

const type = async (text: string) => {
  w!.keys.typeText(text)
  await Bun.sleep(30)
  await w!.settle()
  await Bun.sleep(10)
  await w!.settle()
}

describe("typing /swarm", () => {
  test("calls onSwarm and never sends a prompt", async () => {
    let swarms = 0
    w = await mountWindow({ onSwarm: () => { swarms++ } })
    await type("/swarm")
    w.keys.pressEnter()
    await Bun.sleep(30)
    await w.settle()
    expect(swarms).toBe(1)
    const prompts = w.calls.filter((c) => c.method === "session.prompt")
    expect(prompts).toEqual([])
  })

  test("with an argument the line goes to the session as typed", async () => {
    let swarms = 0
    w = await mountWindow({ onSwarm: () => { swarms++ } })
    await type("/swarm status of the conformation swarm")
    w.keys.pressEnter()
    await Bun.sleep(30)
    await w.settle()
    expect(swarms).toBe(0)
    expect(w.calls.some((c) => c.method === "session.prompt")).toBe(true)
  })

  test("/swarm is completed by the window's own commands", async () => {
    w = await mountWindow({ onSwarm: () => {} })
    await type("/sw")
    const frame = w.t.captureCharFrame()
    expect(frame).toContain("/swarm")
  })
})
