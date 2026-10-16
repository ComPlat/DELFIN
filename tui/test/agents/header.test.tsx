/**
 * The agent view's head says its keys at its right: what opens a session,
 * what stops one (twice), what goes back -- from the keymap, so a rebound
 * key shows; as many as fit beside the counts, the most needed first.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { headerKeys } from "../../src/agents/view"
import { Keymap } from "../../src/keymap/keymap"
import { width } from "../../src/text"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
})

const now = Date.now() / 1000

async function openView(cols: number) {
  w = await mountWindow({ width: cols, height: 24 })
  w.store.push({ type: "session.summary", session: "s1", ts: 0,
                 summary: { sessionId: "s1", title: "this window", workspace: "/w", state: "idle", pending: 0,
                            createdAt: 900, lastActivity: now } } as never)
  await w.settle()
  w.keys.pressArrow("left")
  await w.settle()
  return w.t.captureCharFrame().split("\n").find((l) => l.includes("◆ AGENTS")) ?? ""
}

describe("the agent view's head", () => {
  test("names the keys that open, stop and go back", async () => {
    const head = await openView(110)
    expect(head).toContain("◆ AGENTS  1 session")
    expect(head).toMatch(/→\/enter open · ctrl\+x ×2 stop · esc back\s*$/)
  })

  test("on a narrow terminal keeps what fits, opening and stopping first", async () => {
    const head = await openView(60)
    expect(head).toContain("◆ AGENTS  1 session")
    expect(head).toContain("→/enter open · ctrl+x ×2 stop")
    expect(width(head.trimEnd())).toBeLessThanOrEqual(60)
  })

  test("says a rebound key", () => {
    const km = new Keymap({ agentsStop: ["ctrl+k"], agentsOpen: ["return"] })
    expect(headerKeys(km, 200)).toBe("enter open · ctrl+k ×2 stop · esc back")
    expect(headerKeys(km, 30)).toBe("enter open · ctrl+k ×2 stop")
    expect(headerKeys(km, 5)).toBe("")
  })
})
