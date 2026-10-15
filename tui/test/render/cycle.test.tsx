/**
 * shift+tab's ladder in the window: default → accept edits → plan →
 * bypass → default, as in the REPL. Bypass is stepped over unless the
 * session says the ladder may stop there (the managed policy allows it,
 * isolation holds each command); a session that never said so is stepped
 * over too. The first step onto it in a session waits for y; any other
 * key -- shift+tab again included -- leaves the posture where it was.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { POSTURES, nextPosture } from "../../src/app"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

async function at(approval: string, bypass: { blocked?: string | null; confirmed?: boolean } = {},
                  answerPauseMs = 0) {
  const w = open = await mountWindow({ height: 30, answerPauseMs })
  w.store.update((d) => {
    const s = d.sessions[w.sid]!
    s.approval = approval
    s.state = "idle"
    s.bypassBlocked = bypass.blocked === undefined ? "" : bypass.blocked
    s.bypassConfirmed = Boolean(bypass.confirmed)
  })
  await w.settle()
  const modes = () => w.calls.filter((c) => c.method === "session.setMode").map((c) => c.params)
  const tab = async () => {
    w.keys.pressTab({ shift: true })
    await w.settle()
  }
  return { w, modes, tab }
}

describe("the ladder", () => {
  test("its order, with and without bypass", () => {
    expect([...POSTURES]).toEqual(["default", "acceptEdits", "plan", "bypassPermissions"])
    expect(nextPosture("default")).toBe("acceptEdits")
    expect(nextPosture("acceptEdits")).toBe("plan")
    expect(nextPosture("plan")).toBe("default")
    expect(nextPosture("plan", true)).toBe("bypassPermissions")
    expect(nextPosture("bypassPermissions", true)).toBe("default")
    expect(nextPosture("diff_approval", true)).toBe("default")
  })
})

describe("shift+tab in the window", () => {
  test("default steps to accept edits, accept edits to plan", async () => {
    const { modes, tab } = await at("default")
    await tab()
    expect(modes()).toEqual([{ sessionId: "s1", mode: "acceptEdits" }])
  })

  test("from plan, bypass waits for y: nothing is switched until then", async () => {
    const { w, modes, tab } = await at("plan")
    await tab()
    expect(modes()).toEqual([])
    expect(w.footer().join("\n")).toContain("y switches, any other key stays in plan")
    await w.keys.typeText("y")
    await w.settle()
    expect(modes()).toEqual([{ sessionId: "s1", mode: "bypassPermissions", confirmed: true }])
    expect(w.footer().join("\n")).not.toContain("y switches")
  })

  test("an unanswered offer stays on the previous posture; shift+tab again moves past bypass", async () => {
    const { w, modes, tab } = await at("plan")
    await tab()
    await w.keys.typeText("x")
    await w.settle()
    expect(modes()).toEqual([])
    expect(w.footer().join("\n")).not.toContain("y switches")
    expect(w.footer().join("\n")).toContain("❯ x")         // the key did what it does
    await tab()
    expect(modes()).toEqual([])
    await tab()                                            // the offer open: on past bypass
    expect(modes()).toEqual([{ sessionId: "s1", mode: "default" }])
  })

  test("a y typed sooner than the answer pause is typing, not an answer", async () => {
    const { w, modes, tab } = await at("plan", {}, 10_000)
    await tab()
    await w.keys.typeText("y")
    await w.settle()
    expect(modes()).toEqual([])
  })

  test("bypass the policy or the missing cage rules out is stepped over", async () => {
    for (const blocked of ["the administrator's managed policy disallows it", null]) {
      const { modes, tab } = await at("plan", { blocked })
      await tab()
      expect(modes()).toEqual([{ sessionId: "s1", mode: "default" }])
      await open?.close()
      open = null
    }
  })

  test("once confirmed in the session, the ladder passes through bypass without asking", async () => {
    const { w, modes, tab } = await at("plan", { confirmed: true })
    await tab()
    expect(modes()).toEqual([{ sessionId: "s1", mode: "bypassPermissions" }])
    expect(w.footer().join("\n")).not.toContain("y switches")
  })
})
