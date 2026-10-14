/**
 * The banner and the status line, as the REPL shows them: the banner's
 * safety lines kept and its wordmark box made the window's head; the turn
 * line counting this turn's tokens and cost, estimating output from what
 * streamed until the provider counts it; a wait line between turns.
 */
import { afterEach, expect, test } from "bun:test"
import { turnUsage } from "../../src/footer/footer"
import { bannerOf } from "../../src/scrollback/present"
import { mountWindow } from "../harness"

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

const ev = (type: string, fields: Record<string, unknown>) => ({ type, session: "s1", ts: 0, ...fields }) as never

test("the REPL's banner: its wordmark becomes the head, its safety lines stay", () => {
  const banner = [
    "╭────────╮", "│ DELFIN │", "╰────────╯", "",
    "delfin-agent · kit/kit.deepseek-v4-flash · solo",
    "workspace  ~/project  (git main · 2 uncommitted)",
    "approval   default",
    "isolation  on (bubblewrap) — a command the agent runs writes only inside the workspace (--no-isolate)",
    "session    25016b8f",
    "esc interrupt · shift+tab approval mode · /help · ctrl+d exit",
  ]
  expect(bannerOf({ banner, provider: "kit", model: "kit.deepseek-v4-flash", title: "parser fix" })).toEqual({
    head: "kit/kit.deepseek-v4-flash · solo · parser fix",
    lines: banner.slice(5, -1),
  })
  expect(bannerOf({ banner: [], provider: "fake", model: "m", mode: "solo" }).head).toBe("fake/m · solo")
})

test("this turn's tokens and cost, not the conversation's", () => {
  const turn = { id: "t2", source: "user", startedAt: 0, label: "", streamed: 0,
                 base: { inputTokens: 1000, outputTokens: 200, costUsd: 0.01 } }
  expect(turnUsage(turn, { inputTokens: 1000, outputTokens: 200, costUsd: 0.01 })).toEqual(["↑0 ↓…"])
  expect(turnUsage({ ...turn, streamed: 400 }, { inputTokens: 2500, outputTokens: 200, costUsd: 0.01 }))
    .toEqual(["↑1.5k ↓~100"])
  expect(turnUsage(turn, { inputTokens: 2500, outputTokens: 260, costUsd: 0.0125 })).toEqual(["↑1.5k ↓60", "$0.0025"])
})

test("the footer's turn line counts from the turn's start", async () => {
  const w = open = await mountWindow({ width: 100 })
  w.store.push(ev("usage", { inputTokens: 5000, outputTokens: 900, costUsd: 0 }))   // the session so far
  w.store.push(ev("turn.started", { seq: 1, turnId: "t2", source: "user" }))
  w.store.push(ev("text.started", { turnId: "t2", partId: "p1" }))
  w.store.push(ev("text.delta", { partId: "p1", delta: "x".repeat(80) }))
  await w.settle()
  expect(w.footer().join("\n")).toContain("↑0 ↓~20")
  w.store.push(ev("usage", { inputTokens: 6200, outputTokens: 950, costUsd: 0 }))
  await w.settle()
  expect(w.footer().join("\n")).toContain("↑1.2k ↓50")
})

test("a wait between turns is shown, and goes when the answer comes", async () => {
  const w = open = await mountWindow({ width: 100 })
  w.store.push(ev("wait", { text: "the endpoint is down; trying again in 30 s" }))
  await w.settle()
  expect(w.footer().join("\n")).toContain("… the endpoint is down; trying again in 30 s")
  w.store.push(ev("turn.started", { seq: 1, turnId: "t1", source: "harness" }))
  w.store.push(ev("text.delta", { partId: "p1", delta: "back" }))
  await w.settle()
  expect(w.footer().join("\n")).not.toContain("the endpoint is down")
})

test("a turn ended by a cancel says it stopped", async () => {
  const w = open = await mountWindow({ width: 100 })
  w.store.push(ev("turn.started", { seq: 1, turnId: "t1", source: "user" }))
  w.store.push(ev("turn.ended", { seq: 2, turnId: "t1", durationMs: 20100, inputTokens: 0, outputTokens: 0,
                                  costUsd: 0, toolCalls: 0, stopped: true }))
  await w.settle()
  const text = w.t.externalOutput.takeText()
  expect(text).toContain("◼ stopped · 20.1s")
  expect(text).not.toContain("✓ 20.1s")
})
