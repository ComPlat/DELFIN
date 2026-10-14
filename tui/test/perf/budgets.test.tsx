/**
 * The performance budgets (docs/TUI_PLAN.md, phase 7), headless:
 *
 * - a keystroke reaches the frame in under 16 ms (median), never over 100;
 * - a session of 1,000 rows is written again (a resize, ctrl+l) in under 1 s.
 *
 * Launch to footer (300 ms with the service running) is measured live:
 * test/e2e/check_startup_tmux.py. Timings on a loaded machine vary, so the
 * bounds are what a person notices, not the best case.
 */
import { afterEach, expect, test } from "bun:test"
import { itemWriter } from "../../src/scrollback/rows"
import { DARK } from "../../src/theme/theme"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await w?.close()
  w = null
})

const median = (xs: number[]) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)]!

test("a keystroke reaches the frame in under 16 ms", async () => {
  w = await mountWindow({ width: 100, height: 30 })
  await w.settle()
  const times: number[] = []
  let typed = ""
  for (const ch of "the quick brown fox jumps over the lazy dog") {
    const t0 = performance.now()
    w.keys.pressKey(ch)
    typed += ch
    await w.t.renderOnce()
    if (!w.t.captureCharFrame().includes(`❯ ${typed}`)) {
      // the textarea takes the key a tick later on some paths
      await Bun.sleep(0)
      await w.t.renderOnce()
    }
    times.push(performance.now() - t0)
    expect(w.t.captureCharFrame()).toContain(`❯ ${typed}`)
  }
  console.log(`keystroke to frame: median ${median(times).toFixed(1)} ms, max ${Math.max(...times).toFixed(1)} ms`)
  expect(median(times)).toBeLessThan(16)
  expect(Math.max(...times)).toBeLessThan(100)
})

test("a session of 1,000 rows is written again in under 1 s", async () => {
  w = await mountWindow({ width: 100, height: 30 })
  const sid = w.sid
  let seq = 0
  const push = (type: string, fields: Record<string, unknown>) =>
    w!.store.push({ type, session: sid, ts: 0, seq: ++seq, ...fields } as never)
  // 32 turns of about 33 rows each: a prompt, a tool call, an answer.
  for (let n = 1; n <= 32; n++) {
    push("user.message", { text: `turn ${n}`, source: "user" })
    push("turn.started", { turnId: `t${n}`, source: "user" })
    push("tool.started", { turnId: `t${n}`, callId: `c${n}`, name: "read_file", input: { path: "README.md" } })
    push("tool.ended", { turnId: `t${n}`, callId: `c${n}`, name: "read_file", output: "# demo\n", meta: { ok: true } })
    push("text.ended", { turnId: `t${n}`, partId: `p${n}`,
                         text: Array.from({ length: 14 }, (_, i) => `Line ${i + 1} of the answer to turn ${n}.`).join("\n\n") })
    push("turn.ended", { turnId: `t${n}`, durationMs: 400 })
  }
  await w.settle()
  await w.view.queue.drain()
  const rows = w.t.externalOutput.takeText().split("\n").length
  expect(rows).toBeGreaterThan(1000)
  // What a replay costs: every row written again, in order (the reset of
  // the terminal before it needs a real one; check_resize_tmux.py has it).
  const write = itemWriter(w.t.renderer, () => DARK)
  const t0 = performance.now()
  for (const item of w.view.transcript.items) await write(item)
  const took = performance.now() - t0
  console.log(`writing ${rows} rows again: ${took.toFixed(0)} ms`)
  expect(w.t.externalOutput.takeText()).toContain("Line 14 of the answer to turn 32.")
  expect(took).toBeLessThan(1000)
})
