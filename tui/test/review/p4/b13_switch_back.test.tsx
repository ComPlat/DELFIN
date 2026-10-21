/**
 * B13: switching back to a session this window showed before. switchTo
 * (main.tsx) prints a divider and attaches; the attach's replay is pushed
 * through the store, whose seen-seq filter drops every event the window
 * saw before -- so after "now: <A>" only A's banner comes, none of its
 * turns. The "divider + replay" promise holds only for a session never
 * shown in this window. Driven here exactly as switchTo drives it.
 */
import { afterEach, expect, test } from "bun:test"
import { ensureSession } from "../../../src/store/state"
import { withReplay } from "../../../src/store/store"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const info = (sid: string, title: string) => ({ type: "session.info", session: sid, ts: 0, seq: 0, sessionId: sid,
  workspace: "/w", provider: "kit", model: "m", mode: "solo", approval: "default", title, banner: [], notices: [] })
const turn = (sid: string, base: number, q: string, a: string) => [
  { type: "user.message", session: sid, ts: 0, seq: base + 1, text: q, source: "user" },
  { type: "turn.started", session: sid, ts: 0, seq: base + 2, turnId: `t${base}`, source: "user" },
  { type: "text.ended", session: sid, ts: 0, seq: base + 3, turnId: `t${base}`, partId: `p${base}`, text: a },
  { type: "turn.ended", session: sid, ts: 0, seq: base + 4, turnId: `t${base}`, durationMs: 100, inputTokens: 1, outputTokens: 1, costUsd: 0, toolCalls: 0 },
]
const A = [info("A", "session-A"), ...turn("A", 0, "question one to A", "Answer one from A.")]
const B = [info("B", "session-B"), ...turn("B", 0, "question to B", "Answer from B.")]

async function switchTo(id: string, events: any[]) {
  // main.tsx switchTo, minus the requests: unsee -> attach -> register -> show(divider) -> replay
  w!.store.unsee(id)
  await withReplay(w!.store, async () => ({ session: { title: `session-${id}` }, events }), (attached: any) => {
    w!.store.update((d) => { ensureSession(d, id, attached.session) })
    w!.view.show(id, { divider: `now: ${attached.session.title}` })
    return attached.events
  })
  await w!.settle()
}

test("A, then B, then back to A", async () => {
  w = await mount({ show: false })
  await switchTo("A", A)
  await switchTo("B", B)
  w.t.externalOutput.takeText()
  await switchTo("A", A)                         // the same replay A's attach returns
  const after = w.t.externalOutput.takeText()
  console.log(after.split("\n").filter((l) => l.trim()).join("\n"))
  expect(after).toContain("now: session-A")
  expect(after).toContain("question one to A")
})
