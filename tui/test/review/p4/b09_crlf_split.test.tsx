/**
 * B09: an answer with Windows line ends whose "\r\n" is split between two
 * deltas. AnswerBlocks normalises each delta on its own ("\r" -> "\n",
 * then "\n"), so the streamed text has extra newlines, text.ended's whole
 * text does not start with what was committed, and the presenter commits
 * the answer AGAIN after a false "arrived while the window was
 * reconnecting" note.
 */
import { afterEach, expect, test } from "bun:test"
import { mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

test("CRLF split across deltas", async () => {
  w = await mount()
  const full = "First paragraph of the answer.\r\n\r\nSecond paragraph.\r\n\r\nThird."
  const deltas = ["First paragraph of the answer.\r", "\n\r\nSecond paragraph.\r", "\n\r\nThird."]
  w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
  w.store.push({ type: "text.started", session: "s1", ts: 0, turnId: "t1", partId: "p1" } as never)
  for (const d of deltas) w.store.push({ type: "text.delta", session: "s1", ts: 0, partId: "p1", delta: d } as never)
  w.store.push({ type: "text.ended", session: "s1", ts: 0, seq: 2, turnId: "t1", partId: "p1", text: full } as never)
  await w.settle()
  const text = w.t.externalOutput.takeText()
  console.log(text.split("\n").filter((l) => l.trim()).join("\n"))
  expect(text.split("First paragraph").length - 1).toBe(1)
  expect(text).not.toContain("reconnecting")
})
