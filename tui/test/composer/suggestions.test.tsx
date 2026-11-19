/**
 * The session's suggestions, as the dashboard offers them in its box: the
 * empty prompt shows the first; Tab takes it, Tab again the next. Not
 * while a turn runs or a question is open.
 */
import { afterEach, expect, test } from "bun:test"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const info = (suggestions: string[]) => ({
  type: "session.info", session: "s1", ts: 0, seq: 0, sessionId: "s1", workspace: "/w",
  provider: "fake", model: "m", mode: "solo", approval: "default", banner: [], notices: [], suggestions,
})
const promptLine = () => w!.footer().find((l) => l.startsWith("❯ ")) ?? ""
const tab = async () => {
  w!.keys.pressTab()
  await Bun.sleep(5)
  await w!.settle()
}

test("the first suggestion waits in the empty prompt; tab takes it, tab again the next", async () => {
  w = await mountWindow()
  w.store.push(info(["Run the tests", "Write the summary"]) as never)
  await w.settle()
  expect(w.footer().join("\n")).toContain("Run the tests  (tab takes it, again for the next)")
  await tab()
  expect(promptLine()).toContain("❯ Run the tests")
  await tab()
  expect(promptLine()).toContain("❯ Write the summary")
  await tab()
  expect(promptLine()).toContain("❯ Run the tests")
})

test("a typed line is never replaced, and an update brings new ones", async () => {
  w = await mountWindow()
  w.store.push(info(["Run the tests"]) as never)
  await w.settle()
  await w.keys.typeText("my own")
  await w.settle()
  await tab()
  expect(promptLine()).toContain("❯ my own")
  w.store.push({ type: "session.updated", session: "s1", ts: 0, seq: 1, suggestions: ["Next thing"] } as never)
  await w.settle()
  expect(w.store.state.sessions["s1"]!.suggestions).toEqual(["Next thing"])
})

test("not while a turn runs", async () => {
  w = await mountWindow()
  w.store.push(info(["Run the tests"]) as never)
  w.store.push({ type: "turn.started", session: "s1", ts: 0, seq: 1, turnId: "t1", source: "user" } as never)
  await w.settle()
  expect(w.footer().join("\n")).not.toContain("tab takes it")
  await tab()
  expect(promptLine()).not.toContain("Run the tests")
})
