/**
 * Prompts typed while the window's first session starts are held, said to
 * be, and go out once it is up -- every one, in the order typed; the note
 * that they would goes once they have.
 */
import { afterEach, expect, test } from "bun:test"
import { ensureSession } from "../../src/store/state"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

test("prompts typed while the session starts all go out, in order, and the note goes with them", async () => {
  w = await mountWindow({ show: false })
  w.store.update((d) => { d.starting = true })
  await w.settle()
  for (const n of [1, 2, 3]) {
    await w.keys.typeText(`hello number ${n}`)
    w.keys.pressEnter()
    await w.settle()
  }
  expect(w.footer().join("\n")).toContain("the session is starting")
  expect(w.calls.filter((c) => c.method === "session.prompt")).toEqual([])
  w.store.update((d) => {
    ensureSession(d, "s1", {})
    d.current = "s1"
    d.starting = false
  })
  await Bun.sleep(20)
  await w.settle()
  expect(w.calls.filter((c) => c.method === "session.prompt").map((c) => c.params.text))
    .toEqual(["hello number 1", "hello number 2", "hello number 3"])
  expect(w.footer().join("\n")).not.toContain("the session is starting")
})
