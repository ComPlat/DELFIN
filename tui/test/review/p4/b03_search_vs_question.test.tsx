/**
 * B03: a ctrl+r search left open when a question opens. ctrl+c is meant to
 * work everywhere (app.tsx header): with the search open it is skipped by
 * the ladder (`!composer.search()`), and with a question open the search
 * never sees it either -- ctrl+c does nothing, the turn cannot be ended.
 */
import { afterEach, expect, test } from "bun:test"
import { FIXTURE, mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const edit = FIXTURE.events.find((e: any) => e.type === "permission.asked" && e.kind === "confirm")
const upToEdit = FIXTURE.events.slice(0, FIXTURE.events.indexOf(edit) + 1)

test("ctrl+c with a ctrl+r search open and a question open", async () => {
  w = await mount({ history: [{ text: "run orca opt" }] })
  for (const e of upToEdit.slice(0, -1)) w.store.push(e)   // a turn runs
  await w.settle()
  w.keys.pressKey("r", { ctrl: true })                     // the user searches the history...
  await w.settle()
  w.store.push(edit)                                       // ...and a question opens
  await w.settle()
  // The question closes the search (the draft comes back), and takes the keys.
  expect(w.footer().join("\n")).not.toContain("search history")
  w.keys.pressKey("c", { ctrl: true })                     // end the turn
  await w.settle()
  w.keys.pressKey("c", { ctrl: true })
  await w.settle()
  const cancels = w.calls.filter((c) => c.method === "session.cancel").length
  console.log("session.cancel calls:", cancels, "replies:", JSON.stringify(w.replies()))
  console.log(w.footer().filter((l) => l.trim()).join("\n"))
  expect(cancels).toBeGreaterThan(0)
})
