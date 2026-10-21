/**
 * The agent view's alternate screen: entering waits for what is queued,
 * holds what is queued after, and the way back puts the footer where it
 * was and writes what was held. A quick open-close-open ends open.
 */
import { afterEach, expect, test } from "bun:test"
import { createTestRenderer } from "@opentui/core/testing"
import { createAgentScreen } from "../../src/agents/screen"
import { CommitQueue } from "../../src/scrollback/queue"

let setup: Awaited<ReturnType<typeof createTestRenderer>> | null = null
afterEach(async () => {
  await setup?.renderer.destroy()
  setup = null
})

async function mount() {
  setup = await createTestRenderer({ width: 60, height: 20, screenMode: "split-footer", footerHeight: 4,
                                     externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false })
  const r = setup.renderer
  const queue = new CommitQueue({ requestLive: () => r.requestLive(), dropLive: () => r.dropLive() })
  return { r, queue, screen: createAgentScreen(r, queue),
           offset: () => (r as unknown as { renderOffset: number }).renderOffset }
}

test("what is queued before is written first, what comes after waits for the way back", async () => {
  const { r, queue, screen, offset } = await mount()
  await setup!.renderOnce()
  const before = offset()
  const order: string[] = []
  void queue.run(async () => {
    await Bun.sleep(30)
    order.push("before")
  })
  await screen.enter()
  expect(order).toEqual(["before"])
  void queue.run(() => { order.push("meanwhile") })
  expect(screen.active).toBe(true)
  expect(r.screenMode).toBe("alternate-screen")
  await Bun.sleep(20)
  expect(order).toEqual(["before"])                 // held while the view is up
  screen.leave()
  await queue.drain()
  expect(order).toEqual(["before", "meanwhile"])
  expect(r.screenMode).toBe("split-footer")
  expect(screen.active).toBe(false)
  await setup!.renderOnce()
  expect(offset()).toBe(before)
})

test("open, close, open in a row ends on the alternate screen", async () => {
  const { r, screen } = await mount()
  const first = screen.enter()
  screen.leave()
  const again = screen.enter()
  await first
  await again
  expect(r.screenMode).toBe("alternate-screen")
  expect(screen.active).toBe(true)
  screen.leave()
  await Bun.sleep(10)
  expect(r.screenMode).toBe("split-footer")
})
