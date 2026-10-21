/**
 * The window after resizes, on OpenTUI's test renderer: a size that comes
 * back to the one the scrollback was last written at is written again,
 * because the window was drawn at the other size in between (the tmux
 * resize storm of test/e2e/check_resize_tmux.py, in two steps).
 */
import { afterEach, expect, test } from "bun:test"
import { FIXTURE, mountWindow } from "../harness"

type Window = Awaited<ReturnType<typeof mountWindow>>
let open: Window | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

/** Past the session view's settle time (250 ms), then drawn. */
async function settled(w: Window) {
  await Bun.sleep(400)
  await w.settle()
}

test("a resize that ends at the width of the last replay writes the session again", async () => {
  const w = (open = await mountWindow({ width: 100, height: 30 }))
  for (const e of FIXTURE.events.slice(0, 30)) w.store.push(e)
  await w.settle()
  expect(w.view.transcript.replayCount).toBe(0)

  w.t.resize(80, 30)
  await settled(w)
  expect(w.view.transcript.replayCount).toBe(1)

  // Wider, then back to 80 before it settles: drawn at 120 meanwhile.
  w.t.resize(120, 30)
  w.t.resize(80, 30)
  await settled(w)
  expect(w.view.transcript.replayCount).toBe(2)

  // Nothing changed since: nothing to write.
  await settled(w)
  expect(w.view.transcript.replayCount).toBe(2)
})
