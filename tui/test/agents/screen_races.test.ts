/**
 * The agent view's alternate screen (src/agents/screen.ts) on
 * OpenTUI's own renderer (the test renderer is the real CliRenderer).
 *
 * Found in review 8, each a lost or a stale row before the fix:
 *
 * A  rows committed just before ← (not yet presented by a frame) are
 *    counted in the seed: the footer comes back below them.
 * B  the way back while a frame is being drawn is deferred by OpenTUI;
 *    the held queue is released at once, and a replay queued meanwhile
 *    (a resize while the view was up) throws.
 * C  ← esc ← while that deferred way back is pending: the seed is taken
 *    on the alternate screen (renderOffset 0); the next way back puts the
 *    footer at the top.
 */
import { afterEach, expect, test } from "bun:test"
import { TextRenderable } from "@opentui/core"
import { createTestRenderer } from "@opentui/core/testing"
import { createAgentScreen } from "../../src/agents/screen"
import { CommitQueue } from "../../src/scrollback/queue"
import { commitSurface } from "../../src/scrollback/writers"

let setup: Awaited<ReturnType<typeof createTestRenderer>> | null = null
afterEach(async () => {
  await setup?.renderer.destroy()
  setup = null
})

async function mount() {
  setup = await createTestRenderer({ width: 60, height: 20, screenMode: "split-footer", footerHeight: 4,
                                     externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false })
  const r = setup.renderer
  const errors: string[] = []
  const queue = new CommitQueue({ requestLive: () => r.requestLive(), dropLive: () => r.dropLive() },
                                (e) => errors.push(String(e)))
  const offset = () => (r as unknown as { renderOffset: number }).renderOffset
  const rows = (n: number, tag: string) => (ctx: any) => ({
    root: new TextRenderable(ctx, { content: Array.from({ length: n }, (_, i) => `${tag} ${i + 1}`).join("\n") }),
  })
  const commit = (n: number, tag: string) => queue.run(() => commitSurface(r, rows(n, tag), 10))
  return { r, queue, errors, screen: createAgentScreen(r, queue), offset, commit }
}

async function frames(n = 3) {
  for (let i = 0; i < n; i++) {
    await setup!.renderOnce()
    await Bun.sleep(5)
  }
}

test("control: 3 rows committed and presented move the footer down by 3", async () => {
  const { offset, commit } = await mount()
  await frames()
  const start = offset()
  await commit(3, "row")
  await frames()
  expect(offset()).toBe(start + 3)
})

test("A: ← right after a commit (not yet presented): the footer must come back below it", async () => {
  const { screen, offset, commit } = await mount()
  await frames()
  const start = offset()
  void commit(3, "row")          // queued; written (enqueued) before the view goes up
  await screen.enter()           // no frame between the commit and the seed
  await frames()                 // the commit is presented, then the alternate screen
  screen.leave()
  await frames()
  console.log(`A: start ${start}, after the way back ${offset()} (want ${start + 3})`)
  expect(offset()).toBe(start + 3)
})

test("B: esc during a frame: a replay queued while the view was up must not fail", async () => {
  const { r, queue, errors, screen, commit } = await mount()
  await frames()
  await commit(3, "row")
  await frames()
  await screen.enter()
  await frames()
  // A resize while the view is up queues a replay (transcript.replay),
  // whose first step, resetSplitFooterForReplay, throws unless the
  // renderer's own _screenMode is split-footer.
  const R = r as unknown as { _screenMode: string; rendering: boolean }
  let modeWhenReplayRuns = ""
  void queue.run(() => { modeWhenReplayRuns = R._screenMode })
  const frame = setup!.renderOnce()   // a frame in flight (renderer.rendering)
  const renderingAtEsc = R.rendering
  screen.leave()
  await frame
  await queue.drain()
  await frames()
  console.log(`B: rendering at esc ${renderingAtEsc}; mode when the replay task ran: ${modeWhenReplayRuns}`)
  expect(modeWhenReplayRuns).toBe("split-footer")
})

test("C: ← esc ← with the way back still pending: the seed must not be the alternate screen's 0", async () => {
  const { r, screen, offset, commit } = await mount()
  await frames()
  await commit(5, "row")
  await frames()
  const start = offset()
  await screen.enter()
  await frames()
  const R = r as unknown as { _screenMode: string; rendering: boolean; renderOffset: number }
  const frame = setup!.renderOnce()   // esc arrives while a frame is drawn
  screen.leave()
  for (let i = 0; i < 5; i++) await Promise.resolve()
  console.log(`C: after esc: rendering ${R.rendering}, _screenMode ${R._screenMode}, renderOffset ${R.renderOffset}`)
  const again = screen.enter()        // ← again before that frame ended
  await frame
  await again
  await frames()
  screen.leave()
  await frames()
  console.log(`C: start ${start}, after ← esc ← esc ${offset()}`)
  expect(offset()).toBe(start)
})
