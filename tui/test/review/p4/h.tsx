/** Review 6 harness: the real App + session view on OpenTUI's test renderer,
 *  wired as main.tsx wires them, from the pristine export of 4aef3b0a. */
import { createTestRenderer } from "@opentui/core/testing"
import { render } from "@opentui/solid"
import { App, type Requester } from "../../../src/app"
import { History, type Entry } from "../../../src/composer/history"
import { Keymap } from "../../../src/keymap/keymap"
import { startSessionView } from "../../../src/session"
import { ensureSession } from "../../../src/store/state"
import { createAppStore } from "../../../src/store/store"
import { DARK } from "../../../src/theme/theme"

export const FIXTURE = await Bun.file(`${import.meta.dir}/../../fixtures/fake-conversation.json`).json()

export interface Call { method: string; params: Record<string, any> }

export async function mount(options: {
  width?: number; height?: number; client?: Requester; sid?: string; show?: boolean
  answerPauseMs?: number; history?: Entry[]; onResume?: (id: string, ws: string) => void
  turn?: boolean
} = {}) {
  const t = await createTestRenderer({
    width: options.width ?? 100, height: options.height ?? 30, screenMode: "split-footer", footerHeight: 4,
    externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false,
  })
  const store = createAppStore()
  const theme = () => DARK
  const view = startSessionView(t.renderer, store, theme, { showThinking: false, replayTurns: 20 })
  const calls: Call[] = []
  const client = options.client ?? {
    request: async (method: string, params: Record<string, any> = {}) => {
      calls.push({ method, params })
      return { applied: true, message: "", accepted: true }
    },
  }
  let quits = 0
  const sid = options.sid ?? "s1"
  if (options.show !== false) {
    store.update((d) => {
      const s = ensureSession(d, sid, FIXTURE.session)
      s.state = "idle"
      if (options.turn) s.turn = { id: "t1", source: "user", startedAt: Date.now(), label: "thinking…",
                                   base: { inputTokens: 0, outputTokens: 0, costUsd: 0 }, streamed: 0 }
    })
    view.show(sid)
  }
  const history = new History(null, options.history ?? [])
  await render(() => (
    <App renderer={t.renderer} store={store} keymap={new Keymap()} theme={theme} client={client}
         history={history} sentAs={view.sentAs}
         onQuit={() => { quits++ }} onRedraw={() => {}} onExpand={() => view.expandLast()}
         onLines={(title, lines) => view.lines(title, lines)}
         onResume={options.onResume}
         answerPauseMs={options.answerPauseMs ?? 0} />
  ), t.renderer)
  const settle = async () => {
    store.flush()
    await view.queue.drain()
    await t.flush()
    await t.renderOnce()
    await Bun.sleep(5)
    await t.renderOnce()
  }
  const footer = () => {
    const frame = t.captureCharFrame().split("\n")
    return frame.slice(frame.length - t.renderer.footerHeight - 1)
  }
  const frame = () => t.captureCharFrame()
  const escape = async () => {
    t.mockInput.pressEscape()
    await Bun.sleep(40)
    await settle()
  }
  const close = async () => {
    view.stop()
    await t.renderer.destroy()
  }
  return { t, keys: t.mockInput, store, view, calls, settle, footer, frame, escape, close, history,
           quits: () => quits, sid,
           prompts: () => calls.filter((c) => c.method === "session.prompt").map((c) => c.params),
           replies: () => calls.filter((c) => c.method === "permission.reply").map((c) => c.params) }
}
