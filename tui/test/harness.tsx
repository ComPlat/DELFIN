/**
 * The window, headless: the app, the store and the session view wired the
 * way main.tsx wires them, on OpenTUI's test renderer, with a fake service
 * that records what it is asked. Shared by the tests that drive the
 * window with events and keys.
 */
import { createTestRenderer } from "@opentui/core/testing"
import { render } from "@opentui/solid"
import { App, type Requester } from "../src/app"
import { Keymap } from "../src/keymap/keymap"
import { startSessionView } from "../src/session"
import { ensureSession } from "../src/store/state"
import { createAppStore } from "../src/store/store"
import { DARK } from "../src/theme/theme"

export const FIXTURE = await Bun.file(`${import.meta.dir}/fixtures/fake-conversation.json`).json()

export interface Call {
  method: string
  params: Record<string, any>
}

export async function mountWindow(options: {
  width?: number
  height?: number
  client?: Requester
  sid?: string
  show?: boolean
  answerPauseMs?: number
} = {}) {
  const t = await createTestRenderer({
    width: options.width ?? 100, height: options.height ?? 30, screenMode: "split-footer", footerHeight: 4,
    externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false,
  })
  const store = createAppStore()
  const theme = () => DARK
  const view = startSessionView(t.renderer, store, theme, { showThinking: false, replayLimit: 400 })
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
    store.update((d) => { ensureSession(d, sid, FIXTURE.session) })
    view.show(sid)
  }
  await render(() => (
    <App renderer={t.renderer} store={store} keymap={new Keymap()} theme={theme} client={client}
         onQuit={() => { quits++ }} onRedraw={() => {}} answerPauseMs={options.answerPauseMs ?? 0} />
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
  const close = async () => {
    view.stop()
    await t.renderer.destroy()
  }
  return { t, store, view, calls, settle, footer, close, quits: () => quits, sid,
           replies: () => calls.filter((c) => c.method === "permission.reply").map((c) => c.params) }
}
