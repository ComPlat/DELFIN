/**
 * The window, headless: the app, the store and the session view wired the
 * way main.tsx wires them, on OpenTUI's test renderer, with a fake service
 * that records what it is asked. Shared by the tests that drive the
 * window with events and keys.
 */
import { createTestRenderer } from "@opentui/core/testing"
import { render } from "@opentui/solid"
import { Pins } from "../src/agents/pins"
import { App, type Requester } from "../src/app"
import { Keymap } from "../src/keymap/keymap"
import { createScreenItems } from "../src/scrollback/screenview"
import { startSessionView } from "../src/session"
import { ensureSession } from "../src/store/state"
import { createAppStore } from "../src/store/store"
import { DARK, type Theme } from "../src/theme/theme"
import { ungutter } from "./gutter"

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
  onResume?: (id: string, workspace: string) => void
  onNewSession?: () => Promise<string | null>
  pins?: Pins
  /** The whole terminal (the window's default) instead of inline. */
  full?: boolean
  /** The theme, as a signal (dark by default). */
  theme?: () => Theme
  /** Inline under GNU screen: a resize (or ctrl+o) cannot write the
   *  session again. */
  divider?: boolean
} = {}) {
  const t = await createTestRenderer({
    width: options.width ?? 100, height: options.height ?? 30,
    screenMode: options.full ? "alternate-screen" : "split-footer", footerHeight: 4,
    externalOutputMode: options.full ? "passthrough" : "capture-stdout", useMouse: Boolean(options.full),
    // As main.tsx: buttons and drags only, and no click-to-focus.
    enableMouseMovement: false, autoFocus: false,
    exitOnCtrlC: false,
  })
  const store = createAppStore()
  const theme = options.theme ?? (() => DARK)
  const screen = options.full ? createScreenItems() : undefined
  const view = startSessionView(t.renderer, store, theme, { showThinking: false, replayTurns: 20, screen,
                                                           dividerOnResize: Boolean(options.divider) })
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
         onQuit={() => { quits++ }} onRedraw={() => {}} onExpand={() => view.expandLast()}
         details={view.details} canShowDetails={() => view.canShowDetails} onDetails={() => view.toggleDetails()}
         onLines={(title, lines) => view.lines(title, lines)}
         onResume={options.onResume} onNewSession={options.onNewSession} pins={options.pins}
         answerPauseMs={options.answerPauseMs ?? 0} screen={screen} />
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
    return frame.slice(frame.length - t.renderer.footerHeight - 1).map(ungutter)
  }
  const close = async () => {
    view.stop()
    await t.renderer.destroy()
  }
  // A lone Esc is held 20 ms by the input parser (it may start a sequence).
  const escape = async () => {
    t.mockInput.pressEscape()
    await Bun.sleep(40)
    await settle()
  }
  /** The whole screen's lines, without the gutter (full screen). */
  const screenLines = () => t.captureCharFrame().split("\n").map(ungutter)
  return { t, keys: t.mockInput, store, view, calls, settle, footer, close, escape, quits: () => quits, sid,
           screen, screenLines,
           replies: () => calls.filter((c) => c.method === "permission.reply").map((c) => c.params) }
}
