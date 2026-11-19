/**
 * The masked login field: a cluster's password and OTP, typed into the
 * window's own state and sent in exactly one cluster.login, which the
 * agent service hands to ctun and nowhere else (clusters_feed.login).
 *
 *   typing      into the focused field, each character drawn as •
 *   backspace   deletes one          ctrl+u      clears the field
 *   tab ↑ ↓     the other field      paste       into the focused field
 *   enter       the next field; on the last one, log in (not with a
 *               field empty)
 *   esc         cancel: both fields cleared, nothing sent
 *
 * Opened from /clusters (enter on a logged-out cluster) or from a login
 * item of NEEDS YOU (its cluster's row asked for first). While it is open
 * it has every key (app.tsx routes them here before anything else) and
 * the prompt has no focus: nothing typed here reaches the composer, the
 * input history, a session's prompt, a question's answer or an attention
 * item's. What is typed lives in this module's closure, not in a store a
 * view draws -- the view knows how many characters each field holds --
 * and is let go the moment the request goes out. It is in no note, log or
 * line on the screen: an error text that held it would show *** instead.
 * A login may take a minute or two; "logging in…" says it is on its way.
 */
import { decodePasteBytes, stripAnsiSequences, type KeyEvent } from "@opentui/core"
import { usePaste } from "@opentui/solid"
import { For, Show, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import type { Keymap } from "../keymap/keymap"
import { fit } from "../text"
import type { Theme } from "../theme/theme"
import { clusterRows, refusal, type ClusterRow } from "./view"

export type Field = "password" | "otp"

export interface LoginState {
  open: boolean
  name: string
  target: string
  fields: Field[]
  focus: number
  /** How many characters each field holds: what the view draws, as •. */
  lengths: Record<Field, number>
  /** The cluster's row is being asked for (a login item names only it). */
  looking: boolean
  /** cluster.login is on its way. */
  busy: boolean
  error: string
}

export interface LoginDeps {
  request(method: string, params: Record<string, unknown>): Promise<any>
  keymap: Keymap
  say(text: string): void
  /** A login went through (an open /clusters asks again). */
  loggedIn(name: string): void
}

export interface Login {
  state: LoginState
  /** Open for a cluster of clusters.list. */
  open(row: ClusterRow): void
  /** Open for a cluster known by its name (a login item's): its row is
   *  asked for; `typed` is what was typed for it already. */
  openFor(name: string, typed?: string): void
  close(): void
  handleKey(key: KeyEvent): boolean
  /** Pasted text, into the focused field (line breaks left out). */
  paste(text: string): void
  /** Rows it takes in the footer. */
  height(): number
}

const LABEL: Record<Field, string> = { password: "password", otp: "OTP" }

const fieldsOf = (row: ClusterRow): Field[] =>
  [...(row.requiresPassword ? ["password" as const] : []), ...(row.requiresOtp ? ["otp" as const] : [])]

/** `text` with every secret in it masked: the service names none, and a
 *  window does not rely on it. */
function masked(text: string, secrets: string[]): string {
  return secrets.reduce((out, secret) => out.split(secret).join("***"), text)
}

export function createLogin(deps: LoginDeps): Login {
  const [state, setState] = createStore<LoginState>({
    open: false, name: "", target: "", fields: [], focus: 0, lengths: { password: 0, otp: 0 },
    looking: false, busy: false, error: "",
  })
  /** What is typed: here, and nowhere else until it is sent. */
  let typed: Record<Field, string> = { password: "", otp: "" }
  /** Each opening's: an answer for one before is not this one's. */
  let generation = 0

  const forget = () => {
    typed = { password: "", otp: "" }
    setState("lengths", { password: 0, otp: 0 })
  }
  const set = (field: Field, value: string) => {
    typed[field] = value
    setState("lengths", field, [...value].length)
  }
  const put = (text: string) => {
    const field = state.fields[state.focus]
    if (field && !state.busy && text) set(field, typed[field] + text)
  }

  const open = (row: ClusterRow) => {
    generation++
    forget()
    setState({ open: true, name: row.name, target: row.target, fields: fieldsOf(row), focus: 0, looking: false,
               busy: false, error: "" })
  }
  const close = () => {
    generation++
    forget()
    setState({ open: false, name: "", target: "", fields: [], focus: 0, looking: false, busy: false, error: "" })
  }
  /** Closed, and why in the message row. */
  const end = (why: string) => {
    close()
    deps.say(why)
  }

  const openFor = (name: string, first = "") => {
    const mine = ++generation
    forget()
    // What is typed before the row is in goes into its first field.
    setState({ open: true, name, target: "", fields: ["password"], focus: 0, looking: true, busy: false, error: "" })
    put(first)
    void deps.request("clusters.list", {}).then((got) => {
      if (mine !== generation) return
      const row = clusterRows(got?.clusters).find((r) => r.name === name)
      if (!row) return end(`${name || "that cluster"} is not one of DELFIN's clusters (remote.clusters does not name it)`)
      if (row.state === "refused") return end(refusal(row))
      if (row.state === "logged_in") return end(`already logged in to ${name}`)
      const early = typed.password
      const fields = fieldsOf(row)
      forget()
      setState({ looking: false, target: row.target, fields })
      if (fields[0]) set(fields[0], early)
    }, (error) => {
      if (mine === generation) end(`${name}: ${String((error as Error)?.message ?? error)}`)
    })
  }

  /** Enter on the last field: the one request, with what was typed --
   *  let go here at once. An empty field sends nothing (a login that
   *  fails may count against the account). */
  const submit = () => {
    const empty = state.fields.findIndex((field) => !typed[field])
    if (empty >= 0) {
      setState({ focus: empty, error: `the ${LABEL[state.fields[empty]!]} is empty` })
      return
    }
    const mine = generation
    const name = state.name
    const params: Record<string, unknown> = { name }
    for (const field of state.fields) params[field] = typed[field]
    const secrets = state.fields.map((field) => typed[field]).filter(Boolean)
    forget()
    setState({ busy: true, error: "", focus: 0 })
    const failed = (why: string) => {
      const said = masked(why, secrets)
      if (mine === generation) setState({ busy: false, error: said, focus: 0 })
      else deps.say(`not logged in to ${name}: ${said}`)
    }
    void deps.request("cluster.login", params).then((got) => {
      if (!got?.ok) return failed(String(got?.error || got?.code || "the login did not go through"))
      if (mine === generation) close()
      deps.say(got.opened === false ? `${name} was logged in already` : `logged in to ${name}`)
      deps.loggedIn(name)
    }, (error) => failed(String((error as Error)?.message ?? error)))
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    if (key.name === "escape" || (key.ctrl && key.name === "c")) {
      close()
      return true
    }
    if (state.busy) return true            // logging in: the keys wait
    const n = state.fields.length
    if (key.name === "tab" || key.name === "up" || key.name === "down") {
      if (n > 1) setState({ focus: (state.focus + (key.name === "up" || key.shift ? n - 1 : 1)) % n })
      return true
    }
    if (deps.keymap.matches("submit", key)) {
      if (state.looking) return true
      if (state.focus < n - 1) setState({ focus: state.focus + 1 })
      else submit()
      return true
    }
    const field = state.fields[state.focus]
    if (!field) return true
    if (key.name === "backspace") set(field, [...typed[field]].slice(0, -1).join(""))
    else if (key.ctrl && key.name === "u") set(field, "")
    else if (!key.ctrl && !key.meta && key.sequence && !/[\x00-\x1f\x7f]/.test(key.sequence)) put(key.sequence)
    return true                            // nothing else leaves the field
  }

  return {
    state, open, openFor, close, handleKey,
    paste: (text) => put(text.replace(/[\x00-\x1f\x7f]/g, "")),
    // The border, the cluster's line, a line a field (or one saying none
    // is needed) and the status line.
    height: () => state.open ? 4 + Math.max(1, state.fields.length) : 0,
  }
}

/** The hints line while the field is open: what its keys do now. */
export function loginHelp(login: Login): string {
  const s = login.state
  if (s.busy) return "esc closes it: the login goes on, and a note says how it went"
  const last = s.focus >= s.fields.length - 1
  return [last ? "enter logs in" : "enter next field", s.fields.length > 1 ? "tab other field" : "",
          s.fields.length ? "ctrl+u clears" : "", "esc cancels"].filter(Boolean).join(" · ")
}

export function LoginView(props: { login: Login; theme: () => Theme; width: () => number }): JSX.Element {
  const t = props.theme
  const l = props.login
  const inner = () => Math.max(10, props.width() - 4)
  // A paste goes into the focused field (the prompt has neither the keys
  // nor pastes while this is open).
  usePaste((event) => {
    if (!l.state.open) return
    event.preventDefault()
    l.paste(stripAnsiSequences(decodePasteBytes(event.bytes)))
  })
  const status = () => l.state.looking ? `looking up ${l.state.name}…`
    : l.state.busy ? "logging in… (it can take a minute or two)"
    : l.state.error || "typed here, it goes to ctun only: never to a session, a model or the history"
  return (
    <box flexShrink={0} flexDirection="column" border borderStyle="rounded" borderColor={t().accent}
         title={` log in to ${fit(l.state.name, Math.max(4, inner() - 12))} `} titleAlignment="left"
         paddingLeft={1} paddingRight={1}>
      <text fg={t().muted} wrapMode="none">{fit([l.state.name, l.state.target].filter(Boolean).join(" · "), inner())}</text>
      <Show when={l.state.fields.length} fallback={
        <text fg={t().faint} wrapMode="none">{fit("nothing to type: enter logs in", inner())}</text>
      }>
        <For each={l.state.fields}>
          {(field, i) => {
            const on = () => i() === l.state.focus && !l.state.busy
            return (
              <text wrapMode="none">
                <span style={{ fg: on() ? t().accent : t().faint }}>{on() ? "› " : "  "}</span>
                <span style={{ fg: on() ? t().bright : t().muted }}>{LABEL[field].padEnd(10)}</span>
                <span style={{ fg: t().bright }}>{fit("•".repeat(l.state.lengths[field]), inner() - 13)}</span>
                <span style={{ fg: t().accent }}>{on() ? "▏" : ""}</span>
              </text>
            )
          }}
        </For>
      </Show>
      <text fg={l.state.error && !l.state.busy && !l.state.looking ? t().err : t().faint} wrapMode="none">
        {fit(status(), inner())}
      </text>
    </box>
  )
}
