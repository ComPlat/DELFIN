/**
 * The add form of /clusters: a cluster for DELFIN, written by the agent
 * service (cluster.add: ctun's entry and remote.clusters) only for what a
 * person sent from here.
 *
 *   1  the cluster   ↑ ↓ one of the clusters ctun knows (clusters.presets),
 *                    or "another cluster" by its login host; enter picks it
 *   2  its fields    name, your user name there (with the preset's hint),
 *                    the login host (another cluster only), the partition
 *                    DELFIN's jobs go to, the work dir and install root;
 *                    tab ↑ ↓ move, enter the next field and on the last one
 *                    adds it, ctrl+u clears, esc goes back to step 1
 *   3  its host key  a preset's published key that the server presents is
 *                    pinned by the service at once ("matches the published
 *                    key"); any other key is shown here, and y pins exactly
 *                    that one (cluster.pin), n leaves it unpinned -- DELFIN
 *                    then refuses the cluster, and enter on its row asks
 *                    for the key again
 *
 * Then the login field opens for it (login.tsx). An agent's proposal (a
 * NEEDS YOU item, cluster_propose) opens the form at step 2 filled in: the
 * person checks every field and sends it, or esc sends nothing. While it
 * is open the form has every key (app.tsx routes them here first).
 */
import { decodePasteBytes, stripAnsiSequences, type KeyEvent } from "@opentui/core"
import { usePaste } from "@opentui/solid"
import { For, Match, Switch, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import type { Keymap } from "../keymap/keymap"
import { oneLine } from "../store/state"
import { fit } from "../text"
import type { Theme } from "../theme/theme"

/** A cluster ctun knows, as clusters.presets has it. */
export interface Preset {
  name: string
  title: string
  host: string
  requiresOtp: boolean
  userHint: string
  fingerprint: string
  fingerprintSource: string
  partition: string
}

/** What an agent proposed (a cluster_add item's fields). */
export interface Proposal {
  name?: string
  preset?: string
  user?: string
  host?: string
  partition?: string
}

export type FieldKey = "name" | "user" | "host" | "partition" | "workDir" | "installRoot"

export interface FormField {
  key: FieldKey
  label: string
  value: string
  hint: string
}

export interface HostKeyAsk {
  /** The key the server presents now, to pin once the person says y. */
  offered: string
  /** The preset's published key, when the server presents another one. */
  expected: string
  source: string
  error: string
}

export type AddStep = "preset" | "fields" | "hostkey"

export interface AddState {
  open: boolean
  step: AddStep
  presets: Preset[]
  /** The presets are being asked for. */
  looking: boolean
  /** The picked entry of step 1: a preset's index, presets.length for another cluster. */
  chosen: number
  fields: FormField[]
  focus: number
  /** A request is on its way (cluster.add, .hostkey, .pin). */
  busy: boolean
  error: string
  /** Filled in from an agent's proposal. */
  proposed: boolean
  /** The cluster that was added (steps 3 on). */
  name: string
  hostKey: HostKeyAsk | null
}

export interface AddDeps {
  request(method: string, params: Record<string, unknown>): Promise<any>
  keymap: Keymap
  say(text: string): void
  /** The cluster is added (and its key pinned or not): /clusters asks
   *  again, and with a pinned key its login field opens. */
  added(name: string, pinned: boolean): void
}

export interface AddForm {
  state: AddState
  /** Step 1, or step 2 filled in with a proposal. */
  open(proposal?: Proposal): void
  /** Step 3 for a cluster that is added already but whose key is not pinned. */
  openPin(name: string): void
  close(): void
  handleKey(key: KeyEvent): boolean
  paste(text: string): void
  /** Rows it takes in the footer. */
  height(): number
}

export const DEFAULT_WORK_DIR = "~/delfin/jobs"
export const DEFAULT_INSTALL_ROOT = "~/delfin/installs"

const str = (v: unknown, max: number) => oneLine(typeof v === "string" ? v : v == null ? "" : String(v), max)

/** The presets clusters.presets sent, each field as its type says. */
export function presetsOf(raw: unknown): Preset[] {
  if (!Array.isArray(raw)) return []
  return (raw as Record<string, unknown>[]).filter((p) => p && typeof p === "object" && str(p.name, 64)).map((p) => ({
    name: str(p.name, 64), title: str(p.title, 120) || str(p.name, 64), host: str(p.host, 120),
    requiresOtp: Boolean(p.requiresOtp), userHint: str(p.userHint, 160),
    fingerprint: str(p.fingerprint, 120), fingerprintSource: str(p.fingerprintSource, 200),
    partition: str(p.partition, 64),
  }))
}

export function fieldsFor(preset: Preset | undefined, proposal: Proposal = {}): FormField[] {
  const out: FormField[] = [
    { key: "name", label: "name", value: proposal.name || preset?.name || "",
      hint: "what DELFIN and ctun call it" },
    { key: "user", label: "user", value: proposal.user || "", hint: preset?.userHint || "your user name on the cluster" },
  ]
  // Another cluster, or a preset without a public host (IntNano): the host is asked.
  if (!preset || !preset.host) out.push({ key: "host", label: "login host", value: proposal.host || "", hint: "e.g. login.cluster.example.org" })
  out.push(
    { key: "partition", label: "partition", value: proposal.partition || preset?.partition || "",
      hint: "where DELFIN's jobs go; empty: the cluster's default" },
    { key: "workDir", label: "work dir", value: DEFAULT_WORK_DIR, hint: "the jobs' folders, on the cluster" },
    { key: "installRoot", label: "install root", value: DEFAULT_INSTALL_ROOT, hint: "DELFIN's installs, on the cluster" },
  )
  return out
}

export function createAddForm(deps: AddDeps): AddForm {
  const [state, setState] = createStore<AddState>({
    open: false, step: "preset", presets: [], looking: false, chosen: 0, fields: [], focus: 0,
    busy: false, error: "", proposed: false, name: "", hostKey: null,
  })
  /** Each opening's: an answer for one before is not this one's. */
  let generation = 0

  const preset = (): Preset | undefined => state.presets[state.chosen]

  const close = () => {
    generation++
    setState({ open: false, step: "preset", fields: [], focus: 0, busy: false, error: "", proposed: false,
               name: "", hostKey: null, looking: false })
  }

  const open = (proposal?: Proposal) => {
    const mine = ++generation
    setState({ open: true, step: "preset", presets: [], looking: true, chosen: 0, fields: [], focus: 0,
               busy: false, error: "", proposed: Boolean(proposal), name: "", hostKey: null })
    void deps.request("clusters.presets", {}).then((got) => {
      if (mine !== generation) return
      const presets = presetsOf(got?.presets)
      setState({ presets, looking: false, error: str(got?.error, 300) })
      if (!proposal) return
      const at = presets.findIndex((p) => p.name === proposal.preset)
      const chosen = at >= 0 ? at : presets.length
      if (proposal.preset && at < 0) setState({ error: `the proposed preset ${proposal.preset} is not one ctun knows: another cluster` })
      setState({ chosen, step: "fields", fields: fieldsFor(presets[at], proposal), focus: 0 })
    }, (error) => {
      if (mine === generation) setState({ looking: false, error: str((error as Error)?.message ?? error, 300) })
    })
  }

  /** The key the server presents, asked again (step 3). */
  const askKey = (name: string) => {
    const mine = generation
    setState({ busy: true, error: "" })
    void deps.request("cluster.hostkey", { name }).then((got) => {
      if (mine !== generation) return
      if (got?.pinned) {
        close()
        deps.say(`${name}'s key is pinned already`)
        return deps.added(name, true)
      }
      setState({ busy: false, hostKey: { offered: str(got?.offered, 120), expected: "", source: "",
                                         error: str(got?.error, 300) } })
    }, (error) => {
      if (mine === generation) setState({ busy: false, hostKey: { offered: "", expected: "", source: "",
                                                                  error: str((error as Error)?.message ?? error, 300) } })
    })
  }

  const openPin = (name: string) => {
    ++generation
    setState({ open: true, step: "hostkey", presets: [], looking: false, chosen: 0, fields: [], focus: 0,
               busy: false, error: "", proposed: false, name, hostKey: null })
    askKey(name)
  }

  const pick = () => {
    if (state.looking) return
    setState({ step: "fields", fields: fieldsFor(preset()), focus: 0, error: "" })
  }

  const value = (key: FieldKey) => state.fields.find((f) => f.key === key)?.value.trim() ?? ""

  const submit = () => {
    const empty = state.fields.findIndex((f) => (f.key === "name" || f.key === "user" || f.key === "host") && !f.value.trim())
    if (empty >= 0) return setState({ focus: empty, error: `the ${state.fields[empty]!.label} is empty` })
    const mine = generation
    const p = preset()
    const params: Record<string, unknown> = { name: value("name"), user: value("user") }
    if (p) params.preset = p.name
    if (!p || !p.host) params.host = value("host")
    for (const key of ["partition", "workDir", "installRoot"] as const) if (value(key)) params[key] = value(key)
    setState({ busy: true, error: "" })
    void deps.request("cluster.add", params).then((got) => {
      if (mine !== generation) return
      if (!got?.ok) return setState({ busy: false, error: str(got?.error || got?.code || "it was not added", 300) })
      const name = str(got.name, 64) || String(params.name)
      const key = got.hostKey && typeof got.hostKey === "object" ? got.hostKey : {}
      if (key.pinned) {
        close()
        deps.say(`added ${name}: the server presents the published key${key.source ? ` (${str(key.source, 120)})` : ""}, pinned`)
        return deps.added(name, true)
      }
      setState({ busy: false, step: "hostkey", name,
                 hostKey: { offered: str(key.offered, 120), expected: str(key.expected, 120),
                            source: str(key.source, 200), error: str(key.error, 300) } })
    }, (error) => {
      if (mine === generation) setState({ busy: false, error: str((error as Error)?.message ?? error, 300) })
    })
  }

  const pin = () => {
    const name = state.name
    const offered = state.hostKey?.offered
    if (!offered) return
    const mine = generation
    setState({ busy: true, error: "" })
    void deps.request("cluster.pin", { name, fingerprint: offered }).then((got) => {
      if (mine !== generation) return
      if (!got?.ok) return setState({ busy: false, error: str(got?.error || got?.code || "not pinned", 300) })
      close()
      deps.say(`pinned ${name}'s key ${offered}`)
      deps.added(name, true)
    }, (error) => {
      if (mine === generation) setState({ busy: false, error: str((error as Error)?.message ?? error, 300) })
    })
  }

  /** n or esc at step 3: the cluster stays added, its key unpinned. */
  const leaveUnpinned = () => {
    const name = state.name
    close()
    deps.say(`${name} is added, its key not pinned: DELFIN does not use it until it is (enter on its row asks again)`)
    deps.added(name, false)
  }

  const setField = (text: string) => {
    const f = state.fields[state.focus]
    if (f) setState("fields", state.focus, "value", text.slice(0, 200))
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    const plain = !key.ctrl && !key.meta
    const letter = plain ? (key.sequence ?? "").toLowerCase() : ""
    if (key.ctrl && key.name === "c") {
      if (state.step === "hostkey" && state.name) leaveUnpinned()
      else close()
      return true
    }
    if (state.busy) return true            // a request is on its way: the keys wait
    if (state.step === "preset") {
      const n = state.presets.length + 1
      if (key.name === "escape") close()
      else if (key.name === "up" || key.name === "down") setState({ chosen: (state.chosen + (key.name === "up" ? n - 1 : 1)) % n })
      else if (deps.keymap.matches("submit", key)) pick()
      return true
    }
    if (state.step === "hostkey") {
      if (letter === "y" && state.hostKey?.offered) pin()
      else if (letter === "r") askKey(state.name)
      else if (letter === "n" || key.name === "escape") leaveUnpinned()
      return true
    }
    // the fields
    const n = state.fields.length
    if (key.name === "escape") {
      if (state.proposed) close()
      else setState({ step: "preset", error: "" })
      return true
    }
    if (key.name === "tab" || key.name === "up" || key.name === "down") {
      const back = key.name === "up" || (key.name === "tab" && key.shift)
      setState({ focus: (state.focus + (back ? n - 1 : 1)) % n })
      return true
    }
    if (deps.keymap.matches("submit", key)) {
      if (state.focus < n - 1) setState({ focus: state.focus + 1 })
      else submit()
      return true
    }
    const f = state.fields[state.focus]
    if (!f) return true
    if (key.name === "backspace") setField([...f.value].slice(0, -1).join(""))
    else if (key.ctrl && key.name === "u") setField("")
    else if (plain && key.sequence && !/[\x00-\x1f\x7f]/.test(key.sequence)) setField(f.value + key.sequence)
    return true                            // nothing else leaves the form
  }

  return {
    state, open, openPin, close, handleKey,
    paste: (text) => {
      if (state.open && state.step === "fields" && !state.busy) {
        const f = state.fields[state.focus]
        if (f) setField(f.value + text.replace(/[\x00-\x1f\x7f]/g, ""))
      }
    },
    // The border, the heading, the step's lines and the status line.
    height: () => !state.open ? 0
      : 4 + (state.step === "preset" ? Math.max(1, state.presets.length + 1)
        : state.step === "fields" ? state.fields.length : 3),
  }
}

/** The hints line while the form is open: what its keys do now. */
export function addHelp(form: AddForm): string {
  const s = form.state
  if (s.busy) return "asking the agent service…"
  if (s.step === "preset") return "↑↓ select · enter pick · esc cancels"
  if (s.step === "hostkey") {
    return [s.hostKey?.offered ? "y pins this key" : "", "r asks again", "n leaves it unpinned"].filter(Boolean).join(" · ")
  }
  const last = s.focus >= s.fields.length - 1
  return [last ? "enter adds it" : "enter next field", "tab ↑↓ move", "ctrl+u clears",
          s.proposed ? "esc sends nothing" : "esc back"].join(" · ")
}

export function AddView(props: { form: AddForm; theme: () => Theme; width: () => number }): JSX.Element {
  const t = props.theme
  const f = props.form
  const inner = () => Math.max(10, props.width() - 4)
  usePaste((event) => {
    if (!f.state.open) return
    event.preventDefault()
    f.paste(stripAnsiSequences(decodePasteBytes(event.bytes)))
  })
  const title = () => f.state.step === "hostkey" ? ` the host key of ${fit(f.state.name, 30)} ` : " add a cluster "
  const heading = () => {
    const s = f.state
    if (s.step === "preset") return s.looking ? "asking ctun which clusters it knows…" : "which cluster?"
    if (s.step === "hostkey") return `${s.name} is added: does its server present the right key?`
    const p = s.presets[s.chosen]
    const which = p ? `${p.title} · ${p.host}` : "another cluster, by its login host"
    return (s.proposed ? "proposed by an agent: check every field · " : "") + which
  }
  const status = () => {
    const s = f.state
    if (s.busy) return s.step === "hostkey" ? "asking the server for its key…" : "adding it… (ctun's config and DELFIN's settings)"
    if (s.error) return s.error
    if (s.step === "preset") return "ctun's entry and DELFIN's are written only when you send the form"
    if (s.step === "fields") return s.fields[s.focus]?.hint ?? ""
    return ""
  }
  const keyLines = (): { text: string; look: "text" | "warn" | "err" | "muted" }[] => {
    const k = f.state.hostKey
    if (!k) return [{ text: "asking the server for its key…", look: "muted" }]
    if (k.error) return [{ text: k.error, look: "err" }, { text: "r asks again; n leaves it unpinned for now", look: "muted" }]
    const out: { text: string; look: "text" | "warn" | "err" | "muted" }[] = [{ text: `presented  ${k.offered}`, look: "text" }]
    if (k.expected) {
      out.push({ text: `published  ${k.expected}`, look: "err" })
      out.push({ text: "NOT the key its operators publish: pin it only if you know why it changed", look: "err" })
    } else {
      out.push({ text: "no published key to compare with: compare it with what the operators say", look: "warn" })
    }
    return out
  }
  return (
    <box flexShrink={0} flexDirection="column" border borderStyle="rounded" borderColor={t().accent}
         title={title()} titleAlignment="left" paddingLeft={1} paddingRight={1}>
      <text fg={t().muted} wrapMode="none">{fit(heading(), inner())}</text>
      <Switch>
        <Match when={f.state.step === "preset"}>
          <For each={[...f.state.presets.map((p) => `${p.title}  ·  ${p.host || "its host given"}${p.fingerprint ? "  · published key" : ""}`),
                      "another cluster  ·  by its login host"]}>
            {(line, i) => {
              const on = () => i() === f.state.chosen && !f.state.looking
              return (
                <text wrapMode="none">
                  <span style={{ fg: on() ? t().accent : t().faint }}>{on() ? "› " : "  "}</span>
                  <span style={{ fg: on() ? t().bright : t().text, bold: on() }}>{fit(line, inner() - 2)}</span>
                </text>
              )
            }}
          </For>
        </Match>
        <Match when={f.state.step === "fields"}>
          <For each={f.state.fields}>
            {(field, i) => {
              const on = () => i() === f.state.focus && !f.state.busy
              return (
                <text wrapMode="none">
                  <span style={{ fg: on() ? t().accent : t().faint }}>{on() ? "› " : "  "}</span>
                  <span style={{ fg: on() ? t().bright : t().muted }}>{field.label.padEnd(13)}</span>
                  <span style={{ fg: t().bright }}>{fit(field.value, inner() - 16)}</span>
                  <span style={{ fg: t().accent }}>{on() ? "▏" : ""}</span>
                </text>
              )
            }}
          </For>
        </Match>
        <Match when={f.state.step === "hostkey"}>
          <For each={keyLines()}>
            {(line) => (
              <text wrapMode="none" fg={line.look === "err" ? t().err : line.look === "warn" ? t().warn
                                         : line.look === "muted" ? t().faint : t().bright}>
                {fit(line.text, inner())}
              </text>
            )}
          </For>
        </Match>
      </Switch>
      <text fg={f.state.error && !f.state.busy ? t().err : t().faint} wrapMode="none">{fit(status(), inner())}</text>
    </box>
  )
}
