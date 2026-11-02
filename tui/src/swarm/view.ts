/**
 * The /swarm view (SW2b): manage swarms from the terminal UI, over the
 * service's swarm.* requests (swarm_rpc.py). Keyboard-first, following the
 * agents view's conventions (arrows, return, esc; the keymap owns keys).
 *
 *   goal · packages as cards (phase, members with health)
 *   WAITS ON YOU: the triage queue, one row per question; answer with a
 *   confirmation step, never a bare key.
 *
 * The new-swarm flow is the gate that matters: goal -> validate -> plan
 * (dry run) -> shown and held -> start only on an explicit confirmation.
 * An invalid spec never reaches a start; a refused start is shown, not
 * swallowed.
 */
import type { KeyEvent } from "@opentui/core"
import { createStore } from "solid-js/store"
import type { Keymap } from "../keymap/keymap"
import { width } from "../text"

/** One status payload (swarm.status / SW1's contract). */
export interface SwarmStatus {
  schema_version: number
  goal: string
  packages: Record<string, unknown>[]
  waiting_for_user: Record<string, unknown>[]
}

/** A triage item as swarm.triage.list returns it: the fields the view
 *  reads, over the raw record the service gave (the answer goes with the
 *  RAW record; the service re-derives the verdict there). */
export interface TriageItem extends Record<string, unknown> {
  item_id: string
  package?: string
  question?: string
  raw: Record<string, unknown>
}

export interface SwarmPlanResult {
  goal: string
  members?: Record<string, unknown>[]
  dry_run_text: string
}

export interface SwarmStartResult {
  dry_run: boolean
  goal: string
  started: Record<string, unknown>[]
  dry_run_text?: string | null
  refused?: string | null
}

export interface SwarmViewDeps {
  keymap: Keymap
  now: () => number
  askStatus(): Promise<SwarmStatus>
  askValidate(specYaml: string): Promise<{ valid: boolean; errors: string[] }>
  askStart(specYaml: string, confirm: boolean): Promise<SwarmStartResult>
  askTriage(): Promise<Record<string, unknown>[]>
  /** The answer goes with the item the service gave: it re-derives the
   *  verdict there (swarm_rpc.py passes item + answer only). */
  askAnswer(item: Record<string, unknown>, answer: string): Promise<{ answered: boolean }>
  note(text: string): void
  draft(): string
}

export interface SwarmRow {
  id: string
  name: string
  phase: string
  health: string
  members: Record<string, unknown>[]
}

export interface SwarmViewState {
  open: boolean
  selected: number
  queueOpen: boolean
  queueSelected: number
  status: SwarmStatus | null
  queue: TriageItem[]
  loading: boolean
  /** The new-swarm flow's held plan: a start waits on the user. */
  pendingStart: { specYaml: string; plan: SwarmStartResult } | null
  /** A queue row taken for an answer; it goes out only after confirm. */
  pendingAnswer: { item: TriageItem } | null
}

/** The members' healths, counted: "2 working · 1 waiting · 1 failed". */
export function memberHealth(members: Record<string, unknown>[]): string {
  const counts = new Map<string, number>()
  for (const m of members) {
    const health = String(m.health ?? "unknown")
    counts.set(health, (counts.get(health) ?? 0) + 1)
  }
  if (!counts.size) return "no members"
  return [...counts.entries()].map(([health, n]) => `${n} ${health}`).join(" · ")
}

export function createSwarmView(deps: SwarmViewDeps) {
  const [state, setState] = createStore<SwarmViewState>({
    open: false, selected: 0, queueOpen: false, queueSelected: 0,
    status: null, queue: [], loading: false, pendingStart: null, pendingAnswer: null,
  })

  const rows = (): SwarmRow[] =>
    (state.status?.packages ?? []).map((p) => ({
      id: String(p.id ?? ""),
      name: String(p.name ?? p.id ?? "package"),
      phase: String(p.phase ?? "unknown"),
      health: memberHealth(Array.isArray(p.members) ? p.members : []),
      members: Array.isArray(p.members) ? p.members : [],
    }))

  const waits = (): TriageItem[] =>
    (state.queue ?? []).map((q) => ({
      item_id: String(q.item_id ?? ""),
      package: q.package === undefined ? undefined : String(q.package ?? ""),
      question: q.question === undefined ? undefined : String(q.question ?? ""),
      raw: (q.raw ?? q) as Record<string, unknown>,
    }))

  const refresh = async (): Promise<void> => {
    setState({ loading: true })
    try {
      const st = await deps.askStatus()
      const items = (await deps.askTriage().catch(() => [] as Record<string, unknown>[]))
        .map((raw) => ({ ...raw, item_id: String(raw.item_id ?? ""), raw }))
      setState({ status: st, queue: items, open: true, loading: false })
    } catch (error) {
      setState({ loading: false })
      throw error
    }
  }

  const open = async (): Promise<void> => {
    await refresh()
  }
  const close = (): void => {
    setState({ open: false })
  }

  const select = (list: unknown[], at: number, n = 1): number =>
    list.length ? ((at + n) % list.length + list.length) % list.length : 0

  /** Return/enter on a queue row: hold the item, ask for the answer; the
   *  answer itself goes out through confirmAnswer. */
  const takeAnswer = (item: TriageItem): boolean => {
    setState({ pendingAnswer: { item } })
    deps.note(`answer for ${item.package ?? item.item_id}: type the answer (or yes / no) and confirm`)
    return true
  }

  /** The direct answer of an item (the agents view's answerItem path):
   *  it holds the item and asks for the confirmation step. */
  const answer = async (item: TriageItem, text: string): Promise<boolean> => {
    setState({ pendingAnswer: { item } })
    return confirmAnswer(text)
  }

  const confirmAnswer = async (answer: string): Promise<boolean> => {
    const held = state.pendingAnswer
    if (!held) {
      deps.note("nothing is held for an answer")
      return false
    }
    if (!answer.trim()) return false
    const done = await deps.askAnswer(held.item.raw, answer.trim()).catch((error) => {
      deps.note(`not answered: ${String((error as Error)?.message ?? error)}`)
      return { answered: false }
    })
    if (done.answered) {
      setState({ pendingAnswer: null, queue: state.queue.filter((q) => q.raw !== held.item.raw) })
      deps.note(`answered ${held.item.package ?? held.item.item_id}`)
      return true
    }
    return false
  }

  const newSwarm = async (goalText: string): Promise<boolean> => {
    // The goal is one line, as one quoted scalar: a typed line must never
    // add spec structure (a "goal: x\npackages: [...]" line would change
    // the spec's shape). Newlines fold to spaces, quotes are doubled in
    // the double-quoted form.
    const goal = goalText.trim().replace(/\s+/g, " ")
    if (!goal) return false
    const specYaml = `goal: "${goal.replace(/"/g, '\\"')}"\n`
    const check = await deps.askValidate(specYaml).catch(() => ({ valid: false, errors: ["validation failed"] }))
    if (!check.valid) {
      deps.note(`the spec is not valid: ${check.errors.join("; ")}`)
      return false
    }
    const plan = await deps.askStart(specYaml, false).catch((error) => {
      deps.note(`no plan: ${String((error as Error)?.message ?? error)}`)
      return null
    }) as SwarmStartResult | null
    if (!plan) return false
    // Held, never started: the service's dry-run plan is shown and the
    // user confirms; only confirmStart sends confirm=true.
    setState({ pendingStart: { specYaml, plan } })
    deps.note(plan.dry_run_text || `plan for ${plan.goal}: confirm to start`)
    return true
  }

  /** The explicit confirmation. Without a held plan it refuses. */
  const confirmStart = async (): Promise<boolean> => {
    const held = state.pendingStart
    if (!held) {
      deps.note("no plan is held: /swarm new first")
      return false
    }
    const result = await deps.askStart(held.specYaml, true).catch((error) => {
      deps.note(`not started: ${String((error as Error)?.message ?? error)}`)
      return null
    })
    if (!result) return false
    setState({ pendingStart: null })
    if (result.refused) {
      deps.note(`the service refused the start: ${result.refused}`)
      return false
    }
    deps.note(`swarm started: ${result.goal} (${result.started.length} members)`)
    void refresh()
    return true
  }

  /** The flow's spec edit: the offered editor transforms the text (the
   *  window opens its real $EDITOR; tests pass a transform). */
  const editSpec = async (text: string, transform?: (t: string) => string): Promise<string> => {
    if (!transform) return text
    return transform(text)
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    const km = deps.keymap
    const typing = deps.draft().length > 0
    if (key.name === "escape" && !typing) {
      if (state.pendingAnswer) setState({ pendingAnswer: null })
      else if (state.pendingStart) setState({ pendingStart: null })
      else if (state.queueOpen) setState({ queueOpen: false })
      else close()
      return true
    }
    // The confirmation keys, while a plan or an answer is held (the
    // agents view's question keys): y confirms the HELD thing only,
    // n lets it go; any other key is neither. Not while typing.
    if (!typing && state.pendingStart && key.name === "n") {
      setState({ pendingStart: null })
      deps.note("the plan is let go: /swarm new starts again")
      return true
    }
    if (!typing && state.pendingStart && key.name === "y") {
      void confirmStart()
      return true
    }
    if (!typing && state.pendingAnswer && key.name === "n") {
      setState({ pendingAnswer: null })
      deps.note("the answer is let go")
      return true
    }
    if (km.matches("agentsUp", key) || km.matches("agentsDown", key)) {
      if (typing) return true
      const up = km.matches("agentsUp", key)
      if (state.queueOpen) {
        setState({ queueSelected: select(waits(), state.queueSelected, up ? -1 : 1) })
      } else {
        setState({ selected: select(rows(), state.selected, up ? -1 : 1) })
      }
      return true
    }
    if (km.matches("agentsOpen", key) && !typing) {
      if (state.queueOpen) {
        const item = waits()[state.queueSelected]
        if (item) return takeAnswer(item)
        return true
      }
      // On a package row: its members' detail goes into the note.
      const row = rows()[state.selected]
      if (row) deps.note(`${row.name} · ${row.phase} · ${row.health}`)
      return true
    }
    return false
  }

  const toggleQueue = (): boolean => {
    if (!state.open) return false
    setState({ queueOpen: !state.queueOpen })
    return true
  }

  return { state, rows, waits, refresh, open, close, handleKey, newSwarm, confirmStart, confirmAnswer, answer,
           editSpec, toggleQueue }
}
