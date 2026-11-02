/**
 * The /swarm view (SW2b): package cards with phase, members with health,
 * the "waits on you" queue with the dialog text and approve/deny, and a
 * new-swarm flow that asks before it starts anything. The view lives in
 * tui/src/swarm/ and follows the agents view's keymap conventions.
 *
 * Red control: all cases failed on the tree without tui/src/swarm/
 * (module missing). They go green with the view in place.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { createSwarmView, type SwarmStatus } from "../../src/swarm/view"
import { Keymap } from "../../src/keymap/keymap"

const now = 1700000000

/** A status payload as the service's swarm.status returns it (the shape of
 *  SW1's status contract, phase 0). */
function status(packages = 2, waiting = 1): SwarmStatus {
  const goal = "swarm goal: measure convergence"
  const pkgs = Array.from({ length: packages }, (_, i) => ({
    id: `p${i}`, name: `package ${i}`, phase: "builder", members: [
      { id: `m${i}a`, name: "s-a", health: "working", activity: "editing" },
      { id: `m${i}b`, name: "s-b", health: "waiting", activity: "a question" },
    ],
  }))
  const items = Array.from({ length: waiting }, (_, i) => ({
    item_id: `w${i}`, package: `p${i}`, question: `confirm package ${i}?`,
  }))
  return { schema_version: 1, goal, packages: pkgs, waiting_for_user: items }
}

function makeView(st = status()) {
  const drafts: string[] = []
  const requests: { method: string; params: Record<string, unknown> }[] = []
  const view = createSwarmView({
    keymap: new Keymap(),
    now: () => now,
    askStatus: async () => { requests.push({ method: "swarm.status", params: {} }); return st },
    askValidate: async (spec) => { requests.push({ method: "swarm.validate", params: { specYaml: spec } });
      return { valid: true, errors: [] } },
    askStart: async (spec, confirm) => { requests.push({ method: "swarm.start",
      params: { specYaml: spec, dryRun: !confirm, confirm: Boolean(confirm) } });
      return confirm ? { dry_run: false, goal: "goal", started: [], refused: null }
        : { dry_run: true, goal: "goal", started: [], dry_run_text: "plan: 2 members", refused: null } },
    askTriage: async () => { requests.push({ method: "swarm.triage.list", params: {} }); return st.waiting_for_user },
    askAnswer: async (item, answer) => { requests.push({ method: "swarm.triage.answer",
      params: { item, answer } }); return { answered: true, item_id: String(item.item_id) } },
    note: (text) => { drafts.push(text) },
    draft: () => "",
  })
  return { view, requests, notes: drafts }
}

afterEach(() => { /* stateless: the view holds no globals */ })

describe("the /swarm view's data", () => {
  test("packages become rows with phase and member health", async () => {
    const { view } = makeView()
    await view.refresh()
    expect(view.state.open).toBe(true)
    const rows = view.rows()
    expect(rows.length).toBe(2)
    expect(rows[0]!.name).toBe("package 0")
    expect(rows[0]!.phase).toBe("builder")
    expect(rows[0]!.health).toBe("1 working · 1 waiting")
  })

  test("the waits-on-you queue carries the question and the answers", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    const queue = view.waits()
    expect(queue.length).toBe(1)
    expect(queue[0]!.question).toContain("confirm package 0")
    expect(requests.some((r) => r.method === "swarm.triage.list")).toBe(true)
  })

  test("an answered item goes out through swarm.triage.answer with the item itself", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    const [item] = view.waits()
    const done = await view.answer(item!, "yes")
    expect(done).toBe(true)
    const call = requests.find((r) => r.method === "swarm.triage.answer")!
    expect(call.params.item).toEqual(item!.raw)
    expect(call.params.answer).toBe("yes")
  })
})

describe("the new-swarm flow", () => {
  test("dry-run first: a start without confirmation sends dryRun and asks again", async () => {
    const { view, requests, notes } = makeView()
    await view.refresh()
    await view.newSwarm("goal: do the thing")
    // The plan came back; the flow shows it and waits, it has not started.
    expect(requests.filter((r) => r.method === "swarm.start").length).toBe(1)
    const start = requests.find((r) => r.method === "swarm.start")!
    expect(start.params.dryRun).toBe(true)
    expect(start.params.confirm).not.toBe(true)
    expect(notes.join("\n")).toContain("plan: 2 members")
    expect(view.state.pendingStart).toBeTruthy()
  })

  test("start only after an explicit confirmation", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    await view.newSwarm("goal: do the thing")
    const started = await view.confirmStart()
    expect(started).toBe(true)
    const start = requests.filter((r) => r.method === "swarm.start").pop()!
    expect(start.params.dryRun).toBe(false)
    expect(start.params.confirm).toBe(true)
  })

  test("an invalid spec never reaches a start", async () => {
    // Validation fails -> no plan, no start; the notes say why.
    const notes: string[] = []
    const starts: unknown[] = []
    const bad = createSwarmView({
      keymap: new Keymap(), now: () => now,
      askStatus: async () => status(),
      askValidate: async () => ({ valid: false, errors: ["members: empty"] }),
      askStart: async () => { starts.push("start"); return { dry_run: true, goal: "x", started: [] } },
      askTriage: async () => [],
      askAnswer: async () => ({ answered: true }),
      note: (t) => notes.push(t), draft: () => "",
    }) as never as { state: { pendingStart: unknown }; newSwarm(g: string): Promise<boolean> }
    await bad.newSwarm("goal: broken")
    expect(bad.state.pendingStart).toBeFalsy()
    expect(starts.length).toBe(0)
    expect(notes.join("\n")).toContain("not valid")
  })

  test("the flow edits the spec in an editor when one is offered", async () => {
    const { view } = makeView()
    await view.refresh()
    const edited = await view.editSpec("goal: first\n", (text) => text + "members: []\n")
    expect(edited).toContain("members: []")
    void view
  })
})

describe("keyboard-first navigation", () => {
  test("up and down move over packages and queue rows, esc closes", async () => {
    const { view } = makeView()
    await view.refresh()
    const km = new Keymap()
    const down = { name: "down", sequence: "down", ctrl: false, meta: false, shift: false }
    const up = { name: "up", sequence: "up", ctrl: false, meta: false, shift: false }
    view.handleKey(down as never)
    expect(view.state.selected).toBe(1)
    view.handleKey(up as never)
    expect(view.state.selected).toBe(0)
    view.handleKey({ name: "escape", sequence: "esc", ctrl: false, meta: false, shift: false } as never)
    expect(view.state.open).toBe(false)
    void km
  })

  test("enter on a waits-on-you row answers it, once confirmed", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    const right = { name: "return", sequence: "return", ctrl: false, meta: false, shift: false }
    view.toggleQueue()
    view.handleKey(right as never)
    // No answer sent without the confirmation step: the view asks first.
    expect(requests.some((r) => r.method === "swarm.triage.answer")).toBe(false)
    expect(view.state.pendingAnswer).toBeTruthy()
    const sent = await view.confirmAnswer("yes")
    expect(sent).toBe(true)
    expect(requests.some((r) => r.method === "swarm.triage.answer")).toBe(true)
  })

  test("a goal that is YAML-shaped or multi-line cannot change the spec's shape", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    // A typed line carrying YAML keys must not become spec structure.
    await view.newSwarm("goal: steal\npackages: [x]")
    const plan = requests.find((r) => r.method === "swarm.start")
    expect(plan).toBeTruthy()
    expect(plan!.params.specYaml).not.toContain("\npackages")
    expect(String(plan!.params.specYaml).split("\n").length).toBe(2)  // goal + trailing newline
    // An empty goal makes no plan at all.
    const empty = await makeView().view.newSwarm("   ")
    expect(empty).toBe(false)
    expect(requests.filter((r) => r.method === "swarm.start").length).toBe(1)
  })

  test("y confirms a held plan, n lets it go, another key is neither", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    await view.newSwarm("goal: do the thing")
    expect(view.state.pendingStart).toBeTruthy()
    const y = { name: "y", sequence: "y", ctrl: false, meta: false, shift: false }
    const n = { name: "n", sequence: "n", ctrl: false, meta: false, shift: false }
    const other = { name: "z", sequence: "z", ctrl: false, meta: false, shift: false }
    // z is neither: the plan stays held.
    view.handleKey(other as never)
    expect(view.state.pendingStart).toBeTruthy()
    view.handleKey(n as never)
    expect(view.state.pendingStart).toBe(null)
    expect(requests.filter((r) => r.method === "swarm.start").length).toBe(1)  // the dry run only
    await view.newSwarm("goal: again")
    view.handleKey(y as never)
    // y is the explicit confirmation: the held spec goes out confirmed.
    const start = requests.filter((r) => r.method === "swarm.start").pop()!
    expect(start.params.confirm).toBe(true)
    expect(start.params.dryRun).toBe(false)
  })

  test("y does nothing while an ANSWER is held: the answer is typed, not keyed", async () => {
    const { view, requests } = makeView()
    await view.refresh()
    view.toggleQueue()
    view.handleKey({ name: "return", sequence: "return", ctrl: false, meta: false, shift: false } as never)
    expect(view.state.pendingAnswer).toBeTruthy()
    view.handleKey({ name: "y", sequence: "y", ctrl: false, meta: false, shift: false } as never)
    // No answer went out: the text is typed and enter sends it.
    expect(requests.some((r) => r.method === "swarm.triage.answer")).toBe(false)
    expect(view.state.pendingAnswer).toBeTruthy()
  })
})
