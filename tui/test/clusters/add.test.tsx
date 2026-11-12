/**
 * Adding and removing a cluster in /clusters (clusters/add.tsx, view.tsx):
 * a opens the add form -- the clusters ctun knows, or another one by its
 * host -- and enter on its last field sends exactly one cluster.add with
 * what was typed. A published key the server presents comes back pinned
 * and the login field opens; any other key is shown, and only y pins
 * exactly that one (cluster.pin). x asks first, and only y removes. An
 * agent's proposal (a NEEDS YOU item) opens the form filled in, is never
 * answered, and sends nothing until the person does. Rendered headless
 * through the real App against a fake agent service.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { mkdtempSync, rmSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"
import { presetsOf } from "../../src/clusters/add"
import { clusterRows } from "../../src/clusters/view"
import { History } from "../../src/composer/history"
import type { AttentionItem } from "../../src/protocol/generated"
import { mountWindow, type Call } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
const dirs: string[] = []
afterEach(async () => {
  await w?.close()
  w = null
  for (const d of dirs.splice(0)) rmSync(d, { recursive: true, force: true })
})

const PUBLISHED = "SHA256:5mZYEpKigwK5ibBMHRrh3WIkOtCqomJW6H7OMbPk3ec"
const OTHER = "SHA256:another-key-0000000000000000000000000000000"
const PRESETS = [
  { name: "bwuni3", title: "bwUniCluster 3.0 (KIT / bwHPC)", host: "uc3.scc.kit.edu", requiresOtp: true,
    requiresPassword: true, userHint: "KIT users: 'ka_' + your KIT username, e.g. ka_ab1234",
    budget: { limit: 200, unit: "jobh" }, fingerprint: PUBLISHED,
    fingerprintSource: "https://wiki.bwhpc.de/e/BwUniCluster3.0/Login", partition: "cpu" },
  { name: "kcist", title: "KCIST (KIT)", host: "login2.kcist.kit.edu", requiresOtp: false, requiresPassword: true,
    userHint: "", budget: null, fingerprint: null, fingerprintSource: null, partition: "small" },
]
const row = (extra: Record<string, unknown>) => ({
  name: "kcist", target: "me@login2.kcist.kit.edu", state: "logged_out", problems: [], requiresPassword: true,
  requiresOtp: false, workDir: "~/delfin/jobs", defaultPartition: "small", budget: null, jobs: null, orca: null,
  error: null, ...extra,
})

/** The fake service's clusters: kcist, and bwuni3 once it was added. */
const listed = (calls: Call[]) => [row({}), ...(calls.some((c) => c.method === "cluster.add")
  ? [row({ name: "bwuni3", target: "ka_ab1234@uc3.scc.kit.edu", requiresOtp: true })] : [])]

async function mount(options: { version?: number; clusters?: (calls: Call[]) => Record<string, unknown>[];
                                add?: (params: Record<string, any>) => unknown;
                                hostkey?: () => unknown; pin?: (params: Record<string, any>) => unknown } = {}) {
  const calls: Call[] = []
  const dir = mkdtempSync(join(tmpdir(), "delfin-add-"))
  dirs.push(dir)
  const history = new History(join(dir, "tui_history.jsonl"))
  const client = {
    request: async (method: string, params: Record<string, any> = {}) => {
      calls.push({ method, params })
      if (method === "clusters.list") return { clusters: options.clusters?.(calls) ?? [row({})] }
      if (method === "clusters.presets") return { presets: PRESETS, error: null }
      if (method === "cluster.add") {
        return options.add?.(params) ?? { ok: true, code: null, error: null, name: params.name,
                                          hostKey: { pinned: true, fingerprint: PUBLISHED, source: PRESETS[0]!.fingerprintSource,
                                                     offered: null, expected: PUBLISHED, error: null } }
      }
      if (method === "cluster.hostkey") return options.hostkey?.() ?? { offered: OTHER, pinned: false, error: null }
      if (method === "cluster.pin") return options.pin?.(params) ?? { ok: true, code: null, error: null }
      if (method === "cluster.remove") return { ok: true, code: null, error: null, loggedOut: false, keysRemoved: 1 }
      if (method === "cluster.login") return { ok: true, opened: true, code: null, error: null }
      if (method === "sessions.list") return { sessions: [{ sessionId: "s1", state: "idle" }], elsewhere: [] }
      if (method === "session.replay") return { events: [] }
      if (method === "attention.answer" || method === "attention.dismiss") return { ok: true }
      return { accepted: true, applied: true, message: "" }
    },
  }
  w = await mountWindow({ width: 110, height: 34, full: true, client, history })
  w.store.update((d) => { d.serviceVersion = options.version ?? 7 })
  await w.settle()
  return { calls, history }
}

const screen = () => w!.screenLines().join("\n")
const called = (calls: Call[], method: string) => calls.filter((c) => c.method === method).map((c) => c.params)
const line = async (text: string) => {
  await w!.keys.typeText(text)
  w!.keys.pressEscape()
  await Bun.sleep(40)
  w!.keys.pressEnter()
  await w!.settle()
}
const settled = async () => {
  for (let i = 0; i < 3; i++) await w!.settle()
}
const key = async (k: string) => {
  w!.keys.pressKey(k)
  await settled()
}
const type = async (text: string) => {
  await w!.keys.typeText(text)
  await w!.settle()
}
const enter = async () => {
  w!.keys.pressEnter()
  await settled()
}
const down = async () => {
  w!.keys.pressArrow("down")
  await w!.settle()
}

describe("reading what the service sends", () => {
  test("the presets and a row's ORCA, each field as its type says", () => {
    expect(presetsOf([...PRESETS, { title: "no name" }, null]).map((p) => [p.name, p.fingerprint, p.partition]))
      .toEqual([["bwuni3", PUBLISHED, "cpu"], ["kcist", "", "small"]])
    const [r] = clusterRows([row({ state: "logged_in", orca: { found: false, path: "", hint: "put your ORCA into the home" } })])
    expect(r!.orca).toEqual({ found: false, path: "", hint: "put your ORCA into the home" })
  })
})

describe("the add form", () => {
  test("a opens it with the clusters ctun knows; enter picks one and fills in what it knows", async () => {
    const { calls } = await mount()
    await line("/clusters")
    await settled()
    expect(screen()).toContain("a add")
    await key("a")
    let text = screen()
    expect(called(calls, "clusters.presets")).toHaveLength(1)
    expect(text).toContain("add a cluster")
    expect(text).toContain("› bwUniCluster 3.0 (KIT / bwHPC)  ·  uc3.scc.kit.edu  · published key")
    expect(text).toContain("another cluster  ·  by its login host")
    await enter()
    text = screen()
    expect(text).toMatch(/› name\s+bwuni3▏/)
    expect(text).toMatch(/partition\s+cpu/)
    expect(text).toMatch(/work dir\s+~\/delfin\/jobs/)
    expect(text).toMatch(/install root\s+~\/delfin\/installs/)
    await down()
    expect(screen()).toContain("KIT users: 'ka_' + your KIT username")      // the user field's hint
    expect(called(calls, "cluster.add")).toEqual([])
  })

  test("enter on the last field sends one cluster.add; a published key that matches is pinned, then the login field", async () => {
    const { calls } = await mount({ clusters: listed })
    await line("/clusters")
    await settled()
    await key("a")
    await enter()                                    // bwuni3
    await down()
    await type("ka_ab1234")
    for (let i = 0; i < 4; i++) await enter()        // user -> partition -> work dir -> install root -> add
    await settled()
    expect(called(calls, "cluster.add")).toEqual([{ name: "bwuni3", user: "ka_ab1234", preset: "bwuni3",
                                                    partition: "cpu", workDir: "~/delfin/jobs",
                                                    installRoot: "~/delfin/installs" }])
    const text = screen()
    expect(text).toContain("log in to bwuni3")      // the login field, for the cluster just added
    expect(text).toContain("ka_ab1234@uc3.scc.kit.edu")
    expect(called(calls, "cluster.pin")).toEqual([])
    expect(called(calls, "clusters.list").length).toBeGreaterThanOrEqual(2)
  })

  test("another key is shown beside the published one, and y pins exactly the key shown", async () => {
    const { calls } = await mount({ clusters: listed, add: (p) => ({ ok: true, code: null, error: null, name: p.name,
      hostKey: { pinned: false, fingerprint: null, source: PRESETS[0]!.fingerprintSource, offered: OTHER,
                 expected: PUBLISHED, error: null } }) })
    await line("/clusters")
    await settled()
    await key("a")
    await enter()
    await down()
    await type("ka_ab1234")
    for (let i = 0; i < 4; i++) await enter()
    const text = screen()
    expect(text).toContain(`presented  ${OTHER}`)
    expect(text).toContain(`published  ${PUBLISHED}`)
    expect(text).toContain("NOT the key its operators publish")
    expect(called(calls, "cluster.pin")).toEqual([])
    await key("y")
    expect(called(calls, "cluster.pin")).toEqual([{ name: "bwuni3", fingerprint: OTHER }])
    expect(screen()).toContain("log in to bwuni3")
  })

  test("another cluster is given by its host; n leaves its key unpinned, and enter on its row asks again", async () => {
    let added = false
    const { calls } = await mount({
      clusters: () => added ? [row({ name: "lab", target: "me@login.lab.example", state: "refused",
                                     problems: ["ctun has not pinned its host key; check the fingerprint"] })] : [],
      add: (p) => {
        added = true
        return { ok: true, code: null, error: null, name: p.name,
                 hostKey: { pinned: false, fingerprint: null, source: null, offered: OTHER, expected: null, error: null } }
      },
    })
    await line("/clusters")
    await settled()
    expect(screen()).toContain("no clusters yet: a adds one")
    await key("a")
    await down()
    await down()                                     // past kcist: another cluster
    await enter()
    await type("lab")
    await enter()
    await type("me")
    await enter()
    await type("login.lab.example")
    for (let i = 0; i < 4; i++) await enter()
    expect(called(calls, "cluster.add")).toEqual([{ name: "lab", user: "me", host: "login.lab.example",
                                                    workDir: "~/delfin/jobs", installRoot: "~/delfin/installs" }])
    expect(screen()).toContain("no published key to compare with")
    await key("n")
    expect(called(calls, "cluster.pin")).toEqual([])
    expect(screen()).toContain("lab is added, its key not pinned")
    await settled()
    expect(screen()).toContain("enter pin its key")
    await enter()                                    // on the refused row: step 3, the key asked again
    expect(called(calls, "cluster.hostkey")).toEqual([{ name: "lab" }])
    expect(screen()).toContain(`presented  ${OTHER}`)
    await key("y")
    expect(called(calls, "cluster.pin")).toEqual([{ name: "lab", fingerprint: OTHER }])
  })

  test("an empty field or the service's no keeps the form open and says why", async () => {
    const { calls } = await mount({ add: () => ({ ok: false, code: "invalid", name: "kcist",
                                                  error: "partition 'bad name' is not a name Slurm would accept", hostKey: null }) })
    await line("/clusters")
    await settled()
    await key("a")
    await down()
    await enter()                                    // kcist
    for (let i = 0; i < 5; i++) await enter()        // the user is empty
    expect(called(calls, "cluster.add")).toEqual([])
    expect(screen()).toContain("the user is empty")
    await type("me")
    for (let i = 0; i < 4; i++) await enter()
    expect(called(calls, "cluster.add")).toHaveLength(1)
    expect(screen()).toContain("is not a name Slurm would accept")
    expect(screen()).toContain("add a cluster")      // still open
  })

  test("esc goes back to the clusters, then closes the form; nothing was sent", async () => {
    const { calls } = await mount()
    await line("/clusters")
    await settled()
    await key("a")
    await enter()
    w!.keys.pressEscape()
    await Bun.sleep(40)
    await settled()
    expect(screen()).toContain("which cluster?")
    w!.keys.pressEscape()
    await Bun.sleep(40)
    await settled()
    expect(screen()).not.toContain("add a cluster")
    expect(screen()).toContain("◆ CLUSTERS")
    expect(called(calls, "cluster.add")).toEqual([])
  })
})

describe("removing a cluster", () => {
  test("x asks first; only y removes it, in one cluster.remove", async () => {
    const { calls } = await mount()
    await line("/clusters")
    await settled()
    await key("x")
    expect(screen()).toContain("remove kcist from DELFIN and ctun? y removes")
    await key("q")
    expect(called(calls, "cluster.remove")).toEqual([])
    expect(screen()).toContain("kept kcist")
    await key("x")
    await key("y")
    expect(called(calls, "cluster.remove")).toEqual([{ name: "kcist" }])
    expect(screen()).toContain("removed kcist from DELFIN and ctun")
  })
})

test("with a service before version 7, a and x say to restart it and send nothing", async () => {
  const { calls } = await mount({ version: 6 })
  await line("/clusters")
  await settled()
  expect(screen()).not.toContain("a add")
  await key("a")
  expect(screen()).toContain("restart the agent service to add or remove clusters here")
  await key("x")
  expect(called(calls, "clusters.presets")).toEqual([])
  expect(called(calls, "cluster.remove")).toEqual([])
})

test("a logged-in cluster without ORCA says so under its row", async () => {
  await mount({ clusters: () => [row({ state: "logged_in", orca: { found: false, path: "",
                                                                   hint: "put your ORCA into the cluster's home" } })] })
  await line("/clusters")
  await settled()
  expect(screen()).toContain("no ORCA where DELFIN looks: put your ORCA into the cluster's home")
})

describe("an agent's proposal in NEEDS YOU", () => {
  const proposal: AttentionItem = {
    id: "att-p", kind: "addCluster", group: "needs", blocking: false, name: "add the cluster bwuni3",
    text: "an agent proposes bwuni3", detail: "An agent proposes bwuni3 (bwuni3) for DELFIN.", workspace: "/w",
    host: "", sessionId: "s1", sessionTitle: "", sessionOpen: false, sessionSaved: false, createdAt: Date.now() / 1000 - 30,
    options: [], job: null, cluster: "bwuni3",
    proposal: { name: "bwuni3", preset: "bwuni3", user: "ka_ab1234", host: "", partition: "cpu" },
  }

  test("enter opens the form filled in over /clusters; it is checked and sent by the person, never answered", async () => {
    const { calls } = await mount()
    w!.store.push({ type: "attention.updated", session: "", ts: 0, items: [proposal] } as never)
    await w!.settle()
    w!.keys.pressArrow("left")
    await w!.settle()
    for (let i = 0; i < 6 && !screen().split("\n").some((l) => l.startsWith(" › ") && l.includes("add the cluster bwuni3")); i++) {
      await down()
    }
    w!.keys.pressKey(" ")
    await w!.settle()
    expect(screen()).toContain("→ opens the add form filled in for bwuni3")
    await enter()
    const text = screen()
    expect(text).toContain("proposed by an agent: check every field")
    expect(text).toMatch(/user\s+ka_ab1234/)
    expect(text).toMatch(/partition\s+cpu/)
    expect(called(calls, "cluster.add")).toEqual([])
    expect(called(calls, "attention.answer")).toEqual([])
    for (let i = 0; i < 5; i++) await enter()
    expect(called(calls, "cluster.add")).toEqual([{ name: "bwuni3", user: "ka_ab1234", preset: "bwuni3", partition: "cpu",
                                                    workDir: "~/delfin/jobs", installRoot: "~/delfin/installs" }])
  })

  test("esc on a proposal sends nothing", async () => {
    const { calls } = await mount()
    w!.store.push({ type: "attention.updated", session: "", ts: 0, items: [proposal] } as never)
    await w!.settle()
    w!.keys.pressArrow("left")
    await w!.settle()
    for (let i = 0; i < 6 && !screen().split("\n").some((l) => l.startsWith(" › ") && l.includes("add the cluster bwuni3")); i++) {
      await down()
    }
    await enter()
    expect(screen()).toContain("proposed by an agent")
    w!.keys.pressEscape()
    await Bun.sleep(40)
    await settled()
    expect(screen()).not.toContain("proposed by an agent")
    expect(called(calls, "cluster.add")).toEqual([])
    expect(called(calls, "attention.answer")).toEqual([])
  })
})
