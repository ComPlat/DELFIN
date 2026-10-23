/**
 * The TUI's protocol types are generated from the Python definition
 * (delfin/agent/service/protocol.py). This fails when the committed file
 * and the definition drift apart: run `bun scripts/gen-protocol.ts`.
 */
import { describe, expect, test } from "bun:test"
import { GENERATED, generate, loadSchema, python } from "../../scripts/gen-protocol"
import { DURABLE_EVENTS, ERRORS, PROTOCOL_VERSION } from "../../src/protocol/generated"

const pythonFound = Bun.spawnSync([python(), "-c", "import delfin"], {
  cwd: `${import.meta.dir}/../../..`,
  env: { ...process.env, PYTHONPATH: `${import.meta.dir}/../../..` },
}).exitCode === 0

describe.skipIf(!pythonFound)("generated protocol types", () => {
  test("match the Python definition", async () => {
    expect(await Bun.file(GENERATED).text()).toBe(generate(loadSchema()))
  })

  test("carry the version, the durable events and the error codes", () => {
    const schema = loadSchema()
    expect(PROTOCOL_VERSION).toBe(schema.protocolVersion)
    expect([...DURABLE_EVENTS].sort()).toEqual(schema.durableEvents)
    expect(ERRORS.NOT_ALLOWED).toBe(-32003)
    expect(ERRORS.UNKNOWN_SESSION).toBe(-32004)
  })
})
