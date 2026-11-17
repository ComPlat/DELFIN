/**
 * The window talks only to a socket in a plain 0700 directory of its own
 * uid: one reached through a symlinked directory, or in a directory others
 * may write, may be planted and is refused before anything is sent. The
 * control: the same listener in a proper directory is connected to.
 */
import { expect, test } from "bun:test"
import { chmodSync, mkdtempSync, symlinkSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"
import { ServiceConnection, checkSocketDir } from "../../src/transport/rpc"

const HELLO = '{"protocolVersion":1,"serviceVersion":4,"pid":1,"host":"h","mayAct":true,"sessions":[]}'

function listener(path: string) {
  let opened = 0
  const s = Bun.listen({
    unix: path,
    socket: {
      open() { opened++ },
      data(socket, data) {
        for (const line of Buffer.from(data).toString().split("\n").filter(Boolean)) {
          const msg = JSON.parse(line)
          socket.write(`{"jsonrpc":"2.0","id":${msg.id},"result":${HELLO}}\n`)
        }
      },
    },
  })
  return { stop: () => s.stop(true), opened: () => opened }
}

test("a socket behind a symlinked directory is refused; the real one is used", async () => {
  const real = mkdtempSync(join(tmpdir(), "tui-sd-"))
  const link = join(mkdtempSync(join(tmpdir(), "tui-sl-")), "rt")
  symlinkSync(real, link)
  const planted = listener(join(real, "agent.sock"))
  const refused = await new ServiceConnection(join(link, "agent.sock"), [20]).connect().catch((e) => e)
  expect(String(refused?.message)).toContain("not a plain directory")
  expect(planted.opened()).toBe(0)
  // Control: the same listener, through its own directory.
  const conn = new ServiceConnection(join(real, "agent.sock"), [20])
  await conn.connect()
  expect(planted.opened()).toBe(1)
  conn.close()
  planted.stop()
})

test("a directory others may enter is refused", () => {
  const open = mkdtempSync(join(tmpdir(), "tui-so-"))
  chmodSync(open, 0o755)
  expect(() => checkSocketDir(join(open, "agent.sock"))).toThrow(/open to others/)
  chmodSync(open, 0o700)
  expect(() => checkSocketDir(join(open, "agent.sock"))).not.toThrow()
})
