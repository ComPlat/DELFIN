/**
 * JSON-RPC over the agent service's Unix socket: one JSON object per line.
 *
 * `RpcClient` speaks the framing and matches answers to requests; it is
 * transport-agnostic so a test hands it a fake. `ServiceConnection` opens
 * the socket with Bun, says hello, and reconnects when the service goes
 * away and comes back (an upgrade, a restart): requests fail while it is
 * down, events stop, and `onUp` runs again once it is back, so the app can
 * attach to its sessions again. A service that vanishes without a word (it
 * crashed) is started again (`startService`); one that says it stops on
 * purpose (`service.stopping`) is let go: the connection closes for good.
 */
import type { AnyEvent, Requests, Results } from "../protocol/generated"
import { ERRORS, PROTOCOL_VERSION } from "../protocol/generated"

export class RpcError extends Error {
  constructor(readonly code: number, message: string) {
    super(message)
  }
}

/** The connection to the service is gone (not a session's worker). */
export const CONNECTION_LOST = -1

/** What carries lines both ways. */
export interface Wire {
  write(line: string): void
  close(): void
}

type Pending = { resolve: (value: any) => void; reject: (error: Error) => void; method: string }

export class RpcClient {
  private nextId = 0
  private pending = new Map<number, Pending>()
  private buffer = ""
  // One per connection: a lost one may have stopped inside a character.
  private decoder = new TextDecoder()
  private readonly listeners = new Set<(event: AnyEvent) => void>()

  constructor(private wire: Wire | null = null) {}

  attach(wire: Wire): void {
    this.wire = wire
  }

  onEvent(listener: (event: AnyEvent) => void): () => void {
    this.listeners.add(listener)
    return () => this.listeners.delete(listener)
  }

  request<M extends keyof Requests & keyof Results>(method: M, params: Requests[M],
                                                     options?: { timeoutMs?: number }): Promise<Results[M]>
  request(method: string, params?: Record<string, unknown>, options?: { timeoutMs?: number }): Promise<any>
  request(method: string, params: Record<string, unknown> = {},
          options: { timeoutMs?: number } = {}): Promise<any> {
    if (!this.wire) return Promise.reject(new RpcError(CONNECTION_LOST, "not connected"))
    const id = ++this.nextId
    return new Promise((resolve, reject) => {
      let timer: ReturnType<typeof setTimeout> | undefined
      const done = (fn: () => void) => {
        if (timer) clearTimeout(timer)
        fn()
      }
      this.pending.set(id, {
        resolve: (value) => done(() => resolve(value)),
        reject: (error) => done(() => reject(error)),
        method,
      })
      if (options.timeoutMs) {
        timer = setTimeout(() => {
          if (!this.pending.delete(id)) return
          reject(new RpcError(CONNECTION_LOST, `${method}: no answer within ${Math.round(options.timeoutMs! / 1000)} s`))
        }, options.timeoutMs)
      }
      try {
        this.wire!.write(JSON.stringify({ jsonrpc: "2.0", id, method, params }) + "\n")
      } catch (error) {
        this.pending.delete(id)
        done(() => reject(new RpcError(CONNECTION_LOST, `could not send ${method}: ${(error as Error).message}`)))
      }
    })
  }

  /** Bytes from the wire, in order; whole lines are handled as they complete. */
  feed(chunk: Uint8Array | string): void {
    this.buffer += typeof chunk === "string" ? chunk : this.decoder.decode(chunk, { stream: true })
    let nl: number
    while ((nl = this.buffer.indexOf("\n")) >= 0) {
      const line = this.buffer.slice(0, nl).trim()
      this.buffer = this.buffer.slice(nl + 1)
      if (line) this.handle(line)
    }
  }

  private handle(line: string): void {
    let message: any
    try {
      message = JSON.parse(line)
    } catch {
      return
    }
    if (!message || typeof message !== "object") return
    if (message.method === "event") {
      for (const listener of this.listeners) {
        try {
          listener(message.params as AnyEvent)
        } catch (error) {
          console.error("event listener failed:", error)
        }
      }
      return
    }
    const waiting = this.pending.get(message.id)
    if (!waiting) return
    this.pending.delete(message.id)
    if (message.error) {
      waiting.reject(new RpcError(Number(message.error.code), String(message.error.message)))
    } else {
      waiting.resolve(message.result)
    }
  }

  /** The wire is gone: every request still waiting fails. */
  lost(reason: string): void {
    this.wire = null
    this.buffer = ""
    this.decoder = new TextDecoder()
    for (const [id, waiting] of this.pending) {
      this.pending.delete(id)
      waiting.reject(new RpcError(CONNECTION_LOST, `${waiting.method}: ${reason}`))
    }
  }
}

export type ConnectionState = "connecting" | "up" | "reconnecting" | "closed"

export interface ServiceHello {
  protocolVersion: number
  serviceVersion: number
  pid: number
  host: string
  mayAct: boolean
  sessions: Record<string, unknown>[]
}

/** How long a service gets to answer hello. */
export const HELLO_TIMEOUT_MS = 10_000

/** Failed reconnects before a vanished service is started again, and how
 *  many more before the next try if that did not bring one up. */
export const START_AFTER_ATTEMPTS = 2
export const START_AGAIN_EVERY = 10

export interface ConnectionOptions {
  /** Start the next service (one vanished without a word); true when one
   *  runs. Never called after a stop on purpose. */
  startService?: () => Promise<boolean>
}

export class ServiceConnection {
  readonly rpc = new RpcClient()
  state: ConnectionState = "connecting"
  hello: ServiceHello | null = null
  /** Why the connection closed for good, when it did on its own. */
  problem = ""
  /** It closed because the service stopped on purpose (not a crash). */
  stoppedOnPurpose = false
  private closing = false
  private socket: any = null
  /** What the service said before it went: it stops on purpose. */
  private stopping: { reason: string; restart: boolean } | null = null
  private readonly stateListeners = new Set<(state: ConnectionState) => void>()
  private readonly upListeners = new Set<(hello: ServiceHello) => void>()

  constructor(readonly path: string, private readonly retryMs: number[] = [200, 500, 1000, 2000, 5000],
              private readonly options: ConnectionOptions = {}) {
    this.rpc.onEvent((event) => {
      const e = event as { type?: string; reason?: string; restart?: boolean }
      if (e.type === "service.stopping") this.stopping = { reason: String(e.reason ?? "stop"), restart: Boolean(e.restart) }
    })
  }

  onState(listener: (state: ConnectionState) => void): () => void {
    this.stateListeners.add(listener)
    return () => this.stateListeners.delete(listener)
  }

  /** Runs after every successful hello: the first connect and each reconnect. */
  onUp(listener: (hello: ServiceHello) => void): () => void {
    this.upListeners.add(listener)
    return () => this.upListeners.delete(listener)
  }

  private setState(state: ConnectionState): void {
    // Said once per change: a listener that closes the connection on
    // "closed" must not hear "closed" again (it would never end).
    if (this.state === state) return
    this.state = state
    for (const listener of this.stateListeners) listener(state)
  }

  /** Connect and say hello; resolves with the hello, rejects only when the
   *  first attempt fails (a window with no service has nothing to show). */
  async connect(): Promise<ServiceHello> {
    await this.open()
    return this.greet()
  }

  private async open(): Promise<void> {
    const rpc = this.rpc
    const encoder = new TextEncoder()
    // Bytes the socket did not take yet; written on drain, in order.
    let unsent = new Uint8Array(0)
    const self = this
    const flush = (socket: any) => {
      while (unsent.length) {
        const wrote = socket.write(unsent)
        if (wrote <= 0) return
        unsent = unsent.subarray(wrote)
      }
    }
    this.socket = await Bun.connect({
      unix: this.path,
      socket: {
        data(_socket, data) {
          rpc.feed(data)
        },
        drain(socket) {
          flush(socket)
        },
        close() {
          self.lost("the agent service closed the connection")
        },
        error(_socket, error) {
          self.lost(`the connection failed: ${error.message}`)
        },
      },
    })
    const socket = this.socket
    rpc.attach({
      write(line: string) {
        const bytes = encoder.encode(line)
        if (unsent.length) {
          const joined = new Uint8Array(unsent.length + bytes.length)
          joined.set(unsent)
          joined.set(bytes, unsent.length)
          unsent = joined
          return
        }
        unsent = bytes
        flush(socket)
      },
      close() {
        socket.end()
      },
    })
  }

  private async greet(): Promise<ServiceHello> {
    const hello = await this.rpc.request("hello", { protocolVersion: PROTOCOL_VERSION, client: "delfin-tui" },
                                         { timeoutMs: HELLO_TIMEOUT_MS })
    this.hello = hello as ServiceHello
    this.setState("up")
    for (const listener of this.upListeners) listener(this.hello)
    return this.hello
  }

  private lost(reason: string): void {
    // Whatever waits on this connection fails now -- also a hello sent
    // while reconnecting, which would otherwise wait forever.
    this.rpc.lost(reason)
    this.socket = null
    if (this.state === "reconnecting" || this.state === "closed") return
    if (this.closing) {
      this.setState("closed")
      return
    }
    const stopped = this.stopping
    this.stopping = null
    if (stopped && !stopped.restart) {
      // Stopped on purpose (`service stop`, an emergency stop): starting
      // another would undo it. The sessions are saved.
      this.problem = `the agent service stopped (${stopped.reason}); the sessions are saved, `
        + "and `delfin-agent tui` offers them again"
      this.stoppedOnPurpose = true
      this.closing = true
      this.setState("closed")
      return
    }
    this.setState("reconnecting")
    void this.reconnect()
  }

  private async reconnect(): Promise<void> {
    for (let attempt = 0; !this.closing; attempt++) {
      const wait = this.retryMs[Math.min(attempt, this.retryMs.length - 1)]!
      await Bun.sleep(wait)
      if (this.closing) return
      // Still not back: start the next one (the launcher waits for one
      // that is starting or stopping, and starts one only when none is).
      const failed = attempt
      if (this.options.startService && failed >= START_AFTER_ATTEMPTS
          && (failed - START_AFTER_ATTEMPTS) % START_AGAIN_EVERY === 0) {
        await this.options.startService().catch(() => false)
        if (this.closing) return
      }
      try {
        await this.open()
        await this.greet()
        return
      } catch (error) {
        // The socket of a failed attempt is closed here, or every attempt
        // would leave one open on the service.
        const socket = this.socket
        this.socket = null
        try {
          socket?.end()
        } catch {
          // already gone
        }
        this.rpc.lost("still reconnecting")
        if (error instanceof RpcError && error.code === ERRORS.VERSION_MISMATCH) {
          // The service was updated under this window: retrying cannot help.
          this.problem = `the agent service now speaks another protocol (${error.message}); `
            + "start the window again"
          this.closing = true
          this.setState("closed")
          return
        }
      }
    }
  }

  close(): void {
    this.closing = true
    this.setState("closed")
    try {
      this.socket?.end()
    } catch {
      // already gone
    }
  }
}
