"""The coordinator-edge interface.

An edge node hosts drivers and its own safety kernel; a coordinator plans
and schedules across modules and reaches an edge only through
:meth:`EdgeNode.handle`. Safety holds on the edge alone: its kernel, its
dead-man and its e-stop work with the coordinator gone.

Envelope (JSON)::

    {"v": 1, "node": "edge_synth", "op": "set", "args": {...},
     "actor": {"person": "alice", "session": "run-7", "coordinator": "c1"},
     "manifest": "<sha256 the coordinator planned against>",
     "nonce": "<128-bit hex>", "ts": 1767225600.0,
     "mac": "<HMAC-SHA256 over the canonical body, per-node key>"}

The edge refuses, before any kernel call: another node's envelope; a bad
MAC; a timestamp outside +-``WINDOW_S``; a nonce it has seen (kept for
``2 * WINDOW_S``); a manifest checksum other than the one pinned on the
edge (a manifest altered in transit, or a coordinator planning against a
stale version). The per-node key lives on the edge in
``<lab_dir>/nodes/<node>.key`` (0600) and with the coordinator; an edge's
key verifies nothing on another edge.

Transport. Phase 2 calls ``handle`` in-process (the tests, a single-host
lab). Phase 3 carries the same envelope over mutual TLS (each node and the
coordinator with an X.509 certificate from the installation's CA); the MAC
stays as a second check, so a TLS misconfiguration alone does not open the
kernel. Plain MQTT is not a command channel.

Operations: ``devices``, ``read``, ``set``, ``plan_check``, ``heartbeat``,
``release``, ``lock_module``, ``unlock_module``, ``estop``. ``estop`` is the
one operation accepted without a manifest match: a stop must never be
refused for a version mismatch.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import time
from pathlib import Path
from typing import Callable, Optional

from . import audit as _audit
from . import store as _store

WINDOW_S = 30.0


def _canonical(body: dict) -> bytes:
    return json.dumps(body, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def key_path(node: str) -> Path:
    return _store.lab_dir() / "nodes" / f"{node}.key"


def provision_key(node: str) -> bytes:
    """Create (once) the per-node key; the admin's step on enrolment."""
    p = key_path(node)
    p.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if not p.exists():
        fd = os.open(str(p), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w") as fh:
            fh.write(secrets.token_hex(32))
    return bytes.fromhex(p.read_text().strip())


def load_key(node: str) -> bytes:
    return bytes.fromhex(key_path(node).read_text().strip())


def seal(body: dict, key: bytes) -> dict:
    env = dict(body)
    env["mac"] = hmac.new(key, _canonical(body), hashlib.sha256).hexdigest()
    return env


class EdgeClient:
    """The coordinator's side for one node."""

    def __init__(self, node: str, key: bytes, *, transport: Callable,
                 manifest_checksum: str, coordinator: str = "coordinator",
                 clock: Callable[[], float] = time.time):
        self.node = node
        self.key = key
        self.transport = transport
        self.checksum = manifest_checksum
        self.coordinator = coordinator
        self.clock = clock

    def call(self, op: str, actor: dict, **args) -> dict:
        body = {"v": 1, "node": self.node, "op": op, "args": args,
                "actor": {**actor, "coordinator": self.coordinator},
                "manifest": self.checksum, "nonce": secrets.token_hex(16),
                "ts": self.clock()}
        return self.transport(seal(body, self.key))


class EdgeNode:
    def __init__(self, node: str, *, kernel=None, key: Optional[bytes] = None,
                 clock: Callable[[], float] = time.time):
        from .safety import Kernel
        self.node = node
        self.clock = clock
        self.kernel = kernel or Kernel.load(node=node, clock=clock)
        self.key = key if key is not None else load_key(node)

    def _reject(self, why: str) -> dict:
        _audit.alert("lab_refused", f"edge {self.node}: {why}")
        _audit.record("refused", node=self.node, code="envelope", reason=why)
        return {"ok": False, "code": "envelope", "reason": why}

    def handle(self, env: dict) -> dict:
        from .safety import Actor, pinned
        if not isinstance(env, dict) or env.get("v") != 1:
            return self._reject("not a v1 envelope")
        if env.get("node") != self.node:
            return self._reject(f"envelope for node {env.get('node')!r}")
        body = {k: v for k, v in env.items() if k != "mac"}
        want = hmac.new(self.key, _canonical(body), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(want, str(env.get("mac") or "")):
            return self._reject("bad MAC (not signed with this node's key)")
        now = self.clock()
        try:
            ts = float(env.get("ts"))
        except (TypeError, ValueError):
            return self._reject("no timestamp")
        if abs(now - ts) > WINDOW_S:
            return self._reject(f"timestamp {now - ts:+.0f} s outside the "
                                f"+-{WINDOW_S:g} s window (replay or clock)")
        nonce = str(env.get("nonce") or "")
        if len(nonce) < 16:
            return self._reject("no nonce")
        with _store.state() as doc:
            seen = doc["nonces"].setdefault(self.node, {})
            for n, t in list(seen.items()):
                if now - t > 2 * WINDOW_S:
                    seen.pop(n, None)
            replay = nonce in seen
            if not replay:
                seen[nonce] = now
        if replay:
            return self._reject("nonce already used (replay refused)")
        op = env.get("op")
        args = env.get("args") or {}
        a = env.get("actor") or {}
        actor = Actor(str(a.get("person") or ""), str(a.get("session") or ""),
                      0, f"coordinator:{a.get('coordinator', '')}",
                      str(a.get("coordinator") or ""))
        if not actor.person or not actor.session:
            return self._reject("actor without person or session")
        if op == "estop":
            from . import estop as _estop
            return _estop.trigger(module=args.get("module"),
                                  reason=args.get("reason", "coordinator"),
                                  by=actor.person, kernel=self.kernel)
        if env.get("manifest") != pinned().get("checksum") or (
                self.kernel.manifest is None
                or env.get("manifest") != self.kernel.manifest.checksum):
            _audit.alert("lab_manifest_refused",
                         f"edge {self.node}: coordinator manifest "
                         f"{str(env.get('manifest'))[:12]} is not the pinned one")
            return {"ok": False, "code": "manifest", "reason":
                    "the coordinator's manifest checksum is not the one "
                    "pinned on this edge"}
        k = self.kernel
        if op == "devices":
            return {"ok": True, "devices": k.devices()}
        if op == "read":
            return k.read(args.get("device"), actor)
        if op == "set":
            return k.set(actor, args.get("device"), args.get("channel"),
                         args.get("value"), args.get("unit"),
                         args.get("approval"))
        if op == "plan_check":
            return k.plan_check(actor, args.get("steps") or [])
        if op == "heartbeat":
            return {"ok": True, "renewed": k.heartbeat(actor)}
        if op == "release":
            return {"ok": True, "released": k.release(actor,
                                                      args.get("device"))}
        if op in ("lock_module", "unlock_module"):
            return self._module_lock(op, actor, str(args.get("module") or ""))
        return self._reject(f"unknown op {op!r}")

    def _module_lock(self, op: str, actor, module: str) -> dict:
        m = self.kernel.manifest
        if m is None or module not in m.modules \
                or m.modules[module].node != self.node:
            return {"ok": False, "code": "unknown_module",
                    "reason": f"module {module!r} is not hosted on {self.node}"}
        with _store.state() as doc:
            holder = doc["module_locks"].get(module)
            mine = holder and (holder.get("person"), holder.get(
                "session")) == actor.key()
            if op == "unlock_module":
                if mine:
                    doc["module_locks"].pop(module, None)
                return {"ok": bool(mine)}
            if holder and not mine:
                return {"ok": False, "code": "module_locked", "reason":
                        f"module {module} is held by {holder.get('person')} "
                        f"(run {holder.get('session')}, coordinator "
                        f"{holder.get('coordinator')})"}
            doc["module_locks"][module] = {**actor.as_dict(),
                                           "since": self.clock()}
        return {"ok": True, "module": module}
