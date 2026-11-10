"""The safety kernel: the only path from a request to a driver's write.

Input: an actor (person, session, pid, host) and a request (device,
channel, value, unit, approval id). Output: an outcome dict, ``{"ok":
True, ...}`` with the readback, or ``{"ok": False, "code": ..., "reason":
...}`` with the exact reason. Nothing else in DELFIN calls a driver's
``write``; the tools, the CLI, the edge interface and the workflow runner
all come through :meth:`Kernel.set`.

Order of checks for a write (each refusal is audited and returned):

1. e-stop (global or the device's module) and a latched device;
2. the manifest is pinned and unchanged (``pin``);
3. device known, hosted on this node, channel writable;
4. level: ``forbidden`` is never executed;
5. value: a state from the declared list; a number with a unit of the
   right dimension (converted), inside ``min``/``max``;
6. dead-man sweep (an owner that died is put safe before anyone else acts);
7. lock: held by another (person, session) means refused, never taken;
   a module locked by another actor (a coordinator's run) means refused;
8. rate: ``min_interval_s`` since the last accepted write of the channel,
   ``max_step`` against the current readback (unknown means refused);
9. interlocks: every interlock naming this channel, evaluated with the
   proposed value and FRESH readings of everything else. A reading older
   than its device's ``stale_after_s`` is UNKNOWN; UNKNOWN fails closed;
10. approval for ``hazardous``: a single-use grant bound to this action;
11. one transaction: the lock taken, the approval consumed, the dead-man
    obligation recorded; then ``intent`` audited, the driver called,
    ``result`` audited, the device read back and ``readback`` audited. A
    readback outside ``tolerance`` puts the device in its safe state,
    latches it and alerts.

The kernel does not read the agent's permission mode. That is the point:
bypass, all_free, an accepted dialog, a prompt or a file the model read
change nothing here.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional

from . import approvals as _approvals
from . import audit as _audit
from . import expr as _expr
from . import manifest as _manifest
from . import store as _store
from . import units as _units
from .drivers import make as _make_driver


@dataclass(frozen=True)
class Actor:
    person: str
    session: str
    pid: int = 0
    host: str = ""
    coordinator: str = ""

    @classmethod
    def here(cls, session: str, person: Optional[str] = None) -> "Actor":
        return cls(person or _store.person(), str(session), os.getpid(),
                   _store.host())

    def key(self) -> tuple:
        return (self.person, self.session)

    def as_dict(self) -> dict:
        return {"person": self.person, "session": self.session,
                "pid": self.pid, "host": self.host,
                "coordinator": self.coordinator}


def _same(holder: dict, actor: Actor) -> bool:
    return (holder.get("person"), holder.get("session")) == actor.key()


def refuse(code: str, reason: str, **extra) -> dict:
    out = {"ok": False, "code": code, "reason": reason}
    out.update(extra)
    return out


# -- manifest pin -------------------------------------------------------------

def pin_path() -> Path:
    return _store.lab_dir() / "manifest.pin.json"


def pin(manifest: _manifest.Manifest, *, by: Optional[str] = None) -> dict:
    """Pin the manifest's checksum: the person's step after reviewing it."""
    from delfin.agent.approval_guard import inside_agent_session
    if inside_agent_session():
        raise PermissionError("a manifest is pinned by a person, not by a "
                              "process an agent session started")
    rec = {"checksum": manifest.checksum, "version": manifest.version,
           "installation": manifest.installation,
           "by": by or _store.person(), "ts": _store.now()}
    with _store.locked_json("manifest.pin.json") as doc:
        doc.clear()
        doc.update(rec)
    _audit.record("manifest_pinned", **rec)
    return rec


def pinned() -> dict:
    return _store.read_json(pin_path(), {})


# -- the kernel -----------------------------------------------------------------

class Kernel:
    def __init__(self, manifest: Optional[_manifest.Manifest], *,
                 node: Optional[str] = None,
                 clock: Callable[[], float] = time.time,
                 load_error: str = ""):
        self.manifest = manifest
        self.node = node
        self.clock = clock
        self.load_error = load_error
        self._drivers: dict = {}

    @classmethod
    def load(cls, *, node: Optional[str] = None,
             clock: Callable[[], float] = time.time) -> "Kernel":
        path = _store.manifest_path()
        try:
            m = _manifest.load(path)
            return cls(m, node=node, clock=clock)
        except _manifest.ManifestError as exc:
            return cls(None, node=node, clock=clock,
                       load_error=f"manifest refused: {exc}")
        except Exception as exc:
            return cls(None, node=node, clock=clock,
                       load_error=f"manifest unreadable: {exc}")

    # -- helpers ---------------------------------------------------------------

    def driver(self, device):
        drv = self._drivers.get(device.driver)
        if drv is None:
            drv = _make_driver(device.driver)
            if hasattr(drv, "clock"):
                drv.clock = self.clock
            self._drivers[device.driver] = drv
        return drv

    def _manifest_problem(self) -> Optional[str]:
        if self.manifest is None:
            return self.load_error or "no lab manifest on this installation"
        p = pinned()
        if not p.get("checksum"):
            return ("the manifest is not pinned; a person reviews it and runs "
                    "`delfin-agent lab manifest pin`")
        if p["checksum"] != self.manifest.checksum:
            return ("the manifest changed since it was pinned (checksum "
                    f"{self.manifest.checksum[:12]} != pinned "
                    f"{p['checksum'][:12]}); every write is refused until a "
                    "person reviews and pins it again")
        return None

    def _device(self, device_id: str):
        m = self.manifest
        d = m.devices.get(str(device_id or "")) if m else None
        if d is None:
            known = ", ".join(sorted(m.devices)) if m else ""
            return None, refuse("unknown_device",
                                f"no device {device_id!r} in the manifest; "
                                f"known: {known}")
        if self.node and d.node != self.node:
            return None, refuse("wrong_node", f"{d.id} is hosted on node "
                                f"{d.node}, not {self.node}")
        return d, None

    def _fresh(self, device) -> dict:
        """``{channel: (value, unit)}`` or UNKNOWN per channel, read now."""
        try:
            got = self.driver(device).read(device)
        except Exception as exc:
            return {"__error__": str(exc)}
        now = self.clock()
        out = {}
        for ch, r in got.items():
            age = now - float(r.ts)
            out[ch] = _expr.UNKNOWN if age > device.stale_after_s \
                else (r.value, r.unit)
            out.setdefault("__age__", {})[ch] = round(age, 3)
        return out

    def _refused(self, actor: Optional[Actor], device_id: str, out: dict,
                 **ctx) -> dict:
        _audit.record("refused", device=device_id, code=out["code"],
                      reason=out["reason"],
                      **(actor.as_dict() if actor else {}), **ctx)
        if out["code"] not in ("approval_required",):
            _audit.alert("lab_refused", f"{device_id}: {out['reason']}")
        return out

    def _normalise(self, device, channel, value, unit):
        """(value in the channel's unit, error dict or None)."""
        if channel.is_state:
            v = _manifest.state_str(value)
            if v not in channel.states:
                return None, refuse("bad_state", f"{device.id}.{channel.name} "
                                    f"takes one of {', '.join(channel.states)}"
                                    f", not {value!r}")
            return v, None
        if isinstance(value, bool):
            return None, refuse("bad_value", "a number is required")
        try:
            num = float(value)
        except (TypeError, ValueError):
            return None, refuse("bad_value", f"{value!r} is not a number")
        if num != num or num in (float("inf"), float("-inf")):
            return None, refuse("bad_value", "not a finite number")
        try:
            got = _units.canonical(unit)
            num = _units.convert(num, got, channel.unit)
        except _units.UnitError as exc:
            return None, refuse("bad_unit", f"{device.id}.{channel.name} is "
                                f"in {channel.unit}: {exc}")
        if not channel.min <= num <= channel.max:
            return None, refuse("out_of_range",
                                f"{device.id}.{channel.name} = {num:g} "
                                f"{channel.unit} is outside the declared "
                                f"range {channel.min:g}..{channel.max:g} "
                                f"{channel.unit}")
        return num, None

    def _interlocks(self, device, channel, proposed, overrides=None):
        """None when every interlock naming the channel holds, else refusal."""
        target = _expr.Ref(device.id, channel.name)
        overrides = dict(overrides or {})
        overrides[target] = (proposed, channel.unit)
        cache: dict = {}
        unknown_why: dict = {}

        def lookup(ref):
            if ref in overrides:
                return overrides[ref]
            dev = self.manifest.devices.get(ref.device)
            if dev is None:
                return _expr.UNKNOWN
            if ref.device not in cache:
                cache[ref.device] = self._fresh(dev)
            got = cache[ref.device]
            if "__error__" in got:
                unknown_why[ref.name] = f"unreadable ({got['__error__']})"
                return _expr.UNKNOWN
            val = got.get(ref.channel, _expr.UNKNOWN)
            if val is _expr.UNKNOWN:
                age = got.get("__age__", {}).get(ref.channel)
                unknown_why[ref.name] = (
                    f"reading is {age:g} s old (stale after "
                    f"{dev.stale_after_s:g} s)" if age is not None
                    else "no reading")
            return val

        for il in self.manifest.interlocks():
            when, req = _expr.parse(il.when), _expr.parse(il.requires)
            if target not in _expr.refs(when) | _expr.refs(req):
                continue
            w = _expr.evaluate(when, lookup)
            if w is False:
                continue
            r = _expr.evaluate(req, lookup)
            if r is True:
                continue
            if r is _expr.UNKNOWN or w is _expr.UNKNOWN:
                why = "; ".join(f"{k}: {v}" for k, v in sorted(
                    unknown_why.items())) or "a reading is unknown"
                return refuse("interlock_unknown",
                              f"interlock '{il.text()}' cannot be confirmed "
                              f"({why}); unknown fails closed")
            seen = []
            for ref in sorted(_expr.refs(req), key=lambda x: x.name):
                v = lookup(ref)
                seen.append(f"{ref.name}="
                            f"{'unknown' if v is _expr.UNKNOWN else v[0]}")
            return refuse("interlock", f"interlock '{il.text()}' is not met "
                          f"({', '.join(seen)})")
        return None

    @staticmethod
    def _is_safe_value(device, channel, value) -> bool:
        if channel.name not in device.safe_state:
            return False
        sv = device.safe_state[channel.name]
        if channel.is_state:
            return str(value) == str(sv)
        return abs(float(value) - float(sv)) <= 1e-9

    def _hazardous(self, device) -> bool:
        return any(self.manifest.level(device, c) == "hazardous"
                   for c in device.writes.values())

    def _at_safe(self, device, values: dict) -> bool:
        for ch, sv in device.safe_state.items():
            v = values.get(ch)
            if v is _expr.UNKNOWN or v is None:
                return False
            v = v[0]
            spec = device.writes[ch]
            if spec.is_state:
                if str(v) != str(sv):
                    return False
            elif abs(float(v) - float(sv)) > max(spec.tolerance, 1e-9):
                return False
        return True

    # -- public ----------------------------------------------------------------

    def devices(self) -> list[dict]:
        if self.manifest is None:
            return []
        snap = _store.snapshot()
        out = []
        for d in self.manifest.devices.values():
            if self.node and d.node != self.node:
                continue
            fresh = self._fresh(d)
            state = {}
            for ch, v in fresh.items():
                if ch.startswith("__"):
                    continue
                state[ch] = "unknown (stale)" if v is _expr.UNKNOWN else v[0]
            out.append({
                "id": d.id, "module": d.module,
                "reads": {k: c.unit or ("text" if not c.states else
                                        list(c.states))
                          for k, c in d.reads.items()},
                "writes": {k: ({"states": list(c.states)} if c.is_state else
                               {"unit": c.unit, "min": c.min, "max": c.max})
                           | {"level": self.manifest.level(d, c)}
                           for k, c in d.writes.items()},
                "state": state if "__error__" not in fresh
                else {"error": fresh["__error__"]},
                "lock": (snap["locks"].get(d.id) or {}).get("person"),
                "latched": (snap["latched"].get(d.id) or {}).get("reason"),
            })
        return out

    def estopped(self, device=None, snap=None) -> Optional[str]:
        snap = snap or _store.snapshot()
        g = snap["estop"].get("global")
        if g:
            return f"global e-stop by {g.get('by')}: {g.get('reason', '')}"
        if device is not None:
            m = snap["estop"]["modules"].get(device.module)
            if m:
                return (f"e-stop of module {device.module} by {m.get('by')}: "
                        f"{m.get('reason', '')}")
        return None

    def read(self, device_id: str, actor: Optional[Actor] = None) -> dict:
        problem = None if self.manifest else (
            self.load_error or "no lab manifest on this installation")
        if problem:
            return refuse("no_manifest", problem)
        d, err = self._device(device_id)
        if err:
            return err
        if actor is not None:
            self.heartbeat(actor)
        try:
            got = self.driver(d).read(d)
        except Exception as exc:
            return refuse("driver_error", f"{d.id}: {exc}")
        now = self.clock()
        readings = {}
        for ch, r in got.items():
            age = now - float(r.ts)
            readings[ch] = {"value": r.value, "unit": r.unit,
                            "age_s": round(age, 3),
                            "fresh": age <= d.stale_after_s}
        return {"ok": True, "device": d.id, "readings": readings}

    def set(self, actor: Actor, device_id: str, channel_name: str, value: Any,
            unit: Optional[str] = None, approval: Optional[str] = None) -> dict:
        now = self.clock()
        problem = self._manifest_problem()
        if problem:
            _audit.alert("lab_manifest_refused", problem)
            return self._refused(actor, str(device_id),
                                 refuse("manifest", problem))
        d, err = self._device(device_id)
        if err:
            return self._refused(actor, str(device_id), err)
        snap = _store.snapshot()
        stop = self.estopped(d, snap)
        if stop:
            return self._refused(actor, d.id, refuse(
                "estop", f"{stop}; every write is refused until a person "
                "runs `delfin-agent lab estop --clear`"))
        latched = snap["latched"].get(d.id)
        if latched:
            return self._refused(actor, d.id, refuse(
                "latched", f"{d.id} was put in its safe state "
                f"({latched.get('reason')}); a person clears it with "
                f"`delfin-agent lab clear {d.id}`"))
        ch = d.writes.get(str(channel_name or ""))
        if ch is None:
            return self._refused(actor, d.id, refuse(
                "not_writable", f"{d.id} has no writable channel "
                f"{channel_name!r}; writable: {', '.join(sorted(d.writes))}"))
        level = self.manifest.level(d, ch)
        if level == "forbidden":
            return self._refused(actor, d.id, refuse(
                "forbidden", f"{d.id}.{ch.name} is forbidden in the manifest; "
                "it is never executed, with or without approval"))
        num, err = self._normalise(d, ch, value, unit)
        if err:
            return self._refused(actor, d.id, err)
        if level == "hazardous" and self._is_safe_value(d, ch, num):
            # A write to the channel's own safe value moves the device toward
            # its safe state; it is audited like any write, never held for a
            # person (a cool-down must not wait for an approval).
            level = "benign"
        # 6. dead-man sweep before judging locks.
        from . import watchdog as _watchdog
        _watchdog.sweep(self, now=now)
        snap = _store.snapshot()
        if snap["latched"].get(d.id):
            return self._refused(actor, d.id, refuse(
                "latched", f"{d.id} was just put in its safe state by the "
                "dead-man; a person clears it"))
        holder = snap["locks"].get(d.id)
        if holder and not _same(holder, actor):
            return self._refused(actor, d.id, refuse(
                "locked", f"{d.id} is controlled by {holder.get('person')} "
                f"(session {holder.get('session')}); one session controls a "
                "device at a time, and a lock is never taken over"))
        mlock = snap["module_locks"].get(d.module)
        if mlock and not _same(mlock, actor):
            return self._refused(actor, d.id, refuse(
                "module_locked", f"module {d.module} is held by "
                f"{mlock.get('person')} (run {mlock.get('session')})"))
        # 8. rate.
        key = f"{d.id}.{ch.name}"
        last = snap["last_set"].get(key)
        if last and ch.min_interval_s and now - last["ts"] < ch.min_interval_s:
            return self._refused(actor, d.id, refuse(
                "rate", f"{key}: at most one write per {ch.min_interval_s:g} s"))
        fresh = self._fresh(d)
        if ch.max_step is not None and not ch.is_state:
            cur = fresh.get(ch.readback or ch.name, _expr.UNKNOWN)
            if cur is _expr.UNKNOWN or "__error__" in fresh:
                return self._refused(actor, d.id, refuse(
                    "rate", f"{key}: the current value is unknown, so the "
                    "step size cannot be checked"))
            if abs(num - float(cur[0])) > ch.max_step:
                return self._refused(actor, d.id, refuse(
                    "rate", f"{key}: a step from {float(cur[0]):g} to "
                    f"{num:g} {ch.unit} exceeds max_step {ch.max_step:g}"))
        # 9. interlocks.
        err = self._interlocks(d, ch, num)
        if err:
            return self._refused(actor, d.id, err)
        # 10.-11. approval, lock, obligation: one transaction.
        dig = _approvals.digest(installation=self.manifest.installation,
                                checksum=self.manifest.checksum, device=d.id,
                                channel=ch.name, value=num,
                                person=actor.person, session=actor.session)
        with _store.state() as doc:
            holder = doc["locks"].get(d.id)
            if holder and not _same(holder, actor):
                doc_out = refuse("locked", f"{d.id} was taken by "
                                 f"{holder.get('person')} a moment ago")
            elif self.estopped(d, doc):
                doc_out = refuse("estop", "e-stop arrived during the check")
            elif level == "hazardous" and not approval:
                rec = _approvals.file_request(
                    doc, digest_=dig, person=actor.person,
                    session=actor.session, device=d.id, channel=ch.name,
                    value=num, unit=ch.unit, now=now)
                doc_out = refuse(
                    "approval_required",
                    f"{key} = {num if ch.is_state else format(num, 'g')}{' ' + ch.unit if ch.unit else ''} is "
                    f"hazardous: a person approves exactly this action with "
                    f"`delfin-agent lab approve {rec['id']}`, then repeat the "
                    f"call with approval={rec['id']}",
                    approval_request=rec["id"])
            else:
                doc_out = None
                if level == "hazardous":
                    why = _approvals.consume(doc, approval, digest_=dig,
                                             person=actor.person, now=now)
                    if why:
                        doc_out = refuse("approval", why)
                if doc_out is None:
                    lock = doc["locks"].get(d.id) or {
                        **actor.as_dict(), "since": now}
                    lock["heartbeat"] = now
                    doc["locks"][d.id] = lock
                    if self._hazardous(d):
                        doc["active"][d.id] = {**actor.as_dict(),
                                               "since": now, "heartbeat": now}
                    doc["last_set"][key] = {"value": num, "ts": now}
        if doc_out is not None:
            if doc_out["code"] == "approval_required":
                _audit.record("approval_requested",
                              id=doc_out["approval_request"], device=d.id,
                              channel=ch.name, value=num, unit=ch.unit,
                              **actor.as_dict())
            elif doc_out["code"] == "approval":
                _audit.alert("lab_approval_refused", doc_out["reason"])
            return self._refused(actor, d.id, doc_out)
        if level == "hazardous":
            _audit.record("approval_consumed", id=approval, device=d.id,
                          **actor.as_dict())
        drv = self.driver(d)
        if (self._hazardous(d) and d.device_watchdog
                and drv.capabilities.device_watchdog):
            try:
                drv.arm_watchdog(d, d.deadman_s)
            except Exception:
                pass
        _audit.record("intent", device=d.id, channel=ch.name, value=num,
                      unit=ch.unit, level=level, approval=approval,
                      **actor.as_dict())
        try:
            cmd = drv.write(d, ch.name, num)
        except Exception as exc:
            _audit.record("result", device=d.id, channel=ch.name, ok=False,
                          error=str(exc), **actor.as_dict())
            self.safe(d.id, f"write failed: {exc}", latch=True)
            return refuse("driver_error", f"{d.id}: {exc}; the device was put "
                          "in its safe state")
        _audit.record("result", device=d.id, channel=ch.name, ok=True,
                      command=cmd, **actor.as_dict())
        # An e-stop that arrived while the driver was writing: the stop's own
        # safe-state pass may have run before this write landed, so the
        # write is undone here. The flag is set before that pass starts, so
        # every write that began before it is seen by this check.
        late = self.estopped(d)
        if late:
            self.safe(d.id, f"write landed during an e-stop ({late})",
                      latch=False)
            return refuse("estop", f"{late}; the write was undone")
        after = self._fresh(d)
        rb = after.get(ch.readback or ch.name, _expr.UNKNOWN)
        good = rb is not _expr.UNKNOWN and "__error__" not in after and (
            str(rb[0]) == str(num) if ch.is_state
            else abs(float(rb[0]) - num) <= ch.tolerance + 1e-9)
        _audit.record("readback", device=d.id, channel=ch.name,
                      expected=num, got=None if rb is _expr.UNKNOWN else rb[0],
                      ok=good, **actor.as_dict())
        if not good:
            got = "unknown" if rb is _expr.UNKNOWN else rb[0]
            reason = (f"readback mismatch on {key}: commanded {num}, device "
                      f"reports {got}")
            _audit.alert("lab_readback_mismatch", reason)
            self.safe(d.id, reason, latch=True)
            return refuse("readback_mismatch", reason + "; the device was put "
                          "in its safe state and latched")
        if self._at_safe(d, after):
            with _store.state() as doc:
                doc["active"].pop(d.id, None)
            try:
                drv.disarm_watchdog(d)
            except Exception:
                pass
        return {"ok": True, "device": d.id, "channel": ch.name, "value": num,
                "unit": ch.unit, "readback": rb[0], "command": cmd,
                "level": level}

    def heartbeat(self, actor: Actor) -> int:
        """Renew every lock and dead-man obligation of ``actor``."""
        now = self.clock()
        rearm = []
        with _store.state() as doc:
            n = 0
            for dev, lock in doc["locks"].items():
                if _same(lock, actor):
                    lock["heartbeat"] = now
                    n += 1
            for dev, act in doc["active"].items():
                if _same(act, actor):
                    act["heartbeat"] = now
                    rearm.append(dev)
        if self.manifest is not None:
            for dev in rearm:
                d = self.manifest.devices.get(dev)
                if d is None or (self.node and d.node != self.node):
                    continue
                drv = self.driver(d)
                if d.device_watchdog and drv.capabilities.device_watchdog:
                    try:
                        drv.arm_watchdog(d, d.deadman_s)
                    except Exception:
                        pass
        return n

    def release(self, actor: Actor, device_id: Optional[str] = None) -> list:
        """Give up locks of ``actor``; a hazardous device goes safe first."""
        out = []
        snap = _store.snapshot()
        for dev, lock in list(snap["locks"].items()):
            if not _same(lock, actor) or (device_id and dev != device_id):
                continue
            if dev in snap["active"]:
                res = self.safe(dev, f"released by {actor.person}", latch=False)
                if not res.get("ok"):
                    # Still hazardous and the drive failed: keep the lock and
                    # the active guard so the dead-man keeps watching it, and
                    # do not report the device as released.
                    _audit.alert("lab_release",
                                 f"{dev}: not released; safe state does not "
                                 "hold, the guard stays")
                    continue
            with _store.state() as doc:
                doc["locks"].pop(dev, None)
                doc["active"].pop(dev, None)
            _audit.record("released", device=dev, **actor.as_dict())
            out.append(dev)
        return out

    def safe(self, device_id: str, reason: str, *, latch: bool) -> dict:
        """Drive one device to its manifest safe state and verify it."""
        d = self.manifest.devices.get(device_id) if self.manifest else None
        if d is None:
            return {"ok": False, "device": device_id, "reason": "unknown"}
        drv = self.driver(d)
        try:
            failures = drv.safe_state(d)
        except Exception as exc:
            failures = [str(exc)]
        try:
            drv.disarm_watchdog(d)
        except Exception:
            pass
        after = self._fresh(d)
        confirmed = "__error__" not in after and self._at_safe(d, after)
        with _store.state() as doc:
            if confirmed:
                # The hazardous obligation ends only when the readback
                # confirms the device is safe; otherwise the dead-man and
                # the release guard keep watching it.
                doc["active"].pop(d.id, None)
            if latch:
                doc["latched"][d.id] = {"reason": reason, "ts": self.clock()}
                doc["locks"].pop(d.id, None)
        _audit.record("safe_state", device=d.id, reason=reason,
                      confirmed=confirmed, failures=failures or None)
        if latch:
            _audit.record("latched", device=d.id, reason=reason)
        _audit.alert("lab_safe_state", f"{d.id}: {reason}"
                     + ("" if confirmed else " -- safe state NOT confirmed by "
                        "readback"))
        return {"ok": confirmed, "device": d.id, "failures": failures}

    def plan_check(self, actor: Optional[Actor], steps: list) -> dict:
        """Dry run: every step judged as the kernel would, nothing written."""
        if self.manifest is None:
            return refuse("no_manifest", self.load_error or "no lab manifest")
        problem = self._manifest_problem()
        overrides: dict = {}
        results = []
        ok_all = problem is None
        for i, step in enumerate(steps or []):
            step = step if isinstance(step, dict) else {}
            dev_id = step.get("device")
            d, err = self._device(dev_id)
            verdict: dict
            if err:
                verdict = err
            else:
                ch = d.writes.get(str(step.get("channel") or ""))
                if ch is None:
                    verdict = refuse("not_writable", f"{d.id} has no writable "
                                     f"channel {step.get('channel')!r}")
                elif self.manifest.level(d, ch) == "forbidden":
                    verdict = refuse("forbidden", f"{d.id}.{ch.name} is "
                                     "forbidden")
                else:
                    num, err = self._normalise(d, ch, step.get("value"),
                                               step.get("unit"))
                    if err:
                        verdict = err
                    else:
                        err = self._interlocks(d, ch, num, overrides)
                        if err:
                            verdict = err
                        else:
                            overrides[_expr.Ref(d.id, ch.name)] = (num, ch.unit)
                            level = self.manifest.level(d, ch)
                            verdict = {"ok": True, "level": level,
                                       "needs_approval": level == "hazardous",
                                       "value": num, "unit": ch.unit}
                stop = self.estopped(d)
                if stop and verdict.get("ok"):
                    verdict = refuse("estop", stop)
            verdict = {"step": i, **verdict}
            ok_all = ok_all and bool(verdict.get("ok"))
            results.append(verdict)
        return {"ok": ok_all, "manifest_problem": problem, "steps": results,
                "touched_hardware": False}
