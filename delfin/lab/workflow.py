"""Workflows across modules: steps, sample hand-over, resume.

Input: a workflow document and one :class:`EdgeClient` per node. Output: a
run directory ``<lab_dir>/runs/<run_id>/`` with ``checkpoint.json`` (what is
done, what waits, where every sample is) and, when the run ends, the run
record (``record.py``).

Workflow document::

    {"id": "wf-demo", "samples": {"S-001": {"container": "vial-7",
                                            "module": "synthesis",
                                            "position": "S1"}},
     "steps": [
       {"id": "stir", "op": "set", "device": "stirrer", "channel": "power",
        "value": "on"},
       {"id": "heat", "op": "set", "device": "hotplate", "channel":
        "setpoint", "value": 60, "unit": "degC"},
       {"id": "hot", "op": "wait", "device": "hotplate_probe", "channel":
        "temperature", "cmp": ">=", "value": 55, "unit": "degC",
        "timeout_s": 900},
       {"id": "move", "op": "transfer", "sample": "S-001", "from":
        "synthesis", "to": "analysis", "position": "A1"},
       {"id": "uv", "op": "measure", "device": "uvvis", "sample": "S-001",
        "reader": "reader"}]}

Steps run in order; ``after: [ids]`` is checked (a step whose predecessors
are not done is not started), which is the graph form a parallel scheduler
(phase 5) will use. Semantics per op:

* ``set``: one kernel write through the device's edge. ``approval_required``
  pauses the run (``waiting_approval``, the request id in the checkpoint);
  ``resume`` after the person approved repeats the step with that id.
* ``wait``: polls ``read`` until the comparison holds or ``timeout_s``.
* ``transfer``: the registry moves the sample. Refused when the sample is
  not at ``from``, the position is not declared in ``to`` or already holds
  a sample. The physical move is a robot's or a person's; the next step at
  the destination checks the container's barcode.
* ``measure``: reads the reader's barcode at the module; a barcode other
  than the sample's container stops the run (sample mix-up). Then
  ``run=measure``, the signal is recorded as a measurement, ``run=idle``.

The coordinator takes the module lock of every module the workflow uses
for the whole run, so a second coordinator (or an agent session) cannot
interleave writes. A failed step stops the run (``failed``) and keeps the
locks; ``resume`` continues from the failed step, never re-running a
completed one. The edge's own kernel has already put a faulted device in
its safe state; the coordinator does not undo that.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable, Optional

from . import audit as _audit
from . import store as _store
from . import units as _units


def run_dir(run_id: str) -> Path:
    d = _store.ensure_lab_dir() / "runs" / run_id
    d.mkdir(parents=True, exist_ok=True, mode=0o700)
    return d


def _save(run_id: str, cp: dict) -> None:
    path = run_dir(run_id) / "checkpoint.json"
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(cp, indent=1, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def load_checkpoint(run_id: str) -> Optional[dict]:
    path = run_dir(run_id) / "checkpoint.json"
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


class Coordinator:
    def __init__(self, manifest, clients: dict, *, person: str,
                 name: str = "coordinator",
                 clock: Callable[[], float] = time.time,
                 sleep: Callable[[float], None] = time.sleep):
        self.manifest = manifest
        self.clients = clients          # node id -> EdgeClient
        self.person = person
        self.name = name
        self.clock = clock
        self.sleep = sleep

    def _client_for_device(self, device_id: str):
        d = self.manifest.devices.get(device_id)
        if d is None:
            return None, None
        return self.clients.get(d.node), d

    def _modules(self, wf: dict) -> list[str]:
        mods = []
        for st in wf.get("steps") or []:
            for m in (st.get("module"), st.get("from"), st.get("to")):
                if m and m not in mods:
                    mods.append(m)
            dev = self.manifest.devices.get(st.get("device") or "")
            if dev and dev.module not in mods:
                mods.append(dev.module)
        return mods

    def run(self, wf: dict, run_id: str) -> dict:
        actor = {"person": self.person, "session": run_id}
        cp = load_checkpoint(run_id) or {
            "run_id": run_id, "workflow": wf, "status": "new", "done": [],
            "samples": json.loads(json.dumps(wf.get("samples") or {})),
            "log": [], "measurements": [], "deviations": [], "pending": None,
            "started": self.clock(), "coordinator": self.name,
            "person": self.person,
            "manifest": {"version": self.manifest.version,
                         "checksum": self.manifest.checksum}}
        if cp["status"] == "completed":
            return cp
        for mod in self._modules(wf):
            node = self.manifest.modules[mod].node
            res = self.clients[node].call("lock_module", actor, module=mod)
            if not res.get("ok"):
                cp["status"] = "refused"
                cp["deviations"].append({"step": None, "reason":
                                         res.get("reason"), "ts": self.clock()})
                _save(run_id, cp)
                return cp
        cp["status"] = "running"
        _save(run_id, cp)
        for st in wf.get("steps") or []:
            sid = st["id"]
            if sid in cp["done"]:
                continue
            missing = [a for a in st.get("after") or [] if a not in cp["done"]]
            if missing:
                return self._fail(cp, sid, f"predecessors not done: {missing}")
            outcome = self._step(st, actor, cp)
            cp["log"].append({"step": sid, "op": st.get("op"), "ts":
                              self.clock(), "outcome": outcome})
            _audit.record("run_step", run=run_id, step=sid,
                          ok=bool(outcome.get("ok")),
                          code=outcome.get("code"), person=self.person,
                          session=run_id)
            if outcome.get("code") == "approval_required":
                cp["status"] = "waiting_approval"
                cp["pending"] = {"step": sid,
                                 "approval": outcome.get("approval_request")}
                _save(run_id, cp)
                return cp
            if not outcome.get("ok"):
                return self._fail(cp, sid, outcome.get("reason", "failed"))
            cp["done"].append(sid)
            cp["pending"] = None
            _save(run_id, cp)
            self._heartbeat(actor)
        cp["status"] = "completed"
        cp["ended"] = self.clock()
        for mod in self._modules(wf):
            node = self.manifest.modules[mod].node
            self.clients[node].call("unlock_module", actor, module=mod)
        _save(run_id, cp)
        from . import record as _record
        _record.write(cp, self.manifest)
        return cp

    def abandon(self, run_id: str) -> dict:
        """End a failed or waiting run: release its device and module locks
        (a hazardous device goes to its safe state first) and mark it
        ``abandoned``. A person's decision; nothing resumes it."""
        cp = load_checkpoint(run_id)
        if cp is None:
            return {"ok": False, "reason": f"no run {run_id}"}
        actor = {"person": self.person, "session": run_id}
        for client in self.clients.values():
            client.call("release", actor)
        for mod in self._modules(cp.get("workflow") or {}):
            node = self.manifest.modules[mod].node
            self.clients[node].call("unlock_module", actor, module=mod)
        cp["status"] = "abandoned"
        cp["ended"] = self.clock()
        _save(run_id, cp)
        from . import record as _record
        _record.write(cp, self.manifest)
        return cp

    def _fail(self, cp: dict, sid: str, reason: str) -> dict:
        cp["status"] = "failed"
        cp["failed_step"] = sid
        cp["deviations"].append({"step": sid, "reason": reason,
                                 "ts": self.clock()})
        _save(cp["run_id"], cp)
        from . import record as _record
        _record.write(cp, self.manifest)      # a failure is data too
        return cp

    def _heartbeat(self, actor: dict) -> None:
        for c in self.clients.values():
            try:
                c.call("heartbeat", actor)
            except Exception:
                pass

    def _step(self, st: dict, actor: dict, cp: dict) -> dict:
        op = st.get("op")
        if op == "set":
            client, d = self._client_for_device(st.get("device", ""))
            if client is None:
                return {"ok": False, "code": "unknown_device",
                        "reason": f"no edge for device {st.get('device')!r}"}
            approval = None
            pend = cp.get("pending") or {}
            if pend.get("step") == st["id"]:
                approval = pend.get("approval")
            return client.call("set", actor, device=st["device"],
                               channel=st.get("channel"),
                               value=st.get("value"), unit=st.get("unit"),
                               approval=approval)
        if op == "wait":
            client, d = self._client_for_device(st.get("device", ""))
            if client is None:
                return {"ok": False, "reason": "unknown device"}
            deadline = self.clock() + float(st.get("timeout_s", 600))
            while True:
                res = client.call("read", actor, device=st["device"])
                r = (res.get("readings") or {}).get(st.get("channel"))
                if r and r.get("fresh"):
                    want = float(st["value"])
                    if st.get("unit") and r.get("unit"):
                        want = _units.convert(want, st["unit"], r["unit"])
                    v = float(r["value"])
                    cmpop = st.get("cmp", ">=")
                    if {">=": v >= want, "<=": v <= want, ">": v > want,
                            "<": v < want}.get(cmpop, False):
                        return {"ok": True, "value": v, "unit": r.get("unit")}
                if self.clock() >= deadline:
                    return {"ok": False, "code": "timeout", "reason":
                            f"{st['device']}.{st.get('channel')} did not reach "
                            f"{st.get('cmp', '>=')} {st['value']} in time"}
                self._heartbeat(actor)
                self.sleep(float(st.get("poll_s", 5)))
        if op == "transfer":
            sample = cp["samples"].get(st.get("sample"))
            if sample is None:
                return {"ok": False, "code": "sample", "reason":
                        f"unknown sample {st.get('sample')!r}"}
            if sample.get("module") != st.get("from"):
                return {"ok": False, "code": "sample", "reason":
                        f"sample {st['sample']} is in {sample.get('module')}, "
                        f"not {st.get('from')}"}
            dest = self.manifest.modules.get(st.get("to"))
            if dest is None or st.get("position") not in dest.positions:
                return {"ok": False, "code": "sample", "reason":
                        f"position {st.get('position')!r} is not declared in "
                        f"module {st.get('to')!r}"}
            taken = [k for k, s in cp["samples"].items()
                     if s.get("module") == st["to"]
                     and s.get("position") == st["position"]]
            if taken:
                return {"ok": False, "code": "sample", "reason":
                        f"{st['to']}:{st['position']} already holds {taken}"}
            sample["history"] = sample.get("history", []) + [
                {"from": sample.get("module"), "position": sample.get(
                    "position"), "ts": self.clock()}]
            sample["module"], sample["position"] = st["to"], st["position"]
            return {"ok": True, "sample": st["sample"], "at":
                    f"{st['to']}:{st['position']}"}
        if op == "measure":
            sid = st.get("sample")
            sample = cp["samples"].get(sid)
            client, d = self._client_for_device(st.get("device", ""))
            if sample is None or client is None:
                return {"ok": False, "reason": "unknown sample or device"}
            if sample.get("module") != d.module:
                return {"ok": False, "code": "sample", "reason":
                        f"sample {sid} is registered in {sample.get('module')}"
                        f", not at {d.module}"}
            reader = st.get("reader") or st["device"]
            res = client.call("read", actor, device=reader)
            bc = (res.get("readings") or {}).get("barcode") or {}
            if not bc.get("fresh") or bc.get("value") != sample.get("container"):
                return {"ok": False, "code": "sample_mixup", "reason":
                        f"sample mix-up: {d.module} reads container "
                        f"{bc.get('value')!r} (fresh={bc.get('fresh')}), "
                        f"the registry expects {sample.get('container')!r} "
                        f"for {sid}"}
            res = client.call("set", actor, device=st["device"],
                              channel="run", value="measure")
            if not res.get("ok"):
                return res
            got = client.call("read", actor, device=st["device"])
            sig = (got.get("readings") or {}).get("signal") or {}
            client.call("set", actor, device=st["device"], channel="run",
                        value="idle")
            m = {"sample": sid, "container": sample.get("container"),
                 "device": st["device"], "quantity": "signal",
                 "value": sig.get("value"), "unit": sig.get("unit"),
                 "ts": self.clock()}
            cp["measurements"].append(m)
            return {"ok": True, "measurement": m}
        return {"ok": False, "code": "bad_step", "reason": f"unknown op {op!r}"}
