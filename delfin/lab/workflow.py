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
import logging
import time
from pathlib import Path
from typing import Callable, Optional

from delfin.agent import state_file

from . import audit as _audit
from . import expr as _expr
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
    """The run's checkpoint, or None when there is none.

    Damaged and absent are NOT the same: a damaged file is kept aside and a
    :class:`~delfin.agent.state_file.StateFileUnreadable` is raised, so a
    resume cannot silently rebuild a run's history over it (samples are on
    the bench; the file is the only record of where).  A caller that only
    reads -- ``abandon`` -- catches this and answers, it does not write.
    """
    path = run_dir(run_id) / "checkpoint.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        aside = state_file.quarantine(path)
        _report_checkpoint_loss(run_id, exc, aside)
        raise state_file.StateFileUnreadable(
            f"checkpoint.json of run {run_id} is unreadable ({exc}); the "
            f"file was kept at {aside}") from None


def _report_checkpoint_loss(run_id, exc, aside) -> None:
    logging.getLogger(__name__).warning(
        "The checkpoint of lab run %s is unreadable (%s); it was moved to %s "
        "and kept -- a checkpoint (done steps, sample positions, approvals) "
        "is not reconstructable from the workflow document.",
        run_id, exc, aside)


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

    def _carrier_for(self, module_id: str, position: str):
        """(client, device) of the device that HOLDS a container at ``position``.

        A hand-over has to be made true of the world, and the only device that
        can make it true is one that (a) declares this position in its
        ``address`` and (b) has a writable ``container`` slot -- a rack slot.
        Returns ``(None, None)`` when the module has no such device, which the
        caller reports as a registry-only hand-over rather than silently
        claiming a move. More than one such device is a mis-wired manifest, and
        the first in document order is used so the outcome is deterministic
        rather than dependent on dict order at run time.
        """
        mod = self.manifest.modules.get(module_id)
        if mod is None:
            return None, None
        for dev_id in mod.devices:
            dev = self.manifest.devices.get(dev_id)
            if dev is None or (dev.address or {}).get("position") != position:
                continue
            if "container" not in dev.writes:
                continue
            client, d = self._client_for_device(dev_id)
            if client is not None:
                return client, d
        return None, None

    def run(self, wf: dict, run_id: str, *, parallel: int = 1) -> dict:
        """Run (or resume) a workflow.

        ``parallel=1`` (the default) keeps the original contract: the steps in
        document order, one at a time. ``parallel=N`` runs up to N steps at
        once, with two hard rules:

        * a step starts only when every ``after`` predecessor is done, and
        * two steps never hold one device at the same time.

        The second rule is what keeps the kernel's actor model intact. The
        kernel takes a per-device lock for the actor that writes, and refuses a
        second actor, so concurrent steps must be given different devices --
        which is exactly the shape of a parallel campaign (four reactors, four
        pumps, four instruments).
        """
        actor = {"person": self.person, "session": run_id}
        try:
            cp = load_checkpoint(run_id)
        except state_file.StateFileUnreadable as exc:
            # A damaged checkpoint must never be rebuilt over: it holds the
            # physical truth of the run (done steps, sample positions).  The
            # damaged file is already aside; the failure is recorded so the
            # run does not claim 'completed' over a history it lost.
            logging.getLogger(__name__).error(
                "Refusing to resume lab run %s over a damaged checkpoint: %s",
                run_id, exc)
            failed = {
                "run_id": run_id, "workflow": wf, "status": "failed",
                "done": [], "samples": {}, "log": [], "measurements": [],
                "deviations": [{"step": None, "reason": str(exc),
                                "ts": self.clock()}],
                "pending": None, "started": self.clock(),
                "coordinator": self.name, "person": self.person,
                "manifest": {"version": self.manifest.version,
                             "checksum": self.manifest.checksum},
            }
            _save(run_id, failed)
            from . import record as _record
            _record.write(failed, self.manifest)   # a loss is data too
            return failed
        cp = cp or {
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
        if parallel and parallel > 1:
            return self._run_parallel(wf, cp, actor, int(parallel))
        for st in wf.get("steps") or []:
            sid = st["id"]
            if sid in cp["done"]:
                continue
            missing = [a for a in st.get("after") or [] if a not in cp["done"]]
            if missing:
                return self._fail(cp, sid, f"predecessors not done: {missing}")
            out = self._dispatch_step(st, cp, actor)
            if out is not None:
                return out
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

    def _dispatch_step(self, st: dict, cp: dict, actor: dict):
        """Run one step and fold its outcome into the checkpoint.

        Returns ``None`` when the run may go on, or the finished checkpoint
        when the run stops here (a failed step, or a step waiting for a
        person's approval). The sequential loop and the parallel scheduler
        both come through this function, so they cannot drift apart on what
        counts as done, what a refusal does, or what is written down.
        """
        sid = st["id"]
        outcome = self._step(st, actor, cp)
        cp["log"].append({"step": sid, "op": st.get("op"), "ts":
                          self.clock(), "outcome": outcome})
        _audit.record("run_step", run=cp["run_id"], step=sid,
                      ok=bool(outcome.get("ok")),
                      code=outcome.get("code"), person=self.person,
                      session=cp["run_id"])
        if outcome.get("code") == "approval_required":
            cp["status"] = "waiting_approval"
            cp["pending"] = {"step": sid,
                             "approval": outcome.get("approval_request")}
            _save(cp["run_id"], cp)
            return cp
        if not outcome.get("ok"):
            return self._fail(cp, sid, outcome.get("reason", "failed"))
        cp["done"].append(sid)
        # Only this step's own request may be cleared: in a parallel run another
        # worker can file a request while this one succeeds, and clearing it
        # would leave the run waiting on an approval id nobody can read.
        if (cp.get("pending") or {}).get("step") == sid:
            cp["pending"] = None
        _save(cp["run_id"], cp)
        return None

    def _step_devices(self, st: dict) -> set:
        """Every device a step will touch: the device, its reader and, for a
        ``set``, the devices its interlocks may read.

        Two steps that share a device must not run at once: the kernel takes
        one actor's lock per device and refuses the next, so a scheduler that
        ignored this would turn a legitimate campaign into a stream of
        ``locked`` refusals. Interlock references are included because reading
        a device is cheap but a write to one the other step is reading is the
        race worth preventing.
        """
        out = set()
        for key in ("device", "reader"):
            if st.get(key):
                out.add(str(st[key]))
        dev = self.manifest.devices.get(str(st.get("device") or ""))
        if dev is not None and st.get("op") == "set":
            ch = dev.writes.get(str(st.get("channel") or ""))
            if ch is not None:
                for il in dev.interlocks:
                    for text in (il.when, il.requires):
                        try:
                            tree = _expr.parse(text)
                        except (_expr.ExprError, _units.UnitError):
                            continue
                        out.update(r.device for r in _expr.refs(tree))
        return out

    def _run_parallel(self, wf: dict, cp: dict, actor: dict, workers: int):
        """Steps whose predecessors are done, up to ``workers`` at a time, and
        never two on one device."""
        steps = [st for st in (wf.get("steps") or [])
                 if st.get("id") not in cp["done"]]
        if not steps:
            return self._finish(cp)
        import threading

        lock = threading.Lock()          # guards the checkpoint
        # Writes are serialised; dwells and waits are not. A kernel write is a
        # handful of milliseconds (and on real hardware the bus is one-at-a-
        # time anyway), while the time a campaign spends is in the dwells and
        # the waits. Overlapping the first buys nothing and races the driver's
        # own store; overlapping the second is the whole point. A 4-lane
        # campaign therefore costs one lane's dwells plus 4x its writes.
        # A re-entrant lock: a worker holds it across a write, and the
        # heartbeat it takes on the way in must not deadlock against it.
        write_lock = threading.RLock()
        done = set(cp["done"])
        busy: set = set()                # devices a live worker holds
        stop_state: list = []            # the checkpoint the run ended on
        results: dict = {}               # step id -> outcome, None while live
        live: list = []                  # threads this scheduler started

        def worker(st):
            # The grant this step was resumed with. Read from the per-step map
            # and not from `cp["pending"]`: in a parallel pass several steps
            # can park, and `pending` keeps only the last one -- so a worker
            # resumed with its own grant would otherwise hand the kernel a
            # request id belonging to a different step.
            grant = (cp.get("waiting") or {}).get(st["id"])
            try:
                if st.get("op") in ("set", "measure", "transfer"):
                    with write_lock:
                        outcome = self._step(st, actor, cp, approval=grant)
                else:
                    # hold / wait: the part that takes the time
                    outcome = self._step(st, actor, cp)
            except Exception as exc:     # a driver error is a step failure
                outcome = {"ok": False, "code": "error", "reason": str(exc)}
            with lock:
                results[st["id"]] = outcome
                sid = st["id"]
                cp["log"].append({"step": sid, "op": st.get("op"), "ts":
                                  self.clock(), "outcome": outcome})
                _audit.record("run_step", run=cp["run_id"], step=sid,
                              ok=bool(outcome.get("ok")),
                              code=outcome.get("code"), person=self.person,
                              session=cp["run_id"])
                if outcome.get("code") == "approval_required":
                    cp["status"] = "waiting_approval"
                    cp["pending"] = {"step": sid, "approval":
                                     outcome.get("approval_request")}
                    # Every parked step keeps its own id, so a resume hands
                    # each worker the grant that was filed for it.
                    cp.setdefault("waiting", {})[sid] = \
                        outcome.get("approval_request")
                    _save(cp["run_id"], cp)
                    stop_state.append(cp)
                elif not outcome.get("ok"):
                    stop_state.append(self._fail(
                        cp, sid, outcome.get("reason", "failed")))
                else:
                    cp["done"].append(sid)
                    if (cp.get("pending") or {}).get("step") == sid:
                        cp["pending"] = None
                    (cp.get("waiting") or {}).pop(sid, None)
                    done.add(sid)
                    _save(cp["run_id"], cp)
                busy.difference_update(self._step_devices(st))
                for i, th in enumerate(live):
                    if th.ident == threading.get_ident():
                        del live[i]
                        break

        def startable():
            """The next step whose predecessors are done and whose devices are
            free, in document order."""
            for st in steps:
                if st["id"] in done or st["id"] in results:
                    continue
                if any(a not in done for a in st.get("after") or []):
                    continue
                if self._step_devices(st) & busy:
                    continue
                return st
            return None

        while True:
            with lock:
                if stop_state:
                    return self._join(stop_state[0], live, lock)
            started = 0
            while started < workers:
                with lock:
                    if stop_state:
                        return self._join(stop_state[0], live, lock)
                    st = startable()
                    if st is None:
                        break
                    results[st["id"]] = None       # claimed, not yet finished
                    busy.update(self._step_devices(st))
                # Heartbeat under the write lock: two kernels stamping
                # `now` microseconds apart would let one see a heartbeat dated
                # in its own future, which the dead-man reads as a backwards
                # clock jump and trips on ("no heartbeat for -0 s").
                with write_lock:
                    self._heartbeat(actor)
                t = threading.Thread(target=worker, args=(st,), daemon=True)
                # Registered BEFORE it starts: a fast worker removes itself on
                # its way out, and appending after start() would re-add a
                # thread that has already gone.
                with lock:
                    live.append(t)
                t.start()
                started += 1
            with lock:
                if stop_state:
                    return self._join(stop_state[0], live, lock)
                pending = [st for st in steps
                           if st["id"] not in done and st["id"] not in results
                           or results.get(st["id"]) is None]
                if not pending:
                    return self._finish(cp)
                # Nothing running and nothing startable is a cycle, or an
                # `after` that names a step that will never be done.
                if started == 0 and not live:
                    stuck = [st["id"] for st in steps
                             if st["id"] not in done]
                    return self._fail(
                        cp, stuck[0] if stuck else None,
                        f"no step can start: predecessors or device locks "
                        f"never free for {stuck}")
            self.sleep(0.01)

    def _join(self, cp: dict, live: list, lock):
        """Let the workers already in flight finish, then hand back ``cp``.

        The caller holds ``lock``; it is RELEASED before joining, because a
        worker needs it to record its own outcome. A worker blocked on it would
        be waited for forever, with the run's devices in an unknown state.
        """
        lock.release()
        try:
            for t in list(live):
                t.join()
        finally:
            lock.acquire()
        return cp

    def _finish(self, cp: dict):
        actor = {"person": self.person, "session": cp["run_id"]}
        for mod in self._modules(cp.get("workflow") or {}):
            node = self.manifest.modules[mod].node
            self.clients[node].call("unlock_module", actor, module=mod)
        cp["status"] = "completed"
        cp["ended"] = self.clock()
        _save(cp["run_id"], cp)
        from . import record as _record
        _record.write(cp, self.manifest)
        return cp

    def abandon(self, run_id: str) -> dict:
        """End a failed or waiting run: release its device and module locks
        (a hazardous device goes to its safe state first) and mark it
        ``abandoned``. A person's decision; nothing resumes it."""
        try:
            cp = load_checkpoint(run_id)
        except state_file.StateFileUnreadable as exc:
            # A reader, not a writer: answer, do not decide the file's fate
            # (the damaged file was already moved aside by the quarantine).
            return {"ok": False, "reason": str(exc)}
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

    def _step(self, st: dict, actor: dict, cp: dict, approval=None) -> dict:
        op = st.get("op")
        if op == "set":
            client, d = self._client_for_device(st.get("device", ""))
            if client is None:
                return {"ok": False, "code": "unknown_device",
                        "reason": f"no edge for device {st.get('device')!r}"}
            if approval is None:
                # Sequential runs: the checkpoint names the one step that is
                # waiting. A parallel run passes its own grant in explicitly,
                # because here `cp["pending"]` is shared state two workers
                # write to -- and the grant one of them was resumed with must
                # not be replaced by another step's request id.
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
        if op == "hold":
            # A dwell: the reaction (or the column) needs time, nothing needs
            # to be written. Sleeping here without heartbeating would trip the
            # dead-man on every device this run controls (design §5: "a turn
            # does not heartbeat"), so the dwell heartbeats as it goes. It is
            # also the op a parallel run overlaps: four holds side by side cost
            # their duration once, not four times.
            seconds = float(st.get("seconds", 0))
            if not seconds >= 0:
                return {"ok": False, "code": "bad_step",
                        "reason": f"hold needs seconds >= 0, not {seconds!r}"}
            deadline = self.clock() + seconds
            poll = max(0.01, min(float(st.get("poll_s", 5)), seconds or 1.0))
            beats = 0
            while True:
                left = deadline - self.clock()
                if left <= 0:
                    return {"ok": True, "held_s": seconds, "heartbeats": beats}
                self._heartbeat(actor)
                beats += 1
                self.sleep(min(poll, left))
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
            # A hand-over is a physical move, not a bookkeeping one. The
            # registry above says where the sample is; the destination's
            # position must then HOLD it, or the next module's measure step
            # reads an empty rack and reports a mix-up for a campaign that was
            # carried out correctly. The carrier is the device that declares
            # this position and can hold a container (a rack slot); a module
            # whose positions have no carrier is a registry-only hand-over, and
            # says so rather than claiming a move it did not make.
            carrier, cdev = self._carrier_for(st["to"], st["position"])
            if carrier is None:
                moved = {"registry_only": True,
                         "why": f"{st['to']}:{st['position']} has no device that "
                                f"holds a container"}
            else:
                res = carrier.call("set", actor, device=cdev.id,
                                   channel="container",
                                   value=sample.get("container"))
                if not res.get("ok"):
                    return {"ok": False, "code": res.get("code") or "handover",
                            "reason": f"hand-over to {st['to']}:{st['position']}"
                                      f" did not place {sample.get('container')}"
                                      f": {res.get('reason')}"}
                # Read the slot back: the move is the world's, not the report's.
                back = carrier.call("read", actor, device=cdev.id)
                held = ((back.get("readings") or {}).get("container") or {})
                if held.get("value") != sample.get("container"):
                    return {"ok": False, "code": "handover", "reason":
                            f"hand-over to {st['to']}:{st['position']} read "
                            f"back {held.get('value')!r}, not "
                            f"{sample.get('container')!r}"}
                moved = {"placed_on": cdev.id}
            return {"ok": True, "sample": st["sample"],
                    "at": f"{st['to']}:{st['position']}", **moved}
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
            # An instrument serves ONE position. Without this check a rack of
            # four instruments takes whichever barcode it happens to read and
            # files it under the sample being measured: fine when a module
            # holds one sample, silently wrong the moment two run in parallel.
            served = (d.address or {}).get("position")
            if served is None:
                return {"ok": False, "code": "instrument", "reason":
                        f"{d.id} declares no address.position, so it cannot be "
                        f"shown to serve the position {sid} occupies"}
            if served != sample.get("position"):
                return {"ok": False, "code": "position", "reason":
                        f"instrument mismatch: {d.id} serves position {served}"
                        f", but {sid} is at {sample.get('position')} in "
                        f"{d.module}"}
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
