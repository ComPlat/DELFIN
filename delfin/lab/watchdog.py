"""Dead-man and e-stop enforcement, run outside the agent loop.

Input: the shared lab state. Output: devices put in their safe state.

``sweep(kernel)`` (the dead-man) puts a device with a hazardous obligation
in its safe state and latches it when its controlling actor

* sent no heartbeat for the device's ``deadman_s``, or
* is a process on this host that no longer exists (killed, crashed).

It also drops plain locks whose heartbeat is older than ``deadman_s``.

``enforce_estop(kernel)`` re-asserts the safe state of every device under an
active e-stop whose readback is not safe, and alerts when it does not hold.

``serve(role)`` runs either loop (``deadman``, ``estop`` or ``both``) in its
own process: the systemd units ``delfin-lab-watchdog`` and
``delfin-lab-estop`` on a Pi, ``delfin-agent lab daemon`` on a plain
workstation. The kernel also calls ``sweep`` before every write, so a dead
owner is put safe even where no daemon runs; the daemon is what covers the
time between writes.
"""

from __future__ import annotations

import signal
import sys
import time
from typing import Optional

from . import audit as _audit
from . import store as _store


def sweep(kernel, *, now: Optional[float] = None) -> list[str]:
    if kernel.manifest is None:
        return []
    now = kernel.clock() if now is None else now
    here = _store.host()
    snap = _store.snapshot()
    tripped = []
    for dev_id, act in list(snap["active"].items()):
        d = kernel.manifest.devices.get(dev_id)
        if d is None or (kernel.node and d.node != kernel.node):
            continue
        age = now - float(act.get("heartbeat", 0))
        dead = (act.get("host") == here and act.get("pid")
                and not _store.pid_alive(act["pid"]))
        if age <= d.deadman_s and not dead:
            continue
        why = (f"dead-man: owner process {act.get('pid')} is gone" if dead
               else f"dead-man: no heartbeat from {act.get('person')} "
                    f"(session {act.get('session')}) for {age:.0f} s "
                    f"(limit {d.deadman_s:g} s)")
        res = kernel.safe(dev_id, why, latch=True)
        _audit.record("deadman", device=dev_id, reason=why,
                      person=act.get("person"), session=act.get("session"))
        if res.get("ok"):
            tripped.append(dev_id)
        else:
            _audit.alert("lab_deadman", f"{dev_id}: dead-man could not put "
                         "the device in its safe state; it is latched and "
                         "needs manual attention")
    with _store.state() as doc:
        for dev_id, lock in list(doc["locks"].items()):
            d = kernel.manifest.devices.get(dev_id)
            limit = d.deadman_s if d else 30.0
            if dev_id in doc["active"]:
                continue
            if now - float(lock.get("heartbeat", 0)) > limit:
                doc["locks"].pop(dev_id, None)
    return tripped


def enforce_estop(kernel) -> list[str]:
    if kernel.manifest is None:
        return []
    snap = _store.snapshot()
    g = snap["estop"].get("global")
    mods = snap["estop"].get("modules") or {}
    fixed = []
    for d in kernel.manifest.devices.values():
        if kernel.node and d.node != kernel.node:
            continue
        if not (g or mods.get(d.module)) or not d.safe_state:
            continue
        if kernel._at_safe(d, kernel._fresh(d)):
            continue
        res = kernel.safe(d.id, "e-stop re-asserted", latch=False)
        if res.get("ok"):
            fixed.append(d.id)
        else:
            _audit.alert("lab_estop", f"{d.id}: safe state does not hold "
                         "under the e-stop")
    return fixed


def serve(role: str = "both", *, interval: float = 1.0,
          once: bool = False, out=sys.stdout) -> int:
    from .safety import Kernel
    stop = {"flag": False}

    def _term(*_a):
        stop["flag"] = True

    signal.signal(signal.SIGTERM, _term)
    noted = False
    while not stop["flag"]:
        kernel = Kernel.load()
        if kernel.manifest is None:
            if not noted:
                print(f"lab daemon ({role}): {kernel.load_error}; waiting",
                      file=out, flush=True)
                noted = True
        else:
            noted = False
            if role in ("deadman", "both"):
                for dev in sweep(kernel):
                    print(f"dead-man: {dev} put in its safe state", file=out,
                          flush=True)
            if role in ("estop", "both"):
                for dev in enforce_estop(kernel):
                    print(f"e-stop: {dev} safe state re-asserted", file=out,
                          flush=True)
        if once:
            break
        time.sleep(interval)
    return 0
