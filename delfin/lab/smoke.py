"""Deploy smoke check: the simulator lab, one short plan, end to end.

Runs in a scratch directory (``HOME`` and ``DELFIN_LAB_DIR`` both point
there for this process), so it never touches the installation's manifest,
state or trust store. Steps: stirrer and hood on (benign); hotplate
setpoint 60 degC refused without approval, then approved in-process and
written with readback; the hood switched off is refused by the interlock;
e-stop puts everything safe and refuses the next write; clear lifts the
stop and leaves the devices safe. Prints ``smoke lab ok`` or ``smoke lab
FAILED <step>: <detail>`` and exits 0 / 1.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile


def run(lab: str) -> tuple[bool, str]:
    from . import approvals, estop, manifest, store
    from .drivers import simulator
    from .safety import Actor, Kernel, pin
    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(manifest.demo_text(), encoding="utf-8")
    pin(manifest.load(path), by="smoke")
    k = Kernel.load()
    me = Actor.here("smoke", person="smoke")

    def need(cond, step, detail=""):
        if not cond:
            raise AssertionError(f"{step}: {detail}")

    for dev, ch in (("stirrer", "power"), ("fume_hood", "power"),
                    ("hotplate", "power")):
        r = k.set(me, dev, ch, "on")
        need(r.get("ok"), f"set {dev}.{ch}", r.get("reason", ""))
    r = k.set(me, "hotplate", "setpoint", 60, "degC")
    need(r.get("code") == "approval_required", "approval asked", str(r))
    approvals.grant(r["approval_request"], approver="smoke")
    r = k.set(me, "hotplate", "setpoint", 60, "degC",
              approval=r["approval_request"])
    need(r.get("ok") and abs(r["readback"] - 60) < 0.02, "readback", str(r))
    r = k.set(me, "stirrer", "power", "off")
    need(r.get("code") == "interlock", "interlock", str(r))
    out = estop.trigger(reason="smoke", by="smoke", kernel=k)
    need(out["ok"], "estop", str(out))
    need(simulator.raw("hotplate")["values"]["setpoint"] == 20,
         "estop safe state", str(simulator.raw("hotplate")))
    r = k.set(me, "lamp", "power", "on")
    need(r.get("code") == "estop", "write after estop", str(r))
    estop.clear(by="smoke")
    need(simulator.raw("hotplate")["values"]["power"] == "off",
         "clear restarts nothing")
    r = k.set(me, "lamp", "power", "on")
    need(r.get("ok"), "write after clear", str(r))
    return True, ""


def main(argv=None) -> int:
    scratch = tempfile.mkdtemp(prefix="delfin-lab-smoke-")
    saved = {k: os.environ.get(k) for k in ("HOME", "DELFIN_LAB_DIR",
                                            "DELFIN_LAB_MANIFEST")}
    os.environ["HOME"] = scratch
    os.environ["DELFIN_LAB_DIR"] = os.path.join(scratch, "lab")
    os.environ.pop("DELFIN_LAB_MANIFEST", None)
    try:
        run(scratch)
        print("smoke lab ok")
        return 0
    except Exception as exc:
        print(f"smoke lab FAILED {exc}")
        return 1
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
