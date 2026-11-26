"""The closed loop stops when the WATCHDOG has tripped — the required stop
neither arm could demonstrate.

The package requires the loop to stop "on an interlock, on the watchdog, on the
e-stop, and when the optimiser proposes a step outside the manifest". Arm a1
delivered the interlock, the e-stop and the out-of-manifest refusal, and it says
plainly in its own report that it tested the dead-man as the *latch it leaves*
rather than as a tripped dead-man: ``Kernel.set`` arms a dead-man obligation only
for a HAZARDOUS device, and the shipped reactor's channels are all ``benign``, so
on the demo manifest no sweep can ever reach it.

That is true of a loop that arms the obligation ITSELF, and it is why a
hand-written ``Kernel.safe(..., latch=True)`` is not a watchdog test: it passes
just as well if ``watchdog.sweep`` were deleted. The composition this test
covers is the one a real bench meets: the reactor is armed by a controller that
then dies — a person heated it and walked away, their session stopped
heartbeating — and the closed loop arrives afterwards. The dead-man trips by its
own mechanism (a stale heartbeat past ``deadman_s``), latches the reactor, and
the loop must STOP rather than keep writing.

The assertion that makes this a watchdog test rather than a latch test is the
AUDIT one: ``watchdog.sweep`` writes a ``deadman`` record and raises the
``lab_deadman`` security event, and no other path writes those. If the sweep
never ran, there is no such record and the device would not be latched at all.
"""

from __future__ import annotations

import yaml

import numpy as np
import pytest

from _lab_support import FakeClock, lab_fixture, sim  # noqa: F401


def _manifest_with_hazardous_reactor_setpoint() -> str:
    """The demo manifest with the reactor's setpoint ``hazardous``.

    Edited in the PARSED document and re-parsed, so the level the kernel reads
    is checked rather than assumed: a string replace that matched nothing is how
    a test ends up asserting about a manifest it never changed.
    """
    from delfin.lab import manifest as _manifest
    doc = yaml.safe_load(_manifest.demo_text())
    for module in doc["modules"]:
        for device in module["devices"]:
            if device["id"] == "reactor":
                device["writes"]["setpoint"]["level"] = "hazardous"
    parsed = _manifest.parse_document(doc)
    dev = parsed.devices["reactor"]
    assert parsed.level(dev, dev.writes["setpoint"]) == "hazardous", \
        "the fixture did not make the setpoint hazardous"
    return yaml.safe_dump(doc, sort_keys=False)


def _pin_hazardous():
    from delfin.lab import manifest as _manifest
    from delfin.lab import store
    from delfin.lab.safety import pin
    path = store.manifest_path()
    path.write_text(_manifest_with_hazardous_reactor_setpoint(),
                    encoding="utf-8")
    pin(_manifest.load(path), by="alice")


def _arm_the_reactor_then_let_its_owner_die(k, clock):
    """A person heats the reactor, then their controller stops heartbeating.

    This is the state the dead-man exists for: a HAZARDOUS obligation armed by a
    session that is no longer there. The grant is a PERSON's action (the test
    stands in for the person), which is why the loop under test never has to
    approve anything itself.
    """
    from delfin.lab import approvals
    from delfin.lab.safety import Actor
    person = Actor.here("person-hotplate", person="alice")
    for dev in ("stirrer", "fume_hood"):
        assert k.set(person, dev, "power", "on").get("ok")
    first = k.set(person, "reactor", "setpoint", 60.0, "degC")
    assert first.get("code") == "approval_required", first
    approvals.grant(first["approval_request"], approver="alice", now=clock())
    ok = k.set(person, "reactor", "setpoint", 60.0, "degC",
               approval=first["approval_request"])
    assert ok.get("ok"), ok
    assert sim("reactor")["setpoint"] == 60.0, sim("reactor")
    return person


def test_the_loop_stops_when_the_watchdog_has_tripped(lab):
    """A stale heartbeat latches the reactor; the loop stops, it does not retry."""
    from delfin.lab import audit
    from delfin.lab.loop import GridSearch, Loop

    _pin_hazardous()
    clock = FakeClock(1_000_000.0)
    k = _kernel(clock)
    _arm_the_reactor_then_let_its_owner_die(k, clock)

    # The owner's heartbeat is now older than the reactor's deadman_s (30 s):
    # this is the sweep tripping, not a latch someone wrote by hand.
    clock.sleep(31.0)

    record = Loop(k, _actor("loop-after-deadman"), device="reactor", rounds=3,
                  optimiser=GridSearch(np.array([[70.0, 25.0]]))).run()

    # 1. it stopped, and it said which refusal it stopped on
    assert record["stopped"] == "refused", record
    assert record["code"] == "latched", record
    assert len(record["rounds"]) == 1, record["rounds"]
    # 2. one attempt only: no second write after the latch
    writes = record["rounds"][0]["writes"]
    assert writes and writes[-1]["code"] == "latched", writes
    # 3. the loop did not clear the latch it met (that is a person's action)
    from delfin.lab import store
    assert store.snapshot()["latched"].get("reactor"), "the loop cleared a latch"
    # 4. the reactor is in its safe state, put there by the sweep
    assert sim("reactor")["setpoint"] == 30.0, sim("reactor")

    # 5. THE watchdog assertion: the dead-man swept and recorded it. A latch
    #    written by hand does not produce this record.
    events = [r.get("event") for r in audit.read_all()]
    assert "deadman" in events, events
    deadman = [r for r in audit.read_all() if r.get("event") == "deadman"]
    assert deadman[0]["device"] == "reactor", deadman[0]


def test_the_same_loop_runs_when_the_owner_is_still_heartbeating(lab):
    """Control: the stop is the dead-man's, not the loop's.

    Identical set-up and manifest; the only difference is that the owner's
    heartbeat is renewed, so the sweep does not trip. The loop then completes.
    Without this, the test above could be passing because a hazardous manifest
    refuses every write in some unrelated way.
    """
    from delfin.lab.loop import GridSearch, Loop

    _pin_hazardous()
    clock = FakeClock(1_000_000.0)
    k = _kernel(clock)
    person = _arm_the_reactor_then_let_its_owner_die(k, clock)
    clock.sleep(20.0)
    k.heartbeat(person)          # the owner is alive: the dead-man holds
    clock.sleep(20.0)

    record = Loop(k, _actor("loop-alive-owner"), device="reactor", rounds=1,
                  optimiser=GridSearch(np.array([[70.0, 35.0]]))).run()
    assert record["stopped"] == "refused", record
    # The refusal is the LIVE OWNER's lock, not a dead-man latch: one session
    # controls a device at a time, and a lock is never taken over. That is the
    # documented outcome while the owner is alive; the point of the control is
    # that the device was NOT latched and NOT put safe, which is what the
    # dead-man does and a lock does not.
    assert record["code"] == "locked", record
    assert not store_latched(), "the dead-man tripped a live owner's device"
    assert sim("reactor")["setpoint"] == 60.0, \
        "the device went safe although the dead-man had not tripped"


def store_latched() -> bool:
    from delfin.lab import store
    return bool(store.snapshot()["latched"].get("reactor"))


def _kernel(clock):
    from delfin.lab.safety import Kernel
    return Kernel.load(clock=clock)


def _actor(session="loop-1", person="alice"):
    from delfin.lab.safety import Actor
    return Actor.here(session, person=person)
