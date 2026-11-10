"""The lab safety/reporting functions must not say "safe" when the drive failed.

T21 (dead-man), T8 (e-stop) and the coordinator's release from
docs/design/lab-control.md: a device is counted as having been put in its
safe state only when the driver's safe-state write is confirmed by readback.
A driver that fails (simulated "stuck" fault) must not be reported as safe,
and a hazardous device whose safe-state write failed must keep its dead-man
guard rather than being silently abandoned.
"""

from __future__ import annotations

import time
from contextlib import contextmanager

from _lab_support import FakeClock, lab_fixture, sim  # noqa: F401

from delfin.lab import approvals


def _kernel(clock=None):
    from delfin.lab.safety import Kernel
    return Kernel.load(clock=clock or time.time)


def _heat(kernel, actor):
    """Heating above 40 degC: stirrer/hood/hotplate on, setpoint 60 (hazardous,
    approved) -- the actor holds an active hazardous obligation."""
    for dev in ("stirrer", "fume_hood", "hotplate"):
        assert kernel.set(actor, dev, "power", "on")["ok"]
    out = kernel.set(actor, "hotplate", "setpoint", 60, "degC")
    assert out["code"] == "approval_required"
    approvals.grant(out["approval_request"], now=kernel.clock())
    out = kernel.set(actor, "hotplate", "setpoint", 60, "degC",
                     approval=out["approval_request"])
    assert out["ok"], out


@contextmanager
def _active(device):
    from delfin.lab import store
    with store.state() as doc:
        act = doc["active"].get(device, {})
        yield act


@contextmanager
def _world_entry(device):
    from delfin.lab.drivers import simulator
    from delfin.lab import store
    with store.locked_json("sim_world.json") as w:
        yield w["devices"][device]


@contextmanager
def injected(device, fault):
    from delfin.lab.drivers import simulator
    simulator.inject(device, fault)
    try:
        yield
    finally:
        simulator.inject(device, None)


def test_deadman_does_not_say_safe_when_the_drive_failed(lab):
    """The dead-man must not report a device "put in its safe state" whose
    safe-state write the driver failed to perform; it stays latched (guarded)
    but is not claimed safe."""
    from delfin.lab import watchdog
    from delfin.lab.safety import Actor
    clk = FakeClock(time.time())
    k = _kernel(clk)
    actor = Actor.here("s1")
    _heat(k, actor)
    with _active("hotplate") as act:
        act["heartbeat"] = clk() - 10000        # owner stopped heartbeating
    with injected("hotplate", "stuck"):
        tripped = watchdog.sweep(k, now=clk())
        assert "hotplate" not in tripped        # NOT reported as put-safe
        assert sim("hotplate")["power"] == "on"  # in fact still hazardous
        # the latch still guards it: a set is refused until a person clears it
        out = k.set(actor, "hotplate", "setpoint", 20, "degC")
        assert out["code"] == "latched"
    # control (driver healthy): the device genuinely goes safe and is reported
    _clear_latch("hotplate")
    with _active("hotplate") as act:
        act["heartbeat"] = clk() - 10000
    assert watchdog.sweep(k, now=clk()) == ["hotplate"]
    assert sim("hotplate")["power"] == "off"


def test_estop_reassert_does_not_say_fixed_when_the_drive_failed(lab):
    """The e-stop daemon must not report a device "re-asserted" whose safe-state
    write the driver failed to perform (it stays on under the stop)."""
    from delfin.lab import estop, watchdog
    from delfin.lab.safety import Actor
    k = _kernel()
    assert k.set(Actor.here("s1"), "lamp", "power", "on")["ok"]
    estop.trigger(reason="probe")
    # someone turns the lamp back on after the stop, and the driver now faults
    with _world_entry("lamp") as e:
        e["values"]["power"] = "on"
    with injected("lamp", "stuck"):
        fixed = watchdog.enforce_estop(k)
        assert "lamp" not in fixed              # NOT reported as re-asserted
        assert sim("lamp")["power"] == "on"     # still on under the stop
    # control (driver healthy): it is re-asserted and reported fixed
    with _world_entry("lamp") as e:
        e["values"]["power"] = "on"
    fixed = watchdog.enforce_estop(k)
    assert fixed == ["lamp"]
    assert sim("lamp")["power"] == "off"


def test_release_keeps_the_guard_when_the_hazardous_drive_failed(lab):
    """Giving up a hazardous device whose safe-state write the driver failed
    must NOT drop its lock+active guard: the device stays hazardous and
    would otherwise be unattended. It is not reported as released."""
    from delfin.lab import store
    from delfin.lab.safety import Actor, Kernel
    clk = FakeClock(time.time())
    k = Kernel.load(clock=clk)
    actor = Actor.here("s1")
    _heat(k, actor)
    with injected("hotplate", "stuck"):
        released = k.release(actor)
        assert "hotplate" not in released       # NOT reported as released
        snap = store.snapshot()
        assert "hotplate" in snap["locks"]      # the guard stays
        assert "hotplate" in snap["active"]
        assert sim("hotplate")["power"] == "on"
    # control (driver healthy): the device goes safe and is released
    released = k.release(actor)
    assert "hotplate" in released
    assert "hotplate" not in store.snapshot()["locks"]
    assert sim("hotplate")["power"] == "off"


def test_release_drops_a_benign_device(lab):
    """Control: a benign held device is released; its lock is dropped and it
    is not driven safe (nothing to undo)."""
    from delfin.lab import store
    from delfin.lab.safety import Actor, Kernel
    k = Kernel.load()
    actor = Actor.here("s1")
    assert k.set(actor, "stirrer", "power", "on")["ok"]
    assert store.snapshot()["locks"].get("stirrer")
    released = k.release(actor)
    assert "stirrer" in released
    assert "stirrer" not in store.snapshot()["locks"]


def _clear_latch(device):
    from delfin.lab import store
    with store.state() as doc:
        doc["latched"].pop(device, None)
