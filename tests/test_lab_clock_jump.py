"""Dead-man must fail safe when the wall clock jumps backward.

A backward jump (NTP correction, manual clock change) makes
``now - heartbeat`` negative, so the dead-man never fires and a hung
owner's hazardous device stays on. A heartbeat timestamp in the future
relative to the sweep clock is impossible without a jump; the device is
put safe instead. Same for the device-side watchdog in the simulator.
"""
import os

from _lab_support import FakeClock, lab_fixture, sim  # noqa: F401  (registers "lab")


def _actor():
    from delfin.lab import safety, store
    return safety.Actor(person=store.person(), session="s1",
                        pid=os.getpid(), host=store.host())


def _hold(k, actor, lab, dev="hotplate", value=60, clock=None):
    r = k.set(actor, dev, "setpoint", value, unit="degC")
    if not r.get("ok") and r.get("approval_request"):
        from delfin.lab import approvals
        req = r["approval_request"]
        approvals.grant(req, approver=actor.person,
                        now=clock.t if clock else None)
        r = k.set(actor, dev, "setpoint", value, unit="degC", approval=req)
    assert r.get("ok"), r
    return r


def test_backward_clock_jump_still_trips_the_deadman(lab):
    from delfin.lab import safety, store, watchdog
    from _lab_support import prepare_heating

    clock = FakeClock(1000.0)
    k = safety.Kernel.load(clock=clock)
    actor = _actor()
    prepare_heating(lab)
    _hold(k, actor, lab, clock=clock)
    with store.state() as doc:
        hb = float(doc["active"]["hotplate"]["heartbeat"])
    # limit = deadman_s: a heartbeat 300 s in the future is far past it

    # the owner hangs (stops heartbeating) and the clock jumps BACKWARD
    # 300 s: the heartbeat now sits in the future.
    k2 = safety.Kernel.load(clock=lambda: hb - 300.0)
    tripped = watchdog.sweep(k2)
    assert tripped == ["hotplate"], \
        f"dead-man must fire on a future heartbeat, tripped={tripped}"
    assert sim("hotplate").get("power") == "off"

    # control: a fresh hold with the clock consistent is untouched
    _hold(k, actor, lab, clock=clock)
    assert watchdog.sweep(k) == []


def test_backward_clock_jump_fires_the_device_watchdog(lab):
    from delfin.lab import safety
    from delfin.lab.drivers import simulator

    man = safety.Kernel.load().manifest
    clock = FakeClock(1000.0)
    k = safety.Kernel.load(clock=clock)
    d = man.devices["hotplate"]
    sim_d = k.driver(d)
    # hotplate heating, device watchdog armed at clock + deadman_s
    from _lab_support import prepare_heating
    prepare_heating(lab)
    _hold(k, _actor(), lab, clock=clock)
    assert sim_d.arm_watchdog(d, d.deadman_s)
    with simulator._world() as w:
        deadline = w["devices"]["hotplate"]["watchdog"]
    assert deadline > clock.t, "control: deadline is in the future"

    # the clock jumps BACKWARD past the deadline: the device must fail
    # safe rather than wait for a wall time it has already passed.
    sim_d.clock = lambda: deadline - 10.0
    r = sim_d.read(d)
    vals = r.get("values") or sim("hotplate")
    assert vals.get("power") == "off", \
        f"device watchdog must fire on a backward jump, values={vals}"
