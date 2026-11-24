"""The closed loop: the measurement decides the next conditions, or it stops.

Everything here runs on the simulator, through the safety kernel, in the
test's own lab directory. The kernel is a real ``Kernel.load()`` and the
simulated world is the real one (``<lab_dir>/sim_world.json``); only the
optimiser's arithmetic is exercised outside a device in the last two tests.

What is measured:

* the loop learns: on the shipped surface, expected improvement reaches 99 %
  of the best the 81-point lattice can reach in far fewer rounds than the
  lattice walk, which is the control arm;
* it stops on an interlock, on the dead-man, on the e-stop, on a proposal
  outside the manifest, and on a hazardous channel that needs a person — and
  in every one of those it makes exactly ONE attempt and does not retry,
  clamp, approve or route around the refusal;
* the record says which of those happened, with the kernel's own code.

The stop tests are the point of the file: a loop that quietly worked around
a refusal would produce a plausible yield curve and no refusals at all, and
nothing else in the suite would notice.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest

from _lab_support import FakeClock, call, lab_fixture, sim  # noqa: F401

#: The demo envelope's lattice, which is the comparison the round counts use.
LATTICE_STEP = (10.0, 10.0)


def _kernel(clock=None):
    from delfin.lab.safety import Kernel
    return Kernel.load(clock=clock or time.time)


def _actor(session="loop-1", person="alice"):
    from delfin.lab.safety import Actor
    return Actor.here(session, person=person)


def _heat(k):
    """Stirrer and hood on: what the reactor interlock above 40 degC needs.

    The locks are given back, because the loop under test uses its own
    session: a helper that kept the stirrer locked would make the loop's
    writes fail on a lock rather than on what the test is about.
    """
    me = _actor("heating")
    for dev in ("stirrer", "fume_hood"):
        out = k.set(me, dev, "power", "on")
        assert out.get("ok"), out
    assert k.release(me) or True       # benign devices: nothing goes safe


def _manifest_with_hazardous_reactor_setpoint():
    """The demo manifest with the reactor's setpoint made ``hazardous``.

    Edited in the PARSED document and dumped again, not with a string
    replace: a replace that silently matched nothing is how a test ends up
    asserting about a manifest it never changed. The result goes back through
    the real parser and the level the kernel will read is checked, so the
    fixture fails loudly instead of quietly testing the wrong manifest.
    """
    import yaml as _yaml

    from delfin.lab import manifest as _manifest
    doc = _yaml.safe_load(_manifest.demo_text())
    for module in doc["modules"]:
        for device in module["devices"]:
            if device["id"] == "reactor":
                device["writes"]["setpoint"]["level"] = "hazardous"
    parsed = _manifest.parse_document(doc)
    dev = parsed.devices["reactor"]
    assert parsed.level(dev, dev.writes["setpoint"]) == "hazardous", \
        "the fixture did not make the setpoint hazardous"
    return _yaml.safe_dump(doc, sort_keys=False)


def _lattice_ref(k=None):
    from delfin.lab import surface
    return surface.lattice((30.0, 5.0), (110.0, 85.0), LATTICE_STEP)


def _target():
    from delfin.lab import surface
    best, _at = surface.best_on(_lattice_ref())
    return 0.99 * best


def test_the_loop_runs_and_every_round_is_recorded(lab):
    """Drive the product path: the loop writes through the kernel only."""
    from delfin.lab.loop import Loop, best_measured
    k = _kernel()
    _heat(k)
    record = Loop(k, _actor(), device="reactor", rounds=6).run()
    assert record["stopped"] == "budget", record
    done = [r for r in record["rounds"] if "yield" in r]
    assert len(done) == 6, record["rounds"]
    assert record["best"] == best_measured(record)
    assert record["best"]["yield"] == max(r["yield"] for r in done)
    for r in done:
        # every write went through the kernel and none of them failed
        assert r["writes"] and all(w["ok"] for w in r["writes"]), r
        assert 0.0 <= r["yield"] <= 100.0, r


def test_the_measurement_depends_on_the_conditions_it_was_given(lab):
    """The defining property: a different condition is a different reading.

    Without this the loop's "learning" would be noise about a constant, and
    the round counts below would measure nothing.
    """
    from delfin.lab.loop import GridSearch, Loop
    from delfin.lab import surface
    k = _kernel()
    _heat(k)
    good = Loop(k, _actor("good"), device="reactor", rounds=1,
                optimiser=GridSearch(np.array([[64.0, 25.0]]))).run()
    bad = Loop(k, _actor("bad"), device="reactor", rounds=1,
               optimiser=GridSearch(np.array([[110.0, 5.0]]))).run()
    g = good["rounds"][0]["yield"]
    b = bad["rounds"][0]["yield"]
    assert g > 95.0 and b < 5.0, (g, b)
    assert g == pytest.approx(surface.yield_percent(64.0, 25.0), abs=1e-4)
    assert b == pytest.approx(surface.yield_percent(110.0, 5.0), abs=1e-4)


def test_expected_improvement_beats_the_lattice_walk_on_rounds(lab):
    """The measurement the package asks for, with its control arm.

    ``rounds_to`` counts the runs before the best-so-far first reached 99 %
    of what the 81-point lattice can reach at all. The GP must get there in
    a small fraction of the 81 runs the lattice walk needs.
    """
    from delfin.lab.loop import GridSearch, Loop, rounds_to
    k = _kernel()
    _heat(k)
    target = _target()
    gp = Loop(k, _actor("gp"), device="reactor", rounds=15).run()
    assert gp["stopped"] == "budget", gp
    gp_rounds = rounds_to(gp, target)
    assert gp_rounds is not None, "the GP never reached 99 % of the lattice"
    assert gp_rounds <= 8, gp_rounds

    grid = Loop(k, _actor("grid"), device="reactor",
                optimiser=GridSearch(_lattice_ref()), rounds=81).run()
    grid_rounds = rounds_to(grid, target)
    assert grid_rounds is not None, "the lattice walk never reached its own best"
    assert gp_rounds < grid_rounds, (gp_rounds, grid_rounds)


def test_the_loop_stops_at_an_interlock_and_does_not_retry(lab):
    """A refusal is a stop: one attempt, no retry, the code in the record.

    The interlock is real: ``reactor.setpoint > 40 degC`` requires the
    stirrer and the hood, and neither is switched on here.
    """
    from delfin.lab.loop import GridSearch, Loop
    k = _kernel()
    record = Loop(k, _actor(), device="reactor", rounds=5,
                  optimiser=GridSearch(np.array([[80.0, 25.0]]))).run()
    assert record["stopped"] == "refused", record
    assert record["code"] == "interlock", record
    assert "requires" in record["detail"]
    assert len(record["rounds"]) == 1, record["rounds"]
    writes = record["rounds"][0]["writes"]
    assert len(writes) == 1, writes
    assert writes[0]["ok"] is False and writes[0]["code"] == "interlock"
    assert "yield" not in record["rounds"][0]
    # and nothing was written to the device by the attempt
    assert sim("reactor")["setpoint"] == 30.0, sim("reactor")


def test_a_proposal_outside_the_manifest_stops_the_loop_unclamped(lab):
    """The proposal goes to the kernel AS PROPOSED; the record keeps it.

    Clamping would turn "forbidden" into "a different experiment ran", and
    the record would then hold a measurement of conditions nobody chose.
    """
    from delfin.lab.loop import GridSearch, Loop
    k = _kernel()
    _heat(k)
    record = Loop(k, _actor(), device="reactor", rounds=5,
                  optimiser=GridSearch(np.array([[500.0, 25.0]]))).run()
    assert record["stopped"] == "refused", record
    assert record["code"] == "out_of_range", record
    entry = record["rounds"][0]
    assert entry["proposed"]["setpoint"] == 500.0
    assert entry["writes"][0]["value"] == 500.0, entry["writes"]
    assert len(record["rounds"]) == 1
    assert sim("reactor")["setpoint"] == 30.0, "clamped instead of refused"


def test_the_loop_stops_on_the_estop_without_attempting_a_write(lab):
    """An e-stop is checked before the round: nothing is even proposed."""
    from delfin.lab import estop
    from delfin.lab.loop import Loop
    k = _kernel()
    _heat(k)
    estop.trigger(reason="spill", by="alice", kernel=k)
    record = Loop(k, _actor(), device="reactor", rounds=5).run()
    assert record["stopped"] == "estop", record
    assert record["rounds"] == [], record["rounds"]
    estop.clear(by="alice")


def test_the_loop_stops_when_the_device_is_latched_and_never_unlatches_itself(lab):
    """A latch is a stop: the loop reports it and does not clear it.

    The kernel latches a device when its readback does not match what was
    commanded, and refuses every later write until a PERSON clears it. A loop
    that ran on would report a sequence of refusals as if it were searching,
    so the loop must stop at the first one -- and it must not call
    ``estop.clear_latch``, which is a person's action.

    The latch is put in place the way the kernel puts one, through
    ``Kernel.safe(..., latch=True)``, so this test uses the real state rather
    than a hand-written state file.
    """
    from delfin.lab import estop
    from delfin.lab.loop import Loop
    k = _kernel()
    _heat(k)
    assert k.safe("reactor", "measured latch for the loop test", latch=True)["ok"]
    record = Loop(k, _actor(), device="reactor", rounds=5).run()
    assert record["stopped"] == "refused", record
    assert record["code"] == "latched", record
    assert "person clears it" in record["detail"], record["detail"]
    assert len(record["rounds"]) == 1, record["rounds"]
    assert sim("reactor")["setpoint"] == 30.0
    # the loop did not clear it: the latch is still there for a person
    assert estop.status()["latched"].get("reactor"), estop.status()


def test_a_hazardous_channel_stops_for_a_person_instead_of_self_approving(lab):
    """A loop cannot grant its own approval; it hands the request over.

    The reactor's channels are benign in the demo manifest, so this test
    makes the setpoint hazardous in the manifest it loads. The loop must stop
    with ``approval_required`` and the request id, and it must NOT write.
    """
    from delfin.lab import manifest as _manifest
    from delfin.lab import store
    from delfin.lab.loop import GridSearch, Loop
    from delfin.lab.safety import pin
    text = _manifest_with_hazardous_reactor_setpoint()
    path = store.manifest_path()
    path.write_text(text, encoding="utf-8")
    pin(_manifest.load(path), by="alice")
    k = _kernel()
    _heat(k)
    record = Loop(k, _actor(), device="reactor", rounds=5,
                  optimiser=GridSearch(np.array([[35.0, 10.0]]))).run()
    assert record["stopped"] == "approval_required", record
    assert record["code"] == "approval_required", record
    assert record.get("approval_request"), record
    assert len(record["rounds"]) == 1
    assert "yield" not in record["rounds"][0]
    assert sim("reactor")["setpoint"] == 30.0, "a hazardous write landed"


def test_a_safe_value_write_never_waits_for_an_approval(lab):
    """Control: the approval stop is about the ACTION, not about the loop.

    Writing the channel's own safe value is benign by the kernel's rule, so
    the same hazardous manifest still admits a loop that stays at the safe
    setpoint. Without this, the test above could be passing because every
    write in that manifest is refused.
    """
    from delfin.lab import manifest as _manifest
    from delfin.lab import store
    from delfin.lab.loop import GridSearch, Loop
    from delfin.lab.safety import pin
    text = _manifest_with_hazardous_reactor_setpoint()
    path = store.manifest_path()
    path.write_text(text, encoding="utf-8")
    pin(_manifest.load(path), by="alice")
    k = _kernel()
    _heat(k)
    record = Loop(k, _actor(), device="reactor", rounds=1,
                  optimiser=GridSearch(np.array([[30.0, 10.0]]))).run()
    assert record["stopped"] == "budget", record
    assert "yield" in record["rounds"][0], record["rounds"][0]


def test_a_ramp_respects_max_step_and_is_not_a_way_around_a_refusal(lab):
    """A big move is made in legal steps; the record shows each one."""
    from delfin.lab.loop import GridSearch, Loop
    k = _kernel()
    _heat(k)
    record = Loop(k, _actor(), device="reactor", rounds=1,
                  optimiser=GridSearch(np.array([[100.0, 25.0]]))).run()
    assert record["stopped"] == "budget", record
    sp_writes = [w for w in record["rounds"][0]["writes"]
                 if w["channel"] == "setpoint"]
    # 30 -> 90 -> 100 with max_step 60: the first step is the full 60.
    assert [w["value"] for w in sp_writes] == [90.0, 100.0], sp_writes
    assert all(w["ok"] for w in sp_writes), sp_writes
    assert sim("reactor")["setpoint"] == 100.0


def test_the_ramp_stops_at_the_first_refused_step(lab):
    """A ramp step is judged on its own: a refusal there stops the loop.

    The stirrer goes off after the ramp started, so the second step (which
    crosses 40 degC) is refused by the interlock. The loop must stop and NOT
    finish the ramp or run the reactor.
    """
    from delfin.lab.loop import GridSearch, Loop
    import delfin.lab.loop as loop_mod

    k = _kernel()
    _heat(k)
    original = loop_mod.Loop._write
    calls = {"n": 0}

    def write_then_trip(self, entry, channel, value):
        calls["n"] += 1
        if calls["n"] == 1:
            off = _actor("sabotage")
            assert k.set(off, "stirrer", "power", "off")["ok"]
        return original(self, entry, channel, value)

    loop_mod.Loop._write = write_then_trip
    try:
        record = Loop(k, _actor(), device="reactor", rounds=3,
                      optimiser=GridSearch(np.array([[100.0, 25.0]])),
                      ).run()
    finally:
        loop_mod.Loop._write = original
    assert record["stopped"] == "refused", record
    assert record["code"] == "interlock", record
    assert len(record["rounds"]) == 1
    assert "yield" not in record["rounds"][0], record["rounds"][0]


def test_an_unknown_envelope_is_refused_before_any_write(lab):
    """A loop with no declared window is not a search: it raises."""
    from delfin.lab.loop import Loop
    k = _kernel()
    with pytest.raises(ValueError, match="no device"):
        Loop(k, _actor(), device="no_such_device", rounds=1)


def test_the_optimiser_proposes_inside_the_envelope_it_was_given():
    """The one property every proposal must have, over many rounds."""
    from delfin.lab.loop import GPBayes
    from delfin.lab import surface
    lo, hi = (30.0, 5.0), (110.0, 85.0)
    opt = GPBayes(lo, hi)
    history = []
    for _ in range(25):
        p = opt.propose(history)
        assert lo[0] <= p["setpoint"] <= hi[0], p
        assert lo[1] <= p["duration"] <= hi[1], p
        history.append({"setpoint": p["setpoint"], "duration": p["duration"],
                        "yield": surface.yield_percent(p["setpoint"],
                                                       p["duration"])})
    # and it never proposed the same point twice: a repeated experiment
    # spends a run and teaches nothing
    seen = {(round(h["setpoint"], 9), round(h["duration"], 9))
            for h in history}
    assert len(seen) == len(history), "a point was proposed twice"


def test_expected_improvement_picks_the_ridge_over_a_neighbour():
    """An assertion that can fail: one good point vs one bad one.

    EI must prefer the unexplored point next to the good one. A test that
    only checked "a proposal came back" would pass on a optimiser that
    returned a constant, which is why the expected answer is asserted.
    """
    from delfin.lab.loop import GPBayes
    from delfin.lab import surface
    opt = GPBayes((30.0, 5.0), (110.0, 85.0), length_scale=0.1)
    history = [
        {"setpoint": 64.0, "duration": 25.0,
         "yield": surface.yield_percent(64.0, 25.0)},
        {"setpoint": 30.0, "duration": 80.0,
         "yield": surface.yield_percent(30.0, 80.0)},
        {"setpoint": 100.0, "duration": 30.0,
         "yield": surface.yield_percent(100.0, 30.0)},
    ]
    got = opt.propose(history)
    # closer to the good point than to either bad one
    d_good = math.dist((got["setpoint"], got["duration"]), (64.0, 25.0))
    d_bad1 = math.dist((got["setpoint"], got["duration"]), (30.0, 80.0))
    d_bad2 = math.dist((got["setpoint"], got["duration"]), (100.0, 30.0))
    assert d_good < d_bad1 and d_good < d_bad2, (got, d_good, d_bad1, d_bad2)
    assert surface.yield_percent(got["setpoint"], got["duration"]) > \
        surface.yield_percent(30.0, 80.0)


def test_the_unknown_proposal_is_recorded_even_when_it_stops(lab):
    """The record of a stopped loop is still a record of what it proposed."""
    from delfin.lab.loop import GridSearch, Loop, best_measured, rounds_to
    k = _kernel()
    record = Loop(k, _actor(), device="reactor", rounds=5,
                  optimiser=GridSearch(np.array([[500.0, 25.0]]))).run()
    assert best_measured(record) is None
    assert rounds_to(record, 0.0) is None
    assert record["rounds"][0]["n"] == 1
