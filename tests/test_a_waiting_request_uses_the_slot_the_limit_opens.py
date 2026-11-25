"""A request already waiting uses the slots the shared limit opens for it.

Two properties of one wait, both measured on the untouched base commit
``41fbf18e2`` before any change, with real child processes and one shared
gate directory (``limit`` slots, one free slot opening while a request
waits):

1. **The limit moves while a request waits.** A peer's congestion halving
   lifts, or the user raises ``agent.endpoint_max_concurrency``, and the
   shared limit grows. On the base a waiting request had already sized its
   slot pool at entry, so it could not see the slot the new limit created:
   it sat out the rest of its 900 s deadline beside a slot nobody held.
   Measured: ``SlotWaitTimeout`` after the deadline on the base, served in
   well under a second on the fix.
2. **A slot that opens while a request waits is used.** This is the
   property the defect broke and the one a "fix" can break in the other
   direction: a queue that only serves strictly in arrival order leaves a
   request waiting behind a waiter that is alive but not asking (a paused
   terminal, or a process the OS descheduled), and the free slot goes
   unused. The base passes this one -- its poll is unordered but it does
   take a free slot -- so a change that turns it red has traded a lottery
   for a queue that can wedge the endpoint, and that is worse than the
   defect it set out to fix. It is asserted here so the trade cannot be
   made silently.

The children are real processes, not threads: the gate is a set of flocked
files shared by every DELFIN process of the user, and threads of one process
share a clock and a lock table, so a thread test would prove nothing about
it. The module path is printed by the child, so a run of the wrong tree
cannot pass as evidence.

Run: pytest -q tests/test_a_waiting_request_uses_the_slot_the_limit_opens.py
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

from conftest import child_env
from delfin.agent import endpoint_gate as G


URL = "https://waitslot.example.org/api/v1"
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The wait must be answered by the slot opening, not by the deadline. The
#: number is the gate's, not the host's: the slot opens at ~0.4 s and the
#: child's deadline is 20 s, so a served request is two orders of magnitude
#: inside it and a starved one cannot reach it by being slow.
_SERVED_WITHIN_S = 8.0

#: How long the waiting child is allowed to wait for the slot the raised
#: limit opens. Well above the 0.5 s refresh, far below the product's 900 s,
#: so a starved request reports a timeout instead of hanging the suite.
_WAITER_DEADLINE_S = 12.0


_CHILD = r"""
import json, os, sys, time
sys.path.insert(0, {repo!r})
os.makedirs({runtime!r}, exist_ok=True)
os.environ["DELFIN_RUNTIME_DIR"] = {runtime!r}
os.environ["DELFIN_ENDPOINT_MAX_CONCURRENCY"] = {max_conc!r}
os.environ.pop("DELFIN_ENDPOINT_PER_MODEL_LIMIT", None)
import delfin
from delfin.agent.endpoint_gate import EndpointGate

holder = {holder!r}
if holder:
    # Holds the one slot so the request below has to wait for it.
    h = EndpointGate({url!r}, "m")
    h.acquire(timeout_s=10.0)
    # The directory the gate actually used: the parent must write the shared
    # limit into THIS one, not into a path it reconstructed (its own
    # resolution of the runtime dir can differ from a child's).
    print(json.dumps({{"role": "holder", "module": delfin.__file__,
                       "gate_dir": str(h._dir())}}), flush=True)
    time.sleep({hold!r})
    h.release()
    sys.exit(0)

time.sleep({delay!r})
g = EndpointGate({url!r}, "m")
t0 = time.monotonic()
try:
    idx = g.acquire(timeout_s={deadline!r})
except Exception as exc:                      # half the measurement
    print(json.dumps({{"role": "waiter", "module": delfin.__file__,
                       "error": type(exc).__name__,
                       "waited": round(time.monotonic() - t0, 3)}}),
          flush=True)
    sys.exit(0)
print(json.dumps({{"role": "waiter", "module": delfin.__file__,
                   "slot": idx, "waited": round(time.monotonic() - t0, 3)}}),
      flush=True)
g.release()
"""


def _start(env, *, runtime, holder, hold, delay, deadline, max_conc="1"):
    src = _CHILD.format(repo=_REPO, runtime=str(runtime), url=URL,
                        holder=holder, hold=hold, delay=delay,
                        deadline=deadline, max_conc=max_conc)
    return subprocess.Popen([sys.executable, "-c", src], env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)


def _records(proc):
    out, err = proc.communicate(timeout=120)
    recs = [json.loads(line) for line in out.splitlines()
            if line.startswith("{")]
    assert recs, f"the child reported nothing; stderr was:\n{err}"
    return recs


def _first_record(proc):
    """The first JSON line a still-running child has printed."""
    line = proc.stdout.readline()
    assert line.startswith("{"), f"the child said {line!r} instead"
    return json.loads(line)


@pytest.fixture
def env(tmp_path):
    runtime = tmp_path / "runtime"
    runtime.mkdir()
    environ = child_env(tmp_path)
    environ["DELFIN_RUNTIME_DIR"] = str(runtime)
    # Two slots CONFIGURED, one in force: that is a gate whose peers halved
    # the limit after congestion, which recovers on its own. The configured
    # maximum is what clamps the shared state, so a test that raised the
    # state past it would be measuring a setting no process agreed to.
    environ["DELFIN_ENDPOINT_MAX_CONCURRENCY"] = "2"
    environ["DELFIN_ENDPOINT_PER_MODEL_LIMIT"] = ""
    return environ, runtime


def test_a_slot_the_shared_limit_opens_mid_wait_is_used(env):
    """The limit grows while a request waits: it must take the new slot.

    The holder keeps the ONE slot the whole time, so the only slot the
    waiting request can ever get is the one the raised limit creates. Red on
    the base: the waiting request had sized its pool at entry (the single
    slot the limit allowed then), so the second slot was outside its scan
    and it timed out beside a slot nobody held. Measured on the base: the
    holder kept slot 0, the limit rose to 2 at 0.6 s, and the waiter still
    reported ``SlotWaitTimeout``; with the fix it was served in ~0.5 s.
    """
    environ, runtime = env
    holder = _start(environ, runtime=runtime, holder=True, hold=15.0,
                    delay=0.0, deadline=30.0)
    try:
        # Read the holder's report while it is still holding: it names the
        # gate directory the children share, so the limit below is written
        # where the waiter is actually looking.
        time.sleep(0.6)
        gate_dir = _first_record(holder)["gate_dir"]
        gdir = Path(gate_dir)
        assert gdir.is_dir(), f"the holder reported no gate dir: {gate_dir}"
        # The position after congestion: the limit is halved in the shared
        # state, then recovers. Two slots are CONFIGURED, so the state is
        # what is in force -- a single slot until it recovers.
        G._write_state(gdir, {"limit": 1})
        waiter = _start(environ, runtime=runtime, holder=False, hold=0.0,
                        delay=0.0, deadline=_WAITER_DEADLINE_S,
                        max_conc="2")
        time.sleep(0.6)                   # it is waiting now
        # Recover, exactly as a peer finishing its quiet period does.
        G._write_state(gdir, {"limit": 2})
        assert G.EndpointGate(URL, "m").current_limit(gdir) == 2, (
            "the raised limit did not take effect, so this test would "
            "measure nothing")
        recs = _records(waiter)
    finally:
        holder.kill()
        holder.wait(timeout=60)
    rec = recs[-1]
    assert rec.get("error") is None, (
        f"the waiting request did not use the slot the raised limit opened: "
        f"{rec.get('error')} after {rec.get('waited')} s with the holder "
        f"still in slot 0 and a second slot free")
    assert rec["slot"] == 1, (
        f"served on slot {rec['slot']}, but slot 0 was held the whole time "
        f"-- the measurement did not run as intended")
    assert rec["waited"] <= _SERVED_WITHIN_S, (
        f"served only after {rec['waited']} s, which is the deadline "
        f"answering rather than the slot")


def test_a_slot_freed_while_waiting_is_taken(env):
    """The slot a holder releases must go to the request already waiting.

    This is the control case: the base passes it, so it fails only if a
    change stops serving a slot that is free -- which is what an arrival
    order enforced as a hard barrier does to a waiter that is alive but not
    asking. Asserted so that trade cannot be made silently.
    """
    environ, runtime = env
    holder = _start(environ, runtime=runtime, holder=True, hold=1.0,
                    delay=0.0, deadline=20.0)
    time.sleep(0.4)
    waiter = _start(environ, runtime=runtime, holder=False, hold=0.0,
                    delay=0.0, deadline=20.0)
    recs = _records(waiter)
    holder.wait(timeout=60)
    rec = recs[-1]
    assert rec.get("error") is None, (
        f"a free slot went unused while a request waited: {rec.get('error')}")
    assert rec["waited"] <= _SERVED_WITHIN_S, (
        f"served only after {rec['waited']} s, which is the deadline "
        f"answering rather than the freed slot")


def test_a_waiter_that_never_asks_again_does_not_hold_the_slot(
        env, tmp_path):
    """A waiter that is frozen while queued must not deny the slot to others.

    Not a defect test: the base passes it (its poll has no queue). It is the
    guard against the wrong fix. A strictly arrival-ordered queue lets a
    request that is *alive but not asking* -- a paused terminal, a process
    the OS descheduled -- hold the free slot for everyone behind it, which
    turns 268 occasional timeouts into a certain one. Measured on such a
    queue: the two later requests timed out after their whole deadline while
    the slot had been free since the holder released it.
    """
    environ, runtime = env
    holder = _start(environ, runtime=runtime, holder=True, hold=1.5,
                    delay=0.0, deadline=20.0)
    time.sleep(0.4)
    frozen = _start(environ, runtime=runtime, holder=False, hold=0.0,
                    delay=0.0, deadline=20.0)      # queues behind the holder
    time.sleep(0.4)
    frozen.send_signal(signal.SIGSTOP)             # alive, not asking
    try:
        later = _start(environ, runtime=runtime, holder=False, hold=0.0,
                       delay=0.0, deadline=6.0)
        recs = _records(later)
        rec = recs[-1]
        assert rec.get("error") is None, (
            f"a free slot was held by a waiter that had stopped asking, so "
            f"the next request timed out ({rec.get('error')} after "
            f"{rec.get('waited')} s)")
    finally:
        frozen.send_signal(signal.SIGCONT)
        _records(frozen)
        holder.wait(timeout=60)
