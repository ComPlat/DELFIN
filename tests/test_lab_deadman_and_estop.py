"""Dead-man and emergency stop hold without the agent loop.

T6 (hung agent), T7 (killed agent, crashed host) and T20 (e-stop while the
agent is busy) from docs/design/lab-control.md section 3. The killed agent
and the e-stop run in REAL separate processes: a SIGKILLed child that held
a heating device, and ``delfin-agent lab estop`` started as its own
program while a thread keeps writing through the tool.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import textwrap
import threading
import time
from pathlib import Path

import pytest

from _lab_support import (FakeClock, approve_and_set, call, lab,  # noqa: F401
                          prepare_heating, sim)

REPO = Path(__file__).resolve().parents[1]


def _kernel(clock=None):
    from delfin.lab.safety import Kernel
    return Kernel.load(clock=clock or time.time)


def _heat(kernel, actor, value=60):
    from delfin.lab import approvals
    for dev in ("stirrer", "fume_hood", "hotplate"):
        assert kernel.set(actor, dev, "power", "on")["ok"]
    r = kernel.set(actor, "hotplate", "setpoint", value, "degC")
    approvals.grant(r["approval_request"], now=kernel.clock())
    out = kernel.set(actor, "hotplate", "setpoint", value, "degC",
                     approval=r["approval_request"])
    assert out["ok"], out


def test_missed_heartbeat_goes_safe(lab):
    from delfin.lab import watchdog
    from delfin.lab.safety import Actor
    clk = FakeClock(time.time())
    k = _kernel(clk)
    me = Actor.here("s1")
    _heat(k, me)
    clk.sleep(20)
    k.heartbeat(me)                       # control: a heartbeat in time
    clk.sleep(20)
    assert watchdog.sweep(k) == []
    assert sim("hotplate")["setpoint"] == 60
    clk.sleep(31)                         # the agent hangs: no heartbeat
    assert watchdog.sweep(k) == ["hotplate"]
    assert sim("hotplate") == {"setpoint": 20, "power": "off"}
    out = k.set(me, "hotplate", "power", "on")
    assert out["code"] == "latched"       # stays down until a person clears


def test_killed_agent_mid_action_goes_safe(lab):
    child = textwrap.dedent("""
        import sys, time
        from delfin.lab import approvals
        from delfin.lab.safety import Actor, Kernel
        k = Kernel.load()
        me = Actor.here("child-session")
        for dev in ("stirrer", "fume_hood", "hotplate"):
            assert k.set(me, dev, "power", "on")["ok"]
        r = k.set(me, "hotplate", "setpoint", 60, "degC")
        approvals.grant(r["approval_request"])
        assert k.set(me, "hotplate", "setpoint", 60, "degC",
                     approval=r["approval_request"])["ok"]
        print("heating", flush=True)
        time.sleep(600)
    """)
    proc = subprocess.Popen([sys.executable, "-c", child], cwd=str(REPO),
                            stdout=subprocess.PIPE, text=True,
                            env=dict(os.environ))
    try:
        assert proc.stdout.readline().strip() == "heating"
        from delfin.lab import watchdog
        k = _kernel()
        assert watchdog.sweep(k) == []    # control: owner alive, in time
        assert sim("hotplate")["setpoint"] == 60
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(timeout=10)
        assert watchdog.sweep(k) == ["hotplate"]
        assert sim("hotplate") == {"setpoint": 20, "power": "off"}
    finally:
        if proc.poll() is None:
            proc.kill()


def test_device_side_watchdog_holds_without_any_process(lab):
    """A crashed host: nothing runs the sweep. The device's own timer
    (Shelly toggle_after, here simulated) puts it safe by itself."""
    from delfin.lab.drivers.simulator import SimulatorDriver
    from delfin.lab.safety import Actor
    clk = FakeClock(time.time())
    k = _kernel(clk)
    _heat(k, Actor.here("s1"))
    hp = k.manifest.devices["hotplate"]
    probe = SimulatorDriver(clock=clk)
    clk.sleep(29)
    assert probe.read(hp)["setpoint"].value == 60     # control: armed, in time
    clk.sleep(2)                                       # 31 s, no heartbeat
    assert probe.read(hp)["setpoint"].value == 20
    assert probe.read(hp)["power"].value == "off"


def test_estop_while_the_agent_is_busy(lab):
    prepare_heating(lab)
    stop = threading.Event()
    seen = []

    def busy_agent():
        speed = 100
        while not stop.is_set():
            speed = 200 if speed == 100 else 100
            seen.append(call(lab, "lab_set", {
                "device": "stirrer", "channel": "speed", "value": speed,
                "unit": "rpm"}))

    t = threading.Thread(target=busy_agent)
    t.start()
    try:
        while not any(r.get("ok") for r in seen):
            time.sleep(0.05)
        res = subprocess.run(
            [sys.executable, "-m", "delfin.agent.cli", "lab", "estop",
             "--reason", "test"], cwd=str(REPO), capture_output=True,
            text=True, timeout=120, env=dict(os.environ))
        assert res.returncode == 0, res.stderr + res.stdout
        n = len(seen)
        deadline = time.time() + 10
        while len(seen) < n + 2 and time.time() < deadline:
            time.sleep(0.05)
    finally:
        stop.set()
        t.join(timeout=30)
    assert sim("stirrer") == {"speed": 0.0, "power": "off"}
    assert sim("hotplate")["power"] == "off"
    assert seen[-1].get("code") == "estop"


def test_writes_blocked_after_estop_until_cleared(lab):
    from delfin.lab import estop
    out = call(lab, "lab_estop", {"reason": "spill"})
    assert out["ok"] and "lamp" in out["stopped"]
    out = call(lab, "lab_set", {"device": "lamp", "channel": "power",
                                "value": "on"})
    assert out["code"] == "estop"
    estop.clear()
    assert sim("lamp")["power"] == "off"          # clearing restarts nothing
    assert call(lab, "lab_set", {"device": "lamp", "channel": "power",
                                 "value": "on"})["ok"]


def test_estop_clear_refused_inside_an_agent_session(lab, monkeypatch):
    from delfin.lab import estop
    call(lab, "lab_estop", {})
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "12345")
    with pytest.raises(estop.EstopError):
        estop.clear()
    res = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "lab",
                          "estop", "--clear"], cwd=str(REPO),
                         capture_output=True, text=True, timeout=120,
                         env=dict(os.environ))
    assert res.returncode != 0 and "person" in res.stderr
    monkeypatch.delenv("DELFIN_AGENT_SESSION")
    assert estop.clear()["cleared"] is True         # control


def test_module_estop_leaves_other_module_running(lab):
    assert call(lab, "lab_set", {"device": "lamp", "channel": "power",
                                 "value": "on"})["ok"]
    out = call(lab, "lab_estop", {"module": "synthesis"})
    assert out["scope"] == "synthesis"
    assert sim("lamp")["power"] == "on"              # analysis untouched
    assert call(lab, "lab_set", {"device": "lamp", "channel": "power",
                                 "value": "off"})["ok"]
    out = call(lab, "lab_set", {"device": "stirrer", "channel": "power",
                                "value": "on"})
    assert out["code"] == "estop"


def test_the_estop_daemon_reasserts_the_safe_state(lab):
    from delfin.lab import estop, store, watchdog
    from delfin.lab.drivers import simulator
    call(lab, "lab_estop", {})
    # something outside DELFIN switches the lamp back on
    with store.locked_json("sim_world.json") as w:
        w["devices"]["lamp"]["values"]["power"] = "on"
    assert watchdog.enforce_estop(_kernel()) == ["lamp"]
    assert simulator.raw("lamp")["values"]["power"] == "off"
    estop.clear()
    with store.locked_json("sim_world.json") as w:
        w["devices"]["lamp"]["values"]["power"] = "on"
    assert watchdog.enforce_estop(_kernel()) == []   # control: no stop, no action
