"""Coordinator and edges: envelopes, partitions, modules, sample hand-over.

T8, T13 (envelope replay), T16-T19 of docs/design/lab-control.md, and the
end-to-end test of phase 2: a simulated synthesis module heats a sample
(with a person's approval), hands it to the analysis module, which checks
the container's barcode and measures it; the run record and the Chemotion
export are written. Two edges share one lab directory here (one host);
each edge's kernel accepts only its own node's devices.
"""

from __future__ import annotations

import copy
import json
import time

import pytest

from _lab_support import FakeClock, lab  # noqa: F401


def _setup(clock):
    from delfin.lab import edge
    from delfin.lab.safety import Kernel
    k = Kernel.load()
    edges, clients = {}, {}
    for node in ("edge_synth", "edge_analysis"):
        key = edge.provision_key(node)
        edges[node] = edge.EdgeNode(node, key=key, clock=clock)
        clients[node] = edge.EdgeClient(node, key, transport=edges[node].handle,
                                        manifest_checksum=k.manifest.checksum,
                                        coordinator="c1", clock=clock)
    return k.manifest, edges, clients


WF = {
    "id": "wf-demo",
    "samples": {"S-001": {"container": "vial-7", "module": "synthesis",
                          "position": "S1"}},
    "steps": [
        {"id": "stir", "op": "set", "device": "stirrer", "channel": "power",
         "value": "on"},
        {"id": "hood", "op": "set", "device": "fume_hood", "channel": "power",
         "value": "on"},
        {"id": "hp_on", "op": "set", "device": "hotplate", "channel": "power",
         "value": "on"},
        {"id": "heat", "op": "set", "device": "hotplate", "channel":
         "setpoint", "value": 60, "unit": "degC"},
        {"id": "hot", "op": "wait", "device": "hotplate_probe", "channel":
         "temperature", "cmp": ">=", "value": 55, "unit": "degC",
         "timeout_s": 900, "poll_s": 5},
        {"id": "cool", "op": "set", "device": "hotplate", "channel":
         "setpoint", "value": 20, "unit": "degC", "after": ["hot"]},
        {"id": "move", "op": "transfer", "sample": "S-001", "from":
         "synthesis", "to": "analysis", "position": "A1"},
        {"id": "uv", "op": "measure", "device": "uvvis", "sample": "S-001",
         "reader": "reader"},
    ],
}


def _coordinator(manifest, clients, clk, person="alice", name="c1"):
    from delfin.lab.workflow import Coordinator
    return Coordinator(manifest, clients, person=person, name=name,
                       clock=clk, sleep=clk.sleep)


def _approve_pending(cp, clk):
    from delfin.lab import approvals
    approvals.grant(cp["pending"]["approval"], now=clk())


def _place_vial(clients, barcode="vial-7"):
    from delfin.lab.drivers import simulator
    clients["edge_analysis"].call("read", {"person": "alice",
                                           "session": "probe"},
                                  device="reader")
    clients["edge_analysis"].call("read", {"person": "alice",
                                           "session": "probe"},
                                  device="uvvis")
    simulator.place("reader", barcode=barcode)
    simulator.place("uvvis", signal_next=0.42)


def test_two_module_workflow_with_sample_handover(lab):
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    _place_vial(clients)
    coord = _coordinator(manifest, clients, clk)
    cp = coord.run(WF, "run-1")
    assert cp["status"] == "waiting_approval"
    assert cp["done"] == ["stir", "hood", "hp_on"]
    _approve_pending(cp, clk)                         # the person's step
    cp = coord.run(WF, "run-1")
    assert cp["status"] == "completed", cp["deviations"]
    assert cp["samples"]["S-001"]["module"] == "analysis"
    assert cp["measurements"][0]["value"] == 0.42
    assert cp["measurements"][0]["container"] == "vial-7"
    from delfin.lab.workflow import run_dir
    rec = json.loads((run_dir("run-1") / "record.json").read_text())
    assert rec["schema"] == "delfin.lab.run/1"
    assert rec["manifest"]["checksum"] == manifest.checksum
    assert rec["fair"]["identifier"].endswith(":run-1")
    assert any(r.get("event") == "intent" for r in rec["device_log"])
    chem = json.loads((run_dir("run-1") / "chemotion_export.json").read_text())
    assert chem["samples"][0]["analyses"][0]["value"] == 0.42
    assert chem["unverified"]


def test_sample_mixup_stops_the_run(lab):
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    _place_vial(clients, barcode="vial-9")            # the wrong vial arrived
    coord = _coordinator(manifest, clients, clk)
    cp = coord.run(WF, "run-mix")
    _approve_pending(cp, clk)
    cp = coord.run(WF, "run-mix")
    assert cp["status"] == "failed" and cp["failed_step"] == "uv"
    assert "mix-up" in cp["deviations"][-1]["reason"]
    assert cp["measurements"] == []
    # a transfer from the wrong module is refused before anything moves
    wf = copy.deepcopy(WF)
    wf["steps"] = [{"id": "m", "op": "transfer", "sample": "S-001",
                    "from": "analysis", "to": "synthesis", "position": "S2"}]
    cp = _coordinator(manifest, clients, clk).run(wf, "run-wrong")
    assert cp["status"] == "failed" and "is in synthesis" in \
        cp["deviations"][-1]["reason"]


def test_resume_after_failure_from_last_safe_step(lab):
    from delfin.lab.drivers import simulator
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    _place_vial(clients)
    simulator.inject("uvvis", "offline")
    coord = _coordinator(manifest, clients, clk)
    cp = coord.run(WF, "run-r")
    _approve_pending(cp, clk)
    cp = coord.run(WF, "run-r")
    assert cp["status"] == "failed" and cp["failed_step"] == "uv"
    done_before = list(cp["done"])
    simulator.inject("uvvis", None)
    from delfin.lab import estop
    estop.clear_latch("uvvis")                       # the person's step
    intents_before = _intents("run-r")
    cp = coord.run(WF, "run-r")
    assert cp["status"] == "completed"
    assert cp["done"][:len(done_before)] == done_before
    # nothing completed before the failure was written again
    new = _intents("run-r")[len(intents_before):]
    assert {r["device"] for r in new} == {"uvvis"}


def _intents(run_id):
    from delfin.lab import audit
    return [r for r in audit.read_all()
            if r.get("event") == "intent" and r.get("session") == run_id]


def test_two_coordinators_cannot_hold_one_module(lab):
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    first = _coordinator(manifest, clients, clk, name="c1").run(WF, "run-a")
    assert first["status"] == "waiting_approval"
    second = _coordinator(manifest, clients, clk, person="bob",
                          name="c2").run(WF, "run-b")
    assert second["status"] == "refused"
    assert "held by alice" in second["deviations"][0]["reason"]
    # the agent session cannot write into a module a run holds either
    from delfin.lab.safety import Actor, Kernel
    out = Kernel.load().set(Actor.here("agent-1"), "lamp", "power", "on")
    assert out["code"] == "module_locked"


def test_edge_refuses_bad_mac_replay_and_wrong_manifest(lab):
    from delfin.lab import edge
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    node = edges["edge_analysis"]
    actor = {"person": "alice", "session": "x", "coordinator": "c1"}
    key = edge.load_key("edge_analysis")
    body = {"v": 1, "node": "edge_analysis", "op": "set",
            "args": {"device": "lamp", "channel": "power", "value": "on"},
            "actor": actor, "manifest": manifest.checksum,
            "nonce": "a" * 32, "ts": clk()}
    env = edge.seal(body, key)
    assert node.handle(env)["ok"]                    # control
    out = node.handle(env)                           # the same envelope again
    assert out["code"] == "envelope" and "replay" in out["reason"]
    late = edge.seal({**body, "nonce": "b" * 32, "ts": clk() - 120}, key)
    assert "window" in node.handle(late)["reason"]
    forged = dict(edge.seal({**body, "nonce": "c" * 32}, key))
    forged["args"] = {**forged["args"], "value": "off"}
    assert "MAC" in node.handle(forged)["reason"]
    stale = edge.seal({**body, "nonce": "d" * 32, "manifest": "0" * 64}, key)
    out = node.handle(stale)
    assert out["code"] == "manifest"
    # an e-stop is never refused for a manifest mismatch
    stop = edge.seal({**body, "op": "estop", "args": {"module": "analysis"},
                      "nonce": "e" * 32, "manifest": "0" * 64}, key)
    assert node.handle(stop)["scope"] == "analysis"


def test_an_edge_key_is_not_accepted_by_another_edge(lab):
    from delfin.lab import edge
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    synth_key = edge.load_key("edge_synth")
    body = {"v": 1, "node": "edge_analysis", "op": "set",
            "args": {"device": "lamp", "channel": "power", "value": "on"},
            "actor": {"person": "mallory", "session": "x"},
            "manifest": manifest.checksum, "nonce": "f" * 32, "ts": clk()}
    out = edges["edge_analysis"].handle(edge.seal(body, synth_key))
    assert out["code"] == "envelope" and "MAC" in out["reason"]
    # and an edge refuses devices it does not host, even when signed right
    good = edge.seal({**body, "node": "edge_synth", "args": {
        "device": "lamp", "channel": "power", "value": "on"}}, synth_key)
    assert edges["edge_synth"].handle(good)["code"] == "wrong_node"


def test_partition_goes_safe_on_the_edge(lab):
    from delfin.lab import approvals, watchdog
    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    c = clients["edge_synth"]
    actor = {"person": "alice", "session": "run-p"}
    for dev in ("stirrer", "fume_hood", "hotplate"):
        assert c.call("set", actor, device=dev, channel="power",
                      value="on")["ok"]
    r = c.call("set", actor, device="hotplate", channel="setpoint", value=60,
               unit="degC")
    approvals.grant(r["approval_request"], now=clk())
    assert c.call("set", actor, device="hotplate", channel="setpoint",
                  value=60, unit="degC",
                  approval=r["approval_request"])["ok"]
    k = edges["edge_synth"].kernel
    clk.sleep(25)
    assert c.call("heartbeat", actor)["renewed"] >= 1
    clk.sleep(25)
    assert watchdog.sweep(k) == []                   # control: heartbeats came
    clk.sleep(10)                                    # the network is gone
    assert watchdog.sweep(k) == ["hotplate"]
    from delfin.lab.drivers import simulator
    assert simulator.raw("hotplate")["values"]["setpoint"] == 20
