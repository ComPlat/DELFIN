"""A ChemASAP-shaped campaign: four modules, four samples, in parallel.

What this file pins, and why each assertion can fail:

* the shipped manifest really has the four modules, real units and ranges, a
  hazardous channel and a forbidden one -- a manifest that lost its approvals
  or its ``forbidden`` level would stop refusing;
* ``delfin-agent lab run`` is a subcommand and reaches the coordinator
  (before this change nothing in delfin/ called ``Coordinator`` at all);
* four dwells of the same length overlap under ``--parallel 4`` while the
  same document runs them one at a time under the default -- measured as
  wall time, so a scheduler that ignored ``parallel`` cannot pass;
* two samples in one module are measured on the instrument that serves each
  of their positions, and a mis-wired one is refused instead of filing the
  wrong reading under the right sample;
* the kernel refuses an out-of-range, an interlocked and an unapproved
  hazardous step **in bypass mode**, through the tool the model calls.
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import pytest
import yaml

import _lab_support  # noqa: F401  (the fixture lives there)
from _lab_support import call, lab_fixture  # noqa: F401

EXAMPLES = Path(__file__).resolve().parents[1] / "delfin" / "lab" / "examples"
MANIFEST = EXAMPLES / "chemasap_campaign.yaml"
WORKFLOW = EXAMPLES / "chemasap_campaign_workflow.yaml"

MODULES = ("storage", "synthesis", "purification", "analysis")


# --------------------------------------------------------------------------
# the manifest
# --------------------------------------------------------------------------
def _install_campaign(tmp_path):
    """The campaign manifest, pinned, with an edge key per node."""
    from delfin.lab import edge, manifest, store
    from delfin.lab.safety import pin
    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MANIFEST.read_text(encoding="utf-8"), encoding="utf-8")
    man = manifest.load(path)
    pin(man, by="alice")
    for node in man.nodes:
        edge.provision_key(node)
    return man


def test_the_campaign_models_the_four_modules_with_positions_and_units():
    from delfin.lab import manifest
    man = manifest.load(MANIFEST)
    assert tuple(m for m in MODULES if m in man.modules) == MODULES
    for module in MODULES:
        # Every module has to hold four samples at once, or four samples
        # cannot be in it together.
        assert len(man.modules[module].positions) >= 4, module
        assert set(man.modules[module].positions).isdisjoint({"", None})
    # Real units and ranges, not free numbers.
    setpoint = man.devices["reactor_1"].writes["setpoint"]
    assert (setpoint.unit, setpoint.min, setpoint.max) == ("degC", 20, 120)
    flow = man.devices["column_pump_1"].writes["flow"]
    assert (flow.unit, flow.max) == ("mL/min", 10)
    # Both approval levels that make the kernel refuse something are present.
    levels = {man.level(d, c) for d in man.devices.values()
              for c in d.writes.values()}
    assert {"benign", "hazardous", "forbidden"} <= levels
    # The thermocouple is a different device from the heater (design §T15:
    # an interlock should name an independent sensor).
    assert "reactor_probe_1" in man.devices
    assert man.devices["reactor_probe_1"].module == "synthesis"


def test_every_measure_device_declares_the_position_it_serves():
    """A measure step is only meaningful if the instrument says which
    position it reads; without it the coordinator must refuse (see the
    workflow test below)."""
    from delfin.lab import manifest
    man = manifest.load(MANIFEST)
    for name, dev in man.devices.items():
        if "barcode" in dev.reads:
            assert (dev.address or {}).get("position"), name


def test_the_shipped_workflow_is_a_four_sample_campaign_through_all_four():
    wf = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    assert len(wf["samples"]) == 4
    touched = set()
    for st in wf["steps"]:
        for key in ("from", "to"):
            if st.get(key):
                touched.add(st[key])
    assert set(MODULES) <= touched
    # Every sample is transferred through every module and measured once.
    for sid in wf["samples"]:
        moves = [st for st in wf["steps"]
                 if st["op"] == "transfer" and st["sample"] == sid]
        assert len(moves) == 3, sid
        assert [m["to"] for m in moves] == ["synthesis", "purification",
                                            "analysis"]
        assert len([st for st in wf["steps"]
                    if st["op"] == "measure" and st["sample"] == sid]) == 1
    # Dependencies are resolvable and no id is repeated.
    ids = [st["id"] for st in wf["steps"]]
    assert len(ids) == len(set(ids))
    for st in wf["steps"]:
        for a in st.get("after") or []:
            assert a in ids, (st["id"], a)


# --------------------------------------------------------------------------
# the product entry point
# --------------------------------------------------------------------------
def test_the_cli_offers_lab_run_and_reaches_the_coordinator():
    from delfin.lab import cli
    parser = cli.build_parser()
    subs = next(a.choices for a in parser._actions
                if getattr(a, "choices", None) and a.dest == "cmd")
    assert "run" in subs
    args = parser.parse_args(["run", str(WORKFLOW), "--parallel", "4"])
    assert args.cmd == "run" and args.parallel == 4 and args.workflow


def test_the_cli_refuses_a_workflow_with_a_missing_step_id(tmp_path):
    """A refusal in the loader, reached through the argument the CLI takes."""
    from delfin.lab import cli
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump({"id": "x", "steps": [{"op": "hold"}]}),
                   encoding="utf-8")
    with pytest.raises(ValueError, match="unique id"):
        cli._load_workflow(str(bad))


# --------------------------------------------------------------------------
# parallelism, measured as wall time through cli.main
# --------------------------------------------------------------------------
def _lanes(par: int, dwell: float = 4.0) -> dict:
    steps = []
    for i in range(1, 5):
        steps += [
            {"id": f"open-{i}", "op": "set", "device": f"storage_drawer_{i}",
             "channel": "power", "value": "open"},
            {"id": f"dwell-{i}", "op": "hold", "seconds": dwell,
             "after": [f"open-{i}"]},
            {"id": f"shut-{i}", "op": "set", "device": f"storage_drawer_{i}",
             "channel": "power", "value": "closed", "after": [f"dwell-{i}"]},
        ]
    return {"id": f"lanes-{par}", "samples": {}, "steps": steps}


def _time_a_run(tmp_path, par: int, dwell: float = 4.0) -> float:
    from delfin.lab import cli, workflow
    _install_campaign(tmp_path)
    wf = _lanes(par, dwell)
    path = tmp_path / f"lanes-{par}.yaml"
    path.write_text(yaml.safe_dump(wf, sort_keys=False), encoding="utf-8")
    t0 = time.time()
    rc = cli.main(["run", str(path), "--parallel", str(par),
                   "--run-id", f"lanes-{par}"])
    wall = time.time() - t0
    cp = workflow.load_checkpoint(f"lanes-{par}")
    assert rc == 0 and cp["status"] == "completed", (rc, cp.get("deviations"))
    # All four lanes really ran: 4 x (open, dwell, shut).
    assert len(cp["done"]) == 12
    return wall


def test_four_lanes_overlap_under_parallel_and_not_under_the_default(
        lab, tmp_path_factory):
    """The measurement, not the mechanism: four 4 s dwells cost about one
    dwell at --parallel 4 and about four at the default. A scheduler that
    ignored ``parallel`` gives the same number twice and fails the second
    assertion; one that ran nothing gives a tiny number and fails the first.

    Each run gets its own lab: the kernel locks a device for the session that
    wrote it and never hands it to another, so the second run on the same lab
    would be refused for the first run's lock -- right, but not what is being
    measured here.
    """
    dirs = []
    for par in (4, 1):
        d = tmp_path_factory.mktemp(f"dwell-{par}")
        (d / "home").mkdir()
        (d / "lab").mkdir()
        dirs.append(d)
    saved = {k: os.environ.get(k) for k in ("HOME", "DELFIN_LAB_DIR")}

    def at(d):
        os.environ["HOME"] = str(d / "home")
        os.environ["DELFIN_LAB_DIR"] = str(d / "lab")

    try:
        at(dirs[0])
        par = _time_a_run(dirs[0], 4)
        at(dirs[1])
        seq = _time_a_run(dirs[1], 1)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    # Parallel: one dwell's worth plus overhead, comfortably below two.
    assert par < 2 * 4.0, f"parallel run took {par:.1f}s for 4 x 4s dwells"
    # The control: the same document, one at a time, costs all four.
    assert seq > 3 * 4.0, f"sequential run took {seq:.1f}s for 4 x 4s dwells"
    assert seq > 2 * par, f"parallel {par:.1f}s vs sequential {seq:.1f}s"


# --------------------------------------------------------------------------
# sample hand-over: the right instrument, the right sample
# --------------------------------------------------------------------------
def _two_samples_in_analysis(tmp_path):
    from delfin.lab import edge, manifest, store, workflow
    from delfin.lab.drivers import simulator
    from delfin.lab.safety import Kernel
    from _lab_support import FakeClock
    _install_campaign(tmp_path)
    man = manifest.load(store.manifest_path())
    clk = FakeClock(time.time())
    clients = {}
    for node in man.nodes:
        key = edge.load_key(node)
        clients[node] = edge.EdgeClient(
            node, key, transport=edge.EdgeNode(node, key=key, clock=clk).handle,
            manifest_checksum=man.checksum, coordinator="c1", clock=clk)
    k = Kernel.load(clock=clk)
    for n in (1, 2):
        for dev in (f"rack_scanner_{n}", f"hplc_{n}"):
            k.read(dev)
            simulator.place(dev, barcode=f"vial-{n}")
        simulator.place(f"hplc_{n}", signal_next=0.30 + 0.10 * n)
    return man, clients, workflow.Coordinator(man, clients, person="alice",
                                              name="c1", clock=clk,
                                              sleep=clk.sleep)


def test_two_samples_in_one_module_are_measured_on_their_own_instruments(
        lab, tmp_path):
    man, clients, coord = _two_samples_in_analysis(tmp_path)
    cp = {"samples": {"S-001": {"container": "vial-1", "module": "analysis",
                                "position": "A1"},
                      "S-002": {"container": "vial-2", "module": "analysis",
                                "position": "A2"}},
          "measurements": []}
    actor = {"person": "alice", "session": "t"}
    out = {}
    for sid, dev in (("S-001", "hplc_1"), ("S-002", "hplc_2")):
        r = coord._step({"id": f"m-{sid}", "op": "measure", "device": dev,
                         "sample": sid, "reader": f"rack_scanner_{dev[-1]}"},
                        actor, cp)
        assert r.get("ok"), (sid, r)
        out[sid] = r["measurement"]["value"]
    # The two readings are distinct, so neither sample got the other's vial.
    assert out == {"S-001": pytest.approx(0.40), "S-002": pytest.approx(0.50)}


def test_an_instrument_serving_another_position_is_refused(lab, tmp_path):
    """The guard: S-002 is at A2, but it is asked for on A1's instrument. The
    reading would be vial-1's, filed under S-002."""
    man, clients, coord = _two_samples_in_analysis(tmp_path)
    cp = {"samples": {"S-002": {"container": "vial-2", "module": "analysis",
                                "position": "A2"}},
          "measurements": []}
    r = coord._step({"id": "m", "op": "measure", "device": "hplc_1",
                     "sample": "S-002", "reader": "rack_scanner_1"},
                    {"person": "alice", "session": "t"}, cp)
    assert r.get("ok") is False and r.get("code") == "position", r
    assert "serves position A1" in r["reason"]
    assert cp["measurements"] == []


def test_an_instrument_without_a_declared_position_is_refused(lab, tmp_path):
    """No address.position means the coordinator cannot show the instrument
    serves that sample's position at all, so it refuses rather than accept a
    barcode that only happens to match."""
    from delfin.lab import edge, manifest, store, workflow
    from delfin.lab.drivers import simulator
    from delfin.lab.safety import Kernel
    from _lab_support import FakeClock
    doc = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    for mod in doc["modules"]:
        for dev in mod["devices"]:
            dev.pop("address", None)
    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
    man = manifest.parse_document(doc)
    clk = FakeClock(time.time())
    clients = {}
    for node in man.nodes:
        key = edge.provision_key(node)
        clients[node] = edge.EdgeClient(
            node, key, transport=edge.EdgeNode(node, key=key, clock=clk).handle,
            manifest_checksum=man.checksum, coordinator="c1", clock=clk)
    k = Kernel.load(clock=clk)
    k.read("hplc_1")
    k.read("rack_scanner_1")
    simulator.place("hplc_1", barcode="vial-1")
    simulator.place("rack_scanner_1", barcode="vial-1")
    coord = workflow.Coordinator(man, clients, person="alice", name="c1",
                                 clock=clk, sleep=clk.sleep)
    cp = {"samples": {"S-001": {"container": "vial-1", "module": "analysis",
                                "position": "A1"}},
          "measurements": []}
    r = coord._step({"id": "m", "op": "measure", "device": "hplc_1",
                     "sample": "S-001", "reader": "rack_scanner_1"},
                    {"person": "alice", "session": "t"}, cp)
    assert r.get("ok") is False and r.get("code") == "instrument", r
    assert "declares no address.position" in r["reason"]
    assert cp["measurements"] == []


# --------------------------------------------------------------------------
# the safety kernel, in bypass mode, through the tool the model calls
# --------------------------------------------------------------------------
def test_bypass_mode_does_not_move_any_limit(lab):
    """Out of range, interlocked and unapproved-hazardous are all refused
    while the permission mode is bypassPermissions -- the kernel does not
    read the mode (design doc §4)."""
    _install_campaign(lab)
    # The interlock is about heating a reactor without stirring or extraction,
    # and the approval is about the hazardous write itself. Both are needed
    # for each refusal to be about the thing it names.
    r = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                              "value": 121, "unit": "degC"})
    assert r.get("code") == "out_of_range", r

    r = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                              "value": 80, "unit": "degC"})
    assert r.get("code") == "interlock", r

    for dev in ("stirrer_1", "fume_hood"):
        assert call(lab, "lab_set", {"device": dev, "channel": "power",
                                     "value": "on"}).get("ok")
    r = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                              "value": 80, "unit": "degC"})
    assert r.get("code") == "approval_required", r
    assert r.get("approval_request"), r
    # ... and the mode was bypass the whole time.
    assert r.get("mode") is None or r.get("mode") == "bypassPermissions"

    r = call(lab, "lab_set", {"device": "mains_bypass", "channel": "power",
                              "value": "on"})
    assert r.get("code") == "forbidden", r

    # A benign write in the same mode does go through: the refusals above are
    # about the limits, not about bypass refusing everything.
    r = call(lab, "lab_set", {"device": "storage_drawer_1", "channel": "power",
                              "value": "open"})
    assert r.get("ok"), r


def test_the_pump_flow_is_hazardous_and_interlocked(lab):
    _install_campaign(lab)
    r = call(lab, "lab_set", {"device": "column_pump_1", "channel": "flow",
                              "value": 2, "unit": "mL/min"})
    assert r.get("code") == "interlock", r          # valve still closed
    assert call(lab, "lab_set", {"device": "column_valve_1", "channel": "power",
                                 "value": "on"}).get("ok")
    r = call(lab, "lab_set", {"device": "column_pump_1", "channel": "flow",
                              "value": 2, "unit": "mL/min"})
    assert r.get("code") == "approval_required", r


# --------------------------------------------------------------------------
# through the agent's own lab tools, on the simulator driver
# --------------------------------------------------------------------------
def test_the_agent_tools_list_and_read_the_campaign_devices(lab):
    """``lab_devices`` / ``lab_read`` against the campaign manifest, through
    the executor a model's call goes through (api_client._DocToolExecutor)."""
    _install_campaign(lab)
    out = call(lab, "lab_devices", {})
    ids = {d["id"] for d in out["devices"]}
    assert {"reactor_1", "column_pump_1", "hplc_1", "storage_drawer_1"} <= ids
    r1 = next(d for d in out["devices"] if d["id"] == "reactor_1")
    assert r1["writes"]["setpoint"]["level"] == "hazardous"
    out = call(lab, "lab_read", {"device": "reactor_probe_1"})
    assert out["readings"]["temperature"]["unit"] == "degC"


def test_the_agent_can_drive_one_lane_but_not_past_the_limits(lab):
    """The four-step lane a model would run, on the simulator driver: drawer,
    stirrer, hood, then an approved heat that the interlock still guards."""
    _install_campaign(lab)
    assert call(lab, "lab_set", {"device": "storage_drawer_1", "channel":
                                 "power", "value": "open"}).get("ok")
    assert call(lab, "lab_set", {"device": "stirrer_1", "channel": "power",
                                 "value": "on"}).get("ok")
    assert call(lab, "lab_set", {"device": "fume_hood", "channel": "power",
                                 "value": "on"}).get("ok")
    # Approved by a person, then repeated with the grant: the write lands.
    asked = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                                  "value": 80, "unit": "degC"})
    assert asked.get("code") == "approval_required", asked
    grant = asked["approval_request"]
    from delfin.lab import approvals
    approvals.grant(grant, approver="alice")
    done = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                                 "value": 80, "unit": "degC",
                                 "approval": grant})
    assert done.get("ok") and done["readback"] == pytest.approx(80), done
    # The grant is single use: the same id again is refused.
    again = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                                  "value": 80, "unit": "degC",
                                  "approval": grant})
    assert again.get("ok") is False and again.get("code") == "approval", again


def test_a_grant_for_one_reactor_does_not_heat_another(lab):
    """Approval is bound to device, channel, value and session. A person
    approving reactor_1 at 80 degC must not thereby arm reactor_2."""
    _install_campaign(lab)
    for i in (1, 2):
        assert call(lab, "lab_set", {"device": f"stirrer_{i}", "channel":
                                     "power", "value": "on"}).get("ok")
    assert call(lab, "lab_set", {"device": "fume_hood", "channel": "power",
                                 "value": "on"}).get("ok")
    asked = call(lab, "lab_set", {"device": "reactor_1", "channel": "setpoint",
                                  "value": 80, "unit": "degC"})
    from delfin.lab import approvals
    approvals.grant(asked["approval_request"], approver="alice")
    out = call(lab, "lab_set", {"device": "reactor_2", "channel": "setpoint",
                                "value": 80, "unit": "degC",
                                "approval": asked["approval_request"]})
    assert out.get("ok") is False, out
    assert out.get("code") in ("approval", "approval_required"), out
