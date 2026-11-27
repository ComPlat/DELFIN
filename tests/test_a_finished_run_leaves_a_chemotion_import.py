"""A finished run leaves one file a person imports into Chemotion.

The run reaches this through the coordinator (``record.write`` is called by
``Coordinator.run``), not by calling the exporter. What is asserted is the
file: an import bundle that carries the record, the reaction with its
conditions, one sample whose analysis holds its own data file, the audit
chain the run produced, and the FAIR metadata -- in the shapes the ELN
persists (``container.children``: analyses -> datasets -> attachments).
"""

from __future__ import annotations

import hashlib
import json

import pytest

from _lab_support import FakeClock, lab_fixture  # noqa: F401

from test_lab_edge_and_workflow import WF, _approve_pending, _place_vial, \
    _setup


@pytest.fixture(name="finished_run")
def finished_run(lab):
    """The demo campaign, run to completion through the coordinator."""
    import time

    from delfin.lab.workflow import run_dir

    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    _place_vial(clients)
    from delfin.lab.workflow import Coordinator
    coord = Coordinator(manifest, clients, person="alice", name="c1",
                        clock=clk, sleep=clk.sleep)
    cp = coord.run(WF, "run-1")
    assert cp["status"] == "waiting_approval"
    _approve_pending(cp, clk)
    cp = coord.run(WF, "run-1")
    assert cp["status"] == "completed", cp["deviations"]
    d = run_dir("run-1")
    return {
        "dir": d,
        "record": json.loads((d / "record.json").read_text()),
        "export": json.loads((d / "chemotion_import.json").read_text()),
        "bundle": d / "chemotion_import.zip",
    }


def test_a_finished_run_leaves_an_import_bundle(finished_run):
    """The file exists, is a zip, and its manifest lists every member."""
    import zipfile

    p = finished_run["bundle"]
    assert p.is_file(), "a finished run must leave the import file itself"
    with zipfile.ZipFile(p) as z:
        names = z.namelist()
        assert "manifest.json" in names
        man = json.loads(z.read("manifest.json"))
    assert man["format"] == "delfin.lab.chemotion-import/1"
    assert sorted(man["members"]) == sorted(names)


def test_the_bundle_carries_the_whole_record_not_a_summary(finished_run):
    """Every field the record holds reaches the plan element's body.

    The body is where a reader of the ELN opens the run, so a summary there
    (status, digest, deviations -- what the old export wrote) loses the
    audit chain, the actor and the manifest pin on the element itself.
    """
    import zipfile

    rec = finished_run["record"]
    exp = finished_run["export"]
    body = exp["element"]["body"]
    for field in ("deviations", "device_log", "plan_digest", "manifest",
                  "actor", "fair", "started", "ended", "installation",
                  "measurements", "samples", "steps", "status"):
        assert field in rec, field
        assert field in body, f"record field {field!r} never reaches the body"
    # the same record also rides whole as the bundle's record member, so a
    # reader who wants the bytes rather than the plan uses that
    with zipfile.ZipFile(finished_run["bundle"]) as z:
        assert json.loads(z.read("attachments/record.json")) == rec


def test_the_run_arrives_as_a_reaction_with_its_conditions(finished_run):
    """A reaction, its conditions, and the sample in the reaction's role."""
    exp = finished_run["export"]
    assert len(exp["reactions"]) == 1
    rx = exp["reactions"][0]
    assert rx["name"] == "Automated run run-1"
    # a condition is a set step: what the run HELD a device at. The wait on
    # the probe and the measure are not conditions, the transfer is where the
    # sample went.
    conds = {c["step"]: c for c in rx["conditions"]}
    # Five set steps SUCCEEDED. The refused attempt at "heat" is a step of its
    # own in the record (code approval_required), so an unfiltered list has
    # six entries and one of them names no device at all -- this count is what
    # fails when the refusal is counted as a condition.
    assert set(conds) == {"stir", "hood", "hp_on", "heat", "cool"}
    assert len(rx["conditions"]) == 5
    assert conds["hp_on"]["device"] == "hotplate"
    assert conds["hp_on"]["value"] == "on"
    # the setpoint carries the unit it was set in, which "on" has none of
    assert conds["heat"]["value"] == 60
    assert conds["heat"]["unit"] == "degC"
    # a condition always names the device it held: the refused attempt the
    # kernel logs as its own step carries device None and value None, and
    # would arrive here as a condition of nothing if it were not filtered
    for c in rx["conditions"]:
        assert c["device"], c
        assert c["value"] is not None, c
    # the sample is the reaction's product, with the module it ended in
    assert [p["sample"] for p in rx["products"]] == ["S-001"]
    assert rx["products"][0]["container"] == "vial-7"


def test_an_analysis_carries_its_data_file_in_the_shape_the_eln_persists(
        finished_run):
    """analysis -> dataset -> attachment is the tree the server walks."""
    exp = finished_run["export"]
    sample = [s for s in exp["samples"] if s["name"] == "S-001"][0]
    analyses = sample["container"]["children"]
    assert len(analyses) == 1
    an = analyses[0]
    assert an["container_type"] == "analysis"
    assert an["is_new"] is True and an["attachments"] == []
    dataset = an["children"][0]
    assert dataset["container_type"] == "dataset"
    # the measurement values live in the dataset's metadata ...
    md = dataset["extended_metadata"]
    assert md["delfin"]["value"] == 0.42
    assert md["delfin"]["unit"] == "au"
    assert md["delfin"]["quantity"] == "signal"
    # ... and the measurement itself is a file attached to that dataset
    files = dataset["attachments"]
    assert len(files) == 1
    f = files[0]
    assert f["is_new"] is True and f["filename"].endswith(".json")
    assert f["checksum"] == hashlib.md5(f["raw"].encode("utf-8")).hexdigest()


def test_provenance_names_the_pin_the_chain_and_who_approved(finished_run):
    """Manifest pin, the audit chain head, and the approval, all present."""
    rec = finished_run["record"]
    exp = finished_run["export"]
    prov = exp["provenance"]
    assert prov["manifest"]["checksum"] == rec["manifest"]["checksum"]
    assert prov["plan_digest"] == rec["plan_digest"]
    assert prov["installation"] == rec["installation"]
    # every audit line of the run, and the chain intact across them
    assert len(prov["audit"]) == len(rec["device_log"])
    assert prov["audit_chain_intact"] is True
    assert prov["record_sha256"] == hashlib.sha256(
        json.dumps(rec, indent=1, sort_keys=True).encode("utf-8")).hexdigest()
    # who approved what: the one hazardous write of this run
    assert len(prov["approvals"]) == 1
    ap = prov["approvals"][0]
    assert ap["device"] == "hotplate" and ap["value"] == 60.0
    assert ap["requested_by"] == "alice"
    assert ap["approved_by"] == "alice"
    assert ap["consumed_at"], "an approval the run used must show it was used"


def test_the_bundle_states_its_fair_metadata(finished_run):
    """Findable/accessible/interoperable/reusable fields a repository asks."""
    fair = finished_run["export"]["fair"]
    assert fair["identifier"] == "delfin-lab-run:demo-lab:run-1"
    assert fair["creator"] == "alice"
    assert fair["created"]
    assert fair["provenance"]["manifest_checksum"]
    # the licence is the record's own slot, reported as unset rather than
    # invented -- an unset licence is what a repository must refuse
    assert "license" in fair and fair["license"] is None
    assert fair["license_set"] is False
    assert fair["metadata_standard"] == "delfin.lab.chemotion-import/1"


def test_the_import_says_what_it_could_not_verify(finished_run):
    """Nothing unverified is presented as verified."""
    unv = finished_run["export"]["unverified"]
    assert unv, "the honesty list must stay non-empty"
    # the request shapes verified in the ELN source are no longer claimed
    # unverified; the file-import layout still is
    joined = " ".join(unv)
    assert "container" not in joined.lower() or "children" not in joined.lower()
    assert any("import" in u.lower() or "zip" in u.lower() for u in unv)


def test_two_samples_of_one_quantity_keep_their_own_file_names():
    """A file an analysis carries must be tellable from its twin's.

    Two samples measured for the same quantity are two analyses in two
    containers. If their data files are named "<run>-<quantity>.json" both
    are called the same thing, and a reader holding one cannot say which
    sample it belongs to.
    """
    from delfin.lab import chemotion

    rec = {
        "schema": "delfin.lab.run/1", "run_id": "run-2s",
        "status": "completed", "installation": "demo-lab",
        "manifest": {"version": 1, "checksum": "abc"},
        "plan": {"id": "wf-2s"}, "plan_digest": "d",
        "actor": {"person": "alice"}, "started": None, "ended": None,
        "steps": [], "deviations": [], "device_log": [],
        "measurements": [
            {"sample": "S-001", "device": "uvvis", "quantity": "signal",
             "value": 0.42, "unit": "au", "container": "vial-7", "ts": 1},
            {"sample": "S-002", "device": "uvvis", "quantity": "signal",
             "value": 0.11, "unit": "au", "container": "vial-8", "ts": 2},
        ],
        "samples": {
            "S-001": {"container": "vial-7", "module": "a", "position": "A1"},
            "S-002": {"container": "vial-8", "module": "a", "position": "A2"},
        },
        "fair": {"identifier": "delfin-lab-run:demo-lab:run-2s",
                 "creator": "alice", "license": None, "created": None,
                 "provenance": {"software": "DELFIN", "version": "x",
                                "manifest_checksum": "abc"}},
    }
    model = chemotion.build_import(rec)
    names = [d["attachments"][0]["filename"]
             for s in model["samples"]
             for a in s["container"]["children"]
             for d in a["children"]]
    assert len(names) == 2
    assert len(set(names)) == 2, f"two analyses share a file name: {names}"
    # and each name still says which run and which quantity it is
    for n in names:
        assert n.startswith("run-2s-")
        assert n.endswith("-signal.json")


def test_a_repeat_measurement_keeps_the_data_file_the_analysis_names(finished_run):
    """Every data file the container tree names must be IN the bundle.

    a1's two-sample test asserts the tree's names are distinct, which a
    repeat measurement satisfies too: two analyses of ONE sample and ONE
    quantity were both named ``<run>-<sample>-<quantity>.json``, while the
    bundle's per-measurement members carry a ``-<i>`` index. So the analysis
    named a file that the bundle did not contain -- the link the whole
    "analyses with their data files" requirement is about, broken, with a
    green test over it.

    This asserts the LINK, not the distinctness: whatever the tree names, the
    person importing the bundle must be able to open it.
    """
    import zipfile

    d = finished_run["dir"]
    with zipfile.ZipFile(finished_run["bundle"]) as z:
        members = set(z.namelist())
        model = json.loads(z.read("chemotion_import.json"))

    claimed = [at["filename"]
               for s in model["samples"]
               for an in s["container"]["children"]
               for ds in an["children"]
               for at in ds["attachments"]]
    assert claimed, "the tree names no data file at all"
    missing = [c for c in claimed if f"attachments/{c}" not in members]
    assert not missing, (
        f"an analysis names a data file that is not in the bundle: {missing}; "
        f"bundle has {sorted(m for m in members if 'signal' in m)}")


def test_two_measurements_of_one_sample_get_two_files():
    """The case above, without a campaign: one sample, one quantity, twice.

    Red on the arm before the fix for the named reason -- the two analyses
    both name ``run-2m-S-001-signal.json`` and the bundle has ``-1``/``-2``.
    """
    from delfin.lab import chemotion

    rec = {
        "schema": "delfin.lab.run/1", "run_id": "run-2m",
        "status": "completed", "installation": "demo-lab",
        "manifest": {"version": 1, "checksum": "abc"},
        "plan": {"id": "wf"}, "plan_digest": "d",
        "actor": {"person": "alice"}, "started": None, "ended": None,
        "steps": [], "deviations": [], "device_log": [],
        "measurements": [
            {"sample": "S-001", "device": "uvvis", "quantity": "signal",
             "value": 0.42, "unit": "au", "container": "vial-7", "ts": 1},
            {"sample": "S-001", "device": "uvvis", "quantity": "signal",
             "value": 0.47, "unit": "au", "container": "vial-7", "ts": 2},
        ],
        "samples": {"S-001": {"container": "vial-7", "module": "a",
                              "position": "A1"}},
        "fair": {"identifier": "id", "creator": "alice", "license": None,
                 "created": None, "provenance": {}},
    }
    model = chemotion.build_import(rec)
    tree = [at["filename"]
            for s in model["samples"]
            for an in s["container"]["children"]
            for ds in an["children"]
            for at in ds["attachments"]]
    bundle = set(chemotion.import_files(model))
    assert len(tree) == 2, tree
    assert len(set(tree)) == 2, f"two analyses share a file name: {tree}"
    missing = [n for n in tree if f"attachments/{n}" not in bundle]
    assert not missing, f"the tree names files the bundle lacks: {missing}"


def test_a_campaign_that_measures_one_sample_twice_links_each_file(lab):
    """The repeat case, driven through ``Coordinator.run`` (the real path).

    ``test_two_measurements_of_one_sample_get_two_files`` proves the mapping;
    this drives the campaign and reads the ZIP that a run actually left, so
    the defect cannot come back through the mapping's callers instead.
    """
    import time
    import zipfile

    from delfin.lab.workflow import Coordinator, run_dir
    from test_lab_edge_and_workflow import (WF, _approve_pending, _place_vial,
                                            _setup)

    # the demo campaign, with its single measurement taken twice
    wf = json.loads(json.dumps(WF))
    measure = [s for s in wf["steps"] if s["op"] == "measure"]
    wf["steps"].append({**measure[0], "id": "uv2"})

    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    _place_vial(clients)
    coord = Coordinator(manifest, clients, person="alice", name="c1",
                        clock=clk, sleep=clk.sleep)
    cp = coord.run(wf, "run-repeat")
    if cp["status"] == "waiting_approval":
        _approve_pending(cp, clk)
        cp = coord.run(wf, "run-repeat")
    assert cp["status"] == "completed", cp["deviations"]

    d = run_dir("run-repeat")
    with zipfile.ZipFile(d / "chemotion_import.zip") as z:
        members = set(z.namelist())
        model = json.loads(z.read("chemotion_import.json"))

    tree = [at["filename"]
            for s in model["samples"]
            for an in s["container"]["children"]
            for ds in an["children"]
            for at in ds["attachments"]]
    assert len(tree) == 2, tree
    assert len(set(tree)) == 2, tree
    missing = [n for n in tree if f"attachments/{n}" not in members]
    assert not missing, f"the tree names files the bundle lacks: {missing}"
