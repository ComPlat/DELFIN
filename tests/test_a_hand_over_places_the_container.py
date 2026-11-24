"""A hand-over moves the container, so the next module can prove what it read.

The campaign in ``delfin/lab/examples/`` carries four samples from storage
through synthesis and purification into analysis, and each measure step refuses
unless the container at its position is the one the registry expects. A
``transfer`` that only rewrote the registry therefore left every measure step
reading an empty rack and reporting a mix-up -- for a campaign that was carried
out correctly. These tests pin the world side of a hand-over.

Each failure has a control: the same step, carried out properly, is accepted, so
a test cannot pass because the whole path is broken.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import _lab_support  # noqa: F401  (the fixture lives there)
from _lab_support import FakeClock, lab_fixture  # noqa: F401

EXAMPLES = Path(__file__).resolve().parents[1] / "delfin" / "lab" / "examples"
MANIFEST = EXAMPLES / "chemasap_campaign.yaml"


def _install(tmp_path):
    """The campaign manifest, pinned, one edge key per node."""
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


def _coordinator(man):
    import time

    from delfin.lab import edge, workflow
    from delfin.lab.safety import Kernel
    clk = FakeClock(time.time())
    clients = {}
    for node in man.nodes:
        key = edge.load_key(node)
        clients[node] = edge.EdgeClient(
            node, key, transport=edge.EdgeNode(node, key=key, clock=clk).handle,
            manifest_checksum=man.checksum, coordinator="c1", clock=clk)
    Kernel.load(clock=clk)
    return workflow.Coordinator(man, clients, person="alice", name="c1",
                                clock=clk, sleep=clk.sleep)


def _one_sample():
    return {"S-001": {"container": "vial-1", "module": "storage",
                      "position": "ST1"}}


def test_a_hand_over_places_the_container_on_the_destination_position(
        lab, tmp_path):
    """The registry says S-001 is at A1; the rack slot must then hold vial-1.

    Before this test could pass, ``transfer`` changed only the registry, so the
    rack slot still read empty and the analyser refused the sample.
    """
    from delfin.lab import workflow
    from delfin.lab.safety import Kernel
    man = _install(tmp_path)
    coord = _coordinator(man)
    Kernel.load().read("rack_scanner_1")
    cp = {"samples": _one_sample(), "measurements": []}
    out = coord._step({"id": "to-rack", "op": "transfer", "sample": "S-001",
                       "from": "storage", "to": "analysis", "position": "A1"},
                      {"person": "alice", "session": "t"}, cp)
    assert out.get("ok"), out
    # The world, not the report: the slot the analyser reads holds the vial.
    assert out.get("placed_on") == "rack_scanner_1", out
    read = coord.clients[man.devices["rack_scanner_1"].node].call(
        "read", {"person": "alice", "session": "t"}, device="rack_scanner_1")
    assert read["readings"]["container"]["value"] == "vial-1", read


def test_the_measure_step_accepts_a_sample_that_was_handed_over(lab, tmp_path):
    """The point of the move: after it, the analyser proves the right vial.

    The control is in the same test: the SAME measure step is refused for a
    sample whose container was never placed, so a green result cannot come from
    the check having been loosened.
    """
    from delfin.lab.safety import Kernel
    man = _install(tmp_path)
    coord = _coordinator(man)
    actor = {"person": "alice", "session": "t"}
    Kernel.load().read("rack_scanner_1")
    Kernel.load().read("hplc_1")
    mv = coord._step({"id": "to-rack", "op": "transfer", "sample": "S-001",
                      "from": "storage", "to": "analysis", "position": "A1"},
                     actor, {"samples": _one_sample(), "measurements": []})
    assert mv.get("ok"), mv

    cp = {"samples": _one_sample(), "measurements": []}
    cp["samples"]["S-001"].update(module="analysis", position="A1")
    ok = coord._step({"id": "m1", "op": "measure", "device": "hplc_1",
                      "sample": "S-001", "reader": "rack_scanner_1"}, actor, cp)
    assert ok.get("ok"), ok
    assert ok["measurement"]["container"] == "vial-1", ok

    # Control: the vial is taken out of the slot; the same step is refused.
    gone = _coordinator(man)
    gone.clients[man.devices["rack_scanner_1"].node].call(
        "set", actor, device="rack_scanner_1", channel="container",
        value="empty")
    bad = gone._step({"id": "m2", "op": "measure", "device": "hplc_1",
                      "sample": "S-001", "reader": "rack_scanner_1"}, actor,
                     {"samples": {**{"S-001": {**_one_sample()["S-001"],
                                                "module": "analysis",
                                                "position": "A1"}}},
                      "measurements": []})
    assert bad.get("ok") is False and bad.get("code") == "sample_mixup", bad


def test_an_empty_container_name_is_refused(lab, tmp_path):
    """A slot set to "" would read back as an empty rack; the kernel refuses it.

    Driven through ``Kernel.set`` rather than the tool, because the tool's own
    argument check rejects a missing value before the kernel is reached -- and
    the question here is what the KERNEL does with an empty name.
    """
    from delfin.lab.safety import Actor, Kernel
    _install(tmp_path)
    k = Kernel.load()
    me = Actor.here("t", person="alice")
    out = k.set(me, "rack_scanner_1", "container", "")
    assert out.get("ok") is False and out.get("code") == "bad_value", out
    # Control: a real barcode is accepted and reads back, so the refusal is
    # about the value and not about the slot being unwritable.
    ok = k.set(me, "rack_scanner_1", "container", "vial-7")
    assert ok.get("ok") and ok.get("readback") == "vial-7", ok


def test_the_shipped_campaign_declares_a_slot_for_every_rack_position():
    """Every analysis position has exactly one device that can hold a container.

    Without one the hand-over into analysis could only rewrite the registry, and
    the measure step could never be satisfied.
    """
    from delfin.lab import manifest
    man = manifest.load(MANIFEST)
    analysis = man.modules["analysis"]
    for position in analysis.positions:
        holders = [man.devices[d] for d in analysis.devices
                   if (man.devices[d].address or {}).get("position") == position
                   and "container" in man.devices[d].writes]
        assert len(holders) == 1, (position, [h.id for h in holders])
        assert holders[0].writes["container"].text, holders[0].id
