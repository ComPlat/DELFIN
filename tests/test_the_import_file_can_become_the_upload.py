"""The import file a person receives must be enough to make the upload.

`record.write` leaves `chemotion_import.json` beside the zip. A person who
takes that twin to the ELN -- or feeds it to `chemotion.request_sequence` --
has nothing else. If the twin drops an attachment's text, the upload cannot be
built from it: `request_sequence` raises `KeyError: 'text'`, which is exactly
what a real run against the base did (`scratch/realrun/no-network-base.jsonl`).

The model here is a real finished run's, taken through `Coordinator.run`, so
the assertions hold on the artefact a person has rather than on a hand-made
dict.
"""

from __future__ import annotations

import json

import pytest

from _lab_support import FakeClock, lab_fixture  # noqa: F401

from test_lab_edge_and_workflow import WF, _approve_pending, _place_vial, _setup


@pytest.fixture(name="finished_run")
def finished_run(lab):
    import time

    from delfin.lab.workflow import Coordinator, run_dir

    clk = FakeClock(time.time())
    manifest, edges, clients = _setup(clk)
    _place_vial(clients)
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
        "twin": json.loads((d / "chemotion_import.json").read_text()),
    }


def test_the_written_twin_can_be_turned_back_into_the_upload(finished_run):
    """The twin alone builds the request sequence, with no other input."""
    from delfin.lab import chemotion

    twin = finished_run["twin"]
    assert twin["attachments"], "a finished run attaches nothing"

    # Every attachment keeps its bytes in the twin on disk.
    for a in twin["attachments"]:
        assert a.get("text"), f"{a['filename']} carries no text in the twin"

    # And the upload can be built from that file, with no model in hand.
    seq = chemotion.request_sequence(twin)
    assert seq, "no requests would be sent"


def test_every_file_the_twin_names_is_a_member_of_the_bundle(finished_run):
    """An attachment the twin names but the zip lacks cannot be uploaded."""
    import zipfile

    twin = finished_run["twin"]
    with zipfile.ZipFile(finished_run["dir"] / "chemotion_import.zip") as z:
        members = set(z.namelist())
    for a in twin["attachments"]:
        assert f"attachments/{a['filename']}" in members, (
            f"{a['filename']} is named in the twin but absent from the bundle")


def test_the_dry_run_writes_a_request_for_every_attachment(finished_run, tmp_path):
    """The reviewable artefact shows every attachment's bytes would be sent.

    The bytes are chunked (5000 B per request) and JSON-escaped in the written
    files, so they are reassembled from the parsed requests rather than
    searched for in the raw text.
    """
    from delfin.lab import chemotion

    twin = finished_run["twin"]
    written = chemotion.dry_run(twin, tmp_path / "requests")
    assert written, "a dry run that writes nothing uploads nothing"

    sent = ""
    for p in written:
        for f in json.loads(p.read_text(encoding="utf-8")).get("files", []):
            sent += f["bytes(utf-8)"]
    for a in twin["attachments"]:
        assert a["text"] in sent, (
            f"{a['filename']}'s bytes appear in no request the dry run sends")
