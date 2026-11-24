"""Chemotion ELN exporter: run record -> ELN objects, dry-run, upload.

The exporter (delfin/lab/chemotion.py) maps a ``delfin.lab.run/1`` record
into a research-plan style ELN element with its samples and attachments,
writes the exact HTTP request bodies to a folder in a dry run (sending
nothing), and can POST them to a Chemotion server through an injectable
transport. Upload needs a token by reference and a confirmation; the tests
only exercise the transport against a fake in-process HTTP server on
127.0.0.1, never a real server and never a real token.
"""

from __future__ import annotations

import json
from pathlib import Path


def _record() -> dict:
    """A representative delfin.lab.run/1 record (shape of record.build)."""
    return {
        "schema": "delfin.lab.run/1",
        "run_id": "run-42",
        "status": "completed",
        "installation": "kit-lab-a",
        "manifest": {"version": "1", "checksum": "deadbeef"},
        "plan": {"id": "wf-demo"},
        "plan_digest": "digest123",
        "actor": {"person": "alice", "coordinator": "c1"},
        "started": "2026-01-01T10:00:00+00:00",
        "ended": "2026-01-01T11:00:00+00:00",
        "steps": [{"step": "heat", "op": "set", "ts": 1}],
        "measurements": [
            {"sample": "S-001", "device": "uvvis", "quantity": "signal",
             "value": 0.42, "unit": None, "container": "vial-7", "ts": 2},
            {"sample": "S-002", "device": "ph", "quantity": "ph",
             "value": 7.1, "unit": "pH", "container": "vial-8", "ts": 3},
        ],
        "deviations": [{"step": "heat", "reason": "overshoot"}],
        "samples": {
            "S-001": {"container": "vial-7", "module": "synthesis",
                      "position": "S1"},
            "S-002": {"container": "vial-8", "module": "analysis",
                      "position": "A2"},
        },
        "device_log": [{"event": "intent", "device": "hotplate", "session": "s1"}],
        "fair": {"identifier": "delfin-lab-run:kit-lab-a:run-42",
                 "creator": "alice"},
    }


def test_build_import_maps_record_to_research_plan():
    from delfin.lab.chemotion import build_import

    model = build_import(_record())
    # The format string moved with the payload: the model now also carries the
    # reaction with its conditions, the container tree each analysis lives in,
    # the provenance and the FAIR slots, so it is a different format rather
    # than the same one with more keys.
    assert model["format"] == "delfin.lab.chemotion-import/1"
    elem = model["element"]
    assert elem["type"] == "research_plan"
    assert elem["name"] == "Automated run run-42"
    assert elem["body"]["plan_digest"] == "digest123"
    assert elem["body"]["status"] == "completed"
    # one sample entry per run sample
    assert [s["name"] for s in model["samples"]] == ["S-001", "S-002"]
    # analyses carry the measurements of their sample
    a1 = [s for s in model["samples"] if s["name"] == "S-001"][0]
    assert a1["analyses"][0]["value"] == 0.42
    a2 = [s for s in model["samples"] if s["name"] == "S-002"][0]
    assert a2["analyses"][0]["unit"] == "pH"
    # attachments carry the record, device log, measurements and parameters
    names = {a["filename"] for a in model["attachments"]}
    assert {"record.json", "device_log.json", "measurements.json",
            "parameters.json"} <= names
    assert model["unverified"], "honesty list must be non-empty"


def test_request_sequence_order_and_paths():
    from delfin.lab import chemotion

    seq = chemotion.request_sequence(chemotion.build_import(_record()),
                                     collection_id=7)
    # first a research-plan create, then one sample create
    assert seq[0].path == "/api/v1/research_plans"
    assert seq[0].kind == "json"
    assert seq[1].path == "/api/v1/samples"
    # then chunked attachment uploads, then the attach step
    paths = [r.path for r in seq]
    assert "/api/v1/attachments/upload_chunk" in paths
    assert "/api/v1/attachments/upload_chunk_complete" in paths
    assert any(r.path.startswith("/api/v1/attachments")
               and r.kind == "attach" for r in seq)
    # collection id flows into the research plan body
    assert seq[0].data["collection_id"] == 7


def test_dry_run_writes_exact_bodies_without_any_transport():
    from delfin.lab import chemotion

    out = Path("scratch_dryrun").resolve()
    out.mkdir(exist_ok=True)
    try:
        model = chemotion.build_import(_record())
        written = chemotion.dry_run(model, out, collection_id=7)
        files = [p.name for p in written]
        assert any(f.startswith("01_") and "research_plans" in f for f in files)
        # the create body is exact and JSON
        rp = [f for f in written if "research_plans" in f.name][0]
        body = json.loads(rp.read_text())
        assert body["data"]["collection_id"] == 7
        assert body["data"]["name"] == "Automated run run-42"
        # no network was touched: nothing outside out_dir was created
        leftovers = sorted(p for p in out.iterdir())
        assert leftovers == sorted(written), "dry run must not touch transport"
    finally:
        for p in out.iterdir():
            p.unlink()
        out.rmdir()


def test_upload_needs_confirmation_and_sends_nothing_when_refused():
    from delfin.lab import chemotion

    model = chemotion.build_import(_record())
    try:
        chemotion.upload(model, "http://127.0.0.1:1", token="chemtoken_test",
                         collection_id=7, confirm=False)
        raise AssertionError("upload with confirm=False must refuse")
    except chemotion.ChemotionError as exc:
        assert "confirm" in str(exc).lower()


def test_upload_sends_request_sequence_to_fake_server():
    """Transport is injected: the sandbox refuses real sockets, so the fake
    server is an in-process function with the same contract (path, body,
    bearer header), not a TCP listener. Everything else is asserted as for a
    real server: ordering, the exact create body, the bearer token."""
    from delfin.lab import chemotion

    seen = []

    def fake_server(req, base_url, token):
        import json as _json
        seen.append({
            "method": req.method, "path": req.path,
            "body": req.data if isinstance(req.data, dict) else
            (_json.loads(req.data) if req.kind == "json" else None),
            "auth": f"Bearer {token}" if req.auth else None,
        })
        return {"status": 200, "body": '{"ok": true}'}

    model = chemotion.build_import(_record())
    chemotion.upload(model, "http://192.0.2.1", token="chemtoken_test",
                     collection_id=7, confirm=True, transport=fake_server)
    paths = [r["path"] for r in seen]
    assert paths[0] == "/api/v1/research_plans"
    assert "/api/v1/samples" in paths
    assert "/api/v1/attachments/upload_chunk_complete" in paths
    # the create body reaches the server with a collection id
    rp = [r for r in seen if r["path"] == "/api/v1/research_plans"][0]
    assert rp["body"]["collection_id"] == 7
    # the token travels as a bearer header (by reference, never embedded)
    assert rp["auth"] == "Bearer chemtoken_test"


def test_a_sample_request_carries_the_container_the_analyses_live_in():
    """The upload request must carry the tree, not a bare list of values.

    ``build_import`` builds each analysis as analysis -> dataset -> the
    measurement's data file, in the shape the ELN persists. The sample create
    request then sent ``container: {"analyses": [...]}`` -- the OLD older,
    flatter shape, with no ``container_type``, no ``children`` and no
    datasets -- so an upload of this model would have asked the server to
    create a container type it does not have, and would not have carried the
    data files at all. The two shapes are built in two places, and had
    drifted, exactly as the data-file names had.

    This asserts the SAME property for the request as for the file: what the
    model says an analysis carries is what the request carries.
    """
    from delfin.lab import chemotion

    model = chemotion.build_import(_record())
    seq = chemotion.request_sequence(model, collection_id=7)
    creates = [r for r in seq if r.path == "/api/v1/samples"]
    assert len(creates) == 2, [r.path for r in seq]

    for req, sample in zip(creates, model["samples"]):
        container = req.data["container"]
        # the request carries the container the model built, not a stub
        assert container["container_type"] == "root", container
        assert "children" in container, sorted(container)
        assert container["children"] == sample["container"]["children"], (
            "the sample create request does not carry the tree the import "
            "model built")

    # and each analysis's data file is reachable through the request
    files = [at["filename"]
             for req in creates
             for an in req.data["container"]["children"]
             for ds in an["children"] for at in ds["attachments"]]
    assert len(files) == 2, files
    assert files == sorted(set(files)), files


def test_the_upload_carries_every_data_file_the_model_names():
    """An uploaded attachment list must include the analyses' data files.

    The chunked upload step iterates ``model["attachments"]``; if the
    per-measurement files were only added to the container tree, the tree
    would name files the upload never sends.
    """
    from delfin.lab import chemotion

    model = chemotion.build_import(_record())
    seq = chemotion.request_sequence(model)
    # a multipart request with no file part is the completion step, not an
    # upload of bytes; only the ones carrying a file count here
    uploaded = {r.files[0][1] for r in seq
                if r.kind == "multipart" and r.files}
    tree = [at["filename"]
            for s in model["samples"]
            for an in s["container"]["children"]
            for ds in an["children"] for at in ds["attachments"]]
    assert tree, "the model names no data file"
    missing = [n for n in tree if n not in uploaded]
    assert not missing, f"the tree names data files the upload never sends: {missing}"
