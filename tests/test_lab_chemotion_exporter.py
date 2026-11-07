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
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from threading import Thread


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
    assert model["format"] == "delfin.lab.chemotion/1"
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
        assert any(f.startswith("1_") and "research_plans" in f for f in files)
        # the create body is exact and JSON
        rp = [f for f in written if "research_plans" in f.name][0]
        body = json.loads(rp.read_text())
        assert body["collection_id"] == 7
        assert body["name"] == "Automated run run-42"
        # no network was touched: nothing outside out_dir was created
        leftovers = sorted(p for p in out.iterdir())
        assert leftovers == sorted(written), "dry run must not touch transport"
    finally:
        for p in out.iterdir():
            p.unlink()
        out.rmdir()


class _Recorder(BaseHTTPRequestHandler):
    """fakes a Chemotion server on 127.0.0.1 and records every request."""

    requests = []

    def log_message(self, *a):
        pass

    def _record(self):
        length = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(length) if length else b""
        _Recorder.requests.append({
            "method": self.command, "path": self.path, "body": body,
            "auth": self.headers.get("Authorization")})
        out = b'{"ok": true}'
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(out)))
        self.end_headers()
        self.wfile.write(out)

    do_GET = do_POST = do_PUT = do_DELETE = _record


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
    from delfin.lab import chemotion

    _Recorder.requests = []
    server = HTTPServer(("127.0.0.1", 0), _Recorder)
    base = f"http://127.0.0.1:{server.server_port}"
    t = Thread(target=server.serve_forever, daemon=True)
    t.start()
    try:
        model = chemotion.build_import(_record())
        chemotion.upload(model, base, token="chemtoken_test",
                         collection_id=7, confirm=True)
        paths = [r["path"] for r in _Recorder.requests]
        assert paths[0] == "/api/v1/research_plans"
        assert "/api/v1/samples" in paths
        assert "/api/v1/attachments/upload_chunk_complete" in paths
        # the create body reaches the server with a collection id
        rp = [r for r in _Recorder.requests
              if r["path"] == "/api/v1/research_plans"][0]
        assert json.loads(rp["body"])["collection_id"] == 7
        # the token travels as a bearer header (by reference, never embedded)
        assert rp["auth"] == "Bearer chemtoken_test"
    finally:
        server.shutdown()
        t.join()
