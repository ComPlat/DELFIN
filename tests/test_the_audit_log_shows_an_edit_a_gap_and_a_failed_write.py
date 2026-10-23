"""The audit log is a hash chain: an edit, a removed line and a cut tail
are found, and a write that fails is said.

Input: records appended with ``audit_log.append`` to a log under
``tmp_path``; the file then edited, shortened or replaced by a link.
Output: ``audit_log.verify`` (and ``delfin-agent audit verify``), the
file's mode and the security events.

Semantics:

* every record carries ``prev``, the SHA-256 of the previous line as it
  was written (``audit_log.GENESIS`` for the first line of a new log;
  the rotated log's last line for the first line after a rotation);
* a head file beside the log (0600) holds the hash and size of the last
  line written, so a tail cut off the log is found as well;
* the log is created 0600 in a 0700 directory and opened with O_APPEND
  and O_NOFOLLOW: a link in its place is not written through;
* a failed write records one blocked ``audit_write_failed`` event per
  process instead of vanishing;
* lines written before the chain existed are counted as legacy, not as
  tampering.

The security-event ring holds 200 events; how many it dropped is counted
and shown, so a flood does not hide what came before it.
"""
from __future__ import annotations

import json
import os
import stat

import pytest

from delfin.agent import audit_log as AL
from delfin.agent import security_events as SE


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    SE.clear()
    monkeypatch.setattr(AL, "_WRITE_FAILURE_NOTED", False, raising=False)
    yield
    SE.clear()

def _now_ts() -> str:
    """A timestamp in the current ISO week. The log rotates weekly when the
    first line is from an earlier week; a fixed date made these lines rotate
    away once that week had passed (main went red on 2026-10-19)."""
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")



def _log(tmp_path):
    return tmp_path / "state" / "audit.log"


def _write(path, n=4):
    for i in range(n):
        AL.append(AL.make_record(tool="bash", decision="auto",
                                 command=f"echo {i}"), log_path=path)


def _lines(path):
    return path.read_text(encoding="utf-8").splitlines()


def test_an_intact_log_verifies(tmp_path):
    path = _log(tmp_path)
    _write(path)
    rep = AL.verify(path)
    assert rep["ok"], rep
    assert rep["chained"] == 4 and rep["legacy"] == 0
    recs = [json.loads(line) for line in _lines(path)]
    assert recs[0]["prev"] == AL.GENESIS
    assert all(len(r["prev"]) == 64 for r in recs)


def test_the_log_is_owner_only_from_the_start(tmp_path):
    path = _log(tmp_path)
    _write(path, 1)
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert stat.S_IMODE(path.parent.stat().st_mode) == 0o700
    head = AL.head_path(path)
    assert head.exists() and stat.S_IMODE(head.stat().st_mode) == 0o600


def test_an_edited_line_is_found(tmp_path):
    path = _log(tmp_path)
    _write(path)
    lines = _lines(path)
    lines[1] = lines[1].replace("echo 1", "echo X")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    rep = AL.verify(path)
    assert not rep["ok"]
    assert any(p["line"] == 3 for p in rep["problems"]), rep


def test_a_removed_line_is_found(tmp_path):
    path = _log(tmp_path)
    _write(path)
    lines = _lines(path)
    del lines[1]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    rep = AL.verify(path)
    assert not rep["ok"] and rep["problems"], rep


def test_a_cut_tail_is_found(tmp_path):
    path = _log(tmp_path)
    _write(path)
    lines = _lines(path)
    path.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
    rep = AL.verify(path)
    assert not rep["ok"]
    assert any("head" in p["kind"] for p in rep["problems"]), rep


def test_lines_from_before_the_chain_are_legacy(tmp_path):
    path = _log(tmp_path)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"ts": _now_ts(), "tool": "bash"})
                    + "\n", encoding="utf-8")
    _write(path, 2)
    rep = AL.verify(path)
    assert rep["ok"], rep
    assert rep["legacy"] == 1 and rep["chained"] == 2


def test_a_rotation_carries_the_chain(tmp_path):
    path = _log(tmp_path)
    path.parent.mkdir(parents=True)
    path.write_text("", encoding="utf-8")
    AL.append({"tool": "bash", "decision": "auto"}, log_path=path)
    # Age the first line into an earlier ISO week; the head follows the
    # rewritten bytes so the old log is intact before the rotation.
    lines = _lines(path)
    rec = json.loads(lines[0])
    rec["ts"] = "2020-01-01T00:00:00Z"
    # The signature covered the old timestamp; this record now stands as
    # one written before the key (tests/test_the_audit_chain_is_keyed.py
    # covers a signed rotation).
    rec.pop("mac", None)
    path.write_text(json.dumps(rec, ensure_ascii=False) + "\n", encoding="utf-8")
    AL._write_head(path, AL._line_hash(path.read_bytes().rstrip(b"\n")),
                   path.stat().st_size)
    AL.append({"tool": "bash", "decision": "auto"}, log_path=path)
    rotated = path.with_name("audit-2020-W01.log")
    assert rotated.exists()
    assert AL.verify(rotated)["ok"]
    rep = AL.verify(path)
    assert rep["ok"], rep
    first = json.loads(_lines(path)[0])
    assert first["prev"] == AL._line_hash(rotated.read_bytes().rstrip(b"\n"))


def test_a_link_in_place_of_the_log_is_not_written_through(tmp_path):
    path = _log(tmp_path)
    path.parent.mkdir(parents=True, mode=0o700)
    target = tmp_path / "elsewhere.txt"
    target.write_text("unchanged\n", encoding="utf-8")
    os.symlink(target, path)
    _write(path, 1)
    assert target.read_text(encoding="utf-8") == "unchanged\n"
    assert "audit_write_failed" in [e.kind for e in SE.recent(50)]


def test_a_failed_write_is_said_once(tmp_path):
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("x", encoding="utf-8")
    path = blocker / "audit.log"
    _write(path, 3)
    events = [e for e in SE.recent(50) if e.kind == "audit_write_failed"]
    assert len(events) == 1 and events[0].blocked


def test_the_cli_verifies(tmp_path, capsys):
    from delfin.agent import cli as agent_cli
    path = _log(tmp_path)
    _write(path)
    assert agent_cli.main(["audit", "verify", "--log", str(path)]) == 0
    assert "intact" in capsys.readouterr().out.lower()
    lines = _lines(path)
    lines[0] = lines[0].replace("echo 0", "echo Y")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert agent_cli.main(["audit", "verify", "--log", str(path)]) == 1


def test_the_event_kinds_have_labels():
    for kind in ("audit_write_failed", "audit_chain_broken"):
        assert kind in SE.known_kinds()


def test_the_event_ring_counts_what_it_dropped():
    for i in range(250):
        SE.record("deny_pattern", "bash", f"x{i}")
    n = SE.counts()
    assert n["total"] == 200 and n["dropped"] == 50
    assert "50" in SE.format_panel_html()
