"""The audit chain is keyed: a rewrite without the account key is found.

Input: records appended with ``audit_log.append`` to a log under
``tmp_path``; the trust store redirected beside it. The log is then
rewritten the way a same-account process without the key would: a line
changed and every later ``prev`` and the head file recomputed so the
unkeyed chain is consistent; every ``mac`` stripped and the chain
recomputed; the tail cut and an earlier head file put back. Output:
``audit_log.verify`` and ``delfin-agent audit verify``, with and without
the key available.

Semantics:

* every record written with the key carries ``mac``, HMAC-SHA256 under
  ``trust_store.derived_key("audit_log")`` of the line without it, as its
  last field;
* the trust store keeps an anchor per log (last line hash and size,
  signed); a log that does not end where the anchor says is reported;
* records of a log written before the key existed are counted as
  unkeyed, not as tampering; an unkeyed record after a keyed one is;
* where the key is not available (inside an agent session) verify still
  checks the hash chain and says the signatures were not checked.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import audit_log as AL
from delfin.agent import security_events as SE

MARKER = "DELFIN_AGENT_SESSION"


@pytest.fixture(autouse=True)
def store(tmp_path, monkeypatch):
    monkeypatch.delenv(MARKER, raising=False)
    from delfin.agent import trust_store
    monkeypatch.setattr(trust_store, "store_dir",
                        lambda: tmp_path / "trust_store")
    trust_store.forget_cached_keys()
    SE.clear()
    monkeypatch.setattr(AL, "_WRITE_FAILURE_NOTED", False, raising=False)
    yield tmp_path / "trust_store"
    trust_store.forget_cached_keys()
    SE.clear()

def _now_ts() -> str:
    """A timestamp in the current ISO week. The log rotates weekly when the
    first line is from an earlier week; a fixed date made these lines rotate
    away once that week had passed (main went red on 2026-10-19)."""
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")



def _log(tmp_path):
    return tmp_path / "state" / "audit.log"


def _write(path, n=4, start=0):
    for i in range(start, start + n):
        AL.append(AL.make_record(tool="bash", decision="auto",
                                 command=f"echo {i}"), log_path=path)


def _lines(path):
    return path.read_bytes().rstrip(b"\n").split(b"\n")


def _rechain(path, lines):
    """What a process without the key can do: make the unkeyed chain and
    the head consistent again after its change."""
    out, prev = [], None
    for raw in lines:
        rec = json.loads(raw)
        if prev is not None:
            rec["prev"] = prev
        raw = json.dumps(rec, ensure_ascii=False).encode()
        out.append(raw)
        prev = AL._line_hash(raw)
    data = b"\n".join(out) + b"\n"
    path.write_bytes(data)
    AL._write_head(path, prev, len(data))


def test_every_record_is_signed_and_verifies(tmp_path):
    path = _log(tmp_path)
    _write(path)
    for raw in _lines(path):
        rec = json.loads(raw)
        assert len(rec["mac"]) == 64
        assert raw.endswith(f', "mac": "{rec["mac"]}"}}'.encode())
    rep = AL.verify(path)
    assert rep["ok"], rep
    assert rep["keyed"] == 4


def test_a_rewrite_that_recomputes_the_unkeyed_chain_is_found(tmp_path):
    path = _log(tmp_path)
    _write(path)
    lines = _lines(path)
    lines[1] = lines[1].replace(b"echo 1", b"echo 9")
    _rechain(path, lines)
    rep = AL.verify(path)
    assert not rep["ok"]
    assert any(p["kind"] == "mac" and p["line"] == 2
               for p in rep["problems"]), rep


def test_stripping_every_signature_is_found(tmp_path):
    path = _log(tmp_path)
    _write(path)
    lines = []
    for raw in _lines(path):
        rec = json.loads(raw)
        rec.pop("mac")
        lines.append(json.dumps(rec, ensure_ascii=False).encode())
    _rechain(path, lines)
    rep = AL.verify(path)
    assert not rep["ok"]
    assert any(p["kind"] == "anchor" for p in rep["problems"]), rep


def test_a_cut_tail_with_an_earlier_head_put_back_is_found(tmp_path):
    path = _log(tmp_path)
    _write(path, 2)
    old_head = AL.head_path(path).read_bytes()
    old_size = path.stat().st_size
    _write(path, 2, start=2)
    with open(path, "r+b") as fh:
        fh.truncate(old_size)
    AL.head_path(path).write_bytes(old_head)
    rep = AL.verify(path)
    assert not rep["ok"]
    assert any(p["kind"] == "anchor" for p in rep["problems"]), rep


def test_an_unkeyed_chain_from_before_the_key_continues(tmp_path, store):
    """Migration: the chain written before this change has prev and no
    mac; it is counted, and the keyed records continue it."""
    path = _log(tmp_path)
    path.parent.mkdir(parents=True)
    prev = AL.GENESIS
    lines = []
    for i in range(3):
        raw = json.dumps({"ts": _now_ts(), "tool": "bash",
                          "command": f"old {i}", "prev": prev}).encode()
        lines.append(raw)
        prev = AL._line_hash(raw)
    data = b"\n".join(lines) + b"\n"
    path.write_bytes(data)
    AL._write_head(path, prev, len(data))
    _write(path, 2)
    rep = AL.verify(path)
    assert rep["ok"], rep
    assert rep["unkeyed"] == 3 and rep["keyed"] == 2


def test_an_unkeyed_record_after_a_keyed_one_is_found(tmp_path,
                                                     monkeypatch):
    path = _log(tmp_path)
    _write(path, 2)
    monkeypatch.setenv(MARKER, "77")
    _write(path, 1, start=2)
    monkeypatch.delenv(MARKER)
    rep = AL.verify(path)
    assert not rep["ok"]
    assert any(p["kind"] == "unkeyed" and p["line"] == 3
               for p in rep["problems"]), rep


def test_verify_without_the_key_checks_the_chain_and_says_so(tmp_path,
                                                            monkeypatch):
    path = _log(tmp_path)
    _write(path)
    monkeypatch.setenv(MARKER, "77")
    rep = AL.verify(path)
    assert rep["ok"], rep
    assert any("not checked" in n for n in rep["notes"]), rep
    lines = _lines(path)
    lines[0] = lines[0].replace(b"echo 0", b"echo 8")
    path.write_bytes(b"\n".join(lines) + b"\n")
    assert not AL.verify(path)["ok"]


def test_the_cli_reports_a_keyless_rewrite(tmp_path, capsys):
    from delfin.agent import cli as agent_cli
    path = _log(tmp_path)
    _write(path)
    assert agent_cli.main(["audit", "verify", "--log", str(path)]) == 0
    out = capsys.readouterr().out
    assert "4 signed" in out, out
    lines = _lines(path)
    lines[2] = lines[2].replace(b"echo 2", b"echo 7")
    _rechain(path, lines)
    assert agent_cli.main(["audit", "verify", "--log", str(path)]) == 1


def test_a_rotation_keeps_the_signatures(tmp_path):
    path = _log(tmp_path)
    _write(path, 1)
    raw = _lines(path)[0]
    rec = json.loads(raw)
    rec.pop("mac")
    rec["ts"] = "2020-01-01T00:00:00Z"
    from delfin.agent import trust_store
    body = json.dumps(rec, ensure_ascii=False)
    tag = trust_store.mac("audit_log", body.encode())
    line = (body[:-1] + f', "mac": "{tag}"}}').encode()
    path.write_bytes(line + b"\n")
    AL._write_head(path, AL._line_hash(line), len(line) + 1)
    AL._write_anchor(path, AL._line_hash(line), len(line) + 1)
    _write(path, 1)
    rotated = path.with_name("audit-2020-W01.log")
    assert rotated.exists()
    reps = AL.verify_all(path.parent)
    assert all(r["ok"] for r in reps), reps
