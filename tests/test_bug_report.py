"""Tests for the configurable bug-report bundler.

Covers archive-path resolution (env > setting > fallback, never a
hard-coded site path), the collision-free report layout, and that the
run configuration the maintainer needs is captured.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import bug_report as br


# ---------------------------------------------------------------------------
# Archive-path resolution — must be configurable, never hard-coded
# ---------------------------------------------------------------------------

def test_env_var_wins(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_BUG_ARCHIVE", str(tmp_path / "envarch"))
    got = br.resolve_archive_dir({"agent": {"bug_archive_dir": "/should/lose"}})
    assert got == (tmp_path / "envarch")


def test_setting_used_when_no_env(monkeypatch):
    monkeypatch.delenv("DELFIN_BUG_ARCHIVE", raising=False)
    got = br.resolve_archive_dir({"agent": {"bug_archive_dir": "/team/AGENT_BUGS"}})
    assert str(got) == "/team/AGENT_BUGS"


def test_local_path_ignores_remote_path(monkeypatch):
    # The transfer remote_path is an SSH target, NOT a local write path —
    # resolve_archive_dir (local staging) must not use it.
    monkeypatch.delenv("DELFIN_BUG_ARCHIVE", raising=False)
    got = br.resolve_archive_dir({"transfer": {"remote_path": "/home/grp/archive"}})
    assert got == br._FALLBACK_DIR


def test_fallback_is_per_user_local(monkeypatch):
    monkeypatch.delenv("DELFIN_BUG_ARCHIVE", raising=False)
    got = br.resolve_archive_dir({})
    assert got == br._FALLBACK_DIR
    assert ".delfin" in str(got)          # local, not a site path


def test_no_site_path_baked_into_source():
    # The cluster path the user mentioned must never be in the code.
    src = (br.__file__)
    text = open(src, encoding="utf-8").read()
    assert "/home/qmchem_all" not in text


# ---------------------------------------------------------------------------
# Report layout
# ---------------------------------------------------------------------------

def _write(tmp_path, **over):
    kw = dict(
        chat_messages=[{"role": "user", "content": "tu X"},
                       {"role": "assistant", "content": "ok", "role_label": "Solo"}],
        mode="solo", provider="kit", model="qwen3.5", effort="high",
        perms="plan", backend="openai", role="solo_agent",
        session_id="abcdef123456", input_tokens=120, output_tokens=80,
        cost_usd=0.0123, description="Agent erfand ein Keyword",
        archive_dir=str(tmp_path),
    )
    kw.update(over)
    return br.write_bug_report(**kw)


def test_writes_json_and_md(tmp_path):
    d = _write(tmp_path)
    assert (d / "report.json").is_file()
    assert (d / "report.md").is_file()
    assert d.parent == tmp_path


def test_metadata_captured(tmp_path):
    d = _write(tmp_path)
    js = json.loads((d / "report.json").read_text())
    assert js["model"] == "qwen3.5"
    assert js["perms"] == "plan"
    assert js["effort"] == "high"
    assert js["provider"] == "kit"
    assert js["mode"] == "solo"
    assert js["schema"] == "delfin-bug-report/1"
    assert js["chat_messages"]                       # conversation included
    md = (d / "report.md").read_text()
    assert "erfand ein Keyword" in md
    assert "## Konversation" in md


def test_reports_do_not_collide(tmp_path):
    d1 = _write(tmp_path)
    d2 = _write(tmp_path)
    assert d1 != d2                                  # unique sub-dirs
    assert d1.name != d2.name


def test_dirname_is_sortable_and_tagged(tmp_path):
    d = _write(tmp_path, session_id="zzzzzzzzzzzz")
    # <UTC-ts>_<user>_<session8>_<rand>
    parts = d.name.split("_")
    assert len(parts) >= 4
    assert parts[0][:8].isdigit()                    # date prefix sorts
    assert "zzzzzzzz" in d.name                       # session tag present


def test_empty_chat_still_writes_without_crash(tmp_path):
    d = _write(tmp_path, chat_messages=[])
    assert (d / "report.json").is_file()
    js = json.loads((d / "report.json").read_text())
    assert js["chat_messages"] == []


# ---------------------------------------------------------------------------
# Remote push
# ---------------------------------------------------------------------------

def test_push_noop_without_transfer_config(tmp_path):
    d = _write(tmp_path)
    ok, msg = br.push_report_to_remote(d, host="", user="", remote_path="")
    assert ok is False
    assert "kein Remote" in msg


def test_push_runs_mkdir_then_rsync(tmp_path, monkeypatch):
    d = _write(tmp_path)
    calls = []

    class _OK:
        returncode = 0
        stderr = ""

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return _OK()

    monkeypatch.setattr(br.subprocess, "run", fake_run)
    ok, where = br.push_report_to_remote(
        d, host="login.cluster", user="ka", remote_path="/home/grp/archive",
        port=22,
    )
    assert ok is True
    assert where == f"/home/grp/archive/AGENT_BUGS/{d.name}"
    assert len(calls) == 2                       # mkdir, then rsync
    # rsync command carries the local report dir as a source
    assert any(str(d) in " ".join(map(str, c)) for c in calls)


def test_debug_fields_captured(tmp_path):
    d = _write(
        tmp_path,
        system_prompt="You are the DELFIN solo agent. Pattern 5: ...",
        error_text="Traceback (most recent call last):\n  ValueError: boom",
        denied_commands=["rm -rf /", "/orca set foo bar"],
    )
    js = json.loads((d / "report.json").read_text())
    assert "DELFIN solo agent" in js["system_prompt"]
    assert "ValueError: boom" in js["error_text"]
    assert "rm -rf /" in js["denied_commands"]
    # and the human report surfaces them
    md = (d / "report.md").read_text()
    assert "## Fehler / Traceback" in md
    assert "## Geblockte Commands" in md
    assert "## System-Prompt" in md


def test_settings_snapshot_has_no_secrets():
    snap = br.settings_snapshot({
        "agent": {"model": "sonnet", "bug_archive_dir": "/x"},
        "runtime": {"backend": "slurm"},
        "transfer": {"host": "h", "user": "u", "remote_path": "/r",
                     "port": 22, "password": "SECRET", "ssh_key": "KEY"},
    })
    assert snap["runtime_backend"] == "slurm"
    assert snap["transfer"]["host"] == "h"
    # secrets must not survive the snapshot
    flat = json.dumps(snap)
    assert "SECRET" not in flat
    assert "KEY" not in flat
    assert "password" not in snap["transfer"]


def test_push_reports_rsync_failure(tmp_path, monkeypatch):
    d = _write(tmp_path)

    class _Mk:
        returncode = 0
        stderr = ""

    class _Fail:
        returncode = 23
        stderr = "rsync: connection refused"

    seq = [_Mk(), _Fail()]
    monkeypatch.setattr(br.subprocess, "run", lambda *a, **k: seq.pop(0))
    ok, msg = br.push_report_to_remote(
        d, host="h", user="u", remote_path="/r", port=22,
    )
    assert ok is False
    assert "rsync fehlgeschlagen" in msg
