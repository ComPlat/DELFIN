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


def test_explicit_setting_beats_remote_path(monkeypatch):
    monkeypatch.delenv("DELFIN_BUG_ARCHIVE", raising=False)
    got = br.resolve_archive_dir({
        "agent": {"bug_archive_dir": "/explicit"},
        "transfer": {"remote_path": "/home/grp/archive"},
    })
    assert str(got) == "/explicit"


def test_derives_from_transfer_remote_path(monkeypatch):
    monkeypatch.delenv("DELFIN_BUG_ARCHIVE", raising=False)
    got = br.resolve_archive_dir({"transfer": {"remote_path": "/home/grp/archive"}})
    assert str(got) == "/home/grp/archive/AGENT_BUGS"


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
