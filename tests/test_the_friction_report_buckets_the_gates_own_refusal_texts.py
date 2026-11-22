"""The friction report buckets the gate's own refusal texts.

The audit log's denial records carry the reason the gate handed to the
model: the full sentence, written by ``audit_log.append`` through the
gate's ``_audit_call`` (``reason = payload["error"]``). That sentence does
not contain the ``" was not run: "`` form the bucket patterns were written
around, so the measured leader (auto-allow misses, 29 of the 437 denials
bucketed by hand in 2026-09) landed in "Other refusals", and the variant
that names the reading alternative landed in the WRITE-without-read bucket
because its hint text mentions ``read_file``. A report that puts the most
common denial in "Other" and a read-refusal among write defects reports
nothing, and it reports nothing silently.

This test drives the report over a log the PRODUCT wrote -- records
appended with ``audit_log.append`` (the writer every gate path uses) --
not over a hand-built fixture file.
"""
from __future__ import annotations

import json
from pathlib import Path

from delfin.agent import audit_log, friction_report

#: The gate's own refusal sentences (the texts ``_audit_call`` records as
#: ``reason``). Unisolated auto-allow form: delfin/agent/api_client.py
#: ``_unisolated_fallback_verdict``; reading-hint form: the suffix
#: ``_reading_alternative_hint`` appends; cage hide: the sandbox's
#: DELFIN_CAGE_HIDE denial, observed live in a session of 2026-11-22.
_AUTO_ALLOW_UNISOLATED = (
    "blocked: 'pip install rdkit' is not on the auto-allow list, and this "
    "session runs without isolation (no bubblewrap, Landlock or Seatbelt "
    "holds a command here) with nobody to ask, so only auto-allowed "
    "commands run.")
_AUTO_ALLOW_READING_HINT = (
    "blocked: 'cat /etc/hosts' is not on the auto-allow list. NOTE: the "
    "READING form of 'cat' is auto-allowed: read_file. If that is the "
    "question you were asking, ask it that way.")
_CAGE_HIDE = (
    "blocked: this command would print '/home/x/.delfin'. read denied: "
    "'/home/x/.delfin' is hidden from this run (DELFIN_CAGE_HIDE).")


def _write_audit_log(tmp_path: Path, records: list[dict]) -> Path:
    """A real audit log: the writer's own append path, chain and all."""
    home = tmp_path / "home"
    log = home / ".delfin" / "audit.log"
    for rec in records:
        audit_log.append(rec, log_path=log)
    # The writer ran: the log exists and parses as JSON lines.
    lines = [json.loads(line) for line in log.read_text().splitlines()]
    assert len(lines) == len(records)
    return home


def test_auto_allow_miss_is_bucketed_not_other(tmp_path):
    home = _write_audit_log(tmp_path, [
        {"tool": "bash", "decision": "denied", "command": "pip install rdkit",
         "reason": _AUTO_ALLOW_UNISOLATED},
    ])
    data = friction_report.collect(since_s=0.0, home=home)
    keys = {b["key"] for b in data["buckets"]}
    assert "allowlist_miss" in keys, keys
    assert "other" not in keys, keys


def test_reading_hint_variant_lands_in_the_allowlist_bucket(tmp_path):
    home = _write_audit_log(tmp_path, [
        {"tool": "bash", "decision": "denied", "command": "cat /etc/hosts",
         "reason": _AUTO_ALLOW_READING_HINT},
    ])
    data = friction_report.collect(since_s=0.0, home=home)
    keys = {b["key"] for b in data["buckets"]}
    assert "allowlist_miss" in keys, keys
    # The hint NAMES read_file as the allowed way; it is not a
    # write-without-read refusal.
    assert "read_before_write" not in keys, keys


def test_cage_hide_read_denial_is_a_read_outside_bucket(tmp_path):
    home = _write_audit_log(tmp_path, [
        {"tool": "bash", "decision": "denied",
         "command": "ls $HOME/.delfin", "reason": _CAGE_HIDE},
    ])
    data = friction_report.collect(since_s=0.0, home=home)
    keys = {b["key"] for b in data["buckets"]}
    assert "read_outside_workspace" in keys, keys
    assert "other" not in keys, keys
