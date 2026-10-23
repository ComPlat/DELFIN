"""A scheduled watchpost run is the deterministic scan, never a model turn.

The first version of `watchpost enable` scheduled an ordinary entry whose
prompt said "watchpost run: scan the account ...". The scheduler daemon
runs ordinary entries as an agent turn: every 30 minutes a model would
have interpreted "scan the account" in its own way, at a token cost, with
a toolbox. The entry now carries a machine marker the daemon runs as the
read-only scan itself, like the benchmark entries.
"""

from __future__ import annotations

from types import SimpleNamespace

from delfin.agent import scheduler_daemon as sd
from delfin.watchpost import scheduler_entry
from delfin.watchpost.model import Finding


def test_the_entry_prompt_carries_the_machine_marker():
    assert scheduler_entry.PROMPT_MARKER == sd.WATCHPOST_ENTRY_PREFIX


def test_a_watchpost_entry_runs_the_scan_and_no_engine(tmp_path, monkeypatch):
    scans = []

    def fake_scan(home, report, mode="report", git_repos=None):
        scans.append(mode)
        return [Finding(check="ssh", severity="alert",
                        path="~/.ssh/authorized_keys", line=2,
                        what="new key with command= option",
                        why="a forced command runs on every login")]

    monkeypatch.setattr("delfin.watchpost.scan.run_scan", fake_scan)
    events = []

    def emit(kind, **kw):
        # As the real one does: a kind the inbox does not take is refused.
        # A stand-in that took any kind let "watchpost_alert" through, and
        # the daemon swallowed the real refusal: no alert ever arrived.
        from delfin.agent.attention import ATTENTION_KINDS
        if kind not in ATTENTION_KINDS:
            raise ValueError(f"unknown attention kind {kind!r}")
        events.append(kind)

    monkeypatch.setattr("delfin.agent.attention.emit_attention", emit)

    def no_engine(*a, **k):
        raise AssertionError("a watchpost entry must not build an engine")

    fire = sd.make_fire_callback(engine_factory=no_engine, log=lambda m: None)
    entry = SimpleNamespace(
        id="wp1", kind="interval", workspace=str(tmp_path), reason="",
        prompt=f"{sd.WATCHPOST_ENTRY_PREFIX} diff scan of the account")
    fire(entry)
    assert scans == ["diff"]
    assert "security_alert" in events
