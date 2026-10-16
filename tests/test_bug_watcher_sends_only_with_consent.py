"""The bug watcher sends other people's transcripts to a model provider
only with its own consent.

A report in the archive is another user's session: their conversation,
system prompt with memories, tool outputs. ``triage_report`` hands its
report.md to an engine, i.e. to whichever provider the maintainer's
settings name. Enabling the watcher (``agent.bug_watcher.enabled``) is a
decision about running a daemon; sending the content to a third party is
a separate processing purpose and a transfer to a processor.

``agent.bug_watcher.send_reports``: ``ask`` (default; unknown values
too) sends only with ``consent=True`` -- the maintainer said yes in a
terminal for this run; ``never`` never sends; ``always`` is the standing
opt-in an unattended daemon needs. A withheld report is not marked
triaged, so it is analysed once consent exists.
"""
from __future__ import annotations

import argparse
import json

import pytest

from delfin.agent import bug_watcher as bw


def _report(archive, name="20260101-000000_u_sess_ab12"):
    d = archive / name
    d.mkdir(parents=True)
    (d / "report.json").write_text(json.dumps({
        "schema": "delfin-bug-report/1", "description": "it crashed",
        "model": "m", "mode": "solo"}), encoding="utf-8")
    (d / "report.md").write_text("# Bug report\n\nsomeone's transcript\n",
                                 encoding="utf-8")
    return d


class _Recorder:
    def __init__(self):
        self.sent = []

    def __call__(self, repo_root, settings=None):
        rec = self

        class _E:
            def stream_response(self, user_message="", on_token=None, **kw):
                rec.sent.append(user_message)
                if on_token:
                    on_token("ROOT CAUSE: x")
                return "ROOT CAUSE: x"
        return _E()


def _s(policy=None):
    cfg = {} if policy is None else {"send_reports": policy}
    return {"agent": {"bug_watcher": cfg}}


def test_an_enabled_watcher_alone_sends_nothing(tmp_path):
    d = _report(tmp_path)
    rec = _Recorder()
    t = bw.triage_report(d, settings=_s(), engine_factory=rec)
    assert rec.sent == []
    assert not (d / "triage.md").exists()        # still pending
    assert bw._attempts(d) == 0
    assert "send_reports" in t.summary


def test_never_sends_nothing_even_with_a_yes(tmp_path):
    d = _report(tmp_path)
    rec = _Recorder()
    bw.triage_report(d, settings=_s("never"), engine_factory=rec,
                     consent=True)
    assert rec.sent == []


@pytest.mark.parametrize("policy", [None, "ask", "maybe"])
def test_ask_sends_after_a_yes(tmp_path, policy):
    d = _report(tmp_path)
    rec = _Recorder()
    bw.triage_report(d, settings=_s(policy), engine_factory=rec,
                     consent=True)
    assert len(rec.sent) == 1 and (d / "triage.md").exists()


def test_always_is_the_opt_in_for_the_daemon(tmp_path, monkeypatch):
    _report(tmp_path, "r1")
    _report(tmp_path, "r2")
    monkeypatch.setattr(bw, "announce", lambda *a, **k: None)
    rec = _Recorder()
    bw.run_once(archive=tmp_path, settings=_s("always"), engine_factory=rec)
    assert len(rec.sent) == 2


def test_an_unattended_run_withholds_and_announces_nothing(
        tmp_path, monkeypatch):
    _report(tmp_path)
    told = []
    monkeypatch.setattr(bw, "announce", lambda t, **k: told.append(t))
    rec = _Recorder()
    out = bw.run_once(archive=tmp_path, settings=_s(), engine_factory=rec)
    assert rec.sent == [] and told == []
    assert len(out) == 1 and out[0].withheld


def _cli(tmp_path, monkeypatch, *, tty, answer="", policy=None):
    import delfin.user_settings as us
    from delfin.agent import cli
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: _s(policy))
    rec = _Recorder()
    monkeypatch.setattr(bw, "_default_engine_factory", rec)
    monkeypatch.setattr("sys.stdin.isatty", lambda: tty, raising=False)
    monkeypatch.setattr("builtins.input", lambda _q="": answer)
    args = argparse.Namespace(bug_action="triage", archive=str(tmp_path),
                              name="")
    rc = cli.cmd_bug(args)
    return rc, rec


def test_bug_triage_without_a_terminal_sends_nothing(
        tmp_path, monkeypatch, capsys):
    _report(tmp_path)
    rc, rec = _cli(tmp_path, monkeypatch, tty=False)
    assert rec.sent == []
    assert "send_reports" in capsys.readouterr().out


def test_bug_triage_asks_and_sends_after_a_yes(tmp_path, monkeypatch):
    _report(tmp_path)
    rc, rec = _cli(tmp_path, monkeypatch, tty=True, answer="y")
    assert len(rec.sent) == 1


def test_bug_triage_a_no_sends_nothing(tmp_path, monkeypatch):
    _report(tmp_path)
    rc, rec = _cli(tmp_path, monkeypatch, tty=True, answer="n")
    assert rec.sent == []
