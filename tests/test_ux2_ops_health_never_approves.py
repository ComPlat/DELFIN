"""UX2: `delfin-agent ops health` never acts on anything.

It is read-only by design (classify + collect read state; only --wait
EMITS attention events). These tests prove the read-only contract at the
two places the system enforces acts:
- ops health is NOT an act in approval_guard._OPS_ACTS (the gate's list
  of ops subcommands that act);
- running the command writes no approval file and no pause flag (the
  snapshot of every file under the redirected HOME is identical before
  and after, attention inbox excepted for --wait, which we don't use).
"""

import json

import pytest

from delfin.agent import approval_guard
from delfin.agent import cli
from delfin.agent import session_health as health_mod
from delfin.agent import session_pause as pause_mod
from delfin.agent.state_paths import RedirectedUserState


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    with RedirectedUserState(tmp_path):
        yield tmp_path


def _tree_snapshot(root):
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            out[str(p.relative_to(root))] = p.read_bytes()
    return out


def test_ops_health_is_not_an_act_in_approval_guard():
    assert "health" not in approval_guard._OPS_ACTS
    # and the acts the guard knows never grew an approve/answer shape
    assert approval_guard._OPS_ACTS <= {"send", "escalate", "pause"}


def test_ops_health_run_writes_nothing_under_home(home, capsys, monkeypatch):
    # seed one pause flag so collect has something to read
    pause_mod.pause("nacht-s29", "brake")
    before = _tree_snapshot(home)
    code = cli.main(["ops", "health", "--json"])
    assert code == 0
    rows = json.loads(capsys.readouterr().out)
    assert isinstance(rows, list)
    after = _tree_snapshot(home)
    # NOTHING changed: no approval file, no pause flag touched, no inbox
    # write (no --wait), no lock file left behind.
    assert after == before


def test_ops_health_collect_never_calls_pause_or_resume(monkeypatch):
    calls = []
    monkeypatch.setattr(pause_mod, "pause",
                        lambda *a, **k: calls.append("pause"))
    monkeypatch.setattr(pause_mod, "resume",
                        lambda *a, **k: calls.append("resume"))
    import delfin.agent.ops_interface as oi
    monkeypatch.setattr(oi, "pending_approvals", lambda: [])
    monkeypatch.setattr(oi, "sessions", lambda: [])
    rows = health_mod.collect()
    assert rows == []
    assert calls == []
