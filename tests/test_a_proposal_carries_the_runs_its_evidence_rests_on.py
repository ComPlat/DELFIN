"""A skill proposal is accepted later than it is made -- its evidence must
travel with it.

accept() re-verifies test evidence against a ledger of green runs. The
ledger of the session that proposed the skill is gone by the time a human
accepts it in the dashboard or the CLI, so evidence that did not carry
its runs could never be accepted there. The proposal now stores the green
runs (with their tree fingerprint) that its test evidence was observed
on; accept() re-verifies those, and evidence_freshness still refuses them
once the code under them has moved.
"""
from __future__ import annotations

import subprocess

import pytest

from delfin.agent import session_end, skill_learning, skill_proposals as sp


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    return tmp_path


def _workspace(root, *, git=False):
    ws = root / "ws"
    (ws / "tests").mkdir(parents=True)
    (ws / "tests" / "test_x.py").write_text("def test_t():\n    pass\n")
    if git:
        for cmd in (["git", "init", "-q"],
                    ["git", "-c", "user.email=t@t", "-c", "user.name=t",
                     "add", "."],
                    ["git", "-c", "user.email=t@t", "-c", "user.name=t",
                     "commit", "-qm", "init"]):
            subprocess.run(cmd, cwd=ws, check=True, capture_output=True)
    return ws


_CLEAN = "# Run the covering tests\n\n1. Run the test file first.\n"


def _green(ws=None):
    run = {"command": "pytest -q tests/test_x.py", "exit_code": 0,
           "status": "ok"}
    if ws is not None:
        from delfin.agent import evidence_freshness
        run["fingerprint"] = dict(evidence_freshness.fingerprint(ws))
    return run


def test_stored_runs_make_a_later_accept_possible(home):
    ws = _workspace(home)
    ev = sp.Evidence(kind="test", ref="tests/test_x.py::test_t",
                     runs=[_green()])
    p = sp.propose("covering-tests", _CLEAN, evidence=[ev], source="t")
    assert p.status == "pending", p.findings
    # a later session: no ledger passed, only the workspace
    target = sp.accept(p.name, by="operator", workspace=ws)
    assert target.is_file() and target.name == "SKILL.md"


def test_evidence_without_runs_is_still_refused(home):
    ws = _workspace(home)
    ev = sp.Evidence(kind="test", ref="tests/test_x.py::test_t")
    p = sp.propose("no-runs", _CLEAN, evidence=[ev], source="t")
    with pytest.raises(ValueError, match="evidence"):
        sp.accept(p.name, by="operator", workspace=ws)
    assert sp.get_proposal(p.name).status == "pending"


def test_stored_runs_go_stale_when_the_code_moves(home):
    ws = _workspace(home, git=True)
    ev = sp.Evidence(kind="test", ref="tests/test_x.py::test_t",
                     runs=[_green(ws)])
    p = sp.propose("moves-under-it", _CLEAN, evidence=[ev], source="t")
    # the code under the evidence changes before anyone accepts it
    (ws / "tests" / "test_x.py").write_text("def test_t():\n    assert 0\n")
    subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t",
                    "commit", "-qam", "change"], cwd=ws, check=True,
                   capture_output=True)
    with pytest.raises(ValueError, match="stale"):
        sp.accept(p.name, by="operator", workspace=ws)


def test_the_runs_survive_the_round_trip_to_disk(home):
    ev = sp.Evidence(kind="test", ref="tests/test_x.py::test_t",
                     runs=[_green()])
    p = sp.propose("round-trip", _CLEAN, evidence=[ev], source="t")
    back = sp.get_proposal(p.name)
    assert back.evidence[0].runs == [_green()]


def test_learning_attaches_only_the_runs_that_cover_the_evidence():
    ev = skill_learning._evidence_class()(kind="test",
                                          ref="tests/test_x.py::test_t")
    runs = [_green(),
            {"command": "pytest -q tests/test_other.py", "exit_code": 0,
             "status": "ok"},
            {"command": "pytest -q tests/test_x.py", "exit_code": 1,
             "status": "ok"}]
    skill_learning._attach_runs(ev, runs)
    assert ev.runs == [_green()]


def test_session_end_hands_the_saved_test_ledger_to_learning(monkeypatch):
    from delfin.agent import session_store
    monkeypatch.setattr(session_store, "load_session",
                        lambda sid: {"evidence": {"tests": [_green()]}}
                        if sid == "abc" else None)
    seen = {}

    def _fake_learn(msgs, *, settings=None, llm=None, _propose=None,
                    runs=None):
        seen["runs"] = runs
        return None

    monkeypatch.setattr(skill_learning, "learn_from_session", _fake_learn)
    session_end.learn_at_session_end([], session_id="abc")
    assert seen["runs"] == [_green()]
    assert session_end._session_test_runs("missing") == []


def test_publishing_re_verifies_the_runs_the_evidence_carries(home):
    from delfin.agent import skill_registry
    ws = _workspace(home)
    ev = {"kind": "test", "ref": "tests/test_x.py::test_t",
          "runs": [_green()]}
    ok, detail = skill_registry._gate_evidence([ev], workspace=ws)
    assert ok, detail
    ok, detail = skill_registry._gate_evidence(
        [{"kind": "test", "ref": "tests/test_x.py::test_t"}], workspace=ws)
    assert not ok
