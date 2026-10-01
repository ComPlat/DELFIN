"""Controls for the skill_registry evidence gate (package 8).

Publishing to the team archive is only reachable with a valid evidence
record; invalid evidence is refused with a reason, never repaired and
waved through. Red on the previous commit (publish_skill had no
``evidence`` parameter and no gate).
"""

from __future__ import annotations

from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import skill_registry as sr  # noqa: E402


def _run_ok(calls):
    def run(cmd):
        calls.append(cmd)
        return 0, "", ""
    return run


def _stamped_green_run(ws: Path) -> list[dict]:
    """Green ledger entry stamped with the tree state (the shape
    api_client's _stamp_new_evidence produces in real operation)."""
    from delfin.agent import evidence_freshness as ef
    run = {"command": "tests/test_skill_registry.py", "exit_code": 0,
           "status": "ok", "passed": 1, "failed": 0}
    return [ef.stamp(run, ws)]


_EV_TEST = {"kind": "test", "ref": "tests/test_skill_registry.py"}
_EV_JOB = {"kind": "job", "ref": "12345"}


# ---------------------------------------------------------------------------
# publish: refused without valid evidence
# ---------------------------------------------------------------------------

def test_publish_without_evidence_refused():
    calls: list = []
    ok, msg = sr.publish_skill("casscf-setup", host="h", user="u",
                               remote_path="/r", run_fn=_run_ok(calls))
    assert ok is False and "evidence" in msg.lower()
    assert calls == []                          # nothing left the machine


def test_publish_with_empty_evidence_refused():
    calls: list = []
    ok, msg = sr.publish_skill("casscf-setup", host="h", user="u",
                               remote_path="/r", evidence=[],
                               run_fn=_run_ok(calls))
    assert ok is False and "evidence" in msg.lower()
    assert calls == []


def test_publish_with_only_invalid_evidence_refused():
    calls: list = []
    # a test ref whose file does not exist -> verify_evidence rejects it
    ok, msg = sr.publish_skill("casscf-setup", host="h", user="u",
                               remote_path="/r",
                               evidence=[{"kind": "test",
                                          "ref": "tests/nope.py::t"}],
                               run_fn=_run_ok(calls))
    assert ok is False and "evidence" in msg.lower()
    assert calls == []


def test_publish_with_invalid_evidence_names_the_reason():
    ok, msg = sr.publish_skill("casscf-setup", host="h", user="u",
                               remote_path="/r",
                               evidence=[{"kind": "test",
                                          "ref": "tests/nope.py::t"}],
                               run_fn=_run_ok([]))
    assert ok is False
    assert "not found" in msg or "no green run" in msg


def test_publish_with_valid_evidence_builds_commands():
    calls: list = []
    ws = Path(__file__).resolve().parents[1]
    ok, msg = sr.publish_skill("casscf-setup", host="h", user="u",
                               remote_path="/r",
                               evidence=[{"kind": "test",
                                          "ref": "tests/test_skill_registry.py::test_publish_requires_remote_config"}],
                               workspace=ws,
                               runs=_stamped_green_run(ws),
                               run_fn=_run_ok(calls))
    assert ok is True, msg
    assert len(calls) == 2                      # mkdir + rsync built


# ---------------------------------------------------------------------------
# chemistry domain needs calc-or-test evidence
# ---------------------------------------------------------------------------

def test_publish_chemistry_with_only_job_evidence_refused():
    calls: list = []
    # A job listing that knows the job AND can confirm it ran: still
    # refused, because chemistry demands a calc-or-test evidence entry and
    # a job alone cannot attest it.
    #
    # The stub carries a live state now. Job evidence stopped accepting
    # mere presence in a listing (work/j1-grounding): without a state this
    # case is refused one rule earlier, for "no state to confirm it ran",
    # and never reaches the chemistry rule it exists to test. Weakening
    # the assertion would have hidden that; strengthening the stub keeps
    # the test testing its own subject.
    class _Job:
        job_id = "12345"
        state = "RUNNING"

    ok, msg = sr.publish_skill("casscf-setup", host="h", user="u",
                               remote_path="/r",
                               domain="chemistry",
                               evidence=[_EV_JOB],
                               list_jobs=lambda: [_Job()],
                               run_fn=_run_ok(calls))
    assert ok is False and "chemistry" in msg.lower()
    assert calls == []


# ---------------------------------------------------------------------------
# pull: proposals, never directly active
# ---------------------------------------------------------------------------

def test_pull_returns_proposals_not_installs(tmp_path, monkeypatch):
    """Pull must route every downloaded skill into the proposal pipeline
    (status "pending", source "team"), not into the active skill dir."""
    def run(cmd):
        target = Path(cmd[-1])
        (target / "gamma.md").write_text("# Gamma")
        return 0, "", ""

    created: list[dict] = []

    def fake_propose(name, text):
        created.append({"name": name})
        return "proposed", name

    import delfin.agent.skill_registry as mod
    monkeypatch.setattr(mod, "_propose_pulled_skill", fake_propose,
                        raising=False)
    ok, results = sr.pull_skills(host="h", user="u", remote_path="/r",
                                 dest_dir=tmp_path, run_fn=run)
    assert ok is True
    assert created and created[0]["name"] == "gamma"
    assert any(r[1] == "proposed" for r in results)
