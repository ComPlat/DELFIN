"""Phase 5, real wiring: propose() must honour skill_safety.check.

Runs the REAL proposal store (delfin.agent.skill_proposals, Paket 1)
against the REAL checker (delfin.agent.skill_safety). The public path a
session actually calls is propose() then accept(): an unsafe text must
come back status="blocked" with the findings stored, accept() must
refuse it, and a clean text must stay pending and be acceptable -- with
evidence that passes the real evidence gate (Paket 8), not a stand-in.
"""
from __future__ import annotations

import pytest

from delfin.agent import skill_safety

@pytest.fixture(scope="module")
def sp():
    from delfin.agent import skill_proposals
    return skill_proposals


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    return tmp_path


@pytest.fixture(autouse=True)
def _unwire():
    yield
    skill_safety.wire(refusal_memory=None)


def _ev(sp):
    return [sp.Evidence(kind="test", ref="tests/test_x.py::test_y",
                        detail="green", verified_at="2026-09-27")]


def test_unsafe_text_comes_back_blocked(sp, home):
    skill_safety.wire(refusal_memory=None)
    p = sp.propose("dangerous-cleanup", "```bash\nrm -rf build/\n```",
                   evidence=_ev(sp), source="nacht-s13")
    assert p.status == "blocked"
    assert p.findings, "blocked must carry the findings that blocked it"
    # the findings are what OUR check reported, stored verbatim
    assert p.findings == skill_safety.check("```bash\nrm -rf build/\n```")


def test_accept_refuses_a_blocked_proposal(sp, home):
    skill_safety.wire(refusal_memory=None)
    p = sp.propose("dangerous-cleanup", "```bash\nrm -rf build/\n```",
                   evidence=_ev(sp), source="nacht-s13")
    assert p.status == "blocked"
    with pytest.raises(ValueError, match="blocked"):
        sp.accept("dangerous-cleanup", by="operator")


def test_accept_still_refuses_after_a_reload(sp, home):
    # the refusal must survive a fresh read from disk, not only the
    # in-memory object: a human coming back later must hit the same wall
    skill_safety.wire(refusal_memory=None)
    sp.propose("dangerous-cleanup", "```bash\nrm -rf build/\n```",
               evidence=_ev(sp), source="nacht-s13")
    fresh = sp.get_proposal("dangerous-cleanup")
    assert fresh is not None and fresh.status == "blocked"
    with pytest.raises(ValueError, match="blocked"):
        sp.accept("dangerous-cleanup", by="operator")


def test_clean_text_is_pending_and_acceptable(sp, home):
    skill_safety.wire(refusal_memory=None)
    p = sp.propose("gate-workflow",
                   "Run the covering tests first:\n\n```bash\n"
                   "pytest -q tests/test_x.py\n```\n",
                   evidence=_ev(sp), source="nacht-s13")
    assert p.status == "pending"
    assert p.findings == []
    # accept() re-verifies the evidence (Paket 8): the cited test file
    # exists in the workspace and a green run of it is on record.
    ws = home / "ws"
    (ws / "tests").mkdir(parents=True)
    (ws / "tests" / "test_x.py").write_text("def test_y():\n    pass\n")
    runs = [{"command": "pytest -q tests/test_x.py", "exit_code": 0,
             "status": "ok"}]
    target = sp.accept("gate-workflow", by="operator", workspace=ws,
                       runs=runs)
    assert target.exists()


def test_a_refused_target_blocked_through_the_public_path(sp, home):
    # wire() + propose(): the refusal memory reaches the store through
    # the public path, not only through direct check() calls
    from delfin.agent.refusal_memory import Refusal, RefusalMemory
    mem = RefusalMemory(entries=[Refusal(
        tool="read_file", target="/etc/hosts",
        reason="system file", time="07:00")])
    skill_safety.wire(refusal_memory=mem)
    p = sp.propose("hosts-reader", "```bash\ncat /etc/hosts\n```",
                   evidence=_ev(sp), source="nacht-s13")
    assert p.status == "blocked"
    assert any("refus" in f.lower() for f in p.findings), p.findings


def test_no_evidence_is_no_proposal(sp, home):
    with pytest.raises(ValueError):
        sp.propose("no-evidence", "# Nothing dangerous\n",
                   evidence=[], source="nacht-s13")
