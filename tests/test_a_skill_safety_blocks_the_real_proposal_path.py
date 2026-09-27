"""Phase 5, real wiring: propose() must honour skill_safety.check.

Paket 1 (nacht-s11, agent/s11-lw11, commits 6eb3d906..4744cba8)
implements skill_proposals.propose/accept against the contract. This
file runs the REAL store against the REAL checker (this branch's
delfin/agent/skill_safety.py) by loading s11's module from its commit
— read-only, nothing of Paket 1 lands on this branch from here.

The public path a session actually calls is propose() then accept():
an unsafe text must come back status="blocked" with the findings
stored, accept() must refuse it, and a clean text must stay pending
and be acceptable. A check that cannot even run blocks (fail closed,
Paket 1's own rule — pinned here from the outside too).
"""
from __future__ import annotations

import subprocess
import sys

import pytest

from delfin.agent import skill_safety

_S11_BRANCH = "agent/s11-lw11"


def _load_skill_proposals():
    """Load delfin.agent.skill_proposals from the s11 commit, as package
    module ``delfin.agent.skill_proposals`` so its relative imports
    (``from . import skill_safety``) resolve against THIS branch."""
    text = subprocess.run(
        ["git", "show", f"{_S11_BRANCH}:delfin/agent/skill_proposals.py"],
        capture_output=True, text=True, check=True).stdout
    # Relative import needs a package context; register under the real
    # package name so ``from . import skill_safety`` finds OUR module.
    import delfin.agent  # noqa: F401
    mod = type(sys)("delfin.agent.skill_proposals")
    mod.__package__ = "delfin.agent"
    mod.__file__ = f"<git show {_S11_BRANCH}:delfin/agent/skill_proposals.py>"
    # dataclasses looks the module up in sys.modules by __module__ DURING
    # exec, so it must be registered before, not after.
    sys.modules["delfin.agent.skill_proposals"] = mod
    exec(compile(text, mod.__file__, "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def sp():
    return _load_skill_proposals()


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
    target = sp.accept("gate-workflow", by="operator")
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
