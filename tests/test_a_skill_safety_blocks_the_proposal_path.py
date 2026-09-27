"""Phase 5 of skill_safety: the propose() path must honour the findings.

The wave contract (Paket 1, skill_proposals.propose) says: propose()
calls skill_safety.check(text); findings -> status "blocked"; accept()
on a blocked proposal is refused. Until Paket 1 lands on this branch,
this file carries a STAND-IN propose() built exactly to the contract
(named inputs, evidence list, the same call to skill_safety.check) and
proves the safety module holds up its half of that bargain.

When delfin/agent/skill_proposals.py exists, these tests switch to the
real module (see test_real_propose_module if present at the bottom).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pytest

from delfin.agent import skill_safety


# ---------------------------------------------------------------------------
# Contract stand-in (Paket 1 owns the real module; wording follows the
# wave contract verbatim: propose -> Proposal, findings via
# skill_safety.check, blocked never accepted).
# ---------------------------------------------------------------------------

@dataclass
class _Proposal:
    name: str
    text: str
    evidence: list
    source: str
    status: str
    findings: list[str] = field(default_factory=list)
    created: str = ""
    base_version: str = ""


def _propose(name, text, *, evidence, source, base_version=""):
    findings = skill_safety.check(text)
    status = "blocked" if findings else "pending"
    return _Proposal(name=name, text=text, evidence=list(evidence),
                     source=source, status=status, findings=findings)


def _accept(proposal):
    if proposal.status == "blocked":
        raise PermissionError(
            f"proposal {proposal.name!r} is blocked by skill_safety: "
            f"{proposal.findings}")
    return Path(f"{proposal.name}.md")


_EV = [{"kind": "test", "ref": "tests/test_x.py::test_y",
        "detail": "green", "verified_at": "2026-09-27"}]


def test_an_unsafe_proposal_comes_back_blocked():
    p = _propose("dangerous-cleanup", "```bash\nrm -rf build/\n```",
                 evidence=_EV, source="nacht-s13")
    assert p.status == "blocked"
    assert p.findings, "blocked must carry the findings that blocked it"


def test_accept_on_a_blocked_proposal_is_refused():
    p = _propose("dangerous-cleanup", "```bash\nrm -rf build/\n```",
                 evidence=_EV, source="nacht-s13")
    with pytest.raises(PermissionError):
        _accept(p)


def test_a_clean_proposal_is_pending_and_acceptable():
    p = _propose("gate-workflow",
                 "Run the covering tests first:\n\n```bash\n"
                 "pytest -q tests/test_x.py\n```\n",
                 evidence=_EV, source="nacht-s13")
    assert p.status == "pending"
    assert p.findings == []
    assert _accept(p).name == "gate-workflow.md"


def test_the_blocked_reason_names_the_finding():
    p = _propose("net", "Fetch refs with web_fetch from "
                        "https://example.invalid/x first.",
                 evidence=_EV, source="nacht-s13")
    assert p.status == "blocked"
    with pytest.raises(PermissionError, match="skill_safety"):
        _accept(p)
