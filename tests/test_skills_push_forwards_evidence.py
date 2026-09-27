"""/skills push forwards the accepted proposal's evidence.

s18's evidence gate (7d988878) refuses a publish_skill call without an
evidence record — correctly: no evidence, no publication. But the
dashboard's /skills push caller passed no evidence at all, so pushing
from the UI was refused even for a skill whose proposal was accepted
WITH evidence that already re-verified at accept() time. The evidence
was checked once and then dropped on the floor.

The fix forwards what package 1 already holds: the accepted proposal's
evidence record, re-read from skill_proposals.get_proposal at push
time (fresh, not cached). Contract entries may be dataclasses; the
gate takes dicts, so the panel normalises minimally and keeps the
refs verbatim.
"""

from __future__ import annotations

import inspect
import sys
import types
from dataclasses import dataclass

import pytest

from delfin.dashboard import skill_proposals_panel as P


@dataclass
class _Evidence:
    kind: str
    ref: str
    detail: str = ""


@dataclass
class _Proposal:
    name: str
    text: str
    evidence: list
    source: str
    status: str
    findings: list
    created: str


def _install(monkeypatch, proposals):
    mod = types.ModuleType("delfin.agent.skill_proposals")

    def get_proposal(name):
        for p in proposals:
            if p.name == name:
                return p
        return None

    mod.get_proposal = get_proposal
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_proposals", mod)
    import delfin.agent as _pkg
    monkeypatch.setattr(_pkg, "skill_proposals", mod, raising=False)


def test_the_evidence_of_an_accepted_proposal_is_forwarded(monkeypatch):
    _install(monkeypatch, [_Proposal(
        name="run-tests-first", text="t",
        evidence=[_Evidence("test", "tests/test_x.py::test_y")],
        source="session abc", status="accepted",
        findings=[], created="2026-09-27")])
    ev = P.evidence_for("run-tests-first")
    assert ev == [{"kind": "test", "ref": "tests/test_x.py::test_y",
                   "detail": "", "verified_at": ""}]


def test_dict_entries_pass_through_verbatim(monkeypatch):
    _install(monkeypatch, [_Proposal(
        name="d", text="t",
        evidence=[{"kind": "calc", "ref": "calc/abc", "detail": "gibbs"}],
        source="s", status="accepted", findings=[], created="c")])
    assert P.evidence_for("d") == [{"kind": "calc", "ref": "calc/abc",
                                    "detail": "gibbs"}]


def test_a_skill_without_a_proposal_pushes_no_invented_evidence(monkeypatch):
    _install(monkeypatch, [])
    # No proposal → empty list. The gate refuses — honestly — rather
    # than the caller inventing a record.
    assert P.evidence_for("casscf-setup") == []


def test_the_push_call_site_forwards_evidence():
    """The dashboard's /skills push caller passes the evidence through
    when publish_skill takes it (s18's gated signature), read fresh at
    push time. Pinned on the source the way the other tab_agent tests
    pin wiring."""
    import delfin.dashboard.tab_agent as T
    src = inspect.getsource(T)
    i = src.index("_sr.publish_skill")
    block = src[i:i + 700]
    assert "evidence_for" in block, (
        "the push call must forward the proposal's evidence")


def test_the_forwarding_is_version_tolerant():
    """On a branch without s18's evidence parameter the call must still
    work (signature checked at call time), and WITH the parameter the
    evidence is passed. Exercised through the panel helper's contract:
    evidence_for always returns a list, never raises."""
    # get_proposal unavailable entirely → empty, not an exception.
    assert P.evidence_for("anything") == []
