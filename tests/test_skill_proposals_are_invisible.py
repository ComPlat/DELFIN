"""A proposal is invisible until a human accepts it.

The Hermes failure mode (arXiv 2608.12851) is activation without review:
this pins the opposite. A skill sitting in ``_proposals`` must not be
discoverable through any public path -- ``discover_skills`` and the
``_session_skills`` wrapper the ``skill`` tool goes through -- neither
while pending nor after rejection, and only acceptance (which moves the
folder into the live skills directory) makes it appear.
"""
import sys
import types

import pytest

from delfin.agent import skill_proposals as sp
from delfin.agent.skills import discover_skills


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    safe = types.ModuleType("delfin.agent.skill_safety")
    safe.check = lambda text: []
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", safe)
    # ``from . import skill_safety`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "skill_safety",
                        safe, raising=False)
    ev = types.ModuleType("delfin.agent.evidence")
    ev.verify_evidence = lambda e, **kw: (True, "verified")
    monkeypatch.setitem(sys.modules, "delfin.agent.evidence", ev)
    # ``from . import evidence`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "evidence",
                        ev, raising=False)
    return tmp_path


def _standin_safety(monkeypatch):
    mod = types.ModuleType("delfin.agent.skill_safety")
    mod.check = lambda text: []
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", mod)
    # ``from . import skill_safety`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "skill_safety",
                        mod, raising=False)


def _propose(home, name="proposed-skill", text="# Proposed\nbody"):
    return sp.propose(name, text,
                      evidence=[sp.Evidence(kind="test",
                                            ref="tests/test_a.py::t")],
                      source="test")


def test_a_pending_proposal_is_not_discoverable(home, monkeypatch):
    _standin_safety(monkeypatch)
    _propose(home)
    names = [s.name for s in discover_skills()]
    assert "proposed-skill" not in names


def test_a_rejected_proposal_is_not_discoverable(home, monkeypatch):
    _standin_safety(monkeypatch)
    _propose(home)
    sp.reject("proposed-skill", reason="not general enough", by="tester")
    names = [s.name for s in discover_skills()]
    assert "proposed-skill" not in names


def test_only_acceptance_makes_a_proposal_a_skill(home, monkeypatch):
    _standin_safety(monkeypatch)
    _propose(home)
    target = sp.accept("proposed-skill", by="tester")
    skills = {s.name: s for s in discover_skills()}
    assert "proposed-skill" in skills
    assert skills["proposed-skill"].body == "# Proposed\nbody"


def test_session_skills_hides_a_pending_proposal(home, monkeypatch):
    """The ``skill`` tool's own view, through the api_client wrapper.

    That wrapper is what both the advertised tool surface and the
    executor read, so proving it here covers the public path a proposal
    could otherwise sneak through.
    """
    _standin_safety(monkeypatch)
    from delfin.agent import api_client
    _propose(home)
    perms = types.SimpleNamespace(workspace=home, skip_skill_discovery=False)
    names = [s.name for s in api_client._session_skills(perms)]
    assert "proposed-skill" not in names
    sp.accept("proposed-skill", by="tester")
    names = [s.name for s in api_client._session_skills(perms)]
    assert "proposed-skill" in names


def test_a_workspace_does_not_change_that(home, monkeypatch, tmp_path):
    """Project-scoped skill dirs must not pick proposals up either."""
    _standin_safety(monkeypatch)
    _propose(home)
    names = [s.name for s in discover_skills(tmp_path)]
    assert "proposed-skill" not in names
