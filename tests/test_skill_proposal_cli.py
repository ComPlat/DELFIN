"""`delfin-agent skills proposals ...` — the human side of a proposal.

A proposal only becomes a skill through a person: the CLI lists what is
waiting, `show` gives the full text, evidence and safety findings (the
preview the decision needs), `accept` activates, `reject --reason` turns
it away. All English output; nothing is deleted or overwritten.
"""
import sys
import types

import pytest

from delfin.agent import skill_proposals as sp


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    mod = types.ModuleType("delfin.agent.skill_safety")
    mod.check = lambda text: []
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", mod)
    return tmp_path


def _propose(name="proposal-a", text="# Heading\nBody text."):
    return sp.propose(name, text,
                      evidence=[sp.Evidence(
                          kind="test", ref="tests/test_a.py::test_b",
                          detail="green on 6eb3d906")],
                      source="test")


def _run_cli(argv):
    from delfin.agent import cli
    ns = cli.build_parser().parse_args(argv)
    return ns.func(ns)


def test_proposals_lists_pending_and_blocked(home, monkeypatch):
    _propose()
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", None)
    _propose("proposal-b", "# H\nblocked one")
    monkeypatch.delitem(sys.modules, "delfin.agent.skill_safety")
    import importlib
    sys.modules.pop("delfin.agent.skill_safety", None)
    # restore the clean stand-in for the remaining assertions
    mod = types.ModuleType("delfin.agent.skill_safety")
    mod.check = lambda text: []
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", mod)

    class FakeOut:
        def __init__(self):
            self.parts = []
        def write(self, s):
            self.parts.append(s)

    out = FakeOut()
    monkeypatch.setattr(sys, "stdout", out)
    assert _run_cli(["skills", "proposals"]) == 0
    text = "".join(out.parts)
    assert "proposal-a" in text and "pending" in text
    assert "proposal-b" in text and "blocked" in text


def test_show_prints_text_evidence_and_findings(home, monkeypatch):
    _propose()
    class FakeOut:
        def __init__(self):
            self.parts = []
        def write(self, s):
            self.parts.append(s)
    out = FakeOut()
    monkeypatch.setattr(sys, "stdout", out)
    assert _run_cli(["skills", "show", "proposal-a"]) == 0
    text = "".join(out.parts)
    assert "Body text." in text
    assert "tests/test_a.py::test_b" in text
    assert "green on 6eb3d906" in text


def test_accept_activates_the_skill(home, monkeypatch):
    _propose()
    assert _run_cli(["skills", "accept", "proposal-a"]) == 0
    from delfin.agent.skills import discover_skills
    assert "proposal-a" in [s.name for s in discover_skills()]


def test_reject_needs_a_reason_and_moves_to_rejected(home, monkeypatch):
    _propose()
    assert _run_cli(["skills", "reject", "proposal-a"]) != 0
    assert _run_cli(["skills", "reject", "proposal-a",
                     "--reason", "too narrow"]) == 0
    assert sp.get_proposal("proposal-a") is None
    assert sp.list_proposals(status="rejected")[0].name == "proposal-a"


def test_repl_skills_proposals_lists_and_show_shows(home, monkeypatch):
    _propose()
    from delfin.agent import repl_commands as rc
    ctx = types.SimpleNamespace(workspace=home)
    res = rc._skills(ctx, "proposals")
    assert "proposal-a" in res.output
    res = rc._skills(ctx, "proposals proposal-a")
    assert "Body text." in res.output
    assert "tests/test_a.py::test_b" in res.output
