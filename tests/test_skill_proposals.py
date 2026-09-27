"""Skill proposals: a self-written skill is only a PROPOSAL.

Invisible to discover_skills until a human accepts it; nothing without
evidence; nothing unsafe becomes a rule; never overwritten, never deleted.
"""
import sys
import types

import pytest

from delfin.agent import skill_proposals as sp


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    return tmp_path


def _fake_safety(monkeypatch, findings):
    mod = types.ModuleType("delfin.agent.skill_safety")
    calls = []
    def check(text):
        calls.append(text)
        return list(findings)
    mod.check = check
    mod.calls = calls
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", mod)
    return mod


def _ev(**kw):
    defaults = {"kind": "test", "ref": "tests/test_x.py::test_y"}
    defaults.update(kw)
    return sp.Evidence(**defaults)


def test_without_evidence_there_is_no_proposal(home):
    with pytest.raises(ValueError):
        sp.propose("noop", "body", evidence=[], source="test")


def test_propose_creates_pending_proposal(home):
    p = sp.propose("draft-csv", "how to draft", evidence=[_ev()], source="test")
    assert p.status == "pending"
    assert p.findings == []
    got = sp.get_proposal("draft-csv")
    assert got is not None and got.text == "how to draft"
    assert sp.list_proposals(status="pending")[0].name == "draft-csv"


def test_safety_findings_block_a_proposal(home, monkeypatch):
    _fake_safety(monkeypatch, ["asks for a permission bypass"])
    p = sp.propose("bad", "bypass", evidence=[_ev()], source="test")
    assert p.status == "blocked"
    assert "permission bypass" in p.findings[0]
    assert sp.get_proposal("bad").status == "blocked"
    with pytest.raises(ValueError):
        sp.accept("bad", by="tester")


def test_absent_safety_module_leaves_proposal_pending(home, monkeypatch):
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", None)
    p = sp.propose("plain", "body", evidence=[_ev()], source="test")
    assert p.status == "pending"


def test_name_conflict_gets_a_new_name_never_overwritten(home):
    a = sp.propose("dup", "first", evidence=[_ev()], source="t")
    b = sp.propose("dup", "second", evidence=[_ev()], source="t")
    assert b.name != a.name
    assert sp.get_proposal(a.name).text == "first"


def test_proposals_dir_is_private_and_files_owner_only(home):
    sp.propose("perms", "body", evidence=[_ev()], source="t")
    d = sp.PROPOSALS_DIR
    assert (d.stat().st_mode & 0o777) == 0o700
    assert (d / "perms" / "SKILL.md").stat().st_mode & 0o777 == 0o600
    assert (d / "perms" / "proposal.json").stat().st_mode & 0o777 == 0o600
    assert d.parent.name == "skills"
