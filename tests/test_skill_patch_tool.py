"""skill_propose_patch through the real tool execution path (Paket 4).

Calls the executor the way DELFIN's chat loop does —
``_DocToolExecutor.execute("skill_propose_patch", ...)``, permissions
and all — and follows the proposal through the accept() branch Paket 1
owns (simulated with a stand-in that calls archive_previous_version,
per the contract): patch -> v2 pending, active skill untouched;
accept -> v2 active, v1 under _versions/.
"""

import json
import pathlib
import sys
import types

import pytest

from delfin.agent import api_client as A
from delfin.agent import skill_patch


@pytest.fixture
def only_this_workspace(tmp_path, monkeypatch):
    """A workspace whose skills are the ONLY ones discoverable."""
    empty = tmp_path / "no_pack"
    empty.mkdir()
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(skill_patch.skills, "_PACK_SKILLS_DIR", empty)
    monkeypatch.setattr(pathlib.Path, "home", classmethod(lambda cls: home))
    ws = tmp_path / "projekt"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    return ws


class StandInProposals:
    """Paket 1 stand-in: propose + accept with the version archive."""

    def __init__(self, home):
        self.home = pathlib.Path(home)
        self.proposals = []
        self.accepted = []

    def propose(self, name, text, *, evidence, source, base_version=""):
        rec = {"name": name, "text": text, "evidence": evidence,
               "source": source, "base_version": base_version,
               "status": "pending"}
        self.proposals.append(rec)
        return types.SimpleNamespace(name=name, status="pending",
                                     **{"base_version": base_version})

    def accept(self, name, *, by):
        """What Paket 1's accept() will do for a patched proposal."""
        rec = next(r for r in self.proposals if r["name"] == name)
        assert rec["status"] == "pending"
        # Version branch: hand the previous text to the archive first.
        src = pathlib.Path(rec["active_source"])
        archived = skill_patch.archive_previous_version(
            src, rec["base_version"])
        # Then the new version becomes the active skill.
        src.write_text(rec["text"], encoding="utf-8")
        rec["status"] = "accepted"
        self.accepted.append((name, archived))
        return archived


@pytest.fixture
def store(only_this_workspace, monkeypatch, tmp_path):
    """A fake delfin.agent.skill_proposals visible to the lazy import."""
    fake = StandInProposals(tmp_path / "home")
    mod = types.ModuleType("delfin.agent.skill_proposals")
    mod.propose = fake.propose
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_proposals", mod)
    # Direct stand-in for _load_proposals-based default use.
    monkeypatch.setattr(skill_patch, "_load_proposals",
                        lambda: types.SimpleNamespace(propose=fake.propose))
    return fake


EV = [{"kind": "test", "ref": "tests/test_skill_patch_tool.py::t"}]


def _write(ws, name, text):
    p = ws / ".delfin" / "skills" / f"{name}.md"
    p.write_text(text, encoding="utf-8")
    return p


def _exec(perms, **arguments):
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    return json.loads(ex.execute(
        "skill_propose_patch", {
            "name": "tune", "old": "step one", "new": "step one, verified",
            "reason": "verification was missing",
            "evidence": EV,
            **arguments,
        }, perms))


def test_patch_via_the_executor_leaves_v1_active_and_proposes_v2(
        only_this_workspace, store):
    src = _write(only_this_workspace, "tune",
                 "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    store.proposals  # stand-in wired
    perms = A.KitToolPermissions(workspace=only_this_workspace)
    out = _exec(perms)
    assert out["status"] == "ok", out
    assert out["new_version"] == "2"
    assert out["base_version"] == "1"
    assert out["proposal_status"] == "pending"
    # The active skill is untouched.
    assert "step one, verified" not in src.read_text(encoding="utf-8")
    assert "version: 1" in src.read_text(encoding="utf-8")
    # The proposal carries the patched text and the base version.
    assert len(store.proposals) == 1
    rec = store.proposals[0]
    assert rec["base_version"] == "1"
    assert "version: 2" in rec["text"]
    assert "step one, verified" in rec["text"]
    rec["active_source"] = src  # what accept() will need


def test_accept_activates_v2_and_archives_v1(
        only_this_workspace, store):
    src = _write(only_this_workspace, "tune",
                 "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    perms = A.KitToolPermissions(workspace=only_this_workspace)
    _exec(perms)
    store.proposals[0]["active_source"] = src
    # accept() as Paket 1 will run it: archive v1, activate v2.
    store.accept("tune", by="human")
    now = src.read_text(encoding="utf-8")
    assert "version: 2" in now
    assert "step one, verified" in now
    archived = store.accepted[0][1]
    assert archived.parent.name == "_versions"
    assert archived.name == "1.md"
    assert "step one\n" in archived.read_text(encoding="utf-8")
    assert "version: 1" in archived.read_text(encoding="utf-8")


def test_the_executor_reports_patch_errors_from_the_engine(
        only_this_workspace, store):
    _write(only_this_workspace, "tune",
           "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    perms = A.KitToolPermissions(workspace=only_this_workspace)
    out = _exec(perms, old="no such line")
    assert "error" in out, out
    assert store.proposals == []


def test_plan_mode_refuses_and_proposes_nothing(
        only_this_workspace, store):
    _write(only_this_workspace, "tune",
           "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    perms = A.KitToolPermissions(workspace=only_this_workspace, mode="plan")
    out = _exec(perms)
    assert out["error"].startswith("plan mode"), out
    assert store.proposals == []
