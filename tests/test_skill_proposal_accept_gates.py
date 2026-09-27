"""accept() gates on evidence and archives the previous version.

The acceptance rule (coordinator decision, wave "Lernen"):
- at least ONE evidence entry must verify ok, AND
- no mandatory evidence entry (chemistry: kind calc or test) may be
  invalid -- a proposal whose only chemistry evidence is invalid is not
  acceptable;
- a proposal with base_version != "" archives the previously active
  skill text BEFORE the new one is activated (package 4,
  archive_previous_version); a SkillPatchError is a clean refusal, not
  a crash;
- a team-pulled proposal (source="team", evidence=[] placeholder in the
  pull path) is not acceptable on empty evidence either -- it must go
  through the same verify_evidence gate as any other proposal.

Until packages 4 (skill_patch) and 8 (evidence) are merged onto this
branch, the tests import them from their branches via a stand-in
installed into sys.modules by the fixture -- the stand-in implements
the exact published signatures.
"""
import sys
import types

import pytest

from delfin.agent import skill_proposals as sp


class _StandinArchive:
    """Stand-in for delfin.agent.skill_patch.archive_previous_version
    (package 4, commit 07f82f22). Records the call, moves nothing:
    the real one is exercised in s14's own tests; here we pin the
    accept() side of the contract."""
    def __init__(self):
        self.calls = []
        self.error = None

    def __call__(self, skill_source, base_version):
        self.calls.append((str(skill_source), base_version))
        if self.error is not None:
            raise _SkillPatchError(self.error)
        from pathlib import Path
        return Path(str(skill_source)) / "_versions" / f"{base_version}.md"


class _SkillPatchError(Exception):
    pass


class _StandinVerify:
    """Stand-in for delfin.agent.evidence.verify_evidence (package 8,
    commit 466c4c82): (ok, detail) per entry, never raises."""
    def __init__(self):
        self.results = {}  # ref -> (ok, detail)

    def __call__(self, evidence, *, workspace=None, runs=None,
                 list_jobs=None):
        ref = getattr(evidence, "ref", "") or ""
        return self.results.get(ref, (True, "verified"))


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    # clean safety checker
    safe = types.ModuleType("delfin.agent.skill_safety")
    safe.check = lambda text: []
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", safe)
    # ``from . import skill_safety`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "skill_safety",
                        safe, raising=False)
    # stand-ins for packages 4 and 8
    arch = _StandinArchive()
    patch_mod = types.ModuleType("delfin.agent.skill_patch")
    patch_mod.archive_previous_version = arch
    patch_mod.SkillPatchError = _SkillPatchError
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_patch", patch_mod)
    # ``from . import skill_patch`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "skill_patch",
                        patch_mod, raising=False)
    ver = _StandinVerify()
    ev_mod = types.ModuleType("delfin.agent.evidence")
    ev_mod.verify_evidence = ver
    monkeypatch.setitem(sys.modules, "delfin.agent.evidence", ev_mod)
    # ``from . import evidence`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "evidence",
                        ev_mod, raising=False)
    return types.SimpleNamespace(root=tmp_path, archive=arch, verify=ver)


def _ev(kind="test", ref="tests/test_x.py::t"):
    return sp.Evidence(kind=kind, ref=ref)


def test_accept_with_all_evidence_ok(home):
    p = sp.propose("ok-skill", "# Fine\n", evidence=[_ev()], source="t")
    target = sp.accept(p.name, by="tester")
    assert target.exists()


def test_accept_refuses_when_no_evidence_verifies(home):
    home.verify.results["tests/test_x.py::t"] = (False, "run not in ledger")
    p = sp.propose("bad-ev", "# Fine\n", evidence=[_ev()], source="t")
    with pytest.raises(ValueError, match="evidence"):
        sp.accept(p.name, by="tester")
    # the proposal is untouched -- a refusal is not a deletion
    assert sp.get_proposal(p.name).status == "pending"


def test_accept_refuses_when_a_mandatory_evidence_is_invalid(home):
    good = _ev(ref="tests/test_ok.py::t")
    bad_calc = _ev(kind="calc", ref="/calc/nonexistent")
    home.verify.results["/calc/nonexistent"] = (False, "no such folder")
    p = sp.propose("mixed-ev", "# Fine\n",
                   evidence=[good, bad_calc], source="t")
    with pytest.raises(ValueError, match="calc"):
        sp.accept(p.name, by="tester")


def test_a_team_pull_without_evidence_is_not_acceptable(home):
    p = sp.propose("team-skill", "# Pulled\n",
                   evidence=[_ev(ref="team-placeholder")], source="team")
    home.verify.results["team-placeholder"] = (False,
                                               "team pull has no evidence")
    with pytest.raises(ValueError, match="evidence"):
        sp.accept(p.name, by="tester")


def test_accept_archives_the_previous_version_before_activating(home):
    p = sp.propose("versioned", "# v2 text\n", evidence=[_ev()], source="t",
                   base_version="1.2")
    # a live skill with that name and version must already exist for the
    # archive step to have something to archive
    skills_root = sp._proposals_dir().parent
    live = skills_root / "versioned"
    live.mkdir(parents=True)
    (live / "SKILL.md").write_text("--- version: 1.2 ---\n# v1\n",
                                   encoding="utf-8")
    target = sp.accept(p.name, by="tester")
    assert home.archive.calls, "the previous version must be archived first"
    skill_source, base_version = home.archive.calls[0]
    assert skill_source == str(live / "SKILL.md")
    assert base_version == "1.2"
    # and the archive happened BEFORE the move: the source path recorded
    # is the live one, and the proposal folder is gone from _proposals
    assert not sp._proposals_dir().joinpath(p.name).exists()


def test_second_patch_on_the_same_skill_keeps_every_record(home):
    # Operator review of 2ea09a02: accepting a SECOND patch on the same
    # skill crashed on os.replace(d, target / "_accepted_proposal")
    # because that directory already existed from the first patch --
    # AFTER the SKILL.md replace had already run, leaving a half state.
    skills_root = sp._proposals_dir().parent
    live = skills_root / "double-patched"
    live.mkdir(parents=True)
    (live / "SKILL.md").write_text("--- version: 1.1 ---\n# v1\n",
                                   encoding="utf-8")
    p2 = sp.propose("double-patched", "# v2\n", evidence=[_ev()],
                    source="t", base_version="1.1")
    sp.accept(p2.name, by="tester")
    p3 = sp.propose("double-patched", "# v3\n", evidence=[_ev()],
                    source="t", base_version="1.2")
    target = sp.accept(p3.name, by="tester")
    # v3 is live
    assert "# v3" in (live / "SKILL.md").read_text(encoding="utf-8")
    # both acceptance records survived, one per version, none overwritten
    records = sorted(live.glob("_accepted_proposals/*/proposal.json"))
    assert len(records) == 2, f"expected 2 records, got {records}"
    # and both proposals are gone from _proposals
    assert not sp._proposals_dir().joinpath(p2.name).exists()
    assert not sp._proposals_dir().joinpath(p3.name).exists()


def test_a_failure_in_the_record_step_leaves_the_live_skill_untouched(home):
    # The second half of the operator finding: a step that fails must
    # never have touched the live SKILL.md. We force the record step to
    # fail by planting a colliding directory where the record would go.
    skills_root = sp._proposals_dir().parent
    live = skills_root / "clash-record"
    live.mkdir(parents=True)
    (live / "SKILL.md").write_text("--- version: 1.0 ---\n# v1\n",
                                   encoding="utf-8")
    before = (live / "SKILL.md").read_text(encoding="utf-8")
    p = sp.propose("clash-record", "# v2\n", evidence=[_ev()],
                   source="t", base_version="1.0")
    # occupied record slot: the record step must notice and refuse,
    # BEFORE the live text is replaced
    clash = live / "_accepted_proposals" / "1.0"
    clash.mkdir(parents=True)
    (clash / "proposal.json").write_text("{}",
                                         encoding="utf-8")
    with pytest.raises(ValueError, match="record"):
        sp.accept(p.name, by="tester")
    after = (live / "SKILL.md").read_text(encoding="utf-8")
    assert after == before, "a failed record step must not touch the skill"
    # the proposal is still pending in _proposals
    assert sp.get_proposal(p.name).status == "pending"


def test_a_skill_patch_error_is_a_clean_refusal(home):
    p = sp.propose("clashy", "# v2\n", evidence=[_ev()], source="t",
                   base_version="1.0")
    live = sp._proposals_dir().parent / "clashy"
    live.mkdir(parents=True)
    (live / "SKILL.md").write_text("# v1\n", encoding="utf-8")
    home.archive.error = "archive file for version 1.0 already exists"
    with pytest.raises(ValueError, match="1.0"):
        sp.accept(p.name, by="tester")
    # nothing was activated
    assert sp.get_proposal(p.name).status == "pending"
