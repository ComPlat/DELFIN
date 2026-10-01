"""skill_patch: proposing a new skill version as a patch (alt -> neu).

The ACTIVE skill is never modified — a patch only creates a proposal
(with base_version = the version being patched). Acceptance and the
_versions/ archive are Paket 1's (skill_proposals.accept) calling
archive_previous_version from here.
"""

import pathlib

import pytest

from delfin.agent import skill_patch
from delfin.agent.skill_patch import (
    SkillPatchError,
    archive_previous_version,
    propose_skill_patch,
)


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


def _write(ws, name, text):
    p = ws / ".delfin" / "skills" / f"{name}.md"
    p.write_text(text, encoding="utf-8")
    return p


class FakeProposals:
    """Stand-in for delfin.agent.skill_proposals (Paket 1, contract)."""

    def __init__(self):
        self.calls = []

    def propose(self, name, text, *, evidence, source, base_version=""):
        self.calls.append({
            "name": name, "text": text, "evidence": evidence,
            "source": source, "base_version": base_version,
        })
        return type("P", (), {"name": name, "status": "pending"})()


EV = [{"kind": "test", "ref": "tests/test_skill_patch.py::t"}]


def test_a_patch_proposes_the_next_version_and_leaves_the_active_skill_alone(
        only_this_workspace):
    src = _write(only_this_workspace, "tune",
                 "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    before = src.read_text(encoding="utf-8")
    fake = FakeProposals()
    out = propose_skill_patch("tune", "step one", "step one, then verify",
                              "verification was missing", EV,
                              workspace=only_this_workspace, proposals=fake)
    assert src.read_text(encoding="utf-8") == before, "active skill changed"
    assert out["base_version"] == "1"
    assert out["new_version"] == "2"
    assert len(fake.calls) == 1
    call = fake.calls[0]
    assert call["base_version"] == "1"
    assert "version: 2" in call["text"]
    assert "step one, then verify" in call["text"]
    assert call["source"] == "skill_patch"


def test_a_patch_on_a_skill_without_a_version_starts_at_two(
        only_this_workspace):
    _write(only_this_workspace, "tune",
           "---\nname: tune\n---\n\n# Tune\n\nstep one\n")
    fake = FakeProposals()
    out = propose_skill_patch("tune", "step one", "step two",
                              "why", EV, workspace=only_this_workspace,
                              proposals=fake)
    assert out["base_version"] == "1"
    assert out["new_version"] == "2"
    assert fake.calls[0]["base_version"] == "1"


def test_a_patch_on_a_missing_skill_is_an_error_not_a_new_proposal(
        only_this_workspace):
    fake = FakeProposals()
    with pytest.raises(SkillPatchError, match="no active skill"):
        propose_skill_patch("ghost", "a", "b", "why", EV,
                            workspace=only_this_workspace, proposals=fake)
    assert fake.calls == []


def test_a_patch_without_evidence_is_refused(only_this_workspace):
    _write(only_this_workspace, "tune",
           "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    fake = FakeProposals()
    with pytest.raises(SkillPatchError, match="[Ee]vidence"):
        propose_skill_patch("tune", "step one", "step two", "why", [],
                            workspace=only_this_workspace, proposals=fake)
    assert fake.calls == []


def test_an_ambiguous_old_string_is_an_error_with_line_numbers(
        only_this_workspace):
    _write(only_this_workspace, "tune",
           "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep\nstep\n")
    fake = FakeProposals()
    with pytest.raises(SkillPatchError, match="twice|2 matches|match"):
        propose_skill_patch("tune", "step", "step prime", "why", EV,
                            workspace=only_this_workspace, proposals=fake)
    assert fake.calls == []


def test_a_non_matching_old_string_names_where_it_almost_matched(
        only_this_workspace):
    _write(only_this_workspace, "tune",
           "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    fake = FakeProposals()
    with pytest.raises(SkillPatchError):
        propose_skill_patch("tune", "step three", "x", "why", EV,
                            workspace=only_this_workspace, proposals=fake)
    assert fake.calls == []


# ---------------------------------------------------------------------------
# archive_previous_version — the versioning branch accept() calls (Paket 1)
# ---------------------------------------------------------------------------

def test_archive_moves_the_previous_text_under_versions_and_never_overwrites(
        only_this_workspace):
    src = _write(only_this_workspace, "tune",
                 "---\nname: tune\nversion: 1\n---\n\n# Tune\n\nstep one\n")
    where = archive_previous_version(src, "1")
    assert where.is_file()
    assert "step one" in where.read_text(encoding="utf-8")
    assert where.parent.name == "_versions"
    # Second archive of the same version must not overwrite.
    src.write_text("---\nname: tune\nversion: 1\n---\n\n# Tune\n\nother\n",
                   encoding="utf-8")
    with pytest.raises(SkillPatchError, match="exists"):
        archive_previous_version(src, "1")
    assert "step one" in where.read_text(encoding="utf-8")


def test_archive_rejects_a_versionless_skill(only_this_workspace):
    src = _write(only_this_workspace, "tune", "# Tune\n\nstep one\n")
    with pytest.raises(SkillPatchError, match="[Vv]ersion"):
        archive_previous_version(src, "")
